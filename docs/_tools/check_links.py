"""Check internal links: the built site, and the Markdown as GitHub shows it.

    python docs/_tools/check_links.py [--site site]

Exit status 0 when every link resolves, 1 otherwise. No network access.

1. **The built site.** Every ``href``/``src`` in every HTML page that is not an
   external URL must name a file in the site, and every ``#fragment`` an
   ``id`` in the page it points at. This covers what ``mkdocs build --strict``
   does not: links inside rendered docstrings, the navigation and Material's
   own markup.

2. **The Markdown, as GitHub renders it.** The pages under ``docs/`` are also
   read on GitHub, and the repository's own README files link into them. For
   every relative link in those files the target must be tracked by git -- so
   this checkout and a fresh clone agree; links into gitignored outputs are
   listed as notes -- and an ``#anchor`` into a Markdown file must match one of its headings under
   GitHub's slug rules (lower case, punctuation dropped, spaces to hyphens,
   repeated hyphens kept). MkDocs slugs collapse repeated hyphens and drop
   non-ASCII letters, so an anchor that works in one place can fail in the
   other; this check holds both. Links to API module pages generated at build
   time (``docs/reference/api/<family>/<module>.md``) exist only in the site
   and are accepted here.

Archived case directories are frozen and are not checked.
"""

from __future__ import annotations

import argparse
import re
import sys
import unicodedata
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import hooks  # noqa: E402

ROOT = hooks.ROOT
ARCHIVED = ("1d", "2d-planar", "2d-wedge", "caseFoamEx", "wake-cylinder")


# --------------------------------------------------------------------------- #
# 1. the built site
# --------------------------------------------------------------------------- #

class _Page(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.ids: set[str] = set()
        self.links: list[tuple[str, str]] = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        for key in ("id", "name"):
            if attrs.get(key) and (key == "id" or tag == "a"):
                self.ids.add(attrs[key])
        for key in ("href", "src"):
            if attrs.get(key) is not None and tag in ("a", "img", "script", "link", "source"):
                self.links.append((tag, attrs[key]))


def _parse(path: Path) -> _Page:
    page = _Page()
    page.feed(path.read_text(encoding="utf-8", errors="replace"))
    return page


def check_site(site: Path) -> list[str]:
    problems = []
    pages = {p.resolve(): _parse(p) for p in site.rglob("*.html")}
    for path, page in pages.items():
        for tag, link in page.links:
            parts = urlsplit(link)
            if parts.scheme or link.startswith(("//", "mailto:", "javascript:", "data:")):
                continue
            if not parts.path and not parts.fragment:
                continue
            if not parts.path:
                target = path
            elif parts.path.startswith("/"):
                # 404.html is served at any URL, so MkDocs writes its links
                # relative to the site root
                target = (site / unquote(parts.path).lstrip("/")).resolve()
            else:
                target = (path.parent / unquote(parts.path)).resolve()
            if target.is_dir():
                target = target / "index.html"
            where = path.relative_to(site.resolve()).as_posix()
            if not target.exists():
                problems.append(f"site/{where}: {tag} '{link}' -> missing {target.name}")
                continue
            if parts.fragment and target.suffix == ".html":
                ids = pages[target].ids if target in pages else _parse(target).ids
                if unquote(parts.fragment) not in ids:
                    problems.append(f"site/{where}: '{link}' -> no id '{parts.fragment}'")
    return problems


# --------------------------------------------------------------------------- #
# 2. Markdown as GitHub shows it
# --------------------------------------------------------------------------- #

_FENCE = re.compile(r"^\s*(`{3,}|~{3,})")
_LINK = re.compile(r"(?<!!)\[[^\]\n]*\]\((<[^>\n]+>|[^)\s]+)(?:\s+\"[^\"]*\")?\)")
_CODE = re.compile(r"(`+)[^\n]*?(?<!`)\1(?!`)")
_HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")


def _prose_lines(text: str):
    fence = None
    for number, line in enumerate(text.splitlines(), 1):
        opening = _FENCE.match(line)
        if fence is None and opening:
            fence = opening.group(1)[0]
            continue
        if fence is not None:
            if line.strip().startswith(fence * 3) and not line.strip().strip(fence):
                fence = None
            continue
        yield number, line


def github_slugs(path: Path) -> set[str]:
    """Heading anchors GitHub generates for a Markdown file."""
    seen: dict[str, int] = {}
    slugs = set()
    for _, line in _prose_lines(path.read_text(encoding="utf-8")):
        match = _HEADING.match(line)
        if not match:
            continue
        text = re.sub(r"\{\s*#[^}]*\}\s*$", "", match.group(2))       # attr_list ids
        text = re.sub(r"!?\[([^\]]*)\]\([^)]*\)", r"\1", text)       # links -> text
        text = text.replace("`", "")
        text = re.sub(r"<[^>]+>", "", text)
        slug = "".join(ch for ch in text.lower()
                       if ch in "-_ " or unicodedata.category(ch)[0] in "LN")
        slug = slug.replace(" ", "-")
        count = seen.get(slug, 0)
        seen[slug] = count + 1
        slugs.add(slug if count == 0 else f"{slug}-{count}")
    return slugs


def markdown_files() -> list[Path]:
    files = [ROOT / "README.md", ROOT / "cases" / "ARCHIVE.md",
             ROOT / "tests" / "regression" / "README.md"]
    for path in sorted((ROOT / "cases").rglob("*.md")):
        rel = path.relative_to(ROOT / "cases").parts
        if rel[0] in ARCHIVED or any(p.startswith(("Cases", "Interface", "results", "runs")) for p in rel):
            continue
        files.append(path)
    for path in sorted((ROOT / "docs").rglob("*.md")):
        if "_tools" not in path.parts:
            files.append(path)
    return [f for f in dict.fromkeys(files) if f.is_file()]


def generated_pages() -> set[Path]:
    pages = set()
    for entries in hooks.api_pages().values():
        for _, page in entries:
            pages.add((ROOT / "docs" / page).resolve())
    return pages


def _ignored(rel: str) -> bool:
    import subprocess
    return subprocess.run(["git", "check-ignore", "-q", "--", rel], cwd=ROOT).returncode == 0


def check_markdown(notes: list[str] | None = None) -> list[str]:
    """Broken links in the repository's Markdown, judged against what git tracks.

    A target must be tracked (or, without git, exist), so the answer is the same
    in this checkout and in a fresh clone. A link into a gitignored path names a
    generated output that exists only where it was produced; it cannot be
    checked here and is appended to ``notes`` instead of failing.
    """
    problems = []
    generated = generated_pages()
    tracked = hooks.tracked_paths()
    slug_cache: dict[Path, set[str]] = {}
    for path in markdown_files():
        text = path.read_text(encoding="utf-8")
        for number, line in _prose_lines(text):
            line = _CODE.sub("", line)
            for match in _LINK.finditer(line):
                target = match.group(1).strip("<>")
                parts = urlsplit(target)
                if parts.scheme or target.startswith(("//", "mailto:")):
                    continue
                where = f"{path.relative_to(ROOT).as_posix()}:{number}"
                if not parts.path:
                    file = path
                else:
                    file = (path.parent / unquote(parts.path)).resolve()
                    if file in generated:
                        continue
                    try:
                        rel = file.relative_to(ROOT).as_posix()
                    except ValueError:
                        problems.append(f"{where}: '{target}' -> outside the repository")
                        continue
                    if tracked is not None:
                        present = rel in tracked[0] or rel in tracked[1]
                    else:
                        present = file.exists()
                    if not present:
                        if tracked is not None and _ignored(rel):
                            if notes is not None:
                                notes.append(f"{where}: '{target}' -> gitignored output, "
                                             f"present only where it was produced")
                            continue
                        problems.append(f"{where}: '{target}' -> no such tracked path")
                        continue
                if parts.fragment and file.suffix == ".md":
                    slugs = slug_cache.setdefault(file, github_slugs(file))
                    if parts.fragment not in slugs:
                        problems.append(f"{where}: '#{parts.fragment}' is not a heading "
                                        f"anchor of {file.name} on GitHub")
    return problems


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--site", type=Path, default=ROOT / "site",
                        help="the built site (default: ./site)")
    parser.add_argument("--markdown-only", action="store_true",
                        help="check only the Markdown, not a built site")
    args = parser.parse_args(argv)
    problems = []
    if not args.markdown_only:
        if not args.site.is_dir():
            print(f"check_links: no site at {args.site}; build it first", file=sys.stderr)
            return 1
        problems += check_site(args.site)
    notes: list[str] = []
    problems += check_markdown(notes)
    if notes:
        print(f"check_links: {len(notes)} link(s) into gitignored outputs, not checkable here:")
        for note in notes:
            print(f"  note: {note}")
    if problems:
        print(f"check_links: {len(problems)} broken link(s)")
        for problem in problems:
            print(f"  {problem}")
        return 1
    print("check_links: every internal link and anchor resolves"
          + ("" if args.markdown_only else f" (site: {args.site})"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

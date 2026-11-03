"""MkDocs hooks for the rarefiedPlume documentation build.

Two jobs, both deterministic and both read-only with respect to the repository:

1. **API pages.** Every module under ``plumetools/`` gets a reference page. A
   module with a hand-written page at ``docs/reference/api/<family>/<name>.md``
   keeps it; every other module gets a generated stub that renders its
   docstrings through mkdocstrings. The modules are found by walking the source
   tree, never by importing it, and the navigation is built from the same walk,
   so an added or removed module changes the site without anyone editing a
   list. ``docs/_tools/check_api.py`` checks the hand-written side.

2. **Repository links.** Pages link to source files with ordinary relative
   paths (``../../cases/cai2012/study.yaml``), which work when the Markdown is
   read on GitHub. Those targets are outside ``docs/``, so MkDocs cannot serve
   them; this hook checks that each one is tracked by git (falling back to
   "exists" where git is unavailable) and rewrites it to the file's GitHub URL.
   A link to a missing or untracked path -- a local result, a build product --
   is a warning, and ``mkdocs build --strict`` turns it into a failure, so a
   link that works only in one checkout fails there rather than in CI.

Set ``DOCS_SOURCE_REF`` to link to a commit or branch other than ``main``.
"""

from __future__ import annotations

import logging
import os
import posixpath
import re
from pathlib import Path

from mkdocs.structure.files import File

log = logging.getLogger("mkdocs.hooks.rarefiedplume")

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "plumetools"
API_DIR = "reference/api"
SOURCE_REF = os.environ.get("DOCS_SOURCE_REF", "main")

# Families, in navigation order. Top-level modules are assigned explicitly.
FAMILIES = {
    "legacy": "Legacy source flow (3d-inflow)",
    "cai2012": "cai2012",
    "cai2012ext": "cai2012ext (Kn extension)",
    "markelov1999": "markelov1999",
    "nozzleplume": "nozzleplume",
    "hybridplume": "hybridplume",
    "virtualcollision": "virtualcollision",
    "rarefaction": "rarefaction",
    "viz": "viz",
    "foamio": "foamio",
    "mesh": "mesh",
    "shared": "Shared workflow modules",
}
TOP_LEVEL = {
    "config": "legacy",
    "constants": "legacy",
    "geometry": "legacy",
    "inflow": "legacy",
    "sourceflow": "legacy",
    "species": "legacy",
    "verify_mesh": "legacy",
    "casemanifest": "shared",
    "iteration": "shared",
}


def discover_modules(package: Path = PACKAGE) -> list[str]:
    """Every module's dotted name, from the files on disk, sorted."""
    names = []
    for path in sorted(package.rglob("*.py")):
        parts = path.relative_to(package.parent).with_suffix("").parts
        if "__pycache__" in parts:
            continue
        if parts[-1] == "__init__":
            parts = parts[:-1]
        names.append(".".join(parts))
    return sorted(names)


def page_of(module: str) -> tuple[str, str] | None:
    """``(family, docs path)`` of a module's page; None for the root package.

    A subpackage's own ``__init__`` is documented on its family index page.
    """
    parts = module.split(".")
    if len(parts) == 1:
        return None
    if len(parts) == 2:
        name = parts[1]
        if name in TOP_LEVEL:
            family = TOP_LEVEL[name]
            return family, f"{API_DIR}/{family}/{name}.md"
        if name in FAMILIES:
            return name, f"{API_DIR}/{name}/index.md"
        raise ValueError(f"{module}: a new top-level module; assign it a family in docs/_tools/hooks.py")
    family = parts[1]
    if family not in FAMILIES:
        raise ValueError(f"{module}: a new subpackage; add it to FAMILIES in docs/_tools/hooks.py")
    return family, f"{API_DIR}/{family}/{'.'.join(parts[2:])}.md"


def source_path(module: str) -> str:
    rel = Path(*module.split("."))
    if (ROOT / rel).is_dir():
        return (rel / "__init__.py").as_posix()
    return rel.with_suffix(".py").as_posix()


def github_url(repo_url: str, rel: str, is_dir: bool = False) -> str:
    kind = "tree" if is_dir else "blob"
    return f"{repo_url.rstrip('/')}/{kind}/{SOURCE_REF}/{rel}"


def _stub(module: str, family: str, repo_url: str) -> str:
    entry = ""
    if module.endswith(".__main__"):
        package = module[: -len(".__main__")]
        entry = (f"\n\nThis is the `python -m {package}` entry point; its commands are "
                 f"described in the [CLI reference](../../cli.md).")
    return (
        f"# `{module}`\n\n"
        f"<!-- generated at build time by docs/_tools/hooks.py; add a page at this path to replace it -->\n\n"
        f'<p class="api-meta">Family: <a href="../">{FAMILIES[family]}</a> · '
        f'Source: <a href="{github_url(repo_url, source_path(module))}"><code>{source_path(module)}</code></a></p>'
        f"{entry}\n\n"
        f"::: {module}\n"
        f"    options:\n"
        f"      show_root_heading: false\n"
    )


def api_pages() -> dict[str, list[tuple[str, str]]]:
    """``{family: [(module, docs path), ...]}`` for every module below a family index."""
    pages: dict[str, list[tuple[str, str]]] = {family: [] for family in FAMILIES}
    for module in discover_modules():
        found = page_of(module)
        if found is None:
            continue
        family, path = found
        if not path.endswith("/index.md"):
            pages[family].append((module, path))
    return pages


# --------------------------------------------------------------------------- #
# hooks
# --------------------------------------------------------------------------- #

def on_config(config, **kwargs):
    """Fill the "API reference" navigation section from the source tree."""
    pages = api_pages()

    def fill(items):
        for item in items:
            if isinstance(item, dict):
                for title, value in item.items():
                    if title == "API reference" and isinstance(value, list):
                        for family, title_ in FAMILIES.items():
                            entries = [f"{API_DIR}/{family}/index.md"]
                            for module, path in pages[family]:
                                short = module.split(".", 2)[-1] if family not in ("legacy", "shared") else module.split(".")[-1]
                                entries.append({short: path})
                            value.append({title_: entries})
                        return True
                    if isinstance(value, list) and fill(value):
                        return True
        return False

    if not fill(config["nav"]):
        log.warning("mkdocs.yml has no 'API reference' section to fill")
    return config


def on_files(files, config, **kwargs):
    """Add a generated page for every module without a hand-written one."""
    repo_url = config.get("repo_url") or ""
    docs_dir = Path(config["docs_dir"])
    for family, entries in api_pages().items():
        for module, path in entries:
            if (docs_dir / path).exists():
                continue
            files.append(File.generated(config, path, content=_stub(module, family, repo_url)))
    return files


_FENCE_OPEN = re.compile(r"^(\s*)(`{3,}|~{3,})")
_LINK = re.compile(r"(?<!!)(\[[^\]\n]*\]\()(<[^>\n]+>|[^)\s]+)((?:\s+\"[^\"]*\")?\))")
_REF_DEF = re.compile(r"^(\s{0,3}\[[^\]]+\]:\s*)(\S+)(.*)$")
_CODE_SPAN = re.compile(r"(`+)[^\n]*?(?<!`)\1(?!`)")


_TRACKED: tuple | None = None


def tracked_paths() -> tuple | None:
    """``(files, directories)`` git tracks, as POSIX paths; None without git.

    A link must name something in the repository, not something that happens to
    exist in one checkout: a local result or build product would pass here and
    break in a fresh clone.
    """
    global _TRACKED
    if _TRACKED is None:
        import subprocess
        try:
            out = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True,
                                 check=True).stdout.decode("utf-8")
        except (OSError, subprocess.CalledProcessError):
            _TRACKED = (None,)
            return None
        files = {p for p in out.split("\0") if p}
        dirs = {str(Path(p).parent.as_posix()) for p in files}
        for d in list(dirs):
            while "/" in d:
                d = d.rsplit("/", 1)[0]
                dirs.add(d)
        _TRACKED = (files, dirs)
    return None if _TRACKED == (None,) else _TRACKED


def _rewrite_target(target: str, page_dir: Path, docs_dir: Path, repo_url: str, where: str) -> str:
    raw = target[1:-1] if target.startswith("<") else target
    if re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*:", raw) or raw.startswith(("#", "/")):
        return target
    path_part, hash_, anchor = raw.partition("#")
    if not path_part:
        return target
    resolved = (page_dir / path_part).resolve()
    try:
        resolved.relative_to(docs_dir)
        return target  # a page or asset of this site: MkDocs validates it
    except ValueError:
        pass
    try:
        rel = resolved.relative_to(ROOT).as_posix()
    except ValueError:
        log.warning(f"{where}: link '{raw}' points outside the repository")
        return target
    tracked = tracked_paths()
    if tracked is not None:
        files, dirs = tracked
        if rel not in files and rel not in dirs:
            why = "is not tracked by git (a local or generated file)" if resolved.exists() else "does not exist"
            log.warning(f"{where}: link '{raw}' names a repository path that {why}: {rel}")
            return target
        is_dir = rel in dirs
    elif not resolved.exists():
        log.warning(f"{where}: link '{raw}' names a repository path that does not exist: {rel}")
        return target
    else:
        is_dir = resolved.is_dir()
    return github_url(repo_url, rel, is_dir) + (hash_ + anchor if hash_ else "")


def rewrite_markdown(markdown: str, page_dir: Path, docs_dir: Path, repo_url: str, where: str) -> str:
    out = []
    fence = None
    buffer: list[str] = []

    def flush():
        if not buffer:
            return
        text = "\n".join(buffer)
        spans = []

        def keep(match):
            spans.append(match.group(0))
            return f"\x00{len(spans) - 1}\x00"

        text = _CODE_SPAN.sub(keep, text)
        text = _LINK.sub(lambda m: m.group(1) + _rewrite_target(m.group(2), page_dir, docs_dir, repo_url, where) + m.group(3), text)
        lines = []
        for line in text.split("\n"):
            ref = _REF_DEF.match(line)
            if ref:
                line = ref.group(1) + _rewrite_target(ref.group(2), page_dir, docs_dir, repo_url, where) + ref.group(3)
            lines.append(line)
        text = re.sub(r"\x00(\d+)\x00", lambda m: spans[int(m.group(1))], "\n".join(lines))
        out.append(text)
        buffer.clear()

    for line in markdown.split("\n"):
        opening = _FENCE_OPEN.match(line)
        if fence is None and opening:
            flush()
            fence = opening.group(2)[0] * len(opening.group(2))
            out.append(line)
        elif fence is not None:
            out.append(line)
            if line.strip().startswith(fence) and line.strip().strip(fence[0]) == "":
                fence = None
        else:
            buffer.append(line)
    flush()
    return "\n".join(out)


def on_page_markdown(markdown, page, config, files, **kwargs):
    docs_dir = Path(config["docs_dir"]).resolve()
    if page.file.abs_src_path:
        page_dir = Path(page.file.abs_src_path).resolve().parent
    else:
        page_dir = docs_dir / posixpath.dirname(page.file.src_uri)
    return rewrite_markdown(markdown, page_dir, docs_dir, config.get("repo_url") or "", page.file.src_uri)

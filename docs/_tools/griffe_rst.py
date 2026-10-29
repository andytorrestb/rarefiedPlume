"""Render the reStructuredText inside plumetools docstrings as Markdown.

A Griffe extension, used only by the documentation build. The package's
docstrings are Google-style sections written in reStructuredText: ``.. math::``
blocks, ``:math:`` and ``:func:`` roles, simple tables, section underlines,
literal blocks and a few admonitions. MkDocs renders Markdown, so without this
step those constructs reach the page as literal markup.

The translation happens in memory, after Griffe has loaded the package
statically and before mkdocstrings renders it. No source file is read through an
import and none is written: the docstrings in ``plumetools/`` are unchanged.

What is translated, and to what:

* ``.. math::`` -> a ``$$`` block, or an ``align``/``gather`` environment when
  the block has several rows (MathJax, through ``pymdownx.arithmatex``).
* ``:math:`x``` -> ``$x$``.
* ``:func:``, ``:class:``, ``:mod:``, ``:meth:``, ``:attr:``, ``:data:``,
  ``:exc:``, ``:obj:`` -> a cross-reference when the target is rendered in this
  site, otherwise inline code. Targets are resolved in the docstring's own
  scope, as Sphinx would.
* ``.. code-block:: lang`` and ``text::`` literal blocks -> fenced code.
* ``.. warning::`` / ``.. note::`` -> admonitions.
* Simple (``=====``) tables -> Markdown tables; grid tables -> preformatted text.
* Section underlines -> headings.

Anything else is left as written. ``docs/development/documentation.md`` lists
what this does not handle.
"""

from __future__ import annotations

import re

import griffe

_DIRECTIVE = re.compile(r"^(\s*)\.\. ([a-z-]+)::\s*(.*)$")
_SIMPLE_BORDER = re.compile(r"^(\s*)(=+)((?: +=+)+)\s*$")
_GRID_BORDER = re.compile(r"^\s*\+[-=+]+\+\s*$")
_UNDERLINE = re.compile(r"^(\s*)([-=~^])\2{2,}\s*$")
_FENCE = re.compile(r"^\s*(```|~~~)")
_ROLE = re.compile(r":(func|meth|class|mod|attr|data|exc|obj):`([^`]+)`")
_MATH_ROLE = re.compile(r":math:`([^`]*)`")
_INLINE = re.compile(
    r"(?P<code2>``(?s:.+?)``)"  # literal text may wrap onto the next line
    r"|(?P<math>:math:`[^`]*`)"
    r"|(?P<role>:(?:func|meth|class|mod|attr|data|exc|obj):`[^`]+`)"
    r"|(?P<code1>(?<!`)`(?!`)[^`\n]+`)"
)
_ADMONITIONS = {"warning": "warning", "note": "note", "seealso": "info", "tip": "tip"}


def _indent_of(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def _block_body(lines: list[str], start: int, indent: int) -> tuple[list[str], int]:
    """The lines after ``start`` that are blank or indented deeper than ``indent``."""
    body = []
    i = start
    while i < len(lines):
        line = lines[i]
        if line.strip() and _indent_of(line) <= indent:
            break
        body.append(line)
        i += 1
    while body and not body[-1].strip():
        body.pop()
        i -= 1
    return body, i


def _dedent(body: list[str]) -> list[str]:
    widths = [_indent_of(line) for line in body if line.strip()]
    cut = min(widths) if widths else 0
    return [line[cut:] if line.strip() else "" for line in body]


def _top_level(tex: str, token: str) -> bool:
    """Whether ``token`` occurs outside every ``\\begin..\\end`` environment."""
    depth = 0
    i = 0
    while i < len(tex):
        if tex.startswith("\\begin{", i):
            depth += 1
            i += 7
        elif tex.startswith("\\end{", i):
            depth -= 1
            i += 5
        elif depth == 0 and tex.startswith(token, i):
            if token == "&" and i and tex[i - 1] == "\\":
                i += 1
                continue
            return True
        else:
            i += 1
    return False


def _math_block(body: list[str], pad: str) -> list[str]:
    """A display block. Arithmatex wraps every block in ``\\[..\\]``, so several
    rows become ``aligned``/``gathered`` (an ``align`` would nest illegally), and
    their ``\\tag``s, which those environments refuse, become visible labels."""
    tex = "\n".join(_dedent(body)).strip()
    rows_present = _top_level(tex, "\\\\")
    aligned = _top_level(tex, "&")
    if rows_present or aligned:
        env = "aligned" if aligned else "gathered"
        tex = re.sub(r"\\tag\{([^}]*)\}", r"\\qquad\\text{(\1)}", tex)
        rows = ["$$", f"\\begin{{{env}}}", *tex.splitlines(), f"\\end{{{env}}}", "$$"]
    else:
        rows = ["$$", *tex.splitlines(), "$$"]
    return ["", *[pad + row if row.strip() else "" for row in rows], ""]


def _simple_table(lines: list[str], start: int) -> tuple[list[str], int] | None:
    """Translate an RST simple table starting at ``lines[start]``; None if it is not one."""
    border = lines[start]
    match = _SIMPLE_BORDER.match(border)
    if not match:
        return None
    pad = match.group(1)
    spans = [(m.start(), m.end()) for m in re.finditer(r"=+", border)]
    borders = [start]
    i = start + 1
    while i < len(lines):
        line = lines[i]
        if _SIMPLE_BORDER.match(line) and line.strip().replace(" ", "") and _same_spans(line, spans):
            borders.append(i)
            # a table ends at a border followed by a blank line or the end
            if len(borders) >= 2 and (i + 1 == len(lines) or not lines[i + 1].strip()):
                break
        elif line.strip() and _indent_of(line) < len(pad):
            return None
        i += 1
    else:
        return None
    end = borders[-1]
    header_rows = []
    if len(borders) == 3:
        header_rows = _table_rows(lines[borders[0] + 1 : borders[1]], spans)
        body_rows = _table_rows(lines[borders[1] + 1 : borders[2]], spans)
    elif len(borders) == 2:
        body_rows = _table_rows(lines[borders[0] + 1 : borders[1]], spans)
    else:
        return None
    ncol = len(spans)
    out = [""]
    head = header_rows[0] if header_rows else [""] * ncol
    out.append(pad + "| " + " | ".join(head) + " |")
    out.append(pad + "|" + "|".join(["---"] * ncol) + "|")
    for row in body_rows:
        out.append(pad + "| " + " | ".join(row) + " |")
    out.append("")
    return out, end + 1


def _same_spans(line: str, spans: list[tuple[int, int]]) -> bool:
    found = [(m.start(), m.end()) for m in re.finditer(r"=+", line)]
    return [s for s, _ in found] == [s for s, _ in spans]


def _table_rows(chunk: list[str], spans: list[tuple[int, int]]) -> list[list[str]]:
    rows: list[list[list[str]]] = []
    for line in chunk:
        if not line.strip():
            continue
        cells = []
        for k, (s, e) in enumerate(spans):
            stop = spans[k + 1][0] if k + 1 < len(spans) else None
            cells.append(line[s:stop].strip() if stop else line[s:].strip())
        if cells[0] or not rows:
            rows.append([[c] if c else [] for c in cells])
        else:
            for k, c in enumerate(cells):
                if c:
                    rows[-1][k].append(c)
    return [[" ".join(parts).replace("|", "\\|") for parts in row] for row in rows]


def _convert_blocks(lines: list[str]) -> list[str]:
    out: list[str] = []
    i = 0
    in_fence = False
    while i < len(lines):
        line = lines[i]
        if _FENCE.match(line):
            in_fence = not in_fence
            out.append(line)
            i += 1
            continue
        if in_fence:
            out.append(line)
            i += 1
            continue

        directive = _DIRECTIVE.match(line)
        if directive:
            pad, name, arg = directive.groups()
            body, nxt = _block_body(lines, i + 1, len(pad))
            # directive options (":caption: x") come first and carry nothing we render
            while body and body[0].strip().startswith(":") and re.match(r"^\s*:[\w-]+:", body[0]):
                body.pop(0)
            if name == "math":
                if arg.strip():
                    body = [pad + "    " + arg.strip(), *body]
                out.extend(_math_block(body, pad))
            elif name in ("code-block", "code", "sourcecode"):
                out.extend(["", f"{pad}```{arg.strip() or 'text'}",
                            *[pad + row if row else "" for row in _dedent(body)],
                            f"{pad}```", ""])
            elif name in _ADMONITIONS:
                inner = _convert_blocks(_dedent(body))
                title = f' "{arg.strip()}"' if arg.strip() else ""
                out.extend(["", f"{pad}!!! {_ADMONITIONS[name]}{title}", ""])
                out.extend(pad + "    " + row if row.strip() else "" for row in inner)
                out.append("")
            else:
                out.append(line)
                out.extend(body)
            i = nxt
            continue

        if _SIMPLE_BORDER.match(line):
            table = _simple_table(lines, i)
            if table:
                rows, i = table
                out.extend(rows)
                continue

        if _GRID_BORDER.match(line):
            pad = " " * _indent_of(line)
            j = i
            while j < len(lines) and lines[j].strip() and lines[j].lstrip()[:1] in "+|":
                j += 1
            out.extend(["", f"{pad}```text", *lines[i:j], f"{pad}```", ""])
            i = j
            continue

        # a section title: text, then an underline at least as long, then a gap
        if (
            i + 1 < len(lines)
            and line.strip()
            and _UNDERLINE.match(lines[i + 1])
            and _indent_of(lines[i + 1]) == _indent_of(line)
            and len(lines[i + 1].strip()) >= len(line.strip())
            and (i == 0 or not lines[i - 1].strip())
        ):
            level = {"=": 2, "-": 2, "~": 3, "^": 4}[lines[i + 1].strip()[0]]
            out.extend([" " * _indent_of(line) + "#" * level + " " + line.strip(), ""])
            i += 2
            continue

        # a literal block: "text::" (or a bare "::"), a gap, then an indented block
        if line.rstrip().endswith("::") and not line.strip().startswith(".."):
            pad = _indent_of(line)
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines) and _indent_of(lines[j]) > pad:
                body, nxt = _block_body(lines, j, pad)
                text = line.rstrip()[:-2].rstrip()
                if text.strip():
                    out.append(text + ":")
                out.extend(["", " " * pad + "```text",
                            *[" " * pad + row if row else "" for row in _dedent(body)],
                            " " * pad + "```", ""])
                i = nxt
                continue

        out.append(line)
        i += 1
    return out


class _Resolver:
    """Resolve a role target in a docstring's own scope, to an object this site renders."""

    def __init__(self, collection, rendered_modules: set[str]):
        self.collection = collection
        self.rendered_modules = rendered_modules

    def _lookup(self, path: str):
        try:
            obj = self.collection[path]
        except (KeyError, ValueError):
            return None
        try:
            return obj.final_target if obj.is_alias else obj
        except Exception:  # an alias pointing outside the package
            return None

    def _rendered(self, obj) -> bool:
        if obj is None:
            return False
        module = obj if obj.is_module else obj.module
        if module.path not in self.rendered_modules:
            return False
        if obj.is_module:
            return True
        rest = obj.path[len(module.path) + 1 :].split(".")
        return not any(part.startswith("_") for part in rest)

    def target(self, scope, name: str) -> str | None:
        name = name.lstrip(".")
        candidates = []
        if name.startswith("plumetools"):
            candidates.append(name)
        head, _, tail = name.partition(".")
        node = scope
        while node is not None:
            try:
                resolved = node.resolve(head)
                candidates.append(resolved + ("." + tail if tail else ""))
            except Exception:
                pass
            if node.is_module:
                candidates.append(f"{node.path}.{name}")
            node = node.parent
        for path in candidates:
            obj = self._lookup(path)
            if self._rendered(obj):
                return obj.path
        return None


def _convert_inline(text: str, scope, resolver: _Resolver) -> str:
    def math(match):
        tex = match.group(1).strip().replace("|", "\\vert ")
        return f"${tex}$"

    def role(match):
        kind, body = match.groups()
        explicit = re.match(r"^(.*?)\s*<([^>]+)>$", body)
        if explicit:
            shown, name = explicit.group(1), explicit.group(2)
        else:
            name = body
            shown = body[1:].rsplit(".", 1)[-1] if body.startswith("~") else body.lstrip(".")
        name = name.lstrip("~")
        path = resolver.target(scope, name) if resolver else None
        return f"[`{shown}`][{path}]" if path else f"`{shown}`"

    # One left-to-right scan, so a role is recognised at its own colon and a
    # code span never starts on a role's closing backtick.
    def token(match):
        if match.group("math") is not None:
            return math(_MATH_ROLE.match(match.group(0)))
        if match.group("role") is not None:
            return role(_ROLE.match(match.group(0)))
        return match.group(0)  # a code span: leave it alone

    return _INLINE.sub(token, text)


def convert(text: str, scope=None, resolver: _Resolver | None = None) -> str:
    """Translate one docstring. ``scope``/``resolver`` enable cross-references."""
    lines = _convert_blocks(text.splitlines())
    joined = "\n".join(lines)
    # inline conversion must not touch fenced blocks
    chunks = re.split(r"(^\s*```.*?^\s*```\s*$)", joined, flags=re.M | re.S)
    for k in range(0, len(chunks), 2):
        chunks[k] = _convert_inline(chunks[k], scope, resolver)
    return re.sub(r"\n{3,}", "\n\n", "".join(chunks)).strip("\n")


class RstToMarkdown(griffe.Extension):
    """Translate every docstring in the loaded package once it is fully loaded."""

    def on_package(self, *, pkg: griffe.Module, loader: griffe.GriffeLoader, **kwargs) -> None:
        if pkg.name != "plumetools":
            return
        rendered = {path for path, obj in _walk_modules(pkg)}
        resolver = _Resolver(loader.modules_collection, rendered)
        for obj in _walk(pkg):
            if obj.docstring is not None:
                obj.docstring.value = convert(obj.docstring.value, obj, resolver)


def _walk_modules(module):
    yield module.path, module
    for member in module.modules.values():
        if not member.is_alias:
            yield from _walk_modules(member)


def _walk(obj):
    yield obj
    for member in list(obj.members.values()):
        if member.is_alias:
            continue
        if member.is_module or member.is_class or member.is_function or member.is_attribute:
            yield from _walk(member)

"""
parse_dck.py
------------
Two parsers for the esclab / TRNSYS workflow:

parse_dck(path)
    Parses a TRNSYS deck (.dck) file and returns a list of dicts, one per
    UNIT, each containing:

        {
            "unit":           int,          # UNIT number
            "type":           int,          # TYPE number
            "name":           str,          # descriptive name on the UNIT line
            "parameters":     list[float],  # values listed under PARAMETERS
            "inputs":         list[tuple],  # [(unit, output_index), ...]
            "initial_inputs": list[float],  # values listed under *** INITIAL INPUT VALUES
        }

parse_component(path)
    Parses an esclab Python component file and returns a dict describing the
    class-level Parameter / Input / Output declarations **in source order**:

        {
            "parameters": list[str],   # names declared as Component.Parameter(...)
            "inputs":     list[str],   # names declared as Component.Input(...)
            "outputs":    list[str],   # names declared as Component.Output(...)
        }

Usage
-----
    from parse_dck import parse_dck, parse_component

    units = parse_dck("predawn_solana.dck")
    comp  = parse_component("../components/flownetwork/TeeOut.py")
"""

from __future__ import annotations

import ast
import re
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _try_numeric(token: str) -> float | str:
    """Return float if *token* looks numeric, otherwise return the raw string."""
    try:
        return float(token)
    except ValueError:
        return token


def _parse_value_line(line: str) -> list[float | str]:
    """
    Strip inline comments (everything after '!') and return all whitespace-
    separated tokens on the line, converted to float where possible.
    """
    content = line.split("!")[0].strip()
    if not content:
        return []
    return [_try_numeric(t) for t in content.split()]


def _parse_input_connection(token: str) -> tuple[int, int]:
    """
    Convert a connection token like '109,1' or '0,0' into (unit, output_idx).
    """
    parts = token.split(",")
    return (int(parts[0]), int(parts[1]))


# ---------------------------------------------------------------------------
# Main parser
# ---------------------------------------------------------------------------

# Matches:  UNIT 117 TYPE 4097   ESOL4097-Weather
_UNIT_RE = re.compile(
    r"^UNIT\s+(\d+)\s+TYPE\s+(\d+)\s*(.*)",
    re.IGNORECASE,
)

# Section header keywords that we care about (count follows on same line)
_SECTION_RE = re.compile(
    r"^(PARAMETERS|INPUTS|LABELS|DERIVATIVES|OUTPUTS)\s+(\d+)",
    re.IGNORECASE,
)


def parse_dck(path: str | Path) -> list[dict[str, Any]]:
    """
    Parse a TRNSYS deck file and return one dict per UNIT block.

    Parameters
    ----------
    path : str or Path
        Path to the .dck file.

    Returns
    -------
    list of dict, each with keys:
        unit, type, name, parameters, inputs, initial_inputs
    """
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()

    units: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None

    # State machine variables
    section: str | None = None      # "PARAMETERS" | "INPUTS" | "LABELS" | None
    items_needed: int = 0           # how many items remain to collect for the section
    reading_initial: bool = False   # True while collecting *** INITIAL INPUT VALUES

    def _flush_section() -> None:
        """Called when we need to switch away from the current section."""
        nonlocal section, items_needed, reading_initial
        section = None
        items_needed = 0
        reading_initial = False

    for raw_line in lines:
        line = raw_line.strip()

        # -----------------------------------------------------------------
        # Skip pure comment lines and blank lines
        # (lines starting with '*' that are NOT the INITIAL INPUT marker)
        # -----------------------------------------------------------------
        if not line:
            continue

        # The INITIAL INPUT VALUES marker
        if line.startswith("*** INITIAL INPUT VALUES"):
            _flush_section()
            if current is not None:
                reading_initial = True
            continue

        # Regular comment / metadata lines
        if line.startswith("*"):
            # A new "*  Model …" comment block signals end of the current unit's
            # initial-input section (if we were reading it) but we don't close
            # `current` yet — we close it when we see the next UNIT line.
            _flush_section()
            continue

        # -----------------------------------------------------------------
        # UNIT line — start a new unit
        # -----------------------------------------------------------------
        m = _UNIT_RE.match(line)
        if m:
            # Save the previous unit if any
            if current is not None:
                units.append(current)
            _flush_section()
            current = {
                "unit": int(m.group(1)),
                "type": int(m.group(2)),
                "name": m.group(3).strip(),
                "parameters": [],
                "inputs": [],
                "initial_inputs": [],
            }
            continue

        # If we haven't seen a UNIT line yet, skip everything
        if current is None:
            continue

        # -----------------------------------------------------------------
        # Section headers inside a UNIT block
        # -----------------------------------------------------------------
        m_sec = _SECTION_RE.match(line)
        if m_sec:
            _flush_section()
            section = m_sec.group(1).upper()
            items_needed = int(m_sec.group(2))
            continue

        # -----------------------------------------------------------------
        # Collect PARAMETERS values
        # -----------------------------------------------------------------
        if section == "PARAMETERS" and items_needed > 0:
            values = _parse_value_line(line)
            for v in values:
                if items_needed <= 0:
                    break
                current["parameters"].append(v)
                items_needed -= 1
            continue

        # -----------------------------------------------------------------
        # Collect INPUTS connections
        # -----------------------------------------------------------------
        if section == "INPUTS" and items_needed > 0:
            # Strip inline comment, then grab all comma-containing tokens
            content = line.split("!")[0].strip()
            tokens = content.split()
            for tok in tokens:
                if items_needed <= 0:
                    break
                if "," in tok:
                    current["inputs"].append(_parse_input_connection(tok))
                    items_needed -= 1
            continue

        # -----------------------------------------------------------------
        # Collect LABELS lines (skip them — not part of our output, but we
        # must consume the right number of lines so the parser stays in sync)
        # -----------------------------------------------------------------
        if section == "LABELS" and items_needed > 0:
            items_needed -= 1
            continue

        # -----------------------------------------------------------------
        # Collect INITIAL INPUT VALUES
        # -----------------------------------------------------------------
        if reading_initial:
            values = _parse_value_line(line)
            if values:
                current["initial_inputs"].extend(values)
            continue

    # Don't forget the last unit
    if current is not None:
        units.append(current)

    return units


# ---------------------------------------------------------------------------
# Component parser
# ---------------------------------------------------------------------------

def parse_component(path: str | Path) -> dict[str, list[str]]:
    """
    Parse an esclab Python component file and return the names of all
    ``Component.Parameter``, ``Component.Input``, and ``Component.Output``
    class-body assignments **in the order they appear in the source**.

    The function walks the AST of the first class definition found in the
    file, collecting assignments whose right-hand side is a call to
    ``Component.Parameter``, ``Component.Input``, or ``Component.Output``
    (case-sensitive, matching the esclab convention).  Only direct class-body
    statements are inspected; nested functions / methods are skipped, so
    local variables that happen to shadow a class attribute are ignored.

    Parameters
    ----------
    path : str or Path
        Absolute or relative path to the ``.py`` component file.

    Returns
    -------
    dict with keys:
        ``"parameters"`` - list[str], names of Component.Parameter(...) attrs  
        ``"inputs"``     - list[str], names of Component.Input(...)     attrs  
        ``"outputs"``    - list[str], names of Component.Output(...)    attrs  

    Raises
    ------
    ValueError
        If no class definition is found in the file.
    """
    source = Path(path).read_text(encoding="utf-8", errors="replace")
    tree = ast.parse(source, filename=str(path))

    # Find the first class node (the component class)
    class_node: ast.ClassDef | None = None
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            class_node = node
            break

    if class_node is None:
        raise ValueError(f"No class definition found in {path!r}")
    
    # Get type number from first line of file
    try:
        # first match in the source using regex
        m = re.search(r"\(Type\s+(\d+)\)", source)
        m = m[1]
    except:
        m = ''

    result: dict[str, list[str]] = {
        "parameters": [],
        "inputs": [],
        "outputs": [],
        "type": m,
    }

    # Map Component.<Attr> call name → result key
    _CALL_MAP = {
        "Parameter": "parameters",
        "Input":     "inputs",
        "Output":    "outputs",
    }

    def _classify_call(node: ast.expr) -> str | None:
        """
        Return 'parameters' | 'inputs' | 'outputs' if *node* is a call to
        ``Component.Parameter / Input / Output``, else None.

        Handles both:
          - ``Component.Parameter(...)``          → ast.Attribute on ast.Name
          - bare ``Parameter(...)`` / ``Input(...)`` / ``Output(...)``
            (uncommon but valid after a star-import or local alias)
        """
        if not isinstance(node, ast.Call):
            return None
        func = node.func
        # Component.Parameter(...)
        if (
            isinstance(func, ast.Attribute)
            and func.attr in _CALL_MAP
            and isinstance(func.value, ast.Name)
            and func.value.id == "Component"
        ):
            return _CALL_MAP[func.attr]
        # Bare Parameter(...) / Input(...) / Output(...)
        if isinstance(func, ast.Name) and func.id in _CALL_MAP:
            return _CALL_MAP[func.id]
        return None

    # Only look at *direct* class-body statements (not inside methods)
    for stmt in class_node.body:
        # Simple assignment:  name = Component.Parameter(...)
        if isinstance(stmt, ast.Assign):
            category = _classify_call(stmt.value)
            if category is None:
                continue
            for target in stmt.targets:
                if isinstance(target, ast.Name):
                    result[category].append(target.id)

        # Annotated assignment:  name: type = Component.Parameter(...)
        elif isinstance(stmt, ast.AnnAssign):
            if stmt.value is None:
                continue
            category = _classify_call(stmt.value)
            if category is None:
                continue
            if isinstance(stmt.target, ast.Name):
                result[category].append(stmt.target.id)

    return result


# ---------------------------------------------------------------------------
# CLI / quick test
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    dck_path = Path(__file__).parent / "predawn_solana.dck"
    units = parse_dck(dck_path)

    # print(f"Parsed {len(units)} units.\n")
    # for u in units:
    #     print(
    #         f"  UNIT {u['unit']:>4d}  TYPE {u['type']:>5d}  "
    #         f"name={u['name']!r:40s}  "
    #         f"params={len(u['parameters'])}  "
    #         f"inputs={len(u['inputs'])}  "
    #         f"init_inputs={len(u['initial_inputs'])}"
    #     )

    tee = parse_component(Path(__file__).parent.parent / "components/flownetwork/TeeOut.py")

    pass

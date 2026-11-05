"""Stage dispatch for pripipe.

Each stage is a function `stage_<name>(cfg, rep) -> None` that writes only
inside `cfg.stage_dir(<name>)` plus a `<stage>.done` marker. The marker is
what makes `--from/--to` and re-runs safe.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

from ..config import RunConfig

# Stage order is defined by config.STAGES; every stage below is a placeholder
# until its implementation lands. Keeping the signatures stable lets the CLI
# and manifest code be written first.
IMPLEMENTED: set[str] = set()


def _not_implemented(stage: str) -> None:
    sys.stderr.write(
        f"[pripipe] stage '{stage}' is not implemented yet; "
        f"run `pripipe validate` to check inputs meanwhile\n"
    )
    raise SystemExit(3)


def stage_iso(cfg: RunConfig, rep) -> None:
    from .iso import stage_iso as _impl

    _impl(cfg, rep)


def stage_rna(cfg: RunConfig, rep) -> None:
    from .rna import stage_rna as _impl

    _impl(cfg, rep)


def stage_assemble(cfg: RunConfig, rep) -> None:
    from .assemble import stage_assemble as _impl

    _impl(cfg, rep)


def stage_catalog(cfg: RunConfig, rep) -> None:
    from .catalog import stage_catalog as _impl

    _impl(cfg, rep)


def stage_anchor(cfg: RunConfig, rep) -> None:
    from .anchor import stage_anchor as _impl

    _impl(cfg, rep)


def stage_quant(cfg: RunConfig, rep) -> None:
    from .quant import stage_quant as _impl

    _impl(cfg, rep)


def stage_report(cfg: RunConfig, rep) -> None:
    from .report import stage_report as _impl

    _impl(cfg, rep)


_DISPATCH = {
    "iso": stage_iso,
    "rna": stage_rna,
    "assemble": stage_assemble,
    "catalog": stage_catalog,
    "anchor": stage_anchor,
    "quant": stage_quant,
    "report": stage_report,
}


def run_stages(cfg: RunConfig, rep, stages: list[str]) -> int:
    started = time.time()
    for name in stages:
        marker = cfg.stage_dir(name) / f"{name}.done"
        if marker.exists():
            sys.stderr.write(f"[pripipe] skip {name} (done)\n")
            continue
        sys.stderr.write(f"[pripipe] === {name} ===\n")
        _DISPATCH[name](cfg, rep)
        marker.write_text(time.strftime("%Y-%m-%dT%H:%M:%S"), encoding="utf-8")
    sys.stderr.write(f"[pripipe] finished in {time.time() - started:.1f}s\n")
    return 0

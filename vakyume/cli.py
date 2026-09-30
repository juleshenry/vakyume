"""``vakyume`` command line.

    vakyume list PROJECT                      # equations and their kinds
    vakyume solve PROJECT 2-1 rho=1.2 D=0.1 v=3 mu=1.8e-5
    vakyume report PROJECT [--only 2-34] [--jobs 4]
    vakyume convert-notes PROJECT             # notes/*.py -> equations/*.toml

The legacy scrape/run/reconstruct/make-cpp commands live in ``vakyume.py``.
"""

from __future__ import annotations

import multiprocessing as mp
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Annotated

import typer

from .equations import Library, dump_toml, parse_notes
from .verify import EquationResult, summarize, verify_equation, write_reports

app = typer.Typer(no_args_is_help=True, add_completion=False)

ProjectArg = Annotated[Path, typer.Argument(exists=True, file_okay=False, help="Project directory")]


@app.command("list")
def list_(project: ProjectArg, kind: str | None = None) -> None:
    """List equations with their parsed kind."""
    lib = Library.load(project)
    for eq in lib.equations:
        if kind and eq.kind != kind:
            continue
        extra = f"  ({eq.relation.error})" if eq.relation.error else ""
        typer.echo(f"{eq.id:10} {eq.kind:12} {eq.text}{extra}")
    typer.echo(dict(Counter(e.kind for e in lib.equations)))


@app.command()
def solve(
    project: ProjectArg,
    eq_id: str,
    assignments: list[str],
    select: str = typer.Option("all", help="unique | all | min | max"),
) -> None:
    """Solve an equation for the one variable not given (NAME=VALUE ...)."""
    lib = Library.load(project)
    eq = lib[eq_id]
    knowns = {}
    for a in assignments:
        name, _, value = a.partition("=")
        knowns[name.strip()] = float(value)
    typer.echo(f"{eq.id}: {eq.text}")
    typer.echo(repr(eq.solve(knowns, select=select)))


@app.command("convert-notes")
def convert_notes(project: ProjectArg, force: bool = False) -> None:
    """One-time conversion of legacy notes/*.py into equations/*.toml."""
    out_dir = project / "equations"
    out_dir.mkdir(exist_ok=True)
    for path in sorted((project / "notes").glob("*.py")):
        if path.name.startswith("_"):
            continue
        target = out_dir / f"{path.stem}.toml"
        if target.exists() and not force:
            typer.echo(f"skip {target} (exists; --force to overwrite)")
            continue
        target.write_text(dump_toml(parse_notes(path)))
        typer.echo(f"wrote {target}")


def _verify_one(args: tuple[str, str, int, int]) -> EquationResult:
    project, eq_id, n, seed = args
    return verify_equation(Library.load(project)[eq_id], n=n, seed=seed)


@app.command()
def report(
    project: ProjectArg,
    only: Annotated[list[str] | None, typer.Option(help="Equation ids to verify (repeatable)")] = None,
    samples: int = typer.Option(12, help="Consistent sample points per equation"),
    seed: int = 0,
    jobs: int = typer.Option(1, help="Worker processes (verification is CPU heavy)"),
    write: bool = typer.Option(True, help="Write docs/STATUS.md and docs/status.json"),
) -> None:
    """Verify every solver against its equation and write an honest status report."""
    lib = Library.load(project)
    eqs = [e for e in lib.equations if not only or e.id in only]
    jobs_args = [(str(project), e.id, samples, seed) for e in eqs]
    t0 = time.time()
    results: list[EquationResult] = []

    def show(r: EquationResult) -> None:
        bad = [f"{t.target}:{t.status}" for t in r.targets if t.status != "pass"]
        typer.echo(f"{r.id:10} {r.status:11} {' '.join(bad)} {r.note[:60] if r.status != 'pass' else ''}")

    if jobs > 1:
        ctx = mp.get_context("fork" if "fork" in mp.get_all_start_methods() else "spawn")
        with ProcessPoolExecutor(max_workers=jobs, mp_context=ctx) as pool:
            for r in pool.map(_verify_one, jobs_args):
                results.append(r)
                show(r)
    else:
        for a in jobs_args:
            r = _verify_one(a)
            results.append(r)
            show(r)

    typer.echo(f"\n{summarize(results)}  ({time.time() - t0:.1f}s)")
    if write and not only:
        md, js = write_reports(results, project, project.name)
        typer.echo(f"wrote {md} and {js}")


if __name__ == "__main__":
    app()

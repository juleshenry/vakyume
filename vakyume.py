import argparse
import os
import sys
from vakyume.pipeline import run_pipeline
from vakyume.cpp_gen import main as run_cpp_gen
from vakyume.llm import ask_llm
from vakyume.reconstruct import reconstruct_cli


def extract_notes_from_pdf(pdf_path, output_dir):
    print(f"Attempting to extract notes from PDF: {pdf_path}")
    # In a real scenario, we'd use a PDF library.
    # Here we'll check for 'pdftotext' utility or just simulate if missing.
    try:
        import subprocess

        result = subprocess.run(
            ["pdftotext", pdf_path, "-"], capture_output=True, text=True
        )
        if result.returncode == 0:
            content = result.stdout
        else:
            print("pdftotext failed. Falling back to simple placeholder extraction.")
            content = (
                f"Simulated content from {pdf_path}. Equation 1-1: p * V = n * R * T"
            )
    except Exception:
        content = f"Simulated content from {pdf_path}. Equation 1-1: p * V = n * R * T"

    system_prompt = (
        "You are a formula extraction tool. You convert textbook text into a specific Python-like format.\n\n"
        "Input: Chapter 1-7 ideal gas law. p is pressure, V is volume. p * V = nRT\n"
        "Output:\n"
        "# 1-7 ideal gas law\n"
        '"""\n'
        "p := pressure\n"
        "V := volume\n"
        '"""\n'
        "p * V = n * R * T\n\n"
        "Input: {content}\n"
        "Output:"
    )

    formatted_formulas = ask_llm(system_prompt.format(content=content), "")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_path = os.path.join(output_dir, "extracted_notes.py")
    with open(output_path, "w") as f:
        f.write(formatted_formulas)
    print(f"Saved extracted notes to {output_path}")


def _ensure_project(args):
    """Create project dirs and handle missing notes (shared by cmd_run/cmd_build)."""
    project_dir = args.project
    notes_dir = os.path.join(project_dir, "notes")

    for sub in ["notes", "shards", "reports", "repair_prompts"]:
        d = os.path.join(project_dir, sub)
        if not os.path.exists(d):
            os.makedirs(d)

    if not os.listdir(notes_dir):
        if getattr(args, "pdf", None):
            extract_notes_from_pdf(args.pdf, notes_dir)
        else:
            print(f"Warning: No files found in {notes_dir} and no --pdf provided.")
            print("OCR stage skipped.")

    return project_dir


def cmd_run(args):
    """Run the full pipeline (shard, verify, repair, certify)."""
    project_dir = _ensure_project(args)

    # Step 3: Run Python Pipeline
    print(f"Starting Vakyume pipeline for project '{project_dir}'...")
    analysis = run_pipeline(
        project_dir=project_dir,
        max_rounds=args.max_rounds,
        overwrite=args.overwrite,
        verbose=args.verbose,
        repair_only=args.repair_only,
    )

    print("\nPipeline Stage: Verification & Certification Complete.")
    print(f"Solved: {len(analysis['solved'])}")
    print(f"Inconsistent: {len(analysis['inconsistent'])}")
    print(f"Failed: {len(analysis['failed'])}")

    # Step 4: Run C++ Generation
    if args.cpp:
        print("\nPipeline Stage: C++ Generation...")
        run_cpp_gen(project_dir=project_dir)


def cmd_reconstruct(args):
    """Reconstruct shards into a single importable Python library."""
    reconstruct_cli(
        project_dir=args.project,
        output=args.output,
        stdout=args.stdout,
    )


def cmd_make_cpp(args):
    """Convert the reconstructed Python library (or certified file) to C++."""
    run_cpp_gen(
        project_dir=args.project,
        input_file=args.input,
    )


def cmd_build(args):
    """Full build: scrape, verify, reconstruct, and generate C++."""
    project_dir = _ensure_project(args)

    # Step 3: Run Python Pipeline (shard, verify, repair, certify)
    print(f"\n=== [build] Scraping & Verification ===")
    print(f"Starting Vakyume pipeline for project '{project_dir}'...")
    analysis = run_pipeline(
        project_dir=project_dir,
        max_rounds=args.max_rounds,
        overwrite=args.overwrite,
        verbose=args.verbose,
        repair_only=False,
    )

    print(f"\nSolved: {len(analysis['solved'])}")
    print(f"Inconsistent: {len(analysis['inconsistent'])}")
    print(f"Failed: {len(analysis['failed'])}")

    # Step 4: Reconstruct shards into single Python library
    print(f"\n=== [build] Reconstruction ===")
    reconstruct_cli(project_dir=project_dir)

    # Step 5: Generate C++
    print(f"\n=== [build] C++ Generation ===")
    run_cpp_gen(project_dir=project_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Vakyume: Vacuum Theory Equation Pipeline"
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # ── run (default pipeline) ───────────────────────────────────────────
    run_parser = subparsers.add_parser(
        "run", help="Run the full pipeline (shard, verify, repair, certify)"
    )
    run_parser.add_argument(
        "project", help="Project directory (e.g., projects/VacuumTheory)"
    )
    run_parser.add_argument("--pdf", help="Path to textbook PDF for OCR stage")
    run_parser.add_argument(
        "--cpp", action="store_true", help="Generate C++ library after certification"
    )
    run_parser.add_argument(
        "--max-rounds", type=int, default=10, help="Maximum repair rounds"
    )
    run_parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing shards"
    )
    run_parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show detailed verification output for all families",
    )
    run_parser.add_argument(
        "--repair-only",
        action="store_true",
        help="Skip shard generation and re-verification of passing families; only repair broken ones",
    )
    run_parser.set_defaults(func=cmd_run)

    # ── reconstruct ──────────────────────────────────────────────────────
    recon_parser = subparsers.add_parser(
        "reconstruct",
        help="Rebuild shards into a single importable Python library",
    )
    recon_parser.add_argument(
        "project", help="Project directory (e.g., projects/VacuumTheory)"
    )
    recon_parser.add_argument(
        "-o",
        "--output",
        help="Output file path (default: <project>/vakyume_lib.py)",
    )
    recon_parser.add_argument(
        "--stdout",
        action="store_true",
        help="Print generated source to stdout instead of writing a file",
    )
    recon_parser.set_defaults(func=cmd_reconstruct)

    # ── make-cpp ─────────────────────────────────────────────────────────
    cpp_parser = subparsers.add_parser(
        "make-cpp",
        help="Convert the Python library to a compiled C++ library",
    )
    cpp_parser.add_argument(
        "project", help="Project directory (e.g., projects/VacuumTheory)"
    )
    cpp_parser.add_argument(
        "-i",
        "--input",
        help="Python source file to convert (default: <project>/vakyume_lib.py or vakyume_certified.py)",
    )
    cpp_parser.set_defaults(func=cmd_make_cpp)

    # ── build (all-in-one) ───────────────────────────────────────────────
    build_parser = subparsers.add_parser(
        "build",
        help="Full build: scrape, verify, reconstruct, and generate C++",
    )
    build_parser.add_argument(
        "project", help="Project directory (e.g., projects/VacuumTheory)"
    )
    build_parser.add_argument("--pdf", help="Path to textbook PDF for OCR stage")
    build_parser.add_argument(
        "--max-rounds", type=int, default=10, help="Maximum repair rounds"
    )
    build_parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing shards"
    )
    build_parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show detailed verification output for all families",
    )
    build_parser.set_defaults(func=cmd_build)

    # ── Parse & dispatch ─────────────────────────────────────────────────
    args = parser.parse_args()

    if not args.command:
        # Backwards compatibility: if no subcommand given, try legacy style
        # (positional project arg with flags).  Re-parse with the run parser.
        if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
            sys.argv.insert(1, "run")
            args = parser.parse_args()
        else:
            parser.print_help()
            sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()

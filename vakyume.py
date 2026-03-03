import argparse
import os
import sys
import shutil
from vakyume.master import run_pipeline
from vakyume.cpp_gen import main as run_cpp_gen
from vakyume.llm import ask_llm


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


def main():
    parser = argparse.ArgumentParser(
        description="Vakyume: Vacuum Theory Equation Pipeline"
    )
    parser.add_argument(
        "project", help="Project directory (e.g., projects/VacuumTheory)"
    )
    parser.add_argument("--pdf", help="Path to textbook PDF for OCR stage")
    parser.add_argument(
        "--cpp", action="store_true", help="Generate C++ library after certification"
    )
    parser.add_argument(
        "--max-rounds", type=int, default=10, help="Maximum repair rounds"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing shards"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show detailed verification output for all families",
    )
    parser.add_argument(
        "--repair-only",
        action="store_true",
        help="Skip shard generation and re-verification of passing families; only repair broken ones",
    )

    args = parser.parse_args()

    project_dir = args.project
    notes_dir = os.path.join(project_dir, "notes")

    # Step 1: Ensure project structure exists
    for sub in ["notes", "shards", "reports", "repair_prompts"]:
        d = os.path.join(project_dir, sub)
        if not os.path.exists(d):
            os.makedirs(d)

    # Step 2: Handle missing notes
    if not os.listdir(notes_dir):
        if args.pdf:
            extract_notes_from_pdf(args.pdf, notes_dir)
        else:
            print(f"Warning: No files found in {notes_dir} and no --pdf provided.")
            print("OCR stage skipped.")

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


if __name__ == "__main__":
    main()

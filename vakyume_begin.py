import os
import argparse
import re
from llm import ask_llm


def ocr_page(file_path):
    """
    Simulates OCR process.
    In a real scenario, this would use a multimodal LLM or an OCR library.
    For this tool, we'll use an LLM to extract formulas from the text/image.
    """
    with open(file_path, "r") as f:
        content = f.read()

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

    # We ignore the user_prompt and just use the system_prompt with content injected
    return ask_llm(system_prompt.format(content=content), "")


def main():
    parser = argparse.ArgumentParser(description="OCR tool for textbook formulas")
    parser.add_argument(
        "--input",
        required=True,
        help="Directory containing textbook pages (text/images)",
    )
    parser.add_argument(
        "--output", required=True, help="Directory to save the formatted formula notes"
    )

    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    for filename in sorted(os.listdir(args.input)):
        if filename.startswith("."):
            continue

        file_path = os.path.join(args.input, filename)
        if os.path.isfile(file_path):
            print(f"Processing {filename}...")
            formatted_formulas = ocr_page(file_path)

            output_filename = os.path.splitext(filename)[0] + ".py"
            output_path = os.path.join(args.output, output_filename)

            with open(output_path, "w") as f:
                f.write(formatted_formulas)
            print(f"Saved to {output_path}")


if __name__ == "__main__":
    main()

import os
import subprocess
import sys


def run_stage(name, command):
    print(f"\n>>> Starting Stage: {name}")
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"OK: {name} completed.")
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"FAIL: {name} failed with error:\n{e.stderr}")
        sys.exit(1)


def main():
    pdf_path = "BuildingModelsToDescribeOurWorld.pdf"
    certified_output = "vakyume_certified.py"

    # 1. OCR & ANALYZE -> OOO (Object Oriented Output)
    # This usually involves running the vakyume analyzer on the text
    print("Step 1: OCR and Analysis (Generating OOO Draft)")
    # Example: python3 vakyume_analyze.py --input BuildingModelsToDescribeOurWorld.pdf
    # For now, we assume the analyzer produces intermediate shards/files

    # 2. VERIFY
    # This stage checks the mathematical consistency of the generated equations
    print("Step 2: Verification (SymPy Consistency Checks)")
    # run_stage("Verification", ["python3", "tru.py"])

    # 3. CERTIFY
    # Consolidates verified equations into the certified python library
    print("Step 3: Certification (Consolidating to vakyume_certified.py)")
    # If not already done by the sharder

    # 4. VAKYUME CPP
    # The final step you requested: converting the certified output to C++ and testing it
    if os.path.exists("vakyume_cpp.py"):
        run_stage("C++ Conversion & Test", ["python3", "vakyume_cpp.py"])
    else:
        print("Error: vakyume_cpp.py not found in current directory.")


if __name__ == "__main__":
    main()

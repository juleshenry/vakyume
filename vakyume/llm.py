import ollama

n = "\n"


def ask_llm(system_prompt: str, user_prompt: str, model: str = "phi3:latest"):
    """Base call to Ollama."""
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": user_prompt},
    ]
    response = ollama.chat(model=model, messages=messages)
    return response["message"]["content"]


def escribir_codigo(
    eqn: str, lang: str = "Python", single_variable=None, header=None, **kwargs
):
    """Refined prompt for Phi-3 to solve vacuum equations from scratch."""

    target = single_variable or "result"
    knowns = ""
    if header and "(" in header:
        knowns = header.split("(")[1].split(")")[0]

    system_prompt = (
        f"You are a mathematical programming assistant. Your task is to solve equations for a specific variable in {lang}.\n"
        f"CONSTRAINTS:\n"
        f"- Output ONLY executable {lang} code.\n"
        f"- No commentary, no markdown formatting, no explanations.\n"
        f"- Use 'math' or 'numpy' for functions.\n"
        f"- IMPORTANT: Never use leading zeros in decimal literals (e.g., use 0.5, not 00.5).\n"
        f"- Ensure the result is returned as a list."
    )

    user_prompt = (
        f"Equation: {eqn}\n"
        f"Solve for: {target}\n"
        f"Known variables: {knowns}\n"
        f"Function Header: {header if header else f'def solve_for_{target}({knowns}):'}\n\n"
        f"Provide the {lang} function implementation now (code only):"
    )

    return ask_llm(system_prompt, user_prompt)


def repair_codigo(
    shard_file, shard_code, error=None, inconsistent_variants=None, scores=None
):
    """Refined prompt for Phi-3 to repair existing code."""

    system_prompt = (
        "You are a mathematical code repair assistant. Your goal is to fix inconsistencies or errors in Python equation solver classes.\n"
        "CONSTRAINTS:\n"
        "- Return ONLY the corrected Python methods.\n"
        "- Do NOT use leading zeros in decimal literals.\n"
        "- Ensure all methods in the class are mathematically equivalent.\n"
        "- Return executable code only, no markdown markers."
    )

    if error:
        user_prompt = (
            f"The Python shard '{shard_file}' failed with error: {error}\n\n"
            "CODE TO FIX:\n"
            f"{shard_code}\n\n"
            "Please fix the syntax and logical errors. Return the corrected methods."
        )
    else:
        user_prompt = (
            f"The Python shard '{shard_file}' has inconsistent variants.\n"
            f"Likely incorrect: {inconsistent_variants}\n"
            f"Scores: {scores}\n\n"
            "CODE TO REPAIR:\n"
            f"{shard_code}\n\n"
            "Fix the inconsistent variants so they all yield the same results. Return the corrected methods."
        )

    return ask_llm(system_prompt, user_prompt)


def extract_code(text):
    """Extracts code snippet from markdown blocks or returns text if no blocks found."""
    if "```python" in text:
        parts = text.split("```python")
        if len(parts) > 1:
            code = parts[1].split("```")[0]
            return code.strip()
    elif "```" in text:
        parts = text.split("```")
        if len(parts) > 1:
            code = parts[1].split("```")[0]
            return code.strip()
    return text.strip()

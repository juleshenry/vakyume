import ollama
import re

n = "\n"


def ask_llm(
    system_prompt: str, user_prompt: str, model: str = "phi3:latest", stream=False
):
    """Base call to Ollama."""
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": user_prompt},
    ]
    response = ollama.chat(
        model=model,
        messages=messages,
        stream=stream,
        options={"num_ctx": 8192, "temperature": 0},
    )
    if stream:
        return response
    # The response object can be a dict or a ChatResponse object depending on the version
    if hasattr(response, "message"):
        return response.message.content
    return response["message"]["content"]


def escribir_codigo(
    eqn: str,
    lang: str = "Python",
    single_variable=None,
    header=None,
    stream=False,
    **kwargs,
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
        f"- Do NOT add any chatty comments, explanations, or docstrings.\n"
        f"- Use 'math' or 'numpy' for functions.\n"
        f"- IMPORTANT: Never use leading zeros in decimal literals (e.g., use 0.5, not 00.5).\n"
        f"- DO NOT remove or rename arguments in the provided function header.\n"
        f"- Ensure the result is returned as a list."
    )

    user_prompt = (
        f"Equation: {eqn}\n"
        f"Solve for: {target}\n"
        f"Known variables: {knowns}\n"
        f"Function Header: {header if header else f'def solve_for_{target}({knowns}):'}\n\n"
        f"Provide the {lang} function implementation now (code only):"
    )

    return ask_llm(system_prompt, user_prompt, stream=stream)


def repair_codigo(
    shard_file,
    shard_code,
    error=None,
    broken_variants=None,
    trusted_variants=None,
    scores=None,
    pyeqn=None,
    stream=False,
    **kwargs,
):
    """Refined prompt for Phi-3 to repair existing code with majority context."""

    system_prompt = (
        "You are a mathematical code repair assistant. Your goal is to fix inconsistencies or errors in Python equation solver classes.\n"
        "CONSTRAINTS:\n"
        "- Return ONLY the corrected Python methods.\n"
        "- Include necessary decorators like '@staticmethod' for each method.\n"
        "- Do NOT use leading zeros in decimal literals.\n"
        "- Do NOT add any chatty comments, explanations, or docstrings.\n"
        "- ALWAYS preserve the '# [.pyeqn] <original_equation>' comment at the start of each method body.\n"
        "- DO NOT remove or rename arguments in the method signature. Keep the signature exactly as provided.\n"
        "- Ensure all methods in the class are mathematically equivalent.\n"
        "- Use the 'trusted' variants as the ground truth if provided.\n"
        "- Return executable code only, no markdown markers."
    )

    context = ""
    if pyeqn:
        context = f"Original Equation: {pyeqn}\n"

    if error:
        user_prompt = (
            f"The Python shard '{shard_file}' failed with error: {error}\n\n"
            f"{context}"
            "CODE TO FIX:\n"
            f"{shard_code}\n\n"
            "Please fix the syntax and logical errors. Return the corrected methods."
        )
    elif trusted_variants:
        user_prompt = (
            f"The Python shard '{shard_file}' has inconsistent variants.\n"
            f"{context}"
            f"TRUSTED (Correct): {trusted_variants}\n"
            f"BROKEN (Incorrect): {broken_variants}\n"
            f"Scores: {scores}\n\n"
            "CODE TO REPAIR:\n"
            f"{shard_code}\n\n"
            f"The variants {trusted_variants} are confirmed correct. Repair the variants {broken_variants} "
            f"to be mathematically equivalent to the trusted ones. DO NOT modify the trusted variants."
        )
    else:
        user_prompt = (
            f"The Python shard '{shard_file}' has inconsistent variants.\n"
            f"{context}"
            f"Inconsistent: {broken_variants}\n"
            f"Scores: {scores}\n\n"
            "CODE TO REPAIR:\n"
            f"{shard_code}\n\n"
            "Fix the inconsistent variants so they all yield the same results. Return the corrected methods."
        )

    return ask_llm(system_prompt, user_prompt, stream=stream)


def extract_code(text):
    """Extracts code snippet from markdown blocks or returns text if no blocks found."""
    # Try to find all code blocks and join them
    pattern = re.compile(r"```(?:python)?\n(.*?)\n```", re.DOTALL)
    matches = pattern.findall(text)
    if matches:
        return "\n\n".join(matches).strip()

    # Fallback for single backticks or no blocks
    if "```" in text:
        parts = text.split("```")
        # If there's at least one pair of ```
        if len(parts) >= 3:
            # Join every odd part (between ```)
            code_parts = [parts[i].strip() for i in range(1, len(parts), 2)]
            return "\n\n".join(code_parts).strip()

    return text.strip()

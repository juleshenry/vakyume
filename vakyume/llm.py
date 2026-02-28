import ollama
import re

n = "\n"


def ask_llm(
    system_prompt: str, user_prompt: str, model: str = "phi3:latest", stream=False
):
    """Base call to Ollama."""
    print(
        f"[INPUT] ask_llm: model={model}, system_prompt_len={len(system_prompt)}, user_prompt_len={len(user_prompt)}"
    )
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
    content = ""
    if hasattr(response, "message"):
        content = response.message.content
    else:
        content = response["message"]["content"]

    print(f"[OUTPUT] ask_llm: content_len={len(content)}")
    return content


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
    mismatches=None,
    pyeqn=None,
    stream=False,
    is_subshard=True,
    **kwargs,
):
    """Refined prompt for Phi-3 to repair existing code with majority context."""
    print(
        f"[INPUT] repair_codigo: shard_file={shard_file}, broken={broken_variants}, trusted={trusted_variants}, is_subshard={is_subshard}"
    )

    if is_subshard:
        system_prompt = (
            "You are a mathematical code repair assistant. Your goal is to fix inconsistencies or errors in Python equation solver functions.\n"
            "CONSTRAINTS:\n"
            "- Return ONLY the corrected Python functions.\n"
            "- Do NOT use classes or decorators like '@staticmethod'.\n"
            "- Do NOT use leading zeros in decimal literals.\n"
            "- Do NOT add any chatty comments, explanations, or docstrings.\n"
            "- ALWAYS preserve the '# [.pyeqn] <original_equation>' comment at the start of each function body.\n"
            "- DO NOT remove or rename arguments in the function signature. Keep the signature exactly as provided.\n"
            "- If SymPy cannot solve it analytically, provide a robust numerical/approximate solution (e.g., using scipy.optimize.newton).\n"
            "- Return executable code only, no markdown markers."
        )
    else:
        system_prompt = (
            "You are a mathematical code repair assistant. Your goal is to fix inconsistencies or errors in Python equation solver classes.\n"
            "CONSTRAINTS:\n"
            "- Return ONLY the corrected Python class and methods.\n"
            "- Include necessary decorators like '@staticmethod' for each method.\n"
            "- Do NOT use leading zeros in decimal literals.\n"
            "- Do NOT add any chatty comments, explanations, or docstrings.\n"
            "- ALWAYS preserve the '# [.pyeqn] <original_equation>' comment at the start of each method body.\n"
            "- DO NOT remove or rename arguments in the method signature. Keep the signature exactly as provided.\n"
            "- Ensure all methods in the class are mathematically equivalent and satisfy the original equation.\n"
            "- Use the 'trusted' variants as the ground truth if provided.\n"
            "- If SymPy cannot solve it analytically, provide a robust numerical/approximate solution (e.g., using scipy.optimize.newton).\n"
            "- Return executable code only, no markdown markers."
        )

    context = ""
    if pyeqn:
        context = f"Original Equation: {pyeqn}\n"

    mismatch_context = ""
    if mismatches:
        mismatch_context = "HARMONY FAILURES (Specific examples where code is wrong):\n"
        for var, trials in mismatches.items():
            if trials:
                # Just show the first few trials to keep prompt size manageable
                for i, trial in enumerate(trials[:2]):
                    mismatch_context += f"Variant {var}, Example {i + 1}:\n"
                    if "inputs" in trial:
                        mismatch_context += f"  Inputs: {trial['inputs']}\n"
                        mismatch_context += f"  Got Output: {trial['output']}\n"
                        for m in trial.get("mismatches", []):
                            if "target" in m:
                                mismatch_context += f"  Mismatch with {m['target']}:"
                                if "expected" in m:
                                    mismatch_context += f" expected {m['expected']},"
                                if "got" in m:
                                    mismatch_context += f" got {m['got']}"
                                mismatch_context += "\n"
                            elif "error" in m:
                                mismatch_context += (
                                    f"  Error in {m['target']}: {m['error']}\n"
                                )
                    elif "error" in trial:
                        mismatch_context += f"  Error: {trial['error']}\n"
        mismatch_context += "\n"

    if error:
        user_prompt = (
            f"The Python shard '{shard_file}' failed with error: {error}\n\n"
            f"{context}"
            f"{mismatch_context}"
            "CODE TO FIX:\n"
            f"{shard_code}\n\n"
            "Please fix the syntax and logical errors. Provide an analytic or robust numerical solution. Return the corrected methods."
        )
    elif trusted_variants:
        user_prompt = (
            f"The Python shard '{shard_file}' has inconsistent variants.\n"
            f"{context}"
            f"{mismatch_context}"
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
            f"{mismatch_context}"
            f"Inconsistent: {broken_variants}\n"
            f"Scores: {scores}\n\n"
            "CODE TO REPAIR:\n"
            f"{shard_code}\n\n"
            "Fix the inconsistent variants so they all yield the same results (analytic or numerical). Return the corrected methods."
        )

    return ask_llm(system_prompt, user_prompt, stream=stream)


def extract_code(text, target_name=None):
    """Extracts ONLY the first valid Python code block or a specific function/class definition."""
    # 1. Try to find a markdown code block
    pattern = re.compile(r"```(?:python)?\n(.*?)\n```", re.DOTALL)
    matches = pattern.findall(text)
    
    # If we have a target_name, look through matches or text for it
    if target_name:
        search_text = "\n\n".join(matches) if matches else text
        lines = search_text.splitlines()
        code_lines = []
        found_target = False
        
        # Look for 'def target_name' or 'class target_name'
        target_pattern = re.compile(rf"^\s*(def|class)\s+{re.escape(target_name)}\b")
        
        for line in lines:
            if target_pattern.match(line):
                found_target = True
            elif found_target and (line.strip().startswith("def ") or line.strip().startswith("class ")):
                # Found the start of another definition, stop here
                break
            
            if found_target:
                code_lines.append(line)
        
        if code_lines:
            return "\n".join(code_lines).strip()

    # Fallback to first match if no target_name or target not found in matches
    if matches:
        return matches[0].strip()

    # 2. Try to find the first 'def' and take everything until next 'def' or end
    lines = text.splitlines()
    code_lines = []
    found_def = False
    for line in lines:
        if line.strip().startswith("def ") or line.strip().startswith("class "):
            if found_def:
                break # Stop at second definition
            found_def = True
        if found_def:
            code_lines.append(line)
    
    if code_lines:
        return "\n".join(code_lines).strip()

    return text.strip()

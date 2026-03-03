TAB = " " * 4
TYPE = ": float"
STD = 1
OUTFILE = "vakyume_2025.py"
MAX_COMP_TIME_SECONDS = 1
FUNKTORZ = "*()/-+"
n = "\n"
COOLDOWN_SECONDS = 2


class IntractableSolution(Exception):
    """Exception raised when only numerical methods are allowed."""

    pass


class OllamaOffline(Exception):
    """Exception raised when the Ollama is offline."""

    pass


class UnsolvedException(Exception):
    """Exception raised when equation solving is unavailable."""

    pass

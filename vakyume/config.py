class IntractableSolution(Exception):
    """Exception raised when only numerical methods are allowed."""

    pass


class OllamaOffline(Exception):
    """Exception raised when the Ollama is offline."""

    pass


class UnsolvedException(Exception):
    """Exception raised when equation solving is unavailable."""

    pass


# Code generation constants
TAB = " " * 4
MAX_COMP_TIME_SECONDS = 1
COOLDOWN_SECONDS = 2
FUNKTORZ = "*()/-+"

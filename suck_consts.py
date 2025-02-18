TAB = " " * 4
TYPE = ": float"
STD = 1
OUTFILE = "vakyume_2025.py"
MAX_COMP_TIME_SECONDS = 1
FUNKTORZ = "*()/-+"
n = "\n"
class IntractableSolution(Exception):
    """Exception raised when only numerical methods are allowed."""
    pass
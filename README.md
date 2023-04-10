# Vakyume
# A Python Library for Vacuum Design

Inspired by the 1986 edition of *Process Vacuum System Design and Operation* by James L. Ryans and Daniel L. Roper, this project offers functions related to engineering vacuums. 

We revisit the history of low pressure devices, recreating the example calculations programatically. Furthermore, a supplementary decorator is developed to calculate arbitrary missing values for a given formula. 

# Solver Decorator

As long as one keyword argument is not given, its value is calculated

E.g. 
```
@solver
def einstein(...):
	"""e = m * c ** 2"""
	return ...

einstein(e = 1000) # returns m, (1000 / 8.98755179 e16), ~1.11265 e -14
```


-[ ] Transcription of Formulae
-[ ] Filling in physics constants
-[ ] Inflation adjusting money calculations in chapter 8
-[ ] Implement solver-decorator
-[ ] Confirm worked examples in book
-[ ] Recreate some graphs in book
-[ ] Create arbitrary solver demo that scrolls through various calculations, for fun!

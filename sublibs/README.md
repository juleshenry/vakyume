# sublibs

Predecessor libraries folded into vakyume, imported with their full git history
(`git log -- sublibs/<name>`). They are kept for reference and are not part of
the vakyume package, test suite, or lint run.

- [`equation_permute`](equation_permute): exploration of equation sets and
  n-order permutations; the `.pyeqn` format and permuter that preceded vakyume's
  solver pipeline. Formerly `juleshenry/equation_permute`.
- [`kwasak`](kwasak): the original `@kwasak` missing-variable decorator that
  `vakyume/kwasak.py` is based on. Formerly `juleshenry/KWASAK`.

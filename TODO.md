# TODO

Tied arrays that are still read wrongly, found on 2026-10-02 while fixing the
tied-hash leftovers (each one gives the right answer on the same data untied):

- `group_by()` on a HoA with tied columns returns `{}` without a word.
- `kruskal_test()`, `aov()` and `oneway_test()` with tied group arrays croak
  ("all groups must contain data", "fewer than 2 complete observations",
  "observation 0 is undefined or non-numeric").
- `lm()` and `glm()` on a HoA with tied columns croak "0 degrees of freedom".
- `hoa2hoh()` with a tied key column croaks "has an undefined value at row 0".
- `binom_test()` on a tied vector croaks "successes is undef".
- `fisher_test()` on a tied AoA (the outer array tied) croaks "each row must
  be an array ref"; `chisq_test()` on the same input is right.

Most of these test a fetched cell before its get magic has run. `frame_untied()`
in `LikeR.xs` with `UNTIE_COLS`/`UNTIE_ROWS` would fix any whose result does not
share the input's columns or rows, at the cost of one copy of the tied part.

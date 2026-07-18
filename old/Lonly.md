Retired README section for the original two-argument `Lonly`.

This documented the old `$$`-prototype `Lonly` (backed by `set_difference`) that
was removed from LikeR.xs when the `Lonly` name was reassigned to the former
`get_unique`. Preserved here verbatim as it appeared in README.md.

---

## Lonly

    my @left_only = Lonly(\@left, \@right);
    my $count     = Lonly(\@left, \@right);

Takes **exactly two** array references and returns the values in the left list
that are absent from the right list. Duplicates collapse, the result keeps
left-list order, and scalar context returns the count. Values are compared by
string form (see `get_union`). A non-array-ref argument, an `undef` element,
or anything other than two references is fatal. Mirrors `List::Compare`'s
`get_Lonly`.

    my @a = (1, 2, 3, 4);
    my @b = (3, 4, 5);
    my @l = Lonly(\@a, \@b);                # (1, 2)

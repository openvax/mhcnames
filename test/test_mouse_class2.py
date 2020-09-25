from nose.tools import eq_

from mhcnames import (
    normalized_string,
    compact_string,
    parse,
    Gene,
)


def test_mouse_class2_gene():
    # H2-IAb
    gene = Gene.get("H2", "AB")
    eq_(parse("H2-IAb"), gene)
    eq_(normalized_string("H2-IAb"), "H2-AB")
    eq_(compact_string("H2-IAb"), "AB")

    # with hyphen in "H-2"
    eq_(parse("H-2-IAb"), gene)
    eq_(normalized_string("H-2-IAb"), "H2-AB")
    eq_(compact_string("H-2-IAb"), "AB")

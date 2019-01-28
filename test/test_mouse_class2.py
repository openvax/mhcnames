from nose.tools import eq_

from mhcnames import (
    normalized_string,
    compact_string,
    parse,
    NamedAllele,
)


def test_mouse_class2_alleles():
    # H2-IAb
    eq_(parse("H2-IAb"), NamedAllele("H2", "IA", "b"))
    eq_(normalized_string("H2-IAb"), "H2-IAb")
    eq_(compact_string("H2-IAb"), "IAb")

    # with hyphen in "H-2"
    eq_(parse("H-2-IAb"), NamedAllele("H2", "IA", "b"))
    eq_(normalized_string("H-2-IAb"), "H2-IAb")
    eq_(compact_string("H-2-IAb"), "IAb")

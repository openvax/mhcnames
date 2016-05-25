from nose.tools import eq_

from mhcnames import (
    normalize_allele_name,
    compact_allele_name,
    parse_allele_name,
    AlleleName,
)

def test_mouse_class2_alleles():
    # H2-IAb
    eq_(parse_allele_name("H2-IAb"),
        AlleleName("H-2", "IA", "", "b"))
    eq_(normalize_allele_name("H2-IAb"), "H-2-IAb")
    eq_(compact_allele_name("H2-IAb"), "IAb")

    # with hyphen in "H-2"
    eq_(parse_allele_name("H-2-IAb"),
        AlleleName("H-2", "IA", "", "b"))
    eq_(normalize_allele_name("H-2-IAb"), "H-2-IAb")
    eq_(compact_allele_name("H-2-IAb"), "IAb")

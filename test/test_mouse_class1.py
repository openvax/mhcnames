from nose.tools import eq_
from mhcnames import (
    parse_allele_name,
    AlleleName,
    compact_allele_name,
    normalize_allele_name
)

def test_mouse_class1_alleles_H2_Kk():
    # H2-Kk
    eq_(parse_allele_name("H2-Kk"),
        AlleleName("H-2", "K", "", "k"))
    eq_(normalize_allele_name("H2-Kk"), "H-2-Kk")
    eq_(compact_allele_name("H-2-Kk"), "Kk")

    # with a hyphen in "H-2"
    eq_(parse_allele_name("H-2-Kk"),
        AlleleName("H-2", "K", "", "k"))
    eq_(normalize_allele_name("H-2-Kk"), "H-2-Kk")
    eq_(compact_allele_name("H-2-Kk"), "Kk")

def test_mouse_class1_alleles_H2_Db():
    # H2-Db
    eq_(parse_allele_name("H2-Db"),
        AlleleName("H-2", "D", "", "b"))
    eq_(normalize_allele_name("H2-Db"), "H-2-Db")
    eq_(compact_allele_name("H2-Db"), "Db")

    # with hyphen in "H-2"
    eq_(parse_allele_name("H-2-Db"),
        AlleleName("H-2", "D", "", "b"))
    eq_(normalize_allele_name("H-2-Db"), "H-2-Db")
    eq_(compact_allele_name("H-2-Db"), "Db")

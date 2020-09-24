from nose.tools import eq_
from mhcnames import (
    parse,
    NamedAllele,
    compact_string,
    normalized_string
)


def test_mouse_class1_alleles_H2_Kk():
    H2Kk = NamedAllele.get("H2", "K", "k")

    eq_(parse("H2-Kk"), H2Kk)
    eq_(normalized_string("H2-Kk"), "H2-Kk")
    eq_(compact_string("H-2-Kk"), "Kk")

    # with a hyphen in "H-2"
    eq_(parse("H-2-Kk"), H2Kk)
    eq_(normalized_string("H-2-Kk"), "H2-Kk")
    eq_(compact_string("H-2-Kk"), "Kk")

def test_mouse_class1_alleles_H2_Db():
    H2Db = NamedAllele.get("H2", "D", "b")

    eq_(parse("H2-Db"), H2Db)
    eq_(normalized_string("H2-Db"), "H2-Db")
    eq_(compact_string("H2-Db"), "Db")

    # with hyphen in "H-2"
    eq_(parse("H-2-Db"), H2Db)
    eq_(normalized_string("H-2-Db"), "H2-Db")
    eq_(compact_string("H-2-Db"), "Db")

def test_H2_Kd_without_seps():
    eq_(parse("H2Kd"), NamedAllele.get("H2", "K", "d"))

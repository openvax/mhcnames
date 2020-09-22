from mhcnames import Haplotype, parse
from nose.tools import eq_


def test_parse_H2r():
    haplotype = parse("H2-r")
    print(haplotype)
    assert isinstance(haplotype, Haplotype)
    eq_(haplotype.normalized_string(), "H2-r")

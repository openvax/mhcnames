from mhcnames import Gene, parse
from nose.tools import eq_


def test_locus_H2K_parsing():
    H2K = Gene("H2", "K")
    eq_(parse("H2K", H2K))
    eq_(parse("H2-K", H2K))
    eq_(parse("H-2-K", H2K))
    eq_(parse("H2-K", H2K))
    eq_(parse("h2k", H2K))
    eq_(parse("h-2k", H2K))
    eq_(parse("h-2-k", H2K))
    eq_(parse("h2-k", H2K))

from mhcnames.locus import Locus
from nose.tools import eq_

def test_locus_H2K_parsing():
    H2K = Locus("H2", "K")
    eq_(Locus.parse("H2K", H2K))
    eq_(Locus.parse("H2-K", H2K))
    eq_(Locus.parse("H-2-K", H2K))
    eq_(Locus.parse("H2-K", H2K))
    eq_(Locus.parse("h2k", H2K))
    eq_(Locus.parse("h-2k", H2K))
    eq_(Locus.parse("h-2-k", H2K))
    eq_(Locus.parse("h2-k", H2K))

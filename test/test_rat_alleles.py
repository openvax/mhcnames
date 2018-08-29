
from nose.tools import eq_
from mhcnames import parse, NamedAllele

def test_rat_class1_alleles():
    eq_(parse("RT1-Bb*u"), NamedAllele("RT1", "Bb", "u"))
    eq_(parse("RT1-Bbu"), NamedAllele("RT1", "Bb", "u"))
    eq_(parse("RT1-Db1*a"), NamedAllele("RT1", "Db1", "a"))
    eq_(parse("RT1-Db1a"), NamedAllele("RT1", "Db1", "a"))
    eq_(parse("RT1-DMa*a"), NamedAllele("RT1", "DMA", "a"))
    eq_(parse("RT1-DMAa"), NamedAllele("RT1", "DMA", "a"))
    eq_(parse("RT1-9.5*f"), NamedAllele("RT1", "9.5", "f"))
    eq_(parse("RT1-9.5f"), NamedAllele("RT1", "9.5", "f"))
    eq_(parse("RT1-M3-1*av1"), NamedAllele("RT1", "M3-1", "av1"))
    eq_(parse("RT1-M3-1av1"), NamedAllele("RT1", "M3-1", "av1"))

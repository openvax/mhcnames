import mhcnames
from nose.util import eq_

def test_A_02_01_01_01():
    result = mhcnames.parse_allele_name("A*02:01:01:01")
    eq_(result.species, "HLA")
    eq_(result.gene, "A")
    eq_(result.supertype, "02")
    eq_(result.nonsyn, "01")
    eq_(result.syn, "01")
    eq_(result.intronic, "01")

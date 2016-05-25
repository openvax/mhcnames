from nose.tools import eq_
from mhcnames import parse_allele_name, AlleleName

def test_dog_class2_allele():
    eq_(parse_allele_name("DLA-DQA1*00101"),
        AlleleName("DLA", "DQA1", "01", "01"))

from nose.tools import eq_
from mhcnames import parse_allele_name, AlleleName

def test_sheep_class1_allele():
    eq_(parse_allele_name("Ovar-N*50001"),
        AlleleName("Ovar", "N", "500", "01"))

def test_sheep_class2_allele():
    eq_(parse_allele_name("Ovar-DRB1*0804"),
        AlleleName("Ovar", "DRB1", "08", "04"))

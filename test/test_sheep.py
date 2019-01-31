from nose.tools import eq_
from mhcnames import parse, FourDigitAllele


def test_sheep_class1_allele():
    eq_(
        parse("Ovar-N*50001"),
        FourDigitAllele("OLA", "N", "500", "01"))


def test_sheep_class2_allele():
    eq_(
        parse("Ovar-DRB1*0804"),
        FourDigitAllele("OLA", "DRB1", "08", "04"))

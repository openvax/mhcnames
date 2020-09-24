from nose.tools import eq_
from mhcnames import parse, FourDigitAllele


def test_sheep_class1_allele():
    eq_(
        parse("Ovar-N*50001"),
        FourDigitAllele.get("Ovar", "N", "500", "01"))

def test_sheep_class1_allele_normalize_species():
    eq_(
        parse(
            "Ovar-N*50001",
            normalize_species_prefix=True),
        FourDigitAllele.get("OLA", "N", "500", "01"))

def test_sheep_class2_allele():
    eq_(
        parse("Ovar-DRB1*0804"),
        FourDigitAllele.get("Ovar", "DRB1", "08", "04"))

def test_sheep_class2_allele_normalize_species():
    eq_(
        parse(
            "Ovar-DRB1*0804",
            normalize_species_prefix=True),
        FourDigitAllele.get("OLA", "DRB1", "08", "04"))

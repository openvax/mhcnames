from nose.tools import eq_
from mhcnames import parse, FourDigitAllele


def test_sheep_class1_allele():
    eq_(
        parse("Ovar-N*50001"),
        FourDigitAllele.get("Ovar", "N", "500", "01"))


def test_sheep_class1_allele_string_no_alias():
    eq_(
        parse("Ovar-N*50001").normalized_string(use_species_alias=False),
        "Ovar-N*500:01")

def test_sheep_class1_allele_string_yes_alias():
    eq_(
        parse("Ovar-N*50001").normalized_string(use_species_alias=True),
        "OLA-N*500:01")

def test_sheep_class2_allele():
    eq_(
        parse("Ovar-DRB1*0804"),
        FourDigitAllele.get("Ovar", "DRB1", "08", "04"))

def test_sheep_class2_allele_string_yes_alias():
    eq_(
        parse(
            "Ovar-DRB1*0804").normalized_string(use_species_alias=True),
        "OLA-DRB1*08:04")


def test_sheep_class2_allele_string_no_alias():
    eq_(
        parse(
            "Ovar-DRB1*0804").normalized_string(use_species_alias=False),
        "Ovar-DRB1*08:04")

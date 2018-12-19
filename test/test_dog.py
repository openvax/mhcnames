from nose.tools import eq_
from mhcnames import parse, FourDigitAllele, compact_string, normalized_string


def test_dog_class2_allele():
    eq_(parse("DLA-DQA1*00101"),
        FourDigitAllele("DLA", "DQA1", "001", "01"))

    eq_(normalized_string("DLA-DQA1*00101"),
        "DLA-DQA1*001:01")


def test_only_2_digits_in_first_allele_field():
    eq_(normalized_string("DLA-DQA1*0101"),
        "DLA-DQA1*001:01")


def test_alternative_species_code():
    eq_(normalized_string("Calu-DQA1*00101"),
        "Calu-DQA1*001:01")


def test_dog_class1_allele():
    expected = FourDigitAllele("DLA", "88", "508", "01")
    parsed = parse("DLA-88*50801")
    eq_(parsed, expected)

    eq_(compact_string("DLA-88*50801"), "88*50801")
    eq_(normalized_string("DLA-88*50801"), "DLA-88*508:01")

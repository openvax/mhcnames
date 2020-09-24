from nose.tools import eq_
from mhcnames import parse, FourDigitAllele, compact_string, normalized_string


def test_parse_dog_class2_allele_dla_dqa1_001_01():
    eq_(parse("DLA-DQA1*00101"),
        FourDigitAllele.get("DLA", "DQA1", "001", "01"))

def test_normalized_string_dog_class2_allele_dla_dqa1_001_01():
    eq_(normalized_string("DLA-DQA1*00101"),
        "DLA-DQA1*001:01")


def test_only_2_digits_in_first_allele_field():
    eq_(normalized_string("DLA-DQA1*0101"),
        "DLA-DQA1*001:01")

def test_species_code_calu_no_alias():
    eq_(
        normalized_string("Calu-DQA1*00101", normalize_species_prefix=False),
        "Calu-DQA1*001:01")

def test_species_code_calu_alias():
    eq_(
        normalized_string("Calu-DQA1*00101", normalize_species_prefix=True),
        "DLA-DQA1*001:01")


def test_parse_dog_class1_allele_dla_88_508_01():
    expected = FourDigitAllele.get("DLA", "88", "508", "01")
    parsed = parse("DLA-88*50801")
    eq_(parsed, expected)

def test_compact_string_dog_class1_allele_dla_88_508_01():
    eq_(compact_string("DLA-88*50801"), "88*50801")

def test_normalized_string_dog_class1_allele_dla_88_508_01():
    eq_(normalized_string("DLA-88*50801"), "DLA-88*508:01")

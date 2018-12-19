from nose.tools import eq_
from mhcnames import (
    normalized_string,
    compact_string,
)

hla_02_01_names = [
    "HLA-A*02:01",
    "HLA-A*0201",
    "A*02:01",
    "A*0201",
    "HLA-A02:01",
    # no punctuation
    "A0201",
    "HLA-A0201",
    "A0201",
    # lower case
    "hla-a*0201",
    "a*0201",
    "a*02:01",
    "a0201"
]


def test_hla_long_names():
    expected = "HLA-A*02:01"
    for name in hla_02_01_names:
        result = normalized_string(name)
        eq_(result, expected)


def test_hla_short_names():
    expected = "A0201"
    for name in hla_02_01_names:
        result = compact_string(name)
        eq_(result, expected)


def test_hla_A2_serotype():
    eq_(compact_string("A2"), "A2")
    eq_(normalized_string("A2"), "HLA-A2")

    eq_(compact_string("HLA-A2"), "A2")
    eq_(normalized_string("HLA-A2"), "HLA-A2")


def test_hla_with_3_digit_allele_code():
    # B*15:120
    eq_(normalized_string("HLA-B*15:120"), "HLA-B*15:120")
    eq_(compact_string("HLA-B*15:120"), "B15120")
    eq_(normalized_string("B15120"), "HLA-B*15:120")
    eq_(compact_string("B15120"), "B15120")

    # A*02*123
    eq_(normalized_string("HLA-A*02:123"), "HLA-A*02:123")
    eq_(compact_string("HLA-A*02:123"), "A02123")
    eq_(normalized_string("A02123"), "HLA-A*02:123")
    eq_(compact_string("A02123"), "A02123")

"""
Pig allele names from:
    Analyzing the genetic characteristics and function of the swine
    leukocyte antigen 2 gene in a Chinese inbreed of pigs
"""

from mhcnames import (
    parse_allele_name,
    AlleleName,
)
from nose.tools import eq_

def test_SLA_1_0101_with_seps():
    eq_(
        parse_allele_name("SLA-1*01:01"),
        AlleleName("SLA", "1", "01", "01"))

def test_SLA_1_HB01():
    eq_(
        parse_allele_name("SLA-1-HB01"),
        AlleleName("SLA", "1", "HB", "01"))

def test_SLA_1_0101_no_seps():
    eq_(
        parse_allele_name("SLA-10101"),
        AlleleName("SLA", "1", "01", "01"))

def test_SLA_2_07we01():
    eq_(parse_allele_name("SLA-2*07we01"),
        AlleleName("SLA", "2", "07we", "01"))

def test_SLA_2_jh01():
    eq_(parse_allele_name("SLA-2*jh01"),
        AlleleName("SLA", "2", "jh", "01"))

def test_SLA_2_w09pt22():
    eq_(parse_allele_name("SLA-2*w09pt22"),
        AlleleName("SLA", "2", "w09pt", "22"))

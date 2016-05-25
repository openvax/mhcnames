from nose.tools import eq_
from mhcnames import normalize_allele_name, compact_allele_name

def test_macaque_alleles():
    allele_name = "Mamu-B*082:02"
    eq_(normalize_allele_name(allele_name), "Mamu-B*82:02")
    eq_(compact_allele_name(allele_name), "B8202")

    # expect 3rd zero in the family "007" to be trimmed in the normalized form
    # of this allele
    allele_name = "Mamu-B*007:02"
    eq_(normalize_allele_name(allele_name), "Mamu-B*07:02")
    eq_(compact_allele_name(allele_name), "B0702")

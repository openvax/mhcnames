from nose.tools import eq_

from mhcnames import (
    normalize_allele_name,
    compact_allele_name,
)

def test_human_class2():
    expected = "HLA-DRB1*01:02"
    expected_compact = "DRB10102"
    for name in ["DRB1_0102",
                 "DRB101:02",
                 "HLA-DRB1_0102",
                 "DRB10102",
                 "DRB1*0102",
                 "HLA-DRB1*0102",
                 "HLA-DRB1*01:02"]:
        eq_(normalize_allele_name(name), expected)
        eq_(compact_allele_name(name), expected_compact)

def test_human_class2_alpha_beta():
    expected = "HLA-DPA1*01:05-DPB1*100:01"
    expected_compact = "DPA10105-DPB110001"
    for name in ["DPA10105-DPB110001",
                 "HLA-DPA1*01:05-DPB1*100:01",
                 "hla-dpa1*0105-dpb1*10001",
                 "dpa1*0105-dpb1*10001",
                 "HLA-DPA1*01:05/DPB1*100:01"]:
        eq_(normalize_allele_name(name), expected)
        eq_(compact_allele_name(name), expected_compact)

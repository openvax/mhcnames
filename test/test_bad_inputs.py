from nose.tools import raises
from mhcnames import normalize_allele_name, AlleleParseError

@raises(AlleleParseError)
def test_extra_text_after_allele():
    normalize_allele_name("HLA-A*02:01 zipper")
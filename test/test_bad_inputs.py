from nose.tools import raises
from mhcnames import normalized_string, ParseError

@raises(ParseError)
def test_extra_text_after_allele():
    normalized_string("HLA-A*02:01 zipper")

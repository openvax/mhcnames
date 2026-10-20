import pytest
from mhcnames import normalize_allele_name, AlleleParseError

def test_extra_text_after_allele():
    with pytest.raises(AlleleParseError):
        normalize_allele_name("HLA-A*02:01 zipper")

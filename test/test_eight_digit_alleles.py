import mhcnames
from nose.tools import eq_


def test_A_02_01_01_01():
    result = mhcnames.parse("A*02:01:01:01")
    eq_(result.species_prefix, "HLA")
    eq_(result.gene_name, "A")
    eq_(result.group_id, "02")
    eq_(result.protein_id, "01")
    eq_(result.coding_sequence_id, "01")
    eq_(result.genomic_sequence_id, "01")

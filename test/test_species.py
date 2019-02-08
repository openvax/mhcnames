from mhcnames import Species
from nose.tools import eq_


def test_HLA():
    for hla in ["HLA", "hla"]:
        species = Species.get(hla)
        assert species is not None, "Could not find Species for '%s'" % hla
        eq_(species.common_species_name, "human")


def test_human():
    for human in ["human", "Human", "HUMAN"]:
        species = Species.get(human)
        assert species is not None, "Could not find Species for '%s'" % human
        eq_(species.common_species_name, "human")
from __future__ import print_function, division, absolute_import

import yaml
from os.path import dirname, join

from .normalizing_dictionary import NormalizingDictionary

package_dir = dirname(__file__)
data_dir = join(package_dir, "data")


def get_path(yaml_filename):
    return join(data_dir, yaml_filename)


def load(
        yaml_filename,
        normalize_first_level_keys=False,
        normalize_second_level_keys=False):
    path = get_path(yaml_filename)
    with open(path, 'r') as f:
        result = yaml.load(f)

    if normalize_first_level_keys:
        result = NormalizingDictionary.from_dict(result)

    if normalize_second_level_keys:
        # turn first layer of values in dictionary into NormalizingDictionary
        result = result.map_values(NormalizingDictionary.from_dict)

    return result


# dictionary mapping group -> scientific name -> info dictionary
species = load("species.yaml")

# Dictionary mapping species prefix to MHC class to genes
# For Class II genes, there is an intermediate key indicating gene group
# Examples:
#   gene_ontology["HLA"]["Ia"]
#   gene_ontology["HLA"]["IIa]["DR"]
gene_ontology = load("gene_ontology.yaml", normalize_first_level_keys=True)

# Dictionary mapping each species prefix to a dictionary from old gene
# names to new ones
gene_aliases = load(
    "gene_aliases.yaml",
    normalize_first_level_keys=True,
    normalize_second_level_keys=True)

# Dictionary mapping each species prefix to a dictionary from
# old/retired/provisional allele names to standard new names
allele_aliases = load(
    "allele_aliases.yaml",
    normalize_first_level_keys=True,
    normalize_second_level_keys=True)

# Dictionary mapping species to haplotype name to list of alleles
haplotypes = load(
    "haplotypes.yaml",
    normalize_first_level_keys=True,
    normalize_second_level_keys=True)

# Dictionary mapping species to serotype name to list of alleles
serotypes = load(
    "serotypes.yaml",
    normalize_first_level_keys=True,
    normalize_second_level_keys=True)

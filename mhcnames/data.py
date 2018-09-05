from __future__ import print_function, division, absolute_import

import yaml
from os.path import dirname, join

package_dir = dirname(__file__)
data_dir = join(package_dir, "data")

def get_path(yaml_filename):
    return join(data_dir, yaml_filename)

def expand_keys_with_aliases(dictionary):
    """
    Fill a dictionary with duplicate entries for "aliases" of
    the keys, which consist of:
        (1) making a key uppercase
        (2) removing dashes from a key
        (3) removing dahses and making remaining characters uppercase
    """
    expanded_dictionary = {}
    for (original_key, value) in dictionary.items():
        # recursively traverse dictionaries
        if isinstance(value, dict):
            value = expand_keys_with_aliases(value)
        keys = {
            original_key,
            original_key.upper(),
            original_key.replace("-", ""),
            original_key.upper().replace("-", "")
        }
        for key in keys:
            expanded_dictionary[key] = value
    return expanded_dictionary


def load(yaml_filename):
    path = get_path(yaml_filename)
    with open(path, 'r') as f:
        raw_dictionary = yaml.load(f)
    expanded_dictionary = expand_keys_with_aliases(raw_dictionary)
    return raw_dictionary, expanded_dictionary

gene_ontology, gene_ontology_with_uppercase_keys = load("gene_ontology.yaml")
serotypes, serotypes_with_uppercase_keys = load("serotypes.yaml")
allele_aliases, allele_aliases_with_uppercase_keys = load("allele_aliases.yaml")
gene_aliases, gene_aliases_with_uppercase_keys = load("gene_aliases.yaml")
haplotypes, haplotypes_with_uppercase_keys = load("haplotypes.yaml")

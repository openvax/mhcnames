from __future__ import print_function, division, absolute_import

import yaml
from os.path import dirname, join

package_dir = dirname(__file__)
data_dir = join(package_dir, "data")

def get_path(yaml_filename):
    return join(data_dir, yaml_filename)

def load(yaml_filename):
    path = get_path(yaml_filename)
    with open(path, 'r') as f:
        return yaml.load(f)

gene_ontology = load("gene_ontology.yaml")
serotypes = load("serotypes.yaml")
allele_aliases = load("allele_aliases.yaml")
gene_aliases = load("gene_aliases.yaml")
haplotypes = load("haplotypes.yaml")
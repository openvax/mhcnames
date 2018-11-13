from __future__ import print_function, division, absolute_import

import yaml
from os.path import dirname, join
from collections import defaultdict

package_dir = dirname(__file__)
data_dir = join(package_dir, "data")


def get_path(yaml_filename):
    return join(data_dir, yaml_filename)


def load(yaml_filename):
    path = get_path(yaml_filename)
    with open(path, 'r') as f:
        return yaml.load(f)

# Dictionary mapping species prefix to MHC class to genes
# For Class II genes, there is an intermediate key indicating gene group
# Examples:
#   gene_ontology["HLA"]["Ia"]
#   gene_ontology["HLA"]["IIa]["DR"]
gene_ontology = load("gene_ontology.yaml")

# Dictionary mapping each species prefix to a dictionary from old gene
# names to new ones
gene_aliases = load("gene_aliases.yaml")

# Dictionary mapping each species prefix to a dictionary from
# old/retired/provisional allele names to standard new names
allele_aliases = load("allele_aliases.yaml")

# Dictionary mapping species to haplotype name to list of alleles
haplotypes = load("haplotypes.yaml")

# Dictionary mapping species to serotype name to list of alleles
serotypes = load("serotypes.yaml")


def expand_with_uppercase_and_no_dash(*names):
    """
    Returns set of names which includes original input, uppercase,
    name with all dashes removed, and dash-free uppercase.
    """
    result_set = set([])
    for name in names:
        result_set.add(name)
        no_dash = name.replace("-", "")
        result_set.add(no_dash)
        result_set.add(name.upper())
        result_set.add(no_dash.upper())
    return result_set


def _create_uppercase_and_no_dash_allele_aliases():
    result = defaultdict(dict)

    for (species, aliases) in allele_aliases.items():
        for (old_name, new_name) in aliases.items():
            for key in expand_with_uppercase_and_no_dash(old_name, new_name):
                result[species][key] = new_name
    return result

allele_aliases_with_uppercase_and_no_dash = \
    _create_uppercase_and_no_dash_allele_aliases()


def _create_expanded_allele_aliases():
    """
    Reorganizes allele_aliases dict to contain
           species prefix -> gene name -> allele name -> normalized allele name
    and creates new aliases by including uppercase, and dash-free versions
    of old allele names.
    """
    result = {}
    for species, original_alias_dict in allele_aliases.items():
        gene_to_alias_dict = defaultdict(dict)
        for (original_old_name, new_name) in original_alias_dict.items():
            old_gene_name, old_name_without_gene = original_old_name.split("*")
            new_gene_name, new_name_without_gene = new_name.split("*")
            assert old_gene_name == new_gene_name, (original_old_name, new_name)
            for alias in expand_with_uppercase_and_no_dash(old_name_without_gene):
                gene_to_alias_dict[old_gene_name][alias] = new_name_without_gene
        result[species] = gene_to_alias_dict
    return result

# dictrionary mapping species prefix -> gene name -> allele alias -> normalized name
species_to_gene_to_allele_aliases = _create_expanded_allele_aliases()


def _create_species_to_genes_dict():
    """
    Flatten the gene ontology dictionary to just map each species
    to a list of genes.
    """
    species_to_genes = {}
    for species, species_ontology_dict in gene_ontology.items():
        species_genes = set([])
        if species_ontology_dict is None:
            raise ValueError("Missing gene ontology for '%s'" % species)
        if not isinstance(species_ontology_dict, dict):
            raise ValueError("Expected ontology dictionary for species '%s' but got: %s" % (
                species, species_ontology_dict))
        for mhc_class, genes_or_dict in species_ontology_dict.items():
            if mhc_class.startswith("II"):
                if not isinstance(genes_or_dict, dict):
                    raise ValueError(
                        "Expected MHC class '%s' of '%s' to contain dictionary" % (
                            mhc_class, species))
                # Class II MHC class entries are dictionaries
                # mapping e.g. "DQ->{DQA, DQB1}"
                for genes in genes_or_dict.values():
                    species_genes.update(genes)
            else:
                species_genes.update(genes_or_dict)
        species_to_genes[species] = species_genes
    return species_to_genes


# dictionary mapping species prefix -> set of gene names
species_to_genes = _create_species_to_genes_dict()


def _create_expanded_gene_aliases():
    """
    Returns species -> alias gene -> normalized gene dictionary
    which includes uppercase and dash-free versions of alias gene
    names, as well as the normalized names.

    Also apply uppercase/nodash expansion to all gene names in
    the ontology dictionary.
    """
    result = defaultdict(dict)
    for species_prefix, species_gene_alias_dict in gene_aliases.items():
        for old_name, new_name in species_gene_alias_dict.items():
            for alias in expand_with_uppercase_and_no_dash(old_name):
                result[species_prefix][alias] = new_name
    # apply uppercase and no-dash string transformations to all genes
    for species, genes in species_to_genes.items():
        for gene in genes:
            for alias in expand_with_uppercase_and_no_dash(gene):
                result[species_prefix][alias] = gene
    return result

# dictionary mapping species prefix -> gene alias -> canonical name
species_to_gene_aliases = _create_expanded_gene_aliases()


def _create_serotype_aliases_dict():
    species_to_alias_to_serotype = defaultdict(dict)
    for (species, species_serotypes) in serotypes.items():
        for serotype_name in species_serotypes:
            aliases = {serotype_name}
            if serotype_name.startswith("Cw"):
                aliases.add("C" + serotype_name[2:])
            if serotype_name.startswith("DPw"):
                aliases.add("DP" + serotype_name[3:])
            aliases.add(serotype_name.upper())
            for alias in aliases:
                species_to_alias_to_serotype[species][alias] = serotype_name
    return species_to_alias_to_serotype

species_to_alias_to_serotype = _create_serotype_aliases_dict()


def get_serotype(species_prefix, serotype_name):
    """
    Returns either None (if serotype doesn't exist) or tuple with following
    entries:
        - normalized species prefix
        - normalized serotype name
        - list of alleles in serotype
    """
    for candidate_species_prefix in expand_with_uppercase_and_no_dash(species_prefix):
        serotype_aliases = species_to_alias_to_serotype.get(candidate_species_prefix, {})
        for serotype_alias in expand_with_uppercase_and_no_dash(serotype_name):
            true_serotype_name = serotype_aliases.get(serotype_alias)
            if true_serotype_name is not None:
                allele_list = serotypes[candidate_species_prefix][true_serotype_name]
                return (
                    candidate_species_prefix,
                    true_serotype_name,
                    allele_list
                )
    return None

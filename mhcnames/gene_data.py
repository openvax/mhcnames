# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import


from .data import gene_ontology, gene_aliases
from .normalizing_dictionary import NormalizingDictionary


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
            if isinstance(genes_or_dict, dict):
                for genes in genes_or_dict.values():
                    species_genes.update(genes)
            else:
                assert isinstance(genes_or_dict, (list, tuple, set))
                species_genes.update(genes_or_dict)
        species_to_genes[species] = species_genes
    return species_to_genes

# dictionary mapping species prefix -> set of gene names
species_prefix_to_genes = _create_species_to_genes_dict()


def _create_expanded_gene_aliases():
    """
    Returns species -> alias gene -> normalized gene dictionary
    which includes uppercase and dash-free versions of alias gene
    names, as well as the normalized names.

    Also apply uppercase/nodash expansion to all gene names in
    the ontology dictionary.
    """
    result = NormalizingDictionary(default_value_fn=NormalizingDictionary)
    for species_prefix, species_gene_alias_dict in gene_aliases.items():
        for old_name, new_name in species_gene_alias_dict.items():
            result[species_prefix][old_name] = new_name

    # apply uppercase and no-dash string transformations to all genes
    for species, genes in species_prefix_to_genes.items():
        for gene in genes:
            result[species][gene] = gene

    return result

# dictionary mapping species prefix -> gene alias -> canonical name
species_prefix_to_gene_alias_to_gene = _create_expanded_gene_aliases()

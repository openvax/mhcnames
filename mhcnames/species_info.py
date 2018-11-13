# Copyright (c) 2018. Mount Sinai School of Medicine
#
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

from .mhc_class import class1_subtypes, class2_subtypes


class SpeciesInfo(object):
    def __init__(
            self,
            prefix,
            gene_ontology,
            gene_aliases={},
            allele_aliases={},
            haplotypes={},
            serotypes={}):
        """
        Parameters
        ----------
        prefix : str

        gene_ontology : dict
            Dictionary containing keys which are a subset of
            {Ia, Ib, Ic, IIa, IIb}. The class I entries map to lists of genes
            and the class II entries map to dictionaries such as:
            {"DR": ["DRA", "DRB1", "DRB3", "DRB5"]}

        gene_aliases : dict
            Dictionary mapping non-canonical gene names to their canonical
            forms.

        allele_aliases : dict
            Dictionary mapping non-canonical allele names to their
            canonical forms.

        haplotypes : dict
            Dictionary mapping haplotype names to lists of alleles (one per gene)

        serotypes : dict
            Dictionary mapping serotype names to lists of alleles
            (all expected to be for the same gene)
        """

        self.prefix = prefix
        self.gene_ontology = gene_ontology
        self.gene_aliases = gene_aliases
        self.allele_aliases = allele_aliases
        self.serotypes = serotypes
        self.haplotypes = haplotypes

        self._genes = None
        self._gene_set = None
        self._expanded_gene_alias_dict = None

    def __str__(self):
        return "SpeciesInfo(prefix='%s')" % (self.prefix,)

    @property
    def genes(self):
        if self._genes is None:
            self._genes = self._collect_genes()
        return self._genes

    @property
    def gene_set(self):
        if self._gene_set is None:
            self._gene_set = set(self.genes)
        return self._gene_set

    @property
    def expanded_gene_alias_dict(self):
        if self._expanded_gene_alias_dict is None:
            self._expanded_gene_alias_dict = self._create_expanded_gene_aliases()
        return self._expanded_gene_alias_dict

    def _collect_genes(self):
        """
        Create a dictionary mapping each species to a set of genes
        """
        all_genes = []
        for class1_category in class1_subtypes:
            class1_genes = self.gene_ontology.get(class1_category, [])
            all_genes.extend(class1_genes)
        for class2_category in class2_subtypes:
            class2_genes_dict = self.gene_ontology.get(class2_category, {})
            for class2_genes in class2_genes_dict.values():
                all_genes.extend(class2_genes)
        return all_genes

    def find_matching_gene_name(self, gene_name):
        """
        Use known aliases and normalized capitlization to infer
        the canonical gene name corresponding to the input.
        """
        if gene_name in self.gene_set:
            return gene_name
        return self.expanded_gene_alias_dict.get(gene_name.upper())

    def normalize_gene_name_if_exists(self, gene_name):
        normalized_name = self.find_matching_gene_name(gene_name)
        if normalized_name:
            return normalized_name
        else:
            return gene_name

    def _create_expanded_gene_aliases(self):
        expanded_aliases = {}
        for gene in self.genes:
            expanded_aliases[gene] = gene
            upper = gene.upper()
            if upper != gene:
                expanded_aliases[upper] = gene

        for original_alias, gene_name in self.gene_aliases.items():
            expanded_alias_set = {
                original_alias,
                original_alias.replace("-", ""),
                original_alias.upper(),
                original_alias.replace("-", "").upper()
            }
            for alias in expanded_alias_set:
                expanded_aliases[alias] = gene_name

        return expanded_aliases

    def get_mhc_class(self, gene_name):
        """
        Returns either one of "Ia", "Ib", "IIa", "IIb" or None
        if species can't be found
        """
        gene_name = self.find_matching_gene_name(gene_name)
        for mhc_class, mhc_class_members in self.gene_ontology.items():
            if mhc_class in class1_subtypes:
                for member_gene in mhc_class_members:
                    if member_gene == gene_name:
                        return mhc_class
            elif mhc_class in class2_subtypes:
                for locus, genes in mhc_class_members.items():
                    if locus == gene_name:
                        return mhc_class
                    for candidate_gene_name in genes:
                        if candidate_gene_name == gene_name:
                            return mhc_class
        return None

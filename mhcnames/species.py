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

from collections import OrderedDict

from pytypes import typechecked

from .parsed_result import ParsedResult
from .data import gene_ontology as raw_gene_ontology_dict
from .data import gene_aliases as raw_gene_aliases_dict
from .data import serotypes as raw_serotypes_dict
from .data import haplotypes as raw_haplotypes_dict
from .data import allele_aliases as raw_allele_aliases_dict
from .normalizing_dictionary import NormalizingDictionary

from .species_data import (
    prefix_to_default_common_name,
    prefix_to_scientific_name,
    common_names_to_scientific_names,
    normalize_species_prefix,
    scientific_name_to_prefixes,
    scientific_name_to_canonical_prefix,
    prefix_to_alias,
)


class Species(ParsedResult):
    """
    Representation of a parsed species prefix such as "HLA", "ELA"
    """
    @typechecked
    def __init__(self, species_prefix : str):
        self.species_prefix = normalize_species_prefix(species_prefix)
        self._genes = None
        self._gene_set = None
        self._expanded_gene_alias_dict = None
        self._serotypes = None
        self._allele_aliases = None
        self._gene_aliases = None
        self._haplotypes = None

    def field_names(self):
        return ("species_prefix",)

    @property
    def prefix(self):
        return self.species_prefix

    @property
    def historic_alias(self):
        return prefix_to_alias.get(self.prefix, self.prefix)

    def normalized_string(self, include_species=True, use_species_alias=True):
        if not include_species:
            return ""
        elif use_species_alias:
            return self.historic_alias
        else:
            return self.prefix

    def compact_string(self, include_species=False, use_species_alias=True):
        return self.normalized_string(
            include_species=include_species, use_species_alias=use_species_alias)

    @classmethod
    def get(cls, species_name):
        """
        Alias for find_matching_species function
        """
        if species_name.__class__ is Species:
            return species_name
        return find_matching_species(species_name)

    def to_record(self):
        return OrderedDict([
            ("species_prefix", self.species_prefix),
            ("species_name", self.common_species_name),
        ])

    @property
    def common_species_name(self):
        """
        Returns common species name associated with MHC species
        prefix.
        """
        return prefix_to_default_common_name.get(self.species_prefix).lower()

    @property
    def scientific_species_name(self):
        return prefix_to_scientific_name.get(self.species_prefix).lower()

    @property
    def all_prefixes(self):
        """
        Returns all prefixes used for this species, including aliases
        such as SLA/Susc
        """
        return scientific_name_to_prefixes.get(self.scientific_species_name, [])

    @property
    def gene_ontology(self):
        """
        Dictionary containing keys which are a subset of
        {I, Ia, Ib, Ic, Id, II, IIa, IIb}.
        The class I entries map to lists of genes and the class II entries map
        to dictionaries such as:
            {"DR": ["DRA", "DRB1", "DRB3", "DRB5"]}
        """
        for prefix in self.all_prefixes:
            if prefix in raw_gene_ontology_dict:
                return raw_gene_ontology_dict[prefix]
        return {}

    @property
    def gene_aliases(self):
        """
        Dictionary mapping non-canonical gene names to their canonical forms.
        """
        if self._gene_aliases is None:
            gene_aliases = NormalizingDictionary()
            for prefix in self.all_prefixes:
                gene_aliases.update(raw_gene_aliases_dict.get(prefix, {}))
            self._gene_aliases = gene_aliases
        return self._gene_aliases

    @property
    def allele_aliases(self):
        """
        Dictionary mapping non-canonical allele names to their canonical forms.
        """
        if self._allele_aliases is None:
            allele_aliases = NormalizingDictionary()
            for prefix in self.all_prefixes:
                allele_aliases.update(raw_allele_aliases_dict.get(prefix, {}))
            self._allele_aliases = allele_aliases
        return self._allele_aliases

    @property
    def serotypes(self):
        """
        Dictionary mapping serotype names to lists of alleles
        (all expected to be for the same gene)
        """
        if self._serotypes is None:
            serotypes = NormalizingDictionary()
            for prefix in self.all_prefixes:
                serotypes.update(raw_serotypes_dict.get(prefix, {}))
            self._serotypes = serotypes
        return self._serotypes

    @property
    def haplotypes(self):
        """
        Dictionary mapping haplotype names to lists of alleles (one per gene)
        """
        if self._haplotypes is None:
            haplotypes = NormalizingDictionary()
            for prefix in self.all_prefixes:
                haplotypes.update(raw_haplotypes_dict(prefix, {}))
            self._haplotypes = haplotypes
        return self._haplotypes

    def genes(self):
        if self._genes is None:
            self._genes = self._collect_genes()
        return self._genes

    def gene_set(self):
        if self._gene_set is None:
            self._gene_set = set(self.genes())
        return self._gene_set

    def expanded_gene_aliases(self):
        if self._expanded_gene_alias_dict is None:
            self._expanded_gene_alias_dict = self._create_expanded_gene_aliases()
        return self._expanded_gene_alias_dict

    def _collect_genes(self):
        """
        Create a dictionary mapping each species to a set of genes
        """
        all_genes = []
        for class1_category in {"I", "Ia", "Ib", "Ic", "Id"}:
            class1_genes = self.gene_ontology.get(class1_category, [])
            all_genes.extend(class1_genes)
        for class2_category in {"II", "IIa", "IIb"}:
            class2_genes_dict = self.gene_ontology.get(class2_category, {})
            if type(class2_genes_dict) is not dict:
                raise ValueError(
                    "Malformed class II gene ontology for '%s', got %s for '%s'" % (
                        self.species_prefix,
                        type(class2_genes_dict),
                        class2_category))
            for class2_genes in class2_genes_dict.values():
                all_genes.extend(class2_genes)
        return all_genes

    def find_matching_gene_name(self, gene_name):
        """
        Use known aliases and normalized capitlization to infer
        the canonical gene name corresponding to the input.
        """
        if gene_name in self.gene_set():
            return gene_name
        return self.expanded_gene_aliases().get(gene_name.upper())

    def normalize_gene_name_if_exists(self, gene_name):
        normalized_name = self.find_matching_gene_name(gene_name)
        if normalized_name:
            return normalized_name
        else:
            return gene_name

    def _create_expanded_gene_aliases(self):
        expanded_aliases = {}
        for gene in self.genes():
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

    def get_mhc_class_of_gene(self, gene_name):
        """
        Parameters
        ----------
        gene_name : str

        Returns either one of "I", "Ia", "Ib", "Ic", "Id", "II", IIa", "IIb"
        or None if species can't be found
        """
        gene_name = self.find_matching_gene_name(gene_name)
        for mhc_class, mhc_class_members in self.gene_ontology.items():
            if mhc_class in {"I", "Ia", "Ib", "Ic", "Id"}:
                for member_gene in mhc_class_members:
                    if member_gene == gene_name:
                        return mhc_class
            elif mhc_class in {"II", "IIa", "IIb"}:
                for locus, genes in mhc_class_members.items():
                    if locus == gene_name:
                        return mhc_class
                    for candidate_gene_name in genes:
                        if candidate_gene_name == gene_name:
                            return mhc_class
        return None


# map prefix strings to Species objects
_species_cache = {}

def find_matching_species(name):
    """
    Given an unnormalized species prefix string,
    returns Species object.
    """
    if name not in _species_cache:
        if name in common_names_to_scientific_names:
            scientific_name = common_names_to_scientific_names[name]
        elif name in prefix_to_scientific_name:
            scientific_name = prefix_to_scientific_name[name]
        else:
            # fall back: assume input is already a scientific name,
            # though this is likely to fail at next step,
            # causing returned value to be None
            scientific_name = name

        prefix = scientific_name_to_canonical_prefix.get(scientific_name)
        if prefix is None:
            species = None
        else:
            species = Species(prefix)
        _species_cache[name] = species
    return _species_cache[name]


def find_matching_species_prefix(name):
    """
    Given an unnormalized species prefix string,
    returns normalized prefix string.
    """
    species = find_matching_species(name)
    if species:
        return species.species_prefix
    else:
        return None


def infer_species_prefix_substring(name):
    """
    Trying to parse prefixes of alleles such as:
        HLA-A
    but also ones with dashes in the species prefix:
        H-2-K
    and also those lacking any dashes such as:
        H2K

     ...we also need to consider that alleles, haplotypes, etc may come
     immediately after the gene:
        H2Kk
        HLA-A0201

    Returns the normalized species prefix and the original string that matched
    it or None.
    """
    # Try parsing a few different substrings to get the species,
    # and then use the species gene list to determine what the gene is in this string
    candidate_species_substrings = [name]

    if "-" in name:
        # if name is "H-2-K" then try parsing "H" and "H-2" as a species
        # prefix
        parts_split_by_dash = name.split("-")
        candidate_species_substrings.extend([
            parts_split_by_dash[0],
            parts_split_by_dash[0] + "-" + parts_split_by_dash[1]
        ])
    for seq in candidate_species_substrings:
        for n in [None, 4, 3, 2]:
            original_prefix = seq[:n]
            normalized_prefix = find_matching_species_prefix(name[:n])
            if normalized_prefix is not None:
                return normalized_prefix, original_prefix
    return None

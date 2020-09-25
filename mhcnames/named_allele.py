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

from pytypes import typechecked

from .gene import Gene
from .parsed_result import ParsedResult

class NamedAllele(ParsedResult):
    """
    Some species, such as mouse (H2) and rats (RT1) do not yet use
    the standard nomenclature format (e.g. Species-Gene*Group:Protein).

    Instead they identify each distinct MHC molecule with a letter, or sometimes
    short sequence. Example: H2-Kk, RT1-9.5*f

    Also, some older swine (SLA) alleles seem to have not been translated into
    the updated nomenclature (e.g. SLA-1-CHANGDA)
    """

    @typechecked
    def __init__(self, gene : Gene, allele_name : str):
        self.gene = gene
        self.allele_name = allele_name

    def field_names(self):
        return ("gene", "allele_name")

    @property
    def species(self):
        return self.gene.species

    @property
    def species_prefix(self):
        return self.species.prefix

    @property
    def gene_name(self):
        return self.gene.name

    @property
    def mhc_class(self):
        return self.gene.mhc_class

    @property
    def is_class1(self):
        return self.gene.is_class1

    @property
    def is_class2(self):
        return self.gene.is_class2

    @classmethod
    def get(cls, species_prefix, gene_name, allele_name):
        gene = Gene.get(species_prefix, gene_name)
        if gene is None:
            return None
        return NamedAllele(gene, allele_name)

    def normalized_string(
            self,
            include_species=True,
            use_species_alias=True):
        return "%s%s" % (
            self.gene.normalized_string(
                include_species=include_species,
                use_species_alias=use_species_alias),
            self.allele_name)

    def compact_string(
            self,
            include_species=False,
            use_species_alias=True):
        return "%s%s" % (
            self.gene.compact_string(
                include_species=include_species,
                use_species_alias=use_species_alias),
            self.allele_name)

    def to_record(self):
        d = self.gene.to_record()
        d["allele"] = self.normalized_string()
        d["is_mutant"] = False

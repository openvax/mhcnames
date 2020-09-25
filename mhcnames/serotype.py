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

from typing import List, Union

from pytypes import typechecked

from .gene import Gene
from .species import Species
from .named_allele import NamedAllele
from .numeric_alleles import FourDigitAllele
from .parsed_result import ParsedResult

class Serotype(ParsedResult):

    @typechecked
    def __init__(
            self, species : Species,
            name : str, alleles : Union[List[FourDigitAllele], List[NamedAllele]]):
        if len(alleles) == 0:
            raise ValueError("Cannot create Serotype without alleles")

        genes = list({allele.gene for allele in alleles})
        if len(genes) != 1:
            raise ValueError(
                "Serotype cannot span multiple genes: %s" % (
                    [gene.name for gene in genes],))
        gene = genes[0]
        if gene.species != species:
            raise ValueError(
                "Species inferred from given alleles (%s) is different from %s" % (
                    gene.species,
                    species,
                ))
        self.species = species
        self.gene = gene
        self.name = name
        self.alleles = alleles

    @property
    def species_prefix(self):
        return self.species.prefix

    @property
    def gene_name(self):
        return self.gene.name

    @classmethod
    def field_names(cls):
        return ("species", "name", "alleles")

    def normalized_string(
            self,
            include_species=True,
            use_species_alias=True):
        if include_species:
            return "%s-%s" % (
                self.species.normalized_string(
                    use_species_alias=use_species_alias),
                self.name)
        else:
            return self.name

    def compact_string(self, include_species=False, use_species_alias=True):
        return self.normalized_string(
            include_species=include_species,
            use_species_alias=use_species_alias)

    def to_record(self):
        d = Gene.to_record(self)
        d["serotype"] = self.normalized_string()
        d["alleles_in_serotype"] = ";".join([
            allele.normalized_string() for allele in self.alleles
        ])
        return d

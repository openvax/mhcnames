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
from .mutation import Mutation
from .named_allele import NamedAllele
from .numeric_alleles import FourDigitAllele
from .parsed_result import ParsedResult


class MutantAllele(ParsedResult):
    @typechecked
    def __init__(
            self,
            original_allele : Union[NamedAllele, FourDigitAllele],
            mutations : List[Mutation]):
        self.original_allele = original_allele
        self.mutations = mutations

    @classmethod
    def field_names(cls):
        return (
            "original_allele",
            "mutations",
        )

    @property
    def species(self):
        return self.original_allele.species

    @property
    def species_prefix(self):
        return self.species.prefix

    @property
    def gene_name(self):
        return self.gene.name

    @property
    def gene(self):
        return self.original_allele.gene

    @property
    def gene_name(self):
        return self.gene.name

    def to_dict(self):
        d = Gene.to_dict(self)
        d["original_allele"] = self.original_allele.to_dict()
        d["original_allele_class"] = self.original_allele_class
        d["mutations"] = [m.to_dict() for m in self.mutations]
        return d

    @classmethod
    def from_dict(cls, d):
        original_allele_class = d.pop("original_allele_class")
        if original_allele_class == "NamedAllele":
            original_allele = NamedAllele.from_dict(d)
        elif original_allele_class == "FourDigitAllele":
            original_allele = FourDigitAllele.from_dict(d)
        else:
            raise ValueError(
                "Unable to create MutantAllele from %s" % original_allele_class)
        mutations = [Mutation.from_dict(m) for m in d["mutations"]]
        return cls(original_allele, mutations)

    def mutation_string(self):
        return " ".join([mut.normalized_string() for mut in self.mutations])

    def normalized_string(self, include_species=True):
        return "%s %s mutant" % (
            self.original_allele.normalized_string(
                include_species=include_species),
            self.mutation_string())

    def compact_string(self, include_species=False):
        return "%s %s mutant" % (
            self.original_allele.compact_string(
                include_species=include_species),
            self.mutation_string())

    def to_record(self):
        d = self.original_allele.to_dict()
        d["allele"] = self.normalized_string()
        d["is_mutant"] = True
        d["mutations"] = self.mutation_string()
        return d

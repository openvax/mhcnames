# Copyright (c) 2018-2019. Mount Sinai School of Medicine
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

from .four_digit_allele import FourDigitAllele
from .mutation import Mutation


class MutantAllele(FourDigitAllele):
    def __init__(self, original_allele, mutations):
        self.original_allele = original_allele
        FourDigitAllele.__init__(
            self,
            species_prefix=original_allele.species_prefix,
            gene_name=original_allele.gene_name,
            group_id=original_allele.group_id,
            protein_id=original_allele.protein_id,
            modifier=original_allele.modifier)
        self.mutations = mutations

    @classmethod
    def field_names(cls):
        return (
            "species_prefix",
            "gene_name",
            "group_id",
            "protein_id",
            "modifier",
            "mutations",
        )

    def to_dict(self):
        d = FourDigitAllele.to_dict(self)
        d["mutations"] = [m.to_dict() for m in self.mutations]
        return d

    @classmethod
    def from_dict(cls, d):
        four_digit_allele = FourDigitAllele.from_dict(d)
        mutations = [Mutation.from_dict(m) for m in d["mutations"]]
        return cls(four_digit_allele, mutations)

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

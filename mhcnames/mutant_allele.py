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

from .locus import Locus

class MutantFourDigitAllele(Locus):
    def __init__(self, original_allele, mutations):
        Locus.__init__(
            self,
            species_prefix=original_allele.species_prefix,
            gene_name=original_allele.gene_name)
        self.original_allele = original_allele
        self.mutations = mutations

    def is_mutant(self):
        return True

    def mutation_string(self):
        return " ".join([mut.normalized_string() for mut in self.mutations])

    def normalized_string(self, include_species_prefix=True):
        return "%s %s mutant" % (
            self.original_allele.normalized_string(
                include_species_prefix=include_species_prefix),
            self.mutation_string())

    def compact_string(self, include_species_prefix=False):
        return "%s %s mutant" % (
            self.original_allele.compact_string(
                include_species_prefix=include_species_prefix),
            self.mutation_string())

    def to_dict(self):
        d = self.original_allele.to_dict()
        d["allele_name"] = self.normalized_string()
        d["mutations"] = self.mutation_string()
        return d

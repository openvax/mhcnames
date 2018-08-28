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

from .four_digit_allele import FourDigitAllele

class MutantFourDigitAllele(FourDigitAllele):
    def __init__(self, original_allele, mutations):
        FourDigitAllele.__init__(
            self,
            original_allele.species_prefix,
            original_allele.gene_name,
            original_allele.group_id,
            original_allele.protein_id,
            original_allele.modifier)
        self.mutations = mutations

    def is_mutant(self):
        return True

    def to_dict(self):
        d = FourDigitAllele.to_dict(self)
        d["mutant_allele"] = self.normalized_string()
        d["mutations"] = ";".join([
            mut.normalized_string() for mut in self.mutations])
        return d

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


class SixDigitAllele(FourDigitAllele):
    def __init__(
            self,
            species_prefix,
            gene_name,
            group_id,
            protein_id,
            coding_sequence_id,
            modifier=None):
        FourDigitAllele.__init__(
            self,
            species_prefix,
            gene_name,
            group_id,
            protein_id,
            modifier)
        self.coding_sequence_id = coding_sequence_id

    @classmethod
    def field_names(cls):
        return (
            "species_prefix",
            "gene_name",
            "group_id",
            "protein_id",
            "coding_sequence_id",
            "modifier"
        )

    def normalized_string(self, include_species=True, include_modifier=True):
        """
        Return allele strings like "HLA-A*02:01"
        """
        result = "%s:%s" % (
            FourDigitAllele.normalized_string(
                self,
                include_species=include_species,
                include_modifier=False),
            self.coding_sequence_id)
        if include_modifier and self.modifier:
            result += self.modifier
        return result

    def compact_string(self, include_species=False):
        """
        Compact representation of a SixDigitAllele, omits the "*" and ":"
        in allele names.
            Normalized: HLA-A*02:01:01
            Compact: HLA-A020101
        """
        return "%s%s" % (
            FourDigitAllele.compact_string(include_species=include_species),
            self.coding_sequence_id)

    def to_six_digit_allele(self):
        """
        Descendant class EightDigitAllele can use this
        to project its fields down to only the subset used by SixDigitAllele
        """
        if self.__class__ is SixDigitAllele:
            return self
        else:
            return SixDigitAllele(
                self.species_prefix,
                self.gene_name,
                self.group_id,
                self.protein_id,
                self.coding_sequence_id,
                self.modifier)

    def to_record(self):
        """
        Returns dictionary with all fields of this allele,
        as well as its representations as a gene, allele group,
        four digit allele, and six digit allele.
        """
        d = FourDigitAllele.to_record(self)
        d["allele"] = d["six_digit_allele"] = self.normalized_string()
        return d

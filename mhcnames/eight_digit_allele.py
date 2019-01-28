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

from .six_digit_allele import SixDigitAllele


class EightDigitAllele(SixDigitAllele):
    def __init__(
            self,
            species_prefix,
            gene_name,
            group_id,
            protein_id,
            coding_sequence_id,
            genomic_sequence_id,
            modifier=None):
        SixDigitAllele.__init__(
            self,
            species_prefix,
            gene_name,
            group_id,
            protein_id,
            coding_sequence_id,
            modifier)
        self.genomic_sequence_id = genomic_sequence_id

    @classmethod
    def field_names(cls):
        return (
            "species_prefix",
            "gene_name",
            "group_id",
            "protein_id",
            "coding_sequence_id",
            "genomic_sequence_id",
            "modifier"
        )

    def normalized_string(self, include_species=True, include_modifier=True):
        """
        Return allele strings like "HLA-A*02:01"
        """
        result = "%s:%s" % (
            SixDigitAllele.normalized_string(
                self,
                include_species=include_species,
                include_modifier=False),
            self.coding_sequence_id)
        if include_modifier and self.modifier:
            result += self.modifier
        return result

    def compact_string(self, include_species=False):
        """
        Compact representation of an EightDigitAllele, omits the "*" and ":"
        in allele names.
            Normalized: HLA-A*02:01:01:01
            Compact: HLA-A02010101
        """
        return "%s%s" % (
            SixDigitAllele.compact_string(include_species=include_species),
            self.genomic_sequence_id)

    @classmethod
    def from_four_digit_allele(cls, four_digit_allele):
        """
        Invent a 'default' eight-digit allele for a four-digit allele
        by append ":01:01" at the end of the name.
        For example, HLA-A*02:01 becomes HLA-A*02:01:01:01
        """
        return cls(
            four_digit_allele.species_prefix,
            four_digit_allele.gene_name,
            four_digit_allele.group_id,
            protein_id=four_digit_allele.protein_id,
            coding_sequence_id="01",
            genomic_sequence_id="01",
            modifier=four_digit_allele.modifier)

    def from_six_digit_allele(self, six_digit_allele):
        """
        Create a 'default' eight-digit allele for a six-digit allele
        by append ":01" at the end of the name.
        For example, HLA-A*02:01:01 becomes HLA-A*02:01:01:01
        """
        return EightDigitAllele(
            six_digit_allele.species_prefix,
            six_digit_allele.gene_name,
            six_digit_allele.group_id,
            six_digit_allele.protein_id,
            six_digit_allele.coding_sequence_id,
            genomic_sequence_id="01",
            modifier=self.modifier)

    def to_record(self):
        """
        Returns dictionary with all fields of this allele,
        as well as its representations as a gene, allele group,
        four digit allele, six digit allele, eight_digit_allele.
        """
        d = SixDigitAllele.to_record(self)
        d["allele"] = d["eight_digit_allele"] = self.normalized_string()
        return d

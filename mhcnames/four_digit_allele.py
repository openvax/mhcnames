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

from .allele_group import AlleleGroup

class FourDigitAllele(AlleleGroup):
    """
    Allele name which specifies a unique protein amino acid sequence
    using this kind of notation: "HLA-A*02:01" or more generally:
            Species-Gene*Group:ProteinID
    """

    def __init__(self, species_prefix, gene_name, group_id, protein_id, modifier=None):
        AlleleGroup.__init__(self, species_prefix, gene_name, group_id)
        self.protein_id = protein_id
        self.modifier = modifier

    @classmethod
    def from_allele_group(cls, allele_group, protein_id, modifier=None):
        return cls(
            species_prefix=allele_group.species_prefix,
            gene_name=allele_group.gene_name,
            group_id=allele_group.group_id,
            protein_id=protein_id,
            modifier=modifier)

    @classmethod
    def from_tuple(cls, t):
        return cls(
            species_prefix=t[0],
            gene_name=t[1],
            group_id=t[2],
            protein_id=t[3],
            modifier=t[6])

    def to_tuple(self):
        return (
            self.species_prefix,
            self.gene_name,
            self.group_id,
            self.protein_id,
            self.modifier,
        )

    def normalized_string(self, include_species=True, include_modifier=True):
        """
        Return allele strings like "HLA-A*02:01"
        """
        result = "%s:%s" % (
            AlleleGroup.normalized_string(self, include_species=include_species),
            self.protein_id)
        if include_modifier and self.modifier:
            result += self.modifier
        return result

    def compact_string(self, include_species=True):
        """
        Compact representation of a FourDigitAllele, omits the "*" and ":"
        in allele names
            Normalized: HLA-A*02:01
            Compact: HLA-A0201
        """
        return "%s%s" % (
            AlleleGroup.compact_string(include_species=include_species),
            self.protein_id)

    def to_four_digit_allele(self):
        """
        Descendant classes use this method to project their fields
        down to a FourDigitAllele object.
        """
        if self.__class__ is FourDigitAllele:
            return self
        else:
            return FourDigitAllele(
                self.species_prefix,
                self.gene_name,
                self.group_id,
                self.protein_id,
                modifier=self.modifier)

    def __str__(self):
        return "%s:%s%s" % (
            AlleleGroup.__str__(self),
            self.nonsyn,
            self.modifier if self.modifier else "")

    def to_dict(self):
        """
        Returns dictionary with all fields of this allele,
        as well as its representations as a gene, allele group,
        and four digit allele.
        """
        d = AlleleGroup.to_dict(self)
        d["allele"] = d["four_digit_allele"] = self.normalized_string()
        d["modifier"] = self.modifier
        return d

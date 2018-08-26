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

from serializable import Serializable

class Locus(Serializable):
    def __init__(self, species_prefix, gene_name):
        self.species_prefix = species_prefix
        self.gene_name = gene_name

    def normalized_string(self, include_species=True):
        if include_species:
            return "%s-%s" % (self.species_prefix, self.gene_name)
        else:
            return self.gene_name

    def to_locus(self):
        """
        For Locus objects this acts as a simple copy but descendant classes
        use this method to project their fields down to a Locus object.
        """
        return Locus(self.species_prefix, self.gene_name)

class AlleleGroup(Locus):
    """
    Representation of a group of closely related alleles,
    specified by a species, a gene, and a group identifier,
    such as "HLA-A*02".

    These are not the same as serotypes (e.g. HLA-A2). Serotypes
    correspond to a group of proteins which can all be identified with
    the same antibody. Allele groups in modern MHC nomenclature
    are closely related but not every allele in a group is also
    in the similarly named serotype (and vice versa).
    """
    def __init__(self, species_prefix, gene_name, group):
        Locus.__init__(self, species_prefix, gene_name)
        self.group = group

    def normalized_string(self, include_species=True):
        return "%s*%s" % (
            Locus.normalized_string(self, include_species=include_species),
            self.group)

    def to_allele_group(self):
        """
        For AlleleGroup objects this acts as a simple copy but descendant
        classes use this method to project their fields down to a AlleleGroup
        object.
        """
        return AlleleGroup(self.species_prefix, self.gene_name, self.group)

class FourDigitAllele(AlleleGroup):
    def __init__(self, species_prefix, gene_name, group, protein_id, modifier=None):
        AlleleGroup.__init__(self, species_prefix, gene_name, group)
        self.protein_id = protein_id
        self.modifier = modifier

    _parse_cache = {}

    @classmethod
    def parse_with_extra_characters(cls, name):
        pass

    @classmethod
    def parse(cls, name, default_species_prefix=None):
        if name in cls._parse_cache:
            return cls._parse_cache[name]

        allele_group, extra_characters = AlleleGroup.parse_with_extra_characters(
            name,
            default_species_prefix=default_species_prefix)

        result = FourDigitAllele(
            species_prefix=allele_group.species_prefix,
            gene_name=allele_group.gene_name,
            group=allele_group.group,
            protein_id=protein_id,
            modifier=modifier)

        cls._parse_cache[name] = result
        return result

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

    def to_four_digit_allele(self):
        """
        On FourDigitAllele objects this acts as a simple copy but descendant
        classes use this method to project their fields down to a FourDigitAllele
        object.
        """
        return FourDigitAllele(
            self.species_prefix,
            self.gene_name,
            self.group,
            self.protein_id,
            modifier=self.modifier)

    def to_six_digit_allele(self):
        """
        Invent a 'default' six-digit allele for a four-digit allele
        by append ":01" at the end of the name.
        For example, HLA-A*02:01 becomes HLA-A*02:01:01
        """
        return SixDigitAllele(
            self.species_prefix,
            self.gene_name,
            self.group,
            protein_id=self.protein_id,
            coding_sequence_id="01",
            modifier=self.modifier)

    def to_eight_digit_allele(self):
        """
        Invent a 'default' eight-digit allele for a four-digit allele
        by append ":01:01" at the end of the name.
        For example, HLA-A*02:01 becomes HLA-A*02:01:01:01
        """
        return EightDigitAllele(
            self.species_prefix,
            self.gene_name,
            self.group,
            protein_id=self.protein_id,
            coding_sequence_id="01",
            genomic_sequence_id="01",
            modifier=self.modifier)

    def __str__(self):
        return "%s:%s%s" % (
            AlleleGroup.__str__(self),
            self.nonsyn,
            self.modifier if self.modifier else "")

class SixDigitAllele(FourDigitAllele):
    def __init__(
            self,
            species_prefix,
            gene_name,
            group,
            protein_id,
            coding_sequence_id,
            modifier=None):
        FourDigitAllele.__init__(
            self,
            species_prefix,
            gene_name,
            group,
            protein_id,
            modifier)
        self.coding_sequence_id = coding_sequence_id

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

    def to_six_digit_allele(self):
        """
        Acts as a copy function for SixDigitAllele functions
        but the descendant class EightDigitAllele can use this
        to project its fields down to only the subset used by SixDigitAllele
        """
        return SixDigitAllele(
            self.species_prefix,
            self.gene_name,
            self.group,
            self.protein_id,
            self.coding_sequence_id,
            self.modifier)

    def to_eight_digit_allele(self):
        """
        Create a 'default' eight-digit allele for a six-digit allele
        by append ":01" at the end of the name.
        For example, HLA-A*02:01:01 becomes HLA-A*02:01:01:01
        """
        return EightDigitAllele(
            self.species_prefix,
            self.gene_name,
            self.group,
            self.protein_id,
            self.coding_sequence_id,
            genomic_sequence_id="01",
            modifier=self.modifier)

class EightDigitAllele(SixDigitAllele):
    def __init__(
            self,
            species_prefix,
            gene_name,
            group,
            protein_id,
            coding_sequence_id,
            genomic_sequence_id,
            modifier=None):
        SixDigitAllele.__init__(
            self,
            species_prefix,
            gene_name,
            group,
            protein_id,
            coding_sequence_id,
            modifier)
        self.genomic_sequence_id = genomic_sequence_id

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

class AlphaBetaPair(Serializable):
    def __init__(self, alpha, beta):
        self.alpha = alpha
        self.beta = beta

    def normalized_string(self):
        return "%s/%s" % (
            self.alpha.normalized_string(include_species=True),
            self.beta.normalized_string(include_species=False))

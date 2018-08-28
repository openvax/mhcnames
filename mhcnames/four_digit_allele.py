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

import re

from .allele_parse_error import AlleleParseError

from .locus import AlleleGroup
from .allele_modifiers import check_for_allele_modifier

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
    def parse_substring(cls, name):
        """
        Parses a sequence like "HLA-A*02:01:02" and returns the corresponding
        FourDigitAllele object for "HLA-A*02:01" and the remaining
        characters ":02"
        """
        if name in cls._parse_cache:
            return cls._parse_cache[name]
        name = name.strip()
        allele_group, extra_characters = AlleleGroup.parse_normalized_substring(name)
        extra_characters = extra_characters.strip()
        regex = "^:?([A-Za-z0-9]+)"
        match_obj = re.match(regex, extra_characters)
        if match_obj:
            raise AlleleParseError("Unable parse '%s' beyond allele group %s" % (
                name, allele_group))
        _, end_index = match_obj.span()
        protein_id = extra_characters[:end_index]
        new_extra_characters = extra_characters[end_index:].strip()
        modifier, new_extra_characters = check_for_allele_modifier(new_extra_characters)
        four_digit_allele = FourDigitAllele.from_allele_group(
            allele_group=allele_group,
            protein_id=protein_id,
            modifier=modifier)
        cls._parse_cache[name] = (four_digit_allele, new_extra_characters)
        return four_digit_allele, new_extra_characters

    @classmethod
    def parse(cls, name):
        """
        Parse four-digit allele from normalized representation such
        as HLA-A*02:01. Does not handle variability arising from
        retired alleles, gene aliases, capitalization or non-standard
        separators.
        """
        four_digit_allele, extra_characters = cls.parse_substring(name)
        if len(extra_characters) > 0:
            raise AlleleParseError(
                "Unable to parse '%s', found extra characters '%s' after %s" % (
                    name,
                    extra_characters,
                    four_digit_allele))
        return four_digit_allele

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
        as well as its representations as a locus, allele group,
        and four digit allele.
        """
        d = AlleleGroup.to_dict(self)
        d["four_digit_allele"] = self.normalized_string()
        d["protein_id"] = self.protein_id
        d["modifier"] = self.modifier
        return d

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

from typing import Union

from pytypes import typechecked

from .gene import Gene
from .parsed_result import ParsedResult


class AlleleGroup(ParsedResult):
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
    @typechecked
    def __init__(self, gene : Gene, group_id : str):
        self.gene = gene
        self.group_id = self._adjust_formatting(group_id)

    @property
    def species(self):
        return self.gene.species

    @property
    def species_prefix(self):
        return self.species.prefix

    @property
    def gene_name(self):
        return self.gene.name

    @classmethod
    def get(cls, species_prefix, gene_name, group_id):
        gene = Gene.get(species_prefix, gene_name)
        if gene is None:
            return None
        return AlleleGroup(gene, group_id)

    def _adjust_formatting(self, group_id):
        if self.species_prefix in {"Mamu", "Mafa", "ELA", "DLA", "Ecqa", "Calu", "BoLA"}:
            if len(group_id) == 2:
                if group_id.isdigit():
                    return "0" + group_id
        return group_id

    @classmethod
    def field_names(cls):
        return (
            "gene",
            "group_id"
        )

    def normalized_string(self, include_species=True, use_species_alias=True):
        return "%s*%s" % (
            self.gene.normalized_string(
                include_species=include_species,
                use_species_alias=use_species_alias),
            self.group_id)

    def compact_string(self, include_species=False, use_species_alias=True):
        """
        Compact representation of an AlleleGroup, omits the "*"
        in an allele group.
            Normalized: HLA-A*02
            Compact: HLA-A02
        """
        gene_string = self.gene.compact_string(
            include_species=include_species,
            use_species_alias=use_species_alias)
        requires_sep = (
            (gene_string[-1].isdigit() and self.group_id[0].isdigit()) or
            (gene_string[-1].isalpha() and self.group_id[0].isalpha())
        )
        if requires_sep:
            return "%s*%s" % (
                gene_string,
                self.group_id)
        else:
            return "%s%s" % (
                gene_string,
                self.group_id)

    def to_allele_group(self):
        """
        For AlleleGroup objects this acts as a simple copy but descendant
        classes use this method to project their fields down to a AlleleGroup
        object.
        """
        return AlleleGroup(self.species_prefix, self.gene_name, self.group_id)

    def to_record(self):
        """
        Returns dictionary with all fields of this allele group,
        as well as its representations as a locus.
        """
        d = self.gene.to_record(self)
        d["allele_group"] = self.normalized_string()
        return d


class TwoDigitAllele(ParsedResult):
    """
    A few species have a single numeric field for identifying unique MHC
    proteins, e.g. Anpl-UAA*01
    """
    @typechecked
    def __init__(
            self,
            gene : Gene,
            protein_id : str,
            modifier : Union[None, str]=None):
        self.gene = gene
        self.protein_id = protein_id
        self.modifier = modifier

    @classmethod
    def field_names(cls):
        return (
            "gene",
            "protein_id",
            "modifier"
        )

    @property
    def species(self):
        return self.allele_group.species

    @property
    def species_prefix(self):
        return self.species.prefix

    @property
    def gene_name(self):
        return self.gene.name

    @classmethod
    def get(cls, species_prefix, gene_name, protein_id):
        gene = Gene.get(species_prefix, gene_name)
        if gene is None:
            return None
        return TwoDigitAllele(gene, protein_id)

    def normalized_string(
            self,
            include_species=True,
            use_species_alias=True,
            include_modifier=True):
        """
        Return allele strings like "Anpl-UAA*01"
        """
        allele_group_str = self.gene.normalized_string(
            include_species=include_species,
            use_species_alias=use_species_alias)
        result = "%s*%s" % (allele_group_str, self.protein_id)
        if include_modifier and self.modifier:
            result += self.modifier
        return result

    def compact_string(self, include_species=False, use_species_alias=True):
        """
        Compact representation of a TwoDigitAllele, omits the "*" and ":"
        in allele names
            Normalized: HLA-A*02:01
            Compact: HLA-A0201
        """
        return "%s%s" % (
            self.gene.compact_string(
                include_species=include_species,
                use_species_alias=use_species_alias),
            self.protein_id)

    def to_record(self):
        """
        Returns dictionary with all fields of this allele,
        as well as its representations as a gene, allele group,
        and four digit allele.
        """
        d = self.gene.to_record()
        d["allele"] = self.normalized_string()
        d["modifier"] = self.modifier
        d["is_mutant"] = False
        return d

class FourDigitAllele(ParsedResult):
    """
    Allele name which specifies a unique protein amino acid sequence
    using this kind of notation: "HLA-A*02:01" or more generally:
            Species-Gene*Group:ProteinID
    """
    @typechecked
    def __init__(
            self,
            allele_group : AlleleGroup,
            protein_id : str,
            modifier : Union[None, str] = None):
        self.allele_group = allele_group
        self.protein_id = protein_id
        self.modifier = modifier

    @classmethod
    def field_names(cls):
        return (
            "allele_group",
            "protein_id",
            "modifier"
        )

    @property
    def species(self):
        return self.allele_group.species

    @property
    def species_prefix(self):
        return self.species.prefix

    @property
    def gene(self):
        return self.allele_group.gene

    @property
    def gene_name(self):
        return self.allele_group.gene_name

    @property
    def group_id(self):
        return self.allele_group.group_id

    @property
    def mhc_class(self):
        return self.gene.mhc_class

    @property
    def is_class1(self):
        return self.gene.is_class1

    @property
    def is_class2(self):
        return self.gene.is_class2

    @classmethod
    def get(
            cls,
            species_prefix,
            gene_name,
            group_id,
            protein_id,
            modifier=None):
        allele_group = AlleleGroup.get(species_prefix, gene_name, group_id)
        if allele_group is None:
            return None
        return FourDigitAllele(allele_group, protein_id, modifier=modifier)

    def normalized_string(
            self,
            include_species=True,
            use_species_alias=True,
            include_modifier=True):
        """
        Return allele strings like "HLA-A*02:01"
        """
        allele_group_str = self.allele_group.normalized_string(
            include_species=include_species,
            use_species_alias=use_species_alias)
        result = "%s:%s" % (allele_group_str, self.protein_id)
        if include_modifier and self.modifier:
            result += self.modifier
        return result

    def compact_string(self, include_species=False, use_species_alias=True):
        """
        Compact representation of a FourDigitAllele, omits the "*" and ":"
        in allele names
            Normalized: HLA-A*02:01
            Compact: HLA-A0201
        """
        return "%s%s" % (
            self.allele_group.compact_string(
                include_species=include_species,
                use_species_alias=use_species_alias),
            self.protein_id)

    def to_record(self):
        """
        Returns dictionary with all fields of this allele,
        as well as its representations as a gene, allele group,
        and four digit allele.
        """
        d = self.allele_group.to_record()
        d["allele"] = d["four_digit_allele"] = self.normalized_string()
        d["modifier"] = self.modifier
        d["is_mutant"] = False
        return d


class SixDigitAllele(ParsedResult):

    @typechecked
    def __init__(
            self,
            four_digit_allele : FourDigitAllele,
            coding_sequence_id : str,
            modifier : Union[None, str]=None):
        self.four_digit_allele = four_digit_allele
        self.coding_sequence_id = coding_sequence_id
        self.modifier = modifier


    @classmethod
    def field_names(cls):
        return (
            "four_digit_allele",
            "coding_sequence_id",
            "modifier"
        )

    @property
    def species(self):
        return self.four_digit_allele.species

    @property
    def species_prefix(self):
        return self.species.prefix

    @property
    def gene(self):
        return self.four_digit_allele.gene

    @property
    def gene_name(self):
        return self.gene.name

    @property
    def allele_group(self):
        return self.four_digit_allele.allele_group

    @property
    def group_id(self):
        return self.allele_group.group_id

    @property
    def protein_id(self):
        return self.four_digit_allele.protein_id


    @property
    def mhc_class(self):
        return self.gene.mhc_class

    @property
    def is_class1(self):
        return self.gene.is_class1

    @property
    def is_class2(self):
        return self.gene.is_class2

    @classmethod
    def get(
            cls,
            species_prefix,
            gene_name,
            group_id,
            protein_id,
            coding_sequence_id,
            modifier=None):
        four_digit_allele = FourDigitAllele.get(
            species_prefix, gene_name, group_id, protein_id, modifier=modifier)
        if four_digit_allele is None:
            return None
        return SixDigitAllele(
            four_digit_allele,
            coding_sequence_id=coding_sequence_id,
            modifier=modifier)

    def normalized_string(self, include_species=True, include_modifier=True):
        """
        Return allele strings like "HLA-A*02:01"
        """
        result = "%s:%s" % (
            self.four_digit_allele.normalized_string(
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
            self.four_digit_allele.compact_string(include_species=include_species),
            self.coding_sequence_id)


    def to_record(self):
        """
        Returns dictionary with all fields of this allele,
        as well as its representations as a gene, allele group,
        four digit allele, and six digit allele.
        """
        d = self.four_digit_allele.to_record()
        d["allele"] = d["six_digit_allele"] = self.normalized_string()
        d["modifier"] = self.modifier
        return d


class EightDigitAllele(ParsedResult):

    @typechecked
    def __init__(
            self,
            six_digit_allele : SixDigitAllele,
            genomic_sequence_id : str,
            modifier : Union[None, str]=None):
        self.six_digit_allele = six_digit_allele
        self.genomic_sequence_id = genomic_sequence_id
        self.modifier = modifier

    @classmethod
    def field_names(cls):
        return (
            "six_digit_allele",
            "genomic_sequence_id",
            "modifier"
        )

    @property
    def species(self):
        return self.six_digit_allele.species

    @property
    def species_prefix(self):
        return self.species.prefix

    @property
    def gene(self):
        return self.six_digit_allele.gene

    @property
    def gene_name(self):
        return self.gene.name

    @property
    def allele_group(self):
        return self.six_digit_allele.allele_group

    @property
    def four_digit_allele(self):
        return self.six_digit_allele.four_digit_allele

    @property
    def group_id(self):
        return self.allele_group.group_id

    @property
    def protein_id(self):
        return self.four_digit_allele.protein_id

    @property
    def coding_sequence_id(self):
        return self.six_digit_allele.coding_sequence_id

    @property
    def mhc_class(self):
        return self.gene.mhc_class

    @property
    def is_class1(self):
        return self.gene.is_class1

    @property
    def is_class2(self):
        return self.gene.is_class2

    @classmethod
    def get(
            cls,
            species_prefix,
            gene_name,
            group_id,
            protein_id,
            coding_sequence_id,
            genomic_sequence_id,
            modifier=None):
        six_digit_allele = SixDigitAllele.get(
            species_prefix,
            gene_name,
            group_id,
            protein_id,
            coding_sequence_id,
            modifier=modifier)
        if six_digit_allele is None:
            return None
        return EightDigitAllele(
            six_digit_allele,
            genomic_sequence_id=genomic_sequence_id,
            modifier=modifier)

    def normalized_string(
            self,
            include_species=True,
            use_species_alias=True,
            include_modifier=True):
        """
        Return allele strings like "HLA-A*02:01"
        """
        result = "%s:%s" % (
            self.six_digit_allele.normalized_string(
                self,
                include_species=include_species,
                use_species_alias=use_species_alias,
                include_modifier=False),
            self.genomic_sequence_id)
        if include_modifier and self.modifier:
            result += self.modifier
        return result

    def compact_string(self, include_species=False, use_species_alias=True):
        """
        Compact representation of an EightDigitAllele, omits the "*" and ":"
        in allele names.
            Normalized: HLA-A*02:01:01:01
            Compact: HLA-A02010101
        """
        return "%s%s" % (
            self.six_digit_allele.compact_string(
                include_species=include_species,
                use_species_alias=use_species_alias),
            self.genomic_sequence_id)

    def to_record(self):
        """
        Returns dictionary with all fields of this allele,
        as well as its representations as a gene, allele group,
        four digit allele, six digit allele, eight_digit_allele.
        """
        d = self.six_digit_allele.to_record()
        d["modifier"] = self.modifier
        d["eight_digit_allele"] = self.normalized_string()
        d["allele"] = d["eight_digit_allele"]
        return d

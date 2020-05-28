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

import re

from .alpha_beta_pair import AlphaBetaPair, infer_class2_alpha_chain
from .allele_parse_error import AlleleParseError
from .parsing_helpers import (
    strip_whitespace_and_trim_outer_quotes,
    strip_whitespace_and_dashes,
    split_on_all_seps,
)
from .data import haplotypes
from .mutation import Mutation
from .four_digit_allele import FourDigitAllele
from .six_digit_allele import SixDigitAllele
from .eight_digit_allele import EightDigitAllele
from .allele_group import AlleleGroup
from .gene import Gene
from .mhc_class import MhcClass
from .species import infer_species_prefix_substring, find_matching_species
from .serotype_data import get_serotype
from .mutant_allele import MutantAllele
from .serotype import Serotype
from .haplotype import Haplotype
from .allele_modifiers import valid_allele_modifiers
from .mhc_class_helpers import normalize_mhc_class_string
from .named_allele import NamedAllele
from .species import Species


compact_gene_and_allele_regex = re.compile("([A-Za-z]+)([0-9\:]+)[A-Z]?")

class Parser(object):
    def __init__(
            self,
            default_species_prefix="HLA",
            use_species_alias=True,
            infer_class2_pairing=False,
            gene_seps="*_-"):
        """
        use_species_alias : bool
            For species which have a newer four-digit code and and older locus
            name (such as "Ecqa" / "ELA"), use the older species prefix in the
            result.

        infer_class2_pairing : bool
            If only alpha or beta chain of Class II MHC is given, try
            to infer the missing pair?

        default_species_prefix : str
            If no species prefix is given, which should should be assumed?

        gene_seps : iterable of str
            Possible separators used after gene names
        """
        self.default_species_prefix = default_species_prefix
        self.use_species_alias = use_species_alias
        self.infer_class2_pairing = infer_class2_pairing
        self.gene_seps = gene_seps

        self._parse_cache = {}
        self._haplotype_cache = {}
        self._serotype_cache = {}


    def parse_species_prefix(self, name):
        """
        Returns tuple with two elements:
            - species prefix string
            - remaining string after species prefix
        """
        inferred_prefix_and_original = infer_species_prefix_substring(
            name,
            use_species_alias=self.use_species_alias)
        if inferred_prefix_and_original is None:
            if self.default_species_prefix is None:
                return (None, name)
            else:
                return (self.default_species_prefix, name)
        else:
            species_prefix, original_prefix = inferred_prefix_and_original
            original_prefix_length = len(original_prefix)
            remaining_string = name[original_prefix_length:]
            return species_prefix, remaining_string


    def get_species_prefix_and_info(self, name):
        """
        Returns tuple with elements:
            - Species
            - species prefix
            - remaining string after species prefix
        """
        (species_prefix, remaining_string) = self.parse_species_prefix(name)

        if species_prefix is None:
            raise AlleleParseError("Unable to infer species for '%s'" % name)

        species = find_matching_species(species_prefix)
        if species is None:
            raise AlleleParseError("Unknown species '%s' in '%s'" % (species_prefix, name))
        return species, species_prefix, remaining_string


    def parse_serotype(self, species_prefix, serotype_name):
        """
        Returns Serotype or None
        """
        key = (species_prefix, serotype_name)
        if key in self._serotype_cache:
            return self._serotype_cache[key]
        t = get_serotype(species_prefix, serotype_name)
        if t is None:
            result = None
        else:
            species_prefix, serotype_name, allele_list = t
            parsed_allele_objects = []
            for allele in allele_list:
                parsed_allele_objects.append(self.parse(allele))
            result = Serotype(species_prefix, serotype_name, parsed_allele_objects)
        self._serotype_cache[key] = result
        return result

    def parse_haplotype(self, species_prefix, haplotype_name):
        """
        Returns Haplotype or None
        """
        key = (species_prefix, haplotype_name)
        if key in self._haplotype_cache:
            return self._haplotype_cache[key]
        haplotype_entries = haplotypes.get(species_prefix, {}).get(haplotype_name)
        if haplotype_entries is None:
            result = None
        else:
            alleles = [
                self.parse("%s-%s" % (species_prefix, allele))
                for allele in haplotype_entries
            ]
            result = Haplotype(species_prefix, haplotype_name, alleles)
        self._haplotype_cache[key] = result
        return result


    def parse_allele_after_species_and_gene_name(
            self,
            original_name,
            species_prefix,
            gene_name,
            str_after_gene,
            allow_three_digits_in_first_field=False,
            allow_three_digits_in_second_field=False):

        if str_after_gene[-1].upper() in valid_allele_modifiers:
            modifier = str_after_gene[-1]
            str_after_gene = str_after_gene[:-1]
        else:
            modifier = None

        parts = split_on_all_seps(str_after_gene)

        parsed_fields = []
        for part in parts:
            if part.isdigit():
                if (allow_three_digits_in_first_field and
                        len(parsed_fields) == 0 and
                        len(part) > 4):
                    parsed_fields.append(part[:3])
                    part = part[3:]
                if (allow_three_digits_in_second_field and
                        len(parsed_fields) == 1 and
                        len(part) > 4):
                    parsed_fields.append(part[3:])
                while part:
                    n_parsed = len(parsed_fields)
                    remaining_length = len(part)
                    if remaining_length == 1:
                        raise AlleleParseError("Unable to parse '%s'" % original_name)
                    if (allow_three_digits_in_first_field and n_parsed == 0 and
                            (remaining_length == 3 or remaining_length > 4)):
                        boundary_index = 3
                    elif (allow_three_digits_in_second_field and n_parsed == 1 and
                            (remaining_length == 3 or remaining_length > 4)):
                        boundary_index = 3
                    else:
                        boundary_index = 2
                    parsed_fields.append(part[:boundary_index])
                    part = part[boundary_index:]
            else:
                parsed_fields.append(part)

        if len(parsed_fields) in {0, 1} and modifier is not None:
                raise AlleleParseError("Unexpected modifier '%s' at end of '%s'" % (
                    modifier, original_name))

        if len(parsed_fields) == 0:
            return Gene(species_prefix, gene_name)
        elif len(parsed_fields) == 1:
            return AlleleGroup(
                species_prefix,
                gene_name,
                parsed_fields[0])
        elif len(parsed_fields) == 2:
            return FourDigitAllele(
                species_prefix,
                gene_name,
                parsed_fields[0],
                parsed_fields[1],
                modifier=modifier)
        elif len(parsed_fields) == 3:
            return SixDigitAllele(
                species_prefix,
                gene_name,
                parsed_fields[0],
                parsed_fields[1],
                parsed_fields[2],
                modifier=modifier)
        elif len(parsed_fields) == 4:
            return EightDigitAllele(
                species_prefix,
                gene_name,
                parsed_fields[0],
                parsed_fields[1],
                parsed_fields[2],
                parsed_fields[3],
                modifier=modifier)




    def parse_gene_if_possible(self, species, name):
        """
        Parse gene such as "A" or "DQB" and return it along with
        remaining string.

        Return None if not possible.
        """
        for n in range(len(name), 0, -1):
            substring = name[:n]

            gene_name = species.find_matching_gene_name(substring)

            if gene_name:
                return gene_name, name[n:]
        return None, name


    def parse_gene_name_from_prefix(
            self,
            species,
            original_name,
            str_after_species):
        gene_name, str_after_gene = self.parse_gene_if_possible(species, str_after_species)


        if gene_name is None:
            for sep in self.gene_seps:
                if str_after_species.count(sep) == 1:
                    gene_name, str_after_gene = str_after_species.split(sep)
                    break
        else:
            for sep in self.gene_seps:
                if str_after_gene.startswith(sep):
                    str_after_gene = str_after_gene[1:]

        if gene_name is None:
            # If the string had neither "*" nor "_" then try to collect the gene
            # name as the non-numerical part at the start of the string.
            match = compact_gene_and_allele_regex.fullmatch(str_after_species)
            if match:
                gene_name, str_after_gene = match.groups()

        if gene_name is None:
            # failed to figure out a gene name, give up
            raise AlleleParseError("Unable to parse '%s'" % original_name)

        return gene_name, str_after_gene


    def parse_murine_gene(self, species, gene_name, str_after_gene):
        return NamedAllele(species.prefix, gene_name, str_after_gene)


    def split_by_hyphen_except_gene_names(self, species, str_after_species):
        """
        Split a string into a list of parts by hyphens except keep
        gene names such as "M3-1" together
        """
        parts = str_after_species.split("-")
        parts_with_merged_gene_names = []
        i = 0
        while i < len(parts):
            first_part = parts[i]
            if i + 1 == len(parts):
                parts_with_merged_gene_names.append(first_part)
                break

            next_part = parts[i + 1]
            combined = "%s-%s" % (first_part, next_part)
            if species.find_matching_gene_name(combined):
                parts_with_merged_gene_names.append(combined)
                i += 2
            else:
                parts_with_merged_gene_names.append(first_part)
                i += 1
        return parts_with_merged_gene_names


    def parse_without_mutation(self, name):
        """
        First test to see if MHC name requires any species-specific special logic.
        If the name doesn't fit any of the special species templates then
        try parsing it as the following kinds of names in this order:
            1) eight digit allele
            2) six digit allele
            3) four digit allele
            4) allele group
            5) locus
            6) MHC class
            7) species
        If none of these succeed, then raise an exception
        """
        species, species_prefix, str_after_species = \
            self.get_species_prefix_and_info(name)

        str_after_species = strip_whitespace_and_dashes(str_after_species)

        if len(str_after_species) == 0:
            return species

        try:
            # try to heuristically split apart the gene name and any allele information
            # when the requires separators are missing
            # Examples which will parse correctly here:
            #   A*0201
            #   A*02:01
            #   A_0101
            #   A_01:01
            #   A-0101
            #   A-01:01
            # However this will not work:
            #   - A_01_01
            gene_name, str_after_gene = self.parse_gene_name_from_prefix(
                species,
                name,
                str_after_species)
            # use the canonical gene name e.g. "A" and not "a"
            gene_name = species.normalize_gene_name_if_exists(gene_name)

            str_after_gene = strip_whitespace_and_dashes(str_after_gene)

            if len(str_after_gene) == 0:
                return Gene(species.prefix, gene_name)

            # if the remaining string is an allele string which has
            # been renamed or deprecated, then get its new/canonical form
            new_allele_name = species.allele_aliases.get("%s*%s" % (gene_name, str_after_gene))
            if new_allele_name:
                gene_name, str_after_gene = new_allele_name.split("*")

            if species.prefix in {"H2", "RT1"}:
                return self.parse_murine_gene(species, gene_name, str_after_gene)

            # only allele names which allow three digits in second field seem to be
            # human class I names such as "HLA-B*15:120",
            # it's otherwise typical to allow three digits in the first field
            allow_three_digits_in_second_field = (
                    species_prefix == "HLA" and gene_name in {"A", "B", "C"}
            )
            allow_three_digits_in_first_field = not allow_three_digits_in_second_field

            return self.parse_allele_after_species_and_gene_name(
                original_name=name,
                species_prefix=species_prefix,
                gene_name=gene_name,
                str_after_gene=str_after_gene,
                allow_three_digits_in_first_field=allow_three_digits_in_first_field,
                allow_three_digits_in_second_field=allow_three_digits_in_second_field)
        except AlleleParseError:
            hyphen_parts = self.split_by_hyphen_except_gene_names(
                species,
                str_after_species)

            if len(hyphen_parts) == 1:
                str_after_species = hyphen_parts[0]
            elif len(hyphen_parts) > 2:
                # if after collapsing gene names we still have a hyphen
                # then we expect it to be a Class II alpha-beta pair
                raise AlleleParseError("Unexpected number of '-' in '%s'" % name)
            elif len(hyphen_parts) == 2:
                # this situation is tricky since it might be either
                # a class II allele pair
                #   e.g. DRA1*01:01-DRB1*01:01
                # or a class I allele where '-' is used instead of '*'
                #   e.g. 1-HB01 (swine allele)
                str_after_species = "-".join(hyphen_parts)
                if species.find_matching_gene_name(hyphen_parts[0]) is None:
                    return self.parse_known_alpha_beta_pair(str_after_species)

            serotype_result = self.parse_serotype(species_prefix, str_after_species)
            if serotype_result is not None:
                return serotype_result

            haplotype_result = self.parse_haplotype(species_prefix, str_after_species)
            if haplotype_result:
                return haplotype_result
            else:
                raise AlleleParseError("Unable to parse '%s'" % name)




    def parse_known_alpha_beta_pair(self, name):
        """
        If a name is known to contain "/" then it's
        expected to be of a format like:
            HLA-DQA*01:01/DQB*01:02

        The species information from the first allele
        is used to guide parsing for the second allele.
        """
        if "/" in name:
            parts = name.split("/")
        else:
            parts = name.split("-")

        if len(parts) == 3:
            default_species_prefix, alpha_string, beta_string = parts
        elif len(parts) == 2:
            alpha_string, beta_string = parts
        else:
            raise AlleleParseError(
                "Expected Class II alpha/beta pairing but got %d allele names in '%s'" % (
                    len(parts),
                    name))
        alpha = self.parse(alpha_string)
        beta = self.parse(beta_string)
        return AlphaBetaPair(alpha, beta)


    def parse_with_mutations(self, name):
        """
        Parameters
        ----------
        name : str

        Returns MutantFourDigitAllele
        """
        parts = name.split()
        result_without_mutation = self.parse_without_mutation(parts[0])

        valid_classes = (FourDigitAllele, NamedAllele)
        if result_without_mutation.__class__ not in valid_classes:
            raise AlleleParseError(
                "Cannot apply mutations in '%s' to %s" % (
                    name,
                    result_without_mutation.__class__.__name__))
        # expect names with spaces to be like "A*02:07 T80M mutant"
        # trim off final commas in case we encounter a list of
        # mutations like: "E152A, R155Y, L156Y mutant"
        mutation_strings = [
            p[:-1] if p.endswith(",") else p
            for p in parts[1:]
            if p.upper() != "MUTANT" and p != ""
        ]
        if len(mutation_strings) == 0:
            raise AlleleParseError(
                "Expected '%s' to have mutations but none found" % name)
        mutations = [
            Mutation.parse(mutation_string)
            for mutation_string in mutation_strings
        ]
        return MutantAllele(result_without_mutation, mutations)


    def parse_with_interior_whitespace(self, name):
        """
        If there's whitespace within an allele description then it's
        either a mutant allele or an error.
        """
        lower = name.lower()
        if "mutant" in lower:
            return self.parse_with_mutations(name)
        parts = lower.split()
        if len(parts) == 2:
            # TODO: parse MHC-Id genes and alleles such as "human CD1a"
            raise AlleleParseError("Gene parsing not yet implemented for '%s'" % name)
        elif len(parts) >= 3:
            if parts[-2] == "class" and parts[-1] in {"1", "2", "i", "ii"}:
                mhc_class_string = normalize_mhc_class_string(parts[-1])
                # Parse MHC classes, haplotypes, or serotypes such as:
                # - "HLA class I"
                # - "H2-b class I"
                # - "ELA-A1 class I"
                # - "H2-r class I"
                # - "BF19 class II"
                unrestricted_string = " ".join(parts[:-2])
                unrestricted_result = self.parse(unrestricted_string)
                if unrestricted_result.__class__ is Haplotype:
                    return unrestricted_result.restrict_mhc_class(mhc_class_string)
                elif unrestricted_result.__class__ is Species:
                    return MhcClass(unrestricted_result.species_prefix, mhc_class_string)
                else:
                    raise AlleleParseError(
                        "Unable to parse '%s' in '%s'" % (
                            unrestricted_string,
                            name))

        raise AlleleParseError("Unexpected whitespace in '%s'" % name)

    def parse(self, name):
        """
        Parse any MHC related string, from gene loci to fully specified 8 digit
        alleles, alpha/beta pairings of Class II MHCs, with expression modifiers
        and the description of point mutations in the molecule.

        Example of the complicated inputs this function can handle:
            HLA-DRA*01:02/DRB1*03:01 Q74R mutant
            "H2-Kb E152A, R155Y, L156Y mutant"
            SLA-1*01:01:01:01
            HLA-DRA*01:01 F54C mutant/DRB1*01:01

        Parameters
        ----------
        name : str
            Raw name of MHC locus or allele

        Returns object with one of the following types:
            - Species
            - MhcClass
            - Gene
            - AlleleGroup
            - FourDigitAllele
            - SixDigitAllele
            - EightDigitAllele
            - MutantFourDigitAllele
            - AlphaBetaPair
            - NamedAllele
            - MutantNamedAllele
        """

        if name in self._parse_cache:
            return self._parse_cache[name]

        trimmed_name = strip_whitespace_and_trim_outer_quotes(name)
        if len(trimmed_name) == 0:
            raise AlleleParseError(
                "Cannot parse empty allele name '%s'" % name)
        if " " in trimmed_name or "\t" in trimmed_name:
            result = self.parse_with_interior_whitespace(trimmed_name)
        elif "/" in trimmed_name:
            # parse paired Class II alleles such as 'DRA1*01:01/DRB1*01:01'
            result = self.parse_known_alpha_beta_pair(trimmed_name)
        else:
            result = self.parse_without_mutation(trimmed_name)

        if self.infer_class2_pairing and result.__class__ is not AlphaBetaPair:
            result = infer_class2_alpha_chain(result)

        self._parse_cache[name] = result
        return result

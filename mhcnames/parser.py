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

from .allele_modifiers import valid_allele_modifiers
from .parse_error import ParseError
from .alpha_beta_pair import AlphaBetaPair, infer_class2_alpha_chain
from .data import haplotypes

from .gene import Gene
from .haplotype import Haplotype
from .mhc_class import MhcClass
from .mhc_class_helpers import normalize_mhc_class_string
from .mutant_allele import MutantAllele
from .mutation import Mutation
from .named_allele import NamedAllele
from .numeric_alleles import (
    TwoDigitAllele,
    AlleleGroup,
    FourDigitAllele,
    SixDigitAllele,
    EightDigitAllele
)
from .parsing_helpers import (
    strip_whitespace_and_trim_outer_quotes,
    strip_whitespace_and_dashes,
    split_allele_fields,
    contains_any_letters
)
from .serotype import Serotype
from .serotype_data import get_serotype
from .species import Species, infer_species_prefix_substring, find_matching_species


# default values for Parser parameters, reused in the 'parse' function below
DEFAULT_SPECIES_PREFIX = "HLA"
NORMALIZE_ALLELE_ALIASES = False
INFER_CLASS2_PAIRING = False
GENE_SEPS = "*_-"

class Parser(object):
    def __init__(
            self,
            default_species_prefix=DEFAULT_SPECIES_PREFIX,
            normalize_allele_aliases=NORMALIZE_ALLELE_ALIASES,
            gene_seps=GENE_SEPS):
        """
        default_species_prefix : str
            In the absence of a species prefix, should we assume we're parsing
            a human ("HLA") allele. Set to None to not have a default species.

        normalize_allele_aliases : bool
            Convert old allele aliases to newer names. For example,
            change "SLA-2*07we01" to "SLA-2*07:03"

        gene_seps : iterable of str
            Possible separators used after gene names
        """
        self.default_species_prefix = default_species_prefix
        self.normalize_allele_aliases = normalize_allele_aliases
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
        inferred_prefix_and_original = infer_species_prefix_substring(name)
        if inferred_prefix_and_original is None:
            if self.default_species_prefix is None:
                return (None, name)
            else:
                return (self.default_species_prefix, name)
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
            raise ParseError("Unable to infer species for '%s'" % name)

        species = Species.get(species_prefix)
        if species is None:
            raise ParseError("Unknown species '%s' in '%s'" % (species_prefix, name))
        return species, species_prefix, remaining_string


    def create_serotype_if_exists(self, species, serotype_name):
        """
        Returns Serotype or None
        """
        key = (species.prefix, serotype_name)
        if key in self._serotype_cache:
            return self._serotype_cache[key]
        serotype_tuple = get_serotype(species.prefix, serotype_name)
        if serotype_tuple is None:
            result = None
        else:
            _, serotype_name, allele_list = serotype_tuple
            parsed_allele_objects = []
            for allele in allele_list:
                parsed_allele_objects.append(self.parse(allele))
            result = Serotype(species, serotype_name, parsed_allele_objects)
        self._serotype_cache[key] = result
        return result

    def create_haplotype_if_exists(self, species, haplotype_name):
        """
        Returns Haplotype or None
        """
        key = (species.prefix, haplotype_name)
        if key in self._haplotype_cache:
            return self._haplotype_cache[key]
        haplotype_entries = haplotypes.get(species.prefix, {}).get(haplotype_name)
        if haplotype_entries is None:
            result = None
        else:
            alleles = [
                self.parse("%s-%s" % (species.prefix, allele))
                for allele in haplotype_entries
            ]
            result = Haplotype(species, haplotype_name, alleles)
        self._haplotype_cache[key] = result
        return result


    def parse_allele_after_species_and_gene_name(
            self,
            original_name,
            species,
            gene,
            str_after_gene,
            allow_three_digits_in_first_field=False,
            allow_three_digits_in_second_field=False):

        parsed_fields = split_allele_fields(
            original_name,
            str_after_gene,
            allow_three_digits_in_first_field,
            allow_three_digits_in_second_field)

        if len(parsed_fields) == 0:
            return gene

        if len(parsed_fields) == 1 and not parsed_fields[0].isdigit():
            return NamedAllele(gene, parsed_fields[0].lower())


        # TODO: distinguish between an AlleleGroup and TwoDigitAllele
        allele_group = AlleleGroup(gene, parsed_fields[0])

        if len(parsed_fields) == 1:
            return allele_group

        modifier = None
        for i, field in enumerate(parsed_fields):
            if i == len(parsed_fields) - 1 and field[-1].upper() in valid_allele_modifiers:
                modifier = field[-1].upper()
                field = field[:-1]

            if not field.isdigit():
                raise ParseError(
                    "Expected all fields to be digits in '%s'" % (
                        original_name,))

        four_digit_allele = FourDigitAllele(
            allele_group,
            protein_id=parsed_fields[1],
            modifier=modifier)

        if len(parsed_fields) == 2:
            return four_digit_allele

        six_digit_allele = SixDigitAllele(
            four_digit_allele,
            coding_sequence_id=parsed_fields[2],
            modifier=modifier)

        if len(parsed_fields) == 3:
            return six_digit_allele

        if len(parsed_fields) == 4:
            return EightDigitAllele(
                    six_digit_allele=six_digit_allele,
                    genomic_sequence_id=parsed_fields[3],
                    modifier=modifier)

        raise ParseError("Too many allele fields in '%s'" % original_name)

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

    compact_gene_and_allele_regex = re.compile("([A-Za-z]+)([0-9\:]+)[A-Z]?")

    def parse_gene_from_prefix(
            self,
            species,
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
            match = Parser.compact_gene_and_allele_regex.fullmatch(str_after_species)
            if match:
                gene_name, str_after_gene = match.groups()

        if gene_name is None:
            return None, str_after_species
        gene = Gene.get(species, gene_name)
        return gene, str_after_gene

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


    def parse_serotype_or_haplotype(self, species, str_after_species, original_name):
        hyphen_parts = self.split_by_hyphen_except_gene_names(
            species,
            str_after_species)

        if len(hyphen_parts) == 1:
            str_after_species = hyphen_parts[0]
        elif len(hyphen_parts) > 2:
            # if after collapsing gene names we still have a hyphen
            # then we expect it to be a Class II alpha-beta pair
            return None
        elif len(hyphen_parts) == 2:
            # this situation is tricky since it might be either
            # a class II allele pair
            #   e.g. DRA1*01:01-DRB1*01:01
            # or a class I allele where '-' is used instead of '*'
            #   e.g. 1-HB01 (swine allele)
            str_after_species = "-".join(hyphen_parts)
            if species.find_matching_gene_name(hyphen_parts[0]) is None:
                return self.parse_known_alpha_beta_pair(str_after_species)

        serotype_result = self.create_serotype_if_exists(species, str_after_species)
        if serotype_result is not None:
            return serotype_result

        haplotype_result = self.create_haplotype_if_exists(species, str_after_species)
        if haplotype_result:
            return haplotype_result

        return None

    def parse_without_mutation(self, name, original_name=None):
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
        if original_name is None:
            original_name = name

        species, species_prefix, str_after_species = \
            self.get_species_prefix_and_info(name)

        str_after_species = strip_whitespace_and_dashes(str_after_species)
        if len(str_after_species) == 0:
            return species

        serotype_or_haplotype = self.parse_serotype_or_haplotype(
            species,
            str_after_species,
            original_name=original_name)

        if serotype_or_haplotype:
            return serotype_or_haplotype

        gene, str_after_gene = self.parse_gene_from_prefix(
            species,
            str_after_species)

        if not gene:
            raise ParseError("Unable to parse MHC gene in '%s'" % original_name)

        str_after_gene = strip_whitespace_and_dashes(str_after_gene)

        if len(str_after_gene) == 0:
            return gene

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
        gene_name = gene.name
        if self.normalize_allele_aliases:
            # if the remaining string is an allele string which has
            # been renamed or deprecated, then get its new/canonical form
            # TODO: make this an optional transformation after parsing
            new_allele_name = species.allele_aliases.get(
                "%s*%s" % (gene_name, str_after_gene))
            if new_allele_name:
                gene_name, str_after_gene = new_allele_name.split("*")

        if species.prefix in {"H2", "RT1"}:
            # mouse or rat alleles
            return NamedAllele(gene, str_after_gene.lower())
        elif species.prefix in {"Susc", "SLA"}:
            # parse e.g. "SLA-1-HB03" or "SLA-3-US#11"
            if str_after_gene[:2] == "HB" or "#" in str_after_gene:
                return NamedAllele(
                    gene,
                    str_after_gene.upper())
            elif contains_any_letters(str_after_gene):
                return NamedAllele(
                    gene,
                    str_after_gene.lower())

        # only allele names which allow three digits in second field seem to be
        # human class I names such as "HLA-B*15:120",
        # it's otherwise typical to allow three digits in the first field
        allow_three_digits_in_second_field = (
                species_prefix == "HLA" and gene_name in {"A", "B", "C"}
        )
        allow_three_digits_in_first_field = not allow_three_digits_in_second_field

        return self.parse_allele_after_species_and_gene_name(
            original_name=original_name,
            species=species,
            gene=gene,
            str_after_gene=str_after_gene,
            allow_three_digits_in_first_field=allow_three_digits_in_first_field,
            allow_three_digits_in_second_field=allow_three_digits_in_second_field)


    def parse_known_alpha_beta_pair(self, name, original_name=None):
        """
        If a name is known to contain "/" then it's
        expected to be of a format like:
            HLA-DQA*01:01/DQB*01:02

        The species information from the first allele
        is used to guide parsing for the second allele.
        """
        if original_name is None:
            original_name = name

        if "/" in name:
            parts = name.split("/")
        else:
            parts = name.split("-")

        if len(parts) == 3:
            default_species_prefix, alpha_string, beta_string = parts
        elif len(parts) == 2:
            alpha_string, beta_string = parts
        else:
            raise ParseError(
                "Expected Class II alpha/beta pairing but got %d allele names in '%s'" % (
                    len(parts),
                    original_name))
        alpha = self.parse(alpha_string, infer_class2_pairing=False)
        beta = self.parse(beta_string, infer_class2_pairing=False)
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
            raise ParseError(
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
            raise ParseError(
                "Expected '%s' to have mutations but none found" % name)
        mutations = [
            Mutation.parse(mutation_string)
            for mutation_string in mutation_strings
        ]
        return MutantAllele(result_without_mutation, mutations)


    def parse_with_interior_whitespace(self, name, original_name=None):
        """
        If there's whitespace within an allele description then it's
        either a mutant allele or an error.
        """
        if original_name is None:
            original_name = name
        lower = name.lower()
        if "mutant" in lower:
            return self.parse_with_mutations(name)
        parts = lower.split()
        if len(parts) == 2:
            species_common_name, gene_name = parts
            gene = Gene.get(species_common_name, gene_name)
            if gene:
                return gene
            else:
                raise ParseError("Failed to parse '%s' as gene in '%s'" % (
                    gene_name, original_name))
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
                    return MhcClass(unrestricted_result, mhc_class_string)
                else:
                    raise ParseError(
                        "Unable to parse '%s' in '%s'" % (
                            unrestricted_string,
                            name))

        raise ParseError("Unexpected whitespace in '%s'" % name)

    def parse(self, name, infer_class2_pairing=INFER_CLASS2_PAIRING):
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

        infer_class2_pairing : bool
            If only alpha or beta chain of Class II MHC is given, try
            to infer the missing pair?

        Returns object with one of the following types:
            - Species
            - MhcClass
            - Gene
            - AlleleGroup
            - TwoDigitAllele
            - FourDigitAllele
            - SixDigitAllele
            - EightDigitAllele
            - MutantAllele
            - AlphaBetaPair
            - NamedAllele
        """
        key = (name, infer_class2_pairing)
        if key in self._parse_cache:
            return self._parse_cache[key]

        trimmed_name = strip_whitespace_and_trim_outer_quotes(name)
        if len(trimmed_name) == 0:
            raise ParseError(
                "Cannot parse empty allele name '%s'" % name)
        if "/" in trimmed_name:
            # parse paired Class II alleles such as 'DRA1*01:01/DRB1*01:01'
            result = self.parse_known_alpha_beta_pair(trimmed_name, original_name=name)
        elif " " in trimmed_name or "\t" in trimmed_name:
            result = self.parse_with_interior_whitespace(trimmed_name, original_name=name)
        else:
            result = self.parse_without_mutation(trimmed_name, original_name=name)

        if infer_class2_pairing:
            result = infer_class2_alpha_chain(result)

        self._parse_cache[key] = result
        return result


_parser_cache = {}

def cached_parser(
        default_species_prefix=DEFAULT_SPECIES_PREFIX,
        normalize_allele_aliases=NORMALIZE_ALLELE_ALIASES,
        gene_seps=GENE_SEPS):
    """
    Construct a Parser instance if this combination of arguments hasn't
    been used before, otherwise retrieve an existing parser.
    """
    key = (
        ("normalize_allele_aliases", normalize_allele_aliases),
        ("default_species_prefix", default_species_prefix),
        ("gene_seps", gene_seps),
    )
    if key in _parser_cache:
        return _parser_cache[key]
    else:
        parser = Parser(
            normalize_allele_aliases=normalize_allele_aliases,
            default_species_prefix=default_species_prefix,
            gene_seps=gene_seps)
        _parser_cache[key] = parser
        return parser

def parse(
        raw_string,
        default_species_prefix=DEFAULT_SPECIES_PREFIX,
        normalize_allele_aliases=NORMALIZE_ALLELE_ALIASES,
        infer_class2_pairing=INFER_CLASS2_PAIRING):
    """
    Parse MHC alleles into a structured representation.

    Parameters
    ----------
    raw_string : str
       String corresponding to allele, locus, or other MHC-related name

    default_species_prefix : str
       By default, parse alleles like "A*02:01" as human but it's possible
       to change this to some other species.

    normalize_allele_aliases : bool

    infer_class2_pairing : bool
       If given only the alpha or beta chain of a Class II allele,
       try to infer the most likely pairing from population frequencies.
    """
    parser = cached_parser(
            normalize_allele_aliases=normalize_allele_aliases,
            default_species_prefix=default_species_prefix)
    return parser.parse(raw_string, infer_class2_pairing=infer_class2_pairing)
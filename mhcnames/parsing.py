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

from .alpha_beta_pair import AlphaBetaPair
from .allele_parse_error import AlleleParseError
from .parsing_helpers import (
    strip_whitespace_and_trim_outer_quotes,
    strip_whitespace_and_dashes
)
from .mutation import Mutation
from .four_digit_allele import FourDigitAllele
from .six_digit_allele import SixDigitAllele
from .eight_digit_allele import EightDigitAllele
from .allele_group import AlleleGroup
from .gene import Gene
from .gene_class import GeneClass
from .species_registry import infer_species_prefix_substring, find_matching_species
from .data import (
    allele_aliases_with_uppercase_and_no_dash,
    get_serotype,
)
from .mutant_allele import MutantAllele
from .serotype import Serotype
from .allele_modifiers import valid_allele_modifiers


def parse_species_prefix(name, default_species_prefix=None):
    """
    Returns tuple with two elements:
        - species prefix string
        - remaining string after species prefix
    """
    inferred_prefix_and_original = infer_species_prefix_substring(name)

    if inferred_prefix_and_original is None:
        if default_species_prefix is None:
            return (None, name)
        else:
            return (default_species_prefix, name)
    else:
        species_prefix, original_prefix = inferred_prefix_and_original
        original_prefix_length = len(original_prefix)
        remaining_string = name[original_prefix_length:]
        return species_prefix, remaining_string


def get_species_prefix_and_info(name, default_species_prefix):
    """
    Returns tuple with elements:
        - Species
        - species prefix
        - remaining string after species prefix
    """
    (species_prefix, remaining_string) = \
        parse_species_prefix(name, default_species_prefix=default_species_prefix)
    if species_prefix is None:
        raise AlleleParseError("Unable to infer species for '%s'" % name)

    species = find_matching_species(species_prefix)
    if species is None:
        raise AlleleParseError("Unknown species '%s' in '%s'" % (species_prefix, name))
    return species, species_prefix, remaining_string


def parse_gene_if_possible(name, species_info):
    """
    Parse gene such as "A" or "DQB" and return it along with
    remaining string.

    Return None if not possible.
    """
    for n in range(len(name), 0, -1):
        substring = name[:n]
        gene_name = species_info.find_matching_gene_name(substring)
        if gene_name:
            return gene_name, name[n:]
    return None, name


def normalize_allele_string(species_prefix, allele_sequence_without_species):
    """
    Look up allele name in a species-specific dictionary of aliases
    and, if it's present, substitute allele name with canonical form.
    """
    trimmed = strip_whitespace_and_dashes(allele_sequence_without_species)
    if "*" in trimmed:
        split_by_star = trimmed.split("*")
        if len(split_by_star) > 2:
            raise AlleleParseError(
                "Unexpected number of '*' characters in sequence '%s' " % (
                    allele_sequence_without_species,))
    upper_seq = trimmed.upper()
    species_allele_alias_dict = \
        allele_aliases_with_uppercase_and_no_dash.get(species_prefix, {})
    new_allele_name = species_allele_alias_dict.get(upper_seq)
    if new_allele_name:
        return new_allele_name
    else:
        return trimmed

_serotype_cache = {}


def get_serotype_if_exists(species_prefix, serotype_name):
    key = (species_prefix, serotype_name)
    if key in _serotype_cache:
        return _serotype_cache[key]
    t = get_serotype(species_prefix, serotype_name)
    if t is None:
        result = None
    else:
        species_prefix, serotype_name, allele_list = t
        parsed_allele_objects = []
        for allele in allele_list:
            parsed_allele_objects.append(
                parse(allele, default_species_prefix=species_prefix))
        result = Serotype(species_prefix, serotype_name, parsed_allele_objects)
    _serotype_cache[key] = result
    return result


"""
def normalize_parsed_object(
        species_info,
        parsed_object):
    if species_info is not None:
        if isinstance(parsed_object, Gene):
            old_gene_name = parsed_object.gene_name
            new_gene_name = species_info.normalize_gene_name_if_exists(old_gene_name)

            if old_gene_name != new_gene_name:
                parsed_object = parsed_object.copy(gene_name=new_gene_name)
    return parsed_object
"""


def split_on_all_seps(seq, seps="_:"):
    """
    Split given string on all separators specified

    For example, 02_01:01 will be split into:
        ["02", "01", "01"]
    """
    string_parts = [seq]
    for sep in seps:
        new_parts = []
        for subseq in string_parts:
            new_parts.extend(subseq.split(sep))
        parts = new_parts
    return parts


def parse_allele_after_species_and_gene_name(
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
            raise AlleleParseError("Unexpected part of allele name '%s' in '%s'" % (
                part, original_name))

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


compact_gene_and_allele_regex = re.compile("([A-Za-z]+)([0-9\:]+)[A-Z]?")


def parse_gene_name_from_prefix(original_name, str_after_species, gene_seps="*_"):
    gene_name = None
    for sep in gene_seps:
        if str_after_species.count(sep) == 1:
            gene_name, str_after_gene = str_after_species.split(sep)
            break
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


def parse_without_mutation(
        name,
        default_species_prefix="HLA",
        gene_seps="*_"):
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
    species, species_prefix, str_after_species = get_species_prefix_and_info(
        name,
        default_species_prefix=default_species_prefix)

    str_after_species = normalize_allele_string(
        species_prefix=species_prefix,
        allele_sequence_without_species=str_after_species)

    if len(str_after_species) == 0:
        return species

    serotype_result = get_serotype_if_exists(species_prefix, str_after_species)

    if serotype_result is not None:
        return serotype_result

    # try to heuristically split apart the gene name and any allele information
    # when the requires separators are missing
    # Examples which will parse correctly here:
    #   - A*0201
    #   - A*02:01
    #   - A_0101
    #   - A_01:01
    # However this will not work:
    #   - A_01_01
    # Also, if no gene separator is used (e.g. A0101) then the parsing
    # continues further down.
    gene_name, str_after_gene = parse_gene_name_from_prefix(
        name,
        str_after_species,
        gene_seps=gene_seps)

    # use the canonical gene name e.g. "A" and not "a"
    gene_name = species.normalize_gene_name_if_exists(gene_name)

    # only allele names which allow three digits in second field seem to be
    # human class I names such as "HLA-B*15:120",
    # it's otherwise typical to allow three digits in the first field
    allow_three_digits_in_second_field = (
        species_prefix == "HLA" and gene_name in {"A", "B", "C"}
    )
    allow_three_digits_in_first_field = not allow_three_digits_in_second_field
    return parse_allele_after_species_and_gene_name(
        original_name=name,
        species_prefix=species_prefix,
        gene_name=gene_name,
        str_after_gene=str_after_gene,
        allow_three_digits_in_first_field=allow_three_digits_in_first_field,
        allow_three_digits_in_second_field=allow_three_digits_in_second_field)


def parse_known_alpha_beta_pair(name, default_species_prefix="HLA"):
    """
    If a name is known to contain "/" then it's
    expected to be of a format like:
        HLA-DQA*01:01/DQB*01:02

    The species information from the first allele
    is used to guide parsing for the second allele.
    """
    parts = name.split("/")
    if len(parts) != 2:
        raise AlleleParseError(
            "Expected Class II alpha/beta pairing but got %d allele names in '%s'" % (
                len(parts),
                name))
    alpha_string, beta_string = parts
    alpha = parse(alpha_string, default_species_prefix=default_species_prefix)
    beta = parse(
        beta_string,
        default_species_prefix=alpha.species_prefix)
    return AlphaBetaPair(alpha, beta)


def parse_with_mutations(
        name,
        default_species_prefix):
    """
    Parameters
    ----------
    name : str

    default_species_prefix : str
        If no species can inferred for an allele then use this, e.g. "HLA"

    Returns MutantFourDigitAllele
    """
    parts = name.split()
    result_without_mutation = parse_without_mutation(
        parts[0],
        default_species_prefix=default_species_prefix)
    if result_without_mutation.__class__ is not FourDigitAllele:
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


def parse_with_interior_whitespace(name, default_species_prefix):
    """
    If there's whitespace within an allele description then it's
    either a mutant allele or an error.
    """
    lower = name.lower()
    if "mutant" in lower:
        return parse_with_mutations(
            name,
            default_species_prefix=default_species_prefix)
    parts = lower.split()
    if len(parts) == 2:
        # TODO: parse MHC-Id genes and alleles such as "human CD1a"
        pass
    elif len(parts) == 3:
        if parts[1] == "class" and parts[2] in {"1", "2", "i", "ii"}:
            # parse haplotypes and MHC classes such as:
            # - "HLA class I"
            # - "H2-b class I"
            # - "ELA-A1 class I"
            # - "H2-r class I"
            # - "BF19 class II"
            return GeneClass(parts[0], parts[2])
    raise AlleleParseError("Unexpected whitespace in '%s'" % name)


_parse_cache = {}


def parse(
        name,
        infer_class2_pairing=False,
        default_species_prefix="HLA"):
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

    default_species_prefix : str
        If no species prefix is given, which should should be assumed?

    Returns object with one of the following types:
        - Species
        - GeneClass
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
    cache_key = (name, infer_class2_pairing, default_species_prefix)
    if cache_key in _parse_cache:
        return _parse_cache[cache_key]

    trimmed_name = strip_whitespace_and_trim_outer_quotes(name)
    if len(trimmed_name) == 0:
        raise ValueError(
            "Cannot parse empty allele name '%s'" % name)
    elif "/" in trimmed_name:
        result = parse_known_alpha_beta_pair(
            trimmed_name,
            default_species_prefix=default_species_prefix)
    elif " " in trimmed_name or "\t" in trimmed_name:
        result = parse_with_interior_whitespace(
            trimmed_name,
            default_species_prefix=default_species_prefix)
    else:
        result = parse_without_mutation(
            trimmed_name,
            default_species_prefix=default_species_prefix)

    if infer_class2_pairing and result.__class__ is not AlphaBetaPair:
        if isinstance(result, Gene) and result.is_class2:
            raise AlleleParseError(
                "Inference of paired alpha/beta pair for %s not yet implemented" % (
                    result,))
    _parse_cache[cache_key] = result
    return result

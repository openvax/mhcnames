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


from .human import is_human, parse_human
from .mouse import is_mouse, parse_mouse
from .rat import is_rat, parse_rat
from .swine import is_swine, parse_swine
from .alpha_beta_pair import AlphaBetaPair
from .allele_parse_error import AlleleParseError
from .parsing_helpers import strip_whitespace_and_trim_outer_quotes
from .mutations import Mutation
from .mutant_allele import MutantAllele
from .allele_group import AlleleGroup
from .locus import Locus
from .four_digit_allele import FourDigitAllele
from .six_digit_allele import SixDigitAllele
from .eight_digit_allele import EightDigitAllele

def parse_without_mutation(name, default_species_prefix="HLA"):
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
    if is_mouse(name) or default_species_prefix == "H2":
        return parse_mouse(name)
    elif is_rat(name):
        return parse_rat(name)
    elif is_swine(name):
        return parse_swine(name)
    elif is_human(name):
        return parse_human(name)

    result_classes = [
        EightDigitAllele,
        SixDigitAllele,
        FourDigitAllele,
        AlleleGroup,
        Locus
    ]
    for result_class in result_classes:
        try:
            return result_class.parse(
                name,
                default_species_prefix=default_species_prefix)
        except:
            pass
    raise AlleleParseError("Unable to parse '%s'" % name)

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
        result_without_mutation,
        mutation_strings,
        default_species_prefix):
    """
    Parameters
    ----------
    original_name : str
        Used for error messages, e.g. "A*02:07 T80M mutant"

    result_without_mutation : FourDigitAllele
        The allele being modified by point mutations in its protein
        sequence

    mutation_strings : list of str
        The mutation descriptors (e.g. "T80M") as well as the optional
        "MUTANT" string used to denote the whole allele as mutated.

    default_species_prefix : str
        If no species can inferred for an allele then use this, e.g. "HLA"

    Returns MutantFourDigitAllele
    """
    parts = name.upper().split()
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
        for p in mutation_strings[1:]
        if p != "MUTANT" and p != ""
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

    TODO: add support for:
        - "HLA class I"
        - "H2-b class I"
        - "ELA-A1 class I"
        - "human CD1a"
        - "H2-r class I"
        - "BF19 class II"
    """
    lower = name.lower()
    if "mutant" in lower:
        return parse_with_mutations(
            name,
            default_species_prefix=default_species_prefix)
    else:
        raise ValueError("Unexpected whitespace in '%s'" % name)

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
        - Locus
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
        raise NotImplementedError(
            "Inference of paired alpha/beta pair for %s not yet implemented" % (
                result,))
    _parse_cache[cache_key] = result
    return result
"""
def parse_allele_name(name, species_prefix=None):
    '''Takes an allele name and splits it into four parts:
        1) species prefix
        2) gene name
        3) allele family
        4) allele code

    If species_prefix is provided, that is used instead of getting the species prefix from the name.
    (And in that case, a species prefix in the name will result in an error being raised.)

    For example, in all of the following inputs:
        "HLA-A*02:01"
        "A0201"
        "A00201"
    The result is a AlleleName object. Example:
        AlleleName(
            species="HLA",  # species prefix
            gene="A",    # gene name
            allele_family="02",   # allele family
            allele_code="01",   # allele code
        )

    The logic for other species mostly resembles the naming system for humans,
    except for mice, rats, and swine, which have archaic nomenclature.
    '''
    original = name
    name = name.strip()

    if len(name) == 0:
        raise ValueError("Can't normalize empty MHC allele name")

    species_from_name, name = split_species_prefix(name)

    if species_prefix:
        if species_from_name:
            raise ValueError(
                ("If a species is passed in, we better not have another "
                 "species in the name itself."))
        species = species_prefix
    else:
        species = species_from_name

    if species in ("H-2", "H2"):
        return parse_mouse_allele_name("H-2-" + name)

    if len(name) == 0:
        raise AlleleParseError("Incomplete MHC allele name: %s" % (original,))
    elif not species:
        # assume that a missing species name means we're dealing with a
        # human HLA allele
        if "-" in name:
            raise AlleleParseError("Can't parse allele name: %s" % original)
        species = "HLA"

    if name[0].upper() == "D":
        if len(name) == 7:
            # sometimes we get very compact names like DRB0101
            gene, name = parse_letters(name, 3)
        else:
            # MHC class II genes like "DQA1" need to be parsed with both
            # letters and numbers
            gene, name = parse_alphanum(name, 4)
        # TODO: make a list of known species/gene pairs, along with
        # gene synonyms. That should significantly imporve on this kind of
        # ad-hoc synonym handling.

        if gene.isalpha():
            # expand e.g. DRA -> DRA1, DQB -> DQB1
            gene = gene + "1"
    elif len(name) == 5:
        # example: SLA-30101
        gene, name = name[0], name[1:]
    elif name[0].isalpha():
        # if there are more separators to come, then assume the gene names
        # can have the form "DQA1"
        gene, name = parse_letters(name)
    elif name[0].isdigit():
        gene, name = parse_numbers(name)
    elif len(name) in (6, 7) and ("*" in name or "-" in name or ":" in name):
        # example: SLA-3*0101 or SLA-3*01:01
        gene, name = parse_alphanum(name)
        _, name = parse_separator(name)
    else:
        raise AlleleParseError(
            "Can't parse gene name from allele: %s" % original)

    if len(gene) == 0:
        raise AlleleParseError("No MHC gene name given in %s" % original)
    if len(name) == 0:
        raise AlleleParseError("Malformed MHC type %s" % original)

    gene = gene.upper()
    # skip initial separator
    sep, name = parse_separator(name)

    if species == "SLA":
        if ":" in name:
            parts = name.split(":")
            if len(parts) != 2:
                raise AlleleParseError(
                    "Unexpected number of ':' characters in '%s'" % original)
            family, name = parts
        elif len(name) < 2:
                raise AlleleParseError("Unable to parse '%s'" % original)
        elif name.isalpha() or len(name) == 2:
            # parse sequences serotypes like SLA-1-HB
            # as shorthand for SLA-1-HB01
            family = name
            name = "01"
        else:
            # the family names for pigs can be weirdly complicated
            # such as 'w13sm' but the alleles still always
            # end with two digits e.g. SLA-2*w13sm20
            family = name[:-2]
            name = name[-2:]
    elif len(name) == 4 or (species == "HLA" and gene in ("A", "B", "C")):
        # If all that's left is e.g. "0201" then only parse the
        # first two digits as the family code. Also, human Class I alleles
        # seem to be exceptional in that they have only 2 digit allele
        # families but 3 digit allele codes
        # (other species like sheep have 3 digits followed by 2 digits)
        family, name = parse_numbers(name, max_len=2)
    else:
        family, name = parse_numbers(name, max_len=3)

    sep, name = parse_separator(name)

    allele_code, rest_of_text = parse_numbers(name)

    rest_of_text = rest_of_text.strip()

    if len(rest_of_text) > 0:
        raise AlleleParseError("The suffix '%s' of '%s' was not parsed" % (
            rest_of_text,
            original))

    if len(family) == 1:
        family = "0" + family
    elif len(family) == 3 and family[0] == "0":
        family = family[1:]
    if len(allele_code) == 0:
        allele_code = "01"
    elif len(allele_code) == 3 and allele_code[0] == "0":
        # normalize HLA-A*02:001 into HLA-A*02:01
        allele_code = allele_code[1:]

    return FourDigitAllele(species, gene, family, allele_code)
"""
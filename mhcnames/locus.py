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

from collections import OrderedDict

from serializable import Serializable

from .species_registry import (
    find_matching_species_prefix,
    find_matching_species_info,
)
from .allele_parse_error import AlleleParseError

class Locus(Serializable):
    def __init__(self, species_prefix, gene_name):
        self.species_prefix = species_prefix
        self.gene_name = gene_name

    _parse_cache = {}

    @classmethod
    def parse_substring(
            cls,
            name,
            default_species_prefix="HLA"):
        """
        Parse locus such as "HLA-A" and return any extra characters
        which follow afterward.
        """
        key = (name, default_species_prefix)
        if key in cls._parse_cache:
            return cls._parse_cache[name]

        num_star_characters = name.count("*")
        if num_star_characters == 1:
            split_index = name.find("*")
            before_star, remaining_string = name[:split_index], name[split_index:]
            locus = Locus.parse(
                before_star,
                default_species_prefix=default_species_prefix)
        elif num_star_characters > 1:
            raise AlleleParseError("Unable to parse locus '%s'" % (name,))
        else num_star_characters == 0:
            # trying to parse prefixes of alleles such as:
            #   HLA-A
            # but also ones with dashes in the species prefix:
            #   H-2-K
            # and also those lacking any dashes such as:
            #   H2K
            # ...and since this is a substring parser, we need to consider
            # that there will also be alleles mixed in, e.g.:
            #  H2Kk
            #  HLA-A0201
            #  YYY-A301 -- where "A3" is a gene on a theoretical species YYY
            for n in [2, 3, 4]:
                candidate_prefix = name[:n]
                species_prefix = find_matching_species_prefix(candidate_prefix)
                if species_prefix is not None:
                    break
            if species_prefix is None:
                species_prefix = default_species_prefix
            if species_prefix is None:
                raise AlleleParseError("Unable to parse '%s'" % name)
            remaining_string = name[n:]
            while remaining_string[0] == "-":
                remaining_string = remaining_string[1:]

            species_info = find_matching_species_info(species_prefix)
            if species_info:
                locus = None
                for n in range(len(remaining_string), 0, -1):
                    substring = remaining_string[:n]
                    gene_name = species_info.find_matching_gene_name(substring)
                    locus = Locus(species_prefix, gene_name)
                    remaining_string = remaining_string[n:]
                    break
            else:
                locus = Locus(species_prefix, remaining_string)
                remaining_string = ""
        result = (locus, remaining_string)
        cls._parse_cache[key] = result
        return result

    @classmethod
    def parse(cls, name, default_species_prefix="HLA"):
        """
        Parse Locus and make sure that there are no extra characters at the end
        of the input sequence.
        """
        result, extra_characters = cls.parse_substring(
            name,
            default_species_prefix=default_species_prefix)
        if len(extra_characters) > 0:
            raise AlleleParseError(
                "Unable to parse '%s', found extra characters '%s' after %s" % (
                    name,
                    extra_characters,
                    result))
        return result

    @property
    def species_info(self):
        return find_matching_species_info(self.species_prefix)

    def get_mhc_class(self):
        return self.species_info.get_mhc_class(self.gene_name)

    def normalized_string(self, include_species=True):
        if include_species:
            return "%s-%s" % (self.species_prefix, self.gene_name)
        else:
            return self.gene_name

    def compact_string(self, include_species=True):
        """
        Compact representation of a Locus, currently same as the
        normalized representation.
        """
        return self.normalized_string(include_species=include_species)

    def to_locus(self):
        """
        Descendant classes use this method to project their fields down
        to a Locus object.
        """
        if self.__class__ is Locus:
            return self
        else:
            return Locus(self.species_prefix, self.gene_name)

    def is_mutant(self):
        return False

    def to_dict(self):
        """
        Returns dictionary with all fields of this locus, its normalized,
        representation and the values of the following methods:
            - is_mutant
            - get_mhc_class
        """
        return OrderedDict([
            ("locus", self.normalized_string()),
            ("species_prefix", self.species_prefix),
            ("gene_name", self.gene_name),
            ("is_mutant", self.is_mutant()),
            ("mhc_class", self.get_mhc_class()),
        ])

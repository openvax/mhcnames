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

from .parse_error import ParseError


def strip_whitespace_and_dashes(s : str):
    while s.startswith("-"):
        s = s[1:]
    while s.startswith(" "):
        s = s[1:]
    while s.endswith("-"):
        s = s[:-1]
    while s.endswith(" "):
        s = s[:-1]
    return s


def strip_whitespace_and_trim_outer_quotes(name : str):
    original_name = name
    name = name.strip()
    while name.startswith('"'):
        if name.endswith('"'):
            name = name[1:-1].strip()
        else:
            raise ParseError(
                "Unbalanced double quotes on allele name: %s" % original_name)
    return name

def split_on_all_seps(seq : str, seps="_:"):
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

def contains_any_letters(s : str):
    """
    Returns True if any characters in the sequence are letters.
    """
    for si in s:
        if si.isalpha():
            return True
    return False

def split_allele_fields(
        original_name,
        str_after_gene,
        allow_three_digits_in_first_field,
        allow_three_digits_in_second_field):
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
                    raise ParseError("Unable to parse '%s'" % original_name)
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
    return parsed_fields
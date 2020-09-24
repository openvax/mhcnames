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

from .allele_parse_error import ParseError


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

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


def _create_regex_for_strip_whitespace_and_dashes():
    optional_space_or_dash = "[-\s]*"
    anything_except_space_or_dash = "A-Za-z0-9\._\*:"
    anything_except_space = anything_except_space_or_dash + "-"
    regex_string_to_strip_spaces_and_dashes = \
        "%s([%s][%s]*[%s])%s" % (
            optional_space_or_dash,
            anything_except_space_or_dash,
            anything_except_space,
            anything_except_space_or_dash,
            optional_space_or_dash
        )
    return re.compile(regex_string_to_strip_spaces_and_dashes)

regex_string_to_strip_spaces_and_dashes = _create_regex_for_strip_whitespace_and_dashes()


def strip_whitespace_and_dashes(s):
    while s.startswith("-"):
        s = s[1:]
    while s.startswith(" "):
        s = s[1:]
    while s.endswith("-"):
        s = s[:-1]
    while s.endswith(" "):
        s = s[:-1]
    return s


def strip_whitespace_and_trim_outer_quotes(name):
    original_name = name
    name = name.strip()
    while name.startswith('"'):
        if name.endswith('"'):
            name = name[1:-1].strip()
        else:
            raise AlleleParseError(
                "Unbalanced double quotes on allele name: %s" % original_name)
    return name

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

from serializable import Serializable


class ParsedResult(Serializable):
    """
    Base class for all parsed objects in mhcnames.
    """
    """
    def _field_name_value_pairs(self):
        results = []
        for field_name in self.field_names():
            field_value = getattr(self, field_name)
            if isinstance(field_value, ParsedResult):
                results.extend(field_value._field_name_value_pairs())
            else:
                results.append((field_name, field_value))
        return results

    def _field_name_string_pairs(self):
        results = []
        for k, v in self._field_name_value_pairs():
            if isinstance(v, str):
                results.append((k, "'%s'" % v))
            else:
                results.append((k, "%s" % v))
        return results
    
    def __str__(self):
        return "%s(%s)" % (
            self.__class__.__name__,
            ", ".join(
                ["%s=%s" % (k, v_str)
                for (k, v_str) in self._field_name_string_pairs()]))
                
    def __repr__(self):
        return str(self)
    """
    def normalized_string(self, include_species=True, use_species_alias=True):
        raise NotImplementedError(
            "%s requires implementation of normalized_string() method" % (
                self.__class__.__name__))

    def compact_string(self, include_species=False, use_species_alias=True):
        """
        Compact representation, defaults to omitting species
        """
        return self.normalized_string(include_species=include_species)

    @classmethod
    def field_names(cls):
        raise NotImplementedError(
            "%s requires implementation of field_names() method" % (
                cls.__name__))

    def __eq__(self, other):
        if self.__class__ is not other.__class__:
            return False
        for field in self.field_names():
            if getattr(self, field) != getattr(other, field):
                return False
        return True

    def __hash__(self):
        return sum(hash(getattr(self, field)) for field in self.field_names())

    def to_record(self):
        raise NotImplementedError(
            "%s requires implementation of to_record() method" % (
                self.__class__.__name__))

    def to_dict(self):
        return dict(zip(self.field_names(), self.to_tuple()))

    def to_tuple(self):
        keys = self.field_names()
        values = [getattr(self, k) for k in keys]
        return tuple(values)

    @classmethod
    def from_tuple(cls, t):
        keys = cls.field_names()
        assert len(keys) == len(t)
        d = dict(zip(keys, t))
        return cls.from_dict(d)

    @classmethod
    def from_dict(cls, d):
        return cls(**d)

    def copy(self, **kwargs):
        """
        Make a copy of this object and update any specified fields.
        """
        field_dict = self.to_dict()
        field_dict.update(kwargs)
        return self.__class__.from_dict(field_dict)

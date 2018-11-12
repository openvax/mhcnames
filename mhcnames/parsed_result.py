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

from serializable import Serializable


class ParsedResult(Serializable):
    """
    Base class for all parsed objects in mhcnames.
    """
    def normalized_string(self, include_species=True):
        raise NotImplementedError(
            "%s requires implementation of normalized_string() method" % (
                self.__class__.__name__))

    def compact_string(self, include_species=False):
        """
        Compact representation, defaults to omitting species
        """
        return self.normalized_string(include_species=include_species)

    @classmethod
    def field_names(cls):
        raise NotImplementedError(
            "%s requires implementation of field_names() method" % (
                cls.__name__))

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
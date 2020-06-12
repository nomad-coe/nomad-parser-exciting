# Copyright 2016-2018 The NOMAD Developers Group
# Copyright 2017-2018 Lorenzo Pardini
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from builtins import object
import xml.sax
from nomadcore.simple_parser import mainFunction, CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion import unit_conversion
import os, sys, json

################################################################
# This is the subparser for the exciting XS output
################################################################

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class EPSParser(object):
    """context for wien2k In2 parser"""

    def __init__(self):
        pass

    def startedParsing(self, path, parser):
        """called when parsing starts"""
        self.parser = parser
        # allows to reset values if the same superContext is used
        # to parse different files
        # self.initialize_values()

    def parseEpsilon(self, epsFile, backend, epsEn, epsilon):
        with open(epsFile) as g:
            while 1:
                s = g.readline()
                if not s: break
                s = s.strip()
                s = s.split()
                if len(s) == 0:
                    break
                else:
                    epsEn[-1].append(float(s[0]))
                    epsilon[-1][0].append(float(s[1]))
                    epsilon[-1][1].append(float(s[2]))

def get_cachingLevelForMetaName(metaInfoEnv, CachingLvl):
    cachingLevelForMetaName = {}
    return cachingLevelForMetaName


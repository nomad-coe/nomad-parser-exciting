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
import setup_paths
import xml.sax
from nomadcore.simple_parser import mainFunction, CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion import unit_conversion
import os, sys, json, logging

###########################################################################
# This is the subparser for the excitons, oscillator strength etc of the XS
###########################################################################


class XSParser(object):

    def __init__(self):
        pass

    def startedParsing(self, path, parser):
        self.parser = parser

    def parseExciton(self, excFile, backend, excNum, excEn, excBindEn, osclStr, transCoeff):
        with open(excFile) as g:
            while 1:
                s = g.readline()
                if not s: break
                s = s.strip()
                s = s.split()
                if s[0] != "#":
                    excNum[-1].append(int(s[0]))
                    excEn[-1].append(float(s[1]))
                    excBindEn[-1].append(float(s[2]))
                    osclStr[-1].append(float(s[3]))
                    transCoeff[-1][0].append(float(s[4]))
                    transCoeff[-1][1].append(float(s[5]))

def get_cachingLevelForMetaName(metaInfoEnv, CachingLvl):
    cachingLevelForMetaName = {}
    return cachingLevelForMetaName


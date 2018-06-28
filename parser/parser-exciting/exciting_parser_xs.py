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

################################################################
# This is the subparser for the exciting GW output
################################################################


class XSParser(object):
    """context for wien2k In2 parser"""

    def __init__(self):
        pass
#        self.excitonEnergys = []
#        self.vertexLabels = []
#        self.vertexNum = 0
#        self.parser = None
#        self.secSingleConfIndex = None

    def startedParsing(self, path, parser):
        """called when parsing starts"""
        self.parser = parser
        # allows to reset values if the same superContext is used to parse different files
#        self.initialize_values()

    def parseExciton(self, excFile, backend, excNum, excEn, excBindEn, osclStr, transCoeff):
#        if dftSingleConfigurationGindex is not None:
#            backend.openNonOverlappingSection("section_calculation_to_calculation_refs")
#            backend.addValue("calculation_to_calculation_ref", dftSingleConfigurationGindex)
#            backend.addValue("calculation_to_calculation_kind", "starting_point")
#            backend.closeNonOverlappingSection("section_calculation_to_calculation_refs")
#        logging.error("GW onClose_section_single_configuration_calculation")
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

#        print("s===",excNum)
#          if len(st) == 3:
 
#           if len(s) >= 40:

#        print("vaaa")
#        self.gmaxvr = float(gmaxvr[0])
#        self.unitCellVol = float(unitCellVol[0])
#        backend.openNonOverlappingSection("section_single_configuration_calculation")
#        if dftSingleConfigurationGindex is not None:
#            backend.openNonOverlappingSection("section_calculation_to_calculation_refs")
#            backend.addValue("calculation_to_calculation_ref", dftSingleConfigurationGindex)
#            backend.addValue("calculation_to_calculation_kind", "starting_point")
#            backend.closeNonOverlappingSection("section_calculation_to_calculation_refs")
#
#        dirPath = os.path.dirname(gwFile)
#        if os.access(os.path.join(dirPath, "EVALQP.DAT"), os.F_OK):
 #           eigvalGWFile = os.path.join(dirPath, "EVALQP.DAT")
#        elif os.access(os.path.join(dirPath, "EVALQP.TXT"), os.F_OK):
#            eigvalGWFile = os.path.join(dirPath, "EVALQP.TXT")
 #       else:
#            pass
#        dosGWFile = os.path.join(dirPath, "TDOS-QP.OUT")
 #       bandCarbGWFile = os.path.join(dirPath, "bandstructure-qp.dat")
 #       bandBorGWFile = os.path.join(dirPath, "BAND-QP.OUT")
 #       vertexGWFile = os.path.join(dirPath, "BANDLINES.OUT")
 #       vertexLabGWFile = os.path.join(dirPath, "bandstructure.xml")
 #       selfCorGWFile = os.path.join(dirPath, "SELFC.DAT")
 #       inputFile = os.path.join(dirPath, "input.xml")
 #       inputgw1File = os.path.join(dirPath, "input-gw.xml")
#        inputgw2File = os.path.join(dirPath, "input_gw.xml")

def get_cachingLevelForMetaName(metaInfoEnv, CachingLvl):
    cachingLevelForMetaName = {}
#                                'section_single_configuration_calculation': CachingLvl
#                               }
#    cachingLevelForMetaName["gw_fundamental_gap"] = CachingLevel.Cache
#    cachingLevelForMetaName["gw_optical_gap"] = CachingLevel.Cache
#    cachingLevelForMetaName["gw_fermi_energy"] = CachingLevel.Cache
    return cachingLevelForMetaName


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
import logging
from nomadcore.simple_parser import mainFunction, CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion import unit_conversion
import os, sys, json
import excitingparser.exciting_parser_input as exciting_parser_input

################################################################
# This is the subparser for the exciting GW output
################################################################


class GWParser(object):
    """context for wien2k In2 parser"""

    def __init__(self):
        self.spinTreat = None
        self.vertexDist = []
        self.vertexLabels = []
        self.vertexNum = 0
        self.parser = None
        self.secSingleConfIndex = None

    def startedParsing(self, path, parser):
        """called when parsing starts"""
        self.parser = parser
        # allows to reset values if the same superContext is used
        # to parse different files
        #s self.initialize_values()

    def parseGW(self, gwFile, backend, dftMethodSectionGindex,
                dftSingleConfigurationGindex, xcName, unitCellVol, gmaxvr):
        #  logging.error("GW onClose_section_single_configuration_calculation")
        self.gmaxvr = float(gmaxvr[0])
        self.unitCellVol = float(unitCellVol[0])
        backend.openNonOverlappingSection("section_single_configuration_calculation")
        if dftSingleConfigurationGindex is not None:
            backend.openNonOverlappingSection("section_calculation_to_calculation_refs")
            backend.addValue("calculation_to_calculation_ref", dftSingleConfigurationGindex)
            backend.addValue("calculation_to_calculation_kind", "starting_point")
            backend.closeNonOverlappingSection("section_calculation_to_calculation_refs")

        dirPath = os.path.dirname(gwFile)
        if os.access(os.path.join(dirPath, "EVALQP.DAT"), os.F_OK):
            eigvalGWFile = os.path.join(dirPath, "EVALQP.DAT")
        elif os.access(os.path.join(dirPath, "EVALQP.TXT"), os.F_OK):
            eigvalGWFile = os.path.join(dirPath, "EVALQP.TXT")
        else:
            pass
        dosGWFile = os.path.join(dirPath, "TDOS-QP.OUT")
        bandCarbGWFile = os.path.join(dirPath, "bandstructure-qp.dat")
        bandBorGWFile = os.path.join(dirPath, "BAND-QP.OUT")
        vertexGWFile = os.path.join(dirPath, "BANDLINES.OUT")
        vertexLabGWFile = os.path.join(dirPath, "bandstructure.xml")
        selfCorGWFile = os.path.join(dirPath, "SELFC.DAT")
        inputFile = os.path.join(dirPath, "input.xml")
        inputgw1File = os.path.join(dirPath, "input-gw.xml")
        inputgw2File = os.path.join(dirPath, "input_gw.xml")

        if os.path.exists(inputFile):
            if os.path.exists(inputgw1File):
                inputgwFile = inputgw1File
            elif os.path.exists(inputgw2File):
                inputgwFile = inputgw2File
            else:
                inputgwFile = inputFile
        elif os.path.exists(inputgw1File):
            inputgwFile = inputgw1File
        elif os.path.exists(inputgw2File):
            inputgwFile = inputgw2File

        if os.path.exists(inputgwFile):
            selfGWSetGIndex = backend.openSection("section_method")
            backend.addValue('electronic_structure_method', "G0W0")
            backend.addValue('gw_starting_point', xcName)
            if dftMethodSectionGindex is not None:
                m2mGindex = backend.openNonOverlappingSection("section_method_to_method_refs")
                backend.addValue("method_to_method_ref", dftMethodSectionGindex)
                backend.addValue("method_to_method_kind", "starting_point")
                backend.closeNonOverlappingSection("section_method_to_method_refs")
                with open(inputgwFile) as f:
                    exciting_parser_input.parseInput(f, backend, self.gmaxvr)
            backend.closeSection("section_method",selfGWSetGIndex)

        if os.path.exists(vertexGWFile):
            with open(vertexGWFile) as g:
                while 1:
                    s = g.readline()
                    if not s: break
                    s = s.strip()
                    s = s.split()
                    if len(s) > 0:
                        if not self.vertexDist:
                            self.vertexDist.append(float(s[0]))
                        elif float(s[0]) != self.vertexDist[-1]:
                            self.vertexDist.append(float(s[0]))
                self.vertexNum = len(self.vertexDist)-1

        if os.path.exists(vertexLabGWFile):
            with open(vertexLabGWFile) as g:
                while 1:
                    s = g.readline()
                    if not s: break
                    s = s.strip()
                    s = s.split()
                    if s[0] == "<vertex":
                        f = s[4].split("\"")
                        self.vertexLabels.append(f[1])

        if os.path.exists(eigvalGWFile):
            eigvalGWGIndex = backend.openSection("section_eigenvalues")
            with open(eigvalGWFile) as g:
                qpGWKpoint=[]
                # For full names of these variables look below where
                # they are added to backaend.
                Vxc = [[],[]]
                Sx = [[],[]]
                Sc = [[],[]]
                qpE = [[],[]]
                Znk = [[],[]]
                fromH = unit_conversion.convert_unit_function("hartree", "J")
                lines_remaining_in_file = os.path.getsize(eigvalGWFile)
                while 1:
                    s = g.readline()
                    lines_remaining_in_file -= len(s)
                    if not s: break
                    s = s.strip()
                    if "k-point" in s.split():
                        qpGWKpoint.append([])
                        for i in range(0, 2):
                            Vxc[i].append([])
                            Sx[i].append([])
                            Sc[i].append([])
                            qpE[i].append([])
                            Znk[i].append([])
                        x,y,z = float(s.split()[3]),float(s.split()[4]),float(s.split()[5])
                        qpGWKpoint[-1].append(x)
                        qpGWKpoint[-1].append(y)
                        qpGWKpoint[-1].append(z)
                    else:
                        s=s.split()
                        if len(s) == 0:
                            continue
                        else:
                            try: not int(s[0])
                            except ValueError:
                                continue
                            if self.spinTreat:
                                pass
                            else:
                                for i in range(0, 2):
                                    try:
                                        qpE[i][-1].append(fromH(float(s[3])))
                                        Sx[i][-1].append(fromH(float(s[4])))
                                        Sc[i][-1].append(fromH(float(s[5])))
                                        Vxc[i][-1].append(fromH(float(s[6])))
                                        Znk[i][-1].append(float(s[9]))
                                    except IndexError:
                                        if not lines_remaining_in_file:
                                            logging.error(
                                                "Last line of GW evalqp file incomplete.")
                                        else:
                                            logging.error(
                                                "Non-last line of GW evalqp file"
                                                " incomplete.")

        backend.addValue("eigenvalues_kpoints", qpGWKpoint)
        backend.addValue("number_of_eigenvalues", len(qpE[0]))
        backend.addValue("number_of_eigenvalues_kpoints", len(qpGWKpoint))
        backend.addValue("eigenvalues_values", qpE)
        backend.addValue("gw_qp_linearization_prefactor", Znk)
        backend.closeSection("section_eigenvalues",eigvalGWGIndex)
        backend.addValue("gw_self_energy_x", Sx)
        backend.addValue("gw_self_energy_c", Sc)
        backend.addValue("gw_xc_potential", Vxc)

        #  ############# DOS  ##############
        if os.path.exists(dosGWFile):
            dosGWGIndex = backend.openSection("section_dos")
            ha_per_joule = unit_conversion.convert_unit(1, "hartree", "J")
            fromH = unit_conversion.convert_unit_function("hartree", "J")
            with open(dosGWFile) as g:
                dosValues = [[],[]]
                dosEnergies = []
                while 1:
                    s = g.readline()
                    if not s: break
                    s = s.strip()
                    s = s.split()
                    ene, value = fromH(float(s[0])), float(s[1])*ha_per_joule
                    dosEnergies.append(ene)
                    if not self.spinTreat:
                        for i in range(0,2):
                            dosValues[i].append(value)
                    else:
                        pass
            backend.addValue("dos_energies", dosEnergies)
            backend.addValue("dos_values", dosValues)
            backend.addValue("number_of_dos_values", len(dosEnergies))
            backend.closeSection("section_dos",dosGWGIndex)

        # ############### BANDSTRUCTURE ###############
        if os.path.exists(bandCarbGWFile):
            bandGWGIndex = backend.openSection("section_k_band")
            fromH = unit_conversion.convert_unit_function("hartree", "J")

            with open(bandCarbGWFile) as g:
                bandEnergies = [[],[]]
                kpoint = []
                dist = []
                Kindex = [0]
                segmK = []
                segmLength = []
                bandEnergiesSegm = []
                bandGWBE = []
                while 1:
                    s = g.readline()
                    if not s: break
                    s = s.strip()
                    s = s.split()
                    if not self.spinTreat:
                        if len(s) == 0:
                            for i in range(0,2):
                                bandEnergies[i].append([])
                        elif s[0] == "#":
                            for i in range(0,2):
                                bandEnergies[i].append([])
                                numBand = int(s[2])
                                numK = int(s[3])
                        elif len(s) > 0:
                            for i in range(0,2):
                                bandEnergies[i][-1].append(fromH(float(s[6])))
                            if int(s[0]) == 1:
                                kpoint.append([])
                                dist.append(float(s[5]))
                                kpoint[-1].append([float(s[2]),float(s[3]),float(s[4])])
                    else:
                        pass

                for i in range(0,2):
                    bandEnergies[i].pop()

                for i in range(1,numK):
                    if dist[i] == dist[i-1]:
                        Kindex.append(i)
                Kindex.append(numK)
                for i in range(0,len(Kindex)-1):
                    segmK.append(dist[Kindex[i]:Kindex[i+1]])

                for i in range(0,len(segmK)):
                    segmLength.append(len(segmK[i]))
                for i in range(0,2):
                    bandEnergiesSegm.append([])
                    for j in range(0,numBand):
                         bandEnergiesSegm[i].append([])
                         for k in range (0,len(Kindex)-1):
                             bandEnergiesSegm[i][j].append(bandEnergies[i][j][Kindex[k]:Kindex[k+1]])
            for i in range(0,len(Kindex)-1):
                bandGWBE.append([])
                for j in range(0,2):
                    bandGWBE[i].append([])
                    for k in range(0,segmLength[i]):
                        bandGWBE[i][j].append([])
                        for l in range(0,numBand):
                            bandGWBE[i][j][-1].append(bandEnergiesSegm[j][l][i][k])

            for i in range(0,len(Kindex)-1):
                bandGWSegmGIndex = backend.openSection("section_k_band_segment")
                backend.addValue("band_energies", bandGWBE[i])
                backend.closeSection("section_k_band_segment",bandGWSegmGIndex)

            backend.closeSection("section_k_band",bandGWGIndex)

        if os.path.exists(bandBorGWFile) and not os.path.exists(bandCarbGWFile):
            bandGWGIndex = backend.openSection("section_k_band")
            fromH = unit_conversion.convert_unit_function("hartree", "J")
            with open(bandBorGWFile) as g:
                bandEnergies = [[[]],[[]]]
                kappa = [[[]],[[]]]
                dist1 = [[]]
                Kindex = [0]
                segmK = []
                segmLength = []
                bandEnergiesSegm = []
                bandGWBE = []
                while 1:
                    s = g.readline()
                    if not s: break
                    s = s.strip()
                    s = s.split()
                    if not self.spinTreat:
                        if len(s) == 0:
                            for i in range(0,2):
                                bandEnergies[i].append([])
                                kappa[i].append([])
                            dist1.append([])
                        elif len(s) > 0:
                            for i in range(0,2):
                                bandEnergies[i][-1].append(fromH(float(s[1])))
                                kappa[i][-1].append(float(s[0]))
                            dist1[-1].append(float(s[0]))
                        numK = len(kappa[0][0])
                        for i in kappa[0][0]:
                            if kappa[0][0].count(i) > 1:
                                kappa[0][0].remove(i)
                    else:
                        pass
                for i in range(0,2):
                    bandEnergies[i].pop()
                numBand = len(bandEnergies[0])
                for i in range(1,numK + self.vertexNum-1):
                    if dist1[0][i] == dist1[0][i-1]:
                        Kindex.append(i)
                Kindex.append(numK + self.vertexNum-1)
                for i in range(0,len(Kindex)-1):
                    segmK.append(dist1[0][Kindex[i]:Kindex[i+1]])

                for i in range(0,len(segmK)):
                    segmLength.append(len(segmK[i]))
                for i in range(0,2):
                    bandEnergiesSegm.append([])
                    for j in range(0,numBand):
                         bandEnergiesSegm[i].append([])
                         for k in range (0,len(Kindex)-1):
                             bandEnergiesSegm[i][j].append(bandEnergies[i][j][Kindex[k]:Kindex[k+1]])
            for i in range(0,len(Kindex)-1):
                bandGWBE.append([])
                for j in range(0,2):
                    bandGWBE[i].append([])
                    for k in range(0,segmLength[i]):
                        bandGWBE[i][j].append([])
                        for l in range(0,numBand):
                            bandGWBE[i][j][-1].append(bandEnergiesSegm[j][l][i][k])

            for i in range(0,len(Kindex)-1):
                bandGWSegmGIndex = backend.openSection("section_k_band_segment")
                backend.addValue("band_energies", bandGWBE[i])
                backend.closeSection("section_k_band_segment",bandGWSegmGIndex)

            backend.closeSection("section_k_band",bandGWGIndex)
        backend.closeNonOverlappingSection("section_single_configuration_calculation")
        # logging.error("done GW onClose_section_single_configuration_calculation")

    # def onOpen_section_method(self, backend, gIndex, section):
    #   fava = section["gw_frequency_number"]
    #   print("fava=",fava)

    # def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
    #   if dftSingleConfigurationGindex is not None:
    #     if self.secSingleConfIndex is None:
    #       self.secSingleConfIndex = gIndex
    #       singleGIndex = backend.openSection("section_single_configuration_calculation")
    #     fermi = section["gw_fermi_energy"]
    #     fundamental = section["gw_fundamental_gap"]
    #     optical = section["gw_optical_gap"]
    #     backend.addValue("gw_fermi_energy", fermi)
    #     backend.addValue("gw_fundamental_gap", fundamental)
    #     backend.addValue("gw_optical_gap", optical)
    #     backend.closeSection("section_single_configuration_calculation",singleGIndex)
    #   else:
    #     singleGIndex = backend.openSection("section_single_configuration_calculation")
    #     fermi = section["gw_fermi_energy"]
    #     fundamental = section["gw_fundamental_gap"]
    #     optical = section["gw_optical_gap"]
    #     backend.addValue("gw_fermi_energy", fermi)
    #     backend.addValue("gw_fundamental_gap", fundamental)
    #     backend.addValue("gw_optical_gap", optical)
    #     backend.closeSection("section_single_configuration_calculation",singleGIndex)

def buildGWMatchers():
    return SM(
    name = 'root',
    weak = True,
    startReStr = "\=\s*Main GW output file\s*\=",
    # sections = ["section_run"],
    subMatchers = [
    SM(
      startReStr = "\-\s*frequency grid\s*\-",
      endReStr = "\-\s*Peak memory estimate \(Mb, per process\)\:\s*\-",
      sections = ["section_method"],
      subMatchers = [
      # SM(r"\s*Type\:\s*\<\s*(?P<gw_dummy>[-a-zA-Z]+)\s*\>"),
        SM(r"\s*(?P<gw_frequency_number>[0-9]+)\s*(?P<gw_frequency_values__hartree>[0-9]\.[0-9]*([E]?[-]?[-0-9]+))\s*(?P<gw_frequency_weights>[0-9]\.[0-9]*([E]?[-]?[-0-9]+))", repeats = True)
    ]),
    SM(
      startReStr = "\-\s*G0W0\s*\-",
      endReStr = "\=\s*GW timing info \(seconds\)\s*\=",
      sections = ["section_single_configuration_calculation"],
      subMatchers = [
        SM(r"\s*Fermi energy\:\s*(?P<gw_fermi_energy__hartree>[-0-9.]+)"),
        SM(r"\s* Fundamental BandGap \(eV\)\:\s*(?P<gw_fundamental_gap__eV>[0-9.]+)"),
        SM(r"\s*Direct BandGap \(eV\)\:\s*(?P<gw_fundamental_gap__eV>[0-9.]+)"),
        SM(r"\s* Optical BandGap \(eV\)\:\s*(?P<gw_optical_gap__eV>[0-9.]+)")
    ])
    ])
def get_cachingLevelForMetaName(metaInfoEnv, CachingLvl):
    cachingLevelForMetaName = {}
    # 'section_single_configuration_calculation': CachingLvl
    # }
    # cachingLevelForMetaName["gw_fundamental_gap"] = CachingLevel.Cache
    # cachingLevelForMetaName["gw_optical_gap"] = CachingLevel.Cache
    # cachingLevelForMetaName["gw_fermi_energy"] = CachingLevel.Cache
    return cachingLevelForMetaName


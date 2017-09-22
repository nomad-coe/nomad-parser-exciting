from builtins import object
import setup_paths
import xml.sax
from nomadcore.simple_parser import mainFunction, CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion import unit_conversion
import os, sys, json

################################################################
# This is the subparser for the exciting GW output
################################################################


class GWContext(object):
    """context for wien2k In2 parser"""

    def __init__(self):
        self.parser = None
        self.spinTreat = None
        self.vertexDist = []
        self.vertexLabels = []
        self.vertexNum = 0

    def initialize_values(self):
        """allows to reset values if the same superContext is used to parse different files"""
        pass

    def startedParsing(self, path, parser):
  
        """called when parsing starts"""
        self.parser = parser
        self.initialize_values()

    def onClose_x_exciting_section_GW(self, backend, gIndex, section):

        dirPath = os.path.dirname(self.parser.fIn.name)
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

        if os.path.exists(inputFile):
            selfGWSetGIndex = backend.openSection("x_exciting_section_GW_settings")
            singularity = 'mpd'
            actype = 'pade'
            npol = 0
            scrtype = "rpa"
            snempty = 0
            coreflag = "all"
            fgrid = "gaule2"
            k1 = 0
            k2 = 0
#            f1 = 0
#            f2 = 0
            s1 = 0
            s2 = 0
#            m1 = 0
#            m2 = 0
#            bc1 = 0
#            bc2 = 0
#            sc1 = 0
#            sc2 = 0
            with open(inputFile) as g:
                i = 0
                while 1:
                    s = g.readline()
                    i += 1
                    if not s: break
                    s = s.strip()
                    s = s.split('=')
                    if s[0] == "<gw": k1 = i
                    if s[0] == "</gw>": k2 = i
#                    if s[0] == "<freqgrid": f1 = i
#                    if s[0] == "</freqgrid>": f2 = i
                    if s[0] == "<selfenergy": s1 = i
                    if s[0] == "</selfenergy>": s2 = i
#                    if s[0] == "<mixbasis": m1 = i
#                    if s[0] == "</mixbasis>": m2 = i
#                    if s[0] == "<barecoul": bc1 = i
#                    if s[0] == "</barecoul>": bc2 = i
#                    if s[0] == "<scrcoul": sc1 = i
#                    if s[0] == "</scrcoul>": sc2 = i
            with open(inputFile) as g:
                i = 0
                while 1:
                    s = g.readline()
                    i += 1
                    if not s: break
                    s = s.strip()
                    s = s.split('=')
                    if (s[0] == "coreflag") and (i >= k1):
                        coreflag = s[1][1:-1]
                    if (s[0] == "singularity") and (i >= k1):
                        freq_conv = s[1][1:-1]
                    if (s[0] == "actype") and (i >= k1):
                        actype = s[1][1:-1]
                    if (s[0] == "npol") and (i >= k1):
                        npol = s[1][1:-1]
                    if (s[0] == "scrtype") and (i >= k1):
                        scrtype = s[1][1:-1]
                    if (s[0] == "nempty") and (i >= s1) and (i <= s2):
                        snempty = s[1][1:-1]
                    if (s[0] == "fgrid") and (i >= k1):
                        fgrid = s[1][1:-1]
            backend.addValue("x_exciting_GW_frequency_grid_type", fgrid)
            backend.addValue("x_exciting_GW_self_energy_c_empty_states", int(snempty))
            backend.addValue("x_exciting_GW_core_treatment", coreflag)
            backend.addValue("x_exciting_GW_self_energy_singularity_treatment", singularity)
            backend.addValue("x_exciting_GW_self_energy_c_analytical_continuation", actype)
            backend.addValue("x_exciting_GW_self_energy_c_number_of_poles", int(npol))
            backend.addValue("x_exciting_GW_screened_Coulomb", scrtype)
            backend.closeSection("x_exciting_section_GW_settings",selfGWSetGIndex)

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
            eigvalGWOutGIndex = backend.openSection("x_exciting_section_GW_output")
            eigvalGWGIndex = backend.openSection("x_exciting_section_GW_qp_eigenvalues")
            with open(eigvalGWFile) as g:
                qpGWKpoint=[]
                Sx = [[],[]]
                Sc = [[],[]]
                qpE = [[],[]]
                Znk = [[],[]]
                fromH = unit_conversion.convert_unit_function("hartree", "J")
                while 1:
                    s = g.readline()
                    if not s: break
                    s = s.strip()
                    if "k-point" in s.split():
                        qpGWKpoint.append([])
                        for i in range(0,2):
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
                                for i in range(0,2):
                                    Sx[i][-1].append(fromH(float(s[4])))
                                    Sc[i][-1].append(fromH(float(s[5])))
                                    qpE[i][-1].append(fromH(float(s[3])))
                                    Znk[i][-1].append(float(s[9]))
        backend.addValue("x_exciting_GW_qp_eigenvalues_kpoints", qpGWKpoint)
        backend.addValue("x_exciting_GW_qp_number_of_eigenvalues", len(qpE[0]))
        backend.addValue("x_exciting_GW_qp_number_of_eigenvalues_kpoints", len(qpGWKpoint))
        backend.addValue("x_exciting_GW_qp_eigenvalues_values", qpE)
        backend.addValue("x_exciting_GW_qp_linearization_prefactor", Znk)
        backend.closeSection("x_exciting_section_GW_qp_eigenvalues",eigvalGWGIndex)

        selfGWGIndex = backend.openSection("x_exciting_section_GW_self_energy")
        backend.addValue("x_exciting_GW_self_energy_x", Sx)
        backend.addValue("x_exciting_GW_self_energy_c", Sc)
        backend.closeSection("x_exciting_section_GW_self_energy",selfGWGIndex)
        backend.closeSection("x_exciting_section_GW_output",eigvalGWOutGIndex)
####################DOS######################

        if os.path.exists(dosGWFile):
            dosGWOutGIndex = backend.openSection("x_exciting_section_GW_output")
            dosGWGIndex = backend.openSection("x_exciting_section_GW_dos")
            fromH = unit_conversion.convert_unit_function("hartree", "J")
            with open(dosGWFile) as g:
                dosValues = [[],[]]
                dosEnergies = []
                while 1:
                    s = g.readline()
                    if not s: break
                    s = s.strip()
                    s = s.split()
                    ene, value = fromH(float(s[0])), float(s[1])
                    dosEnergies.append(ene)
                    if not self.spinTreat:
                        for i in range(0,2):
                            dosValues[i].append(value)
                    else:
                        pass
            backend.addValue("x_exciting_GW_dos_energies", dosEnergies)
            backend.addValue("x_exciting_GW_dos_values", dosValues)
            backend.addValue("x_exciting_GW_number_of_dos_values", len(dosEnergies))
            backend.closeSection("x_exciting_section_GW_dos",dosGWGIndex)        
            backend.closeSection("x_exciting_section_GW_output",dosGWOutGIndex)
##################BANDSTRUCTURE#####################

        if os.path.exists(bandCarbGWFile):
            bandGWOutGIndex = backend.openSection("x_exciting_section_GW_output")
            bandGWGIndex = backend.openSection("x_exciting_section_GW_k_band")
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
                bandGWSegmGIndex = backend.openSection("x_exciting_section_GW_k_band_segment")
                backend.addValue("x_exciting_GW_band_energies", bandGWBE[i])
                backend.closeSection("x_exciting_section_GW_k_band_segment",bandGWSegmGIndex)

            backend.closeSection("x_exciting_section_GW_k_band",bandGWGIndex)
            backend.closeSection("x_exciting_section_GW_output",bandGWOutGIndex)

        if os.path.exists(bandBorGWFile) and not os.path.exists(bandCarbGWFile):
            bandGWOutGIndex = backend.openSection("x_exciting_section_GW_output")
            bandGWGIndex = backend.openSection("x_exciting_section_GW_k_band")
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
                bandGWSegmGIndex = backend.openSection("x_exciting_section_GW_k_band_segment")
                backend.addValue("x_exciting_GW_band_energies", bandGWBE[i])
                backend.closeSection("x_exciting_section_GW_k_band_segment",bandGWSegmGIndex)

            backend.closeSection("x_exciting_section_GW_k_band",bandGWGIndex)
            backend.closeSection("x_exciting_section_GW_output",bandGWOutGIndex)

def buildGWMatchers():
    return SM(
    name = 'root',
    weak = True,
    startReStr = "\*\s*GW input parameters\s*\*",
        sections = ["x_exciting_section_GW","x_exciting_section_GW_settings"],
    subMatchers = [
        SM(r"\s*Number of empty states \(GW\):\s*(?P<x_exciting_GW_polarizability_empty_states>[0-9]+)\s*")
    ])

def get_cachingLevelForMetaName(metaInfoEnv, CachingLvl):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
        CachingLvl: Sets the CachingLevel for the sections k_band, run, and single_configuration_calculation.
            This allows to run the parser without opening new sections.

    Returns:
        Dictionary with metaname as key and caching level as value.
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = {
                               'section_run': CachingLvl,
                               'section_method': CachingLvl
                              }
    return cachingLevelForMetaName


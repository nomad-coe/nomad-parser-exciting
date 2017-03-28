from builtins import object
import setup_paths
import xml.sax
from nomadcore.simple_parser import mainFunction, CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion import unit_conversion
import os, sys, json

################################################################
# This is the subparser for the WIEN2k input file (.in2)
################################################################


class GWContext(object):
    """context for wien2k In2 parser"""

    def __init__(self):
        self.parser = None
        self.spinTreat = None
        self.vertexDist = []
        self.vertexLabels = []

    def initialize_values(self):
        """allows to reset values if the same superContext is used to parse different files"""
        pass

    def startedParsing(self, path, parser):
  
        """called when parsing starts"""
        self.parser = parser
        self.initialize_values()

    def onClose_x_exciting_section_GW_method(self, backend, gIndex, section):

        dirPath = os.path.dirname(self.parser.fIn.name)
        eigvalGWFile = os.path.join(dirPath, "EVALQP.DAT")
        dosGWFile = os.path.join(dirPath, "TDOS-QP.OUT")
        bandGWFile = os.path.join(dirPath, "bandstructure-qp.dat")
        vertexGWFile = os.path.join(dirPath, "BANDLINES.OUT")
        vertexLabGWFile = os.path.join(dirPath, "bandstructure.xml")

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
                        x,y,z,weight = float(s.split()[3]),float(s.split()[4]),float(s.split()[5]),float(s.split()[6])
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
                                    Sx[i][-1].append(float(s[4]))
                                    Sc[i][-1].append(float(s[5]))
                                    qpE[i][-1].append(float(s[3]))
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

####################DOS######################

        if os.path.exists(dosGWFile):
            dosGWGIndex = backend.openSection("x_exciting_section_GW_dos")
            with open(dosGWFile) as g:
                dosValues = [[],[]]
                dosEnergies = []
                while 1:
                    s = g.readline()
                    if not s: break
                    s = s.strip()
                    s = s.split()
                    ene, value = float(s[0]), float(s[1])
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

##################BANDSTRUCTURE#####################

        if os.path.exists(bandGWFile):
            bandGWGIndex = backend.openSection("x_exciting_section_GW_k_band")
            bandGWSegmGIndex = backend.openSection("x_exciting_section_GW_k_band_segment")
            fromH = unit_conversion.convert_unit_function("hartree", "J")

            with open(bandGWFile) as g:
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
#                                ene = fromH(float(s[6]))
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
#                    print("i=",i)
#                    print("dist[i-1]=",dist[i-1])
#                    print("dist[i]=",dist[i])
                    if dist[i] == dist[i-1]:
#                        pass
                        Kindex.append(i)
                Kindex.append(numK)
#                        del dist[i]
                for i in range(0,len(Kindex)-1): 
                    segmK.append(dist[Kindex[i]:Kindex[i+1]])

                for i in range(0,len(segmK)):
#                    print("i=",i)
#                    print("segmK[i]=",segmK[i])
                    segmLength.append(len(segmK[i]))
#                    bandEnergiesSegm.append([])
                for i in range(0,2):
                    bandEnergiesSegm.append([])
                    for j in range(0,numBand):
                         bandEnergiesSegm[i].append([])
                         for k in range (0,len(Kindex)-1):
#                             print("l=",k)
#                             print("len(Kindex)=",len(Kindex))
#                             print("Kindex[l]=",Kindex[k])
#                             print("bandEnergies[Kindex[k]:Kindex[k+1]]=",bandEnergies[Kindex[k]:Kindex[k+1]])
                             bandEnergiesSegm[i][j].append(bandEnergies[i][j][Kindex[k]:Kindex[k+1]])
#                    print("i=",i)
#                    print("Kindex[i]=",Kindex[i])
#                    del dist[Kindex[i]]

#            print("bandEnergies=",bandEnergies)  
#            print("self.vertexDist=",self.vertexDist)
#            print("dist=",dist)
#            print("len(dist)=",len(dist))
#            print("segmK=",segmK)
#            print("segmLength=",segmLength)
#            print("bandEnergiesSegm=",bandEnergiesSegm)
#            print("Kindex=",Kindex)
            for i in range(0,len(Kindex)-1):
                bandGWBE.append([])
                for j in range(0,2):
                    bandGWBE[i].append([])
                    for k in range(0,segmLength[i]):
                        bandGWBE[i][j].append([])
                        for l in range(0,numBand):
                            bandGWBE[i][j][-1].append(bandEnergiesSegm[j][l][i][k])

#            print("bandGWBE=",bandGWBE)
            for i in range(0,len(Kindex)-1):
                backend.addValue("x_exciting_GW_band_energies", bandGWBE[i])

            backend.closeSection("x_exciting_section_GW_k_band_segment",bandGWSegmGIndex)
            backend.closeSection("x_exciting_section_GW_k_band",bandGWGIndex)
#            backend.closeSection("x_exciting_section_GW_k_band_segment",bandGWSegmGIndex)

def buildGWMatchers():
    return SM(
    name = 'root',
    weak = True,
    startReStr = "\*\s*GW input parameters\s*\*",
        sections = ["x_exciting_section_GW_method", "x_exciting_section_GW"],
    subMatchers = [
#        SM(name = 'GWinput',
#          startReStr = r"(?P<x_wien2k_system_nameIn>.*)"),
#          startReStr = "\s*GW taskname:\s*"
        SM(r"\s*(?P<x_exciting_GW_type>[-a-zA-Z0-9]+)\s*-\s*[-a-zA-Z0-9]+\s*run"),
        SM(r"\s*(?P<x_exciting_GW_type>[-a-zA-Z0-9]+)\s*-\s*[-a-zA-Z0-9]+\s*run"),
#        SM(r"\s*(?P<x_wien2k_in2c_switch>[A-Z]+)\s*.*"),
#        SM(r"\s*(?P<x_wien2k_in2c_emin>[-+0-9.]+)\s*(?P<x_wien2k_in2c_ne>[-+0-9.]+)\s*(?P<x_wien2k_in2c_espermin>[-+0-9.]+)\s*(?P<x_wien2k_in2c_esper0>[-+0-9.]+)\s*.*"),
#        SM(r"\s*(?P<x_wien2k_smearing_kind>[A-Z]+)\s*\s*(?P<x_wien2k_smearing_width__rydberg>[-+0-9.]+)\s*.*"),
#        SM(r"\s*(?P<x_wien2k_in2c_gmax>[-+0-9.]+)\s*GMAX")
#        SM(r"\s*GW taskname\:\s*")
     
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


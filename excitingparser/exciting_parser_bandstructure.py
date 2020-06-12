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

import xml.sax
import logging
import numpy as np
from nomadcore.unit_conversion import unit_conversion

class BandHandler(xml.sax.handler.ContentHandler):
    def __init__(self, backend, spinTreat):
        self.backend = backend
        self.bandSectionGIndex = -1
        self.normBandSectionGIndex = -1
        self.inBand = False
        self.energy=[]
        self.energySpin = [[],[]]
        self.bandEnergies = [[],[]]
        self.distance = []
        self.vertexCoord = []
        self.vertexLabels = []
        self.vertexDist = []
        self.spinTreat = spinTreat

    def endDocument(self):
            self.inBand = False
            self.backend.closeSection("section_k_band",self.bandSectionGIndex)
#            self.backend.closeSection("section_k_band_segment",self.normBandSectionGIndex)
            self.bandSectionGIndex = -1
            self.normBandSectionGIndex = -1

    def startElement(self, name, attrs):
        if name == "bandstructure":
            self.bandSectionGIndex = self.backend.openSection("section_k_band")
            self.backend.addValue("band_structure_kind","electronic")
#            self.normBandSectionGIndex = self.backend.openSection("section_k_band_segment")
            self.inBand = True
        elif name == "band":
            self.energy.append([])
            self.distance.append([])
        elif name == "point" and self.inBand:
            fromH = unit_conversion.convert_unit_function("hartree", "J")
            self.energy[-1].append(fromH(float(attrs.getValue('eval'))))
            self.distance[-1].append(float(attrs.getValue('distance')))
        elif name == "vertex" and self.inBand:
            self.vertexCoord.append(attrs.getValue("coord"))
            label = attrs.getValue("label")
            if (label == "Gamma") or (label == "GAMMA") or (label == "gamma"):
                self.vertexLabels.append('\u0393')
            else:
                self.vertexLabels.append(label)
            self.vertexDist.append(attrs.getValue("distance"))

    def endElement(self, name):
        if name == 'bandstructure':
            bandEnergiesBE = []
            vertexDist=[]
            vertexDummy = []
            kBandSegm = []
            bandKpoints = []
            bandEnergies = []
            bandOcc = []
            step = []
            numkPointsPerSegmL = []
            numkPointsPerSemIncr = []
            coordinate = []
            vertexNum = len(self.vertexLabels)
            kmesh = len(self.energy[-1])
            bands = len(self.energy)
            bands2 = int(bands/2)
            kpoints=self.distance[-1]
            vertexDist.append(kpoints[0])

            i=0
            while i < vertexNum:
                coordinate.append(self.vertexCoord[i].split())
                i +=1

            for i in kpoints:
                if kpoints.count(i)>1:
                    index=kpoints.index(i)
                    vertexDist.append(kpoints[index])

            vertexDist.append(kpoints[-1])

            for i in vertexDist:
                if vertexDist.count(i)>1:
                    index=vertexDist.index(i)
                    prova = vertexDist.pop(index)

            i = 0
            while i < vertexNum-1:
                bandKpoints.append([])
                bandEnergies.append([])
                bandOcc.append([])
                step.append([])
                initial = kpoints.index(vertexDist[i])
                final = kpoints.index(vertexDist[i+1])
                kBandSegm.append(kpoints[kpoints.index(vertexDist[i]):kpoints.index(vertexDist[i+1])])
                i +=1

            j = 0
            while j < vertexNum -2:
                bodda = kBandSegm[j+1].pop(0)
                kBandSegm[j].append(bodda)
                j+=1

            kBandSegm[-1].append(kpoints[-1])
            numkPointsPerSemIncr.append(0)

            for i in range(0,vertexNum-1):
                numkPointsPerSegm = len(kBandSegm[i])
                numkPointsPerSegmL.append(numkPointsPerSegm)
                numkPointsPerSemIncr.append(numkPointsPerSemIncr[-1]+numkPointsPerSegmL[-1])
                step[i].append((float(coordinate[i+1][0])-float(coordinate[i][0]))/(numkPointsPerSegm-1))
                step[i].append((float(coordinate[i+1][1])-float(coordinate[i][1]))/(numkPointsPerSegm-1))
                step[i].append((float(coordinate[i+1][2])-float(coordinate[i][2]))/(numkPointsPerSegm-1))
                for j in range(0,numkPointsPerSegm):
                    bandKpoints[i].append([])
                    for k in range(0,3):
                        bandKpoints[i][j].append(float(coordinate[i][k])+j*step[i][k])

            i=0
            while i < vertexNum-1:
                bandEnergiesBE.append([])
                for j in range(0,2):
                    bandEnergiesBE[i].append([])
                    for k in range(0,numkPointsPerSegmL[i]):
                        bandEnergiesBE[i][j].append([])
                i+=1

            for i in range(0,vertexNum):
                self.vertexCoord[i]=self.vertexCoord[i].split()
                for j in range(0,3):
                    self.vertexCoord[i][j] = float(self.vertexCoord[i][j])

            self.eigenSectionGIndex = self.backend.openSection("section_eigenvalues")
            self.backend.addValue("number_of_band_segment_eigenvalues",bands)
            self.backend.closeSection("section_eigenvalues",self.eigenSectionGIndex)

            if not self.spinTreat:
                self.energySpin[0] = self.energy[0:bands]
                self.energySpin[1] = self.energy[0:bands]
                for i in range (0,bands):
                    self.bandEnergies[0].append([])
                    self.bandEnergies[1].append([])
                    for j in range(0,vertexNum-1):
                        self.bandEnergies[0][i].append(self.energySpin[0][i][numkPointsPerSemIncr[j]:numkPointsPerSemIncr[j+1]])
                        self.bandEnergies[1][i].append(self.energySpin[1][i][numkPointsPerSemIncr[j]:numkPointsPerSemIncr[j+1]])
                for i in range (0,vertexNum-1):
                   self.normBandSectionGIndex = self.backend.openSection("section_k_band_segment")
                   for j in range(0,bands):
                       for k in range(0,numkPointsPerSegmL[i]):
                            bandEnergiesBE[i][0][k].append(self.bandEnergies[0][j][i][k])
                            bandEnergiesBE[i][1][k].append(self.bandEnergies[1][j][i][k])
                   self.backend.addValue("band_k_points",bandKpoints[i])
                   self.backend.addValue("band_segm_start_end",self.vertexCoord[i:i+2])
                   self.backend.addValue("number_of_k_points_per_segment",numkPointsPerSegmL[i])
                   self.backend.addValue("band_segm_labels",self.vertexLabels[i:i+2])
                   self.backend.addValue("band_energies",bandEnergiesBE[i])
                   self.backend.closeSection("section_k_band_segment",self.normBandSectionGIndex)
            else: #### check for spin polarized!!!!
                self.energySpin[0] = self.energy[0:bands2]
                self.energySpin[1] = self.energy[bands2:bands]
                for i in range (0,bands2):
                    self.bandEnergies[0].append([])
                    self.bandEnergies[1].append([])
                    for j in range(0,vertexNum-1):
                        self.bandEnergies[0][i].append(self.energySpin[0][i][numkPointsPerSemIncr[j]:numkPointsPerSemIncr[j+1]])
                        self.bandEnergies[1][i].append(self.energySpin[1][i][numkPointsPerSemIncr[j]:numkPointsPerSemIncr[j+1]])
                for i in range (0,vertexNum-1):
                    self.normBandSectionGIndex = self.backend.openSection("section_k_band_segment")
                    for j in range(0,bands2):
                        for k in range(0,numkPointsPerSegmL[i]):
                            bandEnergiesBE[i][0][k].append(self.bandEnergies[0][j][i][k])
                            bandEnergiesBE[i][1][k].append(self.bandEnergies[1][j][i][k])
                    self.backend.addValue("band_k_points",bandKpoints[i])
                    self.backend.addValue("band_segm_start_end",self.vertexCoord[i:i+2])
                    self.backend.addValue("number_of_k_points_per_segment",numkPointsPerSegmL[i])
                    self.backend.addValue("band_segm_labels",self.vertexLabels[i:i+2])
                    self.backend.addValue("band_energies",bandEnergiesBE[i])
                    self.backend.closeSection("section_k_band_segment",self.normBandSectionGIndex)

    def startElementNS(self, name, qname, attrs):
        attrDict={}
        for name in attrs.getNames():
            attrDict[name] = attrs.getValue(name)
        logging.error("start element %s ns %s attr %s", name, qname, attrDict)

    def endElementNS(self, name, qname):
        logging.error("end element %s ns %s", name, qname)

    def characters(self, content):
        pass

def parseBand(inF, backend, spinTreat):
    handler = BandHandler(backend, spinTreat)
    logging.info("will parse")
    xml.sax.parse(inF, handler)
    logging.info("did parse")

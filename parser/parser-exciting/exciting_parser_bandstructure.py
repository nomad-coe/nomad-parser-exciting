import xml.sax
import logging
import numpy as np

class BandHandler(xml.sax.handler.ContentHandler):
    def __init__(self, backend, spinTreat):
        self.backend = backend
        self.bandSectionGIndex = -1
        self.inBand = False
        self.energy=[]
        self.energySpin = [[],[]]
        self.distance = []
        self.vertexCoord = []
        self.vertexLabels = []
        self.vertexDist = []
        self.spinTreat = spinTreat

    def endDocument(self):
            self.inBand = False
            self.backend.closeSection("x_exciting_section_bandstructure",self.bandSectionGIndex)
            self.bandSectionGIndex = -1

    def startElement(self, name, attrs):
        if name == "bandstructure":
            self.bandSectionGIndex = self.backend.openSection("x_exciting_section_bandstructure")
            self.inBand = True
        elif name == "band":
            self.energy.append([])
            self.distance.append([])
        elif name == "point" and self.inBand:
            self.energy[-1].append(float(attrs.getValue('eval')))
            self.distance[-1].append(float(attrs.getValue('distance')))
        elif name == "vertex" and self.inBand:
            self.vertexCoord.append(attrs.getValue("coord"))
            self.vertexLabels.append(attrs.getValue("label"))
            self.vertexDist.append(attrs.getValue("distance"))

    def endElement(self, name):
        if name == 'bandstructure':
            vertexDummy = []
            vertexNum = len(self.vertexLabels)
            kmesh = len(self.energy[-1])
            bands = len(self.energy)
            bands2 = int(bands/2)
            for i in range(0,vertexNum):
                self.vertexCoord[i]=self.vertexCoord[i].split() 
                for j in range(0,3):
                    self.vertexCoord[i][j] = float(self.vertexCoord[i][j])
            self.backend.addValue("x_exciting_band_number_of_kpoints",kmesh)
            self.backend.addValue("x_exciting_band_number_of_vertices",vertexNum)
            self.backend.addValue("x_exciting_band_vertex_labels",self.vertexLabels)
            self.backend.addValue("x_exciting_band_vertex_coordinates", self.vertexCoord)
            self.backend.addValue("x_exciting_band_k_points",self.distance[-1])
            self.backend.addValue("x_exciting_band_structure_kind","electronic")
            if not self.spinTreat:
                self.energySpin[0] = self.energy[0:bands]
                self.energySpin[1] = self.energy[0:bands]
                self.backend.addValue("x_exciting_band_number_of_eigenvalues",bands)
                self.backend.addValue("x_exciting_band_energies",self.energySpin)
#                print("self.energySpin=",self.energySpin)
            else:
                self.energySpin[0] = self.energy[0:bands2]
                self.energySpin[1] = self.energy[bands2:bands]
                self.backend.addValue("x_exciting_band_number_of_eigenvalues",bands2)
                self.backend.addValue("x_exciting_band_energies",self.energySpin)
                
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
    logging.error("will parse")
    xml.sax.parse(inF, handler)
    logging.error("did parse")

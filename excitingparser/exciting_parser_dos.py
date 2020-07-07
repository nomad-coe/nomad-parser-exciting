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
from nomadcore.unit_conversion.unit_conversion import convert_unit_function
from nomadcore.unit_conversion.unit_conversion import convert_unit
from nomadcore.unit_conversion import unit_conversion


class DosHandler(xml.sax.handler.ContentHandler):
    def __init__(self, backend, spinTreat, unitCellVol, fermiEnergy):
        self.backend = backend
        self.dosSectionGIndex = -1
        self.inDos = False
        self.dosProjSectionGIndex = -1
        self.inDosProj = False
        self.spinTreat = spinTreat
        self.totDos = []
        self.totDosSpin = [[],[]]
        self.dosProj = []
        self.dosProjSpin = []
        self.energy = []
        self.energySpin = []
        self.speciesrn = []
        self.atom = []
        self.numDosVal = 0
        self.dosProjDummy = []
        self.dosProjDummy2 = []
        self.unitCellVol = float(unitCellVol[0])
        self.fermiEnergy = fermiEnergy

    def endDocument(self):
        self.inDosProj = False
        self.backend.closeSection("section_dos",self.dosSectionGIndex)
        self.dosSectionGIndex = -1

    def startElement(self, name, attrs):
        ha_per_joule = convert_unit(1, "hartree", "J")
        joule_in_ev = convert_unit(1, "eV", "J")
        fromH = unit_conversion.convert_unit_function("hartree", "J")
        if name == "totaldos":
            self.dosSectionGIndex = self.backend.openSection("section_dos")
            self.inDos = True
        elif name == "partialdos":
            self.speciesrn.append(int(attrs.getValue('speciesrn')))
            self.atom.append(int(attrs.getValue('atom')))
            self.dosProjDummy.append([])
            self.dosProjDummy2.append([])
            if self.speciesrn [-1] == 1 and self.atom[-1] == 1:
                self.dosProjSectionGIndex = self.backend.openSection("section_atom_projected_dos")
                self.inDosProj = True
        elif name == "point":
            if self.inDos:
                self.totDos.append(float(attrs.getValue('dos'))/ha_per_joule)
                self.energy.append(fromH(float(attrs.getValue('e'))))
            elif self.inDosProj:
                self.dosProj.append(float(attrs.getValue('dos'))/ha_per_joule)
                self.energy.append(fromH(float(attrs.getValue('e'))))
        elif name == "diagram":
            if not self.speciesrn: pass
            elif self.speciesrn [-1] == 1 and self.atom[-1] == 1:
                if not self.spinTreat:
                    self.dosProjSpin.append([])
                else:
                    nspin = int(attrs.getValue("nspin"))
                    if nspin == 1:
                        self.dosProjSpin.append([])
        elif name == "limrep":
            for i in range(0,len(self.dosProjSpin)):
                for j in range(0,2):
                    self.dosProjSpin[i].append([])
                    for k in range(0,len(self.speciesrn)):
                        self.dosProjDummy2[k].append([])
                        self.dosProjSpin[i][j].append([])
            for i in range (0,len(self.speciesrn)):
                if not self.spinTreat:
                    self.dosProjDummy[i] = self.dosProj[i*len(self.dosProjSpin)*self.numDosVal:(i+1)*len(self.dosProjSpin)*self.numDosVal]
                else:
                    self.dosProjDummy[i] = self.dosProj[i*int(2*len(self.dosProjSpin))*self.numDosVal:(i+1)*int(2*len(self.dosProjSpin))*self.numDosVal]
            for j in range(0,int(2*len(self.dosProjSpin))):
                for i in range (0,len(self.speciesrn)):
                    self.dosProjDummy2[i][j] = self.dosProjDummy[i][j*self.numDosVal:(j+1)*self.numDosVal]
            if not self.spinTreat:
                for i in range(0,len(self.dosProjSpin)):
                    for j in range(0,len(self.speciesrn)):
                        self.dosProjSpin[i][0][j] = self.dosProjDummy2[j][i]
                        self.dosProjSpin[i][1][j] = self.dosProjDummy2[j][i]
            else:
                for j in range(0,len(self.speciesrn)):
                    for i in range(0,int(2*len(self.dosProjSpin))):
                        if i < len(self.dosProjSpin):
                            self.dosProjSpin[i][0][j] = self.dosProjDummy2[j][i]
                        else:
                            k = int(i - len(self.dosProjSpin))
                            self.dosProjSpin[k][1][j] = self.dosProjDummy2[j][i]
    def endElement(self, name):
        if name == 'totaldos':
            self.inDos = False
            if not self.spinTreat:
                self.numDosVal = len(self.energy)
                self.totDosSpin[0] = self.totDos[0:self.numDosVal]
                self.totDosSpin[1] = self.totDos[0:self.numDosVal]
                self.energySpin = self.energy[0:self.numDosVal]
                self.backend.addValue("dos_values", np.array(self.totDosSpin) * self.unitCellVol)
                self.backend.addValue("dos_energies", self.energySpin + self.fermiEnergy)
                self.backend.addValue("number_of_dos_values", self.numDosVal)
                self.backend.addValue("dos_kind", "electronic")
            else:
                self.numDosVal = int(len(self.energy) / 2)
                self.totDosSpin[0] = self.totDos[0:self.numDosVal]
                self.totDosSpin[1] = self.totDos[self.numDosVal:int(2 * (self.numDosVal))]
                self.energySpin = self.energy[0:self.numDosVal]
                self.backend.addValue("dos_values", np.array(self.totDosSpin) * self.unitCellVol)
                self.backend.addValue("dos_energies", self.energySpin + self.fermiEnergy)
                self.backend.addValue("number_of_dos_values", self.numDosVal)

        elif name == 'partialdos':
            pass
        elif name == 'interstitialdos':
            self.backend.addValue("atom_projected_dos_values_lm",self.dosProjSpin)
            self.backend.addValue("number_of_lm_atom_projected_dos",len(self.dosProjSpin))
            self.backend.addValue("number_of_atom_projected_dos_values",self.numDosVal)
            self.backend.addValue("atom_projected_dos_energies",self.energy[0:self.numDosVal])
            self.backend.addValue("atom_projected_dos_m_kind","spherical")
            self.backend.closeSection("section_atom_projected_dos",self.dosProjSectionGIndex)
            self.dosProjSectionGIndex = -1
            self.inDosProj = False
    def startElementNS(self, name, qname, attrs):
        attrDict={}
        for name in attrs.getNames():
            attrDict[name] = attrs.getValue(name)
        logging.error("start element %s ns %s attr %s", name, qname, attrDict)

    def endElementNS(self, name, qname):
        logging.error("end element %s ns %s", name, qname)

    def characters(self, content):
        pass

def parseDos(inF, backend, spinTreat, unitCellVol, fermiEnergy):
    handler = DosHandler(backend, spinTreat, unitCellVol, fermiEnergy)
    logging.info("will parse")
    xml.sax.parse(inF, handler)
    logging.info("did parse")

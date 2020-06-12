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


class InputHandler(xml.sax.handler.ContentHandler):
    def __init__(self, backend, gmaxvr):
        self.gmaxvr = gmaxvr
        self.backend = backend
        self.inputSectionGIndex = -1
        self.inGWInput = False
        self.freqmax = 1.0
        self.singularity = 'mpd'
        self.actype = 'pade'
        self.npol = 0
        self.scrtype = "rpa"
        self.snempty = 0
        self.pnempty = 0
        self.coreflag = "all"
        self.fgrid = "gaule2"
        self.lmaxmb = 3
        self.epsmb = 0.0001
        self.gmb = 1.0
        self.sciavtype = "isotropic"
        self.cutofftype = "none"
        self.pwm = 2.0
        self.ngridqDum = [1, 1, 1]
        self.ngridq = [1, 1, 1]
        self.nomeg = 16

        self.freqgrid = "none"
        self.selfenergy = "none"
        self.mixbasis = "none"
        self.barecoul = "none"
        self.scrcoul = "none"

    def endDocument(self):
        if self.freqgrid == "none":
            self.backend.addValue("gw_max_frequency", self.freqmax)
            self.backend.addValue("gw_frequency_grid_type", self.fgrid)
            self.backend.addValue("gw_number_of_frequencies", self.nomeg)
        if self.selfenergy == "none":
            self.backend.addValue("gw_self_energy_c_number_of_poles", int(self.npol))
            self.backend.addValue("gw_self_energy_c_number_of_empty_states", int(self.snempty))
            self.backend.addValue("gw_self_energy_singularity_treatment", self.singularity)
            self.backend.addValue("gw_self_energy_c_analytical_continuation", self.actype)
        if self.mixbasis == "none":
            self.backend.addValue("gw_mixed_basis_lmax", self.lmaxmb)
            self.backend.addValue("gw_mixed_basis_tolerance", self.epsmb)
            self.backend.addValue("gw_mixed_basis_gmax", self.gmb*self.gmaxvr)
        if self.barecoul == "none":
            self.backend.addValue("gw_bare_coulomb_gmax", self.pwm*self.gmb*self.gmaxvr)
            self.backend.addValue("gw_bare_coulomb_cutofftype", self.cutofftype)
        if self.scrcoul == "none":
            self.backend.addValue("gw_screened_coulomb_volume_average",self.sciavtype)
            self.backend.addValue("gw_screened_Coulomb", self.scrtype)
        self.backend.addValue("gw_basis_set", "mixed")
        self.backend.addValue("gw_qp_equation_treatment", "linearization")
        for j in range(0, 3):
            # This is causing an error: Invalid literal for int()
            # with base 10.
            self.ngridq[j] = int(self.ngridqDum[j])
        self.backend.addValue("gw_ngridq", self.ngridq)

    def startElement(self, name, attrs):
        fromH = unit_conversion.convert_unit_function("hartree", "J")
        if name == "gw":
            # self.inputSectionGIndex = self.backend.openSection("section_system")
            self.inGWInput = True
            try:
                self.coreflag = attrs.getValue('coreflag')
                self.backend.addValue("gw_core_treatment", self.coreflag)
            except:
                self.coreflag = "all"
                self.backend.addValue("gw_core_treatment", self.coreflag)
            try:
                self.pnempty = attrs.getValue('nempty')
                self.backend.addValue("gw_polarizability_number_of_empty_states", int(self.pnempty))
            except:
                self.pnempty = 0
                self.backend.addValue("gw_polarizability_number_of_empty_states", int(self.pnempty))
            try:
                dummy = attrs.getValue('ngridq')
                self.ngridqDum = dummy.split()
            except:
                self.ngridqDum = [1, 1, 1]
        elif name == "freqgrid":
            self.freqgrid = "freqgrid"
            try:
                self.freqmax = attrs.getValue('freqmax')
                self.backend.addValue("gw_max_frequency", float(self.freqmax))
            except:
                self.freqmax = 1.0
                self.backend.addValue("gw_max_frequency", self.freqmax)
            try:
                self.fgrid = attrs.getValue('fgrid')
                self.backend.addValue("gw_frequency_grid_type", self.fgrid)
            except:
                self.fgrid = "gaule2"
                self.backend.addValue("gw_frequency_grid_type", self.fgrid)
            try:
                self.nomeg = attrs.getValue('nomeg')
                self.backend.addValue("gw_number_of_frequencies", int(self.nomeg))
            except:
                self.nomeg = 16
                self.backend.addValue("gw_number_of_frequencies", self.nomeg)

        elif name == "selfenergy":
            self.selfenergy = "selfenergy"
            try:
                self.npol = attrs.getValue('npol')
                self.backend.addValue("gw_self_energy_c_number_of_poles", int(self.npol))
            except:
                self.npol = 0
                self.backend.addValue("gw_self_energy_c_number_of_poles", self.npol)
            try:
                self.snempty = attrs.getValue('nempty')
                self.backend.addValue("gw_self_energy_c_number_of_empty_states", int(self.snempty))
            except:
                self.snempty = 0
                self.backend.addValue("gw_self_energy_c_number_of_empty_states", self.snempty)
            try:
                self.singularity = attrs.getValue('singularity')
                self.backend.addValue("gw_self_energy_singularity_treatment", self.singularity)
            except:
                self.singularity = 'mpd'
                self.backend.addValue("gw_self_energy_singularity_treatment", self.singularity)
            try:
                self.actype = attrs.getValue('actype')
                self.backend.addValue("gw_self_energy_c_analytical_continuation", self.actype)
            except:
                self.actype = 'pade'
                self.backend.addValue("gw_self_energy_c_analytical_continuation", self.actype)

        elif name == "mixbasis":
            self.mixbasis = "mixbasis"
            try:
                self.lmaxmb = attrs.getValue('lmaxmb')
                self.backend.addValue("gw_mixed_basis_lmax", int(self.lmaxmb))
            except:
                self.lmaxmb = 3
                self.backend.addValue("gw_mixed_basis_lmax", self.lmaxmb)
            try:
                self.epsmb = attrs.getValue('epsmb')
                self.backend.addValue("gw_mixed_basis_tolerance", float(self.epsmb))
            except:
                self.epsmb = 0.0001
                self.backend.addValue("gw_mixed_basis_tolerance", self.epsmb)
            try:
                self.gmb = attrs.getValue('gmb')
                self.backend.addValue("gw_mixed_basis_gmax", float(self.gmb)*self.gmaxvr)
            except:
                self.gmb = 1.0
                self.backend.addValue("gw_mixed_basis_gmax", self.gmb*self.gmaxvr)

        elif name == "barecoul":
            self.barecoul = "barecoul"
            try:
                self.pwm = attrs.getValue('pwm')
                self.backend.addValue("gw_bare_coulomb_gmax", float(self.pwm)*float(self.gmb)*self.gmaxvr)
            except:
                self.pwm = 2.0
                self.backend.addValue("gw_bare_coulomb_gmax", self.pwm*float(self.gmb)*self.gmaxvr)
            try:
                self.cutofftype = attrs.getValue('cutofftype')
                self.backend.addValue("gw_bare_coulomb_cutofftype", self.cutofftype)
            except:
                self.cutofftype = "none"
                self.backend.addValue("gw_bare_coulomb_cutofftype", self.cutofftype)

        elif name == "scrcoul":
            self.scrcoul = "scrcoul"
            try:
                self.sciavtype = attrs.getValue('sciavtype')
                self.backend.addValue("gw_screened_coulomb_volume_average",self.sciavtype)
            except:
                self.sciavtype = "isotropic"
                self.backend.addValue("gw_screened_coulomb_volume_average",self.sciavtype)
            try:
                self.scrtype = attrs.getValue('scrtype')
                self.backend.addValue("gw_screened_Coulomb", self.scrtype)
            except:
                self.scrtype = "rpa"
                self.backend.addValue("gw_screened_Coulomb", self.scrtype)

    def endElement(self, name):
        pass

def parseInput(inF, backend, gmaxvr):
    handler = InputHandler(backend, gmaxvr)
    logging.info("will parse")
    xml.sax.parse(inF, handler)
    logging.info("did parse")

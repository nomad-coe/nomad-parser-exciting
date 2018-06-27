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
    def __init__(self, backend, rgkmax):
        self.rgkmax = rgkmax[0]
        self.rgkmaxScr = rgkmax[0]
#        print("self.rgkmax===",self.rgkmax)
        self.backend = backend
        self.inputSectionGIndex = -1
        self.inXSInput = False
        self.fromH = unit_conversion.convert_unit_function("hartree", "J")
        self.fromB_1 = 1.88972613289*10.**10.
        self.broad = 0.01
        self.gqmax = 0.0
        self.lmaxapw = 10
        self.nemptyXS = 5
        self.nemptyScr = 0
        self.ngridqDum = [1, 1, 1]
        self.ngridkXSDum = [1, 1, 1]
        self.ngridkScrDum = [0, 0, 0]
        self.ngridq = [1,1,1]
        self.ngridkXS = [1,1,1]
        self.ngridkScr = [0,0,0]
        self.rgkmaxXs = 0.0
        self.rgkmaxscr = 0.0
        self.scissor = 0.0
        self.vkloffXS = [0.0, 0.0, 0.0]
        self.vkloffXSDum = [0.0, 0.0, 0.0]
#        self.vkloffScr = [-1.0, -1.0, -1.0]
#        self.vkloffScrDum = "-1.0 -1.0 -1.0"
        self.screening = "none"
        self.bse = "none"
        self.screentype = "full"
############ BSE variables################
        self.aresbse = "True"
        self.bsetype = "singlet"
        self.lmaxdielt = 14
        self.nstlbse = [0, 0, 0, 0]        ######## DA FARE
        self.nstlbseDum = [0, 0, 0, 0]        ######## DA FARE
        self.nstlxas = [0, 0]              ######## DA FARE
        self.nstlxasDum = [0, 0]              ######## DA FARE
        self.rgkmaxBse = rgkmax[0]
        self.sciavbd = "True"
        self.sciavqbd = "False"
        self.sciavqhd = "False"
        self.sciavqwg = "False"
        self.sciavtype = "spherical"
        self.xas = "False"
        self.xasatom = 0
        self.xasedge = "K"
        self.xasspecies = 0

#        self.xstype = "BSE"

    def endDocument(self):
#        pass
#        if self.freqgrid == "none":
#            self.backend.addValue("gw_max_frequency", self.freqmax)
#            self.backend.addValue("gw_frequency_grid_type", self.fgrid)
#            self.backend.addValue("gw_number_of_frequencies", self.nomeg)
 #       self.backend.addValue("gw_basis_set", "mixed")
#        self.backend.addValue("gw_qp_equation_treatment", "linearization")
        for j in range(0,3):

            self.ngridq[j] = int(self.ngridqDum[j])
            self.ngridkXS[j] = int(self.ngridkXSDum[j])
            self.vkloffXS[j] = float(self.vkloffXSDum[j])
            self.ngridkScr[j] = int(self.ngridkScrDum[j])
#            self.nstlbse[j] = int (self.nstlbseDum[j])
#            self.nstlxas[j] = int (self.nstlxasDum[j])
#            self.vkloffScr[j] = float(self.vkloffScrDum[j])

        for j in range(0,4):
            self.nstlbse[j] = int (self.nstlbseDum[j])
        for j in range(0,2):
            self.nstlxas[j] = int (self.nstlxasDum[j])

        self.backend.addValue("x_exciting_xs_ngridq", self.ngridq)        
        self.backend.addValue("x_exciting_xs_ngridk", self.ngridkXS)
        self.backend.addValue("x_exciting_xs_vkloff", self.vkloffXS)
        self.backend.addValue("x_exciting_xs_screening_ngridk", self.ngridkScr)
        self.backend.addValue("x_exciting_xs_bse_number_of_bands", self.nstlbse)
#        self.backend.addValue("x_exciting_xs_bse_xas_number_of_bands", self.nstlxas)

#        for j in range(0,4):
#            self.nstlbse[j] = int (self.nstlbseDum[j])
#        for j in range(0,2):
#            self.nstlxas[j] = int (self.nstlxasDum[j])

#        self.backend.addValue("x_exciting_xs_screeninig_vkloff", self.vkloffScr)

        if self.rgkmaxXs == 0.0:
            self.backend.addValue("x_exciting_xs_rgkmax", self.rgkmax)
        else:
            self.backend.addValue("x_exciting_xs_rgkmax", float(self.rgkmaxXs))

        if self.screening != "none" and self.rgkmaxScr == 0.0:
            self.backend.addValue("x_exciting_xs_screening_rgkmax", self.rgkmax)
        else:
            self.backend.addValue("x_exciting_xs_screening_rgkmax", float(self.rgkmaxScr))

        if self.bse != "none" and self.rgkmaxBse == 0.0:
            self.backend.addValue("x_exciting_xs_bse_rgkmax", self.rgkmax)
        else:
            self.backend.addValue("x_exciting_xs_bse_rgkmax", float(self.rgkmaxBse))

    def startElement(self, name, attrs):
        if name == "xs":
#            self.inputSectionGIndex = self.backend.openSection("section_system")
            self.inXSInput = True
            xstype = attrs.getValue('xstype')
            self.backend.addValue("x_exciting_xs_xstype", xstype)
            self.backend.addValue('x_exciting_electronic_structure_method', xstype)

            try:
                self.broad = attrs.getValue('broad')
                self.backend.addValue("x_exciting_xs_broadening", self.fromH(float(self.broad)))
            except:
                self.backend.addValue("x_exciting_xs_broadening", self.fromH(self.broad))
            try:
                self.gqmax = attrs.getValue('gqmax')
                self.backend.addValue("x_exciting_xs_gqmax", self.fromB_1*(float(self.gqmax)))
            except:
                self.backend.addValue("x_exciting_xs_gqmax", self.fromB_1*(self.gqmax))
            try:
                self.gqmax = attrs.getValue('lmaxapw')
                self.backend.addValue("x_exciting_xs_lmaxapw", int(self.lmaxapw))
            except:
                self.backend.addValue("x_exciting_xs_lmaxapw", self.lmaxapw)
            try:
                self.nemptyXS = attrs.getValue('nempty')
                self.backend.addValue("x_exciting_xs_number_of_empty_states", int(self.nemptyXS))
            except:
                self.backend.addValue("x_exciting_xs_number_of_empty_states", self.nemptyXS)
            try:
                dummy = attrs.getValue('ngridq')
                self.ngridqDum = dummy.split()
            except:
                self.ngridqDum = [1, 1, 1]
            try:
                dummy = attrs.getValue('ngridk')
                self.ngridkXSDum = dummy.split()
            except:
                self.ngridkXSDum = [1, 1, 1]
            try:
                self.rgkmaxXs = attrs.getValue('rgkmax')
            except:
                self.rgkmaxXs = 0.0
            try:
                self.scissor = attrs.getValue('scissor')
                self.backend.addValue("x_exciting_xs_scissor", self.fromH(float(self.scissor)))
            except:
                self.backend.addValue("x_exciting_xs_scissor", self.scissor)
            try:
                dummy = attrs.getValue('vkloff')
                self.vkloffXSDum = dummy.split()
            except:
                self.vkloffXSDum = [0.0, 0.0, 0.0]

        elif name == "screening":
            self.screening = "screening"
            try:
                self.nemptySrc = attrs.getValue('nempty')
                self.backend.addValue("x_exciting_xs_screening_number_of_empty_states", int(self.nemptyScr))
            except:
                self.backend.addValue("x_exciting_xs_screening_number_of_empty_states", self.nemptyScr)
            try:
                dummy = attrs.getValue('ngridk')
                self.ngridkScrDum = dummy.split()
#                print("dummo===",self.ngridkScrDum)
            except:
                self.ngridkScrDum = [0, 0, 0]
            try:
                self.rgkmaxScr = attrs.getValue('rgkmax')
            except:
                self.rgkmaxSrc = 0.0
            try:
                self.screentype = attrs.getValue('screentype')
                self.backend.addValue("x_exciting_xs_screening_type", int(self.screentype))
            except:
                self.backend.addValue("x_exciting_xs_screening_type", self.screentype)

#        self.aresbse = "true"
#        self.bsetype = "singlet"
#        self.lmaxdielt = 14
#        self.nstlbse = [0, 0, 0, 0]
#        self.nstlxas = [0, 0]
#        self.rgkmaxBse = rgkmax[0]
#        self.sciavbd = "true"
#        self.sciavqbd = "false"
#        self.sciavqhd = "false"
#        self.sciavqwg = "false"
#        self.sciavtype = "spherical"
#        self.xas = "false"
#        self.xasatom = 0
#        self.xasedge = "K"
#        self.xasspecies = 0

        elif name == "BSE":
            self.bse = "BSE"
            try:
                self.aresbse = attrs.getValue('aresbse')
                self.backend.addValue("x_exciting_xs_bse_antiresonant", self.aresbse)
            except:
                self.backend.addValue("x_exciting_xs_bse_antiresonant", self.aresbse)
            try:
                self.bsetype = attrs.getValue('bsetype')
                self.backend.addValue("x_exciting_xs_bse_type", self.bsetype)
            except:
                self.backend.addValue("x_exciting_xs_bse_type", self.bsetype)
            try:
                self.lmaxdielt = attrs.getValue('lmaxdielt')
                self.backend.addValue("x_exciting_xs_bse_angular_momentum_cutoff", int(self.lmaxdielt))
            except:
                self.backend.addValue("x_exciting_xs_bse_angular_momentum_cutoff", self.lmaxdielt)
            try:
                self.rgkmaxBse = attrs.getValue('rgkmax')
                self.backend.addValue("x_exciting_xs_bse_rgkmax", float(self.rgkmax))
            except:
                self.backend.addValue("x_exciting_xs_bse_rgkmax", self.rgkmax)
            try:
                self.sciavbd = attrs.getValue('sciavbd')
                self.backend.addValue("x_exciting_xs_bse_sciavbd", self.sciavbd)
            except:
                self.backend.addValue("x_exciting_xs_bse_sciavbd", self.sciavbd)
            try:
                self.sciavqbd = attrs.getValue('sciavqbd')
                self.backend.addValue("x_exciting_xs_bse_sciavqbd", self.sciavqbd)
            except:
                self.backend.addValue("x_exciting_xs_bse_sciavqbd", self.sciavqbd)
            try:
                self.sciavqhd = attrs.getValue('sciavqhd')
                self.backend.addValue("x_exciting_xs_bse_sciavqhd", self.sciavqhd)
            except:
                self.backend.addValue("x_exciting_xs_bse_sciavqhd", self.sciavqhd)
            try:
                self.sciavqwg = attrs.getValue('sciavqwg')
                self.backend.addValue("x_exciting_xs_bse_sciavqwg", self.sciavqwg)
            except:
                self.backend.addValue("x_exciting_xs_bse_sciavqwg", self.sciavqwg)
            try:
                self.sciavtype = attrs.getValue('sciavtype')
                self.backend.addValue("x_exciting_xs_bse_sciavtype", self.sciavtype)
            except:
                self.backend.addValue("x_exciting_xs_bse_sciavtype", self.sciavtype)
            try:
                self.xas = attrs.getValue('xas')
                self.backend.addValue("x_exciting_xs_bse_xas", self.xas)
            except:
                self.backend.addValue("x_exciting_xs_bse_xas", self.xas)
            try:
                self.xasatom = attrs.getValue('xasatom')
                self.backend.addValue("x_exciting_xs_bse_xasatom", self.xasatom)
            except:
                self.backend.addValue("x_exciting_xs_bse_xasatom", self.xasatom)
            try:
                self.xasedge = attrs.getValue('xasedge')
                self.backend.addValue("x_exciting_xs_bse_xasedge", self.xasedge)
            except:
                self.backend.addValue("x_exciting_xs_bse_xasedge", self.xasedge)
            try:
                self.xasspecies = attrs.getValue('xasspecies')
                self.backend.addValue("x_exciting_xs_bse_xasspecies", self.xasspecies)
            except:
                self.backend.addValue("x_exciting_xs_bse_xasspecies", self.xasspecies)
            try:
                dummy = attrs.getValue('nstlbse')
                self.nstlbseDum = dummy.split()
            except:
                self.nstlbseDum = [0, 0, 0, 0]
            try:
                dummy = attrs.getValue('nstlxas')
                self.nstlxasDum = dummy.split()
            except:
                self.nstlxasDum = [0, 0]
#            try:
#                dummy = attrs.getValue('vkloff')
#                self.vkloffSrcDum = dummy.split()lmaxdielt
#            except:
#                self.vkloffSrcDum = "0.0 0.0 0.0"

#            try:
#                self.fgrid = attrs.getValue('fgrid')
#                self.backend.addValue("gw_frequency_grid_type", self.fgrid)
#            except:
#                self.fgrid = "gaule2"
#                self.backend.addValue("gw_frequency_grid_type", self.fgrid)
#            try:
#                self.nomeg = attrs.getValue('nomeg')
#                self.backend.addValue("gw_number_of_frequencies", int(self.nomeg))
#            except:
#                self.nomeg = 16
#                self.backend.addValue("gw_number_of_frequencies", self.nomeg)
#
#        elif name == "selfenergy":
#            self.selfenergy = "selfenergy"
#            try:
 #               self.npol = attrs.getValue('npol')
 #               self.backend.addValue("gw_self_energy_c_number_of_poles", int(self.npol))
 #           except:
 #               self.npol = 0
 #               self.backend.addValue("gw_self_energy_c_number_of_poles", self.npol)
 #           try:
 #               self.snempty = attrs.getValue('nempty')
 #               self.backend.addValue("gw_self_energy_c_number_of_empty_states", int(self.snempty))
 #           except:
 #               self.snempty = 0
 #               self.backend.addValue("gw_self_energy_c_number_of_empty_states", self.snempty)
 #           try:
 #               self.singularity = attrs.getValue('singularity')
#                self.backend.addValue("gw_self_energy_singularity_treatment", self.singularity)
#            except:
#                self.singularity = 'mpd'
#                self.backend.addValue("gw_self_energy_singularity_treatment", self.singularity)
#            try:
#                self.actype = attrs.getValue('actype')
#                self.backend.addValue("gw_self_energy_c_analytical_continuation", self.actype)
#            except:
 #               self.actype = 'pade'
 #               self.backend.addValue("gw_self_energy_c_analytical_continuation", self.actype)
#
#        elif name == "mixbasis":
##            self.mixbasis = "mixbasis"
#            try:
#                self.lmaxmb = attrs.getValue('lmaxmb')
#                self.backend.addValue("gw_mixed_basis_lmax", int(self.lmaxmb))
#            except:
#                self.lmaxmb = 3
#                self.backend.addValue("gw_mixed_basis_lmax", self.lmaxmb)
#            try:
#                self.epsmb = attrs.getValue('epsmb')
#                self.backend.addValue("gw_mixed_basis_tolerance", float(self.epsmb))
#            except:
#                self.epsmb = 0.0001
#                self.backend.addValue("gw_mixed_basis_tolerance", self.epsmb)
#            try:
#                self.gmb = attrs.getValue('gmb')
#                self.backend.addValue("gw_mixed_basis_gmax", float(self.gmb)*self.gmaxvr)
#            except:
#                self.gmb = 1.0
#                self.backend.addValue("gw_mixed_basis_gmax", self.gmb*self.gmaxvr)
#
#        elif name == "barecoul":
#            self.barecoul = "barecoul"
#            try:
#                self.pwm = attrs.getValue('pwm')
#                self.backend.addValue("gw_bare_coulomb_gmax", float(self.pwm)*float(self.gmb)*self.gmaxvr)
#            except:
#                self.pwm = 2.0
#                self.backend.addValue("gw_bare_coulomb_gmax", self.pwm*float(self.gmb)*self.gmaxvr)
#            try:
#                self.cutofftype = attrs.getValue('cutofftype')
#                self.backend.addValue("gw_bare_coulomb_cutofftype", self.cutofftype)
#            except:
#                self.cutofftype = "none"
 #               self.backend.addValue("gw_bare_coulomb_cutofftype", self.cutofftype)
#
#        elif name == "scrcoul":
#            self.scrcoul = "scrcoul"
#            try:
#                self.sciavtype = attrs.getValue('sciavtype')
#                self.backend.addValue("gw_screened_coulomb_volume_average",self.sciavtype)
#            except:
#                self.sciavtype = "isotropic"
#                self.backend.addValue("gw_screened_coulomb_volume_average",self.sciavtype)
#            try:
#                self.scrtype = attrs.getValue('scrtype')
#                self.backend.addValue("gw_screened_Coulomb", self.scrtype)
#            except:
#                self.scrtype = "rpa"
#                self.backend.addValue("gw_screened_Coulomb", self.scrtype)

    def endElement(self, name):
        pass

def parseInput(inF, backend, gmaxvr):
    handler = InputHandler(backend, gmaxvr)
    logging.error("will parse")
    xml.sax.parse(inF, handler)
    logging.error("did parse")

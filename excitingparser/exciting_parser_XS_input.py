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
    def __init__(self, backend, rgkmax):
        # self.xsType = None
        self.rgkmax = rgkmax[0]
        self.rgkmaxScr = rgkmax[0]
        # print("self.rgkmax===",self.rgkmax)
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
        # self.vkloffScr = [-1.0, -1.0, -1.0]
        # self.vkloffScrDum = "-1.0 -1.0 -1.0"
        self.screening = "none"
        self.bse = False
        # self.screentype = "full"

############ BSE variables ################
        self.aresbse = True
        # self.bsetype = "singlet"
        self.lmaxdielt = 14
        self.nstlbse = [0, 0, 0, 0]
        self.nstlbseDum = [0, 0, 0, 0]
        self.nstlxas = [0, 0]
        self.nstlxasDum = [0, 0]
        self.rgkmaxBse = rgkmax[0]
        self.sciavbd = True
        self.sciavqbd = False
        self.sciavqhd = False
        self.sciavqwg = False
        self.sciavtype = "spherical"
        self.xas = False
        self.xasatom = 0
        self.xasedge = "K"
        self.xasspecies = 0
        self.tddft = False
        self.acon = "false"
        self.ahc = "false"
        self.drude = [0.0, 0.0]
        self.fxcbsesplit = 0.00001
        self.fxctype = "RPA"
        self.lmaxalda = 3
        self.nwacont = 0
        self.tetra = False
        self.tetradf = "false"
        # self.xstype = "BSE"

    def endDocument(self):
        #        pass
        #if self.freqgrid == "none":
        #    self.backend.addValue("gw_max_frequency", self.freqmax)
        #    self.backend.addValue("gw_frequency_grid_type", self.fgrid)
        #    self.backend.addValue("gw_number_of_frequencies", self.nomeg)
        #       self.backend.addValue("gw_basis_set", "mixed")
        #self.backend.addValue("gw_qp_equation_treatment", "linearization")

        if self.tetradf == "true":
            self.backend.addValue("x_exciting_xs_tetra", True)
        else:
            self.backend.addValue("x_exciting_xs_tetra", False)

        for j in range(0,3):

            self.ngridq[j] = int(self.ngridqDum[j])
            self.ngridkXS[j] = int(self.ngridkXSDum[j])
            self.vkloffXS[j] = float(self.vkloffXSDum[j])
            self.ngridkScr[j] = int(self.ngridkScrDum[j])
            # self.nstlbse[j] = int (self.nstlbseDum[j])
            # self.nstlxas[j] = int (self.nstlxasDum[j])
            # self.vkloffScr[j] = float(self.vkloffScrDum[j])

        for j in range(0,4):
            self.nstlbse[j] = int (self.nstlbseDum[j])
        for j in range(0,2):
            self.nstlxas[j] = int (self.nstlxasDum[j])

        self.backend.addValue("x_exciting_xs_ngridq", self.ngridq)
        self.backend.addValue("x_exciting_xs_ngridk", self.ngridkXS)
        self.backend.addValue("x_exciting_xs_vkloff", self.vkloffXS)
        self.backend.addValue("x_exciting_xs_screening_ngridk", self.ngridkScr)
        if self.bse == True:
            self.backend.addValue("x_exciting_xs_bse_number_of_bands", self.nstlbse)
        elif self.xas == True:
            self.backend.addValue("x_exciting_xs_bse_xasatom", self.xasatom)
            self.backend.addValue("x_exciting_xs_bse_xasedge", self.xasedge)
            self.backend.addValue("x_exciting_xs_bse_xasspecies", self.xasspecies)
            self.backend.addValue("x_exciting_xs_bse_xas_number_of_bands", self.nstlxas)
        # self.backend.addValue("x_exciting_xs_bse_xas_number_of_bands", self.nstlxas)

        # for j in range(0,4):
        #    self.nstlbse[j] = int (self.nstlbseDum[j])
        # for j in range(0,2):
        #   self.nstlxas[j] = int (self.nstlxasDum[j])
        #   self.backend.addValue("x_exciting_xs_screeninig_vkloff", self.vkloffScr)

        if self.rgkmaxXs == 0.0:
            self.backend.addValue("x_exciting_xs_rgkmax", self.rgkmax)
        else:
            self.backend.addValue("x_exciting_xs_rgkmax", float(self.rgkmaxXs))

        if self.screening != "none" and self.rgkmaxScr == 0.0:
            self.backend.addValue("x_exciting_xs_screening_rgkmax", self.rgkmax)
        else:
            self.backend.addValue("x_exciting_xs_screening_rgkmax", float(self.rgkmaxScr))

        if self.bse == True:
            if self.rgkmaxBse == 0.0:
                self.backend.addValue("x_exciting_xs_bse_rgkmax", self.rgkmax)
            else:
                self.backend.addValue("x_exciting_xs_bse_rgkmax", float(self.rgkmaxBse))

    def startElement(self, name, attrs):
        # xsType = None
        if name == "xs":
            # self.inputSectionGIndex = self.backend.openSection("section_system")
            self.inXSInput = True
            xstype = attrs.getValue('xstype')
            self.backend.addValue("x_exciting_xs_xstype", xstype)
            self.backend.addValue('x_exciting_electronic_structure_method', xstype)
            # print("xstyppe===",xstype)

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
                # print("dummo===",self.ngridkScrDum)
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

            # self.aresbse = "true"
            # self.bsetype = "singlet"
            # self.lmaxdielt = 14
            # self.nstlbse = [0, 0, 0, 0]
            # self.nstlxas = [0, 0]
            # self.rgkmaxBse = rgkmax[0]
            # self.sciavbd = "true"
            # self.sciavqbd = "false"
            # self.sciavqhd = "false"
            # self.sciavqwg = "false"
            # self.sciavtype = "spherical"
            # self.xas = "false"
            # self.xasatom = 0
            # self.xasedge = "K"
            # self.xasspecies = 0

        elif name == "BSE":
            self.bse = True
            # xsType = "BSE"
            try:
                self.aresbse = attrs.getValue('aresbse')
                if self.aresbse == "true":
                    self.backend.addValue("x_exciting_xs_bse_antiresonant", True)
                else:
                    self.backend.addValue("x_exciting_xs_bse_antiresonant", False)
            except:
                self.backend.addValue("x_exciting_xs_bse_antiresonant", True)
            try:
                self.lmaxdielt = attrs.getValue('lmaxdielt')
                self.backend.addValue("x_exciting_xs_bse_angular_momentum_cutoff", int(self.lmaxdielt))
            except:
                self.backend.addValue("x_exciting_xs_bse_angular_momentum_cutoff", self.lmaxdielt)
            try:
                self.rgkmaxBse = attrs.getValue('rgkmax')
                # self.backend.addValue("x_exciting_xs_bse_rgkmax", float(self.rgkmax))
            except:
                pass
                # self.backend.addValue("x_exciting_xs_bse_rgkmax", self.rgkmax)
            try:
                self.sciavbd = attrs.getValue('sciavbd')
                if self.sciavqbd == "true":
                    self.backend.addValue("x_exciting_xs_bse_sciavbd", True)
                else:
                    self.backend.addValue("x_exciting_xs_bse_sciavbd", False)
            except:
                self.backend.addValue("x_exciting_xs_bse_sciavbd", True)
            try:
                self.sciavqbd = attrs.getValue('sciavqbd')
                if self.sciavqbd == "true":
                    self.backend.addValue("x_exciting_xs_bse_sciavqbd", True)
                else:
                    self.backend.addValue("x_exciting_xs_bse_sciavqbd", False)
            except:
                self.backend.addValue("x_exciting_xs_bse_sciavqbd", False)
            try:
                self.sciavqhd = attrs.getValue('sciavqhd')
                if self.sciavqhd == "true":
                    self.backend.addValue("x_exciting_xs_bse_sciavqhd", True)
                else:
                    self.backend.addValue("x_exciting_xs_bse_sciavqhd", False)
            except:
                self.backend.addValue("x_exciting_xs_bse_sciavqhd", False)
            try:
                self.sciavqwg = attrs.getValue('sciavqwg')
                if self.sciavqwg == "true":
                    self.backend.addValue("x_exciting_xs_bse_sciavqwg", True)
                else:
                    self.backend.addValue("x_exciting_xs_bse_sciavqwg", False)
            except:
                self.backend.addValue("x_exciting_xs_bse_sciavqwg", False)
            try:
                self.sciavtype = attrs.getValue('sciavtype')
                self.backend.addValue("x_exciting_xs_bse_sciavtype", self.sciavtype)
            except:
                self.backend.addValue("x_exciting_xs_bse_sciavtype", self.sciavtype)
            try:
                self.xas = attrs.getValue('xas')
                if self.xas == "false":
                    # print("xas===",self.xas)
                    self.backend.addValue("x_exciting_xs_bse_xas", False)
                else:
                    # print("xas===",self.xas)
                    self.backend.addValue("x_exciting_xs_bse_xas", True)
            except:
                # print("xasdefault===",self.xas)
                self.backend.addValue("x_exciting_xs_bse_xas", self.xas)
            try:
                self.xasatom = int(attrs.getValue('xasatom'))
                # self.backend.addValue("x_exciting_xs_bse_xasatom", self.xasatom)
            except:
                pass
                # self.backend.addValue("x_exciting_xs_bse_xasatom", self.xasatom)
            try:
                self.xasedge = attrs.getValue('xasedge')
                # self.backend.addValue("x_exciting_xs_bse_xasedge", self.xasedge)
            except:
                pass
                # self.backend.addValue("x_exciting_xs_bse_xasedge", self.xasedge)
            try:
                self.xasspecies = int(attrs.getValue('xasspecies'))
                # self.backend.addValue("x_exciting_xs_bse_xasspecies", self.xasspecies)
            except:
                pass
                # self.backend.addValue("x_exciting_xs_bse_xasspecies", self.xasspecies)
                # if self.xas == False:
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

        elif name == "tddft":
            self.tddft = True

            try:
                self.acon = attrs.getValue('acont')
                if self.acon == "true":
                    self.backend.addValue("x_exciting_xs_tddft_analytic_continuation", True)
                else:
                    self.backend.addValue("x_exciting_xs_tddft_analytic_continuation", False)
            except:
                    self.backend.addValue("x_exciting_xs_tddft_analytic_continuation", False)
            try:
                self.ahc = attrs.getValue('ahc')
                if self.ahc == "true":
                    self.backend.addValue("x_exciting_xs_tddft_anomalous_Hall_conductivity", True)
                else:
                    self.backend.addValue("x_exciting_xs_tddft_anomalous_Hall_conductivity", False)
            except:
                    self.backend.addValue("x_exciting_xs_tddft_anomalous_Hall_conductivity", False)
            try:
                self.aresdf = attrs.getValue('aresdf')
                if self.aresdf == "true":
                    self.backend.addValue("x_exciting_xs_tddft_anti_resonant_dielectric", True)
                elif self.aresdf == "false":
                    self.backend.addValue("x_exciting_xs_tddft_anti_resonant_dielectric", False)
            except:
                    self.backend.addValue("x_exciting_xs_tddft_anti_resonant_dielectric", True)
            try:
                self.aresfxc = attrs.getValue('aresfxc')
                if self.aresfxc == "true":
                    self.backend.addValue("x_exciting_xs_tddft_anti_resonant_xc_kernel", True)
                elif self.aresfxc == "false":
                    self.backend.addValue("x_exciting_xs_tddft_anti_resonant_xc_kernel", False)
            except:
                    self.backend.addValue("x_exciting_xs_tddft_anti_resonant_xc_kernel", True)
            try:
                dummy = attrs.getValue('drude')
                self.drude = dummy.split()
                for j in range (0,2):
                    self.drude[j] = float(self.drude[j])
                self.backend.addValue("x_exciting_xs_tddft_drude", self.drude)
            except:
                self.backend.addValue("x_exciting_xs_tddft_drude", self.drude)
            try:
                self.fxcbsesplit = attrs.getValue('fxcbsesplit')
                self.backend.addValue("x_exciting_xs_tddft_split_parameter", self.fromH(float(self.fxcbsesplit)))
            except:
                self.backend.addValue("x_exciting_xs_tddft_split_parameter", self.fromH(self.fxcbsesplit))
            try:
                self.fxctype = attrs.getValue('fxctype')
                self.backend.addValue("x_exciting_xs_tddft_xc_kernel", self.fxctype)
            except:
                self.backend.addValue("x_exciting_xs_tddft_xc_kernel", self.fxctype)
            try:
                intraband = attrs.getValue('intraband')
                if intraband == "true":
                    self.backend.addValue("x_exciting_xs_tddft_finite_q_intraband_contribution", True)
                else:
                    self.backend.addValue("x_exciting_xs_tddft_finite_q_intraband_contribution", False)
            except:
                    self.backend.addValue("x_exciting_xs_tddft_finite_q_intraband_contribution", False)
            try:
                kerndiag = attrs.getValue('kerndiag')
                if kerndiag == "true":
                    self.backend.addValue("x_exciting_xs_tddft_diagonal_xc_kernel", True)
                else:
                    self.backend.addValue("x_exciting_xs_tddft_diagonal_xc_kernel", False)
            except:
                    self.backend.addValue("x_exciting_xs_tddft_diagonal_xc_kernel", False)
            try:
                self.lmaxalda = attrs.getValue('lmaxalda')
                self.backend.addValue("x_exciting_xs_tddft_lmax_alda", int(self.lmaxalda))
            except:
                self.backend.addValue("x_exciting_xs_tddft_lmax_alda", self.lmaxalda)
            try:
                mdfqtype = attrs.getValue('mdfqtype')
                if mdfqtype == "0":
                    self.backend.addValue("x_exciting_xs_tddft_macroscopic_dielectric_function_q_treatment", int(mdfqtype))
                else:
                    self.backend.addValue("x_exciting_xs_tddft_macroscopic_dielectric_function_q_treatment", int(mdfqtype))
            except:
                    mdfqtype = 0
                    self.backend.addValue("x_exciting_xs_tddft_macroscopic_dielectric_function_q_treatment", int(mdfqtype))
            try:
                self.nwacont = attrs.getValue('nwacont')
                self.backend.addValue("x_exciting_xs_tddft_analytic_continuation_number_of_intervals", int(self.nwacont))
            except:
                self.backend.addValue("x_exciting_xs_tddft_analytic_continuation_number_of_intervals", self.nwacont)

        elif name == "tetra":
            self.tetra = True
            try:
                self.tetradf = attrs.getValue('tetradf')
                # if self.tetradf == "true":
                #   selfbackend.addValue("x_exciting_xs_tetra", True)
                # else:
                #   self.backend.addValue("x_exciting_xs_tetra", False)
            except:
                self.tetradf == "false"
                # self.backend.addValue("x_exciting_xs_tetra", False)

    def endElement(self, name):
        pass

def parseInput(inF, backend, gmaxvr):
    handler = InputHandler(backend, gmaxvr)
    logging.info("will parse")
    xml.sax.parse(inF, handler)
    logging.info("did parse")

# Copyright 2016-2018 The NOMAD Developers Group
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
#
# Maintainers:
# 2017-2018: Lorenzo Pardini <loren.pard@gmail.com>
# 2019-2020: Cuauhtemoc Salazar

from builtins import object
import numpy as np
from nomadcore.simple_parser import AncillaryParser, CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM, mainFunction
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.caching_backend import CachingLevel
from nomadcore.unit_conversion import unit_conversion
from nomadcore.unit_conversion.unit_conversion import convert_unit_function
import os, sys, json
import excitingparser.exciting_parser_dos as exciting_parser_dos
import excitingparser.exciting_parser_bandstructure as exciting_parser_bandstructure
import excitingparser.exciting_parser_gw as exciting_parser_gw
import excitingparser.exciting_parser_GS_input as exciting_parser_GS_input
import excitingparser.exciting_parser_XS_input as exciting_parser_XS_input
import excitingparser.exciting_parser_xs as exciting_parser_xs
import excitingparser.exciting_parser_eps as exciting_parser_eps
from ase import Atoms
import logging

logger = logging.getLogger("nomad.ExcitingParser")

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class ExcitingParserContext(object):

  def __init__(self):
    self.parser = None
    # self.mainFileUri = sys.argv[1]    #exciting !!!!!!LOCAL HOME!!!!!!!!             OKOKOKOK
    self.mainFileUri = sys.argv[2]  #exciting !!! FOR NOMAD URI nmd:// or sbt -> zip file!!!!!!!!   OKOKKOOK
    self.mainFilePath = None
    self.mainFile = None
    self.volumeCubOpt = False
    self.volumeOpt = False
    self.volumeCubOptIndex = 0

  def initialize_values(self):
    self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv
    self.atom_pos = []
    self.atom_labels = []

  def startedParsing(self, path, parser):
    self.parser=parser
    self.initialize_values()
    # self.atom_pos = []
    # self.atom_labels = []
    self.XSSetGIndex = None
    self.secMethodIndex = None
    self.secSystemIndex = None
    self.secSingleConfIndex = None
    self.spinTreat = None
    self.sim_cell = []
    self.cell_format = ''
    self.secRunIndex = None
    self.unit_cell_vol = 0
    self.xcName = None
    self.gmaxvr = 0
    self.rgkmax = 0
    self.energy_thresh = []
    self.samplingMethod = None
    self.secSamplingMethodIndex = None
    self.geometryForceThreshold = 0
    self.frameSequence = []
    self.samplingGIndex = 0
    self.dummy = 0
    self.i=0
    self.clathrates = False
    self.bsetype = None
    self.screentype = None
    self.tensorComp = []
    self.excitonEnergies = []
    self.xsTetra = False
    self.xsAC = False
    self.xsNAR = False
    self.xstype = None
    self.tddftKernel = None
    # self.xsType = None
    # self.volumeOptIndex = 0

  def onOpen_section_run(self, backend, gIndex, section):
    #  self.i=self.i+1
    #  print("self.i=",self.i)
    curDir = os.getcwd()
    mainFile = self.parser.fIn.fIn.name           #####exciting sbt -> zip filei r from NOMAD URI nmd:// ####     YES ?????????? sure???? check first    OKOKOKKO
    self.mainFilePath = os.path.dirname(mainFile)           #####exciting sbt -> zip file####     YES ?????????? sure???? check first    OKOKOKKO
    #  self.mainFilePath = self.mainFileUri[0:-9]     #####exciting LOCAL HOME #######   YES                      OKOKOKOK
    #  print("self.mainFilePath===",self.mainFilePath)
    os.chdir(self.mainFilePath)
    os.chdir('../')
    #    self.i=self.i+1
    #    print("self.i=",self.i)
    #    self.volumeOptIndex = 0
    #    for root, dirs, files in os.walk('.'):
    #      if root == '.':
    #        for directory in dirs:
    #          if directory[-4:-1]==
    #        print("root=",root)
    #        print("dirs=",dirs)
    #        print("files=",files)

    #####volume optimization for simple cubic systems########
    if 'INFO_VOL' in os.listdir('.'):
      self.volumeOpt = True
      #  print("self.volumeOpt=",self.volumeOpt)
      with open('INFO_VOL') as g:
          while 1:
            s = g.readline()
            if not s: break
            s = s.strip()
    else:
      for files in os.listdir('.'):
      #  print("files=",files)
      #  backend.addValue("x_exciting_dummy2", files)
        if files[0:7] == 'rundir-':
          self.volumeCubOptIndex+=1
    os.chdir(curDir)
    #  curDir = os.getcwd()
    #  dirPath = self.mainFileUri[0:-9]
    #  os.chdir(dirPath)
    #  os.chdir('../')
    #  i = 0
    #  for files in os.listdir('.'):
    #    if files[0:7] == 'rundir-':
    #      i+=1
    #    if i>1:
    #      self.volumeOpt = True
    #      optGindex = backend.openSection("section_method")
    #      print("optGindex=",optGindex)
    #      backend.addValue("x_exciting_volume_optimization", self.volumeOpt)
    #      backend.closeSection("section_method", optGindex)
    #  os.chdir(curDir)

  def onOpen_section_sampling_method(self, backend, gIndex, section):
    self.secSamplingMethodIndex = gIndex
    backend.addValue("sampling_method", "geometry_optimization")
    #  print("self.secSamplingMethodIndex=",self.secSamplingMethodIndex)
    self.samplingMethod = "geometry_optimization"
    #  print("self.samplingMethod=",self.samplingMethod)

  def onOpen_section_system(self, backend, gIndex, section):
    self.secSystemIndex = gIndex
    curDir = os.getcwd()
    mainFile = self.parser.fIn.fIn.name           #####exciting sbt -> zip file or from NOMAD URI nmd:// ####     YES ?????????? sure???? check first    OKOKOKKO
    self.mainFilePath = os.path.dirname(mainFile)           #####exciting sbt -> zip file####     YES ?????????? sure???? check first    OKOKOKKO
    #  self.mainFilePath = self.mainFileUri[0:-9]     #####exciting LOCAL HOME#######   YES                      OKOKOKOK
    #  print("self.mainFilePath===",self.mainFilePath)
    os.chdir(self.mainFilePath)
    #  print("listdir=",os.listdir('.'))
    if 'str.out' in os.listdir('.'):
      self.clathrates = True
      #  print("clathrate_vero===",self.clathrates)
      backend.addValue("x_exciting_clathrates", True)
      clathrate_labels = []
      clathrate_coordinates = []
      with open('str.out') as g:
          while 1:
            s = g.readline()
            if not s: break
            s = s.strip()
            s = s.split()
            if len(s) == 4:
              #  print("s=",s[0],float(s[0]))
              clathrate_coordinates.append([float(s[0]),float(s[1]),float(s[2])])
              clathrate_labels.append(s[3])
      #  print("clathrate_coordinates=",clathrate_coordinates)
      backend.addArrayValues("x_exciting_clathrates_atom_coordinates", np.asarray(clathrate_coordinates))
      backend.addValue("x_exciting_clathrates_atom_labels", clathrate_labels)
    else:
      #  print("clathrate_falso===",self.clathrates)
      backend.addValue("x_exciting_clathrates", False)
      #  backend.addArrayValues("x_exciting_clathrates_atom_coordinates", np.array(clathrate_coordinates))
      #  backend.addValue("x_exciting_clathrates_atom_labels", clathrate_labels)
    os.chdir(curDir)

  def onOpen_section_single_configuration_calculation(self, backend, gIndex, section):
    if self.secSingleConfIndex is None:
      self.secSingleConfIndex = gIndex
    self.frameSequence.append(gIndex)

  def onOpen_section_method(self, backend, gIndex, section):
    if self.secMethodIndex is None:
      self.secMethodIndex = gIndex

  def onOpen_x_exciting_section_geometry_optimization(self, backend, gIndex, section):
    #    """Trigger called when x_abinit_section_dataset is opened.
    #    """
    self.samplingGIndex = backend.openSection("section_sampling_method")

  def onClose_section_run(self, backend, gIndex, section):
    self.secRunIndex = gIndex
    #  backend.addValue("x_exciting_dummy", self.volumeOptIndex)
    ########### VOLUME OPTIMIZATION #########################
    #    curDir = os.getcwd()
    ##    mainFile = self.parser.fIn.fIn.name           #####exciting sbt -> zip file####     YES ?????????? sure???? check first    OKOKOKKO
    ##    dirPath = os.path.dirname(mainFile)           #####exciting sbt -> zip file####     YES ?????????? sure???? check first    OKOKOKKO
    #    dirPath = self.mainFileUri[0:-9]     #####exciting LOCAL HOME or from NOMAD URI nmd://  #######   YES                      OKOKOKOK
    #    os.chdir(dirPath)
    #    os.chdir('../')
    #    self.volumeOptIndex = 0
    #    for files in os.listdir('.'):
    #      if files[0:7] == 'rundir-':
    #        self.volumeOptIndex+=1
    ######################independent from above#################
    if self.volumeCubOptIndex>1:
      self.volumeCubOpt = True
      optGindex = backend.openSection("section_method")
      backend.addValue("x_exciting_volume_optimization", self.volumeCubOpt)
      backend.closeSection("section_method", optGindex)
    #  os.chdir(curDir)

    ####################################TEST############################
    mainFile = self.parser.fIn.fIn.name
    dirPath = os.path.dirname(self.parser.fIn.name)
    #  print("dirPath===",dirPath)
    gw_File = os.path.join(dirPath, "GW_INFO.OUT")
    gwFile = os.path.join(dirPath, "GWINFO.OUT")
    xsFile = os.path.join(dirPath, "INFOXS.OUT")
    for gFile in [gw_File, gwFile]:
      if os.path.exists(gFile):
        gwParser = exciting_parser_gw.GWParser()
        gwParser.parseGW(gFile, backend,
                         dftMethodSectionGindex = self.secMethodIndex,
                         dftSingleConfigurationGindex = self.secSingleConfIndex,
                         xcName = self.xcName,
                         unitCellVol = self.unit_cell_vol,
                         gmaxvr = self.gmaxvr)

        subParser = AncillaryParser(
            fileDescription = exciting_parser_gw.buildGWMatchers(),
            parser = self.parser,
            cachingLevelForMetaName = exciting_parser_gw.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.PreOpenedIgnore),
            superContext = gwParser)
        with open(gFile) as fIn:
            subParser.parseFile(fIn)
        break

    if os.path.exists(xsFile):
      excNum = []
      excEn = []
      excBindEn = []
      osclStr = []
      transCoeff = []
      epsEn = []
      sigmaEn = []
      epsilon = []
      sigma = []
      lossEn = []
      loss = []
      lossDummy = []
      xstype = "TDDFT"
      self.xstype = "TDDFT"
      fromeV = unit_conversion.convert_unit_function("eV", "J")
      #  print("BSE!!!")
      #  if os.path.exists(inputgwFile):
      inputXSFile = os.path.join(dirPath, "input.xml")
      self.XSSetGIndex = backend.openSection("section_method")
      try:  # tmk:
        with open(inputXSFile) as f:
          exciting_parser_XS_input.parseInput(f, backend, self.rgkmax)
      except FileNotFoundError:
        logger.warning("File not found: {}" .format(inputXSFile))
      except Exception as err:
        logger.error("Exception on {}" .format(__file__), exc_info=err)
        #    xstype = section["x_exciting_xs_type"]
        #    print("xstype===",xstype)
        #    print("xsType===",exciting_parser_XS_input.InputHandler.self.xsType)
        #  backend.addValue('x_exciting_electronic_structure_method', "BSE")   ########!!!!!!!!!!! So far, BSE is not a valid value for electronic_structure_method. It must be added! #########
        #    teste = section["x_exciting_xs_tetra"]
        #    print("teste===",teste)
        #  backend.addValue('x_exciting_xs_starting_point', self.xcName)
        #  teste = section["x_exciting_xs_tetra"]
        #  print("teste===",teste)
      if self.secMethodIndex is not None:
        m2mGindex = backend.openNonOverlappingSection("section_method_to_method_refs")
        backend.addValue("method_to_method_ref", self.secMethodIndex)
        backend.addValue("method_to_method_kind", "starting_point")
        backend.closeNonOverlappingSection("section_method_to_method_refs")

        #  xstype = section["x_exciting_xs_type"]
        #  print("xstype===",xstype)
        for files in os.listdir(dirPath):
            if files[0:11] == "EXCITON_BSE":
                xstype = "BSE"
                self.xstype = "BSE"
                dummyBse = files[11:13]
                #  self.tensorComp = files[-6:-4]
                self.tensorComp.append(files[-6:-4])

                if dummyBse == 'si':
                    self.bsetype = 'singlet'
                    #  name = "EXCITON_BSE" + self.bsetype
                elif dummyBse == 'tr':
                    self.bsetype = 'triplet'
                elif dummyBse == 'RP':
                    self.bsetype = 'RPA'
                elif dummyBse == 'IP':
                    self.bsetype = 'IP'
                name = "EXCITON_BSE" + self.bsetype + '_SCR'
                #  nameEps = "EPSILON_BSE" + self.bsetype + '_SCR'
                if files[len(name):len(name)+4] == 'full':
                    self.screentype = 'full'
                elif files[len(name):len(name)+4] == 'diag':
                    self.screentype = 'diag'
                elif files[len(name):len(name)+4] == 'noin':
                    self.screentype = 'noinvdiag'
                elif files[len(name):len(name)+4] == 'long':
                    self.screentype = 'longrange'
          #  else:
          #    xstype = "TDDFT"
          #  print("xstype===",xstype)
          #  teste = section["x_exciting_xs_tetra"]
          #  print("teste===",teste)
        if xstype == "BSE":
            numberOfComponents = len(self.tensorComp)
            #  backend.addValue("x_exciting_xs_screening_type", self.screentype)
            backend.addValue("x_exciting_xs_bse_type", self.bsetype)
      #  print("===", self.tensorComp[0])
      #  with open(inputXSFile) as f:
      #    exciting_parser_XS_input.parseInput(f, backend, self.rgkmax)
      #    teste = section["x_exciting_xs_tetra"]
      #    print("teste===",teste)
      #  teste = section["x_exciting_xs_tetra"]
      #  print("teste===",teste)
      backend.closeSection("section_method",self.XSSetGIndex)
      #  teste = section["x_exciting_xs_tetra"]
      #  print("teste===",teste)
      #  numberOfComponents = len(self.tensorComp)
      #  print("xstype=====",xstype)
      backend.openNonOverlappingSection("section_single_configuration_calculation")
      if self.secSingleConfIndex is not None:
          backend.openNonOverlappingSection("section_calculation_to_calculation_refs")
          backend.addValue("calculation_to_calculation_ref", self.secSingleConfIndex)
          backend.addValue("calculation_to_calculation_kind", "starting_point")
          backend.closeNonOverlappingSection("section_calculation_to_calculation_refs")

      if xstype == "BSE":
          for i in range(numberOfComponents):
              excNum.append([])
              excEn.append([])
              epsEn.append([])
              sigmaEn.append([])
              excBindEn.append([])
              osclStr.append([])
              transCoeff.append([[],[]])
              epsilon.append([[],[]])
              sigma.append([[],[]])
              lossEn.append([])
              loss.append([[],[]])
            #  lossDummy.append([])

              outputXSFile = os.path.join(dirPath, "EXCITON_BSE" + self.bsetype + '_SCR' + self.screentype + "_OC" + self.tensorComp[i] + ".OUT")
              outputEpsFile = os.path.join(dirPath, "EPSILON_BSE" + self.bsetype + '_SCR' + self.screentype + "_OC" + self.tensorComp[i] + ".OUT")
              outputSigmaFile = os.path.join(dirPath, "SIGMA_BSE" + self.bsetype + '_SCR' + self.screentype + "_OC" + self.tensorComp[i] + ".OUT")
              outputLossFile = os.path.join(dirPath, "LOSS_BSE" + self.bsetype + '_SCR' + self.screentype + "_OC" + self.tensorComp[i] + ".OUT")
              # - - - - -
              try:
                with open(outputXSFile) as g:
                    xsParser = exciting_parser_xs.XSParser()
                    xsParser.parseExciton(outputXSFile, backend, excNum, excEn, excBindEn, osclStr, transCoeff) #, dftMethodSectionGindex = self.secMethodIndex,
              except FileNotFoundError:
                logger.warning("File not found: {}" .format(outputXSFile))
              except Exception as err:
                logger.error("Exception on {}" .format(__file__), exc_info=err)
              # - - - - -
              try:
                with open(outputEpsFile) as g:
                    epsParser = exciting_parser_eps.EPSParser()
                    epsParser.parseEpsilon(outputEpsFile, backend, epsEn, epsilon) #, dftMethodSectionGindex = self.secMethodIndex,
              except FileNotFoundError:
                logger.warning("File not found: {}" .format(outputEpsFile))
              except Exception as err:
                logger.error("Exception on Exciting subparser", exc_info=err)
              # - - - - -
              try:
                with open(outputSigmaFile) as g:
                    sigmaParser = exciting_parser_eps.EPSParser()
                    sigmaParser.parseEpsilon(outputSigmaFile, backend, sigmaEn, sigma) #, dftMethodSectionGindex = self.secMethodIndex,
              except FileNotFoundError:
                logger.warning("File not found: {}" .format(outputSigmaFile))
              except Exception as err:
                logger.error("Exception on {}" .format(__file__), exc_info=err)


              # with open(outputLossFile) as g:
              #     lossParser = exciting_parser_eps.EPSParser()
              #     lossParser.parseEpsilon(outputLossFile, backend, lossEn, loss) #, dftMethodSectionGindex = self.secMethodIndex,
              #     dftSingleConfigurationGindex = self.secSingleConfIndex)
          backend.addValue("x_exciting_xs_bse_number_of_components",numberOfComponents)
          backend.addValue("x_exciting_xs_bse_number_of_excitons",len(excNum))
          backend.addValue("x_exciting_xs_bse_number_of_energy_points",len(epsEn))
          backend.addValue("x_exciting_xs_bse_exciton_energies",excEn)
          backend.addValue("x_exciting_xs_bse_exciton_binding_energies",excBindEn)
          backend.addValue("x_exciting_xs_bse_exciton_oscillator_strength",osclStr)
          backend.addValue("x_exciting_xs_bse_exciton_amplitude_re",transCoeff[0])
          backend.addValue("x_exciting_xs_bse_exciton_amplitude_im",transCoeff[1])
          backend.addValue("x_exciting_xs_bse_epsilon_energies",epsEn)
          backend.addValue("x_exciting_xs_bse_epsilon_re",epsilon[0])
          backend.addValue("x_exciting_xs_bse_epsilon_im",epsilon[1])
          backend.addValue("x_exciting_xs_bse_sigma_energies",sigmaEn)
          backend.addValue("x_exciting_xs_bse_sigma_re",sigma[0])
          backend.addValue("x_exciting_xs_bse_sigma_im",sigma[1])
          backend.addValue("x_exciting_xs_bse_loss_energies",lossEn)
          backend.addValue("x_exciting_xs_bse_loss",loss[0])

      if xstype == "TDDFT":
          dielFunctEne = []
          dielFunctLoc=[[],[]]
          dielFunctNoLoc=[[],[]]
          dielTensNoSym = []
          dielTensSym = []
          lossFunctionLoc = []
          lossFunctionNoLoc = []
          qCartesian = []
          QFile = os.path.join(dirPath, "QPOINTS.OUT")
          qLattice = []
          qPlusG = []
          qPlusGCartesian = []
          qPlusGLattice = []
          qpointNumber = 0
          sigmaLoc = [[],[]]
          sigmaNoLoc = [[],[]]
          tensorComp=[]

          for files in os.listdir(dirPath):
              if files[0:9] == "SIGMA_NLF":
                  tensorComp.append(files[-13:-11])
                  #  print("tensorComp===",tensorComp)
          if self.xsTetra and self.xsAC and not self.xsNAR:
              ext = "TET_AC_NAR"
          elif not self.xsTetra and self.xsAC and not self.xsNAR:
              ext = "AC_NAR"
          elif not self.xsTetra and not self.xsAC and not self.xsNAR:
              ext = "NAR"
          elif self.xsTetra and self.xsAC and self.xsNAR:
              ext = "TET_AC"
          elif self.xsTetra and not self.xsAC and self.xsNAR:
              ext = "TET"
          elif self.xsTetra and not self.xsAC and not self.xsNAR:
              ext = "TET_NAR"
          else:
              ext=""
          try:
            with open(QFile) as g:
                while 1:
                    s = g.readline()
                    if not s: break
                    s = s.strip()
                    s = s.split()
                    if not is_number(s[1]):
                        qpointNumber = int(s[0] )
                    else:
                        qPlusGCartesian.append([])
                        qPlusGLattice.append([])
                        qLattice.append([float(s[1]),float(s[2]),float(s[3])])
                        qCartesian.append([float(s[1]),float(s[2]),float(s[3])])
                        qPlusG.append(int(s[7]))
                if self.xsTetra and self.xsAC and not self.xsNAR:
                    ext = "TET_AC_NAR"
                elif not self.xsTetra and self.xsAC and not self.xsNAR:
                    ext = "AC_NAR"
                elif not self.xsTetra and not self.xsAC and not self.xsNAR:
                    ext = "NAR"
                elif self.xsTetra and self.xsAC and self.xsNAR:
                    ext = "TET_AC"
                elif self.xsTetra and not self.xsAC and self.xsNAR:
                    ext = "TET"
                elif self.xsTetra and not self.xsAC and not self.xsNAR:
                    ext = "TET_NAR"
                else:
                    ext=""
          except FileNotFoundError:
            logger.warning("File not found: {}" .format(QFile))
          except Exception as err:
            logger.error("Exception on {}" .format(__file__), exc_info=err)
          # ---- QFile closing
          #  xstype = "BSE"
          #  self.xstype = "BSE"
          #  dummyBse = files[11:13]
          #  self.tensorComp = files[-6:-4]
          # WARNING: the following entries don't exist in Exciting Metainfo!!
          #backend.addValue("x_exciting_xs_tddft_number_of_optical_components",len(tensorComp))
          #backend.addValue("x_exciting_xs_tddft_optical_component",tensorComp)
          #backend.addValue("s",qpointNumber)

          for i in range(qpointNumber):
              dielTensSym.append([[],[]])
              dielTensNoSym.append([[],[]])
              dielFunctLoc[0].append([])
              dielFunctLoc[1].append([])
              dielFunctNoLoc[0].append([])
              dielFunctNoLoc[1].append([])
              lossFunctionLoc.append([])
              lossFunctionNoLoc.append([])
              dielFunctEne.append([])
              sigmaLoc[0].append([])
              sigmaLoc[1].append([])
              sigmaNoLoc[0].append([])
              sigmaNoLoc[1].append([])
              if i < 10:
                  qExt00 = '_QMT00'+ str(i+1)
                  qPlusGFile = os.path.join(dirPath, 'GQPOINTS' + qExt00 + '.OUT')
                  DielNoSymFile = os.path.join(dirPath, 'DIELTENS0_NOSYM' + qExt00 + '.OUT')
                  DielSymFile = os.path.join(dirPath, 'DIELTENS0' + qExt00 + '.OUT')
                  with open(qPlusGFile) as g:
                      while 1:
                          s = g.readline()
                          if not s: break
                          s = s.strip()
                          s = s.split()
                          if is_number(s[1]):
                              qPlusGLattice[i].append([float(s[1]),float(s[2]),float(s[3])])
                              qPlusGCartesian[i].append([float(s[4]),float(s[5]),float(s[6])])
                  with open(DielNoSymFile) as g:
                      while 1:
                          s = g.readline()
                          if not s: break
                          s = s.strip()
                          s = s.split()
                          if s and is_number(s[0]):
                              dielTensNoSym[i][0].append([float(s[0]),float(s[1]),float(s[2])])
                              dielTensNoSym[i][1].append([float(s[3]),float(s[4]),float(s[5])])
                  with open(DielSymFile) as g:
                      while 1:
                          s = g.readline()
                          if not s: break
                          s = s.strip()
                          s = s.split()
                          if s and is_number(s[0]):
                              dielTensSym[i][0].append([float(s[0]),float(s[1]),float(s[2])])
                              dielTensSym[i][1].append([float(s[3]),float(s[4]),float(s[5])])
                  for j in range(len(tensorComp)):
                      dielFunctLoc[0][-1].append([])
                      dielFunctLoc[1][-1].append([])
                      dielFunctNoLoc[0][-1].append([])
                      dielFunctNoLoc[1][-1].append([])
                      dielFunctEne[-1].append([])
                      lossFunctionLoc[-1].append([])
                      lossFunctionNoLoc[-1].append([])
                      sigmaLoc[0][-1].append([])
                      sigmaLoc[1][-1].append([])
                      sigmaNoLoc[0][-1].append([])
                      sigmaNoLoc[1][-1].append([])
                      #  print("j===",j)
                      #  print("ext===",ext)
                      #  print("self.tddftKernel===",self.tddftKernel)
                      #  print("qExt00===",qExt00)
                      #  print("tutto===",'EPSILON_' + ext + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT')

                      if self.tddftKernel is not None:
                        ########
                        try:
                          epsilonLocalField = os.path.join(dirPath, 'EPSILON_' + ext + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT')
                          with open(epsilonLocalField) as g:
                            while 1:
                              s = g.readline()
                              if not s: break
                              s = s.strip()
                              s = s.split()
                              ene, epsRe, epsIm = fromeV(float(s[0])),float(s[1]),float(s[2])
                              dielFunctLoc[0][-1][-1].append(epsRe)
                              dielFunctLoc[1][-1][-1].append(epsIm)
                              dielFunctEne[-1][-1].append(ene)
                          backend.addValue("x_exciting_xs_tddft_number_of_epsilon_values",len(dielFunctEne[0][0]))
                          backend.addValue("x_exciting_xs_tddft_epsilon_energies",dielFunctEne[0][0])
                          backend.addValue("x_exciting_xs_tddft_dielectric_function_local_field",dielFunctLoc)
                        except IOError:
                          logger.error("File not processable: %s" % (epsilonLocalField))

                        #########
                        try:
                          epsilonNoLocalField = os.path.join(dirPath, 'EPSILON_' + ext + 'NLF_' + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT')
                          with open(epsilonNoLocalField) as g:
                            while 1:
                              s = g.readline()
                              if not s: break
                              s = s.strip()
                              s = s.split()
                              epsRe, epsIm = float(s[1]),float(s[2])
                              dielFunctNoLoc[0][-1][-1].append(epsRe)
                              dielFunctNoLoc[1][-1][-1].append(epsIm)
                          backend.addValue("x_exciting_xs_tddft_dielectric_function_no_local_field",dielFunctNoLoc)
                        except IOError:
                          logging.error("File not processable: %s" % (epsilonNoLocalField))

                        #########
                        try:
                          lossFunctionLocalFieldFile = os.path.join(dirPath, 'LOSS_' + ext + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT')
                          with open(lossFunctionLocalFieldFile) as g:
                            while 1:
                              s = g.readline()
                              if not s: break
                              s = s.strip()
                              s = s.split()
                              #print("s===",s)
                              if s and is_number(s[0]):
                                loss = float(s[1])
                                lossFunctionLoc[-1][-1].append(loss)
                          backend.addValue("x_exciting_xs_tddft_loss_function_local_field",lossFunctionLoc)
                        except IOError:
                          logger.error("File not processable: %s" % (lossFunctionLocalFieldFile))

                        #########
                        try:
                          lossFunctionNoLocalFieldFile = os.path.join(dirPath, 'LOSS_' + ext + 'NLF_' + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT')
                          with open(lossFunctionNoLocalFieldFile) as g:
                            while 1:
                              s = g.readline()
                              if not s: break
                              s = s.strip()
                              s = s.split()
                              # print("s===",s)
                              if s and is_number(s[0]):
                                loss = float(s[1])
                                lossFunctionNoLoc[-1][-1].append(loss)
                          backend.addValue("x_exciting_xs_tddft_loss_function_no_local_field",lossFunctionNoLoc)
                        except IOError:
                          logger.error("File not processable: %s" % (lossFunctionNoLocalFieldFile))

                        #########
                        sigmaLocalFieldFile = os.path.join(dirPath, 'SIGMA_' + ext + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT')
                        try:
                          with open(sigmaLocalFieldFile) as g:
                            while 1:
                              s = g.readline()
                              if not s: break
                              s = s.strip()
                              s = s.split()
                              if s and is_number(s[0]):
                                sigmaRe, sigmaIm = float(s[1]),float(s[2])
                                sigmaLoc[0][-1][-1].append(sigmaRe)
                                sigmaLoc[1][-1][-1].append(sigmaIm)
                          backend.addValue("x_exciting_xs_tddft_sigma_local_field",sigmaLoc)
                        except IOError:
                          logger.error("File not processable: %s" % (sigmaLocalFieldFile))
                        #
                        #########
                        sigmaNoLocalFieldFile = os.path.join(dirPath, 'SIGMA_' + ext + 'NLF_' + 'FXC' + self.tddftKernel[0] + '_OC' + tensorComp[j] + qExt00 + '.OUT')
                        try:
                          with open(sigmaNoLocalFieldFile) as g:
                            while 1:
                              s = g.readline()
                              if not s: break
                              s = s.strip()
                              s = s.split()
                              if s and is_number(s[0]):
                                sigmaRe, sigmaIm = float(s[1]),float(s[2])
                                sigmaNoLoc[0][-1][-1].append(sigmaRe)
                                sigmaNoLoc[1][-1][-1].append(sigmaIm)
                          backend.addValue("x_exciting_xs_tddft_sigma_no_local_field",sigmaNoLoc)
                        except IOError:
                          logger.error("File not processable: %s" % (sigmaLocalFieldFile))

          #  dielFunctEne[-1][-1].append(ene)
          #  if s and is_number(s[0]):
          #    dielTensSym[i][0].append([float(s[0]),float(s[1]),float(s[2])])
          #    dielTensSym[i][1].append([float(s[3]),float(s[4]),float(s[5])])
          # print("loss===",lossFunctionLoc)
          # print("ext===",ext)
          # print("dielTensSym===",dielTensSym)
          # print("dielTensNoSym===",dielTensNoSym)
          #    qPlusGLattice[i].append([float(s[1]),float(s[2]),float(s[3])])
          #    qPlusGCartesian[i].append([float(s[4]),float(s[5]),float(s[6])])
          #WARNING: the following entries don't exist in the Exciting Metainfo!!
          #backend.addValue("x_exciting_xs_tddft_dielectric_tensor_sym",dielTensSym)
          #backend.addValue("x_exciting_xs_tddft_dielectric_tensor_no_sym",dielTensNoSym)

      backend.closeNonOverlappingSection("section_single_configuration_calculation")


  def onClose_x_exciting_section_geometry_optimization(self, backend, gIndex, section):
    """Trigger called when x_abinit_section_dataset is closed.
    """
    if len(self.frameSequence) > 1:
      frameGIndex = backend.openSection("section_frame_sequence")
      backend.addValue("geometry_optimization_converged", True)
      backend.closeSection("section_frame_sequence", frameGIndex)
    backend.closeSection("section_sampling_method", self.samplingGIndex)

  def onClose_section_frame_sequence(self, backend, gIndex, section):
    """Trigger called when section_framce_sequence is closed.
    """
    backend.addValue("number_of_frames_in_sequence", len(self.frameSequence))
    backend.addArrayValues("frame_sequence_local_frames_ref", np.array(self.frameSequence))
    backend.addValue("frame_sequence_to_sampling_ref", self.samplingGIndex)

    #  print("self.samplingMethod=",self.samplingMethod)
    if self.samplingMethod == "geometry_optimization":
    #  gix = backend.openSection("section_sampling_method")
    #  backend.addValue("XC_functional_name", xcName)
    #  backend.closeSection("section_sampling_method", gix)
    #  geometryForceThreshold = section["x_exciting_geometry_optimization_threshold_force"]
    #  print("geometryForceThreshold=",self.geometryForceThreshold)
        gi = backend.openSection("section_sampling_method")
        backend.addValue("geometry_optimization_threshold_force", self.geometryForceThreshold)
        backend.closeSection("section_sampling_method", gi)
    else:
        pass

  def onClose_x_exciting_section_lattice_vectors(self, backend, gIndex, section):
    latticeX = section["x_exciting_geometry_lattice_vector_x"]
    latticeY = section["x_exciting_geometry_lattice_vector_y"]
    latticeZ = section["x_exciting_geometry_lattice_vector_z"]
    cell = [[latticeX[0],latticeY[0],latticeZ[0]],
            [latticeX[1],latticeY[1],latticeZ[1]],
            [latticeX[2],latticeY[2],latticeZ[2]]]
    self.sim_cell = cell
    backend.addValue("simulation_cell", cell)

  def onClose_x_exciting_section_reciprocal_lattice_vectors(self, backend, gIndex, section):
    recLatticeX = section["x_exciting_geometry_reciprocal_lattice_vector_x"]
    recLatticeY = section["x_exciting_geometry_reciprocal_lattice_vector_y"]
    recLatticeZ = section["x_exciting_geometry_reciprocal_lattice_vector_z"]
    recCell = [[recLatticeX[0],recLatticeY[0],recLatticeZ[0]],
            [recLatticeX[1],recLatticeY[1],recLatticeZ[1]],
            [recLatticeX[2],recLatticeY[2],recLatticeZ[2]]]
    backend.addValue("x_exciting_simulation_reciprocal_cell", recCell)

  def onClose_x_exciting_section_xc(self, backend, gIndex, section):
    xcNr = section["x_exciting_xc_functional"][0]
    xc_internal_map = {
        2: ['LDA_C_PZ', 'LDA_X_PZ'],
        3: ['LDA_C_PW', 'LDA_X_PZ'],
        4: ['LDA_C_XALPHA'],
        5: ['LDA_C_VBH'],
        20: ['GGA_C_PBE', 'GGA_X_PBE'],
        21: ['GGA_C_PBE', 'GGA_X_PBE_R'],
        22: ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
        26: ['GGA_C_PBE', 'GGA_X_WC'],
        30: ['GGA_C_AM05', 'GGA_C_AM05'],
        300: ['GGA_C_BGCP', 'GGA_X_PBE'],
        406: ['HYB_GGA_XC_PBEH']
        }
    if xcNr == 100:
        dirPath = os.path.dirname(self.parser.fIn.name)
        inputGSFile = os.path.join(dirPath, "input.xml")
        try:
          with open(inputGSFile) as f:
            exciting_parser_GS_input.parseInput(f, backend)
        except FileNotFoundError:
          logger.warning("File not found: {}" .format(inputGSFile))
        except Exception as err:
          logger.error("Exception while processing file {}" .format(inputGSFile), exc_info=err)
    else:
        for xcName in xc_internal_map[xcNr]:
          self.xcName = xcName
          gi = backend.openSection("section_XC_functionals")
          backend.addValue("XC_functional_name", xcName)
          backend.closeSection("section_XC_functionals", gi)

  def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
    # Determine the Fermi energy for this calculation.
    energy_reference_fermi = None
    try:
        energy_reference_fermi = backend.superBackend.get_value('x_exciting_fermi_energy', g_index=gIndex).m
    except KeyError:
        try:
            energy_reference_fermi = backend.superBackend.get_value('gw_fermi_energy', g_index=gIndex).m
        except KeyError:
            pass
    if energy_reference_fermi is not None:
        backend.addArrayValues('energy_reference_fermi', [energy_reference_fermi, energy_reference_fermi])  # always two spin channels

    # logger.error("BASE onClose_section_single_configuration_calculation")
    backend.addValue('single_configuration_to_calculation_method_ref', self.secMethodIndex)
    backend.addValue('single_configuration_calculation_to_system_ref', self.secSystemIndex)
    # print("self.samplingMethod=",self.samplingMethod)

    ####################VOLUME TEST BEGIN################################
    ext_uri = []
    # backend.addValue("x_exciting_dummy=",self.volumeOptIndex)
    # backend.addValue("x_exciting_dummy", self.volumeOptIndex)
    if self.volumeCubOptIndex > 1:
      for j in range(1, self.volumeCubOptIndex):
        if (j<10):
          ext_uri.append(self.mainFilePath[0:-9] + 'rundir-0' + str(j) + '/INFO.OUT')
          # backend.addValue("x_exciting_dummy2", ext_uri[-1])
          # print("ext_uri===",ext_uri)
        else:
          ext_uri.append(self.mainFilePath[0:-9] + 'rundir-' + str(j) + '/INFO.OUT')
          # backend.addValue("x_exciting_dummy2", ext_uri[-1])
          # print("ext_uri===",ext_uri)
      # backend.addArrayValues("x_exciting_dummy2", np.asarray(ext_uri))
    for ref in ext_uri:
      refGindex = backend.openSection("section_calculation_to_calculation_refs")
      backend.addValue("calculation_to_calculation_external_url", ref)
      backend.addValue("calculation_to_calculation_kind", "source_calculation")
      backend.closeSection("section_calculation_to_calculation_refs", refGindex)

    ####################VOLUME TEST END############################

    if self.samplingMethod == "geometry_optimization":
      ivalue = section["x_exciting_geometry_optimization_threshold_force"]
      if ivalue is None:
        # then use default value
        self.geometryForceThreshold = self.geometryForceThreshold
        logger.warning("Found suspicious value for geometry optimization "
          "threshold force. Hint: inspect INFO.OUT, is it complete?")
      else:
        try:
          self.geometryForceThreshold = ivalue[0]
        except:
          raise
    forceX = section["x_exciting_geometry_atom_forces_x"]
    if forceX:
      forceY = section["x_exciting_geometry_atom_forces_y"]
      forceZ = section["x_exciting_geometry_atom_forces_z"]
      atoms = len(forceX)
      atom_geometry_forces = []
      for i in range(0,atoms):
          atom_geometry_forces.append([forceX[i],forceY[i],forceZ[i]])
      backend.addValue("atom_forces",atom_geometry_forces)
      # print("geometryForceThreshold=",geometryForceThreshold)
      # backend.addValue("geometry_optimization_threshold_force", geometryForceThreshold)
    # else:
    #   pass
    #
    ##############TO DO. FIX FORCES#####################
    #    forceX = section["x_exciting_atom_forces_x"]
    #    if forceX:
    #      forceY = section["x_exciting_atom_forces_y"]
    #      forceZ = section["x_exciting_atom_forces_z"]
    #      print("forceX===",forceX)
    #      print("forceY===",forceY)
    #      print("forceZ===",forceZ)
    #      forceCoreX = section["x_exciting_atom_core_forces_x"]
    #      forceCoreY = section["x_exciting_atom_core_forces_y"]
    #      forceCoreZ = section["x_exciting_atom_core_forces_z"]
    #      print("forceCoreX===",forceCoreX)
    ##      print("forceCoreY===",forceCoreY)
    #      print("forceCoreZ===",forceCoreZ)
    #      forceIBSX = section["x_exciting_atom_IBS_forces_x"]
    #      forceIBSY = section["x_exciting_atom_IBS_forces_y"]
    #      forceIBSZ = section["x_exciting_atom_IBS_forces_z"]
    #      print("forceIBSX===",forceIBSX)
    #      print("forceIBSY===",forceIBSY)
    #      print("forceIBSZ===",forceIBSZ)
    #      forceHFX = section["x_exciting_atom_HF_forces_x"]
    #      forceHFY = section["x_exciting_atom_HF_forces_y"]
    #      forceHFZ = section["x_exciting_atom_HF_forces_z"]
    #      fConv = convert_unit_function("hartree/bohr", "N")
    #      atoms = len(forceX)
    #      atom_forces = []
    #      atom_core_forces = []
    #      atom_IBS_forces = []
    #      atom_HF_forces = []
    #      for i in range(0,atoms):
    #        print("atoms===",atoms)
    #        print("i===",i)
    #        print("atom_forces===",atom_forces)
    #        print("forceX[i]===",forceX[i])
    #        print("forceY[i]===",forceY[i])
    #        print("forceZ[i]===",forceZ[i])
    #        print("forceCoreX[i]===",forceCoreX[i])
    #        print("forceCoreY[i]===",forceCoreY[i])
    #        print("forceCoreZ[i]===",forceCoreZ[i])
    #        print("forceIBSX[i]===",forceIBSX[i])
    #        print("forceIBSY[i]===",forceIBSY[i])
    #        print("forceIBSZ[i]===",forceIBSZ[i])
    #        atom_forces.append([fConv(forceX[i]),fConv(forceY[i]),fConv(forceZ[i])])
    #        atom_core_forces.append([fConv(forceCoreX[i]),fConv(forceCoreY[i]),fConv(forceCoreZ[i])])
    #        atom_IBS_forces.append([fConv(forceIBSX[i]),fConv(forceIBSY[i]),fConv(forceIBSZ[i])])
    #        atom_HF_forces.append([fConv(forceHFX[i]),fConv(forceHFY[i]),fConv(forceHFZ[i])])
    #      backend.addValue("atom_forces",atom_forces)
    #      backend.addValue("x_exciting_atom_core_forces",atom_core_forces)
    #      backend.addValue("x_exciting_atom_IBS_forces",atom_IBS_forces)
    #      backend.addValue("x_exciting_atom_HF_forces",atom_HF_forces)
    #    print("atom_forces=",atom_forces)
    #
    dirPath = os.path.dirname(self.parser.fIn.name)
    dosFile = os.path.join(dirPath, "dos.xml")
    bandFile = os.path.join(dirPath, "bandstructure.xml")
    fermiSurfFile = os.path.join(dirPath, "FERMISURF.bxsf")
    eigvalFile = os.path.join(dirPath, "EIGVAL.OUT")
    #    logger.error("done BASE onClose_section_single_configuration_calculation")

    try:
      with open(dosFile) as f:
        exciting_parser_dos.parseDos(f, backend, self.spinTreat, self.unit_cell_vol, energy_reference_fermi)
    except FileNotFoundError:
      logger.warning("File not found: {}" .format(dosFile))
    except Exception as err:
      logger.error("Exception while processing file {}" .format(dosFile), exc_info=err)

    if os.path.exists(bandFile):
      with open(bandFile) as g:
        exciting_parser_bandstructure.parseBand(g, backend, self.spinTreat)
    if os.path.exists(eigvalFile):
      eigvalGIndex = backend.openSection("section_eigenvalues")
      with open(eigvalFile) as g:
          eigvalKpoint=[]
          eigvalVal=[]
          eigvalOcc=[]
          eigvalValSpin = [[],[]]
          eigvalOccSpin = [[],[]]
          fromH = unit_conversion.convert_unit_function("hartree", "J")
          while 1:
            s = g.readline()
            if not s: break
            s = s.strip()
            if len(s) < 20:
              if "nstsv" in s.split():
                 nstsv = int(s.split()[0])
                 nstsv2=int(nstsv/2)
              elif "nkpt" in s.split():
                 nkpt = int(s.split()[0])
              continue
            elif len(s) > 50 and ('k-point,' in s.split()):
              try:
                int(s[0])  # assert this line is not a header string
              except ValueError:
                continue
              else:
                eigvalVal.append([])
                eigvalOcc.append([])
                eigvalKpoint.append([float(x) for x in s.split()[1:4]])

            else:
              try:
                int(s[0])
              except ValueError:
                # ignore headers by neglecting lines that
                # don't start with an integer
                continue
              else:
                # FIXME: EIG files with 'partial chg density' have 7 columns
                # FIXME: and we need to process them
                n, e, occ = s.split()[0:3]
                #---------
                try:  # eigenvalues 'e' could be wrongly formatted
                  enew = float(e)
                except Exception as ex:
                  if 'E' not in e.upper():
                    pieces = e.split('-')
                    if (len(pieces) != 2):
                      raise ex
                    try:
                        pieces = [ float(ii) for ii in pieces ]
                    except:
                        raise ex
                    mantissa, exponent = pieces
                    enew = mantissa * 10**(-1 * exponent)
                    logger.warning("In-house conversion '{}' -> {}" .format(e, enew))
                #---------
                eigvalVal[-1].append(fromH(enew))
                eigvalOcc[-1].append(float(occ))
          if not self.spinTreat:
            for i in range(0,nkpt):
              eigvalValSpin[0].append(eigvalVal[i][0:nstsv])
              eigvalOccSpin[0].append(eigvalOcc[i][0:nstsv])
              eigvalValSpin[1].append(eigvalVal[i][0:nstsv])
              eigvalOccSpin[1].append(eigvalOcc[i][0:nstsv])
            backend.addValue("eigenvalues_values", eigvalValSpin)
            backend.addValue("eigenvalues_occupation", eigvalOccSpin)
          else:
            for i in range(0,nkpt):
              eigvalValSpin[0].append(eigvalVal[i][0:nstsv2])
              eigvalOccSpin[0].append(eigvalOcc[i][0:nstsv2])
              eigvalValSpin[1].append(eigvalVal[i][nstsv2:nstsv])
              eigvalOccSpin[1].append(eigvalOcc[i][nstsv2:nstsv])
            backend.addValue("eigenvalues_values", eigvalValSpin)
            backend.addValue("eigenvalues_occupation", eigvalOccSpin)
          backend.addValue("eigenvalues_kpoints", eigvalKpoint)
          backend.closeSection("section_eigenvalues",eigvalGIndex)

    ########################## Parsing Fermi surface ##################
    if os.path.exists(fermiSurfFile):
      fermiGIndex = backend.openSection("x_exciting_section_fermi_surface")
      with open(fermiSurfFile) as g:
        grid = []
        all_vectors = []
        values = []
        origin = []
        vectors = []
        fermi = 0
        number_of_bands = 0
        mesh_size = 0
        fromH = unit_conversion.convert_unit_function("hartree", "J")
        while 1:
          s = g.readline()
          if not s: break
          s = s.strip()
          st = s.split()
          if len(st) == 3:
            if len(s) >= 40:
              all_vectors.append([])
              i = 0
              while i < 3:
                all_vectors[-1].append(float(st[i]))
                i += 1
            elif st[0] == "Fermi":
              fermi = fromH(float(st[2]))
            else:
              j = 0
              while j < 3:
                grid.append(int(st[j]))
                j += 1
          elif len(st) == 2:
            values.append([])
          elif len(s) >= 12 and len(st) == 1:
            try: float(st[0])
            except ValueError:
              continue
            else:
              values[-1].append(float(st[0]))
          elif len(s) < 5 and len(st) == 1:
            number_of_bands = st[0]
        mesh_size = grid[0]*grid[1]*grid[2]
        origin = all_vectors[0]
        vectors = all_vectors[1:]
        backend.addValue("x_exciting_number_of_bands_fermi_surface", int(number_of_bands))
        backend.addValue("x_exciting_number_of_mesh_points_fermi_surface", int(mesh_size))
        backend.addValue("x_exciting_fermi_energy_fermi_surface", float(fermi))
        backend.addArrayValues("x_exciting_grid_fermi_surface", np.asarray(grid))
        backend.addArrayValues("x_exciting_origin_fermi_surface", np.asarray(origin))
        backend.addArrayValues("x_exciting_vectors_fermi_surface", np.asarray(vectors))
        backend.addArrayValues("x_exciting_values_fermi_surface", np.asarray(values))
        backend.closeSection("x_exciting_section_fermi_surface",fermiGIndex)

  def onClose_x_exciting_section_spin(self, backend, gIndex, section):

    spin = section["x_exciting_spin_treatment"][0]
    spin = spin.strip()
    if spin == "spin-polarised":
      self.spinTreat = True
    else:
      self.spinTreat = False

  def onClose_section_system(self, backend, gIndex, section):

    self.unit_cell_vol = section["x_exciting_unit_cell_volume"]
    self.gmaxvr = section["x_exciting_gmaxvr"]
    self.rgkmax = section["x_exciting_rgkmax"]
    backend.addArrayValues('configuration_periodic_dimensions', np.asarray([True, True, True]))

    # next line fixes normalizer error "no lattice vectors but pbc"
    backend.addArrayValues("lattice_vectors", self.sim_cell)

    self.secSystemDescriptionIndex = gIndex

    if self.atom_pos and self.cell_format[0] == 'cartesian':
      #  print("self.atom_pos=",self.atom_pos)
      backend.addArrayValues('atom_positions', np.asarray(self.atom_pos))
    elif self.atom_pos and self.cell_format[0] == 'lattice':
      # print("aaaself.atom_pos=",self.atom_pos)
      # print("aaaself.atom_labels=",self.atom_labels)
       atoms = Atoms(self.atom_labels, self.atom_pos, cell=[(1, 0, 0),(0, 1, 0),(0, 0, 1)])
       atoms.set_cell(self.sim_cell, scale_atoms=True)
       self.atom_pos = atoms.get_positions()
       backend.addArrayValues('atom_positions', np.asarray(self.atom_pos))
    if self.atom_labels is not None:
      # print("aaaself.atom_labels=",self.atom_labels)
      backend.addArrayValues('atom_labels', np.asarray(self.atom_labels))
    self.atom_labels = []

    excSmearingKind = section["x_exciting_smearing_type"]

    smearing_internal_map = {
        "Gaussian": ['gaussian'],
        "Methfessel-Paxton": ['methfessel-paxton'],
        "Fermi-Dirac": ['fermi'],
        "Extended": ['tetrahedra']
        }

    if self.samplingMethod is not "geometry_optimization":
        for smName in smearing_internal_map[excSmearingKind[0]]:
          backend.addValue("smearing_kind", smName)
    else:
        pass

  def onClose_x_exciting_section_atoms_group(self, backend, gIndex, section):
    # print("start.self.atom_labels=",self.atom_labels)
    fromB = unit_conversion.convert_unit_function("bohr", "m")
    formt = section['x_exciting_atom_position_format']
    if self.samplingMethod is not "geometry_optimization":
      pass
    else:
      if self.atom_pos is not None: self.atom_pos = []
      if self.atom_labels is not None: self.atom_labels = []
    self.cell_format = formt
    pos = [section['x_exciting_geometry_atom_positions_' + i] for i in ['x', 'y', 'z']]
    # print("pos=",pos)
    pl = [len(comp) for comp in pos]
    natom = pl[0]
    if pl[1] != natom or pl[2] != natom:
      raise Exception("invalid number of atoms in various components %s" % pl)
    for i in range(natom):
      # print("nattom=",natom)
      # print("i=",i)
      # print("natom=",natom)
      # print("[pos[0][i]=",pos[0][i])
      # print("[pos[1][i]=",pos[1][i])
      # print("[pos[2][i]=",pos[2][i])
      # print("self.atom_pos=",self.atom_pos)
      if formt[0] == 'cartesian':
        self.atom_pos.append([fromB(pos[0][i]), fromB(pos[1][i]), fromB(pos[2][i])])
      else:
        # print("self.atom_labels_prima=",self.atom_labels)
        # print("self.atom_pos_prima=",self.atom_pos)
        self.atom_pos.append([pos[0][i], pos[1][i], pos[2][i]])
        # print("self.atom_pos_dopo=",self.atom_pos)
        # print("self.atom_labels_dopo=",self.atom_labels)
        # print("natom=",natom)
        # print("section['x_exciting_geometry_atom_labels']=",section['x_exciting_geometry_atom_labels'])
        # print("self.samplingMethod[0]=",self.samplingMethod)
    if self.samplingMethod is not "geometry_optimization":
        # print("prima=",self.atom_labels)
        self.atom_labels = self.atom_labels + (section['x_exciting_geometry_atom_labels'] * natom)
        # print("self.atom_labels=",self.atom_labels)
        # print("section['x_exciting_geometry_atom_labels']=",section['x_exciting_geometry_atom_labels'])
    else:
        self.atom_labels = self.atom_labels + section['x_exciting_geometry_atom_labels']
        # print("self.self.XSSetGIndexatom_labels_dopodopo=",self.atom_labels)

  def onClose_section_method(self, backend, gIndex, section):
    if gIndex == self.XSSetGIndex and self.xstype == "TDDFT":
        # tmk: if 'input.xml' is missing, then the next are 'NoneType'
        try:
          self.tddftKernel = section["x_exciting_xs_tddft_xc_kernel"]
        except:
          pass
        try:
          self.xsTetra = section["x_exciting_xs_tetra"][0]
        except:
          pass
        try:
          self.xsAC = section["x_exciting_xs_tddft_analytic_continuation"][0]
        except:
          pass
        try:
          self.xsNAR = section["x_exciting_xs_tddft_anti_resonant_dielectric"][0]
        except:
          pass
    if gIndex == self.secMethodIndex:
      backend.addValue('electronic_structure_method', "DFT")
      try:
        energy_thresh = section["x_exciting_scf_threshold_energy_change"][0]
      except:
        pass
      try:
        potential_thresh = section["x_exciting_scf_threshold_potential_change_list"][0]
      except:
        pass
      try:
        charge_thresh = section["x_exciting_scf_threshold_charge_change_list"][0]
      except:
        pass
      if section["x_exciting_scf_threshold_force_change_list"]:
        force_thresh = section["x_exciting_scf_threshold_force_change_list"][0]
        try:
          backend.addValue('x_exciting_scf_threshold_force_change', force_thresh)
        except:
          pass
      try:
        backend.addValue('scf_threshold_energy_change', energy_thresh)
      except:
        pass
      try:
        backend.addValue('x_exciting_scf_threshold_potential_change', potential_thresh)
      except:
        pass
      try:
        backend.addValue('x_exciting_scf_threshold_charge_change', charge_thresh)
      except:
        pass
      ##########BELOW VOLUME OPTIMIZATION######################
      #    if self.volumeOptIndex>1:
      #      self.volumeOpt = True
      #      optGindex = backend.openSection("section_method")
      #      backend.addValue("x_exciting_volume_optimization", self.volumeOpt)
      #      backend.closeSection("section_method", optGindex)
      #
      #    ext_uri = []
      #    backend.addValue("x_exciting_dummy", self.volumeOptIndex)
      #    if self.volumeOptIndex > 1:
      #      for j in range(1, self.volumeOptIndex):
      #        if (j<10):
      #          ext_uri.append(self.mainFilePath[0:-9] + 'rundir-0' + str(j) + '/INFO.OUT')
      #        else:
      #          ext_uri.append(self.mainFilePath[0:-9] + 'rundir-' + str(j) + '/INFO.OUT')
      #    for ref in ext_uri:
      #      refGindex = backend.openSection("section_calculation_to_calculation_refs")
      #      backend.addValue("calculation_to_calculation_external_url", ref)
      #      backend.addValue("calculation_to_calculation_kind", "source_calculation")
      #      backend.closeSection("section_calculation_to_calculation_refs", refGindex)

mainFileDescription = \
    SM(name = "root matcher",
       startReStr = "",
       weak = True,
       subMatchers = [
         SM(name = "header",
         startReStr = r"\s*(\||\+|\*)\s*EXCITING\s*(?P<program_version>[-a-zA-Z0-9]+)\s*started\s*=",
         fixedStartValues={'program_name': 'exciting', 'program_basis_set_type': '(L)APW+lo' },
            sections = ["section_run", "section_method"],
         subMatchers = [
	   SM(name = 'input',
              startReStr = r"(\||\+|\*)\sStarting initialization",
              endReStr = r"(\||\+|\*)\sEnding initialization",
              sections = ['section_system'],
              subMatchers = [
                SM(startReStr = r"\sLattice vectors \(cartesian\) :",
                sections = ["x_exciting_section_lattice_vectors"],
                subMatchers = [

    SM(startReStr = r"\s*(?P<x_exciting_geometry_lattice_vector_x__bohr>[-+0-9.]+)\s+(?P<x_exciting_geometry_lattice_vector_y__bohr>[-+0-9.]+)\s+(?P<x_exciting_geometry_lattice_vector_z__bohr>[-+0-9.]+)", repeats = True)
                ]),
                SM(startReStr = r"\sReciprocal lattice vectors \(cartesian\) :",
                sections = ["x_exciting_section_reciprocal_lattice_vectors"],
                subMatchers = [

    SM(startReStr = r"\s*(?P<x_exciting_geometry_reciprocal_lattice_vector_x__bohr_1>[-+0-9.]+)\s+(?P<x_exciting_geometry_reciprocal_lattice_vector_y__bohr_1>[-+0-9.]+)\s+(?P<x_exciting_geometry_reciprocal_lattice_vector_z__bohr_1>[-+0-9.]+)", repeats = True)
                ]),
    SM(r"\s*Unit cell volume\s*:\s*(?P<x_exciting_unit_cell_volume__bohr3>[-0-9.]+)"),
    SM(r"\s*Brillouin zone volume\s*:\s*(?P<x_exciting_brillouin_zone_volume__bohr_3>[-0-9.]+)"),
    SM(r"\s*Species\s*:\s*[0-9]\s*\((?P<x_exciting_geometry_atom_labels>[-a-zA-Z0-9]+)\)", repeats = True,
      sections = ["x_exciting_section_atoms_group"],
       subMatchers = [
        SM(r"\s*muffin-tin radius\s*:\s*(?P<x_exciting_muffin_tin_radius__bohr>[-0-9.]+)"),
        SM(r"\s*# of radial points in muffin-tin\s*:\s*(?P<x_exciting_muffin_tin_points>[-0-9.]+)"),
        SM(startReStr = r"\s*atomic positions\s*\((?P<x_exciting_atom_position_format>[-a-zA-Z]+)\)\s*:\s*",
           endReStr = r"\s*magnetic fields\s*",
           subMatchers = [
                    SM(r"\s*(?P<x_exciting_geometry_atom_number>[+0-9]+)\s*:\s*(?P<x_exciting_geometry_atom_positions_x>[-+]?[0-9.]+)\s*(?P<x_exciting_geometry_atom_positions_y>[-+]?[0-9.]+)\s*(?P<x_exciting_geometry_atom_positions_z>[-+]?[0-9.]+)", repeats = True)
         ]) #,
        # SM(startReStr = r"\s*magnetic fields\s*\((?P<x_exciting_magnetic_field_format>[-a-zA-Z]+)\)\s*:\s*",
        #   subMatchers = [
        #     SM(r"\s*(?P<x_exciting_MT_external_magnetic_field_atom_number>[+0-9]+)\s*:\s*(?P<x_exciting_MT_external_magnetic_field_x>[-+0-9.]+)\s*(?P<x_exciting_MT_external_magnetic_field_y>[-+0-9.]+)\s*(?P<x_exciting_MT_external_magnetic_field_z>[-+0-9.]+)", repeats = True)
        #  ])
    ]),
    SM(r"\s*Total number of atoms per unit cell\s*:\s*(?P<x_exciting_number_of_atoms>[-0-9.]+)"),
    SM(r"\s*Spin treatment\s*:\s*(?P<x_exciting_spin_treatment>[-a-zA-Z\s*]+)",
       sections = ["x_exciting_section_spin"]),
    SM(r"\s*k-point grid\s*:\s*(?P<x_exciting_number_kpoint_x>[-0-9.]+)\s+(?P<x_exciting_number_kpoint_y>[-0-9.]+)\s+(?P<x_exciting_number_kpoint_z>[-0-9.]+)"),
    SM(r"\s*k-point offset\s*:\s*(?P<x_exciting_kpoint_offset_x>[-0-9.]+)\s+(?P<x_exciting_kpoint_offset_y>[-0-9.]+)\s+(?P<x_exciting_kpoint_offset_z>[-0-9.]+)"),
    SM(r"\s*Total number of k-points\s*:\s*(?P<x_exciting_number_kpoints>[-0-9.]+)"),
    SM(r"\s*R\^MT_min \* \|G\+k\|_max \(rgkmax\)\s*:\s*(?P<x_exciting_rgkmax>[-0-9.]+)"),
    SM(r"\s*Maximum \|G\+k\| for APW functions\s*:\s*(?P<x_exciting_gkmax__bohr_1>[-0-9.]+)"),
    SM(r"\s*Maximum \|G\| for potential and density\s*:\s*(?P<x_exciting_gmaxvr__bohr_1>[-0-9.]+)"),
    SM(r"\s*G-vector grid sizes\s*:\s*(?P<x_exciting_gvector_size_x>[-0-9.]+)\s+(?P<x_exciting_gvector_size_y>[-0-9.]+)\s+(?P<x_exciting_gvector_size_z>[-0-9.]+)"),
    SM(r"\s*Total number of G-vectors\s*:\s*(?P<x_exciting_gvector_total>[-0-9.]+)"),
    SM(startReStr = r"\s*Maximum angular momentum used for\s*",
        subMatchers = [
          SM(r"\s*APW functions\s*:\s*(?P<x_exciting_lmaxapw>[-0-9.]+)")
        ]),
    SM(r"\s*Total nuclear charge\s*:\s*(?P<x_exciting_nuclear_charge>[-0-9.]+)"),
    SM(r"\s*Total electronic charge\s*:\s*(?P<x_exciting_electronic_charge>[-0-9.]+)"),
    SM(r"\s*Total core charge\s*:\s*(?P<x_exciting_core_charge_initial>[-0-9.]+)"),
    SM(r"\s*Total valence charge\s*:\s*(?P<x_exciting_valence_charge_initial>[-0-9.]+)"),
    SM(r"\s*Effective Wigner radius, r_s\s*:\s*(?P<x_exciting_wigner_radius>[-0-9.]+)"),
    SM(r"\s*Number of empty states\s*:\s*(?P<x_exciting_empty_states>[-0-9.]+)"),
    SM(r"\s*Total number of valence states\s*:\s*(?P<x_exciting_valence_states>[-0-9.]+)"),
    SM(r"\s*Maximum Hamiltonian size\s*:\s*(?P<x_exciting_hamiltonian_size>[-0-9.]+)"),
    SM(r"\s*Maximum number of plane-waves\s*:\s*(?P<x_exciting_pw>[-0-9.]+)"),
    SM(r"\s*Total number of local-orbitals\s*:\s*(?P<x_exciting_lo>[-0-9.]+)"),
    SM(startReStr = r"\s*Exchange-correlation type\s*:\s*(?P<x_exciting_xc_functional>[-0-9.]+)",
       sections = ['x_exciting_section_xc']),
    SM(r"\s*Smearing scheme\s*:\s*(?P<x_exciting_smearing_type>[-a-zA-Z0-9]+)"),
    SM(r"\s*Smearing width\s*:\s*(?P<smearing_width__hartree>[-0-9.]+)"),
    SM(r"\s*Using\s*(?P<x_exciting_potential_mixing>[-a-zA-Z\s*]+)\s*potential mixing")
    ]),
            SM(name = "single configuration iteration",
              startReStr = r"(\||\+|\*)\s*Self-consistent loop started\s*\+",
              sections = ["section_single_configuration_calculation"],
              repeats = True,
              subMatchers = [
                SM(name = "scfi totE",
                 startReStr =r"(\||\+|\*)\s*SCF iteration number\s*:",
                  sections = ["section_scf_iteration"],
                  repeats = True,
                  subMatchers = [
                   SM(r"\s*Total energy\s*:\s*(?P<energy_total_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Fermi energy\s*:\s*(?P<x_exciting_fermi_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Kinetic energy\s*:\s*(?P<electronic_kinetic_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Coulomb energy\s*:\s*(?P<x_exciting_coulomb_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Exchange energy\s*:\s*(?P<x_exciting_exchange_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Correlation energy\s*:\s*(?P<x_exciting_correlation_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Sum of eigenvalues\s*:\s*(?P<energy_sum_eigenvalues_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Effective potential energy\s*:\s*(?P<x_exciting_effective_potential_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Coulomb potential energy\s*:\s*(?P<x_exciting_coulomb_potential_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*xc potential energy\s*:\s*(?P<energy_XC_potential_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Hartree energy\s*:\s*(?P<x_exciting_hartree_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Electron-nuclear energy\s*:\s*(?P<x_exciting_electron_nuclear_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Nuclear-nuclear energy\s*:\s*(?P<x_exciting_nuclear_nuclear_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Madelung energy\s*:\s*(?P<x_exciting_madelung_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Core-electron kinetic energy\s*:\s*(?P<x_exciting_core_electron_kinetic_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Absolute change in total energy   (target)\s*:\s*(?P<energy_change_scf_iteration__hartree>[-0-9.]+)\s*(\s*(?P<scf_threshold_energy_change__hartree>[-0-9.]+))"),
                   SM(r"\s*DOS at Fermi energy \(states\/Ha\/cell\)\s*:\s*(?P<x_exciting_dos_fermi_scf_iteration__hartree_1>[-0-9.]+)"),
                   SM(r"\s*core\s*:\s*(?P<x_exciting_core_charge_scf_iteration>[-0-9.]+)"),
                   SM(r"\s*core leakage\s*:\s*(?P<x_exciting_core_leakage_scf_iteration>[-0-9.]+)"),
                   SM(r"\s*valence\s*:\s*(?P<x_exciting_valence_charge_scf_iteration>[-0-9.]+)"),
                   SM(r"\s*interstitial\s*:\s*(?P<x_exciting_interstitial_charge_scf_iteration>[-0-9.]+)"),
                   SM(r"\s*total charge in muffin-tins\s*:\s*(?P<x_exciting_total_MT_charge_scf_iteration>[-0-9.]+)"),
                   SM(r"\s*Estimated fundamental gap\s*:\s*(?P<x_exciting_gap_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Wall time \(seconds\)\s*:\s*(?P<x_exciting_time_scf_iteration>[-0-9.]+)"),
                   SM(r"\s*RMS change in effective potential \(target\)\s*:\s*(?P<x_exciting_effective_potential_convergence_scf_iteration__hartree>[0-9]+\.[0-9]*([E]?[-]?[0-9]+))\s*\(\s*(?P<x_exciting_scf_threshold_potential_change_list__hartree>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\)"),
                   SM(r"\s*Absolute change in total energy\s*\(target\)\s*:\s*(?P<x_exciting_energy_convergence_scf_iteration>[0-9]+\.[0-9]*([E]?[-]?[0-9]+))\s*\(\s*(?P<x_exciting_scf_threshold_energy_change__hartree>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\)"),
                   SM(r"\s*Absolute change in total energy\s*\(target\)\s*\:\s*\s*\(\s*[0-9]+\.[0-9]+\s*[0-9]+\.[0-9]+\s*[0-9]+\.[0-9]+\s*\)\s*avg =\s*[0-9]+\.[0-9]+\s*\(\s*(?P<x_exciting_scf_threshold_energy_change__hartree>[0-9.]+E?[-][0-9]+)\s*\)"),
                   SM(r"\s*Charge distance\s*\(target\)\s*:\s*(?P<x_exciting_charge_convergence_scf_iteration>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\s*\(\s*(?P<x_exciting_scf_threshold_charge_change_list>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\)"),
                   SM(r"\s*Abs. change in max-nonIBS-force\s*\(target\)\s*:\s*(?P<x_exciting_force_convergence_scf_iteration>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\s*\(\s*(?P<x_exciting_scf_threshold_force_change_list>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\)")
                  ]),
                SM(name="final_quantities",
                  startReStr = r"(\||\+|\*) Convergence targets achieved. Performing final SCF iteration\s*\+",
                  endReStr = r"(\||\+|\*) Self-consistent loop stopped\s*\+",
                   subMatchers = [
                     SM(r"\s*Total energy\s*:\s*(?P<energy_total__hartree>[-0-9.]+)"),
                     SM(r"\s*Fermi energy\s*:\s*(?P<x_exciting_fermi_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Kinetic energy\s*:\s*(?P<electronic_kinetic_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Coulomb energy\s*:\s*(?P<x_exciting_coulomb_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Exchange energy\s*:\s*(?P<x_exciting_exchange_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Correlation energy\s*:\s*(?P<x_exciting_correlation_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Sum of eigenvalues\s*:\s*(?P<energy_sum_eigenvalues__hartree>[-0-9.]+)"),
                     SM(r"\s*Effective potential energy\s*:\s*(?P<x_exciting_effective_potential_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Coulomb potential energy\s*:\s*(?P<x_exciting_coulomb_potential_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*xc potential energy\s*:\s*(?P<energy_XC_potential__hartree>[-0-9.]+)"),
                     SM(r"\s*Hartree energy\s*:\s*(?P<x_exciting_hartree_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Electron-nuclear energy\s*:\s*(?P<x_exciting_electron_nuclear_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Nuclear-nuclear energy\s*:\s*(?P<x_exciting_nuclear_nuclear_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Madelung energy\s*:\s*(?P<x_exciting_madelung_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Core-electron kinetic energy\s*:\s*(?P<x_exciting_core_electron_kinetic_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*DOS at Fermi energy \(states\/Ha\/cell\)\s*:\s*(?P<x_exciting_dos_fermi__hartree_1>[-0-9.]+)"),
                     SM(r"\s*core leakage\s*:\s*(?P<x_exciting_core_leakage>[-0-9.]+)"),
                     SM(r"\s*interstitial\s*:\s*(?P<x_exciting_interstitial_charge>[-0-9.]+)"),
                     SM(r"\s*total charge in muffin-tins\s*:\s*(?P<x_exciting_total_MT_charge>[-0-9.]+)"),
                     SM(r"\s*Estimated fundamental gap\s*:\s*(?P<x_exciting_gap__hartree>[-0-9.]+)")
                   ]) #,
                #                SM(name="final_forces",
                ##                  startReStr = r"(\||\+|\*) Writing atomic positions and forces\s*\-",
                #                  startReStr = r"\s*Total atomic forces including IBS \(cartesian\) \s*:",
                #                  endReStr = r"(\||\+|\*)\s*Groundstate module stopped\s*\*",
                ##                  endReStr = r"\s* Atomic force components including IBS \(cartesian\)\s*:",
                #                  floating = True,
                #                   subMatchers = [
                ##                     SM(name="total_forces",
                ##                     startReStr = r"\s*Total atomic forces including IBS \(cartesian\)\s*:",
                #                       SM(r"\s*atom\s*[0-9]+\s*[A-Za-z]+\s*\:\s*(?P<x_exciting_atom_forces_x>[-0-9.]+)\s*(?P<x_exciting_atom_forces_y>[-0-9.]+)\s*(?P<x_exciting_atom_forces_z>[-0-9.]+)",
                #                          repeats = True )
                ######                     subMatchers = [
                ######                     SM(r"\s*atom\s*(?P<x_exciting_store_total_forces>[0-9]+\s*[A-Za-z]+\s*\:+\s*[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+)",
                ######                          repeats = True)
                ######                   ] )
                ##)
                ##                     print ("number atoms=", x_exciting_number_of_atoms)
                ##                     SM(name="force_components",
                ##                     startReStr = r"\s*Atomic force components including IBS \(cartesian\)\s*:",
                ##                     forwardMatch = True,
                ##                     subMatchers = [
                ##                     SM(r"\s*atom\s*(?P<x_exciting_store_total_forces>[0-9]+\s*[A-Za-z]+\s*\:+\s*[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+\s*[A-Za-z]+\s*[A-Za-z]+)", weak = True),
                ##                     SM(r"\s*(?P<x_exciting_store_total_forces>\s*\:+\s*[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+\s*[A-Za-z]+\s*[A-Za-z]+)"),
                ##                     SM(r"\s*(?P<x_exciting_store_total_forces>\s*\:+\s*[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+\s*[A-Za-z]+\s*[A-Za-z]+)")
                ##                     SM(r"\s*(?P<x_exciting_store_total_forces>\s*\:+\s*[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+\s*[A-Za-z]+\s*[A-Za-z]+)"),
                ##                   ]
                ##                    )
                #                   ]),
                #                 SM(name="force_components",
                #                  startReStr = r"\s* Atomic force components including IBS \(cartesian\)\s*:",
                #                  endReStr = r"(\||\+|\*)\s* Groundstate module stopped\s* \*",
                #                  subMatchers = [
                ##                  startReStr = r"\s* Atomic force components including IBS \(cartesian\)\s*:",
                #                   SM(r"\s*atom\s*[0-9]+\s*[A-Za-z]+\s*\:\s*(?P<x_exciting_atom_HF_forces_x>[-0-9.]+)\s*(?P<x_exciting_atom_HF_forces_y>[-0-9.]+)\s*(?P<x_exciting_atom_HF_forces_z>[-0-9.]+)\s*HF force",
                #                     repeats = True,
                #                     floating = True),
                #                   SM(r"\s*\:\s*(?P<x_exciting_atom_core_forces_x>[-0-9.]+)\s*(?P<x_exciting_atom_core_forces_y>[-0-9.]+)\s*(?P<x_exciting_atom_core_forces_z>[-0-9.]+)\s*core correction",
                #                     repeats = True,
                #                     floating = True),
                #                   SM(r"\s*\:\s*(?P<x_exciting_atom_IBS_forces_x>[-0-9.]+)\s*(?P<x_exciting_atom_IBS_forces_y>[-0-9.]+)\s*(?P<x_exciting_atom_IBS_forces_z>[-0-9.]+)\s*IBS correction",
                #                     repeats = True,
                #                     floating = True),
                ##                   SM(r"(?P<x_exciting_store_total_forces>.*)",
                ##                          repeats = True,
                #                ] )
               ]
            ),
            SM(name = "geometry optimization",
              startReStr = r"(\||\+|\*)\s*Structure-optimization module started*\s*\*",
              sections = ["section_sampling_method","x_exciting_section_geometry_optimization"],
              # fixedStartValues={'sampling_method': 'geometry_optimization'},
              # repeats = True,
              subMatchers = [
                   SM(name = "optimization steps",
                   startReStr = r"(\||\+|\*)\s*Optimization step\s*(?P<x_exciting_geometry_optimization_step>[-0-9]+)\s*\(method = (?P<x_exciting_geometry_optimization_method>[A-Za-z]+)\)\s*\-",
                   sections = ["section_single_configuration_calculation"],
                  # SM(r"\s*Output level for this task is set to normal\s*"),
                  # SM(r"(\||\+|\*)\s*Optimization step (?P<x_exciting_geometry_optimization_step>[-0-9]+)\: Initialize optimization\s*\-"),
                   repeats = True,
                   subMatchers = [
                   SM(r"\s*Maximum force magnitude\s*\(target\)\s*:\s*(?P<x_exciting_maximum_force_magnitude__hartree_bohr_1>[0-9]+\.[0-9]*([E]?[-]?[0-9]+))\s*\(\s*(?P<x_exciting_geometry_optimization_threshold_force__hartree_bohr_1>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\)"),
                   SM(r"\s*Total energy at this optimization step\s*:\s*(?P<energy_total__hartree>[-0-9.]+)"),
                   SM(startReStr = r"\s*Atomic positions at this step \s*\((?P<x_exciting_atom_position_format>[-a-zA-Z]+)\)\s*:\s*",
                  # endReStr = r"\s*Total atomic forces including IBS \(cartesian\) \:",
                  # weak = True,
                   sections = ["section_system","x_exciting_section_atoms_group"],
                  # endReStr = r"\s*magnetic fields\s*",
           subMatchers = [
                    SM(r"\s*atom\s*(?P<x_exciting_atom_number>[+0-9]+)\s*(?P<x_exciting_geometry_atom_labels>[A-Za-z]+)\s*\:\s*(?P<x_exciting_geometry_atom_positions_x>[-+0-9.]+)\s*(?P<x_exciting_geometry_atom_positions_y>[-+0-9.]+)\s*(?P<x_exciting_geometry_atom_positions_z>[-+0-9.]+)", repeats = True)
         ]),
                   SM(startReStr = r"\s*Total atomic forces including IBS \(cartesian\) \:",
                    endReStr = r"\s*Time spent in this optimization step\s*\:\s*(?P<x_exciting_geometry_dummy>[+0-9.]+)\s*seconds",
                    subMatchers = [
                    SM(r"\s*atom\s*[+0-9]+\s*[A-Za-z]+\s*\:\s*(?P<x_exciting_geometry_atom_forces_x__hartree_bohr_1>[-+0-9.]+)\s*(?P<x_exciting_geometry_atom_forces_y__hartree_bohr_1>[-+0-9.]+)\s*(?P<x_exciting_geometry_atom_forces_z__hartree_bohr_1>[-+0-9.]+)", repeats = True)
         ])
         ])
           ]
           )
          ])
    ])




parserInfo = {
  "name": "exciting_parser",
  "version": "1.0"
}

cachingLevelForMetaName = {
    "x_exciting_geometry_lattice_vector_x":CachingLevel.Cache,
    "x_exciting_geometry_lattice_vector_y":CachingLevel.Cache,
    "x_exciting_geometry_lattice_vector_z":CachingLevel.Cache,
    "x_exciting_section_lattice_vectors": CachingLevel.Ignore,
    "x_exciting_geometry_reciprocal_lattice_vector_x":CachingLevel.Cache,
    "x_exciting_geometry_reciprocal_lattice_vector_y":CachingLevel.Cache,
    "x_exciting_geometry_reciprocal_lattice_vector_z":CachingLevel.Cache,
    "x_exciting_section_reciprocal_lattice_vectors": CachingLevel.Ignore,
    "x_exciting_atom_forces_x":CachingLevel.Cache,
    "x_exciting_atom_forces_y":CachingLevel.Cache,
    "x_exciting_atom_forces_z":CachingLevel.Cache,
    "x_exciting_atom_HF_forces_x":CachingLevel.Cache,
    "x_exciting_atom_HF_forces_y":CachingLevel.Cache,
    "x_exciting_atom_HF_forces_z":CachingLevel.Cache,
    "x_exciting_atom_core_forces_x":CachingLevel.Cache,
    "x_exciting_atom_core_forces_y":CachingLevel.Cache,
    "x_exciting_atom_core_forces_z":CachingLevel.Cache,
    "x_exciting_atom_IBS_forces_x":CachingLevel.Cache,
    "x_exciting_atom_IBS_forces_y":CachingLevel.Cache,
    "x_exciting_atom_IBS_forces_z":CachingLevel.Cache,
    "x_exciting_geometry_atom_forces_x":CachingLevel.Cache,
    "x_exciting_geometry_atom_forces_y":CachingLevel.Cache,
    "x_exciting_geometry_atom_forces_z":CachingLevel.Cache
 }


class ExcitingParser():
    """ A proper class envolop for running this parser from within python. """
    def __init__(self, backend, **kwargs):
        self.backend_factory = backend

    def parse(self, mainfile):
        from unittest.mock import patch
        logging.info('exciting parser started')
        logging.getLogger('nomadcore').setLevel(logging.WARNING)
        backend = self.backend_factory('exciting.nomadmetainfo.json')
        with patch.object(sys, 'argv', ['<exe>', '--uri', 'nmd://uri', mainfile]):
            mainFunction(
                mainFileDescription,
                None,
                parserInfo,
                cachingLevelForMetaName = cachingLevelForMetaName,
                superContext=ExcitingParserContext(),
                superBackend=backend)
            # print("#"*50 + "\nPARSING ENDED. NORMALIZER FOLLOWS..."+"\n"*5)# tmk:
        return backend

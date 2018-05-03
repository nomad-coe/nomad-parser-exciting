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
import setup_paths
import numpy as np
from nomadcore.simple_parser import AncillaryParser, CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM, mainFunction
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.caching_backend import CachingLevel
from nomadcore.unit_conversion import unit_conversion
from nomadcore.unit_conversion.unit_conversion import convert_unit_function
import os, sys, json, exciting_parser_dos,exciting_parser_bandstructure, exciting_parser_gw, exciting_parser_GS_input
from ase import Atoms
import logging

class ExcitingParserContext(object):

  def __init__(self):
    self.parser = None

  def initialize_values(self):
    self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv
    self.atom_pos = []
    self.atom_labels = []

  def startedParsing(self, path, parser):
    self.parser=parser
    self.initialize_values()
#    self.atom_pos = []
#    self.atom_labels = []
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
    self.energy_thresh = []
    self.samplingMethod = None
    self.secSamplingMethodIndex = None
    self.geometryForceThreshold = 0
    self.frameSequence = []
    self.samplingGIndex = 0
    self.dummy = 0

  def onOpen_section_sampling_method(self, backend, gIndex, section):
    self.secSamplingMethodIndex = gIndex
    backend.addValue("sampling_method", "geometry_optimization")
#    print("self.secSamplingMethodIndex=",self.secSamplingMethodIndex)
    self.samplingMethod = "geometry_optimization"
#    print("self.samplingMethod=",self.samplingMethod)
    
  def onOpen_section_system(self, backend, gIndex, section):
    self.secSystemIndex = gIndex

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

  def onClose_x_exciting_section_geometry_optimization(self, backend, gIndex, section):
    """Trigger called when x_abinit_section_dataset is closed.
    """
#    print("len(self.frameSequence)=",len(self.frameSequence))
    if len(self.frameSequence) > 1:
      frameGIndex = backend.openSection("section_frame_sequence")
#        if section["x_abinit_geometry_optimization_converged"] is not None:
#          if section["x_abinit_geometry_optimization_converged"][-1] == "are converged":
      backend.addValue("geometry_optimization_converged", True)
#          else:
#            backend.addValue("geometry_optimization_converged", False)
      backend.closeSection("section_frame_sequence", frameGIndex)
#    self.frameSequence = []
    backend.closeSection("section_sampling_method", self.samplingGIndex)

  def onClose_section_frame_sequence(self, backend, gIndex, section):
    """Trigger called when section_framce_sequence is closed.
    """
    backend.addValue("number_of_frames_in_sequence", len(self.frameSequence))
    backend.addArrayValues("frame_sequence_local_frames_ref", np.array(self.frameSequence))
    backend.addValue("frame_sequence_to_sampling_ref", self.samplingGIndex)

#    print("self.samplingMethod=",self.samplingMethod)
    if self.samplingMethod == "geometry_optimization":
#        gix = backend.openSection("section_sampling_method")
#        backend.addValue("XC_functional_name", xcName)
#        backend.closeSection("section_sampling_method", gix)
#        geometryForceThreshold = section["x_exciting_geometry_optimization_threshold_force"]
#        print("geometryForceThreshold=",self.geometryForceThreshold)
        gi = backend.openSection("section_sampling_method")
        backend.addValue("geometry_optimization_threshold_force", self.geometryForceThreshold)
        backend.closeSection("section_sampling_method", gi)
    else:
        pass

    mainFile = self.parser.fIn.fIn.name
    dirPath = os.path.dirname(self.parser.fIn.name)
    gw_File = os.path.join(dirPath, "GW_INFO.OUT")
    gwFile = os.path.join(dirPath, "GWINFO.OUT")
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
        with open(inputGSFile) as f:
            exciting_parser_GS_input.parseInput(f, backend)
    else:
        for xcName in xc_internal_map[xcNr]:
          self.xcName = xcName
          gi = backend.openSection("section_XC_functionals")
          backend.addValue("XC_functional_name", xcName)
          backend.closeSection("section_XC_functionals", gi)

  def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
#    logging.error("BASE onClose_section_single_configuration_calculation")
    backend.addValue('single_configuration_to_calculation_method_ref', self.secMethodIndex)
    backend.addValue('single_configuration_calculation_to_system_ref', self.secSystemIndex)
#    print("self.samplingMethod=",self.samplingMethod)
    if self.samplingMethod == "geometry_optimization":
        self.geometryForceThreshold = section["x_exciting_geometry_optimization_threshold_force"][0]
    forceX = section["x_exciting_geometry_atom_forces_x"]
    if forceX:
      forceY = section["x_exciting_geometry_atom_forces_y"]
      forceZ = section["x_exciting_geometry_atom_forces_z"]
      atoms = len(forceX)
      atom_geometry_forces = []
      for i in range(0,atoms):
          atom_geometry_forces.append([forceX[i],forceY[i],forceZ[i]])
      backend.addValue("atom_forces",atom_geometry_forces)
#        print("geometryForceThreshold=",geometryForceThreshold)
#        backend.addValue("geometry_optimization_threshold_force", geometryForceThreshold)
#    else:
#        pass
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
#    logging.error("done BASE onClose_section_single_configuration_calculation")

    if os.path.exists(dosFile):
      with open(dosFile) as f:
        exciting_parser_dos.parseDos(f, backend, self.spinTreat, self.unit_cell_vol)
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
            elif len(s) > 50:
              eigvalVal.append([])
              eigvalOcc.append([])
              eigvalKpoint.append([float(x) for x in s.split()[1:4]])
            else:
              try: int(s[0])
              except ValueError:
                continue
              else:
                n, e, occ = s.split()
                eigvalVal[-1].append(fromH(float(e)))
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

##########################Parsing Fermi surface##################

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
    backend.addArrayValues('configuration_periodic_dimensions', np.asarray([True, True, True]))

    self.secSystemDescriptionIndex = gIndex

    if self.atom_pos and self.cell_format[0] == 'cartesian':
       backend.addArrayValues('atom_positions', np.asarray(self.atom_pos))
    elif self.atom_pos and self.cell_format[0] == 'lattice':
#       print("aaaself.atom_pos=",self.atom_pos)
#       print("aaaself.atom_labels=",self.atom_labels)
       atoms = Atoms(self.atom_labels, self.atom_pos, cell=[(1, 0, 0),(0, 1, 0),(0, 0, 1)])
       atoms.set_cell(self.sim_cell, scale_atoms=True)
       self.atom_pos = atoms.get_positions()
       backend.addArrayValues('atom_positions', np.asarray(self.atom_pos))
    if self.atom_labels is not None:
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
    fromB = unit_conversion.convert_unit_function("bohr", "m")
    formt = section['x_exciting_atom_position_format']
    if self.atom_pos is not None: self.atom_pos = []
    if self.atom_labels is not None: self.atom_labels = []
    self.cell_format = formt
    pos = [section['x_exciting_geometry_atom_positions_' + i] for i in ['x', 'y', 'z']]
#    print("pos=",pos)
    pl = [len(comp) for comp in pos]
    natom = pl[0]
    if pl[1] != natom or pl[2] != natom:
      raise Exception("invalid number of atoms in various components %s" % pl)
    for i in range(natom):
#      print("i=",i)
#      print("natom=",natom)
#      print("[pos[0][i]=",pos[0][i])
#      print("[pos[1][i]=",pos[1][i])
#      print("[pos[2][i]=",pos[2][i])
#      print("self.atom_pos=",self.atom_pos)
      if formt[0] == 'cartesian':
        self.atom_pos.append([fromB(pos[0][i]), fromB(pos[1][i]), fromB(pos[2][i])])
      else:
#        print("self.atom_labels_prima=",self.atom_labels)
#        print("self.atom_pos_prima=",self.atom_pos)
        self.atom_pos.append([pos[0][i], pos[1][i], pos[2][i]])
#        print("self.atom_pos_dopo=",self.atom_pos)
#        print("self.atom_labels_dopo=",self.atom_labels)
#    print("natom=",natom)
#    print("section['x_exciting_geometry_atom_labels']=",section['x_exciting_geometry_atom_labels'])
#    print("self.samplingMethod[0]=",self.samplingMethod)
    if self.samplingMethod is not "geometry_optimization":
        self.atom_labels = self.atom_labels + (section['x_exciting_geometry_atom_labels'] * natom)
    else:
        self.atom_labels = self.atom_labels + section['x_exciting_geometry_atom_labels']
#    print("self.atom_labels_dopodopo=",self.atom_labels)

  def onClose_section_method(self, backend, gIndex, section):
    if gIndex == self.secMethodIndex:
      backend.addValue('electronic_structure_method', "DFT")
      energy_thresh = section["x_exciting_scf_threshold_energy_change"][0]
      potential_thresh = section["x_exciting_scf_threshold_potential_change_list"][0]
      charge_thresh = section["x_exciting_scf_threshold_charge_change_list"][0]
      if section["x_exciting_scf_threshold_force_change_list"]:
        force_thresh = section["x_exciting_scf_threshold_force_change_list"][0]
        backend.addValue('x_exciting_scf_threshold_force_change', force_thresh)
      backend.addValue('scf_threshold_energy_change', energy_thresh)
      backend.addValue('x_exciting_scf_threshold_potential_change', potential_thresh)
      backend.addValue('x_exciting_scf_threshold_charge_change', charge_thresh)
#      backend.addValue('x_exciting_scf_threshold_force_change', force_thresh)

mainFileDescription = \
    SM(name = "root matcher",
       startReStr = "",
       weak = True,
       subMatchers = [
         SM(name = "header",
         startReStr = r"\s*\|\s*EXCITING\s*(?P<program_version>[-a-zA-Z0-9]+)\s*started\s*=",
         fixedStartValues={'program_name': 'exciting', 'program_basis_set_type': '(L)APW+lo' },
            sections = ["section_run", "section_method"],
         subMatchers = [
	   SM(name = 'input',
              startReStr = r"\|\sStarting initialization",
              endReStr = r"\|\sEnding initialization",
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
        SM(r"\s*muffin-tin radius\s*:\s*(?P<x_exciting_muffin_tin_radius__bohr>[-0-9.]+)", repeats = True),
        SM(r"\s*# of radial points in muffin-tin\s*:\s*(?P<x_exciting_muffin_tin_points>[-0-9.]+)", repeats = True),
        SM(startReStr = r"\s*atomic positions\s*\((?P<x_exciting_atom_position_format>[-a-zA-Z]+)\)\s*:\s*",
           endReStr = r"\s*magnetic fields\s*",
           subMatchers = [
                    SM(r"\s*(?P<x_exciting_geometry_atom_number>[+0-9]+)\s*:\s*(?P<x_exciting_geometry_atom_positions_x>[-+0-9.]+)\s*(?P<x_exciting_geometry_atom_positions_y>[-+0-9.]+)\s*(?P<x_exciting_geometry_atom_positions_z>[-+0-9.]+)", repeats = True)
         ]) #,
#        SM(startReStr = r"\s*magnetic fields\s*\((?P<x_exciting_magnetic_field_format>[-a-zA-Z]+)\)\s*:\s*",
#           subMatchers = [
#                    SM(r"\s*(?P<x_exciting_MT_external_magnetic_field_atom_number>[+0-9]+)\s*:\s*(?P<x_exciting_MT_external_magnetic_field_x>[-+0-9.]+)\s*(?P<x_exciting_MT_external_magnetic_field_y>[-+0-9.]+)\s*(?P<x_exciting_MT_external_magnetic_field_z>[-+0-9.]+)", repeats = True)
#         ])
    ]),
    SM(r"\s*Total number of atoms per unit cell\s*:\s*(?P<x_exciting_number_of_atoms>[-0-9.]+)"),
    SM(r"\s*Spin treatment\s*:\s*(?P<x_exciting_spin_treatment>[-a-zA-Z\s*]+)",
       sections = ["x_exciting_section_spin"]),
    SM(r"\s*k-point grid\s*:\s*(?P<x_exciting_number_kpoint_x>[-0-9.]+)\s+(?P<x_exciting_number_kpoint_y>[-0-9.]+)\s+(?P<x_exciting_number_kpoint_z>[-0-9.]+)"),
    SM(r"\s*k-point offset\s*:\s*(?P<x_exciting_kpoint_offset_x>[-0-9.]+)\s+(?P<x_exciting_kpoint_offset_y>[-0-9.]+)\s+(?P<x_exciting_kpoint_offset_z>[-0-9.]+)"),
    SM(r"\s*Total number of k-points\s*:\s*(?P<x_exciting_number_kpoints>[-0-9.]+)"),
    SM(r"\s*R\^MT_min \* \|G\+k\|_max \(rgkmax\)\s*:\s*(?P<x_exciting_rgkmax__bohr>[-0-9.]+)"),
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
              startReStr = r"\|\s*Self-consistent loop started\s*\+",
              sections = ["section_single_configuration_calculation"],
              repeats = True,
              subMatchers = [
                SM(name = "scfi totE",
                 startReStr =r"\|\s*SCF iteration number\s*:",
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
                   SM(r"\s*Charge distance\s*\(target\)\s*:\s*(?P<x_exciting_charge_convergence_scf_iteration>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\s*\(\s*(?P<x_exciting_scf_threshold_charge_change_list>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\)"),
                   SM(r"\s*Abs. change in max-nonIBS-force\s*\(target\)\s*:\s*(?P<x_exciting_force_convergence_scf_iteration>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\s*\(\s*(?P<x_exciting_scf_threshold_force_change_list>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\)")
                  ]),
                SM(name="final_quantities",
                  startReStr = r"\| Convergence targets achieved. Performing final SCF iteration\s*\+",
                  endReStr = r"\| Self-consistent loop stopped\s*\+",
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
##                  startReStr = r"\| Writing atomic positions and forces\s*\-",
#                  startReStr = r"\s*Total atomic forces including IBS \(cartesian\) \s*:",
#                  endReStr = r"\|\s*Groundstate module stopped\s*\*",
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
#                  endReStr = r"\|\s* Groundstate module stopped\s* \*",
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
              startReStr = r"\|\s*Structure-optimization module started*\s*\*",
              sections = ["section_sampling_method","x_exciting_section_geometry_optimization"],
#              fixedStartValues={'sampling_method': 'geometry_optimization'},
#              repeats = True,
              subMatchers = [
                   SM(name = "optimization steps",
                   startReStr = r"\|\s*Optimization step\s*(?P<x_exciting_geometry_optimization_step>[-0-9]+)\s*\(method = (?P<x_exciting_geometry_optimization_method>[A-Za-z]+)\)\s*\-",
                   sections = ["section_single_configuration_calculation"],
#                   SM(r"\s*Output level for this task is set to normal\s*"),
#                   SM(r"\|\s*Optimization step (?P<x_exciting_geometry_optimization_step>[-0-9]+)\: Initialize optimization\s*\-"),
                   repeats = True,
                   subMatchers = [
                   SM(r"\s*Maximum force magnitude\s*\(target\)\s*:\s*(?P<x_exciting_maximum_force_magnitude__hartree_bohr_1>[0-9]+\.[0-9]*([E]?[-]?[0-9]+))\s*\(\s*(?P<x_exciting_geometry_optimization_threshold_force__hartree_bohr_1>[0-9]\.[0-9]*([E]?[-]?[0-9]+))\)"),
                   SM(r"\s*Total energy at this optimization step\s*:\s*(?P<energy_total__hartree>[-0-9.]+)"),
                   SM(startReStr = r"\s*Atomic positions at this step \s*\((?P<x_exciting_atom_position_format>[-a-zA-Z]+)\)\s*:\s*",
#                   endReStr = r"\s*Total atomic forces including IBS \(cartesian\) \:",
#                   weak = True,
                   sections = ["section_system","x_exciting_section_atoms_group"],
#           endReStr = r"\s*magnetic fields\s*",
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

metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/exciting.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath = metaInfoPath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)

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
if __name__ == "__main__":
    mainFunction(mainFileDescription, metaInfoEnv, parserInfo, cachingLevelForMetaName = cachingLevelForMetaName, superContext=ExcitingParserContext())

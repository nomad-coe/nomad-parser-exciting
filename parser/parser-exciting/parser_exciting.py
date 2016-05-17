import setup_paths
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM, mainFunction
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.caching_backend import CachingLevel
import os, sys, json, exciting_parser_dos

class ExcitingParserContext(object):

  def startedParsing(self, path, parser):
    self.parser=parser

  def onClose_exciting_section_lattice_vectors(self, backend, gIndex, section):
    latticeX = section["exciting_geometry_lattice_vector_x"]
    latticeY = section["exciting_geometry_lattice_vector_y"]
    latticeZ = section["exciting_geometry_lattice_vector_z"]
    cell = [[latticeX[0],latticeY[0],latticeZ[0]],
            [latticeX[1],latticeY[1],latticeZ[1]],
            [latticeX[2],latticeY[2],latticeZ[2]]]
    backend.addValue("simulation_cell", cell)

  def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
    dirPath = os.path.dirname(self.parser.fIn.name)
    dosFile = os.path.join(dirPath, "dos.xml")
    if os.path.exists(dosFile):
      with open(dosFile) as f:
        exciting_parser_dos.parseDos(f, backend)


# description of the input


mainFileDescription = \
    SM(name = "root matcher",
       startReStr = "",
       weak = True,
       subMatchers = [
         SM(name = "header",
         startReStr = r"\s*\|\s*EXCITING\s*(?P<program_version>[-a-zA-Z0-9]+)\s*started\s*=",
         fixedStartValues={'program_name': 'exciting', 'program_basis_set_type': '(L)APW+lo' },
         sections = ["section_run"],
         subMatchers = [
	   SM(name = 'input',
              startReStr = r"\|\sStarting initialization",
              endReStr = r"\|\sEnding initialization",
              sections = ['section_system'],
              subMatchers = [
                SM(startReStr = r"\sLattice vectors \(cartesian\) :",
                   endReStr = r"\sReciprocal lattice vectors \(cartesian\) :",
                sections = ["exciting_section_lattice_vectors"],
                subMatchers = [

    SM(startReStr = r"\s*(?P<exciting_geometry_lattice_vector_x__bohr>[-+0-9.]+)\s+(?P<exciting_geometry_lattice_vector_y__bohr>[-+0-9.]+)\s+(?P<exciting_geometry_lattice_vector_z__bohr>[-+0-9.]+)", repeats = True)
                ]),
    SM(r"\s*Unit cell volume\s*:\s*(?P<exciting_unit_cell_volume__bohr3>[-0-9.]+)"),
    SM(r"\s*Brillouin zone volume\s*:\s*(?P<exciting_brillouin_zone_volume__bohr_3>[-0-9.]+)")
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
                   SM(r"\s*Fermi energy\s*:\s*(?P<exciting_fermi_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Kinetic energy\s*:\s*(?P<electronic_kinetic_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Coulomb energy\s*:\s*(?P<exciting_coulomb_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Exchange energy\s*:\s*(?P<exciting_exchange_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Correlation energy\s*:\s*(?P<exciting_correlation_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Sum of eigenvalues\s*:\s*(?P<energy_sum_eigenvalues_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Effective potential energy\s*:\s*(?P<exciting_effective_potential_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Coulomb potential energy\s*:\s*(?P<exciting_coulomb_potential_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*xc potential energy\s*:\s*(?P<energy_XC_potential_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Hartree energy\s*:\s*(?P<exciting_hartree_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Electron-nuclear energy\s*:\s*(?P<exciting_electron_nuclear_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Nuclear-nuclear energy\s*:\s*(?P<exciting_nuclear_nuclear_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Madelung energy\s*:\s*(?P<exciting_madelung_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Core-electron kinetic energy\s*:\s*(?P<exciting_core_electron_kinetic_energy_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Absolute change in total energy   (target)\s*:\s*(?P<energy_change_scf_iteration__hartree>[-0-9.]+)")
                  ]),
                SM(name="final_quantities",
                  startReStr = r"\| Convergence targets achieved. Performing final SCF iteration\s*\+",
                  endReStr = r"\| Self-consistent loop stopped\s*\+",
                   subMatchers = [
                     SM(r"\s*Total energy\s*:\s*(?P<energy_total__hartree>[-0-9.]+)"),
                     SM(r"\s*Fermi energy\s*:\s*(?P<exciting_fermi_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Kinetic energy\s*:\s*(?P<electronic_kinetic_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Coulomb energy\s*:\s*(?P<exciting_coulomb_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Exchange energy\s*:\s*(?P<exciting_exchange_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Correlation energy\s*:\s*(?P<exciting_correlation_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Sum of eigenvalues\s*:\s*(?P<energy_sum_eigenvalues__hartree>[-0-9.]+)"),
                     SM(r"\s*Effective potential energy\s*:\s*(?P<exciting_effective_potential_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Coulomb potential energy\s*:\s*(?P<exciting_coulomb_potential_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*xc potential energy\s*:\s*(?P<energy_XC_potential__hartree>[-0-9.]+)"),
                     SM(r"\s*Hartree energy\s*:\s*(?P<exciting_hartree_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Electron-nuclear energy\s*:\s*(?P<exciting_electron_nuclear_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Nuclear-nuclear energy\s*:\s*(?P<exciting_nuclear_nuclear_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Madelung energy\s*:\s*(?P<exciting_madelung_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*Core-electron kinetic energy\s*:\s*(?P<exciting_core_electron_kinetic_energy__hartree>[-0-9.]+)")
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
                            "exciting_geometry_lattice_vector_x":CachingLevel.Cache,
                            "exciting_geometry_lattice_vector_y":CachingLevel.Cache,
                            "exciting_geometry_lattice_vector_z":CachingLevel.Cache,
                            "exciting_section_lattice_vectors": CachingLevel.Ignore
                          }

if __name__ == "__main__":
#    for name in metaInfoEnv.infoKinds:
#        print 'nameo', name
    mainFunction(mainFileDescription, metaInfoEnv, parserInfo, cachingLevelForMetaName = cachingLevelForMetaName, superContext=ExcitingParserContext())

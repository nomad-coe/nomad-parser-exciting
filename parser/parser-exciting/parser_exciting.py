import setup_paths
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM, mainFunction
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.caching_backend import CachingLevel
import os, sys, json, exciting_parser_dos,exciting_parser_bandstructure

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

  def onClose_exciting_section_xc(self, backend, gIndex, section):
    xcNr = section["exciting_xc_functional"][0]
    xc_internal_map = {
        2: ['LDA_C_PZ', 'LDA_X_PZ'],
        3: ['LDA_C_PW'],
        4: ['LDA_C_XALPHA'],
        5: ['LDA_C_VBH'],
        20: ['GGA_C_PBE'],
        21: ['GGA_X_PBE_R'],
        22: ['GGA_C_PBE_SOL'],
        26: ['GGA_X_WC'],
        30: ['GGA_C_AM05'],
        300: ['GGA_C_BGCP'],
        406: ['HYB_GGA_C_PBEH']
        }
    for xcName in xc_internal_map[xcNr]:
      gi = backend.openSection("section_XC_functionals")
      backend.addValue("XC_functional_name", xcName)
      backend.closeSection("section_XC_functionals", gi)

  def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
    dirPath = os.path.dirname(self.parser.fIn.name)
    dosFile = os.path.join(dirPath, "dos.xml")
    bandFile = os.path.join(dirPath, "bandstructure.xml")
    if os.path.exists(dosFile):
      with open(dosFile) as f:
        exciting_parser_dos.parseDos(f, backend)
    if os.path.exists(bandFile):
      with open(bandFile) as g:
        exciting_parser_bandstructure.parseBand(g, backend)


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
                sections = ["exciting_section_lattice_vectors"],
                subMatchers = [

    SM(startReStr = r"\s*(?P<exciting_geometry_lattice_vector_x__bohr>[-+0-9.]+)\s+(?P<exciting_geometry_lattice_vector_y__bohr>[-+0-9.]+)\s+(?P<exciting_geometry_lattice_vector_z__bohr>[-+0-9.]+)", repeats = True)
                ]),
                SM(startReStr = r"\sReciprocal lattice vectors \(cartesian\) :",
                sections = ["exciting_section_reciprocal_lattice_vectors"],
                subMatchers = [

    SM(startReStr = r"\s*(?P<exciting_geometry_reciprocal_lattice_vector_x__bohr_1>[-+0-9.]+)\s+(?P<exciting_geometry_reciprocal_lattice_vector_y__bohr_1>[-+0-9.]+)\s+(?P<exciting_geometry_reciprocal_lattice_vector_z__bohr_1>[-+0-9.]+)", repeats = True)
                ]),
    SM(r"\s*Unit cell volume\s*:\s*(?P<exciting_unit_cell_volume__bohr3>[-0-9.]+)"),
    SM(r"\s*Brillouin zone volume\s*:\s*(?P<exciting_brillouin_zone_volume__bohr_3>[-0-9.]+)"),
    SM(r"\s*Species\s*:\s*[0-9]\s*\((?P<exciting_geometry_atom_labels>[-a-zA-Z0-9]+)\)", repeats = True,
       subMatchers = [
    SM(r"\s*muffin-tin radius\s*:\s*(?P<exciting_muffin_tin_radius__bohr>[-0-9.]+)"),
    SM(r"\s*# of radial points in muffin-tin\s*:\s*(?P<exciting_muffin_tin_points>[-0-9.]+)"),
    SM(startReStr = r"\s*atomic positions\s*\(lattice\)\s*:\s*",
      subMatchers = [
        SM(r"\s*[0-9]\s*:\s*(?P<exciting_geometry_atom_positions_x__bohr>[-+0-9.]+)\s*(?P<exciting_geometry_atom_positions_y__bohr>[-+0-9.]+)\s*(?P<exciting_geometry_atom_positions_z__bohr>[-+0-9.]+)", repeats = True)
      ])
    ]),
    SM(r"\s*k-point grid\s*:\s*(?P<exciting_number_kpoint_x>[-0-9.]+)\s+(?P<exciting_number_kpoint_y>[-0-9.]+)\s+(?P<exciting_number_kpoint_z>[-0-9.]+)"),
    SM(r"\s*k-point offset\s*:\s*(?P<exciting_kpoint_offset_x>[-0-9.]+)\s+(?P<exciting_kpoint_offset_y>[-0-9.]+)\s+(?P<exciting_kpoint_offset_z>[-0-9.]+)"),
    SM(r"\s*Total number of k-points\s*:\s*(?P<exciting_number_kpoints>[-0-9.]+)"),
    SM(r"\s*R\^MT_min \* \|G\+k\|_max \(rgkmax\)\s*:\s*(?P<exciting_rgkmax__bohr>[-0-9.]+)"),
    SM(r"\s*Maximum \|G\| for potential and density\s*:\s*(?P<exciting_gmaxvr__bohr_1>[-0-9.]+)"),
    SM(r"\s*G-vector grid sizes\s*:\s*(?P<exciting_gvector_size_x>[-0-9.]+)\s+(?P<exciting_gvector_size_y>[-0-9.]+)\s+(?P<exciting_gvector_size_z>[-0-9.]+)"),
    SM(r"\s*Total number of G-vectors\s*:\s*(?P<exciting_gvector_total>[-0-9.]+)"),
    SM(startReStr = r"\s*Maximum angular momentum used for\s*",
        subMatchers = [
          SM(r"\s*APW functions\s*:\s*(?P<exciting_lmaxapw>[-0-9.]+)")
        ]),
    SM(r"\s*Total nuclear charge\s*:\s*(?P<exciting_nuclear_charge>[-0-9.]+)"),
    SM(r"\s*Total electronic charge\s*:\s*(?P<exciting_electronic_charge>[-0-9.]+)"),
    SM(r"\s*Total core charge\s*:\s*(?P<exciting_core_charge>[-0-9.]+)"),
    SM(r"\s*Total valence charge\s*:\s*(?P<exciting_valence_charge>[-0-9.]+)"),
    SM(r"\s*Effective Wigner radius, r_s\s*:\s*(?P<exciting_wigner_radius>[-0-9.]+)"),
    SM(r"\s*Number of empty states\s*:\s*(?P<exciting_empty_states>[-0-9.]+)"),
    SM(r"\s*Total number of valence states\s*:\s*(?P<exciting_valence_states>[-0-9.]+)"),
    SM(r"\s*Maximum Hamiltonian size\s*:\s*(?P<exciting_hamiltonian_size>[-0-9.]+)"),
    SM(r"\s*Maximum number of plane-waves\s*:\s*(?P<exciting_pw>[-0-9.]+)"),
    SM(r"\s*Total number of local-orbitals\s*:\s*(?P<exciting_lo>[-0-9.]+)"),
    SM(r"\s*Smearing scheme\s*:\s*(?P<exciting_smearing_type>[-a-zA-Z0-9]+)"),
    SM(r"\s*Smearing width\s*:\s*(?P<exciting_smearing_width__hartree>[-0-9.]+)"),
    SM(r"\s*Using\s*(?P<exciting_potential_mixing>[-a-zA-Z\s*]+)\s*potential mixing"),
    SM(startReStr = r"\s*Exchange-correlation type\s*:\s*(?P<exciting_xc_functional>[-0-9.]+)",
       sections = ['exciting_section_xc'])
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
                   SM(r"\s*Absolute change in total energy   (target)\s*:\s*(?P<energy_change_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*DOS at Fermi energy \(states\/Ha\/cell\)\s*:\s*(?P<exciting_dos_fermi_scf_iteration__hartree_1>[-0-9.]+)"),
                   SM(r"\s*core leakage\s*:\s*(?P<exciting_core_leakage_scf_iteration>[-0-9.]+)"),
                   SM(r"\s*interstitial\s*:\s*(?P<exciting_interstitial_charge_scf_iteration>[-0-9.]+)"),
                   SM(r"\s*total charge in muffin-tins\s*:\s*(?P<exciting_total_MT_charge_scf_iteration>[-0-9.]+)"),
                   SM(r"\s*Estimated fundamental gap\s*:\s*(?P<exciting_gap_scf_iteration__hartree>[-0-9.]+)"),
                   SM(r"\s*Wall time \(seconds\)\s*:\s*(?P<exciting_time_scf_iteration>[-0-9.]+)"),
                   SM(r"\s*RMS change in effective potential \(target\)\s*:\s*(?P<exciting_effective_potential_convergence_scf_iteration>[0-9]\.[0-9]*([E]?[-]?[0-9]+))"),
                   SM(r"\s*Absolute change in total energy\s*\(target\)\s*:\s*(?P<exciting_energy_convergence_scf_iteration>[0-9]\.[0-9]*([E]?[-]?[0-9]+))"),
                   SM(r"\s*Charge distance\s*\(target\)\s*:\s*(?P<exciting_charge_convergence_scf_iteration>[0-9]\.[0-9]*([E]?[-]?[0-9]+))")
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
                     SM(r"\s*Core-electron kinetic energy\s*:\s*(?P<exciting_core_electron_kinetic_energy__hartree>[-0-9.]+)"),
                     SM(r"\s*DOS at Fermi energy \(states\/Ha\/cell\)\s*:\s*(?P<exciting_dos_fermi__hartree_1>[-0-9.]+)"),
                     SM(r"\s*core leakage\s*:\s*(?P<exciting_core_leakage>[-0-9.]+)"),
                     SM(r"\s*interstitial\s*:\s*(?P<exciting_interstitial_charge>[-0-9.]+)"),
                     SM(r"\s*total charge in muffin-tins\s*:\s*(?P<exciting_total_MT_charge>[-0-9.]+)"),
                     SM(r"\s*Estimated fundamental gap\s*:\s*(?P<exciting_gap__hartree>[-0-9.]+)")
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

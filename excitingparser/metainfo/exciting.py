#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
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
import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import public

m_package = Package(
    name='exciting_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='exciting.nomadmetainfo.json'))


class x_exciting_section_geometry_optimization(MSection):
    '''
    section for geometry optimization
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_geometry_optimization'))


class x_exciting_section_atoms_group(MSection):
    '''
    a group of atoms of the same type
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_atoms_group'))

    x_exciting_geometry_atom_labels = Quantity(
        type=str,
        shape=[],
        description='''
        labels of atom
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_atom_labels'))

    x_exciting_geometry_atom_number = Quantity(
        type=str,
        shape=[],
        description='''
        number to identify the atoms of a species
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_atom_number'))

    x_exciting_atom_number = Quantity(
        type=str,
        shape=[],
        description='''
        number to identify the atoms of a species in the geometry optimization
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_number'))

    x_exciting_atom_label = Quantity(
        type=str,
        shape=[],
        description='''
        labels of atoms in geometry optimization
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_label'))

    x_exciting_MT_external_magnetic_field_atom_number = Quantity(
        type=str,
        shape=[],
        description='''
        number to identify the atoms of a species on which a magnetic field is applied
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_MT_external_magnetic_field_atom_number'))

    x_exciting_geometry_atom_positions = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        atomic positions
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_atom_positions'))

    x_exciting_geometry_atom_positions_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        x component of atomic position
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_atom_positions_x'))

    x_exciting_geometry_atom_positions_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        y component of atomic position
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_atom_positions_y'))

    x_exciting_geometry_atom_positions_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        z component of atomic position
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_atom_positions_z'))

    x_exciting_MT_external_magnetic_field_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        x component of the magnetic field
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_MT_external_magnetic_field_x'))

    x_exciting_MT_external_magnetic_field_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        y component of the magnetic field
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_MT_external_magnetic_field_y'))

    x_exciting_MT_external_magnetic_field_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        z component of the magnetic field
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_MT_external_magnetic_field_z'))

    x_exciting_muffin_tin_points = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        muffin-tin points
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_muffin_tin_points'))

    x_exciting_muffin_tin_radius = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        muffin-tin radius
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_muffin_tin_radius'))

    x_exciting_atom_position_format = Quantity(
        type=str,
        shape=[],
        description='''
        whether the atomic positions are given in cartesian or vector coordinates
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_position_format'))

    x_exciting_magnetic_field_format = Quantity(
        type=str,
        shape=[],
        description='''
        whether the magnetic field is given in cartesian or vector coordinates
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_magnetic_field_format'))


class x_exciting_section_bandstructure(MSection):
    '''
    bandstructure values
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_bandstructure'))

    x_exciting_band_number_of_vertices = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number of vertices along the kpoint path used for the bandstructure plot
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_band_number_of_vertices'))

    x_exciting_band_number_of_kpoints = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number of points along the kpoint path used for the bandstructure plot
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_band_number_of_kpoints'))

    x_exciting_band_vertex_labels = Quantity(
        type=str,
        shape=['x_exciting_band_number_of_vertices'],
        description='''
        labels of the vertices along the kpoint path used for the bandstructure plot
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_band_vertex_labels'))

    x_exciting_band_vertex_coordinates = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_band_number_of_vertices', 3],
        description='''
        coordinates of the vertices along the kpoint path used for the bandstructure plot
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_band_vertex_coordinates'))

    x_exciting_band_structure_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String to specify the kind of band structure (either electronic or vibrational).
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_band_structure_kind'))

    x_exciting_band_number_of_eigenvalues = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number of eigenvalues per k-point
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_band_number_of_eigenvalues'))

    x_exciting_band_k_points = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_band_number_of_kpoints'],
        description='''
        Fractional coordinates of the k points (in the basis of the reciprocal-lattice
        vectors) for which the electronic energy are given.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_band_k_points'))

    x_exciting_band_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'x_exciting_band_number_of_kpoints', 'x_exciting_band_number_of_eigenvalues'],
        unit='joule',
        description='''
        $k$-dependent energies of the electronic band structure.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_band_energies'))

    x_exciting_band_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Bandstructure energy values
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_band_value'))


class x_exciting_section_dos(MSection):
    '''
    dos values
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_dos'))

    x_exciting_dos_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        energy value for a dos point
        ''',
        categories=[public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_dos_energy'))

    x_exciting_dos_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / joule',
        description='''
        Density of states values
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_dos_value'))


class x_exciting_section_fermi_surface(MSection):
    '''
    Fermi surface values
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_fermi_surface'))

    x_exciting_fermi_energy_fermi_surface = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Fermi energy for Fermi surface
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_fermi_energy_fermi_surface'))

    x_exciting_grid_fermi_surface = Quantity(
        type=np.dtype(np.int32),
        shape=[3],
        description='''
        number of points in the mesh to calculate the Fermi surface
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_grid_fermi_surface'))

    x_exciting_number_of_bands_fermi_surface = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of bands for fermi surface
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_number_of_bands_fermi_surface'))

    x_exciting_number_of_mesh_points_fermi_surface = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of mesh points for fermi surface
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_number_of_mesh_points_fermi_surface'))

    x_exciting_origin_fermi_surface = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Origin (in lattice coordinate) of the region where the Fermi surface is calculated
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_origin_fermi_surface'))

    x_exciting_values_fermi_surface = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_number_of_bands_fermi_surface', 'x_exciting_number_of_mesh_points_fermi_surface'],
        unit='joule',
        description='''
        Fermi surface values
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_values_fermi_surface'))

    x_exciting_vectors_fermi_surface = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description='''
        Vectors (in lattice coordinate) defining the region where the Fermi surface is
        calculated
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_vectors_fermi_surface'))


class x_exciting_section_lattice_vectors(MSection):
    '''
    lattice vectors
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_lattice_vectors'))

    x_exciting_geometry_lattice_vector_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        x component of lattice vector
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_lattice_vector_x'))

    x_exciting_geometry_lattice_vector_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        y component of lattice vector
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_lattice_vector_y'))

    x_exciting_geometry_lattice_vector_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        z component of lattice vector
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_lattice_vector_z'))


class x_exciting_section_reciprocal_lattice_vectors(MSection):
    '''
    reciprocal lattice vectors
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_reciprocal_lattice_vectors'))

    x_exciting_geometry_reciprocal_lattice_vector_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / meter',
        description='''
        x component of reciprocal lattice vector
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_reciprocal_lattice_vector_x'))

    x_exciting_geometry_reciprocal_lattice_vector_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / meter',
        description='''
        y component of reciprocal lattice vector
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_reciprocal_lattice_vector_y'))

    x_exciting_geometry_reciprocal_lattice_vector_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / meter',
        description='''
        z component of reciprocal lattice vector
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_reciprocal_lattice_vector_z'))


class x_exciting_section_spin(MSection):
    '''
    section for exciting spin treatment
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_spin'))

    x_exciting_spin_treatment = Quantity(
        type=str,
        shape=[],
        description='''
        Spin treatment
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_spin_treatment'))


class x_exciting_section_xc(MSection):
    '''
    index for exciting functional
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_xc'))

    x_exciting_xc_functional = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        index for exciting functional
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xc_functional'))


class section_single_configuration_calculation(public.section_single_configuration_calculation):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_single_configuration_calculation'))

    x_exciting_atom_forces = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_number_of_atoms', 3],
        unit='newton',
        description='''
        Forces acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_forces'))

    x_exciting_atom_IBS_forces = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_number_of_atoms', 3],
        unit='newton',
        description='''
        IBS correction to the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_IBS_forces'))

    x_exciting_geometry_optimization_method = Quantity(
        type=str,
        shape=[],
        description='''
        Geometry optimization method
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_optimization_method'))

    x_exciting_geometry_optimization_step = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Geometry optimization step
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_optimization_step'))

    x_exciting_atom_core_forces = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_number_of_atoms', 3],
        unit='newton',
        description='''
        core correction to the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_core_forces'))

    x_exciting_atom_HF_forces = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_number_of_atoms', 3],
        unit='newton',
        description='''
        HF correction to the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_HF_forces'))

    x_exciting_atom_IBS_forces_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        x-component of the IBS correction to the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_IBS_forces_x'))

    x_exciting_atom_IBS_forces_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        y-component of the IBS correction to the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_IBS_forces_y'))

    x_exciting_atom_IBS_forces_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        z-component of the IBS correction to the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_IBS_forces_z'))

    x_exciting_atom_core_forces_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        x-component of the core correction to the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_core_forces_x'))

    x_exciting_atom_core_forces_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        y-component of the core correction to the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_core_forces_y'))

    x_exciting_atom_core_forces_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        z-component of the core correction to the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_core_forces_z'))

    x_exciting_atom_HF_forces_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        x-component of the HF Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_HF_forces_x'))

    x_exciting_atom_HF_forces_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        y-component of the HF Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_HF_forces_y'))

    x_exciting_atom_HF_forces_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        z-component of the HF Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_HF_forces_z'))

    x_exciting_atom_forces_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        x-component of the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_forces_x'))

    x_exciting_atom_forces_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        y-component of the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_forces_y'))

    x_exciting_atom_forces_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        z-component of the Force acting on the atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_atom_forces_z'))

    x_exciting_core_electron_kinetic_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Core-electron kinetic energy final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_core_electron_kinetic_energy'))

    x_exciting_core_leakage = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Core leakage
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_core_leakage'))

    x_exciting_correlation_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Correlation energy final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_correlation_energy'))

    x_exciting_coulomb_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Coulomb energy final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_coulomb_energy'))

    x_exciting_coulomb_potential_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Coulomb potential energy final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_coulomb_potential_energy'))

    x_exciting_dos_fermi = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / joule',
        description='''
        DOS at Fermi energy
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_dos_fermi'))

    x_exciting_effective_potential_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Effective potential energy final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_effective_potential_energy'))

    x_exciting_electron_nuclear_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Electron-nuclear energy final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_electron_nuclear_energy'))

    x_exciting_exchange_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Exchange energy final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_exchange_energy'))

    x_exciting_fermi_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Fermi energy final
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_fermi_energy'))

    x_exciting_gap = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Estimated fundamental gap
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_gap'))

    x_exciting_geometry_atom_forces_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        x component of the force acting on the atom
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_atom_forces_x'))

    x_exciting_geometry_atom_forces_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        y component of the force acting on the atom
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_atom_forces_y'))

    x_exciting_geometry_atom_forces_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        z component of the force acting on the atom
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_atom_forces_z'))

    x_exciting_geometry_dummy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        time for scf in geometry optimization
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_dummy'))

    x_exciting_maximum_force_magnitude = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        Maximum force magnitude in geometry optimization
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_maximum_force_magnitude'))

    x_exciting_geometry_optimization_threshold_force = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        Value of threshold for the force modulus as convergence criterion of the
        geometry_optimization_method used in exciting
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_geometry_optimization_threshold_force'))

    x_exciting_hartree_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Hartree energy final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_hartree_energy'))

    x_exciting_interstitial_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Interstitial charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_interstitial_charge'))

    x_exciting_madelung_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Madelung energy final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_madelung_energy'))

    x_exciting_nuclear_nuclear_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Nuclear-nuclear energy final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_nuclear_nuclear_energy'))

    x_exciting_store_total_forces = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing converged atom forces cartesian
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_store_total_forces'))

    x_exciting_total_MT_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Total charge in muffin-tins
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_total_MT_charge'))

    x_exciting_XC_potential = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        XC potential final
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_XC_potential'))

    x_exciting_xs_bse_epsilon_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_energy_points'],
        unit='joule',
        description='''
        Energies of the dielectric function epsilon
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_epsilon_energies'))

    x_exciting_xs_bse_epsilon_im = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_energy_points'],
        description='''
        Imaginary part of the dielectric function epsilon
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_epsilon_im'))

    x_exciting_xs_bse_epsilon_re = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_energy_points'],
        description='''
        Real part of the dielectric function epsilon
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_epsilon_re'))

    x_exciting_xs_bse_exciton_amplitude_im = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_excitons'],
        description='''
        Imaginary part of the transition amplitude
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_exciton_amplitude_im'))

    x_exciting_xs_bse_exciton_amplitude_re = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_excitons'],
        description='''
        Real part of the transition amplitude
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_exciton_amplitude_re'))

    x_exciting_xs_bse_exciton_binding_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_excitons'],
        unit='joule',
        description='''
        Exciton binding energies
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_exciton_binding_energies'))

    x_exciting_xs_bse_exciton_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_excitons'],
        unit='joule',
        description='''
        Exciton energies
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_exciton_energies'))

    x_exciting_xs_bse_exciton_oscillator_strength = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_excitons'],
        description='''
        Exciton oscillator strength
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_exciton_oscillator_strength'))

    x_exciting_xs_bse_loss = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_energy_points'],
        description='''
        Loss function
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_loss'))

    x_exciting_xs_bse_loss_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_energy_points'],
        unit='joule',
        description='''
        Energies of the loss function
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_loss_energies'))

    x_exciting_xs_bse_number_of_components = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of independent components in the dielectric tensor
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_number_of_components'))

    x_exciting_xs_bse_number_of_energy_points = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Energy mesh for the dielectric function, conductivity and loss function
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_number_of_energy_points'))

    x_exciting_xs_bse_number_of_excitons = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Total number of excitons
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_number_of_excitons'))

    x_exciting_xs_bse_sigma_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_energy_points'],
        unit='joule',
        description='''
        Energies of the conductivity function sigma
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_sigma_energies'))

    x_exciting_xs_bse_sigma_im = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_energy_points'],
        description='''
        Imaginary part of the conductivity function sigma
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_sigma_im'))

    x_exciting_xs_bse_sigma_re = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_bse_number_of_components', 'x_exciting_xs_bse_number_of_energy_points'],
        description='''
        Real part of the conductivity function sigma
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_sigma_re'))

    x_exciting_xs_tddft_dielectric_function_local_field = Quantity(
        type=np.dtype(np.float64),
        shape=[2, 'x_exciting_xs_tddft_number_of_q_points', 'x_exciting_xs_tddft_number_of_components', 'x_exciting_xs_tddft_number_of_epsilon_values'],
        description='''
        Dielectric function epsilon including local-field effects
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_dielectric_function_local_field'))

    x_exciting_xs_tddft_dielectric_function_no_local_field = Quantity(
        type=np.dtype(np.float64),
        shape=[2, 'x_exciting_xs_tddft_number_of_q_points', 'x_exciting_xs_tddft_number_of_components', 'x_exciting_xs_tddft_number_of_epsilon_values'],
        description='''
        Dielectric function epsilon without local-field effects
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_dielectric_function_no_local_field'))

    x_exciting_xs_tddft_dielectric_tensor_no_sym = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_tddft_number_of_q_points', 2, 3, 3],
        description='''
        No-symmetrized static dielectric tensor
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_dielectric_tensor_no_sym'))

    x_exciting_xs_tddft_dielectric_tensor_sym = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_tddft_number_of_q_points', 2, 3, 3],
        description='''
        Symmetrized static dielectric tensor
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_dielectric_tensor_sym'))

    x_exciting_xs_tddft_epsilon_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_tddft_number_of_epsilon_values'],
        unit='joule',
        description='''
        Array containing the set of discrete energy values for the dielectric function
        epsilon.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_epsilon_energies'))

    x_exciting_xs_tddft_loss_function_local_field = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_tddft_number_of_q_points', 'x_exciting_xs_tddft_number_of_components', 'x_exciting_xs_tddft_number_of_epsilon_values'],
        description='''
        Loss function including local-field effects
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_loss_function_local_field'))

    x_exciting_xs_tddft_loss_function_no_local_field = Quantity(
        type=np.dtype(np.float64),
        shape=['x_exciting_xs_tddft_number_of_q_points', 'x_exciting_xs_tddft_number_of_components', 'x_exciting_xs_tddft_number_of_epsilon_values'],
        description='''
        Loss function without local-field effects
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_loss_function_no_local_field'))

    x_exciting_xs_tddft_number_of_epsilon_values = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of energy values for the dielectric function epsilon.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_number_of_epsilon_values'))

    x_exciting_xs_tddft_number_of_optical_components = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of independent components in the dielectric function epsilon
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_number_of_optical_components'))

    x_exciting_xs_tddft_number_of_q_points = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of Q points
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_number_of_q_points'))

    x_exciting_xs_tddft_optical_component = Quantity(
        type=str,
        shape=['x_exciting_xs_tddft_number_of_optical_components'],
        description='''
        Value of the independent optical components in the dielectric function epsilon
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_optical_component'))

    x_exciting_xs_tddft_sigma_local_field = Quantity(
        type=np.dtype(np.float64),
        shape=[2, 'x_exciting_xs_tddft_number_of_q_points', 'x_exciting_xs_tddft_number_of_components', 'x_exciting_xs_tddft_number_of_epsilon_values'],
        description='''
        Sigma including local-field effects
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_sigma_local_field'))

    x_exciting_xs_tddft_sigma_no_local_field = Quantity(
        type=np.dtype(np.float64),
        shape=[2, 'x_exciting_xs_tddft_number_of_q_points', 'x_exciting_xs_tddft_number_of_components', 'x_exciting_xs_tddft_number_of_epsilon_values'],
        description='''
        Sigma without local-field effects
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_sigma_no_local_field'))

    x_exciting_section_bandstructure = SubSection(
        sub_section=SectionProxy('x_exciting_section_bandstructure'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_bandstructure'))

    x_exciting_section_dos = SubSection(
        sub_section=SectionProxy('x_exciting_section_dos'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_dos'))

    x_exciting_section_fermi_surface = SubSection(
        sub_section=SectionProxy('x_exciting_section_fermi_surface'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_fermi_surface'))

    x_exciting_section_MT_charge_atom = SubSection(
        sub_section=SectionProxy('x_exciting_section_MT_charge_atom'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_MT_charge_atom')
    )

    x_exciting_section_MT_moment_atom = SubSection(
        sub_section=SectionProxy('x_exciting_section_MT_moment_atom'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_MT_moment_atom')
    )


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_system'))

    x_exciting_brillouin_zone_volume = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / meter ** 3',
        description='''
        Brillouin zone volume
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_brillouin_zone_volume'))

    x_exciting_core_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Core charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_core_charge'))

    x_exciting_core_charge_initial = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Core charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_core_charge_initial'))

    x_exciting_electronic_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Electronic charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_electronic_charge'))

    x_exciting_empty_states = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of empty states
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_empty_states'))

    x_exciting_clathrates_atom_labels = Quantity(
        type=str,
        shape=['number_of_atoms'],
        description='''
        Labels of the atoms in the clathrates.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_clathrates_atom_labels'))

    x_exciting_clathrates_atom_coordinates = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        Ordered list of the atoms coordinates in the clathrates.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_clathrates_atom_coordinates'))

    x_exciting_clathrates = Quantity(
        type=bool,
        shape=[],
        description='''
        It indicates whether the system is a clathrate.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_clathrates'))

    x_exciting_gkmax = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / meter',
        description='''
        Maximum length of |G+k| for APW functions
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_gkmax'))

    x_exciting_gmaxvr = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / meter',
        description='''
        Maximum length of |G|
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_gmaxvr'))

    x_exciting_gvector_size_x = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        G-vector grid size x
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_gvector_size_x'))

    x_exciting_gvector_size_y = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        G-vector grid size y
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_gvector_size_y'))

    x_exciting_gvector_size_z = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        G-vector grid size z
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_gvector_size_z'))

    x_exciting_gvector_size = Quantity(
        type=np.dtype(np.int32),
        shape=[3],
        description='''
        G-vector grid size
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_gvector_size'))

    x_exciting_gvector_total = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        G-vector total
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_gvector_total'))

    x_exciting_hamiltonian_size = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum Hamiltonian size
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_hamiltonian_size'))

    x_exciting_kpoint_offset_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        K-points offset x component
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_kpoint_offset_x'))

    x_exciting_kpoint_offset_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        K-points offset y component
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_kpoint_offset_y'))

    x_exciting_kpoint_offset_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        K-points offset z component
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_kpoint_offset_z'))

    x_exciting_kpoint_offset = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        K-points offset
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_kpoint_offset'))

    x_exciting_lmaxapw = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Angular momentum cut-off for the APW functions
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_lmaxapw'))

    x_exciting_lo = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Total number of local-orbitals
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_lo'))

    x_exciting_nuclear_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Nuclear charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_nuclear_charge'))

    x_exciting_number_kpoint_x = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number k-points x
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_number_kpoint_x'))

    x_exciting_number_kpoint_y = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number k-points y
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_number_kpoint_y'))

    x_exciting_number_kpoint_z = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number k-points z
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_number_kpoint_z'))

    x_exciting_kpoint_grid = Quantity(
        type=np.dtype(np.int32),
        shape=[3],
        description='''
        kpoint grid
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_kpoint_grid'))

    x_exciting_number_kpoints = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number k-points
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_number_kpoints'))

    x_exciting_number_of_atoms = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The number of atoms in the unit cell
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_number_of_atoms'))

    x_exciting_potential_mixing = Quantity(
        type=str,
        shape=[],
        description='''
        Mixing type for potential
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_potential_mixing'))

    x_exciting_pw = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum number of plane-waves
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_pw'))

    x_exciting_rgkmax = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        Radius MT * Gmax
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_rgkmax'))

    x_exciting_species_rtmin = Quantity(
        type=str,
        shape=[],
        description='''
        Chemical species with radius RT * Gmax
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_species_rtmin'))

    x_exciting_simulation_reciprocal_cell = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='meter',
        description='''
        Reciprocal lattice vectors (in Cartesian coordinates) of the simulation cell. The
        first index runs over the $x,y,z$ Cartesian coordinates, and the second index runs
        over the 3 lattice vectors.
        ''',
        categories=[public.configuration_core],
        a_legacy=LegacyDefinition(name='x_exciting_simulation_reciprocal_cell'))

    x_exciting_smearing_type = Quantity(
        type=str,
        shape=[],
        description='''
        Smearing scheme for KS occupancies
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_smearing_type'))

    x_exciting_smearing_width = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Smearing width for KS occupancies
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_smearing_width'))

    x_exciting_unit_cell_volume = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter ** 3',
        description='''
        unit cell volume
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_unit_cell_volume'))

    x_exciting_valence_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Valence charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_valence_charge'))

    x_exciting_valence_charge_initial = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Valence charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_valence_charge_initial'))

    x_exciting_valence_states = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Total number of valence states
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_valence_states'))

    x_exciting_wigner_radius = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        Effective Wigner radius
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_wigner_radius'))

    x_exciting_number_of_bravais_lattice_symmetries = Quantity(
        type=int,
        shape=[],
        description='''
        Number of Bravais lattice symmetries
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_number_of_bravais_lattice_symmetries'))

    x_exciting_number_of_crystal_symmetries = Quantity(
        type=int,
        shape=[],
        description='''
        Number of crystal symmetries
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_number_of_crystal_symmetries'))

    x_exciting_section_atoms_group = SubSection(
        sub_section=SectionProxy('x_exciting_section_atoms_group'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_atoms_group'))

    x_exciting_section_lattice_vectors = SubSection(
        sub_section=SectionProxy('x_exciting_section_lattice_vectors'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_lattice_vectors'))

    x_exciting_section_reciprocal_lattice_vectors = SubSection(
        sub_section=SectionProxy('x_exciting_section_reciprocal_lattice_vectors'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_reciprocal_lattice_vectors'))

    x_exciting_section_spin = SubSection(
        sub_section=SectionProxy('x_exciting_section_spin'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_spin'))

    x_exciting_section_xc = SubSection(
        sub_section=SectionProxy('x_exciting_section_xc'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_xc'))


class x_exciting_section_MT_charge_atom(MSection):
    '''
    atom-resolved charges in muffin tins
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_MT_charge_atom'))

    x_exciting_MT_charge_atom_index = Quantity(
        type=int,
        shape=[],
        description='''
        index of the atom with muffin-tin charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_MT_charge_atom_index'))

    x_exciting_MT_charge_atom_symbol = Quantity(
        type=str,
        shape=[],
        description='''
        chemical symbol of the atom with muffin-tin charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_MT_charge_atom_symbol'))

    x_exciting_MT_charge_atom_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        value of the muffin-tin charge on atom
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_MT_charge_atom_value'))


class x_exciting_section_MT_moment_atom(MSection):
    '''
    atom-resolved moments in muffin tins
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_exciting_section_MT_moment_atom'))

    x_exciting_MT_moment_atom_index = Quantity(
        type=int,
        shape=[],
        description='''
        index of the atom with muffin-tin moment
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_MT_moment_atom_index'))

    x_exciting_MT_moment_atom_symbol = Quantity(
        type=str,
        shape=[],
        description='''
        chemical symbol of the atom with muffin-tin moment
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_MT_moment_atom_symbol'))

    x_exciting_MT_moment_atom_value = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        unit='coulomb * meter',
        description='''
        value of the muffin-tin moment on atom
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_MT_moment_atom_value'))


class section_scf_iteration(public.section_scf_iteration):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_scf_iteration'))

    x_exciting_charge_convergence_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        exciting charge convergence
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_charge_convergence_scf_iteration'))

    x_exciting_core_electron_kinetic_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Core-electron kinetic energy
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_core_electron_kinetic_energy_scf_iteration'))

    x_exciting_core_charge_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Core charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_core_charge_scf_iteration'))

    x_exciting_valence_charge_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Valence charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_valence_charge_scf_iteration'))

    x_exciting_core_leakage_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Core leakage
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_core_leakage_scf_iteration'))

    x_exciting_correlation_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Correlation energy
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_correlation_energy_scf_iteration'))

    x_exciting_coulomb_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Coulomb energy
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_coulomb_energy_scf_iteration'))

    x_exciting_coulomb_potential_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Coulomb potential energy
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_coulomb_potential_energy_scf_iteration'))

    x_exciting_dos_fermi_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / joule',
        description='''
        DOS at Fermi energy
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_dos_fermi_scf_iteration'))

    x_exciting_effective_potential_convergence_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        exciting effective potential convergence
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_effective_potential_convergence_scf_iteration'))

    x_exciting_force_convergence_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        exciting force convergence
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_force_convergence_scf_iteration'))

    x_exciting_effective_potential_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Effective potential energy
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_effective_potential_energy_scf_iteration'))

    x_exciting_electron_nuclear_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Electron-nuclear energy
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_electron_nuclear_energy_scf_iteration'))

    x_exciting_energy_convergence_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        exciting energy convergence
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_energy_convergence_scf_iteration'))

    x_exciting_exchange_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Exchange energy
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_exchange_energy_scf_iteration'))

    x_exciting_fermi_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Fermi energy
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_fermi_energy_scf_iteration'))

    x_exciting_gap_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Estimated fundamental gap
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_gap_scf_iteration'))

    x_exciting_hartree_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Hartree energy
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_hartree_energy_scf_iteration'))

    x_exciting_IBS_force_convergence_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        exciting IBS force convergence
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_IBS_force_convergence_scf_iteration'))

    x_exciting_interstitial_charge_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Interstitial charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_interstitial_charge_scf_iteration'))

    x_exciting_madelung_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Madelung energy
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_madelung_energy_scf_iteration'))

    x_exciting_nuclear_nuclear_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Nuclear-nuclear energy
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_nuclear_nuclear_energy_scf_iteration'))

    x_exciting_time_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='s',
        description='''
        scf iteration time
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_time_scf_iteration'))

    x_exciting_total_MT_charge_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Total charge in muffin-tins
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_total_MT_charge_scf_iteration'))

    x_exciting_total_charge_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Total charge
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_total_charge_scf_iteration'))

    x_exciting_XC_potential_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        XC potential
        ''',
        categories=[public.energy_component, public.energy_value],
        a_legacy=LegacyDefinition(name='x_exciting_XC_potential_scf_iteration'))

    x_exciting_interstitial_moment_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        unit='coulomb * meter',
        description='''
        Interstitial moment
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_interstitial_moment_scf_iteration'))

    x_exciting_total_MT_moment_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        unit='coulomb * meter',
        description='''
        Total moment in muffin-tins
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_total_MT_moment_scf_iteration'))

    x_exciting_total_moment_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        unit='coulomb * meter',
        description='''
        Total moment
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_total_moment_scf_iteration'))

    x_exciting_section_MT_charge_atom_scf_iteration = SubSection(
        sub_section=SectionProxy('x_exciting_section_MT_charge_atom'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_MT_charge_atom_scf_iteration'))

    x_exciting_section_MT_moment_atom_scf_iteration = SubSection(
        sub_section=SectionProxy('x_exciting_section_MT_moment_atom'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_MT_moment_atom_scf_iteration'))


class section_method(public.section_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_method'))

    x_exciting_dummy = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        dummy metadata for debuging purposes
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_dummy'))

    x_exciting_volume_optimization = Quantity(
        type=bool,
        shape=[],
        description='''
        If the volume optimization is performed.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_volume_optimization'))

    x_exciting_scf_threshold_energy_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Specifies the threshold for the x_exciting_energy_total_scf_iteration change
        between two subsequent self-consistent field (SCF) iterations.
        ''',
        categories=[public.settings_scf],
        a_legacy=LegacyDefinition(name='x_exciting_scf_threshold_energy_change'))

    x_exciting_scf_threshold_potential_change_list = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Specifies the threshold for the x_exciting_effective_potential_convergence between
        two subsequent self-consistent field (SCF) iterations.
        ''',
        categories=[public.settings_scf],
        a_legacy=LegacyDefinition(name='x_exciting_scf_threshold_potential_change_list'))

    x_exciting_scf_threshold_potential_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Specifies the threshold for the x_exciting_effective_potential_convergence between
        two subsequent self-consistent field (SCF) iterations.
        ''',
        categories=[public.settings_scf],
        a_legacy=LegacyDefinition(name='x_exciting_scf_threshold_potential_change'))

    x_exciting_scf_threshold_charge_change_list = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Specifies the threshold for the x_exciting_effective_potential_convergence between
        two subsequent self-consistent field (SCF) iterations.
        ''',
        categories=[public.settings_scf],
        a_legacy=LegacyDefinition(name='x_exciting_scf_threshold_charge_change_list'))

    x_exciting_scf_threshold_charge_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Specifies the threshold for the x_exciting_effective_potential_convergence between
        two subsequent self-consistent field (SCF) iterations.
        ''',
        categories=[public.settings_scf],
        a_legacy=LegacyDefinition(name='x_exciting_scf_threshold_charge_change'))

    x_exciting_scf_threshold_force_change_list = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Convergence tolerance for forces (not including IBS contribution) during the SCF
        run.
        ''',
        categories=[public.settings_scf],
        a_legacy=LegacyDefinition(name='x_exciting_scf_threshold_force_change_list'))

    x_exciting_scf_threshold_force_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        Convergence tolerance for forces (not including IBS contribution) during the SCF
        run
        ''',
        categories=[public.settings_scf],
        a_legacy=LegacyDefinition(name='x_exciting_scf_threshold_force_change'))

    x_exciting_electronic_structure_method = Quantity(
        type=str,
        shape=[],
        description='''
        String identifying the used electronic structure method. Allowed values are BSE
        and TDDFT. Temporary metadata to be canceled when BSE and TDDFT will be added to
        public metadata electronic_structure_method
        ''',
        categories=[public.settings_potential_energy_surface, public.settings_XC],
        a_legacy=LegacyDefinition(name='x_exciting_electronic_structure_method'))

    x_exciting_xs_broadening = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Lorentzian broadening applied to the spectra.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_broadening'))

    x_exciting_xs_bse_angular_momentum_cutoff = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Angular momentum cutoff of the spherical harmonics expansion of the dielectric
        matrix.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_angular_momentum_cutoff'))

    x_exciting_xs_bse_antiresonant = Quantity(
        type=bool,
        shape=[],
        description='''
        If the anti-resonant part of the Hamiltonian is used for the BSE spectrum
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_antiresonant'))

    x_exciting_xs_bse_number_of_bands = Quantity(
        type=np.dtype(np.int32),
        shape=[4],
        description='''
        Range of bands included for the BSE calculation. The first pair of numbers
        corresponds to the band index for local orbitals and valence states (counted from
        the lowest eigenenergy), the second pair corresponds to the band index of the
        conduction states (counted from the Fermi level).
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_number_of_bands'))

    x_exciting_xs_bse_rgkmax = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Smallest muffin-tin radius times gkmax.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_rgkmax'))

    x_exciting_xs_bse_sciavbd = Quantity(
        type=bool,
        shape=[],
        description='''
        if the body of the screened Coulomb interaction is averaged (q=0).
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_sciavbd'))

    x_exciting_xs_bse_sciavqbd = Quantity(
        type=bool,
        shape=[],
        description='''
        if the body of the screened Coulomb interaction is averaged (q!=0).
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_sciavqbd'))

    x_exciting_xs_bse_sciavqhd = Quantity(
        type=bool,
        shape=[],
        description='''
        if the head of the screened Coulomb interaction is averaged (q=0).
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_sciavqhd'))

    x_exciting_xs_bse_sciavqwg = Quantity(
        type=bool,
        shape=[],
        description='''
        if the wings of the screened Coulomb interaction is averaged (q=0).
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_sciavqwg'))

    x_exciting_xs_bse_sciavtype = Quantity(
        type=str,
        shape=[],
        description='''
        how the screened Coulomb interaction matrix is averaged
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_sciavtype'))

    x_exciting_xs_bse_type = Quantity(
        type=str,
        shape=[],
        description='''
        which parts of the BSE Hamiltonian is considered. Possible values are IP, RPA,
        singlet, triplet.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_type'))

    x_exciting_xs_bse_xas = Quantity(
        type=bool,
        shape=[],
        description='''
        X-ray absorption spectrum
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_xas'))

    x_exciting_xs_bse_xas_number_of_bands = Quantity(
        type=np.dtype(np.int32),
        shape=[2],
        description='''
        Range of bands included for the BSE calculation. The pair corresponds to the band
        index of the conduction states (counted from the Fermi level).
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_xas_number_of_bands'))

    x_exciting_xs_bse_xasatom = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Atom number for which the XAS is calculated.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_xasatom'))

    x_exciting_xs_bse_xasedge = Quantity(
        type=str,
        shape=[],
        description='''
        Defines the initial states of the XAS calculation.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_xasedge'))

    x_exciting_xs_bse_xasspecies = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Species number for which the XAS is calculated..
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_bse_xasspecies'))

    x_exciting_xs_gqmax = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / meter',
        description='''
        G vector cutoff for the plane-waves matrix elements.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_gqmax'))

    x_exciting_xs_lmaxapw = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Angular momentum cut-off for the APW functions.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_lmaxapw'))

    x_exciting_xs_ngridk = Quantity(
        type=np.dtype(np.int32),
        shape=[3],
        description='''
        k-point mesh for the excited states calculation.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_ngridk'))

    x_exciting_xs_ngridq = Quantity(
        type=np.dtype(np.int32),
        shape=[3],
        description='''
        q-point mesh for the excited states calculation.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_ngridq'))

    x_exciting_xs_number_of_empty_states = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of empty states used in the calculation of the response function.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_number_of_empty_states'))

    x_exciting_xs_rgkmax = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Smallest muffin-tin radius times gkmax.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_rgkmax'))

    x_exciting_xs_scissor = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Scissor operator
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_scissor'))

    x_exciting_xs_screening_ngridk = Quantity(
        type=np.dtype(np.int32),
        shape=[3],
        description='''
        k-point mesh for the calculation of the screening.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_screening_ngridk'))

    x_exciting_xs_screening_number_of_empty_states = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of empty states used in the calculation of the screening.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_screening_number_of_empty_states'))

    x_exciting_xs_screening_rgkmax = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Smallest muffin-tin radius times gkmax.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_screening_rgkmax'))

    x_exciting_xs_screening_type = Quantity(
        type=str,
        shape=[],
        description='''
        Type of screening used in the calculations. POssible values are full, diag,
        noinvdiag, longrange.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_screening_type'))

    x_exciting_xs_screening_vkloff = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        k-point set offset for the screening
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_screening_vkloff'))

    x_exciting_xs_starting_point = Quantity(
        type=str,
        shape=[],
        description='''
        Exchange-correlation functional of the ground-state calculation. See xc_functional
        list at https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-
        functional
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_starting_point'))

    x_exciting_xs_tddft_analytic_continuation = Quantity(
        type=bool,
        shape=[],
        description='''
        Analytic continuation for the calculation of the Kohn-Sham response function
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_analytic_continuation'))

    x_exciting_xs_tddft_analytic_continuation_number_of_intervals = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of energy intervals (on imaginary axis) for analytic continuation.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_analytic_continuation_number_of_intervals'))

    x_exciting_xs_tddft_anomalous_hall_conductivity = Quantity(
        type=bool,
        shape=[],
        description='''
        If the anomalous Hall conductivity term is included in the calculation of the
        dielectric tensor [see PRB 86, 125139 (2012)].
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_anomalous_hall_conductivity'))

    x_exciting_xs_tddft_anti_resonant_dielectric = Quantity(
        type=bool,
        shape=[],
        description='''
        If the anti-resonant part is considered for the calculation of the dielectric
        function
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_anti_resonant_dielectric'))

    x_exciting_xs_tddft_anti_resonant_xc_kernel = Quantity(
        type=bool,
        shape=[],
        description='''
        If the anti-resonant part is considered for the calculation of the MBPT-derived
        xc-kernel
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_anti_resonant_xc_kernel'))

    x_exciting_xs_tddft_diagonal_xc_kernel = Quantity(
        type=bool,
        shape=[],
        description='''
        If true, only the diagonal part of xc-kernel is used.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_diagonal_xc_kernel'))

    x_exciting_xs_tddft_drude = Quantity(
        type=np.dtype(np.float64),
        shape=[2],
        unit='1 / second',
        description='''
        Parameters defining the semiclassical Drude approximation to intraband term. The
        first value determines the plasma frequency ωp and the second the inverse
        relaxation time ωτ: χD0=1/ω ωp^2/(ω+iωτ)
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_drude'))

    x_exciting_xs_tddft_finite_q_intraband_contribution = Quantity(
        type=bool,
        shape=[],
        description='''
        Parameter determining whether the the intraband contribution is included in the
        calculation for the finite q.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_finite_q_intraband_contribution'))

    x_exciting_xs_tddft_lmax_alda = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Angular momentum cutoff for the Rayleigh expansion of the exponential factor for
        ALDA-kernel
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_lmax_alda'))

    x_exciting_xs_tddft_macroscopic_dielectric_function_q_treatment = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Treatment of macroscopic dielectric function for the Q-point outside the Brillouin
        zone. A value of 0 uses the full Q and the (0,0) component of the microscopic
        dielectric matrix is used. A value of 1 means a decomposition Q=q+Gq and the
        (Qq,Qq) component of the microscopic dielectric matrix is used.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_macroscopic_dielectric_function_q_treatment'))

    x_exciting_xs_tddft_split_parameter = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Split parameter for degeneracy in energy differences of MBPT derived xc kernels.
        See A. Marini, Phys. Rev. Lett., 91, (2003) 256402.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_split_parameter'))

    x_exciting_xs_tddft_xc_kernel = Quantity(
        type=str,
        shape=[],
        description='''
        Defines which xc kernel is to be used. Possible choices are: RPA - Random-phase
        approximation kernel (fxc=0); LRCstatic - Long-range correction kernel (fxc =
        -alpha/q^2) with alpha given by alphalrcdyn see S. Botti et al., Phys. Rev. B 69,
        155112 (2004); LRCdyn - Dynamical long-range correction kernel, with alpha anf
        beta give by alphalrcdyn and betalrcdyn, respectively, see S. Botti et al., Phys.
        Rev. B 72, 125203 (2005); ALDA - Adiabatic LDA kernel, with Vxc given by the spin-
        unpolarised exchange-correlation potential corresponding to the Perdew-Wang
        parameterisation of Ceperley-Alder's Monte-Carlo data, see Phys. Rev. B 45, 13244
        (1992) and Phys. Rev. Lett. 45, 566 (1980); MB1 - BSE derived xc kernel. See L.
        Reining et al., Phys. Rev. Lett. 88, 066404 (2002) and A. Marini et al., Phys.
        Rev. Lett. 91, 256402 (2003); BO - Bootstrap kernel, see S. Sharma et al., Phys.
        Rev. Lett. 107, 186401 (2011); BO_SCALAR - Scalar version of the bootstrap kernel;
        see S. Sharma et al., Phys. Rev. Lett. 107, 186401 (2011); RBO - RPA bootstrap
        kernel; see S. Rigamonti et al., Phys. Rev. Lett. 114, 146402 (2015).
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tddft_xc_kernel'))

    x_exciting_xs_tetra = Quantity(
        type=bool,
        shape=[],
        description='''
        Integration method (tetrahedron) used for the k-space integration in the Kohn-Sham
        response function.
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_tetra'))

    x_exciting_xs_vkloff = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        k-point set offset
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_vkloff'))

    x_exciting_xs_xstype = Quantity(
        type=str,
        shape=[],
        description='''
        Type of excited states calculation, BSE or TDDFT
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_xs_xstype'))


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_run'))

    x_exciting_dummy2 = Quantity(
        type=str,
        shape=[],
        description='''
        dummy metadata for debuging purposes
        ''',
        a_legacy=LegacyDefinition(name='x_exciting_dummy2'))

    x_exciting_section_geometry_optimization = SubSection(
        sub_section=SectionProxy('x_exciting_section_geometry_optimization'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_exciting_section_geometry_optimization'))


m_package.__init_metainfo__()

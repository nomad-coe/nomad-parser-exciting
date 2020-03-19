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
        categories=[public.energy_value, public.energy_component],
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
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_correlation_energy'))

    x_exciting_coulomb_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Coulomb energy final
        ''',
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_coulomb_energy'))

    x_exciting_coulomb_potential_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Coulomb potential energy final
        ''',
        categories=[public.energy_value, public.energy_component],
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
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_effective_potential_energy'))

    x_exciting_electron_nuclear_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Electron-nuclear energy final
        ''',
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_electron_nuclear_energy'))

    x_exciting_exchange_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Exchange energy final
        ''',
        categories=[public.energy_value, public.energy_component],
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
        categories=[public.energy_value, public.energy_component],
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
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_madelung_energy'))

    x_exciting_nuclear_nuclear_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Nuclear-nuclear energy final
        ''',
        categories=[public.energy_value, public.energy_component],
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
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_XC_potential'))

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


class section_scf_iteration(public.section_scf_iteration):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_scf_iteration'))

    x_exciting_charge_convergence_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
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
        categories=[public.energy_value, public.energy_component],
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
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_correlation_energy_scf_iteration'))

    x_exciting_coulomb_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Coulomb energy
        ''',
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_coulomb_energy_scf_iteration'))

    x_exciting_coulomb_potential_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Coulomb potential energy
        ''',
        categories=[public.energy_value, public.energy_component],
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
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_effective_potential_energy_scf_iteration'))

    x_exciting_electron_nuclear_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Electron-nuclear energy
        ''',
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_electron_nuclear_energy_scf_iteration'))

    x_exciting_energy_convergence_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
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
        categories=[public.energy_value, public.energy_component],
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
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_hartree_energy_scf_iteration'))

    x_exciting_IBS_force_convergence_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
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
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_madelung_energy_scf_iteration'))

    x_exciting_nuclear_nuclear_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Nuclear-nuclear energy
        ''',
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_nuclear_nuclear_energy_scf_iteration'))

    x_exciting_time_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
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

    x_exciting_XC_potential_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        XC potential
        ''',
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_exciting_XC_potential_scf_iteration'))


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
        description='''
        Convergence tolerance for forces (not including IBS contribution) during the SCF
        run
        ''',
        categories=[public.settings_scf],
        a_legacy=LegacyDefinition(name='x_exciting_scf_threshold_force_change'))


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

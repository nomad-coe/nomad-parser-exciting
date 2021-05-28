#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
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

import pytest
import numpy as np

from nomad.units import ureg
from nomad.datamodel import EntryArchive
from excitingparser.exciting_parser import ExcitingParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return ExcitingParser()


@pytest.fixture(scope='module')
def silicon_gw(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Si_gw/INFO.OUT', archive, None)
    return archive


def test_gs(parser):
    archive = EntryArchive()
    parser.parse('tests/data/C_gs/INFO.OUT', archive, None)

    assert len(archive.section_run) == 1

    sec_run = archive.section_run[0]
    assert sec_run.program_version == 'CARBON'

    sec_method = sec_run.section_method[0]
    assert sec_method.number_of_spin_channels == 1
    assert sec_method.smearing_width == approx(4.35974472e-22)
    assert sec_method.section_XC_functionals[1].XC_functional_name == 'GGA_X_PBE_SOL'
    assert sec_method.x_exciting_scf_threshold_force_change.magnitude == approx(4.11936175e-12)

    sec_system = sec_run.section_system[0]
    assert sec_system.lattice_vectors[0][0].magnitude == approx(1.72297146e-10)
    assert sec_system.atom_positions[1][0].magnitude == approx(8.61485729e-11)
    assert len(sec_system.atom_labels) == 2
    assert sec_system.x_exciting_section_spin[0].x_exciting_spin_treatment == 'spin-unpolarised'
    assert sec_system.x_exciting_section_atoms_group[0].x_exciting_muffin_tin_radius.magnitude == approx(6.87930374e-11)

    sec_scc = sec_run.section_single_configuration_calculation[0]
    assert sec_scc.energy_total.magnitude == approx(-3.30863556e-16)
    assert np.mean(sec_scc.atom_forces) == 0.0
    assert sec_scc.charge_total.magnitude == approx(1.92261196e-18)
    assert sec_scc.energy_reference_fermi.magnitude == approx(2.4422694e-18)
    assert len(sec_scc.section_scf_iteration) == 12
    assert sec_scc.section_scf_iteration[5].x_exciting_valence_charge_scf_iteration.magnitude == approx(1.28174131e-18)
    assert sec_scc.section_scf_iteration[8].x_exciting_exchange_energy_scf_iteration.magnitude == approx(-4.39756926e-17)
    assert sec_scc.section_scf_iteration[11].electronic_kinetic_energy_scf_iteration.magnitude == approx(3.30404896e-16)
    sec_eig = sec_scc.section_eigenvalues[0]
    assert np.shape(sec_eig.eigenvalues_kpoints) == (30, 3)
    assert sec_eig.eigenvalues_values[0][9][4].magnitude == approx(2.74680139e-18)


def test_strucopt(parser):
    archive = EntryArchive()
    parser.parse('tests/data/GaO_strucopt/INFO.OUT', archive, None)

    sec_systems = archive.section_run[0].section_system
    assert len(sec_systems) == 15
    assert sec_systems[0].atom_labels == ['Ga', 'Ga', 'Ga', 'Ga', 'O', 'O', 'O', 'O', 'O', 'O']
    assert sec_systems[0].x_exciting_gkmax.magnitude == approx(1.13383567e+11)
    assert sec_systems[3].atom_positions[1][1].magnitude == approx(3.07695918e-10)
    assert sec_systems[10].atom_positions[-1][0].magnitude == approx(3.67156876e-11)
    assert sec_systems[1].lattice_vectors[2][1].magnitude == approx(sec_systems[13].lattice_vectors[2][1].magnitude)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 15
    assert len(sec_sccs[0].section_scf_iteration) == 19
    assert sec_sccs[0].section_scf_iteration[10].time_scf_iteration.magnitude == approx(431.84)
    assert sec_sccs[0].section_scf_iteration[18].x_exciting_effective_potential_convergence_scf_iteration[0].magnitude == approx(4.62350928e-26)
    assert sec_sccs[3].x_exciting_maximum_force_magnitude.magnitude == approx(1.64771998e-10)
    assert sec_sccs[6].energy_total.magnitude == approx(-3.58415586e-14)
    assert sec_sccs[9].time_calculation.magnitude == approx(724.33)
    assert len(sec_sccs[-1].x_exciting_section_MT_charge_atom) == 10
    assert sec_sccs[-1].x_exciting_fermi_energy.magnitude == approx(1.03200886e-18)


def test_dos_spinpol(parser):
    archive = EntryArchive()
    parser.parse('tests/data/CeO_dos/INFO.OUT', archive, None)

    sec_dos = archive.section_run[0].section_single_configuration_calculation[0].dos_electronic[0]

    assert np.shape(sec_dos.dos_total[0].dos_values) == (500,)
    assert sec_dos.dos_energies[79].magnitude == approx(-1.70772016e-18)
    assert sec_dos.dos_total[0].dos_values[126] == approx(1.8925127e-10)
    assert sec_dos.dos_total[1].dos_values[136] == approx(1.91606129e-11)
    assert sec_dos.dos_energies[240].magnitude == approx(1.09995544e-18)
    assert sec_dos.dos_total[0].dos_values[220] == approx(5.63875821e-10)
    assert sec_dos.dos_total[1].dos_values[78] == approx(4.33359121e-10)

    assert len(sec_dos.dos_atom_projected) == 150
    assert np.shape(sec_dos.dos_atom_projected[149].dos_values) == (500,)
    assert sec_dos.dos_atom_projected[7].dos_values[116] == approx(1.07677456e+16)
    assert sec_dos.dos_atom_projected[123].dos_values[85] == approx(7.81205293e+11)


def test_xs_tddft(parser):
    archive = EntryArchive()
    parser.parse('tests/data/CSi_tddft/INFO_QMT001.OUT', archive, None)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 2

    assert len(sec_sccs[1].x_exciting_xs_tddft_epsilon_energies) == 10001
    assert np.shape(sec_sccs[1].x_exciting_xs_tddft_dielectric_function_local_field) == (2, 1, 3, 10001)
    assert np.shape(sec_sccs[1].x_exciting_xs_tddft_dielectric_function_no_local_field) == (2, 1, 3, 10001)
    assert np.shape(sec_sccs[1].x_exciting_xs_tddft_loss_function_local_field) == (1, 3, 10001)
    assert np.shape(sec_sccs[1].x_exciting_xs_tddft_loss_function_no_local_field) == (1, 3, 10001)
    assert np.shape(sec_sccs[1].x_exciting_xs_tddft_loss_function_no_local_field) == (1, 3, 10001)
    assert np.shape(sec_sccs[1].x_exciting_xs_tddft_sigma_local_field) == (2, 1, 3, 10001)
    assert np.shape(sec_sccs[1].x_exciting_xs_tddft_sigma_no_local_field) == (2, 1, 3, 10001)


def test_xs_bse(parser):
    archive = EntryArchive()
    parser.parse('tests/data/CHN_bse/INFO.OUT', archive, None)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    # gs + ip + singlet + triplet
    assert len(sec_sccs) == 4
    assert np.shape(sec_sccs[1].x_exciting_xs_bse_epsilon_energies) == (3, 10000)
    assert np.shape(sec_sccs[2].x_exciting_xs_bse_epsilon_energies) == (3, 10000)
    assert np.shape(sec_sccs[3].x_exciting_xs_bse_epsilon_energies) == (3, 10000)
    assert np.shape(sec_sccs[1].x_exciting_xs_bse_epsilon_im) == (3, 10000)
    assert np.shape(sec_sccs[2].x_exciting_xs_bse_epsilon_im) == (3, 10000)
    assert np.shape(sec_sccs[3].x_exciting_xs_bse_epsilon_im) == (3, 10000)
    assert np.shape(sec_sccs[1].x_exciting_xs_bse_epsilon_re) == (3, 10000)
    assert np.shape(sec_sccs[2].x_exciting_xs_bse_epsilon_re) == (3, 10000)
    assert np.shape(sec_sccs[3].x_exciting_xs_bse_epsilon_re) == (3, 10000)
    assert np.shape(sec_sccs[1].x_exciting_xs_bse_exciton_amplitude_im) == (3, 100)
    assert np.shape(sec_sccs[2].x_exciting_xs_bse_exciton_amplitude_im) == (3, 100)
    assert np.shape(sec_sccs[3].x_exciting_xs_bse_exciton_amplitude_im) == (3, 100)
    assert np.shape(sec_sccs[1].x_exciting_xs_bse_exciton_amplitude_re) == (3, 100)
    assert np.shape(sec_sccs[2].x_exciting_xs_bse_exciton_amplitude_re) == (3, 100)
    assert np.shape(sec_sccs[3].x_exciting_xs_bse_exciton_amplitude_re) == (3, 100)
    assert np.shape(sec_sccs[1].x_exciting_xs_bse_exciton_binding_energies) == (3, 100)
    assert np.shape(sec_sccs[2].x_exciting_xs_bse_exciton_binding_energies) == (3, 100)
    assert np.shape(sec_sccs[3].x_exciting_xs_bse_exciton_binding_energies) == (3, 100)
    assert np.shape(sec_sccs[1].x_exciting_xs_bse_exciton_energies) == (3, 100)
    assert np.shape(sec_sccs[2].x_exciting_xs_bse_exciton_energies) == (3, 100)
    assert np.shape(sec_sccs[3].x_exciting_xs_bse_exciton_energies) == (3, 100)
    assert np.shape(sec_sccs[1].x_exciting_xs_bse_exciton_oscillator_strength) == (3, 100)
    assert np.shape(sec_sccs[2].x_exciting_xs_bse_exciton_oscillator_strength) == (3, 100)
    assert np.shape(sec_sccs[3].x_exciting_xs_bse_exciton_oscillator_strength) == (3, 100)


def test_gw(silicon_gw):
    """Basic tests for a GW calculation."""
    sec_methods = silicon_gw.section_run[0].section_method
    assert len(sec_methods) == 2
    assert sec_methods[1].electronic_structure_method == 'G0W0'
    assert sec_methods[1].gw_mixed_basis_gmax.magnitude == approx(226767134954.67346)
    assert sec_methods[1].gw_self_energy_singularity_treatment == 'mpb'
    assert sec_methods[1].gw_number_of_frequencies == 32
    assert sec_methods[1].gw_frequency_values[-1].magnitude == approx(8.22665908e-16)

    sec_sccs = silicon_gw.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 2

    # Check GW properties
    assert approx(sec_sccs[1].gw_fermi_energy.magnitude, 1.09865567e-19)
    assert approx(sec_sccs[1].gw_fundamental_gap.magnitude, 3.42913865e-19)
    assert approx(sec_sccs[1].gw_optical_gap.magnitude, 6.45981597e-19)
    assert np.shape(sec_sccs[1].section_eigenvalues[0].eigenvalues_values) == (1, 3, 20)
    assert sec_sccs[1].section_eigenvalues[0].eigenvalues_kpoints[-3][1] == 0.0
    assert sec_sccs[1].section_eigenvalues[0].eigenvalues_values[0][2][9].magnitude == approx(1.769533187849446e-18, abs=1e-20)
    assert sec_sccs[1].section_eigenvalues[0].gw_qp_linearization_prefactor[0][2][9] == approx(0.79935)
    assert sec_sccs[1].gw_self_energy_x[0][2][0].magnitude == approx(-2.855981572623473e-18, abs=1e-20)
    assert sec_sccs[1].gw_self_energy_c[0][2][14].magnitude == approx(-1.0879742954267992e-18, abs=1e-20)
    assert sec_sccs[1].gw_xc_potential[0][2][6].magnitude == approx(-2.1691473890869554e-18, abs=1e-20)


def test_band_gw_silicon(silicon_gw):
    """Tests that the band structure of silicon is parsed correctly from a GW
    calculation.
    """
    sccs = silicon_gw.section_run[-1].section_single_configuration_calculation
    assert len(sccs) == 2
    gaps = [0.446307, 1.2553776]
    for gap_assumed, scc in zip(gaps, sccs):
        assert len(scc.section_k_band) == 1
        band = scc.section_k_band[0]
        segments = band.section_k_band_segment
        energies = [s.band_energies.to(ureg.electron_volt).magnitude for s in segments]
        energies = np.concatenate(energies, axis=1)

        # Check that an energy reference is reported
        energy_reference = scc.energy_reference_fermi
        if energy_reference is None:
            energy_reference = scc.energy_reference_highest_occupied
        assert energy_reference is not None
        energy_reference = energy_reference.to(ureg.electron_volt).magnitude

        # Check that an appropriately sized band gap is found at the given
        # reference energy
        energies = energies.flatten()
        energies.sort()
        lowest_unoccupied_index = np.searchsorted(energies, energy_reference, "right")[0]
        highest_occupied_index = lowest_unoccupied_index - 1
        gap = energies[lowest_unoccupied_index] - energies[highest_occupied_index]
        assert gap == approx(gap_assumed)


def test_dos_gw_silicon(silicon_gw):
    """Tests that the DOS of silicon is parsed correctly from a GW calculation.
    """
    sccs = silicon_gw.section_run[-1].section_single_configuration_calculation
    assert len(sccs) == 2
    gaps = [0.5442277, 1.360569]
    for gap_assumed, scc in zip(gaps, sccs):
        assert len(scc.dos_electronic) == 1
        dos = scc.dos_electronic[0]
        energies = dos.dos_energies.to(ureg.electron_volt).magnitude
        values = np.array([d.dos_values for d in dos.dos_total])

        # Check that an energy reference is reported
        energy_reference = scc.energy_reference_fermi
        if energy_reference is None:
            energy_reference = scc.energy_reference_highest_occupied
        assert energy_reference is not None
        energy_reference = energy_reference.to(ureg.electron_volt).magnitude

        # Check that an appropriately sized band gap is found at the given
        # reference energy
        nonzero = np.unique(values.nonzero())
        energies = energies[nonzero]
        energies.sort()
        lowest_unoccupied_index = np.searchsorted(energies, energy_reference, "right")[0]
        highest_occupied_index = lowest_unoccupied_index - 1
        gap = energies[lowest_unoccupied_index] - energies[highest_occupied_index]
        assert gap == approx(gap_assumed)

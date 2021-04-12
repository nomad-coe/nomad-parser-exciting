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


@pytest.fixture(scope='module')
def parser():
    return ExcitingParser()


@pytest.fixture(scope='module')
def silicon_gw(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Si/INFO.OUT', archive, None)
    return archive


def test_gs(parser):
    archive = EntryArchive()
    parser.parse('tests/data/C_gs/INFO.OUT', archive, None)

    assert len(archive.section_run) == 1

    sec_run = archive.section_run[0]
    assert sec_run.program_version == 'CARBON'

    sec_method = sec_run.section_method[0]
    assert sec_method.number_of_spin_channels == 1
    assert pytest.approx(sec_method.smearing_width, 4.35974472e-22)
    assert sec_method.section_XC_functionals[1].XC_functional_name == 'GGA_X_PBE_SOL'
    assert pytest.approx(sec_method.x_exciting_scf_threshold_force_change.magnitude, 4.11936175e-12)

    sec_system = sec_run.section_system[0]
    assert pytest.approx(sec_system.lattice_vectors[0][0].magnitude, 1.72297146e-10)
    assert pytest.approx(sec_system.atom_positions[1][0].magnitude, 8.61485729e-11)
    assert len(sec_system.atom_labels) == 2
    assert sec_system.x_exciting_section_spin[0].x_exciting_spin_treatment == 'spin-unpolarised'
    assert pytest.approx(sec_system.x_exciting_section_atoms_group[0].x_exciting_muffin_tin_radius.magnitude, 6.87930374e-11)

    sec_scc = sec_run.section_single_configuration_calculation[0]
    assert pytest.approx(sec_scc.energy_total.magnitude, -3.30863556e-16)
    assert np.mean(sec_scc.atom_forces) == 0.0
    assert pytest.approx(sec_scc.charge_total.magnitude, 1.92261196e-18)
    assert pytest.approx(sec_scc.energy_reference_fermi.magnitude, 2.4422694e-18)
    assert len(sec_scc.section_scf_iteration) == 12
    assert pytest.approx(sec_scc.section_scf_iteration[5].x_exciting_valence_charge_scf_iteration.magnitude, 1.28174131e-18)
    assert pytest.approx(sec_scc.section_scf_iteration[8].x_exciting_exchange_energy_scf_iteration.magnitude, -4.39756926e-17)
    assert pytest.approx(sec_scc.section_scf_iteration[11].electronic_kinetic_energy_scf_iteration.magnitude, 3.30404896e-16)
    sec_eig = sec_scc.section_eigenvalues[0]
    assert np.shape(sec_eig.eigenvalues_kpoints) == (30, 3)
    assert pytest.approx(sec_eig.eigenvalues_values[0][9][4].magnitude, 2.74680139e-18)


def test_strucopt(parser):
    archive = EntryArchive()
    parser.parse('tests/data/GaO_strucopt/INFO.OUT', archive, None)

    sec_systems = archive.section_run[0].section_system
    assert len(sec_systems) == 15
    assert sec_systems[0].atom_labels == ['Ga', 'Ga', 'Ga', 'Ga', 'O', 'O', 'O', 'O', 'O', 'O']
    assert pytest.approx(sec_systems[0].x_exciting_gkmax.magnitude, 1.13383567e+11)
    assert pytest.approx(sec_systems[3].atom_positions[1][1].magnitude, 3.07695918e-10)
    assert pytest.approx(sec_systems[10].atom_positions[-1][0].magnitude, 3.67156876e-11)
    assert pytest.approx(sec_systems[1].lattice_vectors[2][1].magnitude, sec_systems[13].lattice_vectors[2][1].magnitude)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 15
    assert len(sec_sccs[0].section_scf_iteration) == 19
    assert pytest.approx(sec_sccs[0].section_scf_iteration[10].time_scf_iteration.magnitude, 431.84)
    assert pytest.approx(sec_sccs[0].section_scf_iteration[18].x_exciting_effective_potential_convergence_scf_iteration.magnitude, 431.84)
    assert pytest.approx(sec_sccs[3].x_exciting_maximum_force_magnitude.magnitude, 1.64771998e-10)
    assert pytest.approx(sec_sccs[6].energy_total.magnitude, -3.58415586e-14)
    assert pytest.approx(sec_sccs[9].time_calculation.magnitude, 724.33)
    assert len(sec_sccs[-1].x_exciting_section_MT_charge_atom) == 10
    assert pytest.approx(sec_sccs[-1].x_exciting_fermi_energy.magnitude, 1.03200886e-18)


def test_gw(parser):
    archive = EntryArchive()
    parser.parse('tests/data/AlP_gw/INFO.OUT', archive, None)

    sec_methods = archive.section_run[0].section_method
    assert len(sec_methods) == 2
    assert sec_methods[1].electronic_structure_method == 'G0W0'
    assert pytest.approx(sec_methods[1].gw_mixed_basis_gmax.magnitude, 3.0235618e+11)
    assert sec_methods[1].gw_self_energy_singularity_treatment == 'mpb'
    assert sec_methods[1].gw_number_of_frequencies == 32
    assert pytest.approx(sec_methods[1].gw_frequency_values[-1].magnitude, 8.22665908e-16)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 2

    # Check DFT DOS
    assert len(sec_sccs[0].section_dos) == 1
    assert sec_sccs[0].section_dos[0].dos_energies[87].magnitude == pytest.approx(-1.42127678e-18)
    assert sec_sccs[0].section_dos[0].dos_values[0][449] == pytest.approx(2.249331609726956e-10)

    # # Check GW DOS
    # assert len(sec_sccs[1].section_dos) == 1
    # assert pytest.approx(sec_sccs[1].section_dos[0].dos_energies[-40].magnitude, 1.83109278e-18,)
    # assert pytest.approx(sec_sccs[1].section_dos[0].dos_values[0][150], 2.737958376377784e-10)

    # Check GW data
    assert pytest.approx(sec_sccs[1].gw_fermi_energy.magnitude, 1.09865567e-19)
    assert pytest.approx(sec_sccs[1].gw_fundamental_gap.magnitude, 3.42913865e-19)
    assert pytest.approx(sec_sccs[1].gw_optical_gap.magnitude, 6.45981597e-19)
    assert np.shape(sec_sccs[1].section_eigenvalues[0].eigenvalues_values) == (1, 126, 20)
    assert sec_sccs[1].section_eigenvalues[0].eigenvalues_kpoints[-3][1] == 0.5
    assert pytest.approx(sec_sccs[1].section_eigenvalues[0].eigenvalues_values[0][99][9].magnitude, 1.79904866e-18)
    assert pytest.approx(sec_sccs[1].section_eigenvalues[0].gw_qp_linearization_prefactor[0][49], 0.80343)
    assert pytest.approx(sec_sccs[1].gw_self_energy_x[0][19][0].magnitude, -9.30556993e-18)
    assert pytest.approx(sec_sccs[1].gw_self_energy_c[0][109][14].magnitude, -8.85987322e-19)
    assert pytest.approx(sec_sccs[1].gw_xc_potential[0][74][6].magnitude, -2.09363661e-18)

    # Check DFT band structure
    assert len(sec_sccs[0].section_k_band) == 1
    assert len(sec_sccs[0].section_k_band[0].section_k_band_segment) == 5
    assert sec_sccs[0].section_k_band[0].section_k_band_segment[3].band_k_points[2][0] == pytest.approx(0.535714)
    assert sec_sccs[0].section_k_band[0].section_k_band_segment[1].band_energies[0][2][1].magnitude == pytest.approx(-1.13686767e-17)

    # # Check GW band structure
    # assert len(sec_sccs[1].section_k_band) == 1
    # assert len(sec_sccs[1].section_k_band[0].section_k_band_segment) == 5
    # assert sec_sccs[1].section_k_band[0].section_k_band_segment[4].band_energies[0][6][-2].magnitude == pytest.approx(1.50016281e-18)


def test_dos_spinpol(parser):
    archive = EntryArchive()
    parser.parse('tests/data/CeO_dos/INFO.OUT', archive, None)

    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert len(sec_scc.section_dos) == 1
    sec_doss = sec_scc.section_dos

    assert np.shape(sec_doss[0].dos_values) == (2, 500)
    assert sec_doss[0].dos_energies[79].magnitude == pytest.approx(-2.98206539e-18)
    assert sec_doss[0].dos_values[0][126] == pytest.approx(1.8925126966797494e-10)
    assert sec_doss[0].dos_values[1][136] == pytest.approx(1.9160612889870427e-11)
    assert sec_doss[0].dos_energies[240].magnitude == pytest.approx(-1.74389789e-19)
    assert sec_doss[0].dos_values[0][220] == pytest.approx(5.638758214688811e-10)
    assert sec_doss[0].dos_values[1][78] == pytest.approx(4.33359120779702e-10)

    sec_pdoss = sec_scc.section_atom_projected_dos
    assert len(sec_pdoss) == 1
    assert np.shape(sec_pdoss[0].atom_projected_dos_values_lm) == (25, 2, 3, 500)
    assert pytest.approx(sec_doss[0].dos_energies[100].magnitude, sec_pdoss[0].atom_projected_dos_energies[100].magnitude)
    assert pytest.approx(sec_pdoss[0].atom_projected_dos_values_lm[1][0][1][116], 4.264786582546528e-13)
    assert pytest.approx(sec_pdoss[0].atom_projected_dos_values_lm[20][1][0][85], 3.094123855740692e-17)


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


def test_band_silicon_gw(silicon_gw):
    """Tests that the band structure of silicon is parsed correctly from a 
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

        # Plot
        # n_spin_channels = energies.shape[0]
        # import matplotlib.pyplot as mpl
        # mpl.plot(energies[0], color="blue")
        # if n_spin_channels == 2:
            # mpl.plot(energies[1], color="orange")
        # mpl.show()

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
        assert gap == pytest.approx(gap_assumed)


def test_dos_silicon_gw(silicon_gw):
    """Tests that the DOS of silicon is parsed correctly.
    """
    sccs = silicon_gw.section_run[-1].section_single_configuration_calculation
    assert len(sccs) == 2
    gaps = [0.5442277, 1.360569]
    for gap_assumed, scc in zip(gaps, sccs):
        assert len(scc.section_dos) == 1
        dos = scc.section_dos[0]
        energies = dos.dos_energies.to(ureg.electron_volt).magnitude
        values = dos.dos_values

        # Check that an energy reference is reported
        energy_reference = scc.energy_reference_fermi
        if energy_reference is None:
            energy_reference = scc.energy_reference_highest_occupied
        assert energy_reference is not None
        energy_reference = energy_reference.to(ureg.electron_volt).magnitude

        # Plot
        # import matplotlib.pyplot as mpl
        # n_spin_channels = values.shape[0]
        # mpl.plot(values[0], energies, color="blue")
        # if n_spin_channels == 2:
            # mpl.plot(values[0], energies, color="orange")
        # mpl.show()

        # Check that an appropriately sized band gap is found at the given
        # reference energy
        nonzero = np.unique(values.nonzero())
        energies = energies[nonzero]
        energies.sort()
        lowest_unoccupied_index = np.searchsorted(energies, energy_reference, "right")[0]
        highest_occupied_index = lowest_unoccupied_index - 1
        gap = energies[lowest_unoccupied_index] - energies[highest_occupied_index]
        assert gap == pytest.approx(gap_assumed)

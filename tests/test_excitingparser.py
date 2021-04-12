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

from nomad.datamodel import EntryArchive
from excitingparser.exciting_parser import ExcitingParser


@pytest.fixture(scope='module')
def parser():
    return ExcitingParser()


def test_gs(parser):
    archive = EntryArchive()
    parser.parse('tests/data/C_gs/INFO.OUT', archive, None)

    assert len(archive.section_run) == 1

    sec_run = archive.section_run[0]
    assert sec_run.program_version == 'CARBON'

    sec_method = sec_run.section_method[0]
    assert sec_method.number_of_spin_channels == 1
    assert sec_method.smearing_width == pytest.approx(4.35974472e-22)
    assert sec_method.section_XC_functionals[1].XC_functional_name == 'GGA_X_PBE_SOL'
    assert sec_method.x_exciting_scf_threshold_force_change.magnitude == pytest.approx(4.11936175e-12)

    sec_system = sec_run.section_system[0]
    assert sec_system.lattice_vectors[0][0].magnitude == pytest.approx(1.72297146e-10)
    assert sec_system.atom_positions[1][0].magnitude == pytest.approx(8.61485729e-11)
    assert len(sec_system.atom_labels) == 2
    assert sec_system.x_exciting_section_spin[0].x_exciting_spin_treatment == 'spin-unpolarised'
    assert sec_system.x_exciting_section_atoms_group[0].x_exciting_muffin_tin_radius.magnitude == pytest.approx(6.87930374e-11)

    sec_scc = sec_run.section_single_configuration_calculation[0]
    assert sec_scc.energy_total.magnitude == pytest.approx(-3.30863556e-16)
    assert np.mean(sec_scc.atom_forces) == 0.0
    assert sec_scc.charge_total.magnitude == pytest.approx(1.92261196e-18)
    assert sec_scc.energy_reference_fermi.magnitude == pytest.approx(2.4422694e-18)
    assert len(sec_scc.section_scf_iteration) == 12
    assert sec_scc.section_scf_iteration[5].x_exciting_valence_charge_scf_iteration.magnitude == pytest.approx(1.28174131e-18)
    assert sec_scc.section_scf_iteration[8].x_exciting_exchange_energy_scf_iteration.magnitude == pytest.approx(-4.39756926e-17)
    assert sec_scc.section_scf_iteration[11].electronic_kinetic_energy_scf_iteration.magnitude == pytest.approx(3.30404896e-16)
    sec_eig = sec_scc.section_eigenvalues[0]
    assert np.shape(sec_eig.eigenvalues_kpoints) == (30, 3)
    assert sec_eig.eigenvalues_values[0][9][4].magnitude == pytest.approx(2.74680139e-18)


def test_strucopt(parser):
    archive = EntryArchive()
    parser.parse('tests/data/GaO_strucopt/INFO.OUT', archive, None)

    sec_systems = archive.section_run[0].section_system
    assert len(sec_systems) == 15
    assert sec_systems[0].atom_labels == ['Ga', 'Ga', 'Ga', 'Ga', 'O', 'O', 'O', 'O', 'O', 'O']
    assert sec_systems[0].x_exciting_gkmax.magnitude == pytest.approx(1.13383567e+11)
    assert sec_systems[3].atom_positions[1][1].magnitude == pytest.approx(3.07695918e-10)
    assert sec_systems[10].atom_positions[-1][0].magnitude == pytest.approx(3.67156876e-11)
    assert sec_systems[1].lattice_vectors[2][1].magnitude == pytest.approx(sec_systems[13].lattice_vectors[2][1].magnitude)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 15
    assert len(sec_sccs[0].section_scf_iteration) == 19
    assert sec_sccs[0].section_scf_iteration[10].time_scf_iteration.magnitude == pytest.approx(431.84)
    assert sec_sccs[0].section_scf_iteration[18].x_exciting_effective_potential_convergence_scf_iteration[0].magnitude == pytest.approx(4.62350928e-26)
    assert sec_sccs[3].x_exciting_maximum_force_magnitude.magnitude == pytest.approx(1.64771998e-10)
    assert sec_sccs[6].energy_total.magnitude == pytest.approx(-3.58415586e-14)
    assert sec_sccs[9].time_calculation.magnitude == pytest.approx(724.33)
    assert len(sec_sccs[-1].x_exciting_section_MT_charge_atom) == 10
    assert sec_sccs[-1].x_exciting_fermi_energy.magnitude == pytest.approx(1.03200886e-18)


def test_gw(parser):
    archive = EntryArchive()
    parser.parse('tests/data/AlP_gw/INFO.OUT', archive, None)

    sec_methods = archive.section_run[0].section_method
    assert len(sec_methods) == 2
    assert sec_methods[1].electronic_structure_method == 'G0W0'
    assert sec_methods[1].gw_mixed_basis_gmax.magnitude == pytest.approx(3.0235618e+11)
    assert sec_methods[1].gw_self_energy_singularity_treatment == 'mpb'
    assert sec_methods[1].gw_number_of_frequencies == 32
    assert sec_methods[1].gw_frequency_values[-1].magnitude == pytest.approx(8.22665908e-16)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 2
    assert len(sec_sccs[0].section_dos) == 2
    assert sec_sccs[0].section_dos[0].dos_energies[87].magnitude == pytest.approx(-1.42127678e-18)
    assert sec_sccs[0].section_dos[1].dos_energies[-40].magnitude == pytest.approx(1.83109278e-18,)
    assert sec_sccs[0].section_dos[0].dos_values[0][449] == pytest.approx(2.249331609726956e-10)
    assert sec_sccs[0].section_dos[1].dos_values[0][150] == pytest.approx(2.737958376377784e-10)
    assert sec_sccs[1].gw_fermi_energy.magnitude == pytest.approx(1.09865567e-19)
    assert sec_sccs[1].gw_fundamental_gap.magnitude == pytest.approx(3.42913865e-19)
    assert sec_sccs[1].gw_optical_gap.magnitude == pytest.approx(6.45981597e-19)
    assert np.shape(sec_sccs[1].section_eigenvalues[0].eigenvalues_values) == (1, 126, 20)
    assert sec_sccs[1].section_eigenvalues[0].eigenvalues_kpoints[-3][1] == 0.5
    assert sec_sccs[1].section_eigenvalues[0].eigenvalues_values[0][99][9].magnitude == pytest.approx(1.79904866e-18)
    assert sec_sccs[1].section_eigenvalues[0].gw_qp_linearization_prefactor[0][49][19] == pytest.approx(0.80343)
    assert sec_sccs[1].gw_self_energy_x[0][19][0].magnitude == pytest.approx(-9.30556993e-18)
    assert sec_sccs[1].gw_self_energy_c[0][109][14].magnitude == pytest.approx(-8.85987322e-19)
    assert sec_sccs[1].gw_xc_potential[0][74][6].magnitude == pytest.approx(-2.09363661e-18)
    assert len(sec_sccs[1].section_k_band) == 2
    assert len(sec_sccs[1].section_k_band[0].section_k_band_segment) == 5
    assert len(sec_sccs[1].section_k_band[1].section_k_band_segment) == 5
    assert sec_sccs[1].section_k_band[0].section_k_band_segment[3].band_k_points[2][0] == 0.535714
    assert sec_sccs[1].section_k_band[0].section_k_band_segment[1].band_energies[0][2][1].magnitude == pytest.approx(-1.13686767e-17)
    assert sec_sccs[1].section_k_band[1].section_k_band_segment[4].band_energies[0][6][-2].magnitude == pytest.approx(1.50016281e-18)


def test_dos_spinpol(parser):
    archive = EntryArchive()
    parser.parse('tests/data/CeO_dos/INFO.OUT', archive, None)

    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert len(sec_scc.section_dos) == 2
    sec_doss = sec_scc.section_dos
    assert np.shape(sec_doss[0].dos_values) == (2, 500)
    assert sec_doss[0].dos_energies[79].magnitude == pytest.approx(-2.98206539e-18)
    assert sec_doss[0].dos_values[0][126] == pytest.approx(1.8925126966797494e-10)
    assert sec_doss[0].dos_values[1][136] == pytest.approx(1.9160612889870427e-11)
    assert sec_doss[1].dos_energies[240].magnitude == pytest.approx(-1.74389789e-19,)
    assert sec_doss[1].dos_values[0][220] == pytest.approx(5.638758214688811e-10)
    assert sec_doss[1].dos_values[1][78] == pytest.approx(-4.33359120779702e-10)

    sec_pdoss = sec_scc.section_atom_projected_dos
    assert len(sec_pdoss) == 1
    assert np.shape(sec_pdoss[0].atom_projected_dos_values_lm) == (25, 2, 3, 500)
    assert sec_doss[0].dos_energies[100].magnitude == pytest.approx(sec_pdoss[0].atom_projected_dos_energies[100].magnitude)
    assert sec_pdoss[0].atom_projected_dos_values_lm[1][0][1][116] == pytest.approx(4.264786582546528e-13)
    assert sec_pdoss[0].atom_projected_dos_values_lm[20][1][0][85] == pytest.approx(3.094123855740692e-17)


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

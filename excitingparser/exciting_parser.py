import pint
import numpy as np
import os
import re
import logging
from xml.etree import ElementTree

from .metainfo import m_env
from nomad.parsing.parser import FairdiParser

from nomad.parsing.text_parser import UnstructuredTextFileParser, Quantity, FileParser
from nomad.datamodel.metainfo.public import section_single_configuration_calculation,\
    section_run, section_scf_iteration, section_system, section_method, section_XC_functionals,\
    section_sampling_method, section_dos, section_atom_projected_dos, section_k_band,\
    section_eigenvalues, section_k_band_segment, section_method_to_method_refs,\
    section_calculation_to_calculation_refs, section_frame_sequence

from .metainfo.exciting import x_exciting_section_MT_charge_atom, x_exciting_section_MT_moment_atom,\
    x_exciting_section_spin, x_exciting_section_xc, x_exciting_section_fermi_surface


class InputXMLParser(FileParser):
    def __init__(self):
        super().__init__(None)
        self._init_parameters()

    def _init_parameters(self):
        self._elements = None

    @property
    def root(self):
        if self._file_handler is None:
            if self.mainfile is None:
                return
            self._file_handler = ElementTree.parse(self.mainfile).getroot()
            self._init_parameters()

        return self._file_handler

    @property
    def elements(self):
        if self._elements is None:
            self._elements = self.root.findall('.//')

        return self._elements

    def parse(self, key):
        if self._results is None:
            self._results = dict()

        if not self.root:
            return

        key_in = key
        key = key.lstrip('/')
        if key.find('/') > 0:
            parent = os.path.dirname(key)
            child = os.path.basename(key)
            elements = self.root.findall(os.path.join('./', parent))
        else:
            elements = self.elements
            child = key

        val = []
        for element in elements:
            if child:
                v = element.attrib.get(child, None)
                if v is None:
                    v = element.findall(child)
                    v = [e.text for e in v]
                if v:
                    val.append(v)

            else:
                val.append(element.attrib)

        if not val:
            return

        def convert(val_in):
            # TODO do conversion also for dictionary elements
            for i in range(len(val_in)):
                if isinstance(val_in[i], dict):
                    for key, val in val_in[i].items():
                        val_in[i][key] = convert([val])[0]
                elif isinstance(val_in[i], str):
                    # exponential formatting
                    re_float = r'(\d+\.\d+)d(\-\d+)'
                    v = re.sub(re_float, r'\1e\2', val_in[i]).split()
                    val_in[i] = v[0] if len(v) == 1 else v
                elif isinstance(val_in[i], list):
                    val_in[i] = convert(val_in[i])

            try:
                if val_in and isinstance(val_in[0], str) and val_in[0] in ['true', 'false']:
                    val = np.array(val_in) == 'true'
                    val = list([bool(v) for v in val])
                else:
                    val = np.array(val_in, dtype=float)
                    if np.all(np.mod(val, 1) == 0):
                        val = np.array(val, dtype=int)
            except Exception:
                val = val_in

            return val

        val = convert(val)
        val = val[0] if len(val) == 1 else val

        self._results[key_in] = val


class GWInfoParser(UnstructuredTextFileParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        self._quantities = []

        def str_to_frequency(val_in):
            val = [v.split() for v in val_in.split('\n')]
            val = np.transpose(np.array([v for v in val if len(v) == 3], float))
            return dict(
                number=np.array(val[0], dtype=int), values=pint.Quantity(val[1], 'hartree'),
                weights=val[2])

        # TODO Read also input parameters here if input_GW.xml does not exist

        self._quantities.append(
            Quantity(
                'frequency_data', r'frequency list:\s*\<\s*#\s*freqs\s*weight\s*>\s*([\d\.Ee\s\-]+)',
                str_operation=str_to_frequency, repeats=False)
        )

        self._quantities.append(
            Quantity(
                'fermi_energy', r'G0W0\s*\-\s*\-+\s*[\s\S]*?Fermi energy\s*\:(\s*[\d\.]+)\s',
                unit='hartree', repeats=False)
        )

        self._quantities.append(
            Quantity(
                'direct_band_gap', r'G0W0\s*\-\s*\-+\s*[\s\S]*?Direct BandGap\s*\((?P<__unit>\w+)\)\s*\:(\s*[\d\.]+)\s',
                repeats=False)
        )

        self._quantities.append(
            Quantity(
                'fundamental_band_gap', r'G0W0\s*\-\s*\-+\s*[\s\S]*?Fundamental BandGap\s*\((?P<__unit>\w+)\)\s*\:(\s*[\d\.]+)\s',
                repeats=False)
        )

        self._quantities.append(
            Quantity(
                'optical_band_gap', r'G0W0\s*\-\s*\-+\s*[\s\S]*?Optical BandGap\s*\((?P<__unit>\w+)\)\s*\:(\s*[\d\.]+)\s',
                repeats=False)
        )


class DataTextFileParser(FileParser):
    def __init__(self, **kwargs):
        self._dtype = kwargs.get('dtype', float)
        super().__init__(None)

    def _init_parameters(self):
        pass

    @property
    def data(self):
        if self._file_handler is None:
            if self.mainfile is None:
                return

            try:
                self._file_handler = np.loadtxt(self.mainfile, dtype=self._dtype)
            except Exception:
                return

            self._init_parameters()
        return self._file_handler


class ExcitingEvalqpParser(UnstructuredTextFileParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        self._quantities = []

        def str_to_eigenvalue(val_in):
            val = val_in.strip().split('\n')
            kpts = np.array(val[0].split(), dtype=float)
            keys = val[1].split()
            eigs = np.transpose(np.array([v.split() for v in val[2:]], dtype=float))
            eigs = {keys[i]: eigs[i] for i in range(len(keys))}
            return [kpts, eigs]

        self._quantities.append(
            Quantity(
                'kpoints_eigenvalues', r'\s*k\-point \#\s*\d+:\s*([\d\s\.\-]+)([ \w\(\)]+\n)([\s\d\.\-Ee]+)',
                str_operation=str_to_eigenvalue))


class BandstructureDatParser(DataTextFileParser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._init_parameters()
        self._nspin = kwargs.get('nspin', None)
        self._energy_unit = kwargs.get('energy_unit', None)

    def _init_parameters(self):
        # TODO make a parent clas for bandstructure dat and xml
        self._nspin = None
        self._nkpts_segment = None
        self._neigs_segment = None
        self._vertices = None
        self._distances = None
        self._band_energies = None
        self._band_k_points = None

    @property
    def band_energies(self):
        if self._band_energies is None:
            if self.data is None:
                return

            data = np.transpose(self.data)
            n_kpoints = int(max(data[1]))
            bands = data[6:]
            bands = np.reshape(bands, (
                self.number_of_spin_channels, self.number_of_band_segment_eigenvalues, n_kpoints))

            self._band_energies = []
            start = 0
            for nkpts_segment in self.number_of_k_points_per_segment:
                end = start + nkpts_segment
                band_energy = np.array([np.transpose(band)[start:end] for band in bands])
                if self._energy_unit:
                    band_energy = pint.Quantity(band_energy, self._energy_unit)
                self._band_energies.append(band_energy)
                start = end

        return self._band_energies

    @property
    def band_k_points(self):
        if self._band_k_points is None:
            data = np.transpose(self.data)
            self._band_k_points = []
            start = 0
            for nkpts_segment in self.number_of_k_points_per_segment:
                end = start + nkpts_segment
                self._band_k_points.append(
                    np.transpose(data[2:5])[start:end])
                start = end

        return self._band_k_points

    @property
    def distances(self):
        if self._distances is None:
            self._distances = np.transpose(self.data)[5]

        return self._distances

    @property
    def number_of_spin_channels(self):
        if self._nspin is None:
            self._nspin = np.shape(np.transpose(self.data))[0] - 6
        return self._nspin

    @property
    def number_of_k_points_per_segment(self):
        if self._nkpts_segment is None:
            self._nkpts_segment = []
            count = 1
            for i in range(1, len(self.distances)):
                if self.distances[i] == self.distances[i - 1]:
                    self._nkpts_segment.append(count)
                    count = 1
                else:
                    count += 1
            self._nkpts_segment.append(count)

        return self._nkpts_segment

    @property
    def number_of_band_segment_eigenvalues(self):
        if self._neigs_segment is None:
            data = np.transpose(self.data)
            self._neigs_segment = int(max(data[0]))
        return self._neigs_segment


class BandOutParser(DataTextFileParser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._init_parameters()
        self._nspin = kwargs.get('nspin', None)
        self._energy_unit = kwargs.get('energy_unit', None)

    def _init_parameters(self):
        self._nspin = None
        self._distances = None
        self._band_energies = None
        self._neigs_segment = None
        self._nkpts_segment = None

    @property
    def band_energies(self):
        if self._band_energies is None:
            data = np.transpose(self.data)
            n_kpoints = np.where(data[0] == data[0][0])[0][1]
            bands = data[1:]
            bands = np.reshape(bands, (
                self.number_of_spin_channels, self.number_of_band_segment_eigenvalues, n_kpoints))

            self._band_energies = []
            start = 0
            for nkpts_segment in self.number_of_k_points_per_segment:
                end = start + nkpts_segment
                band_energy = np.array([np.transpose(band)[start:end] for band in bands])
                if self._energy_unit:
                    band_energy = pint.Quantity(band_energy, self._energy_unit)
                self._band_energies.append(band_energy)
                start = end

        return self._band_energies

    @property
    def distances(self):
        if self._distances is None:
            dist = np.transpose(self.data)[0]
            n_k_points = np.where(dist == dist[0])[0][1]
            self._distances = dist[:n_k_points]

        return self._distances

    @property
    def number_of_spin_channels(self):
        if self._nspin is None:
            self._nspin = np.shape(np.transpose(self.data)[1:])[0]
        return self._nspin

    @property
    def number_of_k_points_per_segment(self):
        if self._nkpts_segment is None:
            self._nkpts_segment = []
            count = 1
            for i in range(1, len(self.distances)):
                if self.distances[i] == self.distances[i - 1]:
                    self._nkpts_segment.append(count)
                    count = 1
                else:
                    count += 1
            self._nkpts_segment.append(count)

        return self._nkpts_segment

    @property
    def number_of_band_segment_eigenvalues(self):
        if self._neigs_segment is None:
            data = np.transpose(self.data)[0]
            self._neigs_segment = len(np.where(data == data[0])[0])
        return self._neigs_segment


class BandstructureXMLParser(FileParser):
    def __init__(self, **kwargs):
        # TODO make a parent class for dos and bandstructure
        super().__init__(None)
        self._distance_key = 'distance'
        self._coord_key = 'coord'
        self._energy_key = 'eval'
        self._vertex_key = 'vertex'
        self._band_key = 'band'
        self._init_parameters()
        self._nspin = kwargs.get('nspin', None)
        self._energy_unit = kwargs.get('energy_unit', None)

    def _init_parameters(self):
        self._nspin = None
        self._nkpts_segment = None
        self._neigs_segment = None
        self._bands = None
        self._vertices = None
        self._distances = None

    @property
    def distances(self):
        if self._distances is None:
            if not self.bands:
                return

            self._distances = [
                point.attrib.get(self._distance_key) for point in self.bands[0]]
            self._distances = np.array(self._distances, dtype=float)

        return self._distances

    @property
    def bands(self):
        if self._bands is None:
            self._bands = self.root.findall('./%s' % self._band_key)
            if not self._bands:
                # check if atom-resolved
                if self.root.findall('./species'):
                    self.logger.error(
                        'atom-resolved bandstructure currently not supported.')
        return self._bands

    @property
    def vertices(self):
        if self._vertices is None:
            self._vertices = self.root.findall('./%s' % self._vertex_key)
        return self._vertices

    @property
    def number_of_spin_channels(self):
        if self._nspin is None:
            self._nspin = 1
        return self._nspin

    @property
    def number_of_k_points_per_segment(self):
        if self._nkpts_segment is None:
            self._nkpts_segment = []
            count = 1
            for i in range(1, len(self.distances)):
                if self.distances[i] == self.distances[i - 1]:
                    self._nkpts_segment .append(count)
                    count = 1
                else:
                    count += 1
            self._nkpts_segment.append(count)

        return self._nkpts_segment

    @property
    def number_of_band_segment_eigenvalues(self):
        if self._neigs_segment is None:
            self._neigs_segment = len(self.bands) // self.number_of_spin_channels
        return self._neigs_segment

    @property
    def root(self):
        if self._file_handler is None:
            self._file_handler = ElementTree.parse(self.mainfile).getroot()
            self._init_parameters()

        return self._file_handler

    def parse(self, key):
        if self._results is None:
            self._results = dict()

        if not self.bands:
            return

        if key == 'band_energies':
            # TODO I am not certain about the format for the spin polarized case
            # I cannot find an example bandstructure file
            # How about atom-resolved bandstructures?
            res = []
            start = 0
            band_energies = np.zeros((
                self.number_of_spin_channels, self.number_of_band_segment_eigenvalues,
                len(self.distances)), dtype=float)

            for i in range(len(self.bands)):
                band_energies[i % self.number_of_spin_channels][i] = np.array(
                    [e.attrib.get(self._energy_key) for e in self.bands[i]])

            for nkpts_segment in self.number_of_k_points_per_segment:
                end = start + nkpts_segment
                band_energy = np.array([
                    np.transpose(energy)[start:end] for energy in band_energies])
                if self._energy_unit is not None:
                    band_energy = pint.Quantity(band_energy, self._energy_unit)
                res.append(band_energy)
                start = end

        elif key == 'band_k_points':
            res = []
            for i in range(len(self.number_of_k_points_per_segment)):
                start = np.array(
                    self.vertices[i].attrib.get(self._coord_key).split(), dtype=float)
                end = np.array(
                    self.vertices[i + 1].attrib.get(self._coord_key).split(), dtype=float)

                res.append(np.linspace(start, end, self.number_of_k_points_per_segment[i]))

        elif key == 'band_segm_labels':
            res = []
            for i in range(len(self.vertices) - 1):
                start = self.vertices[i].attrib.get('label')
                end = self.vertices[i + 1].attrib.get('label')
                res.append([
                    '\u0393' if start.lower() == 'gamma' else start,
                    '\u0393' if end.lower() == 'gamma' else end])

        elif key == 'band_segm_start_end':
            res = []
            for i in range(len(self.number_of_k_points_per_segment)):
                start = self.vertices[i].attrib.get(self._coord_key).split()
                end = self.vertices[i + 1].attrib.get(self._coord_key).split()
                res.append([start, end])

        else:
            res = None

        self._results[key] = res


class DOSXMLParser(FileParser):
    def __init__(self, **kwargs):
        super().__init__(None)
        self._init_parameters()
        self._nspin_key = 'nspin'
        self._totaldos_key = 'totaldos'
        self._partialdos_key = 'partialdos'
        self._diagram_key = 'diagram'
        self._l_key = 'l'
        self._m_key = 'm'
        self._energy_key = 'e'
        self._dos_key = 'dos'
        self._unit_key = 'unit'
        self._energy_unit = kwargs.get('energy_unit', None)

    def _init_parameters(self):
        self._ndos = None
        self._natoms = None
        self._nspin = None
        self._nlm = None
        self._energies = None
        self._total_dos = None
        self._partial_dos = None

    @property
    def energy_unit(self):
        if self._energy_unit is None:
            axis = self.root.find('./axis')
            if axis is None:
                return

            self._energy_unit = axis.attrib.get(self._unit_key).lower()

        return self._energy_unit

    @property
    def number_of_spin_channels(self):
        if self._nspin is None:
            if not self.total_dos:
                return
            self._nspin = len(self.total_dos)

        return self._nspin

    @property
    def number_of_atoms(self):
        if self._natoms is None:
            partial_dos = self.root.findall('./%s' % self._partialdos_key)
            self._natoms = len(partial_dos)

        return self._natoms

    @property
    def number_of_dos(self):
        if self._ndos is None:
            total_dos = self.root.find('./%s/%s' % (self._totaldos_key, self._diagram_key))
            self._ndos = len(total_dos)

        return self._ndos

    @property
    def number_of_lm(self):
        if self._nlm is None:
            if self.partial_dos is None:
                return

            self._nlm = 0
            l_list = set([int(e.attrib.get(self._l_key)) for e in self.partial_dos])
            for li in l_list:
                self._nlm += 2 * li + 1

        return self._nlm

    @property
    def total_dos(self):
        if self._total_dos is None:
            self._total_dos = self.root.findall('./%s/%s' % (self._totaldos_key, self._diagram_key))
        return self._total_dos

    @property
    def partial_dos(self):
        if self._partial_dos is None:
            self._partial_dos = self.root.findall('./%s/%s' % (self._partialdos_key, self._diagram_key))
        return self._partial_dos

    @property
    def energies(self):
        if self._energies is None:
            if self.total_dos is None:
                return

            self._energies = np.array(
                [float(point.attrib.get(self._energy_key)) for point in self.total_dos[0]])

            if self.energy_unit is not None:
                self._energies = pint.Quantity(self._energies, self.energy_unit)

        return self._energies

    @property
    def root(self):
        if self._file_handler is None:
            self._file_handler = ElementTree.parse(self.mainfile).getroot()
            self._init_parameters()

        return self._file_handler

    def _get_dos(self, diagram):
        dos = np.array(
            [point.attrib.get(self._dos_key) for point in diagram], dtype=float)

        return dos

    def parse(self, key):
        if self._results is None:
            self._results = dict()

        if 'total' in key:
            if not self.total_dos:
                return

            res = np.zeros((self.number_of_spin_channels, self.number_of_dos))

            for i in range(len(self.total_dos)):
                spin = self.total_dos[i].attrib.get(self._nspin_key, i)
                res[i] = self._get_dos(self._total_dos[i])

            if self.energy_unit is not None:
                res = pint.Quantity(res, '1/%s' % self.energy_unit)

        elif 'partial' in key:
            if not self.partial_dos:
                return

            res = np.zeros((
                self.number_of_lm, self.number_of_spin_channels, self.number_of_atoms, self.number_of_dos))

            for i in range(len(self.partial_dos)):
                spin = self.partial_dos[i].attrib.get(self._nspin_key, None)
                if spin is None:
                    spin = (i % (self.number_of_spin_channels * self.number_of_lm)) // self.number_of_lm
                else:
                    spin = int(spin) - 1

                val_l = self.partial_dos[i].attrib.get(self._l_key, None)
                val_m = self.partial_dos[i].attrib.get(self._m_key, None)
                if val_l is None or val_m is None:
                    lm = i % self.number_of_lm
                else:
                    lm = int(val_l) ** 2 + int(val_m) + int(val_l)

                atom = i // (self.number_of_lm * self.number_of_spin_channels)

                res[lm][spin][atom] = self._get_dos(self.partial_dos[i])

            if self.energy_unit is not None:
                res = pint.Quantity(res, '1/%s' % self.energy_unit)

        elif key == 'energies':
            return self.energies

        else:
            res = None

        self._results[key] = res


class ExcitingFermiSurfaceBxsfParser(UnstructuredTextFileParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        self._quantities = []

        self._quantities.append(
            Quantity(
                'fermi_energy', r'Fermi Energy:\s*([\d\.]+)\s*', unit='hartree', repeats=False))

        def str_to_band_parameters(val_in):
            val = val_in.strip().split('\n')

            nbands = int(val[0])
            mesh = np.array(val[1].split(), dtype=int)
            origin = np.array(val[2].split(), dtype=float)
            vector = np.array([v.split() for v in val[3:6]], dtype=float)

            return [nbands, mesh, origin, vector]

        self._quantities.append(
            Quantity(
                'band_parameters', r'BANDGRID_3D_BANDS\s*([\d\.\-Ee\s]+)',
                str_operation=str_to_band_parameters, repeats=False))

        self._quantities.append(
            Quantity(
                'fermi_surface', r'BAND:\s*\d+\s*([\d\-\+\.Ee\s]+)\n *E*', unit='hartree')
        )


class ExcitingEigenvalueParser(UnstructuredTextFileParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        self._quantities = []
        self._quantities.append(
            Quantity(
                'k_points', r'\s*\d+\s*([\d\.Ee\- ]+):\s*k\-point'))

        def str_to_eigenvalues(val_in):
            val = val_in[:val_in.rfind('\n \n')].strip()
            val = np.array([v.split() for v in val.split('\n')], dtype=float)
            val = np.transpose(val)
            occs = val[-1]
            eigs = val[-2]

            nspin = 2 if np.count_nonzero(occs > 1.0) == 0 else 1
            data = dict()
            data['occupancies'] = np.reshape(occs, (nspin, len(occs) // nspin))
            data['eigenvalues'] = np.reshape(eigs, (nspin, len(eigs) // nspin))
            return data

        self._quantities.append(
            Quantity(
                'eigenvalues_occupancies', r'\(state\, eigenvalue and occupancy below\)\s*([\d\.Ee\-\s]+?(?:\n *\n))',
                str_operation=str_to_eigenvalues
            )
        )


class ExcitingGWOutParser(UnstructuredTextFileParser):
    def __init__(self, mainfile, logger):
        super().__init__(mainfile, logger=logger)

    def init_quantities(self):
        self._quantities = []


class ExcitingInfoParser(UnstructuredTextFileParser):
    def __init__(self):
        super().__init__(None)

    @staticmethod
    def _re_pattern(head, key, value=r'[Ee\+\d\.\- ]+', tail='\n'):
        return r'%s[\s\S]*?%s\s*\:*\=*\s*(%s)%s' % (head, key, value, tail)

    @staticmethod
    def _str_to_quantity_tolerances(val_in):
        return val_in.strip().replace('(', '').replace(')', '').split()

    def _init_loop_quantities(self, loop_type):
        # we use the same initialization in reading properties either for loop_type
        # scf_iteraction, final_scf_iteration or final_optimization

        # TODO check for header for other exciting versions
        if loop_type == 'scf_iteration':
            header = r'(?:I|SCF i)teration number'
        elif loop_type == 'final_scf_iteration':
            header = 'Convergence targets achieved. Performing final SCF iteration'
        elif loop_type == 'final_optimization':
            header = 'Force convergence target achieved'

        self._quantities.append(
            Quantity(
                'energy_total_%s' % loop_type, self._re_pattern(header, r'[Tt]*otal energy'),
                unit='hartree')
        )

        def str_to_energy_dict(val_in):
            val = val_in.strip().split('\n')
            energies = dict()
            for v in val:
                v = v.split(':')
                if len(v) < 2:
                    continue
                energies[v[0].strip()] = pint.Quantity(float(v[1]), 'hartree')
            return energies

        self._quantities.append(
            Quantity(
                'energy_contributions_%s' % loop_type, self._re_pattern(
                    header, r'(?:Energies|_)', value=r'[\-\s\w\.\:]+?', tail=r'\n *(?:DOS|Density)'),
                str_operation=str_to_energy_dict
            )
        )

        self._quantities.append(
            Quantity(
                'x_exciting_dos_fermi_%s' % loop_type,
                self._re_pattern(
                    header, r'DOS at Fermi energy \(states\/Ha\/cell\)'),
                unit='1/hartree')
        )

        def str_to_atom_properties_dict(val_in):
            unit = None
            if 'charge' in val_in:
                unit = 'elementary_charge'
            elif 'moment' in val_in:
                unit = 'elementary_charge * bohr'

            val = val_in.strip().split('\n')

            properties = dict()
            atom_resolved = []
            species = None
            for v in val:
                v = v.strip().split(':')
                if len(v) < 2:
                    continue

                elif v[0].startswith('species'):
                    species = re.search('([A-Z][a-z]*)', v[-1]).group(1)

                elif v[0].startswith('atom'):
                    v[0] = v[0].split()
                    v[1] = [float(vi) for vi in v[1].split()]
                    v[1] = v[1][0] if len(v[1]) == 1 else v[1]
                    if species is None:
                        species = v[0][2]
                    atom_resolved.append(((species, pint.Quantity(v[1], unit))))

                else:
                    vi = [float(vii) for vii in v[1].split()]
                    vi = vi[0] if len(vi) == 1 else vi
                    properties[v[0].strip()] = pint.Quantity(vi, unit)

            properties['atom_resolved'] = atom_resolved
            return properties

        self._quantities.append(
            Quantity(
                'charge_contributions_%s' % loop_type, self._re_pattern(
                    header, r'(?:Charges|Electron charges)', value=r'[\-\s\w\.\:\(\)]+?',
                    tail=r'\n *[A-Z\+]'),
                str_operation=str_to_atom_properties_dict)
        )

        self._quantities.append(
            Quantity(
                'moment_contributions_%s' % loop_type, self._re_pattern(
                    header, r'Moments', value=r'[\-\s\w\.\:\(\)]+?',
                    tail=r'\n *[A-Z\+]'),
                str_operation=str_to_atom_properties_dict)
        )

        self._miscellaneous_keys_mapping = {
            'x_exciting_gap': (r'Estimated fundamental gap', 'hartree'),
            'time': (r'Wall time \(seconds\)', 's')}

        for name, key_unit in self._miscellaneous_keys_mapping.items():
            self._quantities.append(
                Quantity(
                    '%s_%s' % (name, loop_type), self._re_pattern(
                        header, key_unit[0]), unit=key_unit[1])
            )

        self._convergence_keys_mapping = {
            'x_exciting_effective_potential_convergence': (
                r'RMS change in effective potential \(target\)', 'hartree'),
            'x_exciting_energy_convergence': (
                r'Absolute change in total energy\s*\(target\)', 'hartree'),
            'x_exciting_charge_convergence': (
                r'Charge distance\s*\(target\)', 'elementary_charge'),
            'x_exciting_IBS_force_convergence': (
                r'Abs\. change in max\-nonIBS\-force\s*\(target\)', 'hartree/bohr')}

        for name, key_unit in self._convergence_keys_mapping.items():
            self._quantities.append(
                Quantity(
                    '%s_%s' % (name, loop_type), self._re_pattern(
                        header, key_unit[0], value=r'[\(\)\d\.\-\+Ee ]+'),
                    str_operation=self._str_to_quantity_tolerances, unit=key_unit[1])
            )

    def init_quantities(self):
        def str_to_array(val_in):
            val = [v.split(':')[-1].split() for v in val_in.strip().split('\n')]
            val = val[0] if len(val) == 1 else val
            return val

        self._quantities = [
            Quantity('program_version', r'\s*EXCITING\s*([\w\-\(\)\. ]+)\s*started', repeats=False),
            Quantity(
                'lattice_vectors', r'Lattice vectors\s*[\(cartesian\)]*\s*:\s*([\-0-9\.\s]+)\n',
                str_operation=str_to_array, unit='bohr', repeats=False),
            Quantity(
                'lattice_vectors_reciprocal', r'Reciprocal lattice vectors\s*[\(cartesian\)]*\s*:\s*([\-0-9\.\s]+)\n',
                str_operation=str_to_array, unit='1/bohr', repeats=False),
        ]

        self._system_keys_mapping = {
            'x_exciting_unit_cell_volume': ('Unit cell volume', 'bohr ** 3'),
            'x_exciting_brillouin_zone_volume': ('Brillouin zone volume', '1/bohr ** 3'),
            'x_exciting_number_of_atoms': ('Total number of atoms per unit cell', None),
            'x_exciting_spin_treatment': ('Spin treatment', None),
            'x_exciting_number_of_bravais_lattice_symmetries': ('Number of Bravais lattice symmetries', None),
            'x_exciting_number_of_crystal_symmetries': ('Number of crystal symmetries', None),
            'x_exciting_kpoint_grid': (r'k\-point grid', None),
            'x_exciting_kpoint_offset': (r'k\-point offset', None),
            'x_exciting_number_kpoints': (r'Total number of k\-points', None),
            'x_exciting_rgkmax': (r'R\^MT\_min \* \|G\+k\|\_max \(rgkmax\)', 'bohr'),
            'x_exciting_species_rtmin': (r'Species with R\^MT\_min', None),
            'x_exciting_gkmax': (r'Maximum \|G\+k\| for APW functions', '1/bohr'),
            'x_exciting_gmaxvr': (r'Maximum \|G\| for potential and density', '1/bohr'),
            'x_exciting_gvector_size': (r'G\-vector grid sizes', None),
            'x_exciting_gvector_total': (r'Total number of G\-vectors', None),
            'x_exciting_lmaxapw': (r'   APW functions', None),
            'x_exciting_nuclear_charge': ('Total nuclear charge', 'elementary_charge'),
            'x_exciting_electronic_charge': ('Total electronic charge', 'elementary_charge'),
            'x_exciting_core_charge_initial': ('Total core charge', 'elementary_charge'),
            'x_exciting_valence_charge_initial': ('Total valence charge', 'elementary_charge'),
            'x_exciting_wigner_radius': (r'Effective Wigner radius, \_s', 'bohr'),
            'x_exciting_empty_states': ('Number of empty states', None),
            'x_exciting_valence_states': ('Total number of valence states', None),
            'x_exciting_hamiltonian_size': ('Maximum Hamiltonian size', None),
            'x_exciting_pw': (r'Maximum number of plane\-waves', None),
            'x_exciting_lo': (r'Total number of local\-orbitals', None),
            'x_exciting_xc_functional': (r'Exchange\-correlation type', None)}

        self._method_keys_mapping = {
            'smearing_kind': ('Smearing scheme', None),
            'smearing_width': ('Smearing width', None)}

        for name, key_unit in self._system_keys_mapping.items():
            self._quantities.append(
                Quantity(
                    name, r'%s\s*:\s*([\s\S]*?)\n' % key_unit[0], unit=key_unit[1], repeats=False)
            )

        for name, key_unit in self._method_keys_mapping.items():
            self._quantities.append(
                Quantity(
                    name, r'%s\s*:\s*([\s\S]*?)\n' % key_unit[0], unit=key_unit[1], repeats=False)
            )

        def get_species_prop(val_in):
            val = val_in.strip().split('\n')
            val = [v.split() for v in val]
            prop = dict()
            prop['number'] = int(val[0][0])
            prop['symbol'] = val[0][1]
            prop['file'] = val[0][2]
            prop['name'] = val[0][3]
            prop['nuclear_charge'] = pint.Quantity(float(val[0][4]), 'elementary_charge')
            prop['electronic_charge'] = pint.Quantity(float(val[0][5]), 'elementary_charge')
            prop['atomic_mass'] = pint.Quantity(float(val[0][6]), 'electron_mass')
            prop['muffin_tin_radius'] = pint.Quantity(float(val[0][7]), 'bohr')
            prop['radial_points'] = int(val[0][8])
            prop['positions_format'] = val[0][9].lstrip('(').rstrip(')')

            positions = np.zeros((len(val) - 1, 3), dtype=float)
            for i in range(1, len(val)):
                positions[i - 1] = val[i][2:5]

            prop['positions'] = positions

            return prop

        species_prop = [
            'parameters loaded from', 'name', 'nuclear charge', 'electronic charge',
            'atomic mass', r'muffin\-tin radius', r'[number#]* of radial points in muffin\-tin']

        species_pattern = r'Species\s*:(\s*\d+\s*)\((\w+)\)\s*' + ''.join(
            [r'%s\s*:(\s*[\s\S]*?)\n *' % p for p in species_prop]
        ) + r'\s*atomic positions( \(\w+\))\s*[\, \w\(\)]*:(\s*[0-9\-\:\.\s]+)'

        self._quantities.append(
            Quantity('species', species_pattern, str_operation=get_species_prop)
        )

        self.quantities.append(
            Quantity('potential_mixing', r'Using ([\w ]+) potential mixing')
        )

        # different names for different versions of exciting
        self._energy_keys_mapping = {
            'energy_total': ['Total energy', 'total energy'],
            'x_exciting_fermi_energy': ['Fermi energy', 'Fermi'],
            'electronic_kinetic_energy': ['Kinetic energy', 'electronic kinetic'],
            'x_exciting_coulomb_energy': ['Coulomb energy', 'Coulomb'],
            'x_exciting_exchange_energy': ['Exchange energy', 'exchange'],
            'x_exciting_correlation_energy': ['Correlation energy', 'correlation'],
            'energy_sum_eigenvalues': ['Sum of eigenvalues', 'sum of eigenvalues'],
            'x_exciting_effective_potential_energy': ['Effective potential energy'],
            'x_exciting_coulomb_potential_energy': ['Coulomb potential energy', 'Coulomb potential'],
            'energy_XC_potential': ['xc potential energy', 'xc potential'],
            'x_exciting_hartree_energy': ['Hartree energy', 'Hartree'],
            'x_exciting_electron_nuclear_energy': ['Electron-nuclear energy', 'electron-nuclear '],
            'x_exciting_nuclear_nuclear_energy': ['Nuclear-nuclear energy', 'nuclear-nuclear'],
            'x_exciting_madelung_energy': ['Madelung energy', 'Madelung'],
            'x_exciting_core_electron_kinetic_energy': ['Core-electron kinetic energy', 'core electron kinetic'],
            'x_exciting_dft_d2_dispersion_correction': ['DFT-D2 dispersion correction']
        }

        self._electron_charge_keys_mapping = {
            'x_exciting_core_charge': ['core'],
            'x_exciting_core_leakage': ['core leakage'],
            'x_exciting_valence_charge': ['valence'],
            'x_exciting_interstitial_charge': ['interstitial'],
            'x_exciting_total_MT_charge': ['total charge in muffin-tins', 'total in muffin-tins'],
            'charge_total': ['total charge'],
            'x_exciting_section_MT_charge_atom': ['atom_resolved']
        }

        self._moment_keys_mapping = {
            'x_exciting_interstitial_moment': ['interstitial'],
            'x_exciting_total_MT_moment': ['total moment in muffin-tins'],
            'x_exciting_total_moment': ['total moment'],
            'x_exciting_section_MT_moment_atom': ['atom_resolved']
        }

        self._init_loop_quantities('scf_iteration')
        self._init_loop_quantities('final_scf_iteration')

        def str_to_symbols(val_in):
            return [v.split()[2] for v in val_in.strip().split('\n')]

        self._quantities.append(
            Quantity(
                'positions_format_final_scf_iteration', re_pattern=self._re_pattern(
                    r'Self-consistent loop stopped', r'Atomic positions\s*\(',
                    value=r'[a-z]+', tail=r'\)'))
        )

        self._quantities.append(
            Quantity(
                'atom_symbols_final_scf_iteration', re_pattern=self._re_pattern(
                    r'Self-consistent loop stopped', r'Atomic positions \(\w+\)',
                    value=r'\s*atom[\-\s\w\.\:]*?', tail=r'\n *Total'),
                str_operation=str_to_symbols)
        )

        self._quantities.append(
            Quantity(
                'atom_positions_final_scf_iteration', re_pattern=self._re_pattern(
                    r'Self-consistent loop stopped', r'Atomic positions \(\w+\)',
                    value=r'\s*atom[\-\s\w\.\:]*?', tail=r'\n *Total'),
                str_operation=str_to_array)
        )

        self._quantities.append(
            Quantity(
                'atom_forces_final_scf_iteration', re_pattern=self._re_pattern(
                    r'Self-consistent loop stopped', r'Total atomic forces including IBS \(\w+\)',
                    value=r'\s*atom[\-\s\w\.\:]*?', tail=r'\n *Atomic'),
                str_operation=str_to_array, unit='hartree/bohr')
        )

        self._quantities.append(
            Quantity(
                'atom_positions_optimization', re_pattern=self._re_pattern(
                    r'Optimization step', r'Atomic positions at this step \(\w+\)',
                    value=r'\s*atom[\-\s\w\.\:]*?', tail=r'\n *Total'),
                str_operation=str_to_array)
        )

        self._quantities.append(
            Quantity(
                'atom_forces_optimization', re_pattern=self._re_pattern(
                    r'Optimization step', r'Total atomic forces including IBS \(\w+\)',
                    value=r'\s*atom[\-\s\w\.\:]*?', tail=r'\n *Time'),
                str_operation=str_to_array, unit='hartree/bohr')
        )

        self._quantities.append(
            Quantity(
                'atom_symbols_optimization', re_pattern=self._re_pattern(
                    r'Optimization step', r'Atomic positions at this step \(\w+\)',
                    value=r'\s*atom[\-\s\w\.\:]*?', tail=r'\n *Total'),
                str_operation=str_to_symbols)
        )

        self._quantities.append(
            Quantity(
                'nstep_optimization', re_pattern=self._re_pattern(
                    r'Optimization step', r'', value=r'\d+', tail=r'\s*\(method'))
        )

        self._quantities.append(
            Quantity(
                'method_optimization', re_pattern=self._re_pattern(
                    r'Optimization step', r'method', value=r'\w+', tail=r'\)'))
        )

        self._quantities.append(
            Quantity(
                'n_scf_iterations_optimization', re_pattern=self._re_pattern(
                    r'Optimization step', r'Number of total scf iterations', value=r'\d+',))
        )

        self._quantities.append(
            Quantity(
                'force_convergence_optimization', re_pattern=self._re_pattern(
                    r'Optimization step', r'Maximum force magnitude\s*\(target\)',
                    value=r'[\(\)\d\.\-\+Ee ]+'), str_operation=self._str_to_quantity_tolerances,
                unit='hartree/bohr')
        )

        self._quantities.append(
            Quantity(
                'total_energy_optimization', re_pattern=self._re_pattern(
                    r'Optimization step', r'Total energy at this optimization step'),
                unit='hartree')
        )

        self._quantities.append(
            Quantity(
                'time_optimization', re_pattern=self._re_pattern(
                    r'Optimization step', r'Time spent in this optimization step', tail='seconds'),
                unit='s')
        )

        self._init_loop_quantities('final_optimization')

        self._quantities.append(
            Quantity(
                'atom_positions_final_optimization', re_pattern=self._re_pattern(
                    r'Force convergence target achieved', r'Optimized atomic positions \(\w+\)',
                    value=r'\s*atom[\-\s\w\.\:]*?', tail=r'\n *Total'),
                str_operation=str_to_array)
        )

        self._quantities.append(
            Quantity(
                'atom_forces_final_optimization', re_pattern=self._re_pattern(
                    r'Force convergence target achieved', r'Total atomic forces including IBS \(\w+\)',
                    value=r'\s*atom[\-\s\w\.\:]*?', tail=r'\n *Atomic'),
                str_operation=str_to_array, unit='hartree/bohr')
        )

    def get_atom_labels(self):
        # atom labels do not change for all configurations?
        # first, we get it from the initial ground state configuration
        labels = self.get('atom_symbols_final_scf_iteration', [None])[0]

        if labels is None:
            # we get it by concatenating species symbols
            species = self.get('species', [])
            labels = []
            for specie in species:
                labels += [specie.get('symbol')] * len(specie.get('positions'))
        return labels

    def get_positions_format(self):
        positions_format = self.get('positions_format_final_scf_iteration', [None])[0]

        if positions_format is None:
            species = self.get('species', [])
            for specie in species:
                positions_format = specie.get('positions_format', None)
                if positions_format is not None:
                    break

        return positions_format

    def get_atom_positions(self, loop_type, index):
        positions = self.get('atom_positions_%s' % loop_type, [None] * (index + 1))[index]

        if positions is None and loop_type == 'final_scf_iteration':
            # we get it from initial positions in species
            species = self.get('species', [])
            if species:
                positions = np.vstack([s.get('positions') for s in species])

        if positions is None:
            return

        positions = np.array(positions)
        positions_format = self.get_positions_format()

        if positions_format == 'lattice':
            cell = self.get('lattice_vectors', None)
            if cell is None:
                return
            positions = np.dot(positions, cell.magnitude)

        return pint.Quantity(positions, 'bohr')

    def get_scf_threshold(self, name):
        thresholds = self.get('%s_scf_iteration' % name, None)
        if thresholds is None:
            return

        for i in range(len(thresholds)):
            if thresholds[i] is not None:
                return thresholds[i][1]

    def get_scf_quantity(self, name):
        n_scf = len(self.get('energy_total_scf_iteration', []))
        quantity = self.get('%s_scf_iteration' % name)
        if quantity is None:
            return

        # this is really problematic if some scf steps dont have the quantity
        # the only thing that we can do is to assume that the first steps are the
        # ones with the missing quantity
        if len(quantity) < n_scf:
            quantity = [None] * (n_scf - len(quantity)) + quantity

        return quantity

    def get_xc_functional_name(self):
        xc_functional_map = {
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
            406: ['HYB_GGA_XC_PBEH']}

        xc_functional = self.get('x_exciting_xc_functional', None)
        if xc_functional is None:
            return []

        name = xc_functional_map.get(xc_functional, [])

        return name

    @property
    def n_optimization_steps(self):
        return len(self.get('nstep_optimization', []))

    def get_number_of_spin_channels(self):
        spin_treatment = self.get('x_exciting_spin_treatment', ['spin-unpolarised'])[0]
        n_spin = 1 if spin_treatment.lower() == 'spin-unpolarised' else 2
        return n_spin

    def get_unit_cell_volume(self):
        return self.get('x_exciting_unit_cell_volume')


class ExcitingParser(FairdiParser):
    def __init__(self):
        super().__init__(
            name='parsers/exciting', code_name='exciting', code_homepage='http://exciting-code.org/',
            mainfile_name_re=r'^.*.OUT(\.[^/]*)?$', mainfile_contents_re=(r'EXCITING.*started'))
        self._metainfo_env = m_env

        self.info_parser = ExcitingInfoParser()
        self.dos_parser = DOSXMLParser(energy_unit='hartree')
        self.bandstructure_parser = BandstructureXMLParser(energy_unit='hartree')
        self.eigval_parser = ExcitingEigenvalueParser()
        self.fermisurf_parser = ExcitingFermiSurfaceBxsfParser()
        self.evalqp_parser = ExcitingEvalqpParser()
        self.dos_out_parser = DataTextFileParser()
        self.bandstructure_dat_parser = BandstructureDatParser(energy_unit='hartree')
        self.band_out_parser = BandOutParser(energy_unit='hartree')
        self.info_gw_parser = GWInfoParser()
        self.input_xml_parser = InputXMLParser()
        self.data_xs_parser = DataTextFileParser()
        self.data_clathrate_parser = DataTextFileParser(dtype=str)

    def get_exciting_files(self, default):
        filename = os.path.join(self.info_parser.maindir, default)
        if not os.path.isfile(filename):
            file_ext = default.split('.')[-1]
            mainfile_base = os.path.basename(
                self.info_parser.mainfile).split('.')[0].replace('INFO', '')
            file_base = default.split('.')[0]

            options = [
                f for f in os.listdir(
                    self.info_parser.maindir) if file_base in f and mainfile_base in f]
            options = [f for f in options if f.endswith(file_ext)]
            options.sort()

            filenames = [os.path.join(self.info_parser.maindir, f) for f in options]
        else:
            filenames = [filename]

        filenames = [f for f in filenames if os.access(f, os.F_OK)]
        return filenames

    def _parse_dos(self, sec_scc):
        if self.dos_parser.get('totaldos', None) is None:
            return

        sec_dos = sec_scc.m_create(section_dos)
        sec_dos.dos_kind = 'electronic'
        sec_dos.number_of_dos_values = self.dos_parser.number_of_dos
        sec_dos.dos_energies = self.dos_parser.energies
        totaldos = self.dos_parser.get('totaldos') * self.info_parser.get_unit_cell_volume()
        # metainfo does not unit a unit for dos
        sec_dos.dos_values = totaldos.to('m**3/joule').magnitude

        partialdos = self.dos_parser.get('partialdos')
        if partialdos is None:
            return

        partialdos = partialdos * self.info_parser.get_unit_cell_volume()
        partialdos = partialdos.to('m**3/J').magnitude
        sec_atom_projected_dos = sec_scc.m_create(section_atom_projected_dos)
        sec_atom_projected_dos.atom_projected_dos_m_kind = 'spherical'
        sec_atom_projected_dos.number_of_lm_atom_projected_dos = self.dos_parser.number_of_lm
        sec_atom_projected_dos.atom_projected_dos_energies = self.dos_parser.energies
        sec_atom_projected_dos.atom_projected_dos_values_lm = partialdos

    def _parse_bandstructure(self, sec_scc):
        # we need to set nspin again as this is overwritten when setting mainfile
        self.bandstructure_parser._nspin = self.info_parser.get_number_of_spin_channels()

        band_energies = self.bandstructure_parser.get('band_energies', None)
        if band_energies is None:
            return

        sec_k_band = sec_scc.m_create(section_k_band)
        sec_k_band.band_structure_kind = 'electronic'

        # imho the number of eigenvalues should also be property of the section_k_band
        # since because the band structure data is not necessarily run with the
        # same settings as that of the eigevalues
        if not sec_scc.section_eigenvalues:
            sec_eigenvalues = sec_scc.m_create(section_eigenvalues)
        else:
            sec_eigenvalues = sec_scc.section_eigenvalues[-1]
        sec_eigenvalues.number_of_band_segment_eigenvalues =\
            self.bandstructure_parser.number_of_band_segment_eigenvalues

        band_k_points = self.bandstructure_parser.get('band_k_points')
        nkpts_segment = self.bandstructure_parser.number_of_k_points_per_segment
        band_seg_labels = self.bandstructure_parser.get('band_segm_labels')
        band_seg_start_end = self.bandstructure_parser.get('band_segm_start_end')
        for nb in range(len(band_energies)):
            sec_k_band_segment = sec_k_band.m_create(section_k_band_segment)
            sec_k_band_segment.number_of_k_points_per_segment = nkpts_segment[nb]
            sec_k_band_segment.band_k_points = band_k_points[nb]
            sec_k_band_segment.band_energies = band_energies[nb]
            sec_k_band_segment.band_segm_labels = band_seg_labels[nb]
            sec_k_band_segment.band_segm_start_end = band_seg_start_end[nb]

    def _parse_eigenvalues(self, sec_scc):
        if self.eigval_parser.get('eigenvalues_occupancies', None) is None:
            return

        sec_eigenvalues = sec_scc.m_create(section_eigenvalues)

        def get_data(key):
            data = self.eigval_parser.get('eigenvalues_occupancies')
            res = np.hstack([v[key] for v in data])
            res = res.reshape((len(res), len(data), len(res[0]) // len(data)))

            if key == 'eigenvalues':
                res = pint.Quantity(res, 'hartree')
            return res

        sec_eigenvalues.eigenvalues_values = get_data('eigenvalues')
        sec_eigenvalues.eigenvalues_occupation = get_data('occupancies')
        sec_eigenvalues.eigenvalues_kpoints = self.eigval_parser.get('k_points')

    def _parse_fermisurface(self, sec_scc):
        fermi_surface = self.fermisurf_parser.get('fermi_surface', [None])[0]
        if fermi_surface is None:
            return

        sec_fermisurface = sec_scc.m_create(x_exciting_section_fermi_surface)

        band_parameters = self.fermisurf_parser.get('band_parameters', None)
        if band_parameters is not None:
            sec_fermisurface.x_exciting_number_of_bands_fermi_surface = band_parameters[0]
            sec_fermisurface.x_exciting_number_of_mesh_points_fermi_surface = np.product(band_parameters[1])
            sec_fermisurface.x_exciting_grid_fermi_surface = band_parameters[1]
            sec_fermisurface.x_exciting_origin_fermi_surface = band_parameters[2]
            sec_fermisurface.x_exciting_vectors_fermi_surface = band_parameters[3]

        fermi_energy = self.fermisurf_parser.get('fermi_energy', None)
        if fermi_energy is not None:
            sec_fermisurface.x_exciting_fermi_energy_fermi_surface = fermi_energy

        sec_fermisurface.x_exciting_values_fermi_surface = fermi_surface

    def _parse_evalqp(self, sec_scc):
        data = self.evalqp_parser.get('kpoints_eigenvalues')
        if data is None:
            return

        sec_eigenvalues = sec_scc.m_create(section_eigenvalues)

        def get_data(key):
            if key == 'k_points':
                return np.array([d[0] for d in data])
            elif key == 'Znk':
                return np.array([d[1].get(key, None) for d in data])
            else:
                energy = np.array([d[1].get(key, None) for d in data])
                if None in energy:
                    return
                return pint.Quantity(np.array([d[1].get(key) for d in data]), 'hartree')

        eigs_gw = get_data('E_GW')
        if eigs_gw[0] is None:
            return

        sec_eigenvalues.number_of_eigenvalues = len(eigs_gw[0])
        sec_eigenvalues.number_of_eigenvalues_kpoints = len(eigs_gw)
        sec_eigenvalues.eigenvalues_kpoints = get_data('k_points')
        sec_eigenvalues.eigenvalues_values = eigs_gw
        sec_eigenvalues.gw_qp_linearization_prefactor = get_data('Znk')

        sec_scc.gw_self_energy_x = get_data('Sx')
        self_energy = get_data('Sc')
        if self_energy is None:
            self_energy = get_data('Re(Sc)')
        sec_scc.gw_self_energy_c = self_energy
        sec_scc.gw_xc_potential = get_data('Vxc')

    def _parse_dos_out(self, sec_scc):
        data = self.dos_out_parser.data
        if data is None:
            return

        # TODO I am not sure about format for spin-polarized case! I assume it is
        # energy dos_up dos_down
        nspin = self.info_parser.get_number_of_spin_channels()
        if nspin != len(data) - 1:
            self.logger.error('Found inconsistent number of spin channels in gw dos!')
            return

        sec_dos = sec_scc.m_create(section_dos)
        sec_dos.number_of_dos_values = len(data)

        data = np.transpose(data)
        sec_dos.dos_energies = pint.Quantity(data[0], 'hartree')
        # metainfo does not have unit for dos
        dos = pint.Quantity(data[1:], '1/hartree') * self.info_parser.get_unit_cell_volume()
        dos = dos.to('m**3/joule').magnitude
        sec_dos.dos_values = dos

    def _parse_bandstructure_dat(self, sec_scc):
        self.bandstructure_dat_parser._nspin = self.info_parser.get_number_of_spin_channels()

        band_energies = self.bandstructure_dat_parser.band_energies
        if band_energies is None:
            return

        sec_k_band = sec_scc.m_create(section_k_band)
        sec_k_band.band_structure_kind = 'electronic'

        band_k_points = self.bandstructure_dat_parser.band_k_points
        nkpts_segment = self.bandstructure_dat_parser.number_of_k_points_per_segment
        for nb in range(len(band_energies)):
            sec_k_band_segment = sec_k_band.m_create(section_k_band_segment)
            sec_k_band_segment.number_of_k_points_per_segment = nkpts_segment[nb]
            sec_k_band_segment.band_k_points = band_k_points[nb]
            sec_k_band_segment.band_energies = band_energies[nb]

    def _parse_band_out(self, sec_scc):
        self.band_out_parser._nspin = self.info_parser.get_number_of_spin_channels()

        band_energies = self.band_out_parser.band_energies
        if band_energies is None:
            return

        sec_k_band = sec_scc.m_create(section_k_band)
        sec_k_band.band_structure_kind = 'electronic'

        sec_k_band = sec_scc.m_create(section_k_band)
        sec_k_band.band_structure_kind = 'electronic'

        nkpts_segment = self.band_out_parser.number_of_k_points_per_segment
        for nb in range(len(band_energies)):
            sec_k_band_segment = sec_k_band.m_create(section_k_band_segment)
            sec_k_band_segment.number_of_k_points_per_segment = nkpts_segment[nb]
            sec_k_band_segment.band_energies = band_energies[nb]

    def parse_file(self, name, section):
        # TODO add support for info.xml
        if name.startswith('dos') and name.endswith('xml'):
            parser = self.dos_parser
            parser_function = self._parse_dos
        elif name.startswith('bandstructure') and name.endswith('xml'):
            parser = self.bandstructure_parser
            parser_function = self._parse_bandstructure
        elif name.startswith('EIGVAL') and name.endswith('OUT'):
            parser = self.eigval_parser
            parser_function = self._parse_eigenvalues
        elif (name.startswith('FERMISURF') or name.startswith('FS')) and name.endswith('bxsf'):
            parser = self.fermisurf_parser
            parser_function = self._parse_fermisurface
        elif name.startswith('EVALQP') and (name.endswith('DAT') or name.endswith('TXT')):
            parser = self.evalqp_parser
            parser_function = self._parse_evalqp
        elif name.startswith('TDOS') and name.endswith('OUT'):
            parser = self.dos_out_parser
            parser_function = self._parse_dos_out
        elif name.startswith('bandstructure') and name.endswith('dat'):
            parser = self.bandstructure_dat_parser
            parser_function = self._parse_bandstructure_dat
        elif name.startswith('BAND') and name.endswith('OUT'):
            parser = self.band_out_parser
            parser_function = self._parse_band_out
        elif name.startswith('input') and name.endswith('xml'):
            parser = self.input_xml_parser
            if self._calculation_type == 'gw':
                parser_function = self._parse_input_gw
            elif self._calculation_type == 'xs':
                parser_function = self._parse_input_xs
            else:
                # TODO implement reading of parameters from input.xml for normal calculations
                # in addition to INFO.OUT
                return
        else:
            return

        files = self.get_exciting_files(name)
        if len(files) > 1:
            self.logger.warn('Found multiple files of type %s. Will read all!' % name)

        for n in range(len(files)):
            parser.mainfile = files[n]
            parser_function(section)

        # free up memory
        parser.mainfile = None

    def _parse_input_xs(self, sec_method):
        xstype = self.input_xml_parser.get('xs/xstype', None)
        if xstype is not None:
            sec_method.x_exciting_xs_xstype = xstype
            sec_method.x_exciting_electronic_structure_method = xstype

        sec_method.x_exciting_xs_broadening = self.input_xml_parser.get(
            'xs/broad', 0.01, 'hartree')
        sec_method.x_exciting_xs_gqmax = self.input_xml_parser.get(
            'xs/gqmax', 0.0, '1/bohr')
        sec_method.x_exciting_xs_lmaxapw = self.input_xml_parser.get('xs/lmaxapw', 10)
        sec_method.x_exciting_xs_number_of_empty_states = self.input_xml_parser.get(
            'xs/nempty', 5)
        sec_method.x_exciting_xs_ngridq = self.input_xml_parser.get('xs/ngridq', [1, 1, 1])
        sec_method.x_exciting_xs_ngridk = self.input_xml_parser.get('xs/ngridk', [1, 1, 1])
        rgkmax = self.input_xml_parser.get('xs/rgkmax', None)
        if rgkmax is None:
            rgkmax = self.info_parser.get('x_exciting_rgkmax', 0.)
        sec_method.x_exciting_xs_rgkmax = rgkmax
        sec_method.x_exciting_xs_scissor = self.input_xml_parser.get('xs/scissor', 0.0)
        sec_method.x_exciting_xs_vkloff = self.input_xml_parser.get('xs/vkloff', [0., 0., 0.])

        # TODO I am not certain if screening/BSE are children of xs
        if self.input_xml_parser.get('xs/screening') is not None:
            sec_method.x_exciting_xs_screening_number_of_empty_states = self.input_xml_parser.get(
                'xs/screening/nempty', 0)
            sec_method.x_exciting_xs_screening_ngridk = self.input_xml_parser.get(
                'xs/screening/ngridk', [0, 0, 0])
            rgkmax = self.input_xml_parser.get('xs/screening/rgkmax', None)
            if rgkmax is None:
                rgkmax = self.info_parser.get('x_exciting_rgkmax', 0.)
            sec_method.x_exciting_xs_screening_rgkmax = rgkmax
            sec_method.x_exciting_xs_screening_type = self.input_xml_parser.get(
                'xs/screening/screentype', 'full')

        if self.input_xml_parser.get('xs/BSE') is not None:
            sec_method.x_exciting_xs_bse_antiresonant = self.input_xml_parser.get(
                'xs/BSE/aresbse', True)
            sec_method.x_exciting_xs_bse_angular_momentum_cutoff = self.input_xml_parser.get(
                'xs/BSE/lmaxdielt', 14)
            rgkmax = self.input_xml_parser.get('xs/BSE/rgkmax', None)
            if rgkmax is None:
                rgkmax = self.info_parser.get('x_exciting_rgkmax', 0)

            sec_method.x_exciting_xs_bse_rgkmax = rgkmax
            sec_method.x_exciting_xs_bse_sciavbd = self.input_xml_parser.get(
                'xs/BSE/sciavbd', True)
            sec_method.x_exciting_xs_bse_sciavqbd = self.input_xml_parser.get(
                'xs/BSE/sciavqbd', False)
            sec_method.x_exciting_xs_bse_sciavqhd = self.input_xml_parser.get(
                'xs/BSE/sciavqhd', False)
            sec_method.x_exciting_xs_bse_sciavqwg = self.input_xml_parser.get(
                'xs/BSE/sciavqwg', False)
            sec_method.x_exciting_xs_bse_sciavtype = self.input_xml_parser.get(
                'xs/BSE/sciavtype', 'spherical')
            sec_method.x_exciting_xs_bse_xas = self.input_xml_parser.get(
                'xs/BSE/xas', False)
            sec_method.x_exciting_xs_bse_number_of_bands = self.input_xml_parser.get(
                'xs/BSE/nstlbse', [0, 0, 0, 0])
            if sec_method.x_exciting_xs_bse_xas:
                sec_method.x_exciting_xs_bse_xasatom = self.input_xml_parser.get(
                    'xs/BSE/xasatom', 0)
                sec_method.x_exciting_xs_bse_xasedge = self.input_xml_parser.get(
                    'xs/BSE/xasedge', 'K')
                sec_method.x_exciting_xs_bse_xasspecies = self.input_xml_parser.get(
                    'xs/BSE/xasspecies', 0)
                sec_method.x_exciting_xs_bse_xas_number_of_bands = self.input_xml_parser.get(
                    'xs/BSE/nstlxas', [0, 0])

        if self.input_xml_parser.get('xs/tddft') is not None:
            sec_method.x_exciting_xs_tddft_analytic_continuation = self.input_xml_parser.get(
                'xs/tddft/acont', False)
            sec_method.x_exciting_xs_tddft_anomalous_Hall_conductivity = self.input_xml_parser.get(
                'xs/tddft/ahc', False)
            sec_method.x_exciting_xs_tddft_anti_resonant_dielectric = self.input_xml_parser.get(
                'xs/tddft/aresdf', False)
            sec_method.x_exciting_xs_tddft_anti_resonant_xc_kernel = self.input_xml_parser.get(
                'xs/tddft/aresfxc', True)
            sec_method.x_exciting_xs_tddft_drude = self.input_xml_parser.get(
                'xs/tddft/drude', [0., 0.])
            sec_method.x_exciting_xs_tddft_split_parameter = self.input_xml_parser.get(
                'xs/tddft/fxcbsesplit', 0.00001, 'hartree')
            sec_method.x_exciting_xs_tddft_xc_kernel = self.input_xml_parser.get(
                'xs/tddft/fxctype', 'RPA')
            sec_method.x_exciting_xs_tddft_finite_q_intraband_contribution = self.input_xml_parser.get(
                'xs/tddft/intraband', False)
            sec_method.x_exciting_xs_tddft_diagonal_xc_kernel = self.input_xml_parser.get(
                'xs/tddft/kerndiag', False)
            sec_method.x_exciting_xs_tddft_lmax_alda = self.input_xml_parser.get(
                'xs/tddft/lmaxalda', 3)
            sec_method.x_exciting_xs_tddft_macroscopic_dielectric_function_q_treatment = self.input_xml_parser.get(
                'xs/tddft/mdfqtype', 0)
            sec_method.x_exciting_xs_tddft_analytic_continuation_number_of_intervals = self.input_xml_parser.get(
                'xs/tddft/nwacont', 0)
            sec_method.x_exciting_xs_tetra = self.input_xml_parser.get(
                'xs/tetra/tetradf', False)

    def _parse_xs_bse(self):
        sec_run = self.archive.section_run[-1]

        def get_files(name):
            bse_types = ['IP', 'singlet', 'triplet', 'RPA']
            scr_types = ['full', 'diag', 'noinvdiag', 'longrange']
            bse_files = []
            for bse_type in bse_types:
                for scr_type in scr_types:
                    files = self.get_exciting_files(
                        '%s_BSE%s_SCR%s.OUT' % (name, bse_type, scr_type))
                    bse_files.append(files)

            if len([f for f in bse_files if f]) > 1:
                self.logger.warn('Multiple %s BSE types identified.' % name)

            return bse_files

        def get_data(files):
            data = []
            for f in files:
                self.data_xs_parser.mainfile = f
                if self.data_xs_parser.data is None:
                    continue
                data.append(self.data_xs_parser.data)

            return data

        def parse_exciton(data, sec_scc):
            n_components = len(data)
            data = np.transpose(np.vstack(data))

            sec_scc.x_exciting_xs_bse_number_of_components = n_components
            n_excitons = len(data[0]) // n_components
            sec_scc.x_exciting_xs_bse_number_of_excitons = n_excitons
            sec_scc.x_exciting_xs_bse_exciton_energies = pint.Quantity(
                np.reshape(data[1], (n_components, n_excitons)), 'hartree')
            sec_scc.x_exciting_xs_bse_exciton_binding_energies = pint.Quantity(
                np.reshape(data[2], (n_components, n_excitons)), 'hartree')
            sec_scc.x_exciting_xs_bse_exciton_oscillator_strength = np.reshape(
                data[3], (n_components, n_excitons))
            sec_scc.x_exciting_xs_bse_exciton_amplitude_re = np.reshape(
                data[4], (n_components, n_excitons))
            sec_scc.x_exciting_xs_bse_exciton_amplitude_im = np.reshape(
                data[5], (n_components, n_excitons))

        def parse_epsilon(data, sec_scc):
            n_components = len(data)
            data = np.transpose(np.vstack(data))
            n_epsilon = len(data[0]) // n_components

            sec_scc.x_exciting_xs_bse_number_of_energy_points = n_epsilon
            sec_scc.x_exciting_xs_bse_epsilon_energies = pint.Quantity(
                np.reshape(data[0], (n_components, n_epsilon)), 'hartree')
            sec_scc.x_exciting_xs_bse_epsilon_re = np.reshape(
                data[1], (n_components, n_epsilon))
            sec_scc.x_exciting_xs_bse_epsilon_im = np.reshape(
                data[2], (n_components, n_epsilon))

        def parse_sigma(data, sec_scc):
            n_components = len(data)
            data = np.transpose(np.vstack(data))
            n_sigma = len(data[0]) // n_components

            sec_scc.x_exciting_xs_bse_sigma_energies = pint.Quantity(
                np.reshape(data[0], (n_components, n_sigma)), 'hartree')
            sec_scc.x_exciting_xs_bse_sigma_re = np.reshape(
                data[1], (n_components, n_sigma))
            sec_scc.x_exciting_xs_bse_sigma_im = np.reshape(
                data[2], (n_components, n_sigma))

        def parse_loss(data, sec_scc):
            n_components = len(data)
            data = np.transpose(np.vstack(data))
            n_loss = len(data[0]) // n_components

            sec_scc.x_exciting_xs_bse_loss_energies = pint.Quantity(
                np.reshape(data[0], (n_components, n_loss)), 'hartree')
            sec_scc.x_exciting_xs_bse_loss = np.reshape(
                data[1], (n_components, n_loss))

        # TODO check if format of files are really correct, i.e. columns are supposed
        # to be what they are. What is the fourth column in epsilon which is not parsed?
        sccs = []
        for quantity in ['EXCITON', 'EPSILON', 'SIGMA', 'LOSS']:
            files = get_files(quantity)
            for i in range(len(files)):
                data = get_data(files[i])
                if not data:
                    sccs.append(None)
                    continue
                if quantity == quantity[0]:
                    sec_scc = sec_run.m_create(section_single_configuration_calculation)
                else:
                    sec_scc = sccs[i]
                    if sec_scc is None:
                        # This is the case when there is a mismatch between files
                        self.logger.warn(
                            'Mismatch in EXCITON and %s files' % quantity)
                        sec_scc = sec_run.m_create(section_single_configuration_calculation)
                if quantity == 'EXCITON':
                    parse_exciton(sec_scc, data)
                elif quantity == 'EPSILON':
                    parse_epsilon(sec_scc, data)
                elif quantity == 'SIGMA':
                    parse_sigma(sec_scc, data)
                elif quantity == 'LOSS':
                    parse_loss(sec_scc, data)
                else:
                    pass

    def _parse_xs_tddft(self):
        sec_run = self.archive.section_run[-1]

        fxctype = self.input_xml_parser.get('xs/tddft/fxctype', 'RPA')

        tetradf = self.input_xml_parser.get('xs/tetra/tetradf', None)
        nwacont = self.input_xml_parser.get('xs/tddft/nwacont', None)
        aresdf = self.input_xml_parser.get('xs/tddft/aresdf', True)

        file_ext_list = [
            'TET' if tetradf else None, 'AC' if nwacont else None, 'NAR' if not aresdf else None]
        file_ext = '_'.join([e for e in file_ext_list if e])

        # read q points
        qpoints = self.input_xml_parser.get('xs/qpointset/qpoint')

        def get_data(quantity, ext):
            # all files related to quantity at alll qpoints
            files = self.get_exciting_files('%s_%s%s%s.OUT' % (quantity, file_ext, ext, fxctype))
            data = [[], [], []]
            for i in range(len(qpoints)):
                data_q = []
                files_q = [f for f in files if f.endswith('QMT%s.OUT' % str(i + 1).rjust(3, '0'))]
                for f in files_q:
                    self.data_xs_parser.mainfile = f
                    if self.data_xs_parser.data is None:
                        continue
                    data_q.append(self.data_xs_parser.data)
                if not data_q:
                    continue

                n_components = len(data_q)
                n_epsilon = len(data_q[0]) // n_components
                data_q = np.transpose(np.vstack(data_q))
                data_q = np.reshape(
                    data_q, (len(data_q), n_components, n_epsilon))
                for j in range(len(data)):
                    data[j].append(data_q[j])

            return data

        for quantity in ['EPSILON', 'LOSS', 'SIGMA']:
            for ext in ['FXC', 'NLF_FXC']:
                data = get_data(quantity, ext)
                if not data[0]:
                    continue

                if quantity == 'EPSILON' and ext == 'FXC':
                    sec_scc = sec_run.m_create(section_single_configuration_calculation)
                    sec_scc.x_exciting_xs_tddft_number_of_epsilon_values = len(data[0][0][0])
                    sec_scc.x_exciting_xs_tddft_epsilon_energies = pint.Quantity(
                        data[0][0][0], 'hartree')
                    sec_scc.x_exciting_xs_tddft_dielectric_function_local_field = data[1:]

                elif quantity == 'EPSILON' and ext == 'NLF_FXC':
                    sec_scc.x_exciting_xs_tddft_dielectric_function_no_local_field = data[1:3]

                elif quantity == 'LOSS' and ext == 'FXC':
                    sec_scc.x_exciting_xs_tddft_loss_function_local_field = data[1]

                elif quantity == 'LOSS' and ext == 'NLF_FXC':
                    sec_scc.x_exciting_xs_tddft_loss_function_no_local_field = data[1]

                elif quantity == 'SIGMA' and ext == 'FXC':
                    sec_scc.x_exciting_xs_tddft_sigma_local_field = data[1:3]

                elif quantity == 'SIGMA' and ext == 'NLF_FXC':
                    sec_scc.x_exciting_xs_tddft_sigma_no_local_field = data[1:3]

    def parse_xs(self):
        sec_run = self.archive.section_run[-1]

        xs_info_files = self.get_exciting_files('INFOXS.OUT')

        if not xs_info_files:
            return

        self._calculation_type = 'xs'
        # inconsistency in the naming convention for xs input xml file
        sec_method = sec_run.m_create(section_method)

        sec_method_to_method_refs = sec_method.m_create(section_method_to_method_refs)
        sec_method_ref = self.archive.section_run[-1].section_method[0]
        sec_method_to_method_refs.method_to_method_ref = sec_method_ref
        sec_method_to_method_refs.method_to_method_kind = 'starting_point'

        self.parse_file('input.xml', sec_method)

        # parse properties
        input_file = self.get_exciting_files('input.xml')
        if not input_file:
            return
        self.input_xml_parser.mainfile = input_file[0]
        xstype = self.input_xml_parser.get('xs/xstype', None)
        if xstype.lower() == 'bse':
            self._parse_xs_bse()
        elif xstype.lower() == 'tddft':
            self._parse_xs_tddft()

    def _parse_input_gw(self, sec_method):
        gmaxvr = self.info_parser.get('x_exciting_gmaxvr', 0)
        sec_method.gw_core_treatment = self.input_xml_parser.get(
            'gw/coreflag', 'all')
        sec_method.gw_polarizability_number_of_empty_states = int(
            self.input_xml_parser.get('gw/nempty', 0))
        sec_method.gw_ngridq = self.input_xml_parser.get('gw/ngridq', [1, 1, 1])
        sec_method.gw_basis_set = 'mixed'
        sec_method.gw_qp_equation_treatment = 'linearization'
        sec_method.gw_max_frequency = self.input_xml_parser.get(
            'gw/freqgrid/freqmax', 1.0)
        sec_method.gw_frequency_grid_type = self.input_xml_parser.get(
            'gw/freqgrid/fgrid', 'gaule2')
        sec_method.gw_number_of_frequencies = self.input_xml_parser.get(
            'gw/freqgrid/nomeg', 16)
        sec_method.gw_self_energy_c_number_of_poles = self.input_xml_parser.get(
            'gw/selfenergy/npol', 0)
        sec_method.gw_self_energy_c_number_of_empty_states = self.input_xml_parser.get(
            'gw/selfenergy/nempty', 0)
        sec_method.gw_self_energy_singularity_treatment = self.input_xml_parser.get(
            'gw/selfenergy/singularity', 'mpd')
        sec_method.gw_self_energy_c_analytical_continuation = self.input_xml_parser.get(
            'gw/selfenergy/actype', 'pade')
        sec_method.gw_mixed_basis_lmax = self.input_xml_parser.get(
            'gw/mixbasis/lmaxmb', 3)
        sec_method.gw_mixed_basis_tolerance = self.input_xml_parser.get(
            'gw/mixbasis/epsmb', 0.0001)
        gmb = self.input_xml_parser.get('gw/mixbasis/gmb', 1.0)
        sec_method.gw_mixed_basis_gmax = gmb * self.info_parser.get('x_exciting_gmaxvr')
        pwm = self.input_xml_parser.get('gw/barecoul/pwm', 2.0)
        sec_method.gw_bare_coulomb_gmax = pwm * gmb * gmaxvr
        sec_method.gw_bare_coulomb_cutofftype = self.input_xml_parser.get(
            'gw/barecoul/cutofftype', 'none')
        sec_method.gw_screened_coulomb_volume_average = self.input_xml_parser.get(
            'gw/scrcoul/sciavtype', 'isotropic')
        sec_method.gw_screened_Coulomb = self.input_xml_parser.get(
            'gw/scrcoul/scrtype', 'rpa')

    def parse_gw(self):
        sec_run = self.archive.section_run[-1]

        # two versions of gw info files
        gw_info_files = ['GW_INFO.OUT', 'GWINFO.OUT']
        for f in gw_info_files:
            if self.get_exciting_files(f):
                self._calculation_type = 'gw'
                gw_info_file = f
                break

        if not self._calculation_type == 'gw':
            return

        sec_method = sec_run.m_create(section_method)
        sec_method.electronic_structure_method = 'G0W0'
        xc_functional_name = ' '.join(self.info_parser.get_xc_functional_name())
        sec_method.gw_starting_point = xc_functional_name
        sec_method_to_method_refs = sec_method.m_create(section_method_to_method_refs)
        sec_method_ref = self.archive.section_run[-1].section_method[0]
        sec_method_to_method_refs.method_to_method_ref = sec_method_ref
        sec_method_to_method_refs.method_to_method_kind = 'starting_point'

        # parse input xml file, there seems to be two versions, input_gw.xml and input-gw.xml
        for f in ['input_gw.xml', 'input-gw.xml']:
            self.parse_file(f, sec_method)

        sec_scc = sec_run.m_create(section_single_configuration_calculation)
        sec_calc_to_calc_refs = sec_scc.m_create(section_calculation_to_calculation_refs)
        sec_scc_ref = sec_run.section_single_configuration_calculation[0]
        sec_calc_to_calc_refs.calculation_to_calculation_ref = sec_scc_ref
        sec_calc_to_calc_refs.calculation_to_calculation_kind = 'starting_point'

        # parse properties
        gw_files = [
            'EVALQP.DAT', 'EVALQP.TXT', 'TDOS-QP.OUT', 'bandstructure-qp.dat',
            'BAND-QP.OUT']
        for f in gw_files:
            self.parse_file(f, sec_scc)

        gw_info_files = self.get_exciting_files(gw_info_file)
        if len(gw_info_files) > 1:
            self.logger.warn('Found multiple GW info files, will read only first!')

        self.info_gw_parser.mainfile = gw_info_files[0]

        frequency_data = self.info_gw_parser.get('frequency_data', None)
        if frequency_data is not None:
            number = frequency_data.get('number')
            sec_method.gw_number_of_frequencies = len(number)
            sec_method.gw_frequency_number = number
            sec_method.gw_frequency_values = frequency_data.get('values')
            sec_method.gw_frequency_weights = frequency_data.get('weights')

        fermi_energy = self.info_gw_parser.get('fermi_energy', None)
        if fermi_energy is not None:
            sec_scc.gw_fermi_energy = fermi_energy

        fundamental_band_gap = self.info_gw_parser.get('direct_band_gap', None)
        if fundamental_band_gap is None:
            fundamental_band_gap = self.info_gw_parser.get('fundamental_band_gap', None)
        if fundamental_band_gap is not None:
            sec_scc.gw_fundamental_gap = fundamental_band_gap

        optical_band_gap = self.info_gw_parser.get('optical_band_gap', None)
        if optical_band_gap is not None:
            sec_scc.gw_optical_gap = optical_band_gap

    def parse_miscellaneous(self):
        sec_run = self.archive.section_run[-1]
        sec_sampling_method = sec_run.m_create(section_sampling_method)

        # TODO there should be a sampling method single_point_calculation
        sec_sampling_method.sampling_method = 'geometry_optimization'

        threshold_force = self.info_parser.get('force_convergence_optimization', [None])[0]
        if threshold_force is not None:
            sec_sampling_method.geometry_optimization_threshold_force = threshold_force

        sec_frame_sequence = sec_run.m_create(section_frame_sequence)
        sec_frame_sequence.number_of_frames_in_sequence = max(self.info_parser.n_optimization_steps, 1)
        sec_frame_sequence.frame_sequence_local_frames_ref = sec_run.section_single_configuration_calculation
        sec_frame_sequence.frame_sequence_to_sampling_ref = sec_sampling_method

    def parse_method(self):
        sec_run = self.archive.section_run[-1]
        sec_method = sec_run.m_create(section_method)

        sec_method.electronic_structure_method = 'DFT'

        smearing_kind_map = {
            'Gaussian': 'gaussian', 'Methfessel-Paxton': 'methfessel-paxton',
            'Fermi-Dirac': 'fermi', 'Extended': 'tetrahedra'}

        smearing_kind = self.info_parser.get('smearing_kind', None)
        if smearing_kind is not None:
            if not isinstance(smearing_kind, str):
                smearing_kind = smearing_kind[0]
            smearing_kind = smearing_kind_map[smearing_kind]
            sec_method.smearing_kind = smearing_kind
        smearing_width = self.info_parser.get('smearing_width', None)
        if smearing_width is not None:
            sec_method.smearing_width = smearing_width

        for name in self.info_parser._convergence_keys_mapping.keys():
            threshold = self.info_parser.get_scf_threshold(name)
            if threshold is None:
                continue

            metainfo_name = 'x_exciting_scf_threshold_%s_change' % name.split('_')[-2]
            setattr(sec_method, metainfo_name, threshold)
            # additionally, set threshold to global metainfo. This is killing me!
            if metainfo_name == 'x_exciting_scf_threshold_energy_change':
                setattr(sec_method, metainfo_name.lstrip('x_exciting_'), threshold)

        xc_functional_names = self.info_parser.get_xc_functional_name()
        if not xc_functional_names:
            # get it from input.xml
            input_file = self.get_exciting_files('input.xml')
            for f in input_file:
                self.input_xml_parser.mainfile = f
                correlation = self.input_xml_parser.get('libxc/correlation', None)
                xc_functional_names.append(correlation)
                exchange = self.input_xml_parser.get('libxc/exchange', None)
                xc_functional_names.append(exchange)

        for name in xc_functional_names:
            if name is None:
                continue
            sec_xc_functional = sec_method.m_create(section_XC_functionals)
            sec_xc_functional.XC_functional_name = name

        sec_method.number_of_spin_channels = self.info_parser.get_number_of_spin_channels()

        if self._calculation_type == 'volume_optimization':
            sec_method.x_exciting_volume_optimization = True

    def parse_scc_full(self, loop_type):
        sec_run = self.archive.section_run[-1]

        # total energy
        total_energy = self.info_parser.get('energy_total_%s' % loop_type)
        if total_energy is None:
            # for calculations which does not output the final scf iteration
            # we take the last scf_iteration as the final scf iteration
            # we do not however read all the energy contributions
            if loop_type == 'final_scf_iteration':
                total_energy = [self.info_parser.get('energy_total_scf_iteration')[-1]]

        if total_energy is None:
            return

        if loop_type == 'scf_iteration':
            metainfo_ext = '_scf_iteration'
            sec_scc = sec_run.section_single_configuration_calculation[0]
            section = sec_scc.m_create(section_scf_iteration)
            index = len(sec_run.section_single_configuration_calculation[0].section_scf_iteration) - 1
        else:
            metainfo_ext = ''
            section = sec_run.m_create(section_single_configuration_calculation)
            index = 0

        section.energy_total = total_energy[index]

        # energy contibutions
        energy_contributions = self.info_parser.get(
            'energy_contributions_%s' % loop_type, [{}] * (index + 1))[index]

        for key, names in self.info_parser._energy_keys_mapping.items():
            val = None
            for name in names:
                val = energy_contributions.get(name, None)
                if val:
                    break
            if val is None:
                continue
            setattr(section, key + metainfo_ext, val)
            if key == 'x_exciting_fermi_energy':
                # set it also in the global fermi energy, this is killing me
                # there should only be one in global
                # the metainfo naming is not consistent for scf_iteration
                # and for global it becomes energy_reference_fermi
                key = 'energy_reference_fermi'
                metainfo_ext_fermi = '_iteration' if loop_type == 'scf_iteration' else metainfo_ext
                val = pint.Quantity(
                    [val.magnitude] * self.info_parser.get_number_of_spin_channels(), 'hartree')
                setattr(section, key + metainfo_ext_fermi, val)

        # charge contributions
        charge_contributions = self.info_parser.get(
            'charge_contributions_%s' % loop_type, [{}] * (index + 1))[index]
        for key, names in self.info_parser._electron_charge_keys_mapping.items():
            val = None
            for name in names:
                val = charge_contributions.get(name, None)
                if val is not None:
                    break
            if val is None:
                continue
            if key == 'x_exciting_section_MT_charge_atom':
                for n in range(len(val)):
                    sec_mt_charge_atom = section.m_create(x_exciting_section_MT_charge_atom)
                    sec_mt_charge_atom.x_exciting_MT_charge_atom_index = n + 1
                    sec_mt_charge_atom.x_exciting_MT_charge_atom_symbol = val[n][0]
                    sec_mt_charge_atom.x_exciting_MT_charge_atom_value = val[n][1]
            else:
                setattr(section, key + metainfo_ext, val)

        # moment contributions
        moment_contributions = self.info_parser.get(
            'moment_contributions_%s' % loop_type, [{}] * (index + 1))[index]
        for key, names in self.info_parser._moment_keys_mapping.items():
            val = None
            for name in names:
                val = moment_contributions.get(name, None)
                if val is not None:
                    break
            if val is None:
                continue
            if key == 'x_exciting_section_MT_moment_atom':
                for n in range(len(val)):
                    sec_mt_moment_atom = section.m_create(x_exciting_section_MT_moment_atom)
                    sec_mt_moment_atom.x_exciting_MT_moment_atom_index = n + 1
                    sec_mt_moment_atom.x_exciting_MT_moment_atom_symbol = val[n][0]
                    sec_mt_moment_atom.x_exciting_MT_moment_atom_value = val[n][1]
            else:
                setattr(section, key + metainfo_ext, val)

        # forces
        forces = self.info_parser.get('atom_forces_%s' % loop_type, [None] * (index + 1))[index]
        if forces is not None:
            section.atom_forces = forces

        # convergence values
        for name in self.info_parser._convergence_keys_mapping.keys():
            if loop_type == 'scf_iteration':
                val = self.info_parser.get_scf_quantity(name)
            else:
                val = self.info_parser.get('%s_%s' % (name, loop_type), None)

            if val is None or val[index] is None:
                continue

            setattr(section, '%s_%s' % (name, loop_type), val[index][0])

        # other metainfo
        for name in self.info_parser._miscellaneous_keys_mapping.keys():
            if loop_type == 'scf_iteration':
                val = self.info_parser.get_scf_quantity(name)
            else:
                val = self.info_parser.get('%s_%s' % (name, loop_type), None)

            if val is None or val[index] is None:
                continue

            val = val[index]
            if name == 'time':
                if loop_type == 'scf_iteration':
                    section.time_scf_iteration = val
                else:
                    section.time_calculation = val
            else:
                setattr(section, '%s_%s' % (name, loop_type), val)

        # references to method and system
        if loop_type != 'scf_iteration':
            section.single_configuration_calculation_to_system_ref = sec_run.section_system[-1]
            section.single_configuration_to_calculation_method_ref = sec_run.section_method[-1]

        # only final_scf_iteration contains scf_iteration steps
        if loop_type != 'final_scf_iteration':
            return

        energy_scf_iteration = self.info_parser.get('energy_total_scf_iteration')
        for n in range(len(energy_scf_iteration)):
            self.parse_scc_full('scf_iteration')

    def parse_system_full(self, loop_type):
        sec_run = self.archive.section_run[-1]

        if loop_type == 'optimization':
            index = len(sec_run.section_system) - 1
        else:
            index = 0

        positions = self.info_parser.get_atom_positions(loop_type, index)
        if positions is None:
            return

        sec_system = sec_run.m_create(section_system)

        sec_system.atom_labels = self.info_parser.get_atom_labels()
        sec_system.atom_positions = positions
        sec_system.configuration_periodic_dimensions = [True] * 3
        sec_system.lattice_vectors = self.info_parser.get('lattice_vectors')
        sec_system.simulation_cell = self.info_parser.get('lattice_vectors')
        # I did not add it to x_exciting_simulation_reciprocal_cell because the
        # metainfo has wrong units?
        # sec_system.lattice_vectors_reciprocal = self.info_parser.get('lattice_vectors_reciprocal')

        if loop_type != 'final_scf_iteration':
            return

        for name in self.info_parser._system_keys_mapping.keys():
            val = self.info_parser.get(name, None)
            if val is None:
                continue

            if name == 'x_exciting_spin_treatment':
                sub_sec = sec_system.m_create(x_exciting_section_spin)
                sub_sec.x_exciting_spin_treatment = val
            elif name == 'x_exciting_xc_functional':
                sub_sec = sec_system.m_create(x_exciting_section_xc)
                sub_sec.x_exciting_xc_functional
            elif name == 'x_exciting_species_rtmin':
                setattr(sec_system, name, ' '.join([str(v) for v in val]))
            else:
                setattr(sec_system, name, val)

        # clathrate info
        clathrate_file = self.get_exciting_files('str.out')
        if clathrate_file:
            sec_system.x_exciting_clathrates = True
            self.data_clathrate_parser.mainfile = clathrate_file[0]
            if self.data_clathrate_parser.data:
                data = np.transpose(self.data_clathrate_parser.data)
                sec_system.x_exciting_clathrates_atom_coordinates = np.transpose(
                    np.array(data[:3], dtype=float))
                sec_system.x_exciting_clathrates_atom_labels = list(data[3])
        else:
            sec_system.x_exciting_clathrates = False

    def parse_configurations(self):
        sec_run = self.archive.section_run[-1]

        # Initial ground state configuration
        self.parse_system_full('final_scf_iteration')
        self.parse_scc_full('final_scf_iteration')

        # Add data to scc
        # TODO add support for more output files and properties
        exciting_files = [
            'dos.xml', 'bandstructure.xml', 'EIGVAL.OUT', 'FERMISURF.bxsf', 'FS.bxsf']
        for f in exciting_files:
            self.parse_file(f, sec_run.section_single_configuration_calculation[-1])

        # Configurations at each optimization step
        for n in range(self.info_parser.n_optimization_steps):
            self.parse_system_full('optimization')

            sec_scc = sec_run.m_create(section_single_configuration_calculation)
            sec_scc.energy_total = self.info_parser.get('total_energy_optimization')[n]
            forces = self.info_parser.get('atom_forces_optimization', [None] * (n + 1))[n]
            if forces is not None:
                sec_scc.atom_forces = forces

            method = self.info_parser.get('method_optimization', [None] * (n + 1))[n]
            if method is not None:
                sec_scc.x_exciting_geometry_optimization_method = method

            nstep = self.info_parser.get('nstep_optimization', [None] * (n + 1))[n]
            if nstep is not None:
                sec_scc.x_exciting_geometry_optimization_step = nstep

            force_convergence = self.info_parser.get('force_convergence_optimization', [None] * (n + 1))[n]
            if force_convergence is not None:
                sec_scc.x_exciting_maximum_force_magnitude = force_convergence[0]
                sec_scc.x_exciting_geometry_optimization_threshold_force = force_convergence[1]

        # Final configuration after optimization
        self.parse_system_full('final_optimization')
        self.parse_scc_full('final_optimization')

        # Volume optimizations
        volume_index = 1
        while True:
            info_volume = self.get_exciting_files('run_dir%s/INFO.OUT' % str(volume_index).rjust(2, '0'))
            if not info_volume:
                break
            sec_calc_to_calc_refs = sec_scc.m_create(section_calculation_to_calculation_refs)
            sec_calc_to_calc_refs.calculation_to_calculation_external_url = info_volume[0]
            sec_calc_to_calc_refs.calculation_to_calculation_kind = 'source_calculation'
            self._calculation_type = 'volume_optimization'

    def _init_parsers(self):
        self.info_parser.mainfile = self.filepath
        self.info_parser.logger = self.logger
        self.dos_parser.logger = self.logger
        self.bandstructure_parser.logger = self.logger
        self.eigval_parser.logger = self.logger
        self.fermisurf_parser.logger = self.logger
        self.evalqp_parser.logger = self.logger
        self.dos_out_parser.logger = self.logger
        self.bandstructure_dat_parser.logger = self.logger
        self.band_out_parser.logger = self.logger
        self.info_gw_parser.logger = self.logger
        self.input_xml_parser.logger = self.logger
        self.data_xs_parser.logger = self.logger
        self.data_clathrate_parser.logger = self.logger

    def parse(self, filepath, archive, logger):
        self.filepath = filepath
        self.archive = archive
        self.logger = logger if logger is not None else logging

        self._calculation_type = None

        self._init_parsers()

        sec_run = self.archive.m_create(section_run)

        sec_run.program_name = 'exciting'
        program_version = self.info_parser.get('program_version')
        if isinstance(program_version, list):
            program_version = ' '.join(program_version)
        sec_run.program_version = program_version
        sec_run.program_basis_set_type = '(L)APW+lo'

        # method goes first since reference needed for sec_scc
        self.parse_method()

        self.parse_configurations()

        self.parse_gw()

        self.parse_xs()

        self.parse_miscellaneous()
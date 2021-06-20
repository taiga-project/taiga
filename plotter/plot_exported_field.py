import os
import numpy
import matplotlib.pyplot as plt
from preproc.init_from_efit import CDBReader

shot_number = '17178'
time = '1097'

global_settings = {
    'taiga_directory': 'results/' + shot_number + '_' + time,
    'reference_spline_directory': 'input/fieldGrid/' + shot_number + '_' + time,
    'reference_spline_rad_filename': 'rcord.dat',
    'reference_spline_z_filename': 'zcord.dat'
}

rad = {
    'name': 'B_rad',
    'title': 'B_{R}',
    'pycdb_name': 'B_R',
    'value_range': [-0.2, 0.2],
    'taiga_spline_filename': 'test_spline_brad.dat',
    'taiga_bspline_filename': 'test_bspline_brad.dat',
    'reference_spline_filename': 'brad.dat'
}
rad.update(global_settings)

z = {
    'name': 'B_z',
    'title': 'B_{Z}',
    'pycdb_name': 'B_Z',
    'value_range': [-0.4, 0.4],
    'taiga_spline_filename': 'test_spline_bz.dat',
    'taiga_bspline_filename': 'test_bspline_bz.dat',
    'reference_spline_filename': 'bz.dat'
}
z.update(global_settings)

tor = {
    'name': 'B_tor',
    'title': 'B_{\phi}',
    'pycdb_name': 'B_phi',
    'value_range': [0, 2],
    'taiga_spline_filename': 'test_spline_btor.dat',
    'taiga_bspline_filename': 'test_bspline_btor.dat',
    'reference_spline_filename': 'btor.dat'
}
tor.update(global_settings)

polflux = {
    'name': 'psi_n',
    'title': '\psi_\mathrm{pol}^\mathrm{n}',
    'pycdb_name': 'psi_n',
    'value_range': [0, 2],
    'taiga_spline_filename': 'test_spline_psi.dat',
    'taiga_bspline_filename': 'test_bspline_psi.dat',
    'reference_spline_filename': 'psi2.dat'
}
polflux.update(global_settings)


class MagneticFieldDataSet:
    def __init__(self):
        self.rad = []
        self.z = []
        self.value = []


class MagneticFieldDataFile:
    def __init__(self):
        self.directory = ''
        self.filename = ''
        self.path = ''

    def get_value(self):
        print('Read from: ' + self.path)
        return numpy.genfromtxt(self.path, dtype=None)

    def set_directory(self, directory):
        self.directory = directory
        self.set_path()

    def set_filename(self, filename):
        self.filename = filename
        self.set_path()

    def set_path(self):
        self.path = self.directory + '/' + self.filename


class MagneticFieldComponent:
    def __init__(self, component):
        self.component = component
        self.name = ''
        self.get_name()
        self.title = ''
        self.get_title()
        self.value_range = [0, 1]
        self.set_value_range()
        self.taiga = MagneticFieldDataSet()
        self.set_taiga()
        self.reference = MagneticFieldDataSet()
        self.set_reference()
        self.taiga_bspline = MagneticFieldDataSet()
        self.set_taiga_bspline()
        self.reference_bspline = MagneticFieldDataSet()
        self.set_reference_bspline()

    def get_name(self):
        self.name = self.component['name']

    def get_title(self):
        self.title = self.component['title']

    def set_value_range(self):
        self.value_range = self.component['value_range']

    def set_taiga(self):
        taiga_file = MagneticFieldDataFile()
        taiga_file.set_directory(self.component['taiga_directory'])
        taiga_file.set_filename(self.component['taiga_spline_filename'])
        taiga_file_content = taiga_file.get_value()
        self.taiga.rad = taiga_file_content[:, 0]
        self.taiga.z = taiga_file_content[:, 1]
        self.taiga.value = taiga_file_content[:, 2]

    def set_taiga_bspline(self):
        taiga_file = MagneticFieldDataFile()
        taiga_file.set_directory(self.component['taiga_directory'])
        taiga_file.set_filename(self.component['taiga_bspline_filename'])
        taiga_file_content = taiga_file.get_value()
        self.taiga_bspline.rad = taiga_file_content[:, 0]
        self.taiga_bspline.z = taiga_file_content[:, 1]
        self.taiga_bspline.value = taiga_file_content[:, 2]

    def set_reference(self):
        reference_data_file = MagneticFieldDataFile()
        reference_data_file.set_directory(self.component['reference_spline_directory'])
        reference_data_file.set_filename(self.component['reference_spline_filename'])
        reference_rad_file = MagneticFieldDataFile()
        reference_rad_file.set_directory(self.component['reference_spline_directory'])
        reference_rad_file.set_filename(self.component['reference_spline_rad_filename'])
        reference_z_file = MagneticFieldDataFile()
        reference_z_file.set_directory(self.component['reference_spline_directory'])
        reference_z_file.set_filename(self.component['reference_spline_z_filename'])
        self.reference.rad = reference_rad_file.get_value()
        self.reference.z = reference_z_file.get_value()
        self.reference.value = reference_data_file.get_value().T

    def set_reference_bspline(self):
        cdb = CDBReader(shot_number, time)
        self.reference_bspline.rad = self.reference.rad
        self.reference_bspline.z = self.reference.z
        attribute = self.component['pycdb_name']
        f = getattr(cdb, attribute)
        self.reference_bspline.value = f(self.reference_bspline.z, self.reference_bspline.rad, grid=True)


class PlotMagneticFieldComponent(MagneticFieldComponent):
    def __init__(self, component_name):
        super().__init__(component_name)
        plt.rc('font', family='serif')
        plt.rc('text', usetex=True)
        fig, ((ax_taiga_spline, ax_reference_spline), (ax_taiga_bspline, ax_reference_bspline)) = \
            plt.subplots(2, 2, sharex='all', sharey='all', subplot_kw=dict(aspect='equal'), figsize=[6.4, 9])
        ax_taiga_spline.tricontourf(self.taiga.rad, self.taiga.z, self.taiga.value, levels=self.get_levels())
        ax_taiga_spline.set_xlabel(r'$R$ [m]')
        ax_taiga_spline.set_ylabel(r'$Z$ [m]')
        ax_taiga_spline.set_title('TAIGA (spline)')
        ax_reference_spline.tick_params(direction='out', left=True, right=True, labelleft=True)
        c = ax_reference_spline.contourf(self.reference.rad, self.reference.z, self.reference.value, levels=self.get_levels())
        ax_reference_spline.set_xlabel(r'$R$ [m]')
        ax_reference_spline.tick_params(direction='out', left=True, right=True)
        ax_reference_spline.set_title('Reference (spline)')
        cbar_spline = fig.colorbar(c, ax=(ax_taiga_spline, ax_reference_spline), pad=0.02, fraction=0.034)
        cbar_spline.set_ticks(self.get_tick_levels())

        i = ~numpy.isnan(self.taiga_bspline.value)
        ax_taiga_bspline.tricontourf(self.taiga_bspline.rad[i], self.taiga_bspline.z[i], self.taiga_bspline.value[i], levels=self.get_levels())
        ax_taiga_bspline.set_xlabel(r'$R$ [m]')
        ax_taiga_bspline.set_ylabel(r'$Z$ [m]')
        ax_taiga_bspline.set_title('TAIGA (B-spline)')
        ax_reference_bspline.tick_params(direction='out', left=True, right=True, labelleft=True)
        c = ax_reference_bspline.contourf(self.reference_bspline.rad, self.reference_bspline.z, self.reference_bspline.value, levels=self.get_levels())
        ax_reference_bspline.set_xlabel(r'$R$ [m]')
        ax_reference_bspline.tick_params(direction='out', left=True, right=True)
        ax_reference_bspline.set_title('Reference (B-spline)')
        plt.suptitle(r'$' + self.title + '$ \@ COMPASS \#' + shot_number + ' (' + time + ' ms) ')
        cbar_bspline = fig.colorbar(c, ax=(ax_taiga_bspline, ax_reference_bspline), pad=0.02, fraction=0.034)
        cbar_bspline.set_ticks(self.get_tick_levels())
        plt.savefig(self.name + '_' + shot_number + '_' + time + '.svg')
        plt.show()

    def get_levels(self):
        return numpy.linspace(self.value_range[0], self.value_range[1], 21)

    def get_tick_levels(self):
        return numpy.linspace(self.value_range[0], self.value_range[1], 9)


if __name__ == "__main__":
    PlotMagneticFieldComponent(rad)
    PlotMagneticFieldComponent(z)
    PlotMagneticFieldComponent(tor)
    PlotMagneticFieldComponent(polflux)

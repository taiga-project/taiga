import os
import numpy
import matplotlib.pyplot as plt

shot_number = '17178'
time = '1097'

global_settings = {
    'taiga_directory': '.',
    'reference_directory': 'input/fieldGrid/'+shot_number+'_'+time,
    'reference_rad_filename': 'rcord.dat',
    'reference_z_filename': 'zcord.dat'
}

rad = {
    'name': 'B_rad',
    'value_range': [-0.2, 0.2],
    'taiga_filename': 'exported_fieldR.dat',
    'reference_filename': 'brad.dat'
}
rad.update(global_settings)

z = {
    'name': 'B_z',
    'value_range': [-0.3, 0.3],
    'taiga_filename': 'exported_fieldZ.dat',
    'reference_filename': 'bz.dat'
}
z.update(global_settings)

tor = {
    'name': 'B_tor',
    'value_range': [0, 2],
    'taiga_filename': 'exported_fieldT.dat',
    'reference_filename': 'btor.dat'
}
tor.update(global_settings)


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
        self.value_range = [0, 1]
        self.set_value_range()
        self.taiga = MagneticFieldDataSet()
        self.set_taiga()
        self.reference = MagneticFieldDataSet()
        self.set_reference()

    def get_name(self):
        self.name = self.component['name']

    def set_value_range(self):
        self.value_range = self.component['value_range']

    def set_taiga(self):
        taiga_file = MagneticFieldDataFile()
        taiga_file.set_directory(self.component['taiga_directory'])
        taiga_file.set_filename(self.component['taiga_filename'])
        taiga_file_content = taiga_file.get_value()
        self.taiga.rad = taiga_file_content[:, 0]
        self.taiga.z = taiga_file_content[:, 1]
        self.taiga.value = taiga_file_content[:, 2]

    def set_reference(self):
        reference_data_file = MagneticFieldDataFile()
        reference_data_file.set_directory(self.component['reference_directory'])
        reference_data_file.set_filename(self.component['reference_filename'])
        reference_rad_file = MagneticFieldDataFile()
        reference_rad_file.set_directory(self.component['reference_directory'])
        reference_rad_file.set_filename(self.component['reference_rad_filename'])
        reference_z_file = MagneticFieldDataFile()
        reference_z_file.set_directory(self.component['reference_directory'])
        reference_z_file.set_filename(self.component['reference_z_filename'])
        self.reference.rad = reference_rad_file.get_value()
        self.reference.z = reference_z_file.get_value()
        self.reference.value = reference_data_file.get_value().T


class PlotMagneticFieldComponent(MagneticFieldComponent):
    def __init__(self, component_name):
        super().__init__(component_name)
        fig, (ax_taiga, ax_reference) = plt.subplots(1, 2, sharex='all', sharey='all')
        ax_taiga.tricontourf(self.taiga.rad, self.taiga.z, self.taiga.value, levels=self.get_levels())
        ax_taiga.set_xlabel('R [m]')
        ax_taiga.set_ylabel('Z [m]')
        ax_taiga.set_title('TAIGA')
        ax_reference.tick_params(direction='out', left=True, right=True, labelleft=True)
        ax_reference.contourf(self.reference.rad, self.reference.z, self.reference.value, levels=self.get_levels())
        ax_reference.set_xlabel('R [m]')
        ax_reference.tick_params(direction='out', left=True, right=True, labelright=True)
        ax_reference.set_title('Reference')
        plt.suptitle(self.name+' @ Compass #'+shot_number+' ('+time+' ms) ')
        plt.show()

    def get_levels(self):
        return numpy.linspace(self.value_range[0], self.value_range[1], 20)


if __name__ == "__main__":
    PlotMagneticFieldComponent(rad)
    PlotMagneticFieldComponent(z)
    PlotMagneticFieldComponent(tor)

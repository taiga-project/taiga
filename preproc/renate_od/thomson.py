import h5py

from preproc.renate_od.utils import *


class ThomsonProfiles:
    def __init__(self, thomson_directory, export_directory, shot_number, time, reconstruction_id=1):
        self.thomson_directory = thomson_directory
        self.reconstruction_id = reconstruction_id
        self.shot_number = shot_number
        self.time = time
        self.export_directory = export_directory

        self.time_index = []
        self.z_axis = []
        self.temperature_profile = []
        self.temperature_error_profile = []
        self.density_profile = []
        self.density_error_profile = []
        self.normalised_poloidal_flux_profile = []

        self.read_thomson_database()
        self.density = ProfileManager(export_directory, self.normalised_poloidal_flux_profile,
                                      self.density_profile, self.density_error_profile, 'density')
        self.temperature = ProfileManager(export_directory, self.normalised_poloidal_flux_profile,
                                          self.temperature_profile, self.temperature_error_profile, 'temperature')

    def read_thomson_database(self):
        time_dataset = self.get_dataset('TS_record_time', reconstruction_id=self.reconstruction_id)
        self.z_axis = self.get_dataset('TS_z_axis', reconstruction_id=self.reconstruction_id)
        print('Thomson scattering time and geometry files read successfully from: ' + self.thomson_directory)
        self.time_index = self.get_time_index(time_dataset)
        self.temperature_profile = self.get_profile('Te', reconstruction_id=self.reconstruction_id)
        self.temperature_error_profile = self.get_profile('Te_err', reconstruction_id=self.reconstruction_id)
        self.density_profile = self.get_profile('ne', reconstruction_id=self.reconstruction_id)
        self.density_error_profile = self.get_profile('ne_err', reconstruction_id=self.reconstruction_id)
        self.normalised_poloidal_flux_profile = self.get_profile('psi_n', reconstruction_id=1)
        print('Thomson scattering profile files read successfully from: ' + self.thomson_directory)

    def get_dataset(self, field, reconstruction_id=1):
        file = h5py.File(self.thomson_directory + '/' + field + '.' + str(reconstruction_id) + '.h5', 'r')
        return file[field][()]

    def get_time_index(self, time_dataset):
        return numpy.nanargmin((numpy.abs(time_dataset - int(self.time))))

    def get_profile(self, field, reconstruction_id=1):
        return self.get_dataset(field, reconstruction_id)[self.time_index]

    def plot_profiles(self):
        self.density.plot_profile(self.shot_number, self.time, 1e-19)
        self.temperature.plot_profile(self.shot_number, self.time)

    def export_profiles(self):
        self.density.export_profile()
        self.temperature.export_profile()

import h5py

from utils import *


class ThomsonProfiles:
    def __init__(self, thomson_directory, shot_number, time, reconstruction_id, ts_time_source):
        self.thomson_directory = thomson_directory
        self.reconstruction_id = reconstruction_id
        self.ts_time_source = ts_time_source
        self.shot_number = shot_number
        self.time = time

        self.time_index = []
        self.z_axis = []
        self.temperature_profile = []
        self.temperature_error_profile = []
        self.density_profile = []
        self.density_error_profile = []
        self.normalised_poloidal_flux_profile = []

        self.read_thomson_database()
        self.density = ProfileManager(self.normalised_poloidal_flux_profile,
                                      self.density_profile, self.density_error_profile)
        self.temperature = ProfileManager(self.normalised_poloidal_flux_profile,
                                          self.temperature_profile, self.temperature_error_profile)

    def read_thomson_database(self):
        try:
            time_dataset = self.get_dataset('TS_' + self.ts_time_source + '_time', reconstruction_id=1)
            self.z_axis = self.get_dataset('TS_z_axis', self.reconstruction_id)
            print('Thomson scattering time and geometry files read successfully from: ' + self.thomson_directory)
        except OSError:
            raise OSError('Invalid Thomson scattering file structure! '
                          '\nPlease check the directory tree here:\n' + self.thomson_directory)
        except KeyError:
            raise KeyError(
                'Invalid Thomson scattering data structure!\nExample for a correct structure:\nTe.1.h5\n\tTe')

        self.time_index = self.get_time_index(time_dataset)

        try:
            self.temperature_profile = self.get_profile('Te', reconstruction_id=self.reconstruction_id)
            self.temperature_error_profile = self.get_profile('Te_err', reconstruction_id=self.reconstruction_id)
            self.density_profile = self.get_profile('ne', reconstruction_id=self.reconstruction_id)
            self.density_error_profile = self.get_profile('ne_err', reconstruction_id=self.reconstruction_id)
            self.normalised_poloidal_flux_profile = self.get_profile('psi_n', reconstruction_id=1)
            print('Thomson scattering profile files read successfully from: ' + self.thomson_directory)
        except OSError:
            raise OSError('Invalid Thomson scattering file structure! '
                          '\nPlease check the directory tree here:\n' + self.thomson_directory)
        except KeyError:
            raise KeyError(
                'Invalid Thomson scattering data structure!\nExample for a correct structure:\nTe.1.h5\n\tTe')

    def plot_profiles(self):
        self.density.plot_profile(self.shot_number, self.time, r'$n_e~(\mathrm{m}^{-3})$')
        self.temperature.plot_profile(self.shot_number, self.time, r'$T_e~(\mathrm{keV})$')

    def get_time_index(self, time_dataset):
        return numpy.nanargmin((numpy.abs(time_dataset - int(self.time))))

    def get_dataset(self, field, reconstruction_id):
        file = h5py.File(self.thomson_directory + '/' + field + '.' + str(reconstruction_id) + '.h5')
        return file[field][()]

    def get_profile(self, field, reconstruction_id):
        return self.get_dataset(field, reconstruction_id)[self.time_index]

    def export_profiles(self, path):
        export_directory = path + '/' + self.shot_number + '_' + self.time
        self.density.export_profile(export_directory, 'density')
        self.temperature.export_profile(export_directory, 'temperature')

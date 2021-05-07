import os
import h5py
import matplotlib
import numpy
import scipy
from statsmodels.nonparametric.smoothers_lowess import lowess

from stefanikova2016 import *


class Profiles:
    def __init__(self, shot_number='17178', time='1097', reconstruction_id=2,
                 database_directory='input/cdb', thomson_subdir='THOMSON', efit_subdir='EFITXX',
                 ts_time_source='EFIT'):
        self.thomson_profiles = ThomsonProfiles(shot_number, time, reconstruction_id,
                                                database_directory, thomson_subdir, efit_subdir,
                                                ts_time_source)
        self.distance = []
        self.density = []
        self.temperature = []

    def get_distance(self):
        return self.distance

    def get_density(self):
        return self.density

    def get_temperature(self):
        return self.temperature


class MockedProfiles(Profiles):
    def __init__(self):
        super().__init__()
        self.distance = [0, 0.1, 0.2]
        self.density = [1e19, 1.2e19, 1.5e19]
        self.temperature = [1000, 1500, 2500]


class ThomsonProfiles:
    def __init__(self, shot_number='17178', time='1097', reconstruction_id=2,
                 database_directory='input/cdb', thomson_subdir='THOMSON', efit_subdir='EFITXX',
                 ts_time_source='EFIT'):
        try:
            home_directory = os.environ['TAIGA_HOME']
        except KeyError:
            home_directory = '/home/matyi/work/taiga_local'  # '.'

        self.data_directory = home_directory + '/' + database_directory + '/' + str(shot_number)
        self.thomson_directory = self.data_directory + '/' + thomson_subdir
        self.efit_file = self.data_directory + '/' + efit_subdir + '/' + 'EFITXX.' + str(reconstruction_id) + '.h5'
        self.reconstruction_id = reconstruction_id
        self.ts_time_source = ts_time_source
        self.time = time

        self.time_index = []
        self.temperature_profile = []
        self.temperature_error_profile = []
        self.density_profile = []
        self.density_error_profile = []
        self.normalised_poloidal_flux_profile = []

        self.read_thomson_database()
        self.density_manager = ProfileManager(self.normalised_poloidal_flux_profile,
                                              self.density_profile, self.density_error_profile)
        self.temperature_manager = ProfileManager(self.normalised_poloidal_flux_profile,
                                                  self.temperature_profile, self.temperature_error_profile)
        self.plot_profiles()

    def read_thomson_database(self):
        try:
            time_dataset = self.get_ts_dataset('TS_' + self.ts_time_source + '_time', reconstruction_id=1)
            z_axis = self.get_ts_dataset('TS_z_axis', self.reconstruction_id)
            print('Thomson scattering time and geometry files read successfully from: ' + self.thomson_directory)
        except OSError:
            raise OSError('Invalid Thomson scattering file structure! '
                          '\nPlease check the directory tree here:\n' + self.thomson_directory)
        except KeyError:
            raise KeyError(
                'Invalid Thomson scattering data structure!\nExample for a correct structure:\nTe.1.h5\n\tTe')

        self.time_index = self.get_ts_time_index(time_dataset)

        try:
            self.temperature_profile = self.get_ts_profile('Te', reconstruction_id=self.reconstruction_id)
            self.temperature_error_profile = self.get_ts_profile('Te_err', reconstruction_id=self.reconstruction_id)
            self.density_profile = self.get_ts_profile('ne', reconstruction_id=self.reconstruction_id)
            self.density_error_profile = self.get_ts_profile('ne_err', reconstruction_id=self.reconstruction_id)
            self.normalised_poloidal_flux_profile = self.get_ts_profile('psi_n', reconstruction_id=1)
            print('Thomson scattering profile files read successfully from: ' + self.thomson_directory)
        except OSError:
            raise OSError('Invalid Thomson scattering file structure! '
                          '\nPlease check the directory tree here:\n' + self.thomson_directory)
        except KeyError:
            raise KeyError(
                'Invalid Thomson scattering data structure!\nExample for a correct structure:\nTe.1.h5\n\tTe')

    def plot_profiles(self):
        self.density_manager.plot_profile()
        self.temperature_manager.plot_profile()

    def get_ts_time_index(self, time_dataset):
        return (numpy.abs(time_dataset - int(self.time))).argmin()

    def get_ts_dataset(self, field, reconstruction_id):
        file = h5py.File(self.thomson_directory + '/' + field + '.' + str(reconstruction_id) + '.h5')
        return file[field][()]

    def get_ts_profile(self, field, reconstruction_id):
        return self.get_ts_dataset(field, reconstruction_id)[self.time_index]


def set_negatives_to_zero(values):
    return values * (values > 0)


class ProfileManager:
    def __init__(self, x, y, y_error):
        self.x_raw = x
        self.y_raw = y
        self.y_error_raw = y_error
        self.x = x
        self.y = y
        self.y_error = y_error
        self.y_smooth = []
        self.f_smooth = ValueError

        self.x_fine = []
        self.y_fine = []

        self.get_valid_profile()
        self.get_smoothed_profile()

    def get_valid_profile(self):
        valid_indices = self.get_valid_profile_indices()
        self.x = self.x_raw[valid_indices]
        self.y = self.y_raw[valid_indices]
        self.y_error = self.y_error_raw[valid_indices]

    def get_valid_profile_indices(self):
        return numpy.where(~numpy.isnan(self.y) & (self.x > 0.0) & (self.x < 1.2) &
                           (self.y_error < self.y) & self.filter_sol_outliers())

    def filter_sol_outliers(self):
        return (numpy.fmin.accumulate(self.y) == self.y) | (self.x < 1.02)

    def get_smoothed_profile(self):
        self.smooth_profile()
        self.refine_smoothed_profile()

    def smooth_profile(self):
        smooth = lowess(self.y, self.x, is_sorted=True, frac=0.3, it=0)
        y_smooth = smooth[:, 1]
        self.y_smooth = set_negatives_to_zero(y_smooth)

    def refine_smoothed_profile(self):
        self.interpolate_smoothed_profile()
        self.x_fine = numpy.linspace(self.x[0], 1.2, 1000)
        self.y_fine = self.f_smooth(self.x_fine)

    def interpolate_smoothed_profile(self):
        x_ext = numpy.append(self.x, [self.x[-1] + 0.001, 1.2])
        y_ext = numpy.append(self.y_smooth, [0., 0.])
        self.f_smooth = scipy.interpolate.interp1d(x_ext, y_ext, kind='linear', bounds_error=False)

    def plot_profile(self):
        fig, ax = matplotlib.pyplot.subplots()
        ax.plot(self.x_raw, self.y_raw, '.')
        ax.plot(self.x, self.y, '.')
        ax.fill_between(self.x, self.y - 1 * self.y_error, self.y + 1 * self.y_error, color='gray', alpha=0.4)
        ax.fill_between(self.x, self.y - 3 * self.y_error, self.y + 3 * self.y_error, color='gray', alpha=0.3)
        ax.fill_between(self.x, self.y - 5 * self.y_error, self.y + 5 * self.y_error, color='gray', alpha=0.2)
        ax.plot(self.x_fine, self.y_fine)
        matplotlib.pyplot.show()

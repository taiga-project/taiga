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

        self.read_thomson_database()
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
        self.plot_profile(self.normalised_poloidal_flux_profile, self.temperature_profile, self.temperature_error_profile)
        self.plot_profile(self.normalised_poloidal_flux_profile, self.density_profile, self.density_error_profile)

    def get_ts_time_index(self, time_dataset):
        return (numpy.abs(time_dataset - int(self.time))).argmin()

    def get_ts_dataset(self, field, reconstruction_id):
        file = h5py.File(self.thomson_directory + '/' + field + '.' + str(reconstruction_id) + '.h5')
        return file[field][()]

    def get_ts_profile(self, field, reconstruction_id):
        return self.get_ts_dataset(field, reconstruction_id)[self.time_index]

    def filter_sol_outliers(self, x, y):
        return (numpy.fmin.accumulate(y) == y) | (x < 1.02)

    def get_valid_profile_indices(self, x, y, yerr):
        return numpy.where(~numpy.isnan(y) & (x > 0.0) & (x < 1.2) &
                           (yerr < y) & self.filter_sol_outliers(x, y))

    def get_valid_profile(self, x_in, y_in, yerr_in):
        valid_indices = self.get_valid_profile_indices(x_in, y_in, yerr_in)
        x = x_in[valid_indices]
        y = y_in[valid_indices]
        yerr = yerr_in[valid_indices]
        return x, y, yerr

    def set_negatives_to_zero(self, x):
        return x * (x > 0)

    def smooth_profile(self, x, y):
        smooth = lowess(y, x, is_sorted=True, frac=0.3, it=0)
        y_smooth = smooth[:, 1]
        return self.set_negatives_to_zero(y_smooth)

    def get_smoothed_profile(self, x, y):
        y_smooth_non_zero = self.smooth_profile(x, y)
        return self.refine_smoothed_profile(x, y_smooth_non_zero)

    def refine_smoothed_profile(self, x_in, y_in):
        f = self.interpolate_smoothed_profile(x_in, y_in)
        x = numpy.linspace(x_in[0], 1.2, 1000)
        return x, f(x)

    def interpolate_smoothed_profile(self, x, y):
        x_ext = numpy.append(x, [x[-1] + 0.001, 1.2])
        y_ext = numpy.append(y, [0., 0.])
        return scipy.interpolate.interp1d(x_ext, y_ext, kind='linear')

    def get_fitted_profile(self, x, y):
        fit = scipy.optimize.curve_fit(stefanikova_ped_old, x, y,
                                       p0=[y[0], y[-1], 1, 0.02, y[0] / 1000])
        p = fit[0]
        return x, stefanikova_ped(p, x)

    def plot_profile(self, x_raw, y_raw, yerr_raw):
        x, y, yerr = self.get_valid_profile(x_raw, y_raw, yerr_raw)
        x_smooth, y_smooth = self.get_smoothed_profile(x, y)
        # x_fit, y_fit = get_fitted_profile(x, y_smooth)

        fig, ax = matplotlib.pyplot.subplots()
        ax.plot(x_raw, y_raw, '.')
        ax.plot(x, y, '.')
        ax.fill_between(x, y - 1 * yerr, y + 1 * yerr, color='gray', alpha=0.4)
        ax.fill_between(x, y - 3 * yerr, y + 3 * yerr, color='gray', alpha=0.3)
        ax.fill_between(x, y - 5 * yerr, y + 5 * yerr, color='gray', alpha=0.2)
        ax.plot(x_smooth, y_smooth)
        # ax.plot(x_fit, y_fit)
        matplotlib.pyplot.show()



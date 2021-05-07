import os
import h5py
import matplotlib
import numpy
import scipy
from statsmodels.nonparametric.smoothers_lowess import lowess

from stefanikova2016 import *


def get_ts_time_index(time, time_dataset):
    return (numpy.abs(time_dataset - int(time))).argmin()


def get_ts_dataset(field, thomson_directory, reconstruction_id):
    file = h5py.File(thomson_directory + '/' + field + '.' + str(reconstruction_id) + '.h5')
    return file[field][()]


def get_ts_profile(field, time_index, thomson_directory, reconstruction_id):
    return get_ts_dataset(field, thomson_directory, reconstruction_id)[time_index]


def filter_sol_outliers(x, y):
    return (numpy.fmin.accumulate(y) == y) | (x < 1.02)


def get_valid_profile_indices(x, y, yerr):
    return numpy.where(~numpy.isnan(y) & (x > 0.0) & (x < 1.2) &
                       (yerr < y) & filter_sol_outliers(x, y))


def get_valid_profile(x_in, y_in, yerr_in):
    valid_indices = get_valid_profile_indices(x_in, y_in, yerr_in)
    x = x_in[valid_indices]
    y = y_in[valid_indices]
    yerr = yerr_in[valid_indices]
    return x, y, yerr


def set_negatives_to_zero(x):
    return x*(x > 0)


def smooth_profile(x, y):
    smooth = lowess(y, x, is_sorted=True, frac=0.3, it=0)
    y_smooth = smooth[:, 1]
    return set_negatives_to_zero(y_smooth)


def get_smoothed_profile(x, y):
    y_smooth_non_zero = smooth_profile(x, y)
    return refine_smoothed_profile(x, y_smooth_non_zero)


def refine_smoothed_profile(x_in, y_in):
    f = interpolate_smoothed_profile(x_in, y_in)
    x = numpy.linspace(x_in[0], 1.2, 1000)
    return x, f(x)


def interpolate_smoothed_profile(x, y):
    x_ext = numpy.append(x, [x[-1]+0.001, 1.2])
    y_ext = numpy.append(y, [0., 0.])
    return scipy.interpolate.interp1d(x_ext, y_ext, kind='linear')


def get_fitted_profile(x, y):
    fit = scipy.optimize.curve_fit(stefanikova_ped_old, x, y,
                                   p0=[y[0], y[-1], 1, 0.02, y[0]/1000])
    p = fit[0]
    return x, stefanikova_ped(p, x)


def plot_profile(x_raw, y_raw, yerr_raw):
    x, y, yerr = get_valid_profile(x_raw, y_raw, yerr_raw)
    x_smooth, y_smooth = get_smoothed_profile(x, y)
    #x_fit, y_fit = get_fitted_profile(x, y_smooth)

    fig, ax = matplotlib.pyplot.subplots()
    ax.plot(x_raw, y_raw, '.')
    ax.plot(x, y, '.')
    ax.fill_between(x, y - 1 * yerr, y + 1 * yerr, color='gray', alpha=0.4)
    ax.fill_between(x, y - 3 * yerr, y + 3 * yerr, color='gray', alpha=0.3)
    ax.fill_between(x, y - 5 * yerr, y + 5 * yerr, color='gray', alpha=0.2)
    ax.plot(x_smooth, y_smooth)
    #ax.plot(x_fit, y_fit)
    matplotlib.pyplot.show()


def load_profiles(shot_number='17178', time='1097', reconstruction_id=2,
                  database_directory='input/cdb', thomson_subdir='THOMSON', efit_subdir='EFITXX',
                  ts_time_source='EFIT'):
    try:
        home_directory = os.environ['TAIGA_HOME']
    except KeyError:
        home_directory = '/home/matyi/work/taiga_local'  # '.'

    data_directory = home_directory + '/' + database_directory + '/' + str(shot_number)
    thomson_directory = data_directory + '/' + thomson_subdir
    efit_file = data_directory + '/' + efit_subdir + '/' + 'EFITXX.' + str(reconstruction_id) + '.h5'

    try:
        time_dataset = get_ts_dataset('TS_' + ts_time_source + '_time', thomson_directory, reconstruction_id=1)
        z_axis = get_ts_dataset('TS_z_axis', thomson_directory, reconstruction_id)
        print('Thomson scattering time and geometry files read successfully from: ' + thomson_directory)
    except OSError:
        raise OSError('Invalid Thomson scattering file structure! '
                      '\nPlease check the directory tree here:\n' + thomson_directory)
    except KeyError:
        raise KeyError('Invalid Thomson scattering data structure!\nExample for a correct structure:\nTe.1.h5\n\tTe')

    time_index = get_ts_time_index(time, time_dataset)

    try:
        temperature_profile = get_ts_profile('Te', time_index, thomson_directory, reconstruction_id)
        temperature_error_profile = get_ts_profile('Te_err', time_index, thomson_directory, reconstruction_id)
        density_profile = get_ts_profile('ne', time_index, thomson_directory, reconstruction_id)
        density_error_profile = get_ts_profile('ne_err', time_index, thomson_directory, reconstruction_id)
        normalised_poloidal_flux_profile = get_ts_profile('psi_n', time_index, thomson_directory, reconstruction_id=1)
        print('Thomson scattering profile files read successfully from: ' + thomson_directory)
    except OSError:
        raise OSError('Invalid Thomson scattering file structure! '
                      '\nPlease check the directory tree here:\n' + thomson_directory)
    except KeyError:
        raise KeyError('Invalid Thomson scattering data structure!\nExample for a correct structure:\nTe.1.h5\n\tTe')

    plot_profile(normalised_poloidal_flux_profile, temperature_profile, temperature_error_profile)
    plot_profile(normalised_poloidal_flux_profile, density_profile, density_error_profile)

    distance, density, temperature = load_profiles_mock()
    return distance, density, temperature


def load_profiles_mock():
    distance = [0, 0.1, 0.2]
    density = [1e19, 1.2e19, 1.5e19]
    temperature = [1000, 1500, 2500]
    return distance, density, temperature

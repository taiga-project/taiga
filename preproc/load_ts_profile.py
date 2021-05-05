import os
import numpy
import h5py
import matplotlib
import scipy
from statsmodels.nonparametric.smoothers_lowess import lowess

def get_ts_time_index(time, time_dataset):
    return (numpy.abs(time_dataset - int(time))).argmin()


def get_ts_dataset(field, thomson_directory, reconstruction_id):
    file = h5py.File(thomson_directory + '/' + field + '.' + str(reconstruction_id) + '.h5')
    return file[field][()]


def get_ts_profile(field, time_index, thomson_directory, reconstruction_id):
    return get_ts_dataset(field, thomson_directory, reconstruction_id)[time_index]


def stefanikova_core(r, a_height, a_width, a_exp):
    return a_height * (numpy.exp(-r**2 / a_width) ** a_exp)


def stefanikova_ped(r, b_height, b_sol, b_pos, b_width, b_slope):
    return (b_height - b_sol) / 2 * (mtanh((b_pos - r) / (2 * b_width), b_slope) + 1)


def mtanh(x, b_slope):
    return ((1 + b_slope * x) * scipy.exp(x) - scipy.exp(-x)) / (scipy.exp(x) + scipy.exp(-x))


def get_valid_profile_indices(x, y, yerr):
    return numpy.where(~numpy.isnan(x) & (yerr < y) & (x > 0.1) & (x < 1.1))


def get_valid_profile(x, y, yerr):
    valid_indices = get_valid_profile_indices(x, y, yerr)
    x = x[valid_indices]
    y = y[valid_indices]
    yerr = yerr[valid_indices]
    return x, y, yerr


def add_error_to_raw_profile(x, y, yerr):
    x2 = x * numpy.ones((100, len(x)))
    y2 = y + yerr * numpy.random.randn(100, len(y))
    return x2.flatten(), y2.flatten()


def get_fitted_profile(x, y, yerr):
    x2, y2 = add_error_to_raw_profile(x, y, yerr)
    #fit = scipy.optimize.curve_fit(stefanikova_ped, x, y, sigma=yerr, p0=[y2[0], y2[-1], 1, 0.02, y2[0]/1000])
    fit = scipy.optimize.curve_fit(stefanikova_core, x, y, sigma=3 * yerr, p0=[y2[0], 0.4, 3])
    p = fit[0]
    #return stefanikova_core(x, p[0], p[1], p[2])
    return lowess(y, x, is_sorted=True, frac=0.3, it=0)


def plot_profile(x, y, yerr):
    x, y, yerr = get_valid_profile(x, y, yerr)
    y_fit = get_fitted_profile(x, y, yerr)

    fig, ax = matplotlib.pyplot.subplots()
    ax.plot(x, y, '.')
    ax.fill_between(x, y - 3 * yerr, y + 3 * yerr, color='gray', alpha=0.2)
    ax.plot(x, y_fit)
    matplotlib.pyplot.show()


def load_profiles(shot_number='17178', time='1097', reconstruction_id=1,
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

    plot_profile(numpy.sqrt(normalised_poloidal_flux_profile), temperature_profile, temperature_error_profile)
    plot_profile(numpy.sqrt(normalised_poloidal_flux_profile), density_profile, density_error_profile)

    distance, density, temperature = load_profiles_mock()
    return distance, density, temperature


def load_profiles_mock():
    distance = [0, 0.1, 0.2]
    density = [1e19, 1.2e19, 1.5e19]
    temperature = [1000, 1500, 2500]
    return distance, density, temperature

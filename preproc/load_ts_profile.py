import os
import numpy
import h5py


def get_ts_time_index(time, time_dataset):
    return (numpy.abs(time_dataset - int(time))).argmin()


def get_ts_dataset(field, thomson_directory, reconstruction_id):
    file = h5py.File(thomson_directory + '/' + field + '.' + str(reconstruction_id) + '.h5')
    return file[field][()]


def get_ts_profile(field, time_index, thomson_directory, reconstruction_id):
    return get_ts_dataset(field, thomson_directory, reconstruction_id)[time_index]


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
        time_dataset = get_ts_dataset('TS_' + ts_time_source + '_time', thomson_directory, reconstruction_id)
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
        print('Thomson scattering profile files read successfully from: ' + thomson_directory)
    except OSError:
        raise OSError('Invalid Thomson scattering file structure! '
                      '\nPlease check the directory tree here:\n' + thomson_directory)
    except KeyError:
        raise KeyError('Invalid Thomson scattering data structure!\nExample for a correct structure:\nTe.1.h5\n\tTe')

    distance, density, temperature = load_profiles_mock()
    return distance, density, temperature


def load_profiles_mock():
    distance = [0, 0.1, 0.2]
    density = [1e19, 1.2e19, 1.5e19]
    temperature = [1000, 1500, 2500]
    return distance, density, temperature

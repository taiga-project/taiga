import matplotlib.pyplot

from manager import RenateODManager
from beamlet import set_beamlet
from efit import EFITManager
from utils import *


def export_beamlet_profile(export_root=get_home_directory() + '/input/ionProf/',
                           shot_number='17178', time='1097', species='Li', energy='80'):
    z = 0
    tor = 0
    beamlet_geometry = set_beamlet(z, tor)
    r = RenateODManager(beamlet_geometry, shot_number, time, species, energy)
    radial_coordinate, relative_attenuation = r.get_attenuation_profile()
        
    export_directory = export_root + '/'+shot_number+'_'+time
    try:
        os.mkdir(export_directory)
        print('Create directory and write data to ' + export_directory)
    except FileExistsError:
        print('Write data to ' + export_directory)
    else:
        pass
    print('Save RENATE-OD ionisation profile to: ' + export_directory)
    radial_coordinate.fillna(0).to_csv(export_directory+'/rad.dat', index=False, header=False)
    relative_attenuation.fillna(0).to_csv(export_directory+'/ionyeald.dat', index=False, header=False)
    plot_attenuation_profile(shot_number, time, radial_coordinate, relative_attenuation, export_directory)


def get_lcfs_radial_coordinate(shot_number, time, efit_reconstruction_id=1, database_directory='input/cdb', efit_subdir='EFITXX'):
    data_directory = get_home_directory() + '/' + database_directory + '/' + str(shot_number)
    efit_file = data_directory + '/' + efit_subdir + '/' + 'EFITXX.' + str(efit_reconstruction_id) + '.h5'
    efit = EFITManager(efit_file, time)
    return efit.get_time_sliced_data('output/separatrixGeometry/rmidplaneOut')


def plot_attenuation_profile(shot_number, time, radial_coordinate, relative_attenuation, export_directory='.'):
    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(5, 2)
    ax.plot(radial_coordinate, relative_attenuation, '-', linewidth=2)
    matplotlib.pyplot.minorticks_on()
    matplotlib.pyplot.grid(which='major')
    matplotlib.pyplot.xlabel('$R$ [m]', labelpad=-10.5, loc='right')
    matplotlib.pyplot.ylabel('neutral beam attenuation')
    matplotlib.pyplot.title('COMPASS #' + shot_number + ' (' + time + ' ms)')
    R_LCFS = get_lcfs_radial_coordinate(shot_number, time)
    matplotlib.pyplot.axvline(R_LCFS, c='red', ls='--')
    matplotlib.pyplot.text(R_LCFS+0.005, 0.45, 'LCFS', c='red', fontsize=12)
    matplotlib.pyplot.savefig(export_directory+'/attenuation.pdf')
    matplotlib.pyplot.savefig(export_directory+'/attenuation.svg')
    print('Attenuation profile saved as: '+export_directory+'/attenuation.pdf')
    matplotlib.pyplot.show()


if __name__ == "__main__":
    a_shot_number = '17178'
    a_time = '1097'
    export_beamlet_profile(shot_number=a_shot_number, time=a_time)

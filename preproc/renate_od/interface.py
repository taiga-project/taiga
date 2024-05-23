import matplotlib.pyplot

from renate_od.manager import RenateODManager
from renate_od.beamlet import set_beamlet
from renate_od.efit import EFITManager
from renate_od.utils import *


def get_lcfs_radial_coordinate(time, shot_number, efit_reconstruction_id=1, database_directory='input/cdb', efit_subdir='EFITXX'):
    data_directory = os.path.join(get_home_directory(), database_directory, str(shot_number))
    efit_file = os.path.join(data_directory, efit_subdir, 'EFITXX.' + str(efit_reconstruction_id) + '.h5')
    efit = EFITManager(efit_file, time)
    return efit.get_time_sliced_data('output/separatrixGeometry/rmidplaneOut')


class SetProfiles:
    def __init__(self, shot_number, time, species, energy, is_forced_saved=False):
        self.shot_number = str(int(shot_number))
        self.time = str(int(time))
        self.species = species
        self.energy = str(int(energy))
        self.export_root = os.path.join(get_home_directory(), 'input', 'ionProf')
        self.is_saved = is_forced_saved
        self.radial_coordinate = None
        self.relative_attenuation = None
        self.export_directory = None

        self.init_export_directory()
        if self.is_saved:
            self.get_attenuation()
            self.export_beamlet_profile()
            self.plot_attenuation_profile()

    def set_export_directory(self):
        self.export_directory = os.path.join(self.export_root, self.shot_number + '_' + self.time, self.species, self.energy)

    def init_export_directory(self):
        self.set_export_directory()
        try:
            os.makedirs(self.export_directory)
            self.is_saved = True
            print('Create attenuation profile directory: ' + self.export_directory)
        except OSError:
            print('Existing attenuation profile directory: ' + self.export_directory)

    def get_attenuation(self):
        z = 0
        tor = 0
        beamlet_geometry = set_beamlet(z, tor)
        r = RenateODManager(beamlet_geometry, self.shot_number, self.time, self.species, self.energy)
        self.radial_coordinate, self.relative_attenuation = r.get_attenuation_profile()

    def export_beamlet_profile(self):
        print('Save RENATE-OD ionisation profile to: ' + self.export_directory)
        self.radial_coordinate.fillna(0).to_csv(os.path.join(self.export_directory, 'rad.dat'), index=False, header=False)
        self.relative_attenuation.fillna(0).to_csv(os.path.join(self.export_directory, 'ionyeald.dat'), index=False, header=False)

    def plot_attenuation_profile(self):
        fig, ax = matplotlib.pyplot.subplots()
        fig.set_size_inches(5, 2)
        ax.plot(self.radial_coordinate, self.relative_attenuation, '-', linewidth=2)
        matplotlib.pyplot.minorticks_on()
        matplotlib.pyplot.grid(which='major')
        matplotlib.pyplot.xlabel('$R$ [m]', labelpad=-10.5, loc='right')
        matplotlib.pyplot.ylabel('neutral beam attenuation')
        matplotlib.pyplot.title('COMPASS #' + self.shot_number + ' (' + self.time + ' ms)')
        R_LCFS = get_lcfs_radial_coordinate(time=self.time, shot_number=self.shot_number)
        matplotlib.pyplot.axvline(R_LCFS, c='red', ls='--')
        matplotlib.pyplot.text(R_LCFS+0.005, 0.45, 'LCFS', c='red', fontsize=12)
        matplotlib.pyplot.savefig(os.path.join(self.export_directory, 'attenuation.pdf'))
        matplotlib.pyplot.savefig(os.path.join(self.export_directory, 'attenuation.svg'))
        print('Attenuation profile saved to '+ self.export_directory)
        matplotlib.pyplot.show()


if __name__ == "__main__":
    SetProfiles(17178, 1097, 'Li', 80)

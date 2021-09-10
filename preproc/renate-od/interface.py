import matplotlib.pyplot
from manager import RenateODManager



def set_beamlet(z, tor):
    beamlet_geometry = BeamletGeometry()
    beamlet_geometry.rad = numpy.linspace(0.78, 0.35, 200)
    beamlet_geometry.set_with_value(z, 'z', 'rad')
    beamlet_geometry.set_with_value(tor, 'tor', 'rad')
    return beamlet_geometry


def export_beamlet_profile(export_root=get_home_directory() + '/input/ionProf/',
                           shot_number='17178', time='1097', species='Li', energy='80'):
    z = 0
    tor = 0
    beamlet_geometry = set_beamlet(z, tor)
    RenateODManager r(beamlet_geometry, shot_number, time, species, energy)
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
    radial_coordinate.to_csv(export_directory+'/rad.dat', index=False, header=False)
    relative_attenuation.to_csv(export_directory+'/ionyeald.dat', index=False, header=False)
    plot_attenuation_profile(radial_coordinate, relative_attenuation, export_directory, shot_number, time)


def plot_attenuation_profile(radial_coordinate, relative_attenuation, export_directory='.', shot_number, time):
    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(5, 2)
    ax.plot(radial_coordinate, relative_attenuation, '-', linewidth=2)
    matplotlib.pyplot.minorticks_on()
    matplotlib.pyplot.grid(which='both')
    matplotlib.pyplot.xlabel('$R$ [m]', labelpad=-10.5, loc='right')
    matplotlib.pyplot.ylabel('neutral beam attenuation')
    matplotlib.pyplot.title('COMPASS #' + shot_number + ' (' + time + ' ms)')
    R_LCFS = 0.7143 #hack
    matplotlib.pyplot.axvline(R_LCFS, c='red', ls='--')
    matplotlib.pyplot.text(R_LCFS+0.005, 0.45, 'LCFS', c='red', fontsize=12)
    matplotlib.pyplot.savefig(export_directory+'/attenuation.pdf')
    matplotlib.pyplot.savefig(export_directory+'/attenuation.svg')
    print('Attenuation profile saved as: '+export_directory+'/attenuation.pdf')
    matplotlib.pyplot.show()


def mock_beam(shot_number='17178', time='1097', diameter=5e-3, z_length=3, tor_length=3):
    fig, ax = matplotlib.pyplot.subplots()
    species = 'Li'
    energy = '80'
    for z in numpy.linspace(-diameter/2, diameter/2, z_length):
        if z_length == 1:
            z = 0
        for tor in numpy.linspace(-diameter/2, diameter/2, tor_length):
            if tor_length == 1:
                tor = 0
            beamlet_geometry = set_beamlet(z, tor)
            radial_coordinate, ionisation_degree,  attenuation_degree = \
                calculate_beamlet(beamlet_geometry, shot_number, time, species, energy)
            ax.plot(radial_coordinate, attenuation_degree, '-')

    ax.set_xlabel('R [m]')
    ax.set_ylabel('normalised linear density attenuation')
    ax.set_title('COMPASS #'+shot_number+' ('+time+' ms) ')
    matplotlib.pyplot.xlim(0.6, 0.75)
    matplotlib.pyplot.savefig('attenuation_'+shot_number+'_'+time+'.svg')
    matplotlib.pyplot.show()


if __name__ == "__main__":
    shot_number = '21762'
    time = '1120'
    shot_number = '17679'
    time = '1120'
    shot_number = '17178'
    time = '1097'
    export_beamlet_profile(shot_number=shot_number, time=time)
    #mock_beam()


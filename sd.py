import os
from time import localtime, strftime

from preproc.init_from_efit import CDBManager
from preproc.renate_od.interface import SetProfiles
from preproc.renate_od.utils import get_home_directory


class SD:
    def __init__(self, shot_number, time, species, energy):
        self.shot_number = str(int(shot_number))
        self.time = str(int(time))
        self.species = species
        self.energy = str(int(energy))
        self.runnumber = strftime("%Y%m%d%H%M%S", localtime())
        self.parameter_file = os.path.join(get_home_directory(),
                                           "parameter_" + self.runnumber + ".sh")

        self.write_parameter_file()
        self.preproc()
        self.run_core()
        self.delete_parameter_file()

    def preproc(self):
        CDBManager(self.shot_number, self.time)
        SetProfiles(self.shot_number, self.time, self.species, self.energy)

    def run_core(self):
        os.chdir(get_home_directory())
        command = "./taiga_renate.exe -p=" + self.parameter_file + " -r="+self.runnumber
        print(command)
        os.system(command)
        os.system("pwd")

    def write_parameter_file(self):
        f = open(self.parameter_file, "w")
        f.write("##!/usr/bin/bash\n")
        f.write("shotnumber='" + self.shot_number + "'\n")
        f.write("time='" + self.time + "'\n")
        f.write("beammatter='" + self.species + "'\n")
        f.write("energy=" + self.energy + " #keV\n")
        f.write("toroidal_deflection=0 #degree; + 200 V deflection\n")
        f.write("diameter=5 #mm\n")
        f.write("particles=1000\n")
        f.write("detector='0.6846,0.253,0.0,38,0'\n")
        f.write("electric_field_module=0\n")
        f.write("detector_mask='final'\n")
        f.write("solver='rk'\n")
        f.write("secondary_ionisation=1")
        f.close()

    def delete_parameter_file(self):
        os.remove(self.parameter_file)


if __name__ == "__main__":
    a_shot_number = 17178
    a_time = 1097
    a_species = 'Na'
    an_energy = 80
    SD(a_shot_number, a_time, a_species, an_energy)

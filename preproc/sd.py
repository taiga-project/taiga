from init_from_efit import CDBManager
from renate_od.interface import SetProfiles


if __name__ == "__main__":
    a_shot_number = 17178
    a_time = 1097
    a_species = 'Na'
    an_energy = 80
    CDBManager(a_shot_number, a_time)
    SetProfiles(a_shot_number, a_time, a_species, an_energy)

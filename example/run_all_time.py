import taiga.preproc
import taiga.preproc.renate_od
#from taiga import EFITDataReader, CDBReader, ParseToTaiga
#from taiga import export_beamlet_profile
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             os.pardir, 'preproc')))

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             os.pardir, 'preproc', 'renate_od')))
from init_from_efit import EFITDataReader, CDBReader, ParseToTaiga
from interface import export_beamlet_profile

if __name__ == "__main__":
    a_shot_number = 17178
    a_time = 1097
    e = EFITDataReader(a_shot_number)
    time_dataset = e.get_data('time')
    export_beamlet_profile(shot_number=str(a_shot_number), time=str(a_time))
    cr = CDBReader(17178, 1097)
    ParseToTaiga(cr)

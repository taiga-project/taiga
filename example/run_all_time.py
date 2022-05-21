import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             os.pardir, 'preproc')))

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             os.pardir, 'preproc', 'renate_od')))
from init_from_efit import EFITDataReader, CDBReader, ParseToTaiga
from interface import export_beamlet_profile


def run_all_time(shot_number):
    e = EFITDataReader(shot_number)
    times = e.get_data('time')
    time_slices = (times * 1000).astype(int)
    for time in time_slices:
        export_beamlet_profile(shot_number=str(shot_number), time=str(time))
        cr = CDBReader(shot_number, time)
        ParseToTaiga(cr)


if __name__ == "__main__":
    run_all_time(17178)

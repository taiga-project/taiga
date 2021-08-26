#ifndef INIT_DETECTOR_CUH
#define INIT_DETECTOR_CUH

void init_detector(DetectorProp* shared_detector, DetectorProp *device_detector, ShotProp shot);
void set_detector_geometry(ShotProp shot, TaigaCommons *host_common, TaigaCommons *shared_common, DetectorProp *shared_detector);

#endif //INIT_DETECTOR_CUH

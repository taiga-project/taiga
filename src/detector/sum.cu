__global__ void detector_sum(TaigaGlobals *global, TaigaCommons *common, DetectorProp *detector){
    for (int i = 0; i < detector->number_of_detector_cells; ++i){
        detector->counter[i] = 0;
    }
    for (int idx = 0; idx < global->particle_number; ++idx){
        if ((global->detcellid[idx] > 0) 
        & (global->detcellid[idx] < detector->number_of_detector_cells+1) ){
            detector->counter[global->detcellid[idx]-1]+= global->intensity[idx];
        }
    }
}

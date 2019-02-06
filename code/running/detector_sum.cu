__global__ void detector_sum(int *det, int *detcellid, int N, int Ndg){
    for (int i=0; i<Ndg; i++){
        det[i] = 0;
    }
    for (int idx=0; idx<N; idx++){
        if (detcellid[idx]>0){
            det[detcellid[idx]-1]++;
        }
    }
}

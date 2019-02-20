__global__ void undetected(int *detcellid, int N, double *service_var){
    int undetected_counter = 0;
    for (int idx=0; idx<N; idx++){
        if (detcellid[idx]<0){
            undetected_counter++;
        }
    }
    service_var[1] = (double)undetected_counter/N;
}

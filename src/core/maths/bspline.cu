#include "bspline.cuh"

// https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-coef.html
__device__ void bspline(double *B, double x, int k, int index, double *t){
    int i, j, d;
    for (i=0; i<k; ++i){
        B[i] = 0.0;
    }
    B[k] = 1.0;
    if (index<k) {
        B[0]=1.0;
        return;
    }
    for (d=1; d<=k; ++d){
        if (t[index + 1] == t[index - d + 1]){
            B[k - d] = 0;// B[k - d + 1];
        }else {
            B[k - d] = (t[index + 1] - x) / (t[index + 1] - t[index - d + 1]) * B[k - d + 1];
        }
        for (i=index-d+1; i<index; ++i) {
            j = i - index + k;
            if ((t[i + d] == t[i]) || (t[i + d + 1] == t[i + 1])) {
                B[j] = 0.0;
            } else {
                B[j] = (x - t[i]) / (t[i + d] - t[i]) * B[j] + \
                       (t[i + d + 1] - x) / (t[i + d + 1] - t[i + 1]) * B[j + 1];
            }
        }
        if (t[index + d] == t[index]){
            ;//B[k]= 1.0;
        }else{
            B[k] = (x - t[index]) / (t[index + d] - t[index]) * B[k];
        }
    }
}

__device__ double B(double x, int k, int i, double *t){
    if (k == 0){
        return ((t[i] <= x)&(x < t[i + 1])) ? 1.0 : 0.0;
    }
    double c1, c2;
    if (t[i + k] == t[i]){
        c1 = 0.0;
    }else{
        c1 = (x - t[i]) / (t[i + k] - t[i]) * B(x, k - 1, i, t);
    }
    if (t[i + k + 1] == t[i + 1]){
        c2 = 0.0;
    }else{
        c2 = (t[i + k + 1] - x) / (t[i + k + 1] - t[i + 1]) * B(x, k - 1, i + 1, t);
    }
    return c1 + c2;
}

__device__ double get_rad_from_poloidal(double R, double l_br, double l_bt, double l_r, double l_t){
    double tan_phi = l_t/l_r;
    double cos_phi = l_r/R;
    return l_br * cos_phi - l_bt * tan_phi*cos_phi;
}

__device__ double get_tor_from_poloidal(double R, double l_br, double l_bt, double l_r, double l_t){
    double tan_phi = l_t/l_r;
    double cos_phi = l_r/R;
    return l_br * tan_phi*cos_phi + l_bt * cos_phi;
}

__device__ inline double get_major_radius(double l_r, double l_z){
    return hypot(l_r, l_z);
}

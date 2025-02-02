#include "init/thomson.cuh"
#include "utils/cuda.cuh"

void import_thomson_profiles(ShotProp shot, TaigaCommons *c) {
    c->ts_length = read_vector(&c->ts_psi, "input/tsProf", shot.name, "flux.prof");
    read_vector(&c->ts_density, "input/tsProf", shot.name, "density.prof");
    read_vector(&c->ts_temperature, "input/tsProf", shot.name, "temperature.prof");
}

void set_thomson_profiles(ShotProp shot, TaigaCommons *host_common, TaigaCommons *shared_common) {
    double *shared_ts_psi, *shared_ts_density, *shared_ts_temperature;
    import_thomson_profiles(shot, host_common);
    size_t size_ts = host_common->ts_length * sizeof(double);

    CHECK_ERROR(cudaMalloc((void **) &shared_ts_psi, size_ts));
    CHECK_ERROR(cudaMalloc((void **) &shared_ts_density, size_ts));
    CHECK_ERROR(cudaMalloc((void **) &shared_ts_temperature, size_ts));

    CHECK_ERROR(cudaMemcpy(shared_ts_psi,  host_common->ts_psi, size_ts, cudaMemcpyHostToDevice));
    CHECK_ERROR(cudaMemcpy(shared_ts_density,  host_common->ts_density, size_ts, cudaMemcpyHostToDevice));
    CHECK_ERROR(cudaMemcpy(shared_ts_temperature,  host_common->ts_temperature, size_ts, cudaMemcpyHostToDevice));

    shared_common->ts_length = host_common->ts_length;
    shared_common->ts_psi = shared_ts_psi;
    shared_common->ts_density = shared_ts_density;
    shared_common->ts_temperature = shared_ts_temperature;
}

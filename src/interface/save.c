#include "save.h"
#include "utils/prop.h"
#include "utils/basic_functions.h"
#include "dataio/data_export.h"

void save_trajectories(TaigaGlobals *host_global, RunProp run){
    export_data(host_global->rad,  run.particle_number, run.folder_out, run.runnumber, "t_rad.dat");
    export_data(host_global->z,    run.particle_number, run.folder_out, run.runnumber, "t_z.dat");
    export_data(host_global->tor,  run.particle_number, run.folder_out, run.runnumber, "t_tor.dat");
    export_data(host_global->vrad, run.particle_number, run.folder_out, run.runnumber, "t_vrad.dat");
    export_data(host_global->vz,   run.particle_number, run.folder_out, run.runnumber, "t_vz.dat");
    export_data(host_global->vtor, run.particle_number, run.folder_out, run.runnumber, "t_vtor.dat");
    export_data(host_global->intensity, run.particle_number, run.folder_out, run.runnumber, "t_intensity.dat");
    export_data(host_global->time_of_flight, run.particle_number, run.folder_out, run.runnumber, "t_tof.dat");
}

void save_endpoints(TaigaGlobals *host_global, RunProp run){
    export_data(host_global->rad,  run.particle_number, run.folder_out, run.runnumber, "rad.dat");
    export_data(host_global->z,    run.particle_number, run.folder_out, run.runnumber, "z.dat");
    export_data(host_global->tor,  run.particle_number, run.folder_out, run.runnumber, "tor.dat");
    export_data(host_global->vrad, run.particle_number, run.folder_out, run.runnumber, "vrad.dat");
    export_data(host_global->vz,   run.particle_number, run.folder_out, run.runnumber, "vz.dat");
    export_data(host_global->vtor, run.particle_number, run.folder_out, run.runnumber, "vtor.dat");
    export_data(host_global->intensity, run.particle_number, run.folder_out, run.runnumber, "intensity.dat");
    export_data(host_global->time_of_flight, run.particle_number, run.folder_out, run.runnumber, "tof.dat");

    export_table(run.folder_out, run.runnumber, "coords.dat", run.particle_number,
                 host_global->rad, "R [m]",      host_global->z, "Z [m]",      host_global->tor, "T [m]",
                 host_global->vrad, "v_R [m/s]", host_global->vz, "v_Z [m/s]", host_global->vtor, "v_T [m/s]");
}

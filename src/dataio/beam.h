#ifndef BEAM_H
#define BEAM_H

void load_beam(TaigaGlobals *g, BeamProp *beam, ShotProp *shot, RunProp *run);
void init_ion_profile(char* shotname, BeamProfile* prof);
inline double calculate_speed(double beam_energy, double beam_mass){
    return sqrt(2.0 * beam_energy*1000.0*ELEMENTARY_CHARGE/ beam_mass/ AMU);
}

#endif // BEAM_H

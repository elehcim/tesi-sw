#include <string>
//Parameters
extern double mu, ecc, DT, t0, tf, abs_tol, rel_tol, L, d1, d2;
extern int n_frames, n_cores, n_iterations;
extern std::string field_type;
extern std::string custom_file_name;

//Vectors dimensions
extern int nx, ny, nvx, nvy, ne, n1, n2, n_tot, vis_count;

//Vectors extremes
extern double x_min, x_max, y_min, y_max, vx_min, vx_max, vy_min, vy_max, e_min, e_max, dx, dy, dvx, dvy, de;

//Input values
extern double X_0, Y_0, VX_0, VY_0, E_0;

//Others
extern int flags[5], flags_t[5];
extern bool tracers_flag, distance_flag, matlab_flag, name_flag;
extern double d1, d2;

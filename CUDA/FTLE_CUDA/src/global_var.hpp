//Parameters
extern double mu;
extern double ecc;
extern double DT;
extern int n_frames;

//Vectors dimensions
extern int nx;
extern int ny;
extern int nvx;
extern int nvy;
extern int ne;
extern int n1;
extern int n2;
extern int n_tot;
extern int vis_count;

//Vectors extremes
extern double x_min;
extern double x_max;
extern double y_min;
extern double y_max;
extern double vx_min;
extern double vx_max;
extern double vy_min;
extern double vy_max;
extern double e_min;
extern double e_max;
extern double dx;
extern double dy;
extern double dvx;
extern double dvy;
extern double de;

//Input values
extern double X_0;
extern double Y_0;
extern double VX_0;
extern double VY_0;
extern double E_0;

//Others
extern int flags[5], flags_t[5];
extern bool tracers_flag;

//Tracers' vectors dimension
extern int nxt, nyt, nvxt, nvyt, net, n1t, n2t, n_tot_t;

//Tracers' vectors extremes
extern double x_mint, x_maxt, y_mint, y_maxt, vx_mint, vx_maxt, vy_mint, vy_maxt, e_mint, e_maxt, dxt, dyt, dvxt, dvyt, det;

//Tracers' input values
extern double X_0t, Y_0t, VX_0t, VY_0t, E_0t;

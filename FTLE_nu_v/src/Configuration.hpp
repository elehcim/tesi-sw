
double steps;
double pi=3.141592653589793;

double mu=9.537e-4;
double ecc=0.04839;
double T=5;
double R=0.19;
double nu_min=0;
double nu_max=10;
int n_nu=501;
double v_min=0.45;
double v_max=0.6;
int n_v=501;
// dx=0.01 is ~1/20 Earth-Sun distance
// dvx=0.1 is ~1/20 Earth orbital speed
double dx=0.001;
double dy=0.001;
double dvx=0.001;
double dvy=0.001;
double abs_tol=1e-12;
double rel_tol=1e-12;
double in_step=1e-3;
double t0;
double done_perc;
double v_E=29.783; //km/s


double1d nu_0(n_nu);
double1d v_0(n_v);
double1d X1(4);
double1d X2(4);
double1d X3(4);
double1d X4(4);
double1d X5(4);
double1d X6(4);
double1d X7(4);
double1d X8(4);

FILE *ftle_stream;

double a_jup=778412027; //km
double GM_jup=126711995; //km^3/s^2
double GM_sun=132712439935; //km^3/s^2
double GM=GM_jup+GM_sun;

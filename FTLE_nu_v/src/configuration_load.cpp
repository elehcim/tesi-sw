#include <boost/program_options.hpp>
#include <fstream>
#include <stdio.h>
#include <string>
#include <math.h>

#ifdef _WIN32
#define M_PI 3.141592653589793238462
#endif // WIN32

namespace po = boost::program_options;
double mu, ecc, x_min, x_max, y_min, y_max, vx_min, vx_max, vy_min, vy_max, e_min, e_max, X_0, Y_0, VX_0, VY_0, E_0, DT, t0, tf,
        dx, dy, dvx, dvy, de, d1, d2, L, abs_tol, rel_tol;
int nx, ny, nvx, nvy, ne, n1, n2, n_tot, n_frames, n_cores, n_iterations;
int flags[5], flags_t[5];
int vis_count=0;        //Intero in barba a Calvino!
int in_count=0;
bool distance_flag, matlab_flag, name_flag=0;
std::string field_type;
std::string custom_file_name;
int configuration_load (std::string config_file)
{
    po::options_description options("Allowed options on configuration file");
    options.add_options()
        ("parameters.mu", po::value<double>(),"mu")
        ("parameters.ecc", po::value<double>(),"ecc")
        ("parameters.t0", po::value<double>()->default_value(0.0),"t0")
        ("parameters.tf", po::value<double>()->default_value(2*M_PI),"tf")
        ("parameters.DT", po::value<double>()->default_value(10.0),"DT")
        ("parameters.n_frames", po::value<int>()->default_value(1),"set n_frames")
        ("parameters.abs_tol", po::value<double>()->default_value(1e-6), "Integrator absolute tolerance" )
        ("parameters.rel_tol", po::value<double>()->default_value(1e-6), "Integrator relative tolerance" )
        ("parameters.n_cores", po::value<int>()->default_value(1), "Max number of threads")
        ("parameters.file_name", po::value<std::string>(),"Output file name")
        ("vis.var.n_nu", po::value<int>(), "set n_nu")
        ("vis.var.n_v", po::value<int>(), "set n_v")
        ("vis.var.nu_min", po::value<double>()->default_value(0), "set nu_min")
        ("vis.var.nu_max", po::value<double>(), "set nu_max")
        ("vis.var.nu_min", po::value<double>()->default_value(0), "set v_min")
        ("vis.var.nu_max", po::value<double>(), "set v_max")
        ("in.var.dx", po::value<double>(), "set dx")
        ("in.var.dy", po::value<double>(), "set dy")
        ("in.var.dvx", po::value<double>(), "set dvx")
        ("in.var.dvy", po::value<double>(), "set dvy")
        ("distance.flag", po::value<bool>()->default_value(0), "distance flag")
        ("distance.L", po::value<double>(), "distance between primaries")
        ("distance.d1", po::value<double>(), "distance from first primary")
        ("distance.d2", po::value<double>(), "distance from second primary")
        ;

    po::variables_map vm;

    std::ifstream ifs(config_file.c_str());
    store(parse_config_file(ifs, options), vm);
    notify(vm);

    // Load parameters
    if (vm.count("parameters.mu")) {mu=vm["parameters.mu"].as<double>();}
    else    {std::cout<<"Mass parameter missing!\n";
            return 1;}
    if (vm.count("parameters.ecc")) {ecc=vm["parameters.ecc"].as<double>();}
    else    {std::cout<<"Eccentricity missing!\n";
            return 1;}
    if (vm.count("parameters.t0")){t0=vm["parameters.t0"].as<double>();}
    if (vm.count("parameters.tf")){tf=vm["parameters.tf"].as<double>();}
    if (vm.count("parameters.DT")) {DT=vm["parameters.DT"].as<double>();}
    else if (field_type=="FTLE"){std::cout<<"Integration interval missing!\n";
                                return 1;}
    if (vm.count("parameters.n_frames"))
        {
            n_frames=vm["parameters.n_frames"].as<int>();
            if (n_frames<=0){std::cout<<"Number of frames cannot be 0 or less\n";
                            return 1;}
        }
    else    {std::cout<<"Number of frames missing!\n";
            return 1;}
    if (vm.count("parameters.abs_tol")){abs_tol=vm["parameters.abs_tol"].as<double>();}
    if (abs_tol<1e-20){std::cout<<"Please provide an absolute tolerance higher than 1e-20\n";
                        return 1;}
    if (vm.count("parameters.rel_tol")){rel_tol=vm["parameters.rel_tol"].as<double>();}
    if (rel_tol<1e-20){std::cout<<"Please provide a relative tolerance higher than 1e-20\n";
                        return 1;}
    if (vm.count("parameters.n_cores")){n_cores=vm["parameters.n_cores"].as<int>();}
    if (vm.count("parameters.file_name")){custom_file_name=vm["parameters.file_name"].as<std::string>();
                                        name_flag=1;}

    printf("mu=%.4f ecc=%.4f DT=%.2f n_frames=%i\n",mu,ecc,DT,n_frames);

    // Load vectors dimensions and input values

    if (vm.count("vis_var.n_nu")) {n_nu=vm["vis_var.n_nu"].as<double>();}
    if (vm.count("vis_var.n_v")) {n_nu=vm["vis_var.n_v"].as<double>();}
    if (vm.count("vis_var.n_nu")) {n_nu=vm["vis_var.n_nu"].as<double>();}


        // check number of variables
    if (field_type=="FTLE")
    {
    if (vis_count<2)
    {
        std::cout<<"Too few visualization variables!\n";
        return 1;
    }
    else if (vis_count+in_count!=4)
        {
            std::cout<<"Not enough fixed variables. Sum of the number of visualization variables and fixed variables must be 4\n";
            return 1;
        }

    if (in_count>2)
    {
        std::cout<<"Too many fixed variables!\n";
        return 1;
    }
    }
    if (field_type=="FILE")
    {
        if (vm.count("in.var.x") && vm.count("in.var.y")){std::cout<<"Please specify only one spatial fixed variable\n";
                                                            return 1;}
        if (!(vm.count("in.var.x")) && !(vm.count("in.var.y"))){std::cout<<"Please specify one spatial fixed variable\n";
                                                            return 1;}
    }

        //Create dimension of missing vector
    for (int i=0;i<5;i++)
    {
        if (flags[i]==0)
        {
            switch (i)
            {
                case 0: //not supported
                    nx=ny*nvx*nvy*ne;
                    n_tot=nx;
                    break;
                case 1: //not supported
                    ny=nx*nvx*nvy*ne;
                    n_tot=ny;
                    break;
                case 2:
                    nvx=nx*ny*nvy*ne;
                    n_tot=nvx;
                    break;
                case 3:
                    nvy=nx*ny*nvx*ne;
                    n_tot=nvy;
                    break;
                case 4:
                    ne=nx*ny*nvx*nvy;
                    n_tot=ne;
                    break;
            }
        }
    }

/* Distance */
if (vm.count("distance.flag")){distance_flag=vm["distance.flag"].as<bool>();}

if (distance_flag)
{
    if (!(vm.count("distance.L"))){std::cout<<"Please provide distance between primaries\n";
                                    return 1;}
    if (!(vm.count("distance.d1"))){std::cout<<"Please provide forbidden distance from 1st primary\n";
                                    return 1;}
    if (!(vm.count("distance.d1"))){std::cout<<"Please provide forbidden distance from 2nd primary\n";
                                    return 1;}
}

if (vm.count("distance.L")){L=vm["distance.L"].as<double>();}
if (vm.count("distance.d1")){d1=(vm["distance.d1"].as<double>())/L;}
if (vm.count("distance.d2")){d2=(vm["distance.d2"].as<double>())/L;}
/* Matlab */
if (vm.count("matlab.flag")){matlab_flag=vm["matlab.flag"].as<bool>();}
    return 0;
}

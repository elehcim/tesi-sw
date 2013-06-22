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
int vis_count=0;        //Intero in barba a Calvino! (non ho potuto resistere :-) )
int in_count=0;
bool distance_flag, matlab_flag, name_flag=0;
std::string field_type;
std::string custom_file_name;
int configuration_load ()
{
    po::options_description options("Allowed options");
    options.add_options()
        ("parameters.mu", po::value<double>(),"mu")
        ("parameters.ecc", po::value<double>(),"ecc")
        ("parameters.type", po::value<std::string>()->default_value("FTLE"),"Type of field")
        ("parameters.t0", po::value<double>()->default_value(0.0),"t0")
        ("parameters.tf", po::value<double>()->default_value(2*M_PI),"tf")
        ("parameters.DT", po::value<double>()->default_value(10.0),"DT")
        ("parameters.n_iterations", po::value<int>()->default_value(1), "Number of intersections")
        ("parameters.n_frames", po::value<int>()->default_value(1),"set n_frames")
        ("parameters.abs_tol", po::value<double>()->default_value(1e-6), "Integrator absolute tolerance" )
        ("parameters.rel_tol", po::value<double>()->default_value(1e-6), "Integrator relative tolerance" )
        ("parameters.n_cores", po::value<int>(), "Max number of cores")
        ("parameters.file_name", po::value<std::string>(),"Output file name")
        ("vis.var.nx", po::value<int>(), "set nx")
        ("vis.var.ny", po::value<int>(), "set ny")
        ("vis.var.nvx", po::value<int>(), "set nvx")
        ("vis.var.nvy", po::value<int>(), "set nvy")
        ("vis.var.ne", po::value<int>(), "set ne")
        ("vis.var.x_min", po::value<double>(), "set x_min")
        ("vis.var.x_max", po::value<double>(), "set x_max")
        ("vis.var.y_min", po::value<double>(), "set y_min")
        ("vis.var.y_max", po::value<double>(), "set y_max")
        ("vis.var.vx_min", po::value<double>(), "set vx_min")
        ("vis.var.vx_max", po::value<double>(), "set vx_max")
        ("vis.var.vy_min", po::value<double>(), "set vy_min")
        ("vis.var.vy_max", po::value<double>(), "set vy_max")
        ("vis.var.e_min", po::value<double>(), "set e_min")
        ("vis.var.e_max", po::value<double>(), "set e_max")
        ("in.var.x", po::value<double>(), "set x")
        ("in.var.y", po::value<double>(), "set y")
        ("in.var.vx", po::value<double>(), "set vx")
        ("in.var.vy", po::value<double>(), "set vy")
        ("in.var.e", po::value<double>(), "set energy")
        ("distance.flag", po::value<bool>()->default_value(0), "distance flag")
        ("distance.L", po::value<double>(), "distance between primaries")
        ("distance.d1", po::value<double>(), "distance from first primary")
        ("distance.d2", po::value<double>(), "distance from second primary")
        ("matlab.flag", po::value<bool>()->default_value(0), "if true use matlab to plot results")
        ;

    po::variables_map vm;
    std::string config_file;
    config_file="Configuration.cfg";
    std::ifstream ifs(config_file.c_str());
    store(parse_config_file(ifs, options), vm);
    notify(vm);

    // Load parameters
    if (vm.count("parameters.mu")) {mu=vm["parameters.mu"].as<double>();}
    else    {std::cout<<"Mass parameter missing!";
            return 1;}
    if (vm.count("parameters.ecc")) {ecc=vm["parameters.ecc"].as<double>();}
    else    {std::cout<<"Eccentricity missing!";
            return 1;}
    if (vm.count("parameters.type")){field_type=vm["parameters.type"].as<std::string>();}
    if (vm.count("parameters.t0")){t0=vm["parameters.t0"].as<double>();}
    if (vm.count("parameters.tf")){tf=vm["parameters.tf"].as<double>();}
    if (vm.count("parameters.DT")) {DT=vm["parameters.DT"].as<double>();}
    else if (field_type=="FTLE"){std::cout<<"Integration interval missing!\n";    //FIXME Trovare nome più adatto o eliminare errore mettendo un default
                                return 1;}
    if(field_type=="FILE")
            {if (vm.count("parameters.n_iterations")) {n_iterations=vm["parameters.n_iterations"].as<int>();}}
    if (field_type!="FILE" && field_type!="FTLE")
    {std::cout<<"Unrecognized field type\n";
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
    if (n_frames>1 && matlab_flag){std::cout<<"Automatic Matlab visualization not supported if n_frames>1\n";
                                    return 1;}

    printf("mu=%.2f ecc=%.2f DT=%.2f n_frames=%i\n",mu,ecc,DT,n_frames);

    // Load vectors dimensions and input values

    // nx check and load
    if (vm.count("vis.var.nx"))
        {
            flags[0]=1;
            vis_count++;
            nx=vm["vis.var.nx"].as<int>();
            if (vm.count("vis.var.x_min")){x_min=vm["vis.var.x_min"].as<double>();}
            else    {std::cout<<"x range missing!";
                    return 1;}
            if (vm.count("vis.var.x_max")){x_max=vm["vis.var.x_max"].as<double>();}
            else    {std::cout<<"x range missing!";
                    return 1;}
            dx=(x_max-x_min)/(nx-1);
            if (dx<=0)
                {std::cout<<"x_max must be greater than x_min!\n";
                return 1;}
        }
    else
        {
            if (vm.count("in.var.x"))
                {
                    flags[0]=2;
                    in_count++;
                    X_0=vm["in.var.x"].as<double>();
                    nx=3;
                }
            else   {flags[0]=0;}
        }

    // ny check and load
    if (vm.count("vis.var.ny"))
        {
            flags[1]=1;
            vis_count++;
            ny=vm["vis.var.ny"].as<int>();
            if (vm.count("vis.var.y_min")){y_min=vm["vis.var.y_min"].as<double>();}
            else    {std::cout<<"y range missing!\n";
                    return 1;}
            if (vm.count("vis.var.y_max")){y_max=vm["vis.var.y_max"].as<double>();}
            else    {std::cout<<"y range missing!\n";
                    return 1;}
            dy=(y_max-y_min)/(ny-1);
            if (dy<=0)
                {std::cout<<"y_max must be greater than y_min!\n";
                return 1;}
        }
    else
        {
            if (vm.count("in.var.y"))
                {
                    flags[1]=2;
                    in_count++;
                    Y_0=vm["in.var.y"].as<double>();
                    ny=3;
                }
            else   {flags[1]=0;}
        }

    // nvx check and load
    if (vm.count("vis.var.nvx"))
        {
            flags[2]=1;
            vis_count++;
            nvx=vm["vis.var.nvx"].as<int>();
            n1=nvx;
            if (vm.count("vis.var.vx_min")){vx_min=vm["vis.var.vx_min"].as<double>();}
            else    {std::cout<<"vx range missing!\n";
                    return 1;}
            if (vm.count("vis.var.vx_max")){vx_max=vm["vis.var.vx_max"].as<double>();}
            else    {std::cout<<"vx range missing!\n";
                    return 1;}
            dvx=(vx_max-vx_min)/(nvx-1);
            if (dvx<=0)
                {std::cout<<"vx_max must be greater than vx_min!\n";
                return 1;}
        }
    else
        {
            if (vm.count("in.var.vx"))
                {
                    flags[2]=2;
                    in_count++;
                    VX_0=vm["in.var.vx"].as<double>();
                    nvx=3;
                    n1=nvx;
                }
            else   {flags[2]=0;}
        }

        // nvy check and load
    if (vm.count("vis.var.nvy"))
        {
            flags[3]=1;
            vis_count=vis_count+1;
            nvy=vm["vis.var.nvy"].as<int>();
            if (flags[2]){n2=nvy;}
            else{n1=nvy;}
            if (vm.count("vis.var.vy_min")){vy_min=vm["vis.var.vy_min"].as<double>();}
            else    {std::cout<<"vy range missing!\n";
                    return 1;}
            if (vm.count("vis.var.vy_max")){vy_max=vm["vis.var.vy_max"].as<double>();}
            else    {std::cout<<"vy range missing!\n";
                    return 1;}
            dvy=(vy_max-vy_min)/(nvy-1);
            if (dvy<=0)
                {std::cout<<"vy_max must be greater than vy_min!\n";
                return 1;}
        }
    else
        {
            if (vm.count("in.var.vy"))
                {
                    flags[3]=2;
                    in_count=in_count+1;
                    VY_0=vm["in.var.vy"].as<double>();
                    nvy=3;
                    if (flags[2]){n2=nvy;}
                    else{n1=nvy;}

                }
            else   {flags[3]=0;}
        }

        // ne check and load
    if (vm.count("vis.var.ne"))
        {
            flags[4]=1;
            vis_count=vis_count+1;
            ne=vm["vis.var.ne"].as<int>();
            n2=ne;
            if (vm.count("vis.var.e_min")){e_min=vm["vis.var.e_min"].as<double>();}
            else    {std::cout<<"energy range missing!\n";
                    return 1;}
            if (vm.count("vis.var.e_max")){e_max=vm["vis.var.e_max"].as<double>();}
            else    {std::cout<<"energy range missing!\n";
                    return 1;}
            de=(e_max-e_min)/(ne-1);
            if (de<=0)
                {std::cout<<"e_max must be greater than e_min!\n";
                return 1;}
        }
    else
        {
            if (vm.count("in.var.e"))
                {
                    flags[4]=2;
                    in_count=in_count+1;
                    E_0=vm["in.var.e"].as<double>();
                    ne=3;
                    n2=ne;
                }
            else   {flags[4]=0;}
        }

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

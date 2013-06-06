#include <boost/program_options.hpp>
#include <fstream>
#include <stdio.h>
namespace po = boost::program_options;
double mu, ecc, x_min, x_max, y_min, y_max, vx_min, vx_max, vy_min, vy_max, e_min, e_max, X_0, Y_0, VX_0, VY_0, E_0, DT,
        dx, dy, dvx, dvy, de, x_mint, x_maxt, y_mint, y_maxt, vx_mint, vx_maxt, vy_mint, vy_maxt, e_mint, e_maxt, X_0t, Y_0t,
        VX_0t, VY_0t, E_0t, dxt, dyt, dvxt, dvyt, det;
int nx, ny, nvx, nvy, ne, n1, n2, n_tot, n_frames, nxt, nyt, nvxt, nvyt, net, n_tot_t, n1t, n2t;
int flags[5], flags_t[5];
int vis_count=0;        //Intero in barba a Calvino! (non ho potuto resistere :-) )
int in_count=0;
bool tracers_flag;

int configuration_load ()
{
    po::options_description options("Allowed options");
    options.add_options()
        ("parameters.mu", po::value<double>(),"mu")
        ("parameters.ecc", po::value<double>(),"ecc")
        ("parameters.DT", po::value<double>(),"DT")
        ("parameters.n_frames", po::value<int>(),"set n_frames")
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
        ("tracers.nx", po::value<int>(), "set nxt")
        ("tracers.ny", po::value<int>(), "set nyt")
        ("tracers.nvx", po::value<int>(), "set nvxt")
        ("tracers.nvy", po::value<int>(), "set nvyt")
        ("tracers.ne", po::value<int>(), "set net")
        ("tracers.x_min", po::value<double>(), "set x_min for tracers")
        ("tracers.x_max", po::value<double>(), "set x_max for tracers")
        ("tracers.y_min", po::value<double>(), "set y_min for tracers")
        ("tracers.y_max", po::value<double>(), "set y_max for tracers")
        ("tracers.vx_min", po::value<double>(), "set vx_min for tracers")
        ("tracers.vx_max", po::value<double>(), "set vx_max for tracers")
        ("tracers.e_min", po::value<double>(), "set e_min for tracers")
        ("tracers.e_max", po::value<double>(), "set e_max for tracers")
        ("tracers.x", po::value<double>(), "set x for tracers")
        ("tracers.y", po::value<double>(), "set y for tracers")
        ("tracers.vx", po::value<double>(), "set vx for tracers")
        ("tracers.vy", po::value<double>(), "set vy for tracers")
        ("tracers.e", po::value<double>(), "set energy for tracers")
        ("tracers.flag", po::value<bool>()->default_value(0),"tracers' flag")
        ;

    po::variables_map vm;
    std::string config_file;
    config_file="Configuration.cfg";
    std::ifstream ifs(config_file.c_str());
    store(parse_config_file(ifs, options), vm);
    notify(vm);

    // Load parameters
    if (vm.count("parameters.mu")){mu=vm["parameters.mu"].as<double>();}
    else    {std::cout<<"Mass parameter missing!";
            return 1;}
    if (vm.count("parameters.ecc")){ecc=vm["parameters.ecc"].as<double>();}
    else    {std::cout<<"Eccentricity missing!";
            return 1;}
    if (vm.count("parameters.DT")){DT=vm["parameters.DT"].as<double>();}
    else    {std::cout<<"Integration extremes missing!";    //FIXME Trovare nome più adatto o eliminare errore mettendo un default
            return 1;}
    if (vm.count("parameters.n_frames")){n_frames=vm["parameters.n_frames"].as<int>();}
    else    {std::cout<<"Number of frames missing!";
            return 1;}
    printf("mu=%.2f ecc=%.2f DT=%.2f n_frames=%i\n",mu,ecc,DT,n_frames);
    // Load vectors dimensions and input values

    // nx check and load
    if (vm.count("vis.var.nx"))
        {
            flags[0]=1;
            vis_count=vis_count+1;
            nx=vm["vis.var.nx"].as<int>();
            if (vm.count("vis.var.x_min")){x_min=vm["vis.var.x_min"].as<double>();}
            else    {std::cout<<"x range missing!";
                    return 1;}
            if (vm.count("vis.var.x_max")){x_max=vm["vis.var.x_max"].as<double>();}
            else    {std::cout<<"x range missing!";
                    return 1;}
            dx=(x_max-x_min)/(nx-1);
            if (dx<0)
                {std::cout<<"x_max must be greater than x_min!";
                return 1;}
        }
    else
        {
            if (vm.count("in.var.x"))
                {
                    flags[0]=2;
                    in_count=in_count+1;
                    X_0=vm["in.var.x"].as<double>();
                    nx=3;
                }
            else   {flags[0]=0;}
        }

    // ny check and load
    if (vm.count("vis.var.ny"))
        {
            flags[1]=1;
            vis_count=vis_count+1;
            ny=vm["vis.var.ny"].as<int>();
            if (vm.count("vis.var.y_min")){y_min=vm["vis.var.y_min"].as<double>();}
            else    {std::cout<<"y range missing!";
                    return 1;}
            if (vm.count("vis.var.y_max")){y_max=vm["vis.var.y_max"].as<double>();}
            else    {std::cout<<"y range missing!";
                    return 1;}
            dy=(y_max-y_min)/(ny-1);
            if (dy<0)
                {std::cout<<"y_max must be greater than y_min!";
                return 1;}
        }
    else
        {
            if (vm.count("in.var.y"))
                {
                    flags[1]=2;
                    in_count=in_count+1;
                    Y_0=vm["in.var.y"].as<double>();
                    ny=3;
                }
            else   {flags[1]=0;}
        }

    // nvx check and load
    if (vm.count("vis.var.nvx"))
        {
            flags[2]=1;
            vis_count=vis_count+1;
            nvx=vm["vis.var.nvx"].as<int>();
            n1=nvx;
            if (vm.count("vis.var.vx_min")){vx_min=vm["vis.var.vx_min"].as<double>();}
            else    {std::cout<<"vx range missing!";
                    return 1;}
            if (vm.count("vis.var.vx_max")){vx_max=vm["vis.var.vx_max"].as<double>();}
            else    {std::cout<<"vx range missing!";
                    return 1;}
            dvx=(vx_max-vx_min)/(nvx-1);
            if (dvx<0)
                {std::cout<<"vx_max must be greater than vx_min!";
                return 1;}
        }
    else
        {
            if (vm.count("in.var.vx"))
                {
                    flags[2]=2;
                    in_count=in_count+1;
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
            else    {std::cout<<"vy range missing!";
                    return 1;}
            if (vm.count("vis.var.vy_max")){vy_max=vm["vis.var.vy_max"].as<double>();}
            else    {std::cout<<"vy range missing!";
                    return 1;}
            dvy=(vy_max-vy_min)/(nvy-1);
            if (dvy<0)
                {std::cout<<"vy_max must be greater than vy_min!";
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
            else    {std::cout<<"energy range missing!";
                    return 1;}
            if (vm.count("vis.var.e_max")){e_max=vm["vis.var.e_max"].as<double>();}
            else    {std::cout<<"energy range missing!";
                    return 1;}
            de=(e_max-e_min)/(ne-1);
            if (de<0)
                {std::cout<<"e_max must be greater than e_min!";
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
        /*
    if (vis_count>2)
    {
        std::cout<<"Too many visualization variables!";
        return 1;
    }
    */

    if (vis_count<2)
    {
        std::cout<<"Too few visualization variables!";
        return 1;
    }

    if (in_count>2)
    {
        std::cout<<"Too many input variables!";
        return 1;
    }
/*
    if (in_count<2)
    {
        std::cout<<"Too few input variables!";
        return 1;
    }
*/
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

    /* ---TRACERS--- */
if (vm.count("tracers.flag")){tracers_flag=vm["tracers.flag"].as<bool>();}
//std::cout<<"tracers.flag="<<tracers_flag;
if (tracers_flag){
/* nx for tracers check and load*/
    if (vm.count("tracers.nx"))
        {
            flags_t[0]=1;
            //std::cout<<"flags_t[0]="<<flags_t[0];
            //vis_count=vis_count+1;
            nxt=vm["tracers.nx"].as<int>();
            if (vm.count("tracers.x_min")){x_mint=vm["tracers.x_min"].as<double>();}
            else    {std::cout<<"x range for tracers missing!";
                    return 1;}
            if (vm.count("tracers.x_max")){x_maxt=vm["tracers.x_max"].as<double>();}
            else    {std::cout<<"x range for tracers missing!";
                    return 1;}
            dxt=(x_maxt-x_mint)/(nxt-1);
            if (dxt<0)
                {std::cout<<"error in tracers extremes: x_max must be greater than x_min!";
                return 1;}
        }
    else
        {
            if (vm.count("tracers.x"))
                {
                    flags_t[0]=2;
                    //in_count=in_count+1;
                    X_0t=vm["tracers.x"].as<double>();
                    nxt=1;
                }
            else   {flags_t[0]=0;}
        }
/* ny for tracers check and load*/
    if (vm.count("tracers.ny"))
        {
            flags_t[1]=1;
            //vis_count=vis_count+1;
            nyt=vm["tracers.ny"].as<int>();
            if (vm.count("tracers.y_min")){y_mint=vm["tracers.y_min"].as<double>();}
            else    {std::cout<<"y range for tracers missing!";
                    return 1;}
            if (vm.count("tracers.y_max")){y_maxt=vm["tracers.y_max"].as<double>();}
            else    {std::cout<<"y range for tracers missing!";
                    return 1;}
            dxt=(y_maxt-y_mint)/(nyt-1);
            if (dyt<0)
                {std::cout<<"error in tracers extremes: y_max must be greater than y_min!";
                return 1;}
        }
    else
        {
            if (vm.count("tracers.y"))
                {
                    flags_t[1]=2;
                    //in_count=in_count+1;
                    Y_0t=vm["tracers.y"].as<double>();
                    nyt=1;
                }
            else   {flags_t[1]=0;}
        }

/* nvx for tracers check and load*/
    if (vm.count("tracers.nvx"))
        {
            flags_t[2]=1;
            //vis_count=vis_count+1;
            nvxt=vm["tracers.nvx"].as<int>();
            n1t=nvxt;
            if (vm.count("tracers.vx_min")){vx_mint=vm["tracers.vx_min"].as<double>();}
            else    {std::cout<<"vx range for tracers missing!";
                    return 1;}
            if (vm.count("tracers.vx_max")){vx_maxt=vm["tracers.vx_max"].as<double>();}
            else    {std::cout<<"vx range for tracers missing!";
                    return 1;}
            dvxt=(vx_maxt-vx_mint)/(nvxt-1);
            if (dvxt<0)
                {std::cout<<"error in tracers extremes: vx_max must be greater than vx_min!";
                return 1;}
        }
    else
        {
            if (vm.count("tracers.vx"))
                {
                    flags_t[2]=2;
                    //in_count=in_count+1;
                    VX_0t=vm["tracers.vx"].as<double>();
                    nvxt=1;
                    n1t=nvxt;
                }
            else   {flags_t[2]=0;}
        }

/* nvy for tracers check and load*/
    if (vm.count("tracers.nvy"))
        {
            flags_t[3]=1;
            //vis_count=vis_count+1;
            nvxt=vm["tracers.nvy"].as<int>();
            if (flags_t[2]){n2t=nvyt;}
            else{n1t=nvyt;}
            if (vm.count("tracers.vy_min")){vy_mint=vm["tracers.vy_min"].as<double>();}
            else    {std::cout<<"vy range for tracers missing!";
                    return 1;}
            if (vm.count("tracers.vy_max")){vy_maxt=vm["tracers.vy_max"].as<double>();}
            else    {std::cout<<"vy range for tracers missing!";
                    return 1;}
            dvyt=(vy_maxt-vy_mint)/(nvyt-1);
            if (dvyt<0)
                {std::cout<<"error in tracers extremes: vy_max must be greater than vy_min!";
                return 1;}
        }
    else
        {
            if (vm.count("tracers.vy"))
                {
                    flags_t[3]=2;
                    //in_count=in_count+1;
                    VY_0t=vm["tracers.vy"].as<double>();
                    nvyt=1;
                    if (flags_t[2]){n2t=nvyt;}
                    else{n1t=nvyt;}
                }
            else   {flags_t[3]=0;}
        }

/* ne for tracers check and load*/
    if (vm.count("tracers.ne"))
        {
            flags_t[4]=1;
            //vis_count=vis_count+1;
            net=vm["tracers.ne"].as<int>();
            n2t=net;
            if (vm.count("tracers.e_min")){e_mint=vm["tracers.e_min"].as<double>();}
            else    {std::cout<<"energy range for tracers missing!";
                    return 1;}
            if (vm.count("tracers.e_max")){e_maxt=vm["tracers.e_max"].as<double>();}
            else    {std::cout<<"energy range for tracers missing!";
                    return 1;}
            det=(e_maxt-e_mint)/(net-1);
            if (det<0)
                {std::cout<<"error in tracers extremes: e_max must be greater than e_min!";
                return 1;}
        }
    else
        {
            if (vm.count("tracers.e"))
                {
                    flags_t[4]=2;
                    //in_count=in_count+1;
                    E_0t=vm["tracers.e"].as<double>();
                    net=1;
                    n2t=net;
                }
            else   {flags_t[4]=0;}
        }

    /*Create dimension of missing vector for tracers */
    for (int i=0;i<5;i++)
    {
        if (flags_t[i]==0)
        {
            switch (i)
            {
                case 0: //not supported
                    std::cout<<"Case not supported: please provide a value for both x and y";
                    return 1;
                case 1: //not supported
                    std::cout<<"Case not supported: please provide a value for both x and y";
                    return 1;
                case 2:
                    nvxt=nxt*nyt*nvyt*net;
                    n_tot_t=nvxt;
                    break;
                case 3:
                    nvyt=nxt*nyt*nvxt*net;
                    n_tot_t=nvyt;
                    break;
                case 4:
                    net=nxt*nyt*nvxt*nvyt;
                    n_tot_t=net;
                    break;
            }
        }
    }
}
    //printf(" nx=%i\n ny=%i\n nvx=%i\n nvy=%i\n ne=%i\n",nx,ny,nvx,nvy,ne);
    return 0;
}

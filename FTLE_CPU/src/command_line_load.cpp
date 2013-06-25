#include <boost/program_options.hpp>
#include <fstream>
#include <stdio.h>
#include <string>
#include "global_var.hpp"

namespace po = boost::program_options;
std::string config_file;
std::string command_line_load (int ac, char* av[])
{
po::options_description cmdline_options("Allowed options on command line");
    cmdline_options.add_options()
        ("input-file", po::value<std::string>()->default_value("Configuration.cfg"), "Name of the input file")
        ;

po::variables_map vm;
po::store(po::parse_command_line (ac,av,cmdline_options), vm);
po::notify(vm);

config_file=vm["input-file"].as<std::string>();
return config_file;
}

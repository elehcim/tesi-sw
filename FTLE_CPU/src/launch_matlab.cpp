#include <stdio.h>  // sprintf
#include <stdlib.h> // system
#include <iostream> // cout
char command[500];
char matlab_path[1000]="addpath(genpath('../../Functions_library'));addpath(genpath('../../Plot_library'));"; // FIXME
int launch_matlab(char *file_name)
{
/* TODO: In future start the standalone Matlab executable*/

sprintf(command,"matlab -nodesktop -nosplash -r \"%s cr_fig('file','%s');disp('Press any key');pause;quit\"",matlab_path,file_name);
//std::cout << command <<"\n";
system(command);
return 0;
}

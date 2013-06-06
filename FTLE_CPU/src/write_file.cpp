#include <stdio.h>  // sprintf
#include <iostream> // cout
#include <fstream>
#include "global_var.hpp"
#include "type_definitions.hpp"
FILE *ftle_stream;
/*
Output file structure:

coordinate1   coordinate2   (coordinate3)  (coordinate4)    ftle

Number of coordinates varies with the number of gridded variables. Coordinates always appear in this order: x, y, vx, vy, e.
*/

int write_file(char* file_header, char* file_name, int id, double1d x_0,
                double1d y_0, double1d vx_0, double1d vy_0, double1d e_0, double4d ftle)
{
    ftle_stream=fopen(file_name,"w");
    fprintf(ftle_stream,"%s",file_header);

    if (id==11220 || id==11202 || id==11022)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], ftle[i][j][1][1]);
            }
        }
    }

    if(id==12120 || id==12102)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nvx; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", x_0[i], vx_0[j], ftle[i][1][j][1]);
            }
        }
    }

    if(id==12210)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nvy; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", x_0[i], vy_0[j], ftle[i][1][1][j]);
            }
        }
    }

    if(id==12012)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nvy; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", x_0[i], vy_0[j], ftle[i][1][j][1]);
            }
        }
    }

    if(id==12201 || id==12021)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ne; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", x_0[i], e_0[j], ftle[i][1][1][j]);
            }
        }
    }

    if (id==21120 || id==21102)
    {
        for (int i=0; i<ny; i++)
        {
            for (int j=0; j<nvx; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", y_0[i], vx_0[j], ftle[1][i][j][1]);
            }
        }
    }

    if (id==21012)
    {
        for (int i=0; i<ny; i++)
        {
            for (int j=0; j<nvy; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", y_0[i], vy_0[j], ftle[1][i][j][1]);
            }
        }
    }

    if (id==21210)
    {
        for (int i=0; i<ny; i++)
        {
            for (int j=0; j<nvy; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", y_0[i], vy_0[j], ftle[1][i][1][j]);
            }
        }
    }

    if (id==21201 || id==21021)
    {
        for (int i=0; i<ny; i++)
        {
            for (int j=0; j<ne; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", y_0[i], e_0[j], ftle[1][i][1][j]);
            }
        }
    }

    if (id==22110)
    {
        for (int i=0; i<nvx; i++)
        {
            for (int j=0; j<nvy; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", vx_0[i], vy_0[j], ftle[1][1][i][j]);
            }
        }
    }

    if (id==22101)
    {
        for (int i=0; i<nvx; i++)
        {
            for (int j=0; j<ne; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", vx_0[i], e_0[j], ftle[1][1][i][j]);
            }
        }
    }

    if (id==22011)
    {
        for (int i=0; i<nvy; i++)
        {
            for (int j=0; j<ne; j++)
            {
                fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\n", vx_0[i], e_0[j], ftle[1][1][i][j]);
            }
        }
    }
    //3D Cases
    if (id==11120 || id==11102)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<nvx; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vx_0[k], ftle[i][j][k][1]);
                }
            }
        }
    }

    if (id==11210)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<nvy; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vy_0[k], ftle[i][j][1][k]);
                }
            }
        }
    }

    if (id==11012)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<nvy; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vy_0[k], ftle[i][j][k][1]);
                }
            }
        }
    }

    if (id==11021 || id==11201)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<ne; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], e_0[k], ftle[i][j][1][k]);
                }
            }
        }
    }

    if (id==12110)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nvx; j++)
            {
                for (int k=0; k<nvy; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], vx_0[j], vy_0[k], ftle[i][1][j][k]);
                }
            }
        }
    }

    if (id==12101)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nvx; j++)
            {
                for (int k=0; k<ne; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], vx_0[j], e_0[k], ftle[i][1][j][k]);
                }
            }
        }
    }

    if (id==12011)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nvy; j++)
            {
                for (int k=0; k<ne; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], vy_0[j], e_0[k], ftle[i][1][j][k]);
                }
            }
        }
    }

    if (id==12011)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<nvy; j++)
            {
                for (int k=0; k<ne; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], vy_0[j], e_0[k], ftle[i][1][j][k]);
                }
            }
        }
    }

    if (id==21110)
    {
        for (int i=0; i<ny; i++)
        {
            for (int j=0; j<nvx; j++)
            {
                for (int k=0; k<nvy; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", y_0[i], vx_0[j], vy_0[k], ftle[1][i][j][k]);
                }
            }
        }
    }

    if (id==21101)
    {
        for (int i=0; i<ny; i++)
        {
            for (int j=0; j<nvx; j++)
            {
                for (int k=0; k<ne; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", y_0[i], vx_0[j], e_0[k], ftle[1][i][j][k]);
                }
            }
        }
    }

    if (id==21011)
    {
        for (int i=0; i<ny; i++)
        {
            for (int j=0; j<nvy; j++)
            {
                for (int k=0; k<ne; k++)
                {
                    fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\n", y_0[i], vy_0[j], e_0[k], ftle[1][i][j][k]);
                }
            }
        }
    }

//4D case
    if (id==11110)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<nvx; k++)
                {
                    for (int l=0; l<nvy; l++)
                    {
                        fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vx_0[k], vy_0[l], ftle[i][j][k][l]);
                    }
                }
            }
        }
    }

    if (id==11101)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<nvx; k++)
                {
                    for (int l=0; l<ne; l++)
                    {
                        fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vx_0[k], e_0[l], ftle[i][j][k][l]);
                    }
                }
            }
        }
    }

    if (id==11011)
    {
        for (int i=0; i<nx; i++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int k=0; k<nvy; k++)
                {
                    for (int l=0; l<ne; l++)
                    {
                        fprintf(ftle_stream,"%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n", x_0[i], y_0[j], vy_0[k], e_0[l], ftle[i][j][k][l]);
                    }
                }
            }
        }
    }
fclose(ftle_stream);
return 0;
}


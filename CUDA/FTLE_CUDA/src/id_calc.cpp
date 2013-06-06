#include <math.h>
int id_calc (int *arr)
{
    int s=0;
        for (int i=0; i<5; i++)
        {
            s=s+arr[i]*pow(10,4-i);
        }
    return s;
}

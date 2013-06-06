How to write configuration file:

The FTLE executable, reads input values from configuration.cfg, which must be placed in the same directory as the executable itself.
In order for the program to read values correctly, those MUST be written as follows:
parameter_name=parameter_value

Configuration file is divided into several sections, identified by names in squared brackets.
Each section admits only its own parameters. Below there is a list of all the allowable parameters for each section, including a brief description of their meaning. In brackets it is reported the value type of the parameter.

[parameters]		This section groups the physical parameters of the system along with some computational parameters.
mu=(double)		Mass parameter of the system
ecc=(double)		Eccentricity of the primaries' orbit
DT=(double)		Integration interval: the interval over which the integration is performed in order to compute FTLE field
n_frames=(int)		Number of frames. Must be greater than 0. If it is greater than 1 the output wil be a video, otherwise, a static image will be created.
n_cores=(int)		Number of threads over which the integration will be parallelized. Note that it is independent on actual number of computer's hardware cores, so benefits may not exist if this number is increased over the number of actual hardware cores.

[vis.var]			In this section parameters about the variables to be visualized must be specified. It is necessary to have at least 2 (and no more than 4) variables specified and their field's extremes.
nx=(int)			Resolution of x range. The range between x_min and x_max will be divided into nx points.
x_min=(double)	Lower extreme of the x range.
x_max(double)		Upper extreme of the x range.
ny=(int)			Resolution of y range. The range between y_min and y_max will be divided into ny points.
y_min=(double)	Lower extreme of the y range.
y_max(double)		Upper extreme of the y range.
nvx=(int)			Resolution of vx range. The range between vx_min and vx_max will be divided into nvx points.
vx_min=(double)	Lower extreme of the vx range.
vx_max(double)	Upper extreme of the vx range.
nvy=(int)			Resolution of vy range. The range between vy_min and vy_max will be divided into nvy points.
vy_min=(double)	Lower extreme of the vy range.
vy_max(double)	Upper extreme of the vy range.
ne=(int)			Resolution of energy range. The range between e_min and e_max will be divided into ne points.
e_min=(double)	Lower extreme of the energy range.
e_max(double)		Upper extreme of the energy range.

[in.var]			In this section must be specified values for variables which have to be fixed. The sum of the number of variables specified here and above must not be more than 4.
x=(double)		Fixed value for x.
y=(double)		Fixed value for y.
vx=(double)		Fixed value for vx.
vy=(double)		Fixed value for vy.
e=(double)		Fixed value for energy.

[distance]			This section groups parameters about the distance from the primaries check during integration.
flag=(bool)		If true (1) distance check will be performed; if false (0) it will be ignored.
d1=(double)		Forbidden distance from first primary (dimensionless).
d2=(double)		Forbidden distance from second primary (dimensionless).

[matlab]
flag=(bool)		If true (1) matlab will be automatically launched at the end of the computation in order to plot the results; if false (0) the program is ended without launching matlab.

It is possible to write comments everywhere in the configuration file by prepending #.
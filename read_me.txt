The code is written by Mariia Dmitrieva in 2017

It allows to calculate analytical solution for backscattering of a pulse.
The example code is located in example.m - change the parameters for your sphere and run it.
The code depends on all the other functions.

reflectionNumerical.m - calculate the reflection from the sphere
formfunction.m -  computes form function based on the characterics of a 2 layer spherical object
generate_chirp.m - generates orinial pulse ( linear chirp) with a given characteristics

backscattering.m - calculates backscattering itself based on the form function and initial pulse

all the other functions are a support for the formfunction.m and reflectionNumerical.m 



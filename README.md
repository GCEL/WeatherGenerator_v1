# WeatherGenerator_v1
A simple empirical model for down-scaling daily meteorological information from daily information to hourly time step.

First coded by Oliver Browne, UoE
Subseqent modifications by T. L. Smallman (UoE, t.l.smallman@ed.ac.uk)

This repository contains three files.
1) weather_generator_R_call.r; A control interface file written in R to call the weather generator functions which are written in Fortran
2) wg_interface.f90; The Fortran interface file first called by R, but responsible for calling the specific Fortran code.
3) weather_generator.f90; The weather generator source code itself

Inputs:
Maximum/average/minimum daily temperature (Celcius)
Maximum/average/minimum vapour pressure deficit (kPa)
Mean daily shortwave radiation (W/m2)
Daily total precipitation (mm)

Outputs:
Meteorology downscaled to 24 hourly time steps.

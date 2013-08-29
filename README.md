FTLE Calculator
===============
FTLE calculator for astrodynamic purposes.

Directory Structure
================

    root/
      |-CUDA/             #  CUDA sources (non-performing)
      |-FTLE_CPU/         #  Main program
      |-FTLE_Matlab/      #  Matlab implementation of FTLE calculation
      |-MatlabScripts/    #  Various scripts useful for Astrodynamics
      |-MissionAnalysis/  #  Scripts for various computations for 
                             preliminarlily design a Mission to Jupiter
      |-Plot_ftle/        #  Scripts to visualize FTLE fields from output
                             files generated by FTLE_CPU
      |-README.md
      |-Tracers/          #  Scripts to visualize tracers and trajectories from FTLE fields
      |-Tutorials/        #  Tutorials to build the sources

Matlab Configuration
====================
Open the file `Plot_ftle/folder.m` and modify either line 7 (for Windows user) or 11 (for Linux user) writing the path to the root directory.

After a succesful computation of a ftle field, add to the Matlab path the folders (with subfolders):

  * To view FTLE fields
    - Plot_ftle
  * To play with tracers
    - Plot_ftle
    - Tracers
  * To use the Mission analysis routines (not fully implemented for a general use, used for thesis work)
    - Plot_ftle
    - Tracers
    - MissionAnalysis
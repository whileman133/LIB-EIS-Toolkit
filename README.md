# EIS Toolkit for Lithium-Ion Battery Cells

This repository contains MATLAB and Python code for calculating battery impedance from PDAE models constructed in COMSOL Multiphysics and PyBaMM. Also provided is a closed-form transfer function (TF) model for evaluating the impedance of lithium-ion battery (LIB) cells from the lumped-parameter model specified by Plett and Trimboli.[^1]

## MATLAB Resources

MATLAB code is provided for calculating battery impedance from COMSOL models using two approaches: time-domain simulation (TD) and frequency-domain linear perturbation (LP). Both approaches require the [COMSOL LiveLink for MATLAB](https://www.comsol.com/livelink-for-matlab) and have been tested with COMSOL 6.3. We provide MATLAB utility functions to construct the lumped-parameter LIB model in COMSOL.

## Python Resources

Python code is provided for calculating battery impedance from PyBaMM models using frequency-domain linear perturbation (LP). We provide two lumped-parameter models: one for LIB[^1] 

[^1]: Plett, Gregory L., and M. Scott Trimboli. _Battery Management Systems, Volume III: Physics-Based Methods_. Vol.&nbsp;3. Artech House, 2024.

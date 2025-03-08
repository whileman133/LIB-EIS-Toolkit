# EIS Toolkit for Lithium-Ion Battery Cells

This repository contains MATLAB and Python code for calculating battery impedance from PDAE models constructed in COMSOL Multiphysics and PyBaMM. Also provided is a closed-form transfer function (TF) model for evaluating the impedance of lithium-ion battery (LIB) cells from the lumped-parameter model specified by Plett and Trimboli.[^1]

## MATLAB Resources

MATLAB code is provided for calculating battery impedance from COMSOL models using two approaches: time-domain simulation (TD) and frequency-domain linear perturbation (LP). Both approaches require the [COMSOL LiveLink for MATLAB](https://www.comsol.com/livelink-for-matlab) and have been tested with COMSOL 6.3 and MATLAB R2021b. We provide MATLAB utility functions to construct the Plett and Trimboli model for LIB in COMSOL.

## Python Resources

Python code is provided for calculating battery impedance from PyBaMM models using frequency-domain linear perturbation (LP). We provide two lumped-parameter PyBaMM models: one implementing the Plett and Trimboli LIB model and another implementing the lithium-metal battery (LMB) model specified by Hileman et al.[^2] The code has been tested with Python version 3.12.7 and the requirements specified in `requirements.txt`.

[^1]: Plett, Gregory L., and M. Scott Trimboli. _Battery Management Systems, Volume III: Physics-Based Methods_. Vol.&nbsp;3. Artech House, 2024.

[^2]: Hileman, Wesley A., M. Scott Trimboli, and Gregory L. Plett. "Estimating the Values of the PDE Model Parameters of Rechargeable Lithium-Metal Battery Cells Using Linear Electrochemical Impedance Spectroscopy." ASME Letters in Dynamic Systems and Control 4, no. 4 (2024).

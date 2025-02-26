"""
dempPerturbation.py

Demonstrate perturbation analysis of a PyBaMM model using PyBaMM-EIS.

-- Changelog --
2025.02.21 | Created | WH
"""
import os

import numpy as np
import matplotlib.pyplot as plt
import scipy

import util
from cellparams import load_cell_params
from util import sim_eis

if __name__ == '__main__':
    params = load_cell_params(os.path.join('XLSX_CELLDEFS','cellNMC30.xlsx'))
    Rc = params.const.Rc()

    # Calculate impedance using PyBaMM-EIS (see utility function sim_eis).
    sim = sim_eis(params, np.logspace(-3.2,4,50), 50)
    Zcell = sim.solution + Rc

    # Load COMSOL perturbation data for comparison.
    _, comsolFreq, comsolZcell = util.load_mat_comsol_eis_sim(os.path.join('MAT_SIMDATA', 'simEIS-50pct.mat'))

    plt.style.use(['tableau-colorblind10', './thesisformat-lg.mplstyle'])
    fig, ax = plt.subplots(constrained_layout=True)
    ax.set_box_aspect(1 / scipy.constants.golden)
    plt.plot(comsolZcell.real * 1000, -comsolZcell.imag * 1000, '.', label='COMSOL')
    plt.plot(Zcell.real*1000, -Zcell.imag*1000, 'x', label='PyBaMM')
    plt.axis('equal')
    plt.xlabel(r"$Z_\mathrm{cell}'\;[m\Omega]$")
    plt.ylabel(r"-$Z_\mathrm{cell}''\;[m\Omega]$")
    plt.title(f"Cell Impedance: Perturbation")
    plt.legend()
    plt.savefig('plots/Zcell-Perturbation-PyBaMM-v-COMSOL.eps')
    plt.show()


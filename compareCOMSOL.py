"""
compareCOMSOL.py

Compare output of COMSOL and PyBaMM LIB cell models.

2025.02.03 | Created | Wesley Hileman <whileman@uccs.edu>
"""
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import util
from cellparams import load_cell_params

if __name__ == "__main__":
    #pybamm.set_logging_level("DEBUG")

    params = load_cell_params(os.path.join('XLSX_CELLDEFS','cellNMC30.xlsx'))
    comsol_data = util.load_mat_comsol_sim(os.path.join('MAT_SIMDATA', 'simCC_1C.mat'))

    time = comsol_data.time
    iapp = comsol_data.iapp
    soc0 = 100
    TdegC = 25

    sim = util.sim_profile(params, time, iapp, soc0, TdegC)
    sim.plot()

    soln = sim.solution
    vcell = soln['vcell [V]'].entries
    vcell = vcell - iapp*params.const.Rc()
    thetass = soln['thetass'](x=comsol_data.xthetass)
    thetass_cc = soln['thetass'](x=float(comsol_data.xthetass[-1]))

    Qdis = np.trapz(iapp, time) / 3600
    socf = soc0 - 100 * Qdis / params.const.Q()
    thetaf = params.pos.theta0() + (socf / 100) * (params.pos.theta100() - params.pos.theta0())
    print('theta(end): ', thetaf)
    rmse = np.sqrt(np.mean((comsol_data.vcell - vcell)**2))
    print('Vcell RMSE:', rmse*1000, 'mV')

    plt.subplots()
    plt.plot(time,iapp)
    plt.title(r'iapp')

    # Plot vcell.
    plt.subplots()
    plt.plot(time, vcell, 'r', label='PyBaMM')
    plt.plot(comsol_data.time, comsol_data.vcell, 'k', label='COMSOL')
    plt.xlabel(r'Time, $t$ [sec]')
    plt.ylabel(r'Cell voltage, $v_{cell}$ [V]')
    plt.title(r'Cell Voltage')
    plt.legend()

    # Plot Thetass
    plt.subplots()
    plt.plot(time, thetass_cc, 'r', label='PyBaMM')
    plt.plot(comsol_data.time, comsol_data.thetass[:,-1], 'k', label='COMSOL')
    plt.xlabel(r'Time, $t$ [sec]')
    plt.ylabel(r'Solid-Surface Composition, $\theta_\mathrm{ss}(\tilde{x}=3)$')
    plt.title(r'Solid-Surface Composition')
    plt.legend()

    plt.show()

    # Plot Thetass across electrode versus time.
    # fig3, ax3 = plt.subplots(2, 1, figsize=(6, 7))
    # ax3[1].plot(time, iapp / cell.cst.QAh)
    # line_iapp, = ax3[1].plot(time[0], iapp[0] / cell.cst.QAh, 'ro')
    # ax3[1].set_title("Applied Current")
    # ax3[1].set_xlabel(r"Time, $t$ [sec]")
    # ax3[1].set_ylabel(r"$i_\mathrm{app}$ [C-rate]")
    # line_thetas, = ax3[0].plot(comsol_data.xthetass, thetass[:, 0], 'r', label='PyBaMM')
    # line_spm, = ax3[0].plot(comsol_data.xthetass, comsol_data.thetass[0, :], 'k', label='COMSOL')
    # ax3[0].legend()
    # ax3[0].set_xlim(min(comsol_data.xthetass), max(comsol_data.xthetass))
    # ax3[0].set_ylim(np.min(thetass), np.max(thetass))
    # ax3[0].set_title("Solid-Surface Composition")
    # ax3[0].set_xlabel(r"Linear Position, $\tilde{x}$")
    # ax3[0].set_ylabel(r"Surface Li Composition, $\theta_\mathrm{ss}=c_\mathrm{ss}/c_\mathrm{s,max}$")
    # plt.tight_layout()
    # # Create animation
    # def update_Thetass(frame):
    #     line_iapp.set_data([comsol_data.time[frame]], [comsol_data.iapp[frame] / cell.cst.QAh])
    #     line_thetas.set_data(comsol_data.xthetass, thetass[:, frame])
    #     line_spm.set_data(comsol_data.xthetass, comsol_data.thetass[frame, :])
    #     return (line_iapp, line_thetas, line_spm)
    # ani = animation.FuncAnimation(
    #     fig3, update_Thetass,
    #     frames=thetass.shape[1],
    #     interval=comsol_data.ts, blit=False, repeat=False)
    # plt.show()


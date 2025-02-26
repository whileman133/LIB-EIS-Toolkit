"""
util.py

Project utilities.
"""
import math
import os
from dataclasses import dataclass

import numpy as np
import pybamm
from scipy.interpolate import PchipInterpolator
import scipy.io as sio

import pybammeis
from cellparams.model import CellParams, _Container
from lumped.libcell import LumpedLIBModel


def get_lib_param_values(params: CellParams, model: LumpedLIBModel, TdegC: float = 25.0):
    """
    Construct PyBaMM parameter values for an LIB cell.
    """

    # Collect parameter values.
    return pybamm.ParameterValues({
        "kD": params.const.kD(TdegC=TdegC),
        "psi [V]": params.const.psi(TdegC=TdegC),
        "Q [Ah]": params.const.Q(TdegC=TdegC),
        "pos_sigma [Ohm-1]": params.pos.sigma(TdegC=TdegC),
        "pos_Ds [s-1]": params.pos.Ds(TdegC=TdegC),
        "pos_kappa [Ohm-1]": params.pos.kappa(TdegC=TdegC),
        "pos_qe [Ah]": params.pos.qe(TdegC=TdegC),
        "pos_theta0": params.pos.theta0(TdegC=TdegC),
        "pos_theta100": params.pos.theta100(TdegC=TdegC),
        "pos_Rf [Ohm]": params.pos.Rf(TdegC=TdegC),
        "pos_Rdl [Ohm]": params.pos.Rdl(TdegC=TdegC),
        "pos_Cdl [F]": params.pos.Cdl(TdegC=TdegC),
        "pos_k0 [A]": params.pos.k0(TdegC=TdegC),
        "pos_alpha": params.pos.alpha(TdegC=TdegC),
        "pos_Uss [V]": pybamm.Interpolant(
            params.pos.Uocp.value.x, params.pos.Uocp.value.y, model.var.pos_thetass, interpolator="linear"
        ),
        "sep_kappa [Ohm-1]": params.sep.kappa(TdegC=TdegC),
        "sep_qe [Ah]": params.sep.qe(TdegC=TdegC),
        "neg_sigma [Ohm-1]": params.neg.sigma(TdegC=TdegC),
        "neg_Ds [s-1]": params.neg.Ds(TdegC=TdegC),
        "neg_kappa [Ohm-1]": params.neg.kappa(TdegC=TdegC),
        "neg_qe [Ah]": params.neg.qe(TdegC=TdegC),
        "neg_theta0": params.neg.theta0(TdegC=TdegC),
        "neg_theta100": params.neg.theta100(TdegC=TdegC),
        "neg_Rf [Ohm]": params.neg.Rf(TdegC=TdegC),
        "neg_Rdl [Ohm]": params.neg.Rdl(TdegC=TdegC),
        "neg_Cdl [F]": params.neg.Cdl(TdegC=TdegC),
        "neg_k0 [A]": params.neg.k0(TdegC=TdegC),
        "neg_alpha": params.neg.alpha(TdegC=TdegC),
        "neg_Uss [V]": pybamm.Interpolant(
            params.neg.Uocp.value.x, params.neg.Uocp.value.y, model.var.neg_thetass, interpolator="linear"
        ),
    })


def get_initial(cell, soc0, TdegC):
    """
    Get PyBaMM initial values for a cell.
    """
    p = _Container()
    p.neg_thetas = cell.neg.theta0(TdegC=TdegC) + (soc0/100)*(cell.neg.theta100(TdegC=TdegC) - cell.neg.theta0(TdegC=TdegC))
    p.neg_Uocp = np.interp(p.neg_thetas, cell.neg.Uocp.value.x, cell.neg.Uocp.value.y)
    p.pos_thetas = cell.pos.theta0(TdegC=TdegC) + (soc0 / 100) * (cell.pos.theta100(TdegC=TdegC) - cell.pos.theta0(TdegC=TdegC))
    p.pos_Uocp = np.interp(p.pos_thetas, cell.pos.Uocp.value.x, cell.pos.Uocp.value.y)
    return p


#
# PyBaMM simulation utilities
#


def sim_eis(
    cell: CellParams,
    freq: np.ndarray,
    soc0: float,
    TdegC: float = 25.0,
):
    """
    Solve the PDE model for EIS.

    :param cell: LMB cell model object containing parameter values.
    :param freq: Frequency vector [Hz].
    :param soc0: Initial cell SOC [%].
    :param TdegC: Cell temperature [degC]. DEFAULT 25.
    :return: PyBaMM-EIS Simulation object
    """

    T = TdegC + 273.15

    # Construct PyBaMM model.
    pde_model = LumpedLIBModel()

    # Collect parameter values.
    param0 = get_initial(cell, soc0, TdegC)
    param = get_lib_param_values(cell, pde_model, TdegC)
    param.update({
        "T [K]": T,
        "neg_thetas0": param0.neg_thetas,
        "neg_Uocp0 [V]": param0.neg_Uocp,
        "pos_thetas0": param0.pos_thetas,
        "pos_Uocp0 [V]": param0.pos_Uocp,
        "iapp [A]": 0,
    }, check_already_exists=False)

    # Run simulation.
    sim = pybammeis.EISSimulation(pde_model, parameter_values=param)
    sim.solve(freq)

    return sim


def sim_profile(
    cell: CellParams,
    time: np.ndarray,
    iapp: np.ndarray,
    soc0: float,
    TdegC: float = 25.0,
):
    """
    Solve the PDE model for the specified parameter values and current profile.

    :param cell: LMB cell model object containing parameter values.
    :param time: Vector of time values [s].
    :param iapp: Vector of applied current values [A].
    :param soc0: Initial cell SOC [%].
    :param TdegC: Cell temperature [degC]. DEFAULT 25.
    :return: PyBaMM Simulation object
    """

    T = TdegC + 273.15

    # Construct PyBaMM model.
    pde_model = LumpedLIBModel()

    # Collect parameter values.
    param0 = get_initial(cell, soc0, TdegC)
    param = get_lib_param_values(cell, pde_model, TdegC)
    param.update({
        "T [K]": T,
        "iapp [A]": pybamm.Interpolant(
            time, iapp, pybamm.t, interpolator="linear"
        ),
        "neg_thetas0": param0.neg_thetas,
        "neg_Uocp0 [V]": param0.neg_Uocp,
        "pos_thetas0": param0.pos_thetas,
        "pos_Uocp0 [V]": param0.pos_Uocp,
    }, check_already_exists=False)

    # Run simulation.
    sim = pybamm.Simulation(pde_model, parameter_values=param)
    sim.solve(t_eval=time)

    return sim


def sim_cc(
    cell: CellParams,
    i_galv: float,
    soc0: float,
    socf: float,
    TdegC: float = 25.0,
    ts: float = 1.0,
    t_rest: float = 0,
):
    """
    Solve the PDE model for a constant-current dis/charge profile.

    :param cell:
    :param i_galv: Constant dis/charge current [A].
    :param soc0: Initial SOC [%].
    :param socf: Final SOC [%].
    :param TdegC: Cell temperature [degC]. DEFAULT 25.
    :param ts: Sampling interval [s]. DEFAULT 1.0.
    :return: PyBaMM Simulation object
    """

    i_galv = abs(i_galv)  # ignore sign, can infer from soc0 and socf

    # Total amount of charge to remove from (+) / add to (-) the cell [Ah].
    QdisAh = (soc0 - socf) / 100 * cell.const.Q(TdegC=TdegC)

    # Total amount of time to spend dis/charging at i_galv [s].
    t_cc = 3600 * abs(QdisAh) / i_galv

    # Build time and applied current vectors.
    time = np.arange(0, t_cc + t_rest + ts, ts)
    iapp = np.zeros_like(time)
    iapp[time <= t_cc] = np.sign(soc0 - socf)*i_galv

    # Run the simulation.
    return sim_profile(cell, time, iapp, soc0, TdegC)


def sim_gitt(
    cell: CellParams,
    t_galv: float,
    t_rest: float,
    i_galv: float,
    soc0: float,
    socf: float,
    TdegC: float = 25.0,
    ts: float = 1.0,
):
    """
    Solve the PDE model for a galvanostatic intermittent titration
    technique (GITT) current profile.

    :param cell: LMB cell model object containing parameter values.
    :param t_galv: Dis/charge interval duration [s].
    :param t_rest: Rest interval duration [s].
    :param i_galv: Current magnitude during dis/charge intervals [A].
    :param soc0: Initial SOC of the cell [%].
    :param socf: Final SOC of the cell [%].
    :param TdegC: Temperature of the cell [degC]. DEFAULT 25.
    :param ts: Sampling interval [s]. DEFAULT 1.0.
    :return: PyBaMM Simulation object
    """

    i_galv = abs(i_galv)  # ignore sign, can infer from soc0 and socf

    # Total amount of charge to remove from (+) / add to (-) the cell [Ah].
    QdisAh = (soc0 - socf)/100 * cell.const.Q(TdegC=TdegC)

    # Total amount of time to spend dis/charging at i_galv [s].
    t_dis = 3600 * abs(QdisAh) / i_galv

    # Number of required dis/charge intervals.
    n_dis = math.ceil(t_dis/t_galv)

    # Total duration of the GITT profile.
    t_gitt = n_dis * (t_galv + t_rest)

    # Build time and applied current vectors.
    time = np.arange(0, t_gitt + ts, ts)
    iapp = np.zeros_like(time)
    for k in range(n_dis):
        t_start = k * (t_galv + t_rest)
        ind = np.logical_and(t_start <= time, time <= t_start + t_galv)
        iapp[ind] = np.sign(soc0 - socf) * i_galv

    # Run the simulation.
    return sim_profile(cell, time, iapp, soc0, TdegC)


def sim_drive(
    cell: CellParams,
    cycle_name: str,
    soc0: float,
    socf: float,
    TdegC: float = 25.0,
    ts: float = 1.0,
):
    """
    Solve the PDE model for a dynamic drive current profile.

    :param cell: LMB cell model.
    :param cycle_name: Name of the drive cycle ('UDDS', 'LA92', or 'US06')
    :param soc0: Initial cell soc [%].
    :param socf: Final cell soc [%].
    :param TdegC: Cell temperature. DEFAULT 25.
    :param ts: Sampling interval [s]. DEFAULT 1.0.
    :return:
    """

    # Load drive cycle from file.
    data = np.genfromtxt(os.path.join('CSV_DRIVECYCLES', f"{cycle_name}.csv"), delimiter=',')
    time0 = data[:, 0]
    iapp0 = data[:, 1]

    # Interpolate over unform time variable.
    time = np.arange(np.min(time0), np.max(time0) + ts, ts)
    iapp = PchipInterpolator(time0, iapp0)(time)

    # Rescale applied current.
    QdisAh0 = np.trapz(iapp, time)/3600
    QdisAh = (soc0 - socf)/100 * cell.cst.QAh
    iapp = iapp*QdisAh/QdisAh0

    # Run the simulation.
    return sim_profile(cell, time, iapp, soc0, TdegC)


#
# MATLAB io utilities.
#

@dataclass
class FOMOutput:
    ts: float
    time: np.ndarray
    iapp: np.ndarray
    vcell: np.ndarray
    thetass: np.ndarray       # over pos electrode
    thetae: np.ndarray
    xthetass: np.ndarray
    xthetae: np.ndarray


def load_mat_comsol_sim(filepath):
    """
    Load COMSOL simulation output from .mat file.
    """

    simData = sio.loadmat(filepath)
    time = simData['simData']['time'][0, 0].T[0]
    ts = float(np.mean(np.diff(time)))
    iapp = simData['simData']['Iapp'][0, 0].T[0]
    vcell = simData['simData']['Vcell'][0, 0].T[0]
    thetass = simData['simData']['Thetass'][0, 0]
    thetae = simData['simData']['Thetae'][0, 0]
    xthetass = simData['simData']['xLocs'][0, 0]['Thetass'][0, 0].T[0]
    xthetae = simData['simData']['xLocs'][0, 0]['Thetae'][0, 0].T[0]

    return FOMOutput(ts, time, iapp, vcell, thetass, thetae, xthetass, xthetae)


def load_mat_comsol_eis_sim(filepath):
    """
    Load COMSOL EIS simulation output from .mat file.
    """

    simData = sio.loadmat(filepath)
    socPct = simData['SOCs'][0, 0]*100
    freq = simData['freqs'][0]
    Zcell = simData['Z'].T[0]

    return socPct, freq, Zcell

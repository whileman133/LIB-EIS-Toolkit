"""
lmbcell.py

PyBaMM partial-differential-equation (PDE) model for LMB cell.

2024.10.25 | Created | Wesley Hileman <whileman@uccs.edu>
"""

import numpy as np
import pybamm
from base import BaseLumpedModel


class LumpedLMBModel(BaseLumpedModel):
    """
    Model for a lithium-metal battery cell with eff and pos layers.
    Lumped parameters. Dimensionless geometry. MSMR OCP and kinetics.
    Stoichiometry-dependent Ds.

    NOTE 1: For simplicity, we set the potential reference point to the electrolyte at the surface
            of the Li-metal electrode, i.e., we define phi_e(0) = 0V. This makes phi_s_neg(0) nonzero,
            and it must be solved for to compute the cell voltage, thus phi_s_neg(0) appears as
            a variable of the model.

    NOTE 2: The double-layer models break when Rdl(pos)=0 or Rdl(neg)=0. Use small values instead.

    NOTE 3: You must specify the number of MSMR galleries/reactions for the positive electrode
            when constructing the model.
    """

    def __init__(self, J: int, name="Unnamed lumped-parameter lithium-metal cell model"):
        super().__init__(name)

        # Set number of MSMR galleries/reactions to include in the model.
        if J < 1:
            raise ValueError(f"Number of MSMR galleries (J) must be >= 1. Got {J} instead.")
        self.J = J

        # Initialize model.
        self._init_model_parameters()
        self._init_model_variables()
        self._init_model_equations()
        self._init_model_events()

        # Collect model output variables.
        self.variables = {
            'time [s]': pybamm.t,
            'time [min]': pybamm.t / 60,
            'time [h]': pybamm.t / 3600,
            'iapp [A]': self.param.iapp,
            'iapp [C-rate]': self.param.iapp / self.param.Q,
            'soc': self.var.soc,
            'soc [%]': self.var.soc * 100,
            'pos_thetas_avg': self.var.pos_thetas_avg,
            'vcell [V]': self.var.vcell,
            'thetae': self.var.thetae,
            'FNe [A]': self.var.FNe,
            'phie [V]': self.var.phie,
            'ie [A]': self.var.ie,
            'pos_thetas': self.var.pos_thetas,
            'thetass': self.var.thetass,
            'pos_thetass': self.var.pos_thetass,
            'Uss [V]': self.var.Uss,
            'phis [V]': self.var.phis,
            'is [A]': self.var.is_,
            'ifdl [A]': self.var.ifdl,
            'pos_if [A]': self.var.pos_if,
            'phise [V]': self.var.phise,
            'etas [V]': self.var.etas,
            'neg_phis [V]': self.var.neg_phis,
        }

    @property
    def default_geometry(self):
        return {
            "eff": {self.var.x_eff: {"min": 0, "max": 1}},
            "pos": {self.var.x_pos: {"min": 1, "max": 2}},
            "pos_particle": {self.var.r: {"min": 0, "max": 1}},
        }

    @property
    def default_submesh_types(self):
        return {
            "eff": pybamm.Uniform1DSubMesh,
            "pos": pybamm.Uniform1DSubMesh,
            "pos_particle": pybamm.Uniform1DSubMesh,
        }

    @property
    def default_spatial_methods(self):
        return {
            "eff": pybamm.FiniteVolume(),
            "pos": pybamm.FiniteVolume(),
            "pos_particle": pybamm.FiniteVolume(),
        }

    @property
    def default_var_pts(self):
        return {
            self.var.x_eff: 10,
            self.var.x_pos: 20,
            self.var.r: 20
        }

    @property
    def default_solver(self):
        return pybamm.CasadiSolver(root_tol=1e-3, atol=1e-6, rtol=1e-6)

    @property
    def default_quick_plot_variables(self):
        return [
            "iapp [C-rate]", "soc [%]", "vcell [V]",
            "thetae", "FNe [A]", ["ie [A]", "is [A]"],
            "ifdl [A]", "phie [V]", "thetass",
        ]

    @property
    def default_parameter_values(self):
        if self.J != 2:
            raise ValueError(f"To use default parameter values, use J=2 MSMR galleries.")
        return pybamm.ParameterValues({
            "iapp [A]": 0.269177,
            "vmin [V]": 3.5,
            "vmax [V]": 4.2,
            "T [K]": 25+273.15,
            "W": 0.237798,
            "psi [V]": 0.001739,
            "Q [Ah]": 0.269177,
            "pos_sigma [Ohm-1]": 241.896,
            "pos_Dsref [s-1]": 0.00083125,
            "pos_kappa [Ohm-1]": 5.76159,
            "pos_taue [s]": 25.4882,
            "pos_theta0": 0.999996,
            "pos_theta100": 0.10836198,
            "pos_Rf [Ohm]": 0.047348,
            "pos_Rdl [Ohm]": 0.023674,
            "pos_Cdl [F]": 0.023674,
            "pos_Uocp0 [V]": 4.2,
            "pos_U0_0 [V]": 4.16756,
            "pos_U0_1 [V]": 4.02477,
            "pos_X_0": 0.39669,
            "pos_X_1": 0.60331,
            "pos_omega_0": 1.12446,
            "pos_omega_1": 1.71031,
            "pos_k0_0 [A]": 0.99738674,
            "pos_k0_1 [A]": 0.71850266,
            "pos_beta_0": 0.5,
            "pos_beta_1": 0.5,
            "eff_kappa [Ohm-1]": 34.78626,
            "eff_taue [s]": 0.644645,
            "neg_k0 [A]": 1.895276,
            "neg_beta": 0.5,
            "neg_Rdl [Ohm]": 0.000156245,
            "neg_Cdl [F]": 0.0012672,
        })

    def _init_model_parameters(self):
        # Applied current input.
        self.param.iapp = pybamm.FunctionParameter("iapp [A]", inputs={"t": pybamm.t})

        # Constants.
        self.param.F = pybamm.Parameter("Faraday constant [C.mol-1]")
        self.param.R = pybamm.Parameter("Ideal gas constant [J.K-1.mol-1]")

        # Cell-wide parameters.
        self.param.W = pybamm.Parameter("W")
        self.param.psi = pybamm.Parameter("psi [V]")
        self.param.T = pybamm.Parameter("T [K]")
        self.param.Q = pybamm.Parameter("Q [Ah]")
        self.param.vmin = pybamm.Parameter("vmin [V]")
        self.param.vmax = pybamm.Parameter("vmax [V]")

        # Negative electrode.
        self.param.neg_k0 = pybamm.Parameter("neg_k0 [A]")
        self.param.neg_beta = pybamm.Parameter("neg_beta")
        self.param.neg_Rdl = pybamm.Parameter("neg_Rdl [Ohm]")
        self.param.neg_Cdl = pybamm.Parameter("neg_Cdl [F]")

        # Effective layer.
        self.param.eff_kappa = pybamm.Parameter("eff_kappa [Ohm-1]")
        self.param.eff_taue = pybamm.Parameter("eff_taue [s]")

        # Positive electrode.
        self.param.pos_sigma = pybamm.Parameter("pos_sigma [Ohm-1]")
        self.param.pos_Dsref = pybamm.Parameter("pos_Dsref [s-1]")
        self.param.pos_kappa = pybamm.Parameter("pos_kappa [Ohm-1]")
        self.param.pos_taue = pybamm.Parameter("pos_taue [s]")
        self.param.pos_theta0 = pybamm.Parameter("pos_theta0")
        self.param.pos_theta100 = pybamm.Parameter("pos_theta100")
        self.param.pos_Rdl = pybamm.Parameter("pos_Rdl [Ohm]")
        self.param.pos_Cdl = pybamm.Parameter("pos_Cdl [F]")
        self.param.pos_Rf = pybamm.Parameter("pos_Rf [Ohm]")
        self.param.pos_U0 = [pybamm.Parameter(f"pos_U0_{j} [V]") for j in range(self.J)]
        self.param.pos_X = [pybamm.Parameter(f"pos_X_{j}") for j in range(self.J)]
        self.param.pos_omega = [pybamm.Parameter(f"pos_omega_{j}") for j in range(self.J)]
        self.param.pos_k0 = [pybamm.Parameter(f"pos_k0_{j} [A]") for j in range(self.J)]
        self.param.pos_beta = [pybamm.Parameter(f"pos_beta_{j}") for j in range(self.J)]
        self.param.pos_Uocp0 = pybamm.Parameter("pos_Uocp0 [V]")  # initial OCP of positive electrode

    def _init_model_variables(self):
        #
        # Spatial variables.
        #

        self.var.r = pybamm.SpatialVariable(
            "r",
            domain=["pos_particle"],
            auxiliary_domains={"secondary": "pos"},
            coord_sys="spherical polar",
        )
        self.var.x_eff = pybamm.SpatialVariable("x_eff", domain=["eff"], coord_sys="cartesian")
        self.var.x_pos = pybamm.SpatialVariable("x_pos", domain=["pos"], coord_sys="cartesian")

        #
        # Electrolyte variables.
        #

        # Salt concentration.
        self.var.eff_thetae = pybamm.Variable("eff_thetae", domain="eff")
        self.var.pos_thetae = pybamm.Variable("pos_thetae", domain="pos")
        self.var.thetae = pybamm.concatenation(self.var.eff_thetae, self.var.pos_thetae)

        # Electrolyte potential.
        self.var.eff_phie = pybamm.Variable("eff_phie [V]", domain="eff")
        self.var.pos_phie = pybamm.Variable("pos_phie [V]", domain="pos")
        self.var.phie = pybamm.concatenation(self.var.eff_phie, self.var.pos_phie)

        #
        # Solid variables.
        #

        # Solid potential.
        self.var.neg_phis = pybamm.Variable("neg_phis [V]")
        self.var.eff_phis = pybamm.PrimaryBroadcast(pybamm.Scalar(np.nan), "eff")  # phis d.n.e. in eff
        self.var.pos_phis = pybamm.Variable("pos_phis [V]", domain="pos")
        self.var.phis = pybamm.concatenation(self.var.eff_phis, self.var.pos_phis)

        # Particle potential.
        # ..note: We recast the solid-diffusion equation from thetas to U using the
        #         Baker-Verbrugge composition-dependent solid-diffusivity model.
        self.var.pos_Us = pybamm.Variable("pos_Us [V]", domain="pos_particle", auxiliary_domains={"secondary": "pos"})
        self.var.eff_Uss = pybamm.PrimaryBroadcast(np.nan, "eff")  # Uss d.n.e. in eff
        self.var.pos_Uss = pybamm.surf(self.var.pos_Us)
        self.var.Uss = pybamm.concatenation(self.var.eff_Uss, self.var.pos_Uss)

        #
        # Interface variables.
        #

        # Faradaic plus double-layer current.
        self.var.eff_ifdl = pybamm.PrimaryBroadcast(0, "eff")  # ifdl=0 in eff
        self.var.pos_ifdl = pybamm.Variable("pos_ifdl [A]", domain="pos")
        self.var.ifdl = pybamm.concatenation(self.var.eff_ifdl, self.var.pos_ifdl)

        # Double-layer capacitor voltage.
        self.var.neg_vdl = pybamm.Variable("neg_vdl [V]")
        self.var.eff_vdl = pybamm.PrimaryBroadcast(np.nan, "eff")  # vdl d.n.e. in eff
        self.var.pos_vdl = pybamm.Variable("pos_vdl [V]", domain="pos")
        self.var.vdl = pybamm.concatenation(self.var.eff_vdl, self.var.pos_vdl)

        #
        # Expressions.
        #

        self.var.f = self.param.F / self.param.R / self.param.T

        # MSMR variables.
        self.var.pos_g = []
        self.var.pos_xs = []
        self.var.pos_xss = []
        self.var.pos_dxdU = []
        for j in range(self.J):
            # Collect MSMR parameters for gallery j.
            U0_j = self.param.pos_U0[j]
            X_j = self.param.pos_X[j]
            omega_j = self.param.pos_omega[j]

            # Compute the function g_j(U) for reaction j.
            self.var.pos_g.append(pybamm.exp((self.var.pos_Us - U0_j) * self.var.f / omega_j))

            # Compute partial stoichiometry of gallery j.
            self.var.pos_xs.append(X_j / (1 + self.var.pos_g[j]))
            self.var.pos_xss.append(pybamm.surf(self.var.pos_xs[j]))  # at surface of particles

            # Compute differential capacity (dxj/dU) for gallery j.
            self.var.pos_dxdU.append(
                -self.var.f * X_j * self.var.pos_g[j] / omega_j / (1 + self.var.pos_g[j]) ** 2
            )

        # Accumulate solid stoichiometry, differential capacity, and faradaic current.
        self.var.pos_thetas = sum(self.var.pos_xs)
        self.var.pos_thetass = pybamm.surf(self.var.pos_thetas)
        self.var.pos_dthetas = sum(self.var.pos_dxdU)

        # Compute differential OCP dU/d(theta_s) by reciprocal property.
        self.var.pos_dU = 1 / self.var.pos_dthetas

        self.var.eff_thetass = pybamm.PrimaryBroadcast(np.nan, "eff")  # thetass d.n.e. in eff
        self.var.thetass = pybamm.concatenation(self.var.eff_thetass, self.var.pos_thetass)

    def _init_model_equations(self):

        #
        # Current-overpotential equations.
        #

        # Overpotential (pos).
        self.var.eff_phise = pybamm.PrimaryBroadcast(np.nan, "eff")
        self.var.pos_phise = self.var.pos_phis - self.var.pos_phie
        self.var.phise = pybamm.concatenation(self.var.eff_phise, self.var.pos_phise)
        self.var.eff_etas = pybamm.PrimaryBroadcast(np.nan, "eff")
        self.var.pos_etas = self.var.pos_phise - self.var.pos_Uss - self.param.pos_Rf * self.var.pos_ifdl
        self.var.etas = pybamm.concatenation(self.var.eff_etas, self.var.pos_etas)

        # MSMR faradaic current (pos).
        self.var.pos_i0j = []  # gallery exchange currents
        self.var.pos_ifj = []  # gallery faradaic currents
        for j in range(self.J):
            # Collect MSMR parameters for gallery j.
            U0_j = self.param.pos_U0[j]
            X_j = self.param.pos_X[j]
            omega_j = self.param.pos_omega[j]
            k0_j = self.param.pos_k0[j]
            beta_j = self.param.pos_beta[j]

            # Compute exchange current for gallery j.
            # Standard-form exchange current expression:
            # pos_i0.append(
            #     k0_j / (X_j/2)**omega_j *
            #     pos_xss[j]**(omega_j*beta_j) *
            #     (X_j - pos_xss[j])**(omega_j*(1-beta_j)) *
            #     theta_e_pos**(1 - beta_j)
            # )
            # To avoid the singularity when xj->Xj, we reformulate in terms
            # of the potential. (To derive this equation, substitute xj=Xj/(1+gj) into
            # the original exchange-current expression.)
            self.var.pos_i0j.append(
                k0_j / (X_j / 2) ** omega_j *
                self.var.pos_xss[j] ** omega_j *
                pybamm.exp((1 - beta_j) * (self.var.pos_Uss - U0_j) * self.var.f) *
                self.var.pos_thetae ** (1 - beta_j)
            )

            # Compute faradaic current for gallery j.
            self.var.pos_ifj.append(self.var.pos_i0j[j] * (
                pybamm.exp((1 - beta_j) * self.var.f * self.var.pos_etas)
                - pybamm.exp(-beta_j * self.var.f * self.var.pos_etas)
            ))

        # Total faradaic current (pos).
        self.var.pos_if = sum(self.var.pos_ifj)

        # Double-layer current (pos).
        self.var.pos_idl = (self.var.pos_etas + self.var.pos_Uss - self.var.pos_vdl) / self.param.pos_Rdl

        # Current-overpotential equations (pos).
        self.algebraic[self.var.pos_ifdl] = self.var.pos_if + self.var.pos_idl - self.var.pos_ifdl
        self.rhs[self.var.pos_vdl] = (
            (self.var.pos_etas + self.var.pos_Uss - self.var.pos_vdl)
            / self.param.pos_Rdl / self.param.pos_Cdl
        )
        self.initial_conditions[self.var.pos_ifdl] = 0
        self.initial_conditions[self.var.pos_vdl] = self.param.pos_Uocp0

        # Negative electrode equations.
        # ..note: neg_etas = neg_phis since we define phie(0)=0.
        self.var.neg_idl = (self.var.neg_phis - self.var.neg_vdl) / self.param.neg_Rdl
        self.var.neg_i0 = self.param.neg_k0 * pybamm.boundary_value(self.var.thetae, "left") ** (1 - self.param.neg_beta)
        self.var.neg_if = self.var.neg_i0 * (
            pybamm.exp((1 - self.param.neg_beta) * self.var.f * self.var.neg_phis)
            - pybamm.exp(-self.param.neg_beta * self.var.f * self.var.neg_phis)
        )
        self.algebraic[self.var.neg_phis] = self.var.neg_if + self.var.neg_idl - self.param.iapp
        self.rhs[self.var.neg_vdl] = self.var.neg_idl / self.param.neg_Cdl
        self.initial_conditions[self.var.neg_phis] = 0
        self.initial_conditions[self.var.neg_vdl] = 0

        #
        # Charge conservation in solid.
        #

        self.var.eff_is = pybamm.PrimaryBroadcast(np.nan, "eff")  # is d.n.e. in eff
        self.var.pos_is = -self.param.pos_sigma * pybamm.grad(self.var.pos_phis)
        self.var.is_ = pybamm.concatenation(self.var.eff_is, self.var.pos_is)
        self.algebraic[self.var.pos_phis] = pybamm.div(self.var.pos_is) + self.var.pos_ifdl
        self.boundary_conditions[self.var.pos_phis] = {
            "left":  (pybamm.Scalar(0), "Neumann"),
            "right": (-self.param.iapp / self.param.pos_sigma, "Neumann")
        }
        self.initial_conditions[self.var.pos_phis] = self.param.pos_Uocp0

        #
        # Mass conservation in solid.
        #

        pos_D_coeff = self.var.f * self.param.pos_Dsref * self.var.pos_thetas * (1 - self.var.pos_thetas)
        self.rhs[self.var.pos_Us] = -self.var.pos_dU * pybamm.div(pos_D_coeff * pybamm.grad(self.var.pos_Us))
        self.boundary_conditions[self.var.pos_Us] = {
            "left": (pybamm.Scalar(0), "Neumann"),
            "right": (
                self.var.pos_if * abs(self.param.pos_theta100 - self.param.pos_theta0)
                / 10_800 / self.param.Q / pybamm.surf(pos_D_coeff),
                "Neumann"
            )
        }
        self.initial_conditions[self.var.pos_Us] = self.param.pos_Uocp0

        #
        # Charge conservation in electrolyte.
        #

        self.var.eff_ie = (
            -self.param.eff_kappa * pybamm.grad(self.var.eff_phie)
            + self.param.W * self.param.eff_kappa * self.param.psi * self.param.T * pybamm.grad(self.var.eff_thetae)
        )
        self.var.pos_ie = (
            -self.param.pos_kappa * pybamm.grad(self.var.pos_phie)
            + self.param.W * self.param.pos_kappa * self.param.psi * self.param.T * pybamm.grad(self.var.pos_thetae)
        )
        self.var.ie = pybamm.concatenation(self.var.eff_ie, self.var.pos_ie)
        self.algebraic[self.var.eff_phie] = pybamm.div(self.var.eff_ie) - self.var.eff_ifdl
        self.algebraic[self.var.pos_phie] = pybamm.div(self.var.pos_ie) - self.var.pos_ifdl
        self.boundary_conditions[self.var.eff_phie] = {
            "left": (pybamm.Scalar(0), "Dirichlet"),  # set phi_e=0 at surface of neg electrode
            "right": (
                -self.param.iapp / self.param.eff_kappa
                + self.param.W * self.param.psi * self.param.T * pybamm.boundary_gradient(self.var.eff_thetae, "right"),
                "Neumann"
            ),
        }
        self.boundary_conditions[self.var.pos_phie] = {
            "left":  (pybamm.boundary_value(self.var.eff_phie, "right"), "Dirichlet"),
            "right": (pybamm.Scalar(0), "Neumann")
        }
        self.initial_conditions[self.var.eff_phie] = 0
        self.initial_conditions[self.var.pos_phie] = 0

        #
        # Mass conservation in electrolyte.
        #

        FNLi_eff = -self.param.eff_kappa * self.param.psi * self.param.T * pybamm.grad(self.var.eff_thetae)
        FNLi_pos = -self.param.pos_kappa * self.param.psi * self.param.T * pybamm.grad(self.var.pos_thetae)
        self.rhs[self.var.eff_thetae] = (
            (self.var.eff_ifdl - pybamm.div(FNLi_eff))
            / self.param.eff_taue / self.param.eff_kappa / self.param.psi / self.param.T
        )
        self.rhs[self.var.pos_thetae] = (
            (self.var.pos_ifdl - pybamm.div(FNLi_pos))
            / self.param.pos_taue / self.param.pos_kappa / self.param.psi / self.param.T
        )
        self.boundary_conditions[self.var.eff_thetae] = {
            # Zero salt flux:
            "left": (-self.param.iapp / self.param.eff_kappa / self.param.psi / self.param.T, "Neumann"),
            # Concentration continuity at eff-pos boundary:
            "right": (pybamm.boundary_value(self.var.pos_thetae, "left"), "Dirichlet"),
        }
        self.boundary_conditions[self.var.pos_thetae] = {
            # Flux continuity at eff-pos boundary:
            "left": (
                pybamm.boundary_gradient(self.var.eff_thetae, "right")
                * self.param.eff_kappa / self.param.pos_kappa,
                "Neumann"
            ),
            # Zero salt flux:
            "right": (pybamm.Scalar(0), "Neumann")
        }
        self.initial_conditions[self.var.eff_thetae] = 1.0
        self.initial_conditions[self.var.pos_thetae] = 1.0

        # Compute salt flux.
        self.var.eff_FNe = -(
            self.param.eff_kappa * self.param.psi * self.param.T * pybamm.grad(self.var.eff_thetae)
            + self.var.eff_ie
        )
        self.var.pos_FNe = -(
            self.param.pos_kappa * self.param.psi * self.param.T * pybamm.grad(self.var.pos_thetae)
            + self.var.pos_ie
        )
        self.var.FNe = pybamm.concatenation(self.var.eff_FNe, self.var.pos_FNe)

        #
        # Cell-level quantities.
        #

        # Cell voltage.
        self.var.vcell = pybamm.boundary_value(self.var.pos_phis, "right") - self.var.neg_phis

        # Average stoichiometry (pos).
        self.var.pos_thetas_avg = pybamm.Integral(
            pybamm.Integral(
                self.var.pos_thetas,
                self.var.r
            )/(4*np.pi/3),  # average by volume of unit sphere
            self.var.x_pos
        )

        # Cell state of charge.
        self.var.soc = (
            (self.var.pos_thetas_avg - self.param.pos_theta0)
            / (self.param.pos_theta100 - self.param.pos_theta0)
        )

    def _init_model_events(self):
        self.events += [
            pybamm.Event("Minimum voltage", self.var.vcell - self.param.vmin),
            pybamm.Event("Maximum voltage", self.param.vmax - self.var.vcell),
        ]

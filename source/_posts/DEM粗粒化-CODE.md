---
title: 【代码】DEM 的颗粒流数据粗粒化
date: 2024-08-10 11:01:48
tags: [离散单元法, 粗粒化]
categories:
- [离散单元法]
---
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex/dist/katex.min.css">

基于[这篇](/2024/07/30/DEM粗粒化/)的内容，我用 Python 写了一份对**单种颗粒系统**进行粗粒化的代码。

其具体使用流程为：
1. `obj = CG_mono_3D(...)`: 初始化。
2. `obj.set_domain(...)`: 设置矩形的计算域。
3. `obj.get_CG_data(R, CG_width)`: 对 `R` 点进行粗粒化，其粗粒化的长度为 `CG_width`。

<details>
<summary>DEM_CG_mono.py</summary>

```python
"""
Function to convert DEM simulation data into coarse graining field.

NOTE:
    This function is used for graunlar system which consists of single type of particles.
"""
import typing as _T

import numpy as np
import numba as nb
import numpy.typing as npt

from tqdm import tqdm
from scipy.special import erf
from scipy.integrate import simpson


@nb.jit(forceobj=True)
def CG_Gaussian(r: npt.NDArray[np.double], w: float) -> float:
    """Gaussian Coarse Graining function.

    Parameters
    ====
    r: numpy.NDArray[double]
        Relative position vector, whose shape is `(3,)`.
    w: float
        Course grain width `w`.
                                        
    Return
    ===
    W: float
        Average weight.
    """
    c = 3 * w    # Coarse graining cut-off width
    c_2 = c**2

    # Norm of vector `r`
    r_norm: float = r[0]**2 + r[1]**2 + r[2]**2
    # define weight
    W = 0.0
    V_w = 0.0
    if r_norm <= c_2:
        V_w = np.sqrt(8) * np.pi**1.5 * w**3 * erf(c/w/np.sqrt(2)) - \
            4 * c * w**2 * np.pi * np.exp(-c**2 / w**2 / 2)
        W = np.exp(-r_norm / (2 * w**2)) / V_w
    return W


@nb.jit("f8[:, :](f8, f8, f8[:], f8[:])",
        nogil=True, nopython=True, parallel=True)
def SIGMA_q_k(mass_i: float,
              W_i: float,
              vel_i: npt.NDArray[np.double],
              U: npt.NDArray[np.double]) -> npt.NDArray[np.double]:
    """Calculate part of `sigma_qk` contributed by ball `i`

    Parameters
    ---
    mass_i: float
        Mass of ball i.
    W_i: float
        Weight coefficient of ball i.
    vel_i: numpy.NDArray[double]
        Velocity vector of ball i, whose shape is (3,).
    U: numpy.NDArray[double]
        Coarse grained velocity vector, whose shape is (3,).

    Return
    ---
    res: numpy.NDArray[double]
        3x3 Matrix of `sigma_qk`caused by ball i.
    """
    res = np.zeros((3,3), dtype=np.double)
    for i in nb.prange(3):
        for j in nb.prange(3):
            res[i, j] += mass_i * W_i * \
                (vel_i[i] - U[i]) * (vel_i[j] - U[j])
    return res


@nb.jit("f8[:, :](f8[:], f8[:], f8)",
        nogil=True, nopython=True, parallel=True)
def SIGMA_q_c(F_ij: npt.NDArray[np.double],
              l_ij: npt.NDArray[np.double],
              W_ij: float) -> npt.NDArray[np.double]:
    """Calculate part of `sigma_qk` contributed by contact from i to j.

    Parameters
    ---
    F_ij: numpy.NDArray[double]
        Contact force vector, pointing from i to j.
    l_ij: numpy.NDArray[double]
        Contact branch vector, pointing from i to j.
    W_ij: float
        Integrated weight coefficient, from i to j.
    
    Return
    ---
    res: numpy.NDArray[double]
        3x3 Matrix of `sigma_qc`caused by contact from i to j.
    """
    res = np.zeros((3,3), dtype=np.double)
    for i in nb.prange(3):
        for j in nb.prange(3):
            res[i,j] = F_ij[i] * l_ij[j] * W_ij
    return res


class CG_mono_3D:
    def __init__(self,
                 rho: float, 
                 pos: npt.NDArray[np.double],
                 vel: npt.NDArray[np.double],
                 radius: npt.NDArray[np.double],
                 ball_id: npt.NDArray[np.int_],
                 contact_force: npt.NDArray[np.double],
                 contact_branch: npt.NDArray[np.double],
                 contact_end_ball: npt.NDArray[np.double],
                 in_flow: _T.Optional[npt.NDArray[np.bool_]] = None,
                 dp: _T.Optional[float] = None, *,
                 progress_bar: bool = False
                ) -> None:
        """
        Parameters
        ----
        rho: float
            Density of all particles.
        pos: numpy.NDArray[double]
            Position vector, whose shape is `[particle_num, 3]`.
        vel: numpy.NDArray[double]
            Velocity vector, whose shape is `[particle_num, 3]`.
        radius: numpy.NDArray[double]
            Radius list, whose shape is `(particle_num,)`.
        ball_id: numpy.NDArray[int]
            ID number of all balls.
        contact_force: numpy.NDArray[double]
            Contact force list, whose shape is `[contact_num, 3]`.
        contact_branch: numpy.NDArray[double]
            List of contact branch, whose shape is `[contact_num, 3]`.
        contact_end_ball: numpy.NDArray[double]
            List of contacted balls tuple, whose shape is `[contact_num, 2]`.
        in_flow: numpy.NDArray[np.bool_], optional
            Mask array of ball-in-flow, whose shape is `(particle_num,)`.
        dp: float, optional
            Mean diameter of all particles.
        progress_bar: bool, optional, keyword argument only.
            Switch to control whether show the process. Default is False.

        Note
        ----
        1. Take the `i`-th ball for example, its ball ID is `ball_id[i]`.
           And its radius and velocity are `radius[i]` and `vel[i,:]`, 
           respectively. Note that `ball_id[i]` is not equal to  `i`.
        2. For the `j`-th contact, The ID of two contacted balls are 
           `b1 = contact_end_ball[j, 0]`, and `b2 = contact_end_ball[j ,1]`.
        """
        self.particle_num = pos.shape[0]
        self.contact_num  = contact_force.shape[0]

        # => Assign values

        self._id    = ball_id  # Ball id
        self.pos    = pos      # Position vector
        self.vel    = vel      # Velocity vector
        self.radius = radius   # Ball radius

        # Ball density
        # |=> NOTE: Since this function is used for mono-disperse system,
        #    The density of particles are constant.
        self.rho = rho

        # => Contact properties
        self.c_force  = contact_force     # Contact force
        self.c_branch = contact_branch    # Contact branch
        self.c_ends   = contact_end_ball  # Contact end ball tuple

        # => Optional parameters
        # Ball-in-flow mask array.
        self.inflow = in_flow
        # Mean particle diameter
        if dp is None:
            dp = np.mean(radius) * 2
        else:
            self.dp = dp

        # => Mass array
        self.mass = 4/3 * np.pi * np.power(radius, 3) * rho

        # => Ball-ID position array.
        #    For example, `self._c_end_id[j,0]` corresponds to the first end
        #    ball of `j`-th contact (whose ID is `self.c_ends[j,0]`).
        #    And we have `self.c_ends[j,0] = self._id[self._c_end_id[j,0]]`.
        self._c_end_id = [
            [np.where(ball_id == contact_end_ball[i,0])[0][0],
             np.where(ball_id == contact_end_ball[i,1])[0][0]]
            for i in range(self.contact_num)
        ]
        self._c_end_id = np.array(self._c_end_id, dtype=np.intc)

        # => STATUS of domain config
        self._domain_init = False

        # => Max search width
        self._MAX_width = 6

        # => Option to show progress bar
        self._use_tqdm = progress_bar
        pass

    def set_domain(self,
                   xlim: _T.Tuple[float, float],
                   ylim: _T.Tuple[float, float],
                   zlim: _T.Tuple[float, float],
                   periodic_flag: _T.Tuple[bool, bool, bool]):
        """Function to set cuboid domain of coarse graining.
        
        Parameters
        ----
        xlim: (float, float)
            Range of x-coordinate, contains with `(xMin, xMax)`.
        ylim: (float, float)
            Range of y-coordinate, contains with `(yMin, yMax)`.
        zlim: (float, float)
            Range of z-coordinate, contains with `(zMin, zMax)`.
        periodic_flag: (bool, bool, bool)
            Flag tuple to determinate which axis is periodic.
        """
        xlim = (min(*xlim), max(*xlim))
        ylim = (min(*ylim), max(*ylim))
        zlim = (min(*zlim), max(*zlim))
        self._lims = np.array([xlim, ylim, zlim])

        self._isPeriodic  = periodic_flag
        self._hasPeriodic = np.any(periodic_flag)
    
        # Update state
        self._domain_init = True

    def get_CG_data(self,
                    R: npt.NDArray[np.double],
                    w: float):
        """
        Get coarse grainning data at center node `R`.

        Parameter
        ---
        R: numpy.NDArray[double]
            Center of coarse graining node.
        w: float
            Width of coarse graining.
        
        Return
        ---
        rho: float
            Coarse grained density at node `R`.
        U: numpy.NDArray[double]
            Coarse grained velocity at node `R`.
        sigma: numpy.NDArray[double]
            Coarse grained stress tensor at node `R`.
        """
        assert R.shape == (3,)
        if self._domain_init is False:
            raise Exception("Exception: `self.set_domain` haven't run to set domain.")

        rho, U   = self.get_CG_vel(R, w)
        sigma_qk = self.get_CG_sigma_qk(R, w, U)
        sigma_qc = self.get_CG_sigma_qc(R, w)
        sigma    = sigma_qk + sigma_qc
        return rho, U, sigma

    # ============================ #
    # === Sub-process function === #
    # ============================ #

    def get_CG_vel(self,
                   R: npt.NDArray[np.double],
                   CG_width: float) -> _T.Tuple[float, npt.NDArray[np.double]]:
        """Get velocity and density.

        Parameter
        ---
        R: numpy.NDArray[double]
            Center of coarse graining node.
        CG_width: float
            Width of coarse graining.

        Return
        ---
        rho: float
            Coarse grained density at node `R`.
        U: numpy.NDArray[double]
            Coarse grained velocity at node `R`.
        """
        MAX_width = self._MAX_width * CG_width

        # Variable to store result
        rho, Jx, Jy, Jz = 0.0, 0.0, 0.0, 0.0

        # Start loop
        _ball_iter = range(self.particle_num)
        if self._use_tqdm:
            print("|> [Get coarse graining velocity]")
            _ball_iter = tqdm(_ball_iter, unit="Ball")

        for i in _ball_iter:
            # Skip nodes not in flow group
            if (self.inflow is not None) and (not self.inflow[i]):
                continue

            # Get distance vector
            dR = self._get_distance_vec(R, self.pos[i,:])

            # Skip node which too far from center node `R`.
            if np.linalg.norm(dR) >= MAX_width:
                continue

            # Get weight
            W_i = CG_Gaussian(-dR, CG_width)

            # Assign value
            rho += self.mass[i] * W_i
            Jx  += self.mass[i] * W_i * self.vel[i, 0]
            Jy  += self.mass[i] * W_i * self.vel[i, 1]
            Jz  += self.mass[i] * W_i * self.vel[i, 2]
            pass

        # Get velocity
        U = np.array([Jx / rho, Jy / rho, Jz / rho], dtype=np.double)
        return rho, U

    def get_CG_sigma_qk(self,
                        R: npt.NDArray[np.double],
                        CG_width: float,
                        U: npt.NDArray[np.double]) -> npt.NDArray[np.double]:
        """
        Get stress tensor contributed by kinetic energy.
        
        Parameter
        ---
        R: numpy.NDArray[double]
            Center of coarse graining node.
        CG_width: float
            Width of coarse graining.
        U: numpy.NDArray[double]
            Coarse grained velocity at node `R`.
        
        Return
        ---
        sigma_qk: numpy.NDArray[double]
            Stress tensor contributed by kinetic energy.
        
        Note
        ---
        $\sigma_{\\alpha \\beta}^{q(\mathrm{k})} = \sum_{i \in q} m_i v_{i \\alpha}^{'} v_{i \\beta}^{'} W(\\boldsymbol{r} - \\boldsymbol{r}_i (t))$

        where:
        * $v_{i \\alpha}^{'} = u_{\\alpha}^{q}(\\boldsymbol{r}, t) - u_{i,\\alpha}$ is the $\\alpha$ component of the fluctuation velocity of particle i.
        """
        MAX_width = self._MAX_width * CG_width

        # Result array
        sigma_qk = np.zeros((3,3), dtype=np.double)

        # Start loop
        _ball_iter = range(self.particle_num)
        if self._use_tqdm:
            print("|> [Get stress tensor contributed by kinetic energy]")
            _ball_iter = tqdm(_ball_iter, unit="Ball")

        for i in _ball_iter:
            # Skip nodes not in flow group
            if (self.inflow is not None) and (not self.inflow[i]): continue

            # Get distance vector
            dR = self._get_distance_vec(R, self.pos[i,:])

            # Skip node which too far from center node `R`.
            if np.linalg.norm(dR) >= MAX_width: continue

            W_i = CG_Gaussian(-dR, CG_width)
            sigma_qk += SIGMA_q_k(self.mass[i], W_i, self.vel[i, :], U)
            pass

        return sigma_qk

    def get_CG_sigma_qc(self,
                        R: npt.NDArray[np.double],
                        CG_width: float) -> npt.NDArray[np.double]:
        """
        Get stress tensor contributed by contact.
        
        Parameter
        ---
        R: numpy.NDArray[double]
            Center of coarse graining node.
        CG_width: float
            Width of coarse graining.
        
        Return
        ---
        sigma_qc: numpy.NDArray[double]
            Stress tensor contributed by contact.

        Note
        ---
        $\sigma_{\\alpha \\beta}^{q(\mathrm{ck})} = \sum_{i \in q} \sum\limits_{j \in \bar{Q}}^{j \neq i} F_{ij,\\alpha} l_{ij,\\beta} \int_{0}^{1} W(\\boldsymbol{r} - \\boldsymbol{r}_i + s \\boldsymbol{l}_{ij}) \mathrm{d}s$

        where:
        * $v_{i \\alpha}^{'} = u_{\\alpha}^{q}(\\boldsymbol{r}, t) - u_{i,\\alpha}$ is the $\\alpha$ component of the fluctuation velocity of particle i.
        * $\\bar{Q}$ is the union of all particles classes, $q$.
        * $\\boldsymbol{l}_{ij}$ is as the branch vector between two particles, i and j. (from i to j)
        * $\\boldsymbol{F}_{ij}$ is contact force vector between two particles, i and j. (from i to j)
        """
        MAX_width = self._MAX_width * CG_width

        # Result array
        sigma_qc = np.zeros((3,3), dtype=np.double)

        # Start loop
        _c_iter = range(self.contact_num)
        if self._use_tqdm:
            print("|> [Get stress tensor contributed by contact]")
            _c_iter = tqdm(_c_iter, unit="Contact")

        for i in _c_iter:
            # ID of two balls on the i-th contact
            ball_i1 = self._c_end_id[i, 0]
            ball_i2 = self._c_end_id[i, 1]

            # contact force (from i1 to i2)
            F_i12 = self.c_force[i, :]

            # contact branch (from i1 to i2)
            l_i12 = self.c_branch[i, :]

            # Deal with ball_i1
            dR1 = self._get_distance_vec(R, self.pos[ball_i1, :])
            if np.linalg.norm(dR1) < MAX_width:
                W_ij = self._W_ij_4_sigma_qc(dR1, l_i12, CG_width)
                sigma_qc += SIGMA_q_c(F_i12, l_i12, W_ij)

            # Deal with ball_i2
            dR2 = self._get_distance_vec(R, self.pos[ball_i2, :])
            if np.linalg.norm(dR2) < MAX_width:
                W_ij = self._W_ij_4_sigma_qc(dR2, -l_i12, CG_width)
                sigma_qc += SIGMA_q_c(-F_i12, -l_i12, W_ij)

        return sigma_qc

    # ========================== #
    # === Auxiliary function === #
    # ========================== #

    def _W_ij_4_sigma_qc(self, dR, l_ij, CG_width) -> float:
        """
        Parameters
        ---
        dR: numpy.NDArray[double]
            Relative position vector, whose shape is `(3,)`.
        l_ij: numpy.NDArray[double]
            Contact branch vector, whose shape is `(3,)`.
        CG_width: float
            Width of coarse graining.
        
        Return
        ----
        W_ij: float
            Integrated weight coefficient used in the calculation of
            `sigma_qc`.
        """
        _x = np.linspace(0, 1, 101)
        _y = np.array(
            [CG_Gaussian(-dR + s * l_ij, CG_width) for s in _x]
        )
        W_ij = simpson(_y, x=_x)
        return W_ij

    def _get_distance_vec(self,
                          R: npt.NDArray[np.double],
                          r_i: npt.NDArray[np.double]) -> npt.NDArray[np.double]:
        """Get distance between center of coarse graining node `R` and 
        loaction of particle's center `r_i`.

        Parameter
        ---
        R: numpy.NDArray[double]
            Center of coarse graining node.
        r_i: numpy.NDArray[double]
            Loaction of particle's center.
        
        Return
        ---
        dR: numpy.NDArray[double]
            Distance between `R` and `r_i`, with periodic boundary considered.
            Pointing from `R` to `r_i`.
        """
        dR = r_i - R
        if not self._hasPeriodic: return dR

        # Deal with periodic boundary
        for q in range(3):
            if not self._isPeriodic[q]: continue

            # Deal with periodic boundary
            _r_2L = abs(r_i[q] - self._lims[q,0])
            _r_2R = abs(r_i[q] - self._lims[q,1])
            _R_2L = abs(R[q] - self._lims[q,0])
            _R_2R = abs(R[q] - self._lims[q,1])

            r_near_left = _r_2L < _r_2R
            R_near_left = _R_2L < _R_2R

            r2R = abs(dR[q])
            if R_near_left and (not r_near_left):
                tmp = _r_2R + _R_2L
                if tmp < r2R: dR[q] = -tmp
            elif (not R_near_left) and r_near_left:
                tmp = _r_2L + _R_2R
                if tmp < r2R: dR[q] = tmp
            pass
        return dR

    @property
    def progress_bar(self):
        """Bool flag to control whether show the process by tqdm."""
        return self._use_tqdm

    @progress_bar.setter
    def _set_tqdm(self, value: bool):
        self._use_tqdm = bool(value)
    
    @property
    def search_width_coef(self):
        """Coefficient of search radius, and the search radius is
        `self.search_width_coef * CG_width`.
        """
        return self._MAX_width
    
    @search_width_coef.setter
    def _set_MAX_width(self, value: float):
        if value <= 0:
            raise ValueError("`search_width_coef` should be positive")
        elif value <= 3:
            raise ValueError("`search_width_coef` should be larger than cut off width coefficient (3.0)")
        else:
            self._MAX_width = value

```
</details>

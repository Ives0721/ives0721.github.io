---
title: "【笔记】格子Boltzmann方法对热流体的模拟"
date: 2025-01-05 16:20:31
tags: [格子Boltzmann方法, 流体力学]
categories:
- [流体力学]
- [格子Boltzmann方法]
---
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex/dist/katex.min.css">
<!--Katex style sheet for mobile-->
<style type="text/css">
@media only screen and (max-width: 600px) {
  .katex-display > .katex {
    max-width: 100%;
    overflow-x: auto;
    overflow-y: hidden;
  }
}
</style>

> **宇宙安全声明：**  
> **1.这份文档来自我的一份课程作业, 由于笔者水平有限, 所以在内容上写得会比较粗糙, 或者出现一些错误.如有错误, 敬请谅解并指正.**  
> **2. 后续可能随缘更新？【如果又有空去读什么其他文章的话？】**

# Boltzmann方程和格子Boltzmann方法

## Boltzmann方程

Boltzmann方程源自气体动理论, 它通过建模粒子系统的概率分布函数, 描述其集体行为.速度分布函数 $f\left( \mathbf{x},\mathbf{\xi},t \right)$ 的含义通常表示为：在 $t$ 时刻位于 $\mathbf{x}$ 位置的球体微元 $dV = \left\lbrack \mathbf{x,x +}d\mathbf{x} \right\rbrack$ 内, 粒子速度位于区间 $\lbrack\mathbf{\xi},\mathbf{\xi} + d\mathbf{\xi}\rbrack$ 的分子数量为 $f\left( \mathbf{x},\mathbf{\xi},t \right)d\mathbf{\xi}dV$ .Boltzmann方程可以视作是 $f\left( \mathbf{x},\mathbf{\xi},t \right)$ 的雷诺输运方程, 即：

$$
\frac{\mathrm{D}f}{\mathrm{D}t} = \frac{\partial f}{\partial t} + \mathbf{\xi} \cdot \frac{\partial f}{\partial\mathbf{x}} + \mathbf{a} \cdot \frac{\partial f}{\partial\mathbf{\xi}} = \Omega
\tag{1.1}
$$

其中 $\mathbf{\xi =}\frac{\partial\mathbf{x}}{\partial t}$ 为粒子速度, $\mathbf{a =}\frac{\partial\mathbf{\xi}}{\partial t}\mathbf{=}\frac{\mathbf{F}}{m}$ 为粒子的外力加速度,  $\mathbf{F}$ 为粒子所受外力, $m$ 为粒子质量,  $\Omega$ 为描述分子碰撞作用的项.分布函数 $f$ 的速度矩可导出不同宏观量, 即.

$$
\begin{aligned}
    n = \int f d\mathbf{\xi} ,&\, \rho = m \cdot n,\\
    \rho\mathbf{u} = m\int {f \mathbf{\xi}} d\mathbf{\xi} ,&\, \rho E = m\int {\frac{\xi^{2}}{2} f} d\mathbf{\xi}
\end{aligned}
$$

$n$ 为分子数密度, $\rho$ 为流体密度, $\mathbf{u}$ 为流体速度, $E$为总能的密度.

对分子的碰撞进行建模并非易事.为简化式(1.1)的计算, Bhatnagar等提出了BGK碰撞项.这一思路主张将复杂的碰撞视作当前状态 $f$ 向Maxwell平衡态分布 $f^{(eq)}$ 演化的过程, 即：

$$
\Omega = - \frac{1}{\tau_{0}}\left( f - f^{(eq)} \right)
$$

其中 $\tau_{0}$ 为松弛时间.$D$ 维空间下的Maxwell平衡态分布表示为：

$$
f^{(eq)} = \dfrac{n}{\left( 2\pi R_{g}T \right)^{\frac{D}{2}}}\exp\left\lbrack - \frac{\left( \mathbf{\xi} - \mathbf{u} \right)^{2}}{2R_{g}T} \right\rbrack. 
\tag{1.2}
$$

$R_{g}$ 为气体常数,  $T$ 为温度.

在格子Boltzmann方法(Lattice Boltzmann Method, LBM)中, 使用的分布函数是密度分布函数 $mf$ , 下文为方便起见会将密度分布函数记作 $f$ .密度分布函数的演化方程与式(1.1)一致, 其平衡态分布为：

$$
f^{(eq)} = \frac{\rho}{\left( 2\pi R_{g}T \right)^{\frac{D}{2}}}\exp\left\lbrack - \frac{\left( \mathbf{\xi} - \mathbf{u} \right)^{2}}{2R_{g}T} \right\rbrack.
\tag{1.3}
$$

宏观量则表示为：

$$
\rho = \int f d\mathbf{\xi},\quad
\mathbf{j} = \rho\mathbf{u} = \int {f \mathbf{\xi}}d\mathbf{\xi} ,\quad
\rho E = \int_{}^{}{\frac{\xi^{2}}{2}f}d\mathbf{\xi}.
$$


---

# 参考文献

[^1]:  Chen, Y., Ohashi, H., & Akiyama, M. (1994). **Thermal lattice Bhatnagar-Gross-Krook model without nonlinear deviations in macrodynamic equations**. *Phys. Rev. E*, 50, 2776--2783. [DOI: 10.1103/PhysRevE.50.2776](https://link.aps.org/doi/10.1103/PhysRevE.50.2776)

[^2]:  Chen, Y., Ohashi, H., & Akiyama, M. (1997). **Two-Parameter Thermal Lattice BGK Model with a Controllable Prandtl Number**. *Journal of Scientific Computing*, 12(2), 169--185. [DOI:10.1023/A:1025621832215](https://doi.org/10.1023/A:1025621832215)

[^3]: Nikolaos I. Prasianakis, & Konstantinos Boulouchos (2007). **Lattice Boltzmann Method For Simulation Of Weakly Compressible Flows At Arbitrary Prandtl Number**. *International Journal of Modern Physics C*, 18, 602-609. [DOI:10.1142/S012918310701084X](https://www.worldscientific.com/doi/abs/10.1142/S012918310701084X)

[^4]: Hosseini, S. A., Atif, M., Ansumali, S., & Karlin, I. V. (2023). **Entropic lattice Boltzmann methods: A review**. *Computers & Fluids*, 259, 105884. [DOI:10.1016/j.compfluid.2023.105884](https://doi.org/10.1016/j.compfluid.2023.105884)

[^5]:   Anirudh Jonnalagadda; Atul Sharma; Amit Agrawal (2021). **Single relaxation time entropic lattice Boltzmann methods: A developer's perspective for stable and accurate simulations**. *Computers & Fluids*, 215, 104792. [DOI:10.1016/j.compfluid.2020.104792](https://doi.org/10.1016/j.compfluid.2020.104792)

[^6]: Qian, Y.H (1993). **Simulating thermohydrodynamics with lattice BGK models**. *Journal of Scientific Computing*, 8, 231--242. <DOI:10.1007/BF01060932>

[^7]: Alexander, F., Chen, S., & Sterling, J. (1993). **Lattice Boltzmann thermohydrodynamics**. *Physical Review E*, 47, R2249--R2252. <DOI:10.1103/PhysRevE.47.R2249>

[^8]:  A. Bartoloni, Claudia Battista, Simone Cabasino, P. S. Paolucci, J. Pech, Renata Sarno, Gian Marco Todesco, Mario Torelli, Walter Tross, Piero Vicini, Roberto Benzi, Nicola Cabibbo, Federico Massaioli, & Raffaele Tripiccione (1993). **LBE simulations of Rayleigh-Bénard convection on the APE100 parallel processor**. *International Journal of Modern Physics C*, 04, 993-1006. <DOI:10.1142/S012918319300077X>

[^9]: Xiaoyi He, Shiyi Chen, & Gary D. Doolen (1998). **A Novel Thermal Model for the Lattice Boltzmann Method in Incompressible Limit**. *Journal of Computational Physics*, 146(1), 282-300. [DOI:10.1006/jcph.1998.6057](https://doi.org/10.1006/jcph.1998.6057)

[^10]: Guo, Z., Zheng, C., Shi, B., & Zhao, T. (2007). **Thermal lattice Boltzmann equation for low Mach number flows: Decoupling model**. *Physical Review E*, 75, 036704. [DOI:10.1103/PhysRevE.75.036704](https://doi.org/10.1103/PhysRevE.75.036704)

[^11]:  Dipankar Chatterjee (2009). **An enthalpy-based thermal lattice Boltzmann model for non-isothermal systems**. *Europhysics Letters*, 86(1), 14004. [DOI:10.1209/0295-5075/86/14004](https://iopscience.iop.org/article/10.1209/0295-5075/86/14004)

[^12]:  He, X., & Luo, L.-S. (1997). **Theory of the lattice Boltzmann method: From the Boltzmann equation to the lattice Boltzmann equation**. *Physical Review E*, 56(6), 6811--6817. [DOI:10.1103/PhysRevE.56.6811](https://doi.org/10.1103/PhysRevE.56.6811)

[^13]:  杨鲲 & 单肖文. (2020). **多层速度格子Boltzmann方法进展及展望**. *空气动力学学报*, 40(3), 23--45. [DOI:10.7638/kqdlxxb-2021.0348](https://doi.org/10.7638/kqdlxxb-2021.0348)

[^14]: Koelman, J. M. V. A. (1991). **A Simple Lattice Boltzmann Scheme for Navier-Stokes Fluid Flow**. *Europhysics Letters*, 15(6), 603--607. [DOI:10.1209/0295-5075/15/6/007](https://doi.org/10.1209/0295-5075/15/6/007)

[^15]: Qian, Y. H., D'Humières, D., & Lallemand, P. (1992). **Lattice BGK Models for Navier-Stokes Equation**. *Europhysics Letters,* 17(6), 479--484. [DOI:10.1209/0295-5075/17/6/001](https://doi.org/10.1209/0295-5075/17/6/001)

[^16]: He, X., & Luo, L.-S. (1997). **Lattice Boltzmann Model for the Incompressible Navier--Stokes Equation**. *Journal of Statistical Physics*, 88(3/4), 927--944. [DOI:10.1023/B:JOSS.0000015179.12689.e4](https://doi.org/10.1023/B:JOSS.0000015179.12689.e4)

[^17]: Yong, Wen-an., Luo, Li-Shi (2005). **Nonexistence of H Theorem for Some Lattice Boltzmann Models**. *Journal of Statistical Physics*, 121, 91--103. [DOI:10.1007/s10955-005-5958-9](https://doi.org/10.1007/s10955-005-5958-9)

[^18]: Zaki Abiza, Miguel Chavez, David M. Holman, & Ruddy Brionnaud. **Prediction of Finned Projectile Aerodynamics using a Lattice-Boltzmann Method CFD Solution** (Vol:10, No:05, 2016). World Academy of Science, Engineering and Technology.

[^19]: 邵菲, 韩端锋, 刘强, & 谢伟. (2016). **熵格子Boltzmann方法的亚格子尺度模型**. *中国舰船研究*, 11(3), 43--47. [DOI:10.3969/j.issn.1673-3185.2016.03.008](https://doi.org/10.3969/j.issn.1673-3185.2016.03.008)

[^20]: Frapolli, N., Chikatamarla, S., & Karlin, I. (2014). **Multispeed entropic lattice Boltzmann model for thermal flows**. *Physical Review E*, 90, 043306. [DOI:10.1103/PhysRevE.90.043306](https://link.aps.org/doi/10.1103/PhysRevE.90.043306)

[^21]: N.I. Prasianakis, S.S. Chikatamarla, I.V. Karlin, S. Ansumali, & K. Boulouchos (2006). **Entropic lattice Boltzmann method for simulation of thermal flows**. *Mathematics and Computers in Simulation*, 72(2), 179-183. [DOI:10.1016/j.matcom.2006.05.012](https://doi.org/10.1016/j.matcom.2006.05.012)

[^22]: Chen, Z., Shu, C., & Tan, D. (2017). **Three-dimensional simplified and unconditionally stable lattice Boltzmann method for incompressible isothermal and thermal flows**. *Physics of Fluids*, 29(5), 053601. [DOI:10.1063/1.4983339](https://doi.org/10.1063/1.4983339)

[^23]: Y. Wang, C. Shu, & C.J. Teo (2014). **Thermal lattice Boltzmann flux solver and its application for simulation of incompressible thermal flows**. *Computers & Fluids*, 94, 98-111. [DOI: 10.1016/j.compfluid.2014.02.006](https://www.sciencedirect.com/science/article/pii/S0045793014000619)

[^24]: Sharma, K. V., Straka, R., & Tavares, F. W. (2020). **Current status of Lattice Boltzmann Methods applied to aerodynamic, aeroacoustic, and thermal flows**. *Progress in Aerospace Sciences*, 115, 100616. [DOI:10.1016/j.paerosci.2020.100616](https://doi.org/10.1016/j.paerosci.2020.100616)

[^25]: Li, Q., He, Y. L., Wang, Y., & Tao, W. Q. (2007). **Coupled double-distribution-function lattice Boltzmann method for the compressible Navier-Stokes equations**. *Physical Review E*, 76(5), 056705. [DOI:10.1103/PhysRevE.76.056705](https://doi.org/10.1103/PhysRevE.76.056705)

[^26]: Qu, K., Shu, C., & Chew, Y. (2007). **Alternative method to construct equilibrium distribution functions in lattice-Boltzmann method simulation of inviscid compressible flows at high Mach number**. *Physical Review E*, 75, 036706.

[^27]: Q. Li, Y.L. He, Y. Wang, & G.H. Tang (2009). **Three-dimensional non-free-parameter lattice-Boltzmann model and its application to inviscid compressible flows.** *Physics Letters A*, 373(25), 2101-2108. [DOI:10.1016/j.physleta.2009.04.036](https://doi.org/10.1016/j.physleta.2009.04.036)

[^28]:  Kataoka, T., & Tsutahara, M. (2004). **Lattice Boltzmann model for the compressible Navier-Stokes equations with flexible specific-heat ratio**. *Physical Review E*, 69, 035701. [DOI:10.1103/PhysRevE.69.035701](https://link.aps.org/doi/10.1103/PhysRevE.69.035701)

[^29]: Zheng, L., Shi, B., & Guo, Z. (2008). **Multiple-relaxation-time model for the correct thermohydrodynamic equations**. *Physical Review E*, 78, 026705. [DOI:10.1103/PhysRevE.78.026705](https://link.aps.org/doi/10.1103/PhysRevE.78.026705)

[^30]: Ai-Guo Xu, Guang-Cai Zhang, Yan-Biao Gan, Feng Chen & Xi-Jun Yu (2012). **Lattice Boltzmann modeling and simulation of compressible flows**. *Frontiers of Physics*. 7, 582--600. [DOI:10.1007/s11467-012-0269-5](https://doi.org/10.1007/s11467-012-0269-5)

[^31]: Feng Chen, Aiguo Xu, Guangcai Zhang, Yingjun Li, & Sauro Succi (2010). **Multiple-relaxation-time lattice Boltzmann approach to compressible flows with flexible specific-heat ratio and Prandtl number**. *Europhysics Letters*, 90(5), 54003. [DOI:10.1209/0295-5075/90/54003](https://dx.doi.org/10.1209/0295-5075/90/54003)

[^32]: Gan, Y.-B., Xu, A.-G., Zhang, G.-C., & Li, Y.-J. (2011). **Flux Limiter Lattice Boltzmann Scheme Approach to Compressible Flows with Flexible Specific-Heat Ratio and Prandtl Number**. *Communications in Theoretical Physics*, 56(3), 490--498. [DOI:10.1088/0253-6102/56/3/18](https://doi.org/10.1088/0253-6102/56/3/18)


---
title: "LBM笔记：Exact Difference Method 构建的源项"
date: 2024-09-07 21:38:05
tags: [格子Boltzmann方法]
categories:
- [格子Boltzmann方法, 源项]
---
我在阅读其他文献时，看到有文献 Kupershtokh 的一篇关于 LBM 外力项的文章[^Kupershtokh2004]，他提出的格式在其他文章中又被称为 *Exact difference method（EDM）*。

本想搜索这篇文章的PDF进行阅读，但最后只照着一个网站介绍了部分信息。这份粗略笔记的内容，都是基于其他文献[^Kupershtokh2009][^Li2016][^Li_PRE_2016]中对该方法的介绍而整理的。



# 介绍

## Boltzmann-BGK方程的近似

密度分布函数 $f$ 的 Boltzmann-BGK方程写作：

$$
\frac{\partial f}{\partial t} + \bold{\xi} \cdot \nabla f + \bold{a} \cdot \nabla_{\bold{\xi}} f  = -\frac{1}{\tau_f} \left[ f - f^{(eq)}\right]
\tag{1}
$$

其中：

$$
f^{(eq)} = \frac{\rho}{(2 \pi R_g T)^{D/2}} \exp{\left(- \frac{(\bold{\xi} - \bold{u})^2}{2 R_g T}\right)}
\tag{2}
$$

符号含义分别为
* $\bold{\xi}$ : 粒子速度；
* $\rho$ : 流体密度；
* $\bold{u}$ : 流场宏观速度；
* $R_g$ : 气体常数；
* $T$ : 温度；
* $\bold{a}=\mathrm{d}\bold{u}/\mathrm{d}t$ : 流场加速度。

当**Knudsen数较小**时，Boltzmann方程的体力项可以基于$f^{(eq)}$进行近似，即
$$\nabla_{\bold{\xi}} f \approx \nabla_{\bold{\xi}} f^{(eq)}$$

根据式(2)， $f^{(eq)}$ 是分子随机速度 $\bold{C} = \bold{\xi} - \bold{u}$ 的指数函数，则：
$$
\frac{\partial f^{(eq)}}{\partial \bold{u}} = -\frac{\partial f^{(eq)}}{\partial \bold{\xi}} = -\nabla_{\bold{\xi}} f^{(eq)}
$$

若**假设密度 $\rho$ 为常数**，则 $\frac{\mathrm{d} \rho}{\mathrm{d} t}=0$，那么 $f^{(eq)}$ 的全导数为：

$$
\begin{aligned}
\frac{\mathrm{d} f^{(eq)}}{\mathrm{d} t} &= \frac{\partial f^{(eq)}}{\partial \bold{u}} \cdot \frac{\mathrm{d} \bold{u}}{\mathrm{d} t} + \frac{\partial f^{(eq)}}{\partial \rho} \cdot \frac{\mathrm{d} \rho}{\mathrm{d} t} \\
&\approx \frac{\partial f^{(eq)}}{\partial \bold{u}} \cdot \frac{\mathrm{d} \bold{u}}{\mathrm{d} t} = -\bold{a} \cdot \nabla_{\bold{\xi}} f^{(eq)}
\end{aligned}
$$

也就是说，Kupershtokh 等[^Kupershtokh2004][^Kupershtokh2009]近似的 Boltzmann-BGK方程写作：

$$
\frac{\partial f}{\partial t} + \bold{\xi} \cdot \nabla f = -\frac{1}{\tau_f} \left[ f - f^{(eq)}\right] + \frac{\mathrm{d} f^{(eq)}}{\mathrm{d} t}
\tag{3}
$$

## 近似方程的离散化

将式(3)在$[t,t+\delta_t]$内沿着$\bold{e}_\alpha$方向积分，可得

$$\begin{aligned}
f_\alpha (\bold{x} + \bold{e}_\alpha \delta_t, t+\delta_t) - f_\alpha (\bold{x}, t) =& -\frac{1}{\tau} \left[ f_\alpha (\bold{x}, t) - f_\alpha^{(eq)} (\rho, \bold{u^*}) \right] \\
& + \left[ f_\alpha^{(eq)} (\bold{u^*}(t+\delta_t)) - f_\alpha^{(eq)} (\bold{u^*}(t)) \right] 
\end{aligned}
\tag{4}
$$

其中：

- $\delta_t$ 为时间步长；
- $\bold{x}$ 为网格点位置，$\rho, \bold{u^*}$ 分别为该点在 $t$ 时刻的密度和速度（与宏观速度 $\bold{u}$ 有区别）；
- $f_\alpha, f_\alpha^{(eq)}$ 分别为 $\bold{e}_\alpha$ 方向的**分布函数**和**平衡态**；
- $F_\alpha$ 为外力项，其表达式基于外力 $\bold{F} = \rho \dfrac{\mathrm{d} \bold{u^*}}{\mathrm{d} t} = \rho \bold{a}$ 进行构建。

并且

$$
\begin{aligned}
\bold{u^*}(t) &= \left(\sum_{\alpha} f_\alpha \bold{e}_\alpha \right) / \left(\sum_{\alpha} f_\alpha \right) \\
\bold{u^*}(t+\delta_t) &= \bold{u^*}(t) + \int_{t}^{t+\delta_t} \frac{\mathrm{d} \bold{u^*}}{\mathrm{d} t} \mathrm{d} t \\
&= \bold{u^*}(t) + \int_{t}^{t+\delta_t} \frac{\bold{F}}{\rho} \mathrm{d} t 
\approx \bold{u^*}(t) + \frac{\bold{F}}{\rho} \delta_t
\end{aligned}
\tag{5}
$$

这意味着加速度 $\bold{a}=\bold{F}/\rho$ 在 $[t,t+\delta_t]$ 时间段内被视为常数。

Kupershtokh等在这篇文章[^Kupershtokh2009]说明，在使用该源项时，流体宏观量的计算应表示为：

$$
\rho = \sum_{\alpha} f_\alpha ,\quad \rho \bold{u} = \sum_{i} \bold{e}_\alpha f_\alpha + \frac{\delta_t \bold{F}}{2}
\tag{6}
$$

Li 等[^Li_PRE_2016]通过 Chapman-Enskog 展开证明：该格式向宏观方程引入的误差项为：
$$
\mathbf{Error}_{\mathrm{EDM}} = -\delta_t^2 \nabla\cdot\left( \frac{\bold{FF}}{\rho} \right)
\tag{7}
$$


[^Kupershtokh2004]: Kupershtokh, A. L. [**New method of incorporating a body force term into the lattice Boltzmann equation**](https://www.elibrary.ru/item.asp?id=28981868). in *Proceeding of the 5th international EHD workshop* 241–246 (Poitiers, France, 2004).
[^Kupershtokh2009]: A.L. Kupershtokh, D.A. Medvedev, & D.I. Karpov (2009). **On equations of state in a lattice Boltzmann method**. Computers & Mathematics with Applications, 58(5), 965-974. DOI:10.1016/j.camwa.2009.02.024.
[^Li2016]: Q. Li, K.H. Luo, Q.J. Kang, Y.L. He, Q. Chen, & Q. Liu (2016). **Lattice Boltzmann methods for multiphase flow and phase-change heat transfer**. Progress in Energy and Combustion Science, 52, 62-105. DOI:10.1016/j.pecs.2015.10.001.
[^Li_PRE_2016]: Li, Q., Zhou, P., & Yan, H. (2016). **Revised Chapman-Enskog analysis for a class of forcing schemes in the lattice Boltzmann method**. Physical Review E, 94, 043313. DOI:10.1103/PhysRevE.94.043313.
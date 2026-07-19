---
title: "【笔记】双松弛时间格子 Boltzmann 方法"
date: 2026-07-14 17:28:40
tags: [格子Boltzmann方法]
categories:
- [格子Boltzmann方法]
---
<script defer src="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.js" integrity="sha512-LQNxIMR5rXv7o+b1l8+N1EZMfhG7iFZ9HhnbJkTp4zjNr5Wvst75AqUeFDxeRUa7l5vEDyUiAip//r+EFLLCyA==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>
<script defer src="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/contrib/copy-tex.min.js" integrity="sha512-cQxSkTM4RvFAjdBeBDkrllhYfERwZWjM/oijKfPsmhR0JI2U3fXSlbUaJp5SlgBh/FzYmMyWudRTXLnNj3MbEA==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>
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
.container {
    display: flex;
}
.col, p {
    flex: 1;
    font-family: "Arial";
    text-align: left;
}
</style>

双松弛时间（Two-relaxation-time）格子 Boltzmann 方法是多松弛时间（Multi-Relaxation-Time）LBM 的一种特殊形式。它仅具有两个松弛时间，分为对称和非对称两部分，同时保持了与 LBGK 模型一致的简单算法和计算效率。

# Concepts

## 基本方程

对于 $t$ 时刻的点 $\boldsymbol{x}$ 沿 $\boldsymbol{c}_i$ 方向的分布函数 $f_i(\boldsymbol{x}, t)$，其在 $\Delta t$ 时间步长内的演化为：

$$
\begin{aligned}
f_i(\boldsymbol{x}+\boldsymbol{c}_i \Delta t, t + \Delta t) - f_i(\boldsymbol{x}, t) &= 
-\frac{1}{\tau^{+}} \left[ f_i^{+} - f_i^{(\mathrm{eq})+} \right] \\
&\quad -\frac{1}{\tau^{-}} \left[ f_i^{-} - f_i^{(\mathrm{eq})-} \right] \\
&\quad + \mathbb{F}_i \Delta t
\end{aligned}
\tag{1}
$$

其中：

* $\tau^{+}$ 和 $\tau^{-}$ 分别为两部分的松弛时间（均为正数），其数值由流体运动粘度 $\nu$ 和参数 $\Lambda = (\tau^{+} - \frac{1}{2}) (\tau^{-} - \frac{1}{2})$ 共同决定。d’Humières 等[^Humières2009]提出 $\Lambda$ 控制了 TRT 方程的稳态场。并且得出 $\Lambda=1/4$ 可使 TRT 的运算保持稳定。

$$
\tau^{+} = \frac{\nu}{c_s^2} + \frac{1}{2} ,\quad
\tau^{-} = \frac{2 \Lambda}{2 \tau^{+} - 1} + \frac{1}{2}
$$



* $f_{-i}(\boldsymbol{x}, t)$ 和 $f_{-i}^{(\mathrm{eq})}(\boldsymbol{x}, t)$ 为 $\boldsymbol{c}_{-i}=-\boldsymbol{c}_{i}$ 方向的分布函数和平衡态分布。

* $f_i^{+}$ 和 $f_i^{-}$ 分别为分布函数的**对称**部分和**非对称**部分，即

$$
f_i^{+} = \frac{1}{2}\left( f_i(\boldsymbol{x}, t) + f_{-i}(\boldsymbol{x}, t) \right) ,\quad
f_i^{-} = \frac{1}{2}\left( f_i(\boldsymbol{x}, t) - f_{-i}(\boldsymbol{x}, t) \right)
$$

* 同理， $f_i^{(\mathrm{eq})+}$ 和 $f_i^{(\mathrm{eq})-}$ 分别为平衡态分布 $f_i^{(\mathrm{eq})}$ 的对称部分和非对称部分，即

$$
\begin{aligned}
f^{\mathrm{eq}}_i &= w_i \rho \left[ 1 + \frac{\boldsymbol{c}_i \cdot \boldsymbol{u}}{c_s^2} + \frac{(\boldsymbol{c}_i \cdot \boldsymbol{u})^2}{2 c_s^4} - \frac{u^2}{2 c_s^2} \right] ,\\
f_i^{(\mathrm{eq})+} &= \frac{1}{2}\left( f_i^{(\mathrm{eq})}(\boldsymbol{x}, t) + f_{-i}^{(\mathrm{eq})}(\boldsymbol{x}, t) \right) =w_i \rho \left[ 1 + \frac{(\boldsymbol{c}_i \cdot \boldsymbol{u})^2}{2 c_s^4} - \frac{u^2}{2 c_s^2} \right] ,\\
f_i^{(\mathrm{eq})-} &= \frac{1}{2}\left( f_i^{(\mathrm{eq})}(\boldsymbol{x}, t) - f_{-i}^{(\mathrm{eq})}(\boldsymbol{x}, t) \right) = w_i \rho \frac{\boldsymbol{c}_i \cdot \boldsymbol{u}}{c_s^2} 
\end{aligned}
$$

* $w_i$ 为 $\boldsymbol{c}_{i}$ 方向的权重， $\rho$ 和 $\boldsymbol{u}$ 为宏观密度和速度， $c_s$ 为DnQb模型的格子声速。

* $\mathbb{F}_i$ 为源项 $\boldsymbol{F}$ 在 $i$ 方向上的函数。

## TRT 的外力项

> 下面记 $t_i = {w_i} / {c_s^2}$，以及 $F_{i}^{+} = \frac{1}{2} (\mathbb{F}_i + \mathbb{F}_{-i})$ 和 $F_{i}^{-} = \frac{1}{2} (\mathbb{F}_i - \mathbb{F}_{-i})$ 。

Postma 等[^Postma2020]的研究对几种主流的源项模型进行了研究。不过他原文里列出的源项方程我还有些看不懂，我就照着[Bawazeer等的文章](https://doi.org/10.1016/j.jocs.2023.102034)抄了两个简单的源项：

<!--
// https://www.sciencedirect.com/science/article/pii/S1290072923006129
// https://www.sciencedirect.com/science/article/pii/S0045793020301092
-->

（1）Buick and Greated [^Buick2000]:

$$
\mathbb{F}_i = \left(1 - \frac{1}{2 \tau} \right) t_q \boldsymbol{c}_{i} \cdot
\boldsymbol{F}
$$

转写为正负两项后表示为：

$$
F_{i}^{+} = 0, \quad
F_{i}^{-} = \left(1 - \frac{1}{2 \tau^{-}} \right) t_i \boldsymbol{c}_{i} \cdot
\boldsymbol{F}
$$

（2）Guo et al. [^Guo2002]:

$$
\mathbb{F}_{i} = \left(1 - \frac{1}{2 \tau} \right) w_i \cdot 
\left[\frac{\boldsymbol{c}_i - \boldsymbol{u}}{c_s^2} + \frac{(\boldsymbol{c}_i \cdot \boldsymbol{u}) \boldsymbol{c}_i}{c_s^4}\right] \cdot \boldsymbol{F}
$$

转写为正负两项后表示为：
$$
\begin{aligned}
    F_i^{+} =& \left(1 - \frac{1}{2 \tau^{+}} \right) t_i \cdot \left[-\boldsymbol{u} + \frac{(\boldsymbol{c}_i \cdot \boldsymbol{u}) \boldsymbol{c}_i}{c_s^2}\right] \cdot \boldsymbol{F} \\
    F_i^{-} =& \left(1 - \frac{1}{2 \tau^{-}} \right) t_i \boldsymbol{c}_i \cdot \boldsymbol{F} \\
\end{aligned}
$$ 

宏观密度 $\rho$ 和动量 $\boldsymbol{j}$ 的计算均表示为：

$$
\rho = \sum\limits_{i=0}^{b-1} f_i ,\quad
\boldsymbol{j} = \rho \boldsymbol{u} = \sum\limits_{i=0}^{b-1} f_i \boldsymbol{c}_{i} + \frac{\boldsymbol{F}}{2}
$$

## 单步TRT-LBM的计算流程

(1) 计算 $f_i^{+}$、 $f_i^{-}$ 以及 $f_i^{(\mathrm{eq})+}$、 $f_i^{(\mathrm{eq})-}$。 

(2) 计算碰撞步  
$$
f_i^{(\mathrm{pc})}(\boldsymbol{x}, t) = f_i(\boldsymbol{x}, t)  
-\frac{1}{\tau^{+}} \left[ f_i^{+} - f_i^{(\mathrm{eq})+} \right]
-\frac{1}{\tau^{-}} \left[ f_i^{-} - f_i^{(\mathrm{eq})-} \right] + \mathbb{F}_i \Delta t
\tag{2}
$$

(3) 计算流动步（与常规LBM一致）  

(4) 更新宏观量。  

# TRT-LBM 与其他 LBM 的联系

## TRT 与单松弛 LBGK

将式(1)整理，可得式(3)：
$$
\begin{aligned}
f_i(\boldsymbol{x}+\boldsymbol{c}_i \Delta t, t + \Delta t) - f_i(\boldsymbol{x}, t) &= - \left(\frac{1}{2 \tau^{+}} + \frac{1}{2 \tau^{-}}\right) [f_i(\boldsymbol{x}, t) - f^{(\mathrm{eq})}_i] \\
&\quad -\left(\frac{1}{2 \tau^{+}} - \frac{1}{2 \tau^{-}}\right) [f_{-i}(\boldsymbol{x}, t) - f^{(\mathrm{eq})}_{-i}] \\
&\quad + \mathbb{F}_i \Delta t
\end{aligned}
\tag{3}
$$

所以，当 $\tau^{+} = \tau^{-}$ 时，TRT 方程可以退化为单松弛 LBGK 方程。此时： $\Lambda = \frac{\nu^2}{c_s^4}$，$\tau^{+}$ 也回归到单松弛 LBGK 方程的定义。

## TRT 与 MRT

下文为方便起见，记 $\mathtt{w}_1 = (\frac{1}{\tau^{+}} + \frac{1}{\tau^{-}})/2$ 和 $\mathtt{w}_2 = (\frac{1}{\tau^{+}} - \frac{1}{\tau^{-}})/2$ ，并以D2Q9模型为例。

先定义：
* $\boldsymbol{f} = [f_0(\boldsymbol{x}, t), ..., f_8(\boldsymbol{x}, t)]^\mathrm{T}$ 分布函数向量；
* $\boldsymbol{f}^{(\mathrm{pc})} = [f_0^{(\mathrm{pc})}(\boldsymbol{x}, t), ..., f_8^{(\mathrm{pc})}(\boldsymbol{x}, t)]^\mathrm{T}$ 碰撞后的分布函数向量；  
* $\boldsymbol{f}^{(\mathrm{eq})} = [f_0^{(\mathrm{eq})}(\boldsymbol{x}, t), ..., f_8^{(\mathrm{eq})}(\boldsymbol{x}, t)]^\mathrm{T}$ 平衡态分布函数向量。   

则 TRT 方程的碰撞步在D2Q9模型下的矩阵形式为：
$$
\boldsymbol{f}^{(\mathrm{pc})} = [\bold{H}] (\boldsymbol{f} - \boldsymbol{f}^{(\mathrm{eq})})
\tag{4}
$$

其中矩阵 $[\bold{H}]$ 为：

$$
[\bold{H}] = \begin{bmatrix}
\mathtt{w}_1+\mathtt{w}_2 &&&&&&&\\
& \mathtt{w}_1 & & \mathtt{w}_2 &&&&&\\
&& \mathtt{w}_1 & & \mathtt{w}_2 &&&&\\
& \mathtt{w}_2 & & \mathtt{w}_1 &&&&&\\
&& \mathtt{w}_2 & & \mathtt{w}_1 &&&&\\
&&&&& \mathtt{w}_1 & & \mathtt{w}_2 &\\
&&&&&& \mathtt{w}_1 & & \mathtt{w}_2 \\
&&&&& \mathtt{w}_2 & & \mathtt{w}_1 &\\
&&&&&& \mathtt{w}_2 & & \mathtt{w}_1
\end{bmatrix}
$$

对比 MRT 方程的碰撞步在D2Q9模型下的矩阵形式：
$$
\boldsymbol{f}^{(\mathrm{pc})} = [\bold{M}^{-1}][\bold{S}][\bold{M}] (\boldsymbol{f} - \boldsymbol{f}^{(\mathrm{eq})})
\tag{5}
$$
其中 $[\bold{S}]$ 为松弛系数的对角矩阵。变换矩阵 $[\bold{M}]$ 为：

$$
[\bold{M}] = \begin{bmatrix}
1&1&1&1&1&1&1&1&1\\
-4&-1&-1&-1&-1&2&2&2&2\\
4&-2&-2&-2&-2&1&1&1&1\\
0&1&0&-1&0&1&-1&-1&1\\
0&-2&0&2&0&1&-1&-1&1\\
0&0&1&0&-1&1&1&-1&-1\\
0&0&-2&0&2&1&1&-1&-1\\
0&1&-1&1&-1&0&0&0&0\\
0&0&0&0&0&1&-1&1&-1
\end{bmatrix}
$$

求解方程 $[\bold{H}] = [\bold{M}^{-1}][\bold{S}][\bold{M}]$，可得：

$$
\begin{aligned}
[\bold{S}] =& [\bold{M}][\bold{H}][\bold{M}^{-1}] \\
=&\mathrm{diag}(
    \mathtt{w}_1+\mathtt{w}_2,
    \mathtt{w}_1+\mathtt{w}_2,
    \mathtt{w}_1+\mathtt{w}_2,\\&\quad
    \mathtt{w}_1-\mathtt{w}_2,
    \mathtt{w}_1-\mathtt{w}_2,
    \mathtt{w}_1-\mathtt{w}_2,\\&\quad
    \mathtt{w}_1-\mathtt{w}_2,
    \mathtt{w}_1+\mathtt{w}_2,
    \mathtt{w}_1+\mathtt{w}_2
)
\end{aligned}
$$

可见 TRT 模型就是 MRT 模型的特例，将 MRT 中所使用的诸多松弛系数减少至到 2 个（$\frac{1}{\tau^{+}}$ 和 $\frac{1}{\tau^{-}}$）。

<!--
1. Ginzburg, I., d’Humières, D. & Kuzmin, A. Optimal Stability of Advection-Diffusion Lattice Boltzmann Models with Two Relaxation Times for Positive/Negative Equilibrium. J Stat Phys 139, 1090–1143 (2010). https://doi.org/10.1007/s10955-010-9969-9
-->

# 实例代码

为了练手，我写了一个双周期剪切流的代码。

- 计算域： $[0,L] \times [0,L]$，其中 $L=2\pi$。
- 网格大小： $\delta_x = L/N$ ($N=128$)。
- 雷诺数： $Re = 1 / \nu = 100$。 其中 $\nu$ 为流体运动粘度 (Physical)。
- 时间步 $\delta_t$ ： 用 LBM 中 $\tau^{+}$ 的关系反推。
- 速度场
$$
\begin{cases}
    u_x =& U_0 \tanh (\sigma\, (\frac{1}{4} - |\frac{y}{L} - \frac{1}{2}| )) \\
    u_y =& U_0 \delta \sin (2 \pi \, (\frac{x}{L} + \frac{1}{2}) )
\end{cases}
$$
- 压强场 $p$ ： 求解压力泊松方程。由于这个场是周期的，我就用伪谱法算。

代码我放在了[Github Gist](https://gist.github.com/Ives0721/742d6501df0767afc538fb3bcdff40bb)上，欢迎直接查看

[^Ginzburg2008]: Ginzburg, I., Verhaeghe, F., & d’Humières, D. (2008). Two-Relaxation-Time Lattice Boltzmann Scheme: About Parametrization, Velocity, Pressure and Mixed Boundary Conditions. Communications in Computational Physics, 3(2), 427–478.
[^Guo2002]: Guo, Z., Zheng, C., & Shi, B. (2002). Discrete lattice effects on the forcing term in the lattice Boltzmann method. Phys. Rev. E, 65, 046308.
[^Humières2009]: d’Humières, D., & Ginzburg, I. (2009). Viscosity independent numerical errors for Lattice Boltzmann models: From recurrence equations to “magic” collision numbers. Computers & Mathematics with Applications, 58(5), 823–840. https://doi.org/10.1016/j.camwa.2009.02.008
[^Postma2020]: Postma B, Silva G. Force methods for the two-relaxation-times lattice Boltzmann[J]. Physical Review E, 2020, 102(6): 063307.
[^Buick2000]: Buick, J. M., Greated, C. A. Gravity in a lattice Boltzmann model[J]. Physical Review E, 2000, 61(5): 5307--5320.
---
title: "双松弛时间（Two-relaxation-time）的格子 Boltzmann 方法"
tags: [格子Boltzmann方法]
categories:
- [格子Boltzmann方法]
---


# TRT-LBM

双松弛时间的格子 Boltzmann 方法是多松弛时间（Multi-Relaxation-Time）LBM 的一种特殊形式。它仅具有两个松弛时间，分为对称和非对称两部分，同时保持了与 LBGK 模型一致的简单算法和计算效率。

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

* $\mathbb{F}_i$ 为源项在 $i$ 方向上的函数。

## 源项及宏观量的计算

// https://www.sciencedirect.com/science/article/pii/S1290072923006129
// https://www.sciencedirect.com/science/article/pii/S0045793020301092

这里以二阶矩模型的源项为例：

$$
\mathrm{F}_{i} = w_i \cdot 
\left[\frac{\boldsymbol{c}_i - \boldsymbol{u}}{c_s^2} + \frac{(\boldsymbol{c}_i \cdot \boldsymbol{u}) \boldsymbol{c}_i}{c_s^4}\right] \cdot \boldsymbol{F}
$$

其中 $\boldsymbol{F}$ 为外力加速度。此时宏观密度 $\rho$ 、速度 $\boldsymbol{u}$ 和动量 $\boldsymbol{j}$ 的计算表示为：
$$
\rho = \sum\limits_{i=0}^{b-1} f_i ,\quad
\boldsymbol{j} = \rho \boldsymbol{u} = \frac{\boldsymbol{F}}{2} +  \sum\limits_{i=0}^{b-1} f_i \boldsymbol{c}_{i}
$$

令

$$
\begin{aligned}
F_{i}^{+} &= \frac{1}{2} (\mathrm{F}_i + \mathrm{F}_{-i}), \\
F_{i}^{-} &= \frac{1}{2} (\mathrm{F}_i - \mathrm{F}_{-i})
\end{aligned}
$$



这里以 Guo 等[^Guo2002]的源项为例：

$$
\mathrm{F}_{i} = \left(1 - \frac{1}{2 \tau^{+}} \right) w_i \cdot 
\left[\frac{\boldsymbol{c}_i - \boldsymbol{u}}{c_s^2} + \frac{(\boldsymbol{c}_i \cdot \boldsymbol{u}) \boldsymbol{c}_i}{c_s^4}\right] \cdot \boldsymbol{F}
$$

另一种做法则是参照Ginzburg等[^Ginzburg2008]的写法，这里不做展开。

<!--

另一种做法则是参照Ginzburg等[^Ginzburg2008]的写法，将源项 $\mathrm{F}_i(\boldsymbol{r},t)$ 表示为

$$
\mathrm{F}_i(\boldsymbol{r},t) = \mathrm{F}_i^{+}(\boldsymbol{r},t) + \mathrm{F}_i^{-}(\boldsymbol{r},t)
$$

其中 $\mathrm{F}_i^{+} = \frac{1}{2} (\mathrm{F}_i + \mathrm{F}_{-i})$， $\mathrm{F}_i^{-} = \frac{1}{2} (\mathrm{F}_i - \mathrm{F}_{-i})$。

下面记质量源 $M(\boldsymbol{r},t)$ 和体力 $\boldsymbol{F}(\boldsymbol{r},t)$ 。则宏观密度 $\rho$ 和动量 $\boldsymbol{j}$ 的计算表示为：

$$
\rho = \sum\limits_{i=0}^{b-1} f_i + \frac{M}{2} ,\quad
\boldsymbol{j} = \rho \boldsymbol{u} = \sum\limits_{i=0}^{b-1} f_i \boldsymbol{c}_{i} + \frac{\boldsymbol{F}}{2}
$$
-->


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
&\quad + \mathrm{F}_i \Delta t
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

[^Ginzburg2008]: Ginzburg, I., Verhaeghe, F., & d’Humières, D. (2008). Two-Relaxation-Time Lattice Boltzmann Scheme: About Parametrization, Velocity, Pressure and Mixed Boundary Conditions. Communications in Computational Physics, 3(2), 427–478.
[^Guo2002]: Guo, Z., Zheng, C., & Shi, B. (2002). Discrete lattice effects on the forcing term in the lattice Boltzmann method. Phys. Rev. E, 65, 046308.
[^Humières2009]: d’Humières, D., & Ginzburg, I. (2009). Viscosity independent numerical errors for Lattice Boltzmann models: From recurrence equations to “magic” collision numbers. Computers & Mathematics with Applications, 58(5), 823–840. https://doi.org/10.1016/j.camwa.2009.02.008

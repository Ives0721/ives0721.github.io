---
title: "LBM笔记：Exact Difference Method 构建的源项"
date: 2024-09-07 21:38:05
tags: [格子Boltzmann方法]
categories:
- [格子Boltzmann方法, 源项]
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


> 本文同步发表于[知乎](https://zhuanlan.zhihu.com/p/719081120)和[CSDN](https://blog.csdn.net/weixin_43890806/article/details/142057107)

我在阅读其他文献时，看到有文献提到 Kupershtokh 的一篇关于 LBM 外力项的文章[^Kupershtokh2004]，他提出的格式在其他文章中又被称为 *Exact difference method（EDM）*。

本想搜索这篇文章的PDF进行阅读，但最后只找到[一个网站](https://www.elibrary.ru/item.asp?id=28981868)介绍了部分信息。这份粗略笔记的内容，都是基于其他文献[^Kupershtokh2009][^Li2016][^Li_PRE_2016]中对该方法的介绍而整理的。



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
- $\tau = \tau_f / \delta_t$ 为无量纲松弛时间；
- $\bold{x}$ 为网格点位置，$\rho, \bold{u^*}$ 分别为该点在 $t$ 时刻的密度和速度（与宏观速度 $\bold{u}$ 有区别）；
- $f_\alpha, f_\alpha^{(eq)}$ 分别为 $\bold{e}_\alpha$ 方向的**分布函数**和**平衡态**；
$$
\begin{aligned}
f_\alpha^{(eq)} &= w_\alpha \rho \cdot \left[ 1 + \frac{\bold{e}_\alpha \cdot \bold{u}}{c_s^2} + \frac{(\bold{e}_\alpha \cdot \bold{u})^2}{2 c_s^4} - \frac{|\bold{u}|^2}{2 c_s^2} \right] \\
&= w_\alpha \rho \cdot \left[ 1 + \frac{\bold{e}_\alpha \cdot \bold{u}}{c_s^2} + \frac{1}{2 c_s^4}(\bold{u}\bold{u}):(\bold{e}_\alpha \bold{e}_\alpha - c_s^2 [\bold{I}]) \right] \\
\end{aligned}
\tag{4-1}
$$
这里 $w_\alpha$ 为 $\bold{e}_\alpha$ 的权重系数， $[\bold{I}]$ 为单位矩阵；
两个同样为m*n大小的矩阵 $[\bold{A}]$ 和 $[\bold{B}]$ ，则 $\displaystyle [\bold{A}] : [\bold{B}] = \sum^{m}_{i=1}\sum^{n}_{j=1} A_{ij} B_{ij}$ ；
- $F_\alpha$ 为外力项，其表达式基于外力 $\bold{F} = \rho \dfrac{\mathrm{d} \bold{u^*}}{\mathrm{d} t} = \rho \bold{a}$ 进行构建。

$$
F_\alpha = f_\alpha^{(eq)} (\bold{u^*}(t+\delta_t)) - f_\alpha^{(eq)} (\bold{u^*}(t))
$$

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
\mathbf{Error}_{\mathrm{EDM}} = -\delta_t^2 \nabla\cdot\left( \frac{\bold{FF}}{4 \rho} \right)
\tag{7}
$$


## EDM 源项的展开式

为便于分析，这里采用如下书写形式的平衡态：
$$
f_\alpha^{(eq)} = w_\alpha \rho \cdot \left[ 1 + \frac{\bold{e}_\alpha \cdot \bold{u}}{c_s^2} + \frac{1}{2 c_s^4}(\bold{u}\bold{u}):(\bold{e}_\alpha \bold{e}_\alpha - c_s^2 [\bold{I}]) \right]
$$

将其带入 EDM 源项，记 $\Delta\bold{u} = \frac{\bold{F}}{\rho} \delta_t$，整理得：

$$
\begin{aligned}
F_{\alpha} &= f_\alpha^{(eq)} (\bold{u^*}(t+\delta_t)) - f_\alpha^{(eq)} (\bold{u^*}(t)) \\&= 
\rho w_\alpha \left\{ \frac{\bold{e}_\alpha \cdot \Delta\bold{u}}{c_s^2} + 
\newline \frac{1}{2 c_s^4} \left[ (\bold{u}^* + \Delta\bold{u})(\bold{u}^* + \Delta\bold{u}) \right.\right. \\ &\quad 
\left. - \bold{u}^* \bold{u}^* \right] : \left.(\bold{e}_\alpha \bold{e}_\alpha - c_s^2 [\bold{I}]) \right\}
\end{aligned}
$$
这里涉及矩阵 $(\bold{u}^* + \Delta\bold{u})(\bold{u}^* + \Delta\bold{u})$ 和 $\bold{u}^* \bold{u}^*$ 的相减，由于：
$$
(\bold{u}^* + \Delta\bold{u})(\bold{u}^* + \Delta\bold{u}) = (u^*_{i} + \Delta u_i)(u^*_{j} + \Delta u_j), \quad
\bold{u}^* \bold{u}^* = u^*_{i} u^*_{j}
$$

其中 $u_i^*$ 和 $\Delta u_i = F_i \delta_{t} / \rho$ 分别为 $\bold{u}^*$ 和 $\Delta\bold{u}$ 在 i 轴上的分量，所以：
$$
\begin{aligned}
    (\bold{u}^* + \Delta\bold{u})(\bold{u}^* + \Delta\bold{u}) - \bold{u}^* \bold{u}^{*} &=
    (u^*_{i} + \Delta u_i)(u^*_{j} + \Delta u_j) - u^*_{i} u^*_{j} \\ &=  
    \frac{\delta_t}{\rho} (u_i^* F_j + u_j^* F_i) + \frac{\delta_t^2}{\rho^2} F_i F_j \\ &=
    \frac{\delta_t}{\rho} (u_i^* + \frac{F_i \delta_t}{2 \rho}) F_j + \frac{\delta_t}{\rho} (u_j^* + \frac{F_j \delta_t}{2 \rho}) F_i 
\end{aligned}
$$

这里令 $\bold{v}_{\rm{EDM}} = \bold{u}^{*} + \frac{\bold{F}}{2 \rho}  \delta_{t}$ ,则 EDM 源项的展开式化简为：

$$
F_\alpha = \delta_t w_\alpha \left[ \frac{\bold{e}_\alpha \cdot \bold{F}}{c_s^2} +
\frac{1}{2 c_s^4} \left( \bold{v}_{\rm{EDM}} \bold{F}  + \bold{F} \bold{v}_{\rm{EDM}} \right) : (\bold{e}_\alpha \bold{e}_\alpha - c_s^2 [\bold{I}]) \right]
\tag{8}
$$
上式即为 EDM 源项的展开结果。

<!-- 
## 源项误差的 Chapman-Enskog 分析

这里的推导主要参考 Li 等[^Li_PRE_2016]的文章，并补充了一些自己的笔记。

### 源项的各阶速度矩

根据 LBM 格子张量对称性，有前 5 阶表达式：
$$
\sum_{\alpha} w_{\alpha} = 1 ,\quad
\sum_{\alpha} w_{\alpha} e_{\alpha i} = 0 ,\\
\sum_{\alpha} w_{\alpha} e_{\alpha i} e_{\alpha j} = c_s^2 \delta_{ij} ,\quad
\sum_{\alpha} w_{\alpha} e_{\alpha i} e_{\alpha j} e_{\alpha k} = 0 \\
\sum_{\alpha} w_{\alpha} e_{\alpha i} e_{\alpha j} e_{\alpha k} e_{\alpha l} = c_s^4 (\delta_{ij} \delta_{kl} + \delta_{ik} \delta_{jl} + \delta_{il} \delta_{kj}) \\
\sum_{\alpha} w_{\alpha} e_{\alpha i} e_{\alpha j} e_{\alpha k} e_{\alpha l} e_{\alpha m} = 0
$$

其中 $e_{\alpha i}$ 表示 **D 维向量** $\bold{e}_{\alpha}$ 的第 $i$ 个分量， $w_\alpha$ 为 $\bold{e}_{\alpha}$ 的权重系数。 $\delta_{ij}$ 为克罗内克符号。

所以，对于 EDM 格式的源项 $F_{\alpha,\mathrm{EDM}}$ ，其各阶速度矩为：
$$
\sum_{\alpha} F_{\alpha,\mathrm{EDM}} = 0 ,\quad
\sum_{\alpha} \bold{e}_\alpha F_{\alpha,\mathrm{EDM}} = \bold{F} ,\quad\\
\sum_{\alpha} \bold{e}_\alpha \bold{e}_\alpha F_{\alpha,\mathrm{EDM}} = \delta_t \left( (\bold{u}^* \bold{F} + \bold{F}\bold{u}^*) + \delta_t \frac{\bold{FF}}{\rho} \right)
\tag{9}
$$

技术细节可见附录（A）。

### LBM 的 Chapman-Enskog 展开

令展开参数 $\epsilon$ 是与克努森数同阶的小量，则将分布函数展开为

$$
f_{\alpha} = f_{\alpha}^{(0)} + \epsilon f_{\alpha}^{(1)} + \epsilon^2 f_{\alpha}^{(2)} + ...,
$$

源项展开为 $F_{\alpha} = \epsilon F_{1 \alpha}$，时间和空间偏导展开为

$$
\frac{\partial}{\partial t} = \epsilon \frac{\partial}{\partial t_1} + \epsilon^2 \frac{\partial}{\partial t_2}
,\quad
\frac{\partial}{\partial \bold{x}} = \epsilon \frac{\partial}{\partial \bold{x}_1}.
$$

将式(4)中的 $f_\alpha (\bold{x} + \bold{e}_\alpha \delta_t, t+\delta_t)$ 在 $f_\alpha (\bold{x}, t)$ 处进行 Taylor 展开（至2阶项），并把上述展开式代入，可得 $\epsilon$ 各阶项的方程为：

$$
\begin{aligned}
O(\epsilon^0):\quad & 
    f_{\alpha}^{(0)} = f_{\alpha}^{(eq)} ,\\
O(\epsilon^1):\quad &
    \mathrm{D}_{1 \alpha} f_{\alpha}^{(0)} = -\frac{1}{\tau \delta_t} f_{\alpha}^{(1)} + \frac{F_{1 \alpha}}{\delta_t} ,\\
O(\epsilon^2):\quad &
    \frac{\partial f_{\alpha}^{(0)}}{\partial t_2} + \mathrm{D}_{1 \alpha} f_{\alpha}^{(1)} + \frac{\delta_t}{2} \mathrm{D}_{1 \alpha}^2 f_{\alpha}^{(0)} = -\frac{1}{\tau \delta_t} f_{\alpha}^{(2)}
    .\\
\end{aligned}
$$

其中 $\displaystyle \mathrm{D}_{1 \alpha} = \frac{\partial}{\partial t_1} + \bold{e}_{\alpha} \cdot \frac{\partial}{\partial \bold{x}_1}$。

# 附录
## （A）源项各阶速度矩的推导

这里仅以 $\displaystyle \sum_{\alpha} \bold{e}_\alpha \bold{e}_\alpha F_{\alpha,\mathrm{EDM}}$ 的计算举例来简单说明。

$$
\begin{aligned}
    \sum_{\alpha} \bold{e}_{\alpha} \bold{e}_{\alpha} F_{\alpha,\mathrm{EDM}}
    &= \sum_{\alpha} \bold{e}_\alpha \bold{e}_\alpha \cdot \delta_t w_\alpha \left[ \frac{\bold{e}_\alpha \cdot \bold{F}}{c_s^2} + \right.\\
    &\quad\left.
    \frac{1}{2 c_s^4} \left( \bold{v}_{\rm{EDM}} \bold{F}  + \bold{F} \bold{v}_{\rm{EDM}} \right) : (\bold{e}_\alpha \bold{e}_\alpha - c_s^2 [\bold{I}]) \right] \\
    &= \frac{\delta_t}{c_s^2} \sum_{\alpha} \left[ \bold{e}_\alpha \bold{e}_\alpha (w_\alpha \bold{e}_\alpha \cdot \bold{F}) \right] + \\
    &\quad \frac{\delta_t}{2 c_s^4} \sum_{\alpha} w_{\alpha} \bold{e}_\alpha \bold{e}_\alpha \cdot [\left( \bold{v}_{\rm{EDM}} \bold{F}  + \bold{F} \bold{v}_{\rm{EDM}} \right) : (\bold{e}_\alpha \bold{e}_\alpha - c_s^2 [\bold{I}])]
\end{aligned}
$$

一方面，对于 $\displaystyle \sum_{\alpha} \left[ \bold{e}_\alpha \bold{e}_\alpha (w_\alpha \bold{e}_\alpha \cdot \bold{F}) \right]$ 这个矩阵，它的第 $i$ 行第 $j$ 列写作
$$
\sum_{k} F_{k} \sum_{\alpha}(w_{\alpha} e_{\alpha i} e_{\alpha j} e_{\alpha k}) = \sum_{k} F_{k} \cdot 0 = 0
$$

另一方面，已知 $\bold{v}_{\rm{EDM}} = \bold{u}^{*} + \frac{\bold{F}}{2 \rho}  \delta_{t}$，由于 $\left( \bold{v}_{\rm{EDM}} \bold{F}  + \bold{F} \bold{v}_{\rm{EDM}} \right)$ 为矩阵。所以，下文为方便起见，我们记 

$$[\bold{\mathbb{B}}] = \bold{v}_{\rm{EDM}} \bold{F}  + \bold{F} \bold{v}_{\rm{EDM}} = \bold{u}^{*} \bold{F} + \bold{F} \bold{u}^{*} + \frac{\delta_t}{\rho}\bold{F}\bold{F}$$

则  $[\bold{\mathbb{B}}] : (\bold{e}_\alpha \bold{e}_\alpha - c_s^2 [\bold{I}])$ 是一个标量。
对于矩阵 $\displaystyle [\bold{\mathbb{A}}] = \sum_{\alpha} w_{\alpha} \bold{e}_\alpha \bold{e}_\alpha \cdot \{[\bold{\mathbb{B}}] : (\bold{e}_\alpha \bold{e}_\alpha - c_s^2 [\bold{I}])\}$ ，它的第 $i$ 行第 $j$ 列 $\mathbb{A}_{ij}$ 写作

$$
\begin{aligned}
\mathbb{A}_{ij} &= \sum_{\alpha} w_{\alpha} e_{\alpha i} e_{\alpha j} \left[ \sum_{k} \sum_{l} \mathbb{B}_{kl}  (e_{\alpha k} e_{\alpha l} - c_s^2 \delta_{kl}) \right]
\\&= \sum_{k} \sum_{l} \mathbb{B}_{kl} \sum_{\alpha} \left\{ w_{\alpha} e_{\alpha i} e_{\alpha j} (e_{\alpha k} e_{\alpha l} - c_s^2 \delta_{kl}) \right\}
\\&= \sum_{k} \sum_{l} \mathbb{B}_{kl} \cdot c_s^4 (\delta_{il} \delta_{jk} + \delta_{ik} \delta_{jl})
\\&= (\mathbb{B}_{ij} +\mathbb{B}_{ji}) c_s^4
\end{aligned}
$$

因为 $[\mathbb{B}]$ 是对称矩阵，所以 $\mathbb{A}_{ij} =2 c_s^4 \mathbb{B}_{ij}$。也就有：
$$
\begin{aligned}
\sum_{\alpha} \bold{e}_\alpha \bold{e}_\alpha F_{\alpha,\mathrm{EDM}} &= 
\vec{\bold{0}} + \frac{\delta_t}{2 c_s^4} \sum_{\alpha} w_{\alpha} \bold{e}_\alpha \bold{e}_\alpha \cdot [\left( \bold{v}_{\rm{EDM}} \bold{F}  + \bold{F} \bold{v}_{\rm{EDM}} \right) : (\bold{e}_\alpha \bold{e}_\alpha - c_s^2 [\bold{I}])] 
\\&=  \delta_t \left( \bold{v}_{\rm{EDM}} \bold{F}  + \bold{F} \bold{v}_{\rm{EDM}} \right)
\\&= \delta_t \left( \bold{u}^{*} \bold{F} + \bold{F} \bold{u}^{*} + \frac{\delta_t}{\rho}\bold{F}\bold{F} \right)
\end{aligned}
$$

综上， $\displaystyle \sum_{\alpha} \bold{e}_\alpha \bold{e}_\alpha F_{\alpha,\mathrm{EDM}} = \delta_t \left( \bold{u}^* \bold{F} + \bold{F}\bold{u}^* + \delta_t \frac{\bold{FF}}{\rho} \right)$ 成立。其他低阶矩均采用类似方法进行计算。当然，这也同样可以拓展至其他源项的速度矩计算。
-->

[^Kupershtokh2004]: Kupershtokh, A. L. [**New method of incorporating a body force term into the lattice Boltzmann equation**](https://www.elibrary.ru/item.asp?id=28981868). in *Proceeding of the 5th international EHD workshop* 241–246 (Poitiers, France, 2004).
[^Kupershtokh2009]: A.L. Kupershtokh, D.A. Medvedev, & D.I. Karpov (2009). **On equations of state in a lattice Boltzmann method**. Computers & Mathematics with Applications, 58(5), 965-974. DOI:10.1016/j.camwa.2009.02.024.
[^Li2016]: Q. Li, K.H. Luo, Q.J. Kang, Y.L. He, Q. Chen, & Q. Liu (2016). **Lattice Boltzmann methods for multiphase flow and phase-change heat transfer**. Progress in Energy and Combustion Science, 52, 62-105. DOI:10.1016/j.pecs.2015.10.001.
[^Li_PRE_2016]: Li, Q., Zhou, P., & Yan, H. (2016). **Revised Chapman-Enskog analysis for a class of forcing schemes in the lattice Boltzmann method**. Physical Review E, 94, 043313. DOI:10.1103/PhysRevE.94.043313.
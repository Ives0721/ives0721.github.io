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

## 格子Boltzmann方法

广义上说, LBM是离散Boltzmann方法(Discrete Boltzmann Method, DBM)的一种, 是式(1.1)的一种数值格式.本节介绍从式(1.1)到格子Boltzmann方程的离散思路, 推导中仅考虑外力加速度 $\mathbf{a} =\mathbf{F}/ m \equiv 0$ 的情况.

对式(1.1)沿着特征线 $\mathbf{\xi} = \mathbf{x}/t$ 在 $\Delta t$ 时间段内作积分, 略去 $O({\Delta t}^{2})$ , 得

$$
f( \mathbf{x} + \mathbf{\xi}\Delta t,\mathbf{\xi},t + \Delta t ) - f( \mathbf{x},\mathbf{\xi,}t) = -\frac{1}{\tau} \left\lbrack f(\mathbf{x},\mathbf{\xi},t) - f^{(eq)}(\mathbf{x},\mathbf{\xi},t) \right\rbrack
\tag{1.4}
$$

其中 $\tau = \tau_{0}/\Delta t$ 为无量纲松弛时间.为了数值地计算这一连续方程, 需要将其在 $\mathbf{\xi}$ 的空间中进行离散, 并使其平衡态满足速度矩的约束, 即：

$$
\int {f^{(eq)}\left( \mathbf{x},\mathbf{\xi,}t \right) \cdot \Psi\left( \mathbf{\xi} \right)}d\mathbf{\xi} = 
\sum_{\mathbf{i}} {w_{i}\Psi\left( \mathbf{\xi}_{i} \right)f^{(eq)}\left( \mathbf{x},\mathbf{\xi}_{i},t \right)} 
\tag{1.5}
$$

这里 $\Psi\left( \mathbf{\xi} \right)$ 为 $\mathbf{\xi}$ 的多项式, $\mathbf{\xi}_i$ 为离散速度, $w_{i}$ 为积分的权系数.具体到Boltzmann方程中, 则至少应满足质量和动量守恒：

$$
\rho = \sum_{i}^{}f_{i} = \sum_{i} f_{i}^{(eq)},\, \mathbf{j} = \rho\mathbf{u} = \sum_{i} {\mathbf{\xi}_{i}f_{i}} = \sum_{i} {\mathbf{\xi}_{i}f_{i}^{(eq)}}
$$

以及内能 $e$ 的方程：

$$
\rho e = \frac{1}{2}\sum_{i}^{}{\left( \mathbf{\xi}_{i}\mathbf{- u} \right)^{2}f_{i}} = \frac{1}{2}\sum_{i}^{}{\left( \mathbf{\xi}_{i}\mathbf{- u} \right)^{2}f_{i}^{(eq)}}
$$

其中 $f_{i} = W_{i}f\left( \mathbf{x},\mathbf{\xi}_{i}\mathbf{,}t \right)$ 和 $f_{i}^{(eq)} = W_{i}f^{(eq)}\left( \mathbf{x},\mathbf{\xi}_{i}\mathbf{,}t \right)$ .选定合适的离散速度集 $\{\mathbf{c}_{i}\}(i = 0,1,\ldots,b - 1)$ 及其对应的权重后, Boltzmann-BGK方程即可被数值计算.一般地, *n* 维空间下 b 速度的离散速度模型简称 D*n*Qb 模型.

在低马赫数下, 平衡态可对 $\mathbf{u}$ 进行Taylor展开, 并截断至二阶：

$$
\begin{aligned}
f^{(eq)}\left( \mathbf{x},\mathbf{\xi,}t \right) &= \frac{\rho}{\left( 2\pi R_{g}T \right)^{\frac{D}{2}}} \cdot \exp\left( - \frac{\mathbf{\xi}^{2}}{2R_{g}T} \right) \cdot \exp\left( \frac{2\left( \mathbf{\xi \cdot u} \right) - \mathbf{u}^{2}}{2R_{g}T} \right) \\
&= \frac{\rho\exp\left( - \frac{\mathbf{\xi}^{2}}{2R_{g}T} \right)}{\left( 2\pi R_{g}T \right)^{\frac{D}{2}}} \cdot \left\lbrack 1 + \frac{\mathbf{\xi \cdot u}}{R_{g}T} + \frac{\left( \mathbf{\xi \cdot u} \right)^{2}}{{2\left( R_{g}T \right)}^{2}} - \frac{\mathbf{u}^{2}}{2R_{g}T} \right\rbrack + O\left( \mathbf{u}^{3} \right)
\end{aligned}\tag{1.6}
$$

因此, 当使用 **截断的平衡态** 计算 $\int {f^{(eq)} \cdot \Psi\left( \mathbf{\xi} \right)}d\mathbf{\xi}$ 时, 会导出形如 $\int {\exp\left( - x^{2} \right) \cdot \Psi(x)}dx$ 的积分, 它可以用高斯积分进行数值计算.根据上式的具体情况, 记

$$
I_{m} = \int_{- \infty}^{+ \infty}{\exp\left( - \zeta^{2} \right) \cdot \zeta^{m}} d\zeta
$$

其中 $\zeta = \xi/\sqrt{2R_{g}T}$ .$I_{m}$ 是权函数 $\exp\left(-\zeta^{2}\right)$ 在某一实轴上的 $m$ 阶矩.

这里以D2Q9离散速度集为例介绍权重的计算, 其离散速度集为

$$
\mathbf{c}_{\alpha} = \begin{cases}
  (0,0) & ,\alpha = 0 \\
  c\left( \cos\frac{\pi(\alpha - 1)}{2},\sin\frac{\pi(\alpha - 1)}{2} \right) & ,\alpha = 1\ldots 4 \\
  \sqrt{2}c\left( \cos\frac{\pi(2\alpha - 1)}{4},\sin\frac{\pi(2\alpha - 1)}{4} \right) & ,\alpha = 5\ldots 8 \\
\end{cases}
\tag{1.7}
$$

其中 $c = \Delta x/\Delta t$ .二维情况下的 $\Psi$ 可定义为 $\Psi_{m,n}\left( \mathbf{\xi} \right) = \xi_{x}^{m}\xi_{y}^{n}$, 则式(1.6)被展开为 [^12] ：

$$
\begin{aligned}
I & = \int_{}^{}{f^{(eq)} \cdot \Psi_{m,n}\left( \mathbf{\xi} \right)}d\mathbf{\xi} \\
 & = \frac{\rho}{\pi}\left( 2R_{g}T \right)^{\frac{m + n}{2}}\left\{ \left( 1 - \frac{u^{2}}{2R_{g}T} \right)I_{m}I_{n} \right.\  \\
 & + \frac{2}{\sqrt{2R_{g}T}}\left( u_{x}I_{m + 1}I_{n} + u_{y}I_{m}I_{n + 1} \right) \\
 & \left. \  + \frac{1}{R_{g}T}\left( u_{x}^{2}I_{m + 2}I_{n} + 2u_{x}u_{y}I_{m + 1}I_{n + 1} + u_{y}^{2}I_{m}I_{n + 2} \right) \right\} \\
\end{aligned}\tag{1.8}
$$

对D2Q9模型,  $I_{m}$ 可被替换为三点Gauss-Hermite数值积分：

$$
I_{m} = \int_{- \infty}^{+ \infty}{\exp\left( - \zeta^{2} \right) \cdot \zeta^{m}}d\zeta = \sum_{\alpha = 1}^{3}{\omega_{\alpha}\zeta_{\alpha}^{m}}
$$

> **[NOTE]**  
> 积分点 $\zeta_{i}$ 是Hermite多项式 $H_{3}(x) = 8x^{3} - 12x$ 的零点, 对应权重 $\omega_{i} = \frac{2^{3 - 1} \cdot 3!}{H_{2}\left( \zeta_{i} \right)^{2} \cdot 3^{2}} \cdot \sqrt{\pi}$.  
> 所以积分点和权重为：$\zeta_{1} = - \sqrt{3/2},w_{1} = \sqrt{\pi}/6$；$\zeta_{2} = 0,w_{2} = 2\sqrt{\pi}/3$；$\zeta_{3} = \sqrt{3/2},\omega_{3} = \sqrt{\pi}/6$.

则式(1.8)可被数值离散为 [^12] ：

$$
I = \rho\sum_{i,j = 1}^{3}{\Psi\left( \mathbf{\xi}_{i,j} \right) \cdot \frac{\omega_{i}\omega_{j}}{\pi} \cdot \left\lbrack 1 + \frac{\mathbf{\xi}_{i,j}\mathbf{\cdot u}}{R_{g}T} + \frac{\left( \mathbf{\xi}_{i,j}\mathbf{\cdot u} \right)^{2}}{{2\left( R_{g}T \right)}^{2}} - \frac{\mathbf{u}^{2}}{2R_{g}T} \right\rbrack}
\tag{1.9}
$$

此处 $\mathbf{\xi}_{i,j}\mathbf{=}\left( \xi_{i}\mathbf{,}\xi_{j} \right)\mathbf{=}\sqrt{2R_{g}T}\left( \zeta_{i}\mathbf{,}\zeta_{j} \right)$ .也就是说, 若令 $R_{g}T = c_{s}^{2} = c^{2}/3$ , 则 $\mathbf{\xi}_{i,j}$ 便对应D2Q9模型的 $\mathbf{c}_{\alpha}$ .根据式(1.9),  $\mathbf{c}_{\alpha}$ 方向的平衡态分布 $f_{\alpha}^{(eq)}$ 定义为

$$
f_{\alpha}^{(eq)} = \rho w_{\alpha} \cdot \left\lbrack 1 + \frac{\mathbf{c}_{\alpha}\mathbf{\cdot u}}{c_{s}^{2}} + \frac{\left( \mathbf{c}_{\alpha}\mathbf{\cdot u} \right)^{2}}{2c_{s}^{4}} - \frac{\mathbf{u}^{2}}{2c_{s}^{2}} \right\rbrack
\tag{1.10}
$$

其中

$$
w_{\alpha} = \begin{cases}
    4/9  ,& i = j = 2 ,& \alpha = 0 \\
    1/9  ,& (i = 2)\bigoplus(j = 2) ,& \alpha = 1, \ldots, 4 \\
    1/36 ,& (i \neq 2)\bigwedge(j \neq 2) ,& \alpha = 5, \ldots, 8 \\
\end{cases}
$$

Qian等[^15]已给出常用的DnQb模型参数, 大量学者也给出了其他离散速度模型的推导方式 [^12][^13][^14].

在前文中, $R_{g} T$ 为常数一直是LBGK方程的推导前提, 这使得LBGK模型的最佳适用范围被限制在等温流动或无温流动.LBGK在模拟不可压缩流体时, 格子系统的马赫数应小于 0.15 [^16].而LBGK框架虽然可通过部分手段实现亚音速[^18]和高雷诺数[^19]的模拟, 但仍旧无法直接解决热流问题.

就目前来说, 热流体的LBM模拟可以分为多速LBE模型 [^1] [^2] [^3] [^6] [^7] [^20] [^21] 、双分布函数LBE模型 [^8] [^9] [^10] [^11] 、基于简化LBM的实现[^22]、热流格子Boltzmann通量求解器[^23]、多松弛模型[^29] [^30] [^31] 以及与其他数值方法结合的混合模型等多种形式 [^24] .

# 多速LBE模型

多速LBE模型(即Multi-Speed LBE)是一种完全基于LBM经典算法的改进.这类方法演化方程与常规LBM一致, 在保持算法不变的基础上通过使用比常规LBM更为复杂的离散速度模型还原高阶速度矩.它因只使用单一分布函数,物理描述上更接近事实而受到了广泛关注.并且使用更多的离散速度达到高阶精度的思路亦已经被广泛拓展到各种各样的LBM模拟中. 虽然Qian[^6]和Alexander等[^7]才是早期的提出者, 但本章首先以Chen等 [^1] 的高阶MS-LBGK模型的二维形式为例引入多速LBM的概念, 随后介绍几种不同思路下的MS-LBM改进形式.

## 高阶MS-LBGK模型

Chen等 [^1] 提出的多速LBM的平衡态分布函数定义为

$$
\begin{aligned}
f_{\sigma i}^{(eq)} = & A_{\sigma} + B_{\sigma} \cdot \left( \mathbf{c}_{\sigma i} \cdot \mathbf{u} \right) + C_{\sigma}\left\| \mathbf{u} \right\|^{2} + D_{\sigma} \cdot \left( \mathbf{c}_{\sigma i} \cdot \mathbf{u} \right)^{2} \\
 & + E_{\sigma}\left\| \mathbf{u} \right\|^{2} \cdot \left( \mathbf{c}_{\sigma i} \cdot \mathbf{u} \right) + F_{\sigma} \cdot \left( \mathbf{c}_{\sigma i} \cdot \mathbf{u} \right)^{3} \\
 & + G_{\sigma}\left\| \mathbf{u} \right\|^{2} \cdot \left( \mathbf{c}_{\sigma i} \cdot \mathbf{u} \right)^{2} + H_{\sigma}\left\| \mathbf{u} \right\|^{4} \\
\end{aligned}
\tag{2.1}
$$

其中 $A_{\sigma}$ 到 $H_{\sigma}$ 均为待定参数, 它们均可被表示为内能 $e$ 和密度 $\rho$ 的复合函数, 即 $X_{\sigma} = \left( x_{\sigma 0} + x_{\sigma 1}e + x_{\sigma 2}e^{2} \right)\rho$ .并且 $f_{\sigma i}^{(eq)}$ 的各阶速度矩需要满足下列约束条件：

$$
\sum_{\sigma i} f_{\sigma i}^{(eq)} = \rho,\, 
\sum_{\sigma i} {\mathbf{c}_{\sigma i}f_{\sigma i}^{(eq)}} = \rho\mathbf{u},\, 
\sum_{\sigma i} {\mathbf{c}_{\sigma i}\mathbf{c}_{\sigma i}f_{\sigma i}^{(eq)}} = \rho\mathbf{uu} + p\mathbf{I}
$$

$$
\sum_{\sigma i}^{}{\mathbf{c}_{\sigma i}\mathbf{c}_{\sigma i}\mathbf{c}_{\sigma i}f_{\sigma i}^{(eq)}} = \rho\mathbf{uuu} + p\ \lbrack\mathbf{u\delta}\rbrack
$$

$$
\sum_{\sigma i}^{}{\mathbf{c}_{\sigma i}\mathbf{c}_{\sigma i} \vert{c_{\sigma i}}\vert ^{2}f_{\sigma i}^{(eq)}} = \rho|u|^{2}\mathbf{uu}\mathbf{+}p(D + 4)\mathbf{uu}\mathbf{+}p|u|^{2}\mathbf{I}\mathbf{+}\frac{2(D + 2)pe}{D}
$$

其中压强 $p = \frac{2\rho e}{D}$ , 三阶张量 $\left\lbrack \mathbf{u\delta} \right\rbrack = \ u_{\alpha}\delta_{\beta\gamma} + u_{\beta}\delta_{\alpha\gamma} + u_{\gamma}\delta_{\alpha\beta}$ .

为满足这五个守恒方程的约束, 张量 $\displaystyle \sum_{i} \overset{n}{\overbrace{\mathbf{c}_{\sigma i}\mathbf{c}_{\sigma i}...\mathbf{c}_{\sigma i}}}\,$ 需要在 0~6 阶均为各向同性.因此Chen等 [^1] 使用的离散速度集表示为

$$
\mathbf{c}_{\sigma i} = \mathbf{c}_{pki}^{'} = \mathbf{Perm}\left\{ k\left( \underbrace{\pm 1,\ldots, \pm 1}_{p} \mathbf{,} \overset{D - p}{\overbrace{0,...,0}} \right) \right\}
$$

这里 $\mathbf{Perm}$ 是由 $(*)$ 的所有排列构成的集合, $k$ 为格子的缩放系数, 且 $\sigma = c_{\sigma i}^{2} = k^{2}p$.

因此, 高阶MS-LBGK模型的运动粘度为 $\mu = \frac{2\rho e}{D}\left( \tau - \frac{1}{2} \right)$ , 体积粘度 $\lambda = - \frac{2\mu}{D}$ , 热扩散系数 $\kappa = \frac{(D + 2)\rho e}{D}\left( \tau - \frac{1}{2} \right)$ .对于气体常数为 $R_{g}$ 的单原子分子, 定压和定容比热容分别为 $c_{p} = (D + 2)R_{g}/2$ 和 $c_{v} = D R_{g}/2$ .由于无量纲格子系统中 $R_{g} = 1$ , 该模型的Prandtl数为 $Pr = \mu c_{p}/\kappa = 1$ .这也是MS-LBGK模型这类建模方式的常见问题.

目前, 学界已有数类方式在多速LBE模型中实现对Prandtl数的自由调节, 如：双松弛时间方法[^2]；引入熵格式[^3] [^20] [^21] ；修改网格速度的定义 [^13] [^28] [^32] 等.由于已有文献指出若不局限于常规晶格模型则可获得更高精度的离散格式, 并且大多数方法都涉及对D*n*Qb模型的修改, 因此这里仅简要介绍前两种方式的实现思路.

## 对多速LBE模型的改进

### 基于双松弛时间的实现

基于双松弛时间的方法将碰撞分为由 $f_{i}$ 和 $f_{-i}$ 控制的两个部分, 并分别进行松弛.以Chen等[^2]提出的双松弛模型为例, 其演化方程为(不考虑外力项)：

$$
f_{i}\left( \mathbf{x} + \mathbf{c}_{i}\Delta t,t + \Delta t \right) - f_{i}\left( \mathbf{x},t \right) = \Omega_{i} + \Omega_{\mathbb{i}}
\tag{2.2}
$$

其中

$$
\Omega_{i} = - \frac{1}{\tau_{1}}\left( f_{i} - f_{i}^{(eq)} \right),\, \Omega_{\mathbb{i}} = - \frac{1}{\tau_{2}}\left( f_{- i} - f_{- i}^{(eq)} \right)
$$

下标 $-i$ 表示 $\mathbf{c}_{- i} = - \mathbf{c}_{i}$ 方向.$\tau_{v} = \tau_{1}\tau_{2}/(\tau_{1} + \tau_{2})$ ,  $\tau_{k} = \tau_{1}\tau_{2}/(\tau_{2} - \tau_{1})$ .Chen等[^2]指出其导出的宏观方程组为

$$\frac{D\rho}{Dt} = \rho\nabla \cdot \mathbf{u}$$

$$\frac{\partial\left( \rho\mathbf{u} \right)}{\partial t} + \nabla \cdot \left( \rho\mathbf{uu} \right) = - \nabla p + \nabla \cdot \mathbf{T}_{v}$$

$$\frac{\partial(\rho e)}{\partial t} + \nabla \cdot \left( \rho e\mathbf{u} \right) = p\nabla \cdot \mathbf{u} + \nabla \cdot (\kappa\nabla T) + \mathbf{T}_{k}:\nabla\mathbf{u}$$

其中

$$
\mathbf{T}_{v} = 2\mu_{v}\mathbf{S} + \lambda_{v}\left( \nabla \cdot \mathbf{u} \right)\mathbf{I},\mathbf{\, }\mathbf{T}_{k} = 2\mu_{k}\mathbf{S} + \lambda_{k}\left( \nabla \cdot \mathbf{u} \right)\mathbf{I}
$$

$\mathbf{S}$ 为应变率张量.各输运系数为

$$
\mu_{v} = \frac{2}{D}\rho e\left( \tau_{v} - \frac{1}{2} \right),\, \lambda_{v} = \frac{- 4}{D^{2}}\rho e\left( \tau_{v} - \frac{1}{2} \right)
$$

$$
\mu_{k} = \frac{2}{D}\rho e\left( \tau_{k} - \frac{1}{2} \right),\, \lambda_{k} = \frac{- 4}{D^{2}}\rho e\left( \tau_{k} - \frac{1}{2} \right)
$$

温度 $T = De/2$ , 热传导系数 $\kappa = \rho e(D + 2)\left( \tau_{k} - 1/2 \right)/D$.

综上所述, 该模型的Prandtl数为

$$
Pr = \frac{\mu_{v}c_{p}}{\kappa} = \frac{2\tau_{v} - 1}{2\tau_{k} - 1}
$$

这种方法通过使用两个松弛时间 $\tau_{1},\tau_{2}$ 分别对动量和能量方程进行松弛($\tau_{1},\tau_{2} > \frac{1}{2}$), 缓解LBGK单个松弛时间的参数依赖性.但由于能量方程的粘性应力项 $\mathbf{T}_{k} \neq \mathbf{T}_{v}$ , 因而其能量方程是与动量方程存在冲突的.

### 基于熵函数的实现

熵格子Boltzmann方程(Entropic Lattice Boltzmann Equation, ELBE)是另一种基于传统LBE的修改.这里以单松弛的ELBE为例, 它将LBE的演化方程修改为

$$
f_{i}\left( \mathbf{x} + \mathbf{c}_{i}\Delta t,t + \Delta t \right) - f_{i}\left( \mathbf{x},t \right) = - \alpha\beta\left( f_{i} - f_{i}^{(eq)} \right)
\tag{2.3}
$$

其中 $\beta = 1/(1 + 2\nu/c_{s}^{2})$ , $\nu$ 为运动粘度.$\alpha$ 为方程 $H\left( f + \alpha\left( f - f^{(eq)} \right) \right) = H(f)$ 的解, 熵函数 $H$ 定义为[^4]

$$
H = \sum_{i}^{}{f_{i}\ln\left( \frac{f_{i}}{w_{i}} \right)}
\tag{2.4}
$$

$\alpha$ 的解可通过迭代法或近似公式计算[^5].由于部分常规D*n*Qb模型在使用 $f^{(eq)}$ 的低阶展开(如式(1.10))时无法满足熵方程[^17], 因此该方法通常需要对D*n*Qb模型进行修改, 或搭配多速D*n*Qb模型一同使用.

Pransianakis等[^3] [^21] 在ELBE的总体框架上, 将碰撞项视为从 $f_{i}$ 到中间态 $f_{i}^{*}$ 再到平衡态 $f_{i}^{(eq)}$ 的两步松弛, 即：

$$- \frac{1}{\tau_{1}}\left( f_{i} - f_{i}^{*} \right) - \frac{1}{\tau_{2}}\left( f_{i}^{*} - f_{i}^{(eq)} \right)$$

中间态 $f_{i}^{*}$ 是平衡态的扰动, 即：$f_{i}^{*} = f_{i}^{(eq)} + \delta f_{i}^{*}$.校正量$\delta f_{i}^{*}$ 的计算可基于 "中间态的热通量不变" 假设构造 $\delta f_{i}^{*}$ 的方程组进行求解.松弛时间为 $\tau_{1} = \mu/\rho T_{0}$ 和 $\tau_{2} = 2\kappa/\rho T_{0}$ , 分别控制应力张量和热流量的松弛, $T_{0}$ 是格子系统中的参考温度.因此, Prandtl数为 $\Pr = 4\tau_{1}/\tau_{2}$ .

> **[NOTE]**  
> 需要强调的是, 在Pransianakis等[^3]在文中使用的D2Q9模型里, 虽然特征速度方向与式(1.7)一致, 但各个$\mathbf{c}_{i}$的权重 $W_{i}$ 是温度$T$的函数, 即： $\displaystyle W_{i} = \frac{\left( 1 - T^{2} \right)T}{2(1 - T)}c_{i}^{2}$  
> 
>而平衡态分布写作： $\displaystyle f_{i}^{(eq)} = \rho W_{i}\left\{ 1 + \frac{\mathbf{c}_{i} \cdot \mathbf{j}}{\rho T} + \frac{\mathbf{j}\mathbf{j}^{T}}{2(\rho T)^{2}}:\left\lbrack \mathbf{c}_{i}\mathbf{c}_{i}^{T} - \frac{4T^{2} + c_{i}^{2}(1 - 3T)}{2(1 - T)}\mathbf{I} \right\rbrack \right\}$

Frapolli等[^20]提出了将多速LBE模型与熵LBE框架进行结合的思路.以D2Q25-ZOT模型为例, 其演化方程为：

$$
f_{i}\left( \mathbf{x} + \mathbf{c}_{i}\Delta t,t + \Delta t \right) - f_{i}\left( \mathbf{x},t \right) = \omega_{1}\left( f_{i} - f_{i}^{*} \right) + \omega_{2}\left( f_{i}^{*} - f_{i}^{(eq)} \right)
\tag{2.5}
$$

在不考虑ELBE的校正时, 动力粘度 $\mu$ 和热扩散系数 $\kappa$ 的表达式为：

$$
\mu = \begin{cases}
    \left( \frac{1}{\omega_{1}} - \frac{1}{2} \right)\rho T,& \Pr \leq 1 \\
    \left( \frac{1}{\omega_{2}} - \frac{1}{2} \right)\rho T,& \Pr \geq 1
\end{cases}, \,
\kappa = \begin{cases}
    \left( \frac{1}{\omega_{2}} - \frac{1}{2} \right)\rho Tc_{p},& \Pr \leq 1 \\
    \left( \frac{1}{\omega_{1}} - \frac{1}{2} \right)\rho Tc_{p},& \Pr \geq 1
\end{cases}
$$

Frapolli等[^20]认为, 上述计算中仅需要对涉及 $\mu$ 的部分进行校正.当 $\Pr < 1$ 时, $\omega_{1} = \alpha\beta_{1},\omega_{2} = \beta_{2}$ ；当 $\Pr > 1$ 时, $\omega_{1} = \beta_{1},\omega_{2} = \alpha\beta_{2}$ . $\alpha$ 仍为方程 $H\left( f + \alpha\left( f - f^{(eq)} \right) \right) = H(f)$ 的解.因此, $\mu$ 和 $\kappa$ 的计算变为：

$$
\mu = \left\{ \begin{matrix}
\frac{1}{2}\left( \frac{1}{\beta_{1}} - 1 \right)\rho T,\, \Pr \leq 1 \\
\frac{1}{2}\left( \frac{1}{\beta_{2}} - 1 \right)\rho T,\, \Pr \geq 1 \\
\end{matrix} \right.\ ,\, 
\kappa = \left\{ \begin{matrix}
\left( \frac{1}{\beta_{2}} - \frac{1}{2} \right)\rho T,\, \Pr \leq 1 \\
\left( \frac{1}{\beta_{1}} - \frac{1}{2} \right)\rho T,\, \Pr \geq 1 \\
\end{matrix} \right.\ .
$$

> **[NOTE]**  
> 当 $\Pr = 1$ 时模型退化为LBGK方程, 这里不做讨论.

上述所介绍的基于ELBE的热流模拟能够成功还原Fourier-Navier-Stokes方程组, 实现低马赫数下各种热流的数值模拟, 但也没有实现热流和流场的耦合.此外, 由于 $H$ 函数的引入, 这类方法需要额外求解 $H\left( f + \alpha\left( f - f^{(eq)} \right) \right) = H(f)$ , 在数值计算上导致一定不便.

# 双分布函数LBE模型

双分布函数(double distribution function, DDF)的LBE模型的主要思路是将速度场和温度场分别用两套LBE进行计算.DDF-LBE下的热流分为两类：其一是将温度视作被动标量, 这要求温度变化对流动的影响较小；另一类考虑的温度对流场的影响(如可压缩热流).

## 不考虑黏性热耗散和压缩功的DDF-LBE

### 基本框架

Bartoloni等 [^8] 指出：在压力做功和黏性热耗散项均可被忽略的场景中, 可将流场和温度场分别采用两套分布函数体系进行对流扩散.一般来说, 无源项的DDF-LBE演化方程组为：

$$
f_{i}\left( \mathbf{x} + \mathbf{c}_{i}\Delta t,t + \Delta t \right) - f_{i}\left( \mathbf{x},t \right) = - \frac{1}{\tau_{1}}\left( f_{i} - f_{i}^{(eq)} \right)\tag{3.1}
$$

$$
g_{i}\left( \mathbf{x} + \mathbf{c}_{i}\Delta t,t + \Delta t \right) - g_{i}\left( \mathbf{x},t \right) = - \frac{1}{\tau_{2}}\left( g_{i} - g_{i}^{(eq)} \right)\tag{3.2}
$$

其中 $f_{i}$ 和 $g_{i}$ 分别用于计算速度和温度, 并通过如下速度矩还原宏观物理量

$$
\rho = \sum_{i}^{}f_{i},\,
\mathbf{j} = \rho\mathbf{u} = \sum_{i}^{}{\mathbf{c}_{i}f_{i}},\,
\rho T = \sum_{i}^{}g_{i}
\tag{3.3}
$$

需要强调的是, $f_{i}$ 和 $g_{i}$ 分别使用不同的DnQb模型, 因此它们的格子声速分别为 $c_{s1}$ 和 $c_{s2}$ .两个松弛时间( $\tau_{1},\tau_{2}$ )和运动粘度 $\nu$ 、热扩散系数 $\kappa$ 的关系为：

$$
\tau_{1} = \frac{\nu}{c_{s1}^{2}} + \frac{1}{2},\,
\tau_{2} = \frac{\kappa}{c_{s2}^{2}} + \frac{1}{2}
\tag{3.4}
$$

这种思路仅额外添加了一组分布函数及其平衡态的定义, 并且由于其易于实现, 目前是低速小Mach数场景下热流的一种主流LBM计算方法.

## 含黏性热耗散和压缩功的DDF-LBE

前面介绍的DDF-LBE都是在忽略黏性热耗散和压缩功的前提下, 单独另开一个分布函数 $g_{i}$ 来计算温度 $T$ 这一标量的变化.虽然基于温度的建模方式足够简单, 但其应用场景限制较大.为更全面地计算低马赫数热流中流场和温度的相互作用, 学界已经提出了一些方案.

### 基于能量的DDF-LBE

既然可以为温度$T$构建一套LBE, 那么为其他可以表示能量的宏观量(如内能 $\rho e$ [^9] 、总能 $\rho E = \rho e + \rho u^{2}/2$ [^10] 、总焓 $h = e + p/\rho$ [^11] )构建一套LBE也是可行的.通过模拟能量场的演化, 更准确地描述流场和温度场之间的作用.由于这套思路的建模方式总体上较为相似, 下面仅以基于总能$\rho E$的DDF-LBE为例对这种思路进行介绍.

单位体积内的总能可以由分布函数的速度矩表示为

$$
\rho E = \rho e + \frac{\rho u^{2}}{2} = \rho c_{v}T + \frac{\rho u^{2}}{2} = \int_{}^{}{\frac{1}{2}\xi^{2}f}d\mathbf{\xi}
\tag{3.5}
$$

因此可以将总能的分布函数定义为 $h = \frac{1}{2}\xi^{2}f$ , 则 $\rho E = \int h d\mathbf{\xi}$ .

根据 $f$ 的演化方程(即式(1.1)), 可得到 $h$ 的演化方程为

$$
\frac{\partial h}{\partial t} + \mathbf{\xi} \cdot \nabla f + a \cdot \nabla_{\mathbf{\xi}}h - f\mathbf{\xi \cdot a} = \Omega_{h}
\tag{3.6}
$$

其中 $\Omega_{h} = \frac{1}{2}\xi^{2}\Omega_{f}$ 为总能的碰撞算子.Guo等 [^10] 将 $\Omega_{h}$ 视作机械能和内能两部分作用的组合, 分别采用 $\tau_{f}$ 和 $\tau_{h}$ 控制, $\Omega_{h}$ 则表示为

$$
\Omega_{h} = - \frac{1}{\tau_{h}}\left( h - h^{(eq)} \right) + \frac{Z}{\tau_{hf}}\left( f - f^{(eq)} \right)
\tag{3.7}
$$

其中 $Z = \mathbf{\xi \cdot u} - \frac{u^{2}}{2}$, $\frac{1}{\tau_{hf}} = \frac{1}{\tau_{h}} - \frac{1}{\tau_{f}}$.

多原子分子气体的 $f$ 也受到分子转动和振动的影响, 因此其平衡态 $f^{(eq)}$ 为：

$$
{\widetilde{f}}^{(eq)}\left( \mathbf{x,\xi,}\mathbf{\eta},t \right) = \frac{\rho}{\left( 2\pi R_{g}T \right)^{\frac{D + K}{2}}} \cdot \exp\left\lbrack - \frac{\left( \mathbf{\xi} - \mathbf{u} \right)^{2} + \mathbf{\eta}^{2}}{2R_{g}T} \right\rbrack
\tag{3.8}
$$

其中 $\mathbf{\eta}$ 是一个含有K个元素的内部自由度矢量.引入 $\overline{f} = \int f d\mathbf{\eta}$, $\overline{h}\mathbf{=}\int \frac{\mathbf{\xi}^{2}\mathbf{+}\mathbf{\eta}^{2}}{2}f d\mathbf{\eta\ }$ , 则可得到一组新的输运方程：

$$
\frac{\partial\overline{f}}{\partial t} + \mathbf{\xi} \cdot \nabla\overline{f} + \mathbf{a} \cdot \nabla_{\mathbf{\xi}}\overline{f} = - \frac{1}{\tau_{f}}\left( \overline{f} - {\overline{f}}^{(eq)} \right),
\tag{3.9}
$$

$$
\frac{\partial\overline{h}}{\partial t} + \mathbf{\xi} \cdot \nabla\overline{h} + a \cdot \nabla_{\mathbf{\xi}}\overline{h} - \overline{f}\mathbf{\xi \cdot a} = - \frac{1}{\tau_{h}}\left( \overline{h} - {\overline{h}}^{(eq)} \right) + \frac{Z}{\tau_{hf}}\left( \overline{f} - {\overline{f}}^{(eq)} \right).
\tag{3.10}$$

它们的平衡态分别为

$$
{\overline{f}}^{(eq)} = \int_{}^{}{\widetilde{f}}^{(eq)}d\mathbf{\eta =}\frac{\rho}{\left( 2\pi R_{g}T \right)^{\frac{D}{2}}}\exp\left\lbrack - \frac{\left( \mathbf{\xi} - \mathbf{u} \right)^{2}}{2R_{g}T} \right\rbrack,
\tag{3.11}
$$

$$
{\overline{h}}^{(eq)} = \int_{}^{}{\frac{\mathbf{\xi}^{2}\mathbf{+}\mathbf{\eta}^{2}}{2}{\widetilde{f}}^{(eq)}}d\mathbf{\eta =}\frac{\rho\left( \mathbf{\xi}^{2} + KR_{g}T \right)}{2 \cdot \left( 2\pi R_{g}T \right)^{\frac{D}{2}}}\exp\left\lbrack - \frac{\left( \mathbf{\xi} - \mathbf{u} \right)^{2}}{2R_{g}T} \right\rbrack.
\tag{3.12}
$$

可以看到多原子分子的流场与单原子分子的区别主要体现在能量上, 而 ${\overline{f}}^{(eq)}$ 与式(1.10)一致.流场宏观量表示为：

$$
\rho = \int \overline{f}d\mathbf{\xi},\,
\mathbf{j} =\rho\mathbf{u} = \int {\mathbf{\xi}\overline{f}}d\mathbf{\xi},\,
\rho E = \int \overline{h}d\mathbf{\xi}.
$$

定容比热比 $c_{v}$ 和定压比热比 $c_{p}$ 与内部自由度K相关：

$$
c_{v} = \frac{D + K}{2}R_{g},\,
c_{p} = \frac{D + K + 2}{2}R_{g}.
$$

对单原子分子而言, $\mathbf{\eta = 0\ }$和$K = 0$.

式(3.9)和式(3.10)在$\mathbf{c}_{\mathbf{i}}$方向上的方程为

$$
\frac{\partial{\overline{f}}_{i}}{\partial t} + \mathbf{c}_{\mathbf{i}} \cdot \nabla{\overline{f}}_{i} = - \frac{1}{\tau_{f}}\left( {\overline{f}}_{i} - {\overline{f}}_{i}^{(eq)} \right) + F_{i},
\tag{3.13}
$$

$$
\frac{\partial{\overline{h}}_{i}}{\partial t} + \mathbf{c}_{\mathbf{i}} \cdot \nabla{\overline{h}}_{i} = - \frac{1}{\tau_{h}}\left( {\overline{h}}_{i} - {\overline{h}}_{i}^{(eq)} \right) + \frac{Z_{i}}{\tau_{hf}}\left( {\overline{f}}_{i} - {\overline{f}}_{i}^{(eq)} \right) + q_{i}.
\tag{3.14}
$$

其中 $Z_{i} = \mathbf{c}_{\mathbf{i}}\mathbf{\cdot u} - \frac{u^{2}}{2}$ .将其采用二阶积分格式离散后可得：

$$
{\overline{f}}_{i}\left( \mathbf{x} + \mathbf{c}_{i}\Delta t,t + \Delta t \right) - {\overline{f}}_{i}\left( \mathbf{x},t \right) = - \omega_{f}\left( {\overline{f}}_{i} - {\overline{f}}_{i}^{(eq)} \right) + \left( 1 - \frac{\omega_{f}}{2} \right)F_{i}\Delta t,
\tag{3.15}
$$

$$
{\overline{h}}_{i}\left( \mathbf{x} + \mathbf{c}_{i}\Delta t,t + \Delta t \right) - {\overline{h}}_{i}\left( \mathbf{x},t \right) = - \Omega_{h}\left( {\overline{h}}_{i} - {\overline{h}}_{i}^{(eq)} \right) + \frac{Z_{i}}{\tau_{hf}}\left( {\overline{f}}_{i} - {\overline{f}}_{i}^{(eq)} \right) + q_{i}
.\tag{3.16}
$$

Guo等 [^10] 将平衡态使用Hermite展开对 ${\overline{f}}^{(eq)}$ 和 ${\overline{h}}^{(eq)}$ 进行离散, 从而使低阶速度矩不会因高阶项的舍去而产生影响.由于涉及局部温度 $T$ , 因此展开过程是基于参考温度 $T_0$ 进行的.离散平衡态为：

$$
{\overline{f}}_{i}^{(eq)} = \rho w_{i}\left\lbrack 1 + \frac{\mathbf{c}_{\mathbf{i}}\mathbf{\cdot u}}{R_{g}T_{0}} + \frac{1}{2}\left( \frac{\mathbf{c}_{\mathbf{i}}\mathbf{\cdot u}}{R_{g}T_{0}} \right)^{2} - \frac{u^{2}}{2R_{g}T_{0}} \right\rbrack,
\tag{3.17}
$$

$$
{\overline{h}}_{i}^{(eq)} = p_{0}w_{i}\left\lbrack \frac{\mathbf{c}_{\mathbf{i}}\mathbf{\cdot u}}{R_{g}T_{0}} + \left( \frac{\mathbf{c}_{\mathbf{i}}\mathbf{\cdot u}}{R_{g}T_{0}} \right)^{2} - \frac{u^{2}}{R_{g}T_{0}} + \frac{1}{2}\left( \frac{{\mathbf{c}_{\mathbf{i}}}^{2}}{R_{g}T_{0}} - D \right) \right\rbrack + {\overline{f}}_{i}^{(eq)}E.\tag{3.18}
$$

源项为：

$$
F_{i} = \rho w_{i}\left\lbrack \frac{\mathbf{c}_{\mathbf{i}}\mathbf{\cdot}\mathbf{a}}{R_{g}T_{0}} + \frac{\left( \mathbf{c}_{\mathbf{i}}\mathbf{\cdot u} \right)\left( \mathbf{c}_{\mathbf{i}}\mathbf{\cdot}\mathbf{a} \right)}{\left( R_{g}T_{0} \right)^{2}} - \frac{\mathbf{a}\mathbf{\cdot u}}{R_{g}T_{0}} \right\rbrack,\,
q_{i} = \left( \frac{\rho w_{i}E}{R_{g}T_{0}} + f_{i} \right)\left( \mathbf{c}_{\mathbf{i}}\mathbf{\cdot}\mathbf{a} \right).
\tag{3.19}
$$

各宏观量则为：

$$
\rho = \sum_{i}^{}{\overline{f}}_{i},\,
\mathbf{j} = \rho\mathbf{u} = \sum_{i}^{}{\mathbf{c}_{i}{\overline{f}}_{i}},\, 
\rho E = \sum_{i}^{}{\overline{h}}_{i}.
$$

运动粘度为$\mu = \tau_{f}p_{0}$, 热扩散系数为 $\kappa = c_{v}\tau_{h}p_{0} = \frac{D + K}{2}R_{g}\tau_{h}p_{0}$ .因此该模型的Prandtl数为 $\Pr = \mu c_{p}/\kappa = \gamma\tau_{f}/\tau_{h}$ , 其中 $\gamma$ 为比热比.

> **[NOTE]**
> 对应的参考压强 $p_{0} = \rho R_{g}T_{0}$.

### 耦合DDF-LBE

何雅玲团队在Guo等 [^10] 的理论基础上, 提出了一种耦合的DDF-LBE[^25], 实现了热流和流场的耦合模拟.这里主要介绍耦合DDF-LBE与基于总能的DDF-LBE的区别.

为简化讨论, 考虑无外力项耦合DDF-LBE在 $\mathbf{c}_{i}$ 方向上的演化方程组：

$$
\frac{\partial f_{i}}{\partial t} + \mathbf{c}_{\mathbf{i}} \cdot \nabla f_{i} = - \frac{1}{\tau_{f}}\left( f_{i} - f_{i}^{(eq)} \right),\tag{3.20}
$$

$$
\frac{\partial h_{i}}{\partial t} + \mathbf{c}_{\mathbf{i}} \cdot \nabla h_{i} = - \frac{1}{\tau_{h}}\left( h_{i} - h_{i}^{(eq)} \right) + \frac{\left( \mathbf{c}_{\mathbf{i}}\mathbf{\cdot u} \right)}{\tau_{hf}}\left( f_{i} - f_{i}^{(eq)} \right).\tag{3.21}
$$

式(3.21)中将$Z_{i} = \mathbf{c}_{\mathbf{i}}\mathbf{\cdot u} - \frac{u^{2}}{2}$换为$\left( \mathbf{c}_{\mathbf{i}}\mathbf{\cdot u} \right)$并不会影响Navier-Stokes方程中的能量表示.

式(3.20)和式(3.21)中的两个平衡态分布($f_{\alpha}^{(eq)}$和$h_{\alpha}^{(eq)}$)可以用两类方式进行求解.第一种方式是将其写作式(3.22)和式(3.23)的截断形式：

$$
\begin{aligned}
f_{\alpha}^{(eq)} = & A_{i} + B_{\sigma} \cdot \left( \mathbf{c}_{\alpha} \cdot \mathbf{u} \right) + C_{\alpha}\left\| \mathbf{u} \right\|^{2} + D_{\alpha} \cdot \left( \mathbf{c}_{\alpha} \cdot \mathbf{u} \right)^{2} \\
 & + E_{\alpha}\left\| \mathbf{u} \right\|^{2} \cdot \left( \mathbf{c}_{\alpha} \cdot \mathbf{u} \right) + F_{\alpha} \cdot \left( \mathbf{c}_{\alpha} \cdot \mathbf{u} \right)^{3}
\end{aligned},
\tag{3.22}
$$

$$
h_{\alpha}^{(eq)} = K_{\alpha} + L_{\alpha} \cdot \left( \mathbf{c}_{\alpha} \cdot \mathbf{u} \right) + M_{\alpha}\left\| \mathbf{u} \right\|^{2} + N_{\alpha} \cdot \left( \mathbf{c}_{\alpha} \cdot \mathbf{u} \right)^{2}.
\tag{3.23}
$$

则各系数可通过式(3.11)和式(3.12)的Chapman-Enskog展开得出[^25].第二种方式则采用其他$f^{(eq)}$及其基于拉格朗日多项式的离散形式 $f_{\alpha}^{(eq)}$ .令

$$
h^{(eq)} = \left\lbrack E + \left( \mathbf{\xi} - \mathbf{u} \right) \cdot \mathbf{u} \right\rbrack f^{(eq)} + \omega(\mathbf{\xi},T)
$$

其中 $\omega\left( \mathbf{\xi},T \right)$ 满足

$$
\int {\omega(\mathbf{\xi},T)}d\mathbf{\xi =}0,\, 
\int {\omega(\mathbf{\xi},T)\mathbf{\xi}}d\mathbf{\xi =}0,\, 
\int {\omega\left( \mathbf{\xi},T \right)\mathbf{\xi}\mathbf{\xi}^{T}}d\mathbf{\xi =}pR_{g}T\mathbf{I}
$$

因此可通过构建等价离散形式(式(3.24))：

$$
h_{\alpha}^{(eq)} = \frac{w_{\alpha}p}{c^{2}}R_{g}T + \left\lbrack E + \left( \mathbf{c}_{\alpha} - \mathbf{u} \right) \cdot \mathbf{u} \right\rbrack \cdot f_{\alpha}^{(eq)}
\tag{3.24}
$$

计算总能的离散平衡态.对于二维问题, 基于Qu等提出的圆函数[^26]进行修改：

$$
f^{(eq)} = \begin{cases}
\frac{\rho}{2\pi c} & \parallel \mathbf{\xi} - \mathbf{u} \parallel = c \equiv \sqrt{D(\gamma - 1)e} \text{and} \lambda = e_{p}, \\
0 & \text{otherwise}
\end{cases}
\tag{3.25}
$$

其中 $\lambda$ 为静能,  $e_{p} = \left\lbrack 1 - \frac{D}{2}(\gamma - 1) \right\rbrack e$ .对于三维问题, Li等[^27]采用球面函数表示 $f^{(eq)}$ ：

$$
f^{(eq)} = \begin{cases}
\frac{\rho}{4\pi r^{2}} & \left\| \xi - \mathbf{u} \right\| = \left\| \mathbf{r} \right\| = r, \\
0 & \text{otherwise}. \\
\end{cases}
\tag{3.26}
$$

以Li等[^27]使用的D3Q27模型为例, 可得$\rho r^{2}/3 = p$, 因此对于理想气体有$r = \sqrt{3R_{g}T}$.

此外, Li等[^25]指出：空间离散采用五阶WENO和TVD(total variation diminishing scheme)格式捕捉可压缩流中的不连续性；时间离散可采用二阶(或三阶)IMEX Runge-Kutta格式推进.虽然式(3.20)和式(3.21)同样可被离散为如下形式($t$和$t + \Delta t$和上标表示时间)：

$$
f_{\alpha}^{t + \Delta t} = f_{\alpha}^{t} - \Delta t \cdot \left( \mathbf{c}_{\alpha} \cdot \mathbf{\nabla}f_{\alpha}^{t} \right) + \frac{\Delta t}{\tau_{f}^{t}}\left( f_{\alpha}^{(eq),t} - f_{\alpha}^{t} \right)
\tag{3.27}
$$

$$
h_{\alpha}^{t + \Delta t} = h_{\alpha}^{t} + \Delta t\left\lbrack \frac{1}{\tau_{h}^{t}}\left( h_{\alpha}^{(eq),t} - h_{\alpha}^{t} \right) - \frac{\left( \mathbf{c}_{\alpha} \cdot \mathbf{u}^{t} \right)}{\tau_{hf}^{t}}\left( f_{\alpha}^{(eq),t} - f_{\alpha}^{t} \right) - \left( \mathbf{c}_{\alpha} \cdot \mathbf{\nabla}h_{\alpha}^{t} \right) \right\rbrack
\tag{3.28}
$$

但这种简单的形式仅适用于流场中无激波、无间断的情况, 且在时间上仅有一阶精度.

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


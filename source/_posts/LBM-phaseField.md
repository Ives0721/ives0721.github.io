---
title: 【笔记】格子Boltzmann方法的相场模型
tags:
  - 格子Boltzmann方法
  - 相场模型
categories:
  - - 格子Boltzmann方法
  - - 相场模型
date: 2026-04-03 10:28:49
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
    .note-block {
        margin: 20px;
        padding: 15px;
        border: 3px dashed rgba(94, 94, 255, 1);
        background-color: rgba(242, 242, 255, 0.5);
        border-radius: 15px;
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        font-weight: 500;
        transition: all 0.3s ease;
    }
</style>

![](https://img.shields.io/badge/格子Boltzmann方法-PhaseField-green)

最近自己也是在学相场格子Boltzmann方法，这篇笔记主要是在整理刚学到的内容。

# 控制方程

## 相场和表面张力

相场采用局部守恒性 Allen-Cahn 方程描述[1,2]，其核心量为序参数 $\phi$ 。控制方程为：

$$
\frac{\partial \phi}{\partial t} + \nabla \cdot (\phi \boldsymbol{u}) = 
M_{\phi} (\nabla^2\phi - \nabla\cdot(\lambda\boldsymbol{n}))
$$

其中
* $M_{\phi}$ 为迁移率 (mobility)。
* $\boldsymbol{n}=\nabla\phi/\|\nabla\phi\|$ 是界面法向量。
* $\lambda = 4(1-\phi)\phi/W$。 $W$ 为界面厚度，通常可以取 5 个格子单位。

表面张力通常采用如下形式：

$$
\boldsymbol{F}_s = \mu_\phi \nabla \phi
$$

化学势 $\mu_\phi$ 定义为：

$$
\mu_\phi = 4\beta\phi(\phi-1)(\phi-0.5) - k \nabla^2\phi.
$$

$k$ 和 $\beta$ 取决于界面厚度 $W$ 和表面张力 $\sigma$：

$$
k = 3 \sigma W /2, \quad
\beta = 12 \sigma / W
$$

> <b>注：</b> Open Access的[文献[1]](https://www.sciopen.com/article/10.26804/capi.2019.03.01)非常详尽地介绍了这一部分的推导，这里不做赘述。

## 流场

流场的控制方程组分为两个式子，这里我们简单讨论两种不可压流体的情况。

连续性方程为：

$$
\nabla \cdot \boldsymbol{u} = 0,
$$

Navier-Stokes方程为：

$$
\begin{aligned}
    \frac{\partial (\rho\boldsymbol{u})}{\partial t} &+ \nabla \cdot (\rho \boldsymbol{u} \boldsymbol{u}) = - \nabla p\\
    & + \nabla \cdot [\mu (\nabla\boldsymbol{u} + (\nabla\boldsymbol{u})^\mathrm{T})] + \mathbb{F}
\end{aligned}
$$

其中
* $\rho, \boldsymbol{u}, p$ 分别为流场密度、速度和压强。
* $\mu, \nu$ 分别为流体的动力粘度和运动粘度。 $\mu = \rho\nu$ 。
* 外力 $\mathbb{F} = \boldsymbol{F}_s + \boldsymbol{G}$ 。 $\boldsymbol{F}_s$ 为表面张力， $\boldsymbol{G}$ 为其他体力 (如重力)。


# 描述两相流的双分布函数

## 相场

分布函数 $f_i$ 用于描述相场，其演化方程为：

$$
\begin{aligned}
f_i (\boldsymbol{x} &+ \boldsymbol{c}_i \delta_t, t + \delta_t) - f_i (\boldsymbol{x}, t) =\\ 
&-\frac{1}{\tau_f}
\left[f_i (\boldsymbol{x}, t) - f_i^{eq} (\boldsymbol{x}, t)\right] + 
F_i (\boldsymbol{x}, t) \delta_t 
\end{aligned}
\tag{1}
$$

其中，相场中的平衡态 $f_i^{eq}$ 用的是1阶截断的形式：

$$
f_i^{eq} = \omega_i \phi \left( 1 + \frac{\boldsymbol{c}_i \cdot \boldsymbol{u}}{c_s^2} \right).
$$

式(1)里的作用力项 $F_i$ 为：

$$
F_i = \left( 1 - \frac{1}{2 \tau_f} \right)
\frac{w_i}{c_s^2}
\boldsymbol{c}_i \cdot \left[
    \frac{\partial (\phi \boldsymbol{u})}{\partial t} + c_s^2 \lambda \boldsymbol{n}
\right].
$$

松弛时间 $\tau_f$ 由迁移率 $M_{\phi}$ 计算：

$$
M_{\phi} = (\tau_f - \frac{1}{2}) c_s^2 \delta_t
$$

序参数 (order parameter) $\phi$ 是用于描述多相流的等效流体密度 $\rho$ 的参数。若记液相和气相密度分别为 $\rho_l$ 和 $\rho_g$，则：

$$
\rho = \phi \rho_l + (1 - \phi) \rho_g .
$$

而式(1)的模型中， $\phi$ 可表示为： $\displaystyle\phi=\sum_{i}f_{i}$ .
 
<!-- ![相场模型的粘度函数](LBM_相场粘度函数示意图.png) -->
<img src="LBM_相场粘度函数示意图.png">

混合流体的运动粘度 $\nu$ 同样可以由 $\phi$ 进行分配，如上图所示。具体的计算方式包括：

（1）阶跃函数 $\displaystyle\nu(\phi)=\begin{cases}\nu_g,&\phi < 0.5\\ \nu_l,&\phi \ge 0.5\end{cases}$ 。 其中 $\nu_g$ 和 $\nu_l$ 分别为气相和液相的运动粘度。

（2）线性函数 $\displaystyle\nu(\phi)=\phi\nu_{l}+(1-\phi)\nu_{g}$。

（3）倒数的线性形式：$\displaystyle\frac{1}{\nu(\phi)}=\frac{\phi}{\nu_{l}}+\frac{1-\phi}{\nu_{g}}$。这种形式具有二阶可导的特点。

## 流场

分布函数 $g_i$ 用于描述流场，其演化方程为：

$$
\begin{aligned}
g_i (\boldsymbol{x} &+ \boldsymbol{c}_i \delta_t, t + \delta_t) - g_i (\boldsymbol{x}, t) =\\ 
&-\frac{1}{\tau_g}
\left[g_i (\boldsymbol{x}, t) - g_i^{eq} (\boldsymbol{x}, t)\right] + 
G_i (\boldsymbol{x}, t) \cdot \delta_t 
\end{aligned}
\tag{2}
$$

局部平衡态 $g_i^{eq}$ 需要使用下面的版本，以保证速度场的无散性。

$$
g_i^{eq} = 
\begin{cases}
    \frac{p}{c_s^2} (\omega_i - 1) + \rho s_{i}(\boldsymbol{u}), & i = 0\\
    \frac{p}{c_s^2} \omega_i + \rho s_{i}(\boldsymbol{u}), & i \neq 0\\
\end{cases}
\tag{3}
$$

其中
$$
s_{i}(\boldsymbol{u}) = \left[
    \frac{\boldsymbol{c}_i \cdot \boldsymbol{u}}{c_s^2} + 
    \frac{(\boldsymbol{c}_i \cdot \boldsymbol{u})^2}{2 c_s^4} - 
    \frac{\boldsymbol{u} \cdot \boldsymbol{u}}{2 c_s^2}
\right] \cdot \omega_i
$$

源项为

$$
\begin{aligned}
G_i (\boldsymbol{x}, t) &=
\left(1 - \frac{1}{2 \tau_g}\right) \omega_i
\left[
    \boldsymbol{u} \cdot \nabla \rho + 
    \frac{\boldsymbol{c}_i \cdot \mathbb{F}}{c_s^2} + 
    \frac{\boldsymbol{u} \nabla \rho: (\boldsymbol{c}_i \boldsymbol{c}_i - c_s^2 \boldsymbol{I})}{c_s^2}
\right]\\
&=
\left(1 - \frac{1}{2 \tau_g}\right) \omega_i
\left[
    \frac{\boldsymbol{c}_i \cdot \mathbb{F}}{c_s^2} + 
    \frac{(\rho_l - \rho_g) \boldsymbol{u} \nabla\phi: (\boldsymbol{c}_i \boldsymbol{c}_i)}{c_s^2}
\right]
\end{aligned}
\tag{4}
$$

式(4)用到了 $\nabla \rho = (\rho_l - \rho_g) \nabla\phi$ 进行化简。


所以，宏观量的计算为：

$$
\begin{aligned}
\rho \boldsymbol{u} &= \frac{\delta_t \mathbb{F}}{2} + \sum_{i} \boldsymbol{c}_i g_i ,\\
p &= \frac{c_s^2}{1 - \omega_0} \left[ \sum_{i \neq 0}{g_i} + \frac{\delta_t}{2} \boldsymbol{u} \cdot \nabla \rho + \rho s_{0}(\boldsymbol{u}) \right]\\
&= \frac{c_s^2}{1 - \omega_0} \left[ \sum_{i \neq 0}{g_i} + \frac{\delta_t}{2} (\rho_l - \rho_g)\boldsymbol{u} \cdot \nabla\phi + \rho s_{0}(\boldsymbol{u}) \right].
\end{aligned}
\tag{5}
$$

## 相场梯度

记 $\displaystyle\mathcal{A}=-c_s^2 \tau_f \delta_t$, $\displaystyle\mathcal{B} = M_{\phi} \lambda \delta_t$, 和 $\displaystyle\mathcal{C} = \frac{\delta_t}{2}\phi\boldsymbol{u} + \sum_{i}{\boldsymbol{c}_i (f_i - f_i^{eq})}$，
则：

$$
\|\nabla\phi\| = \frac{-\|\mathcal{C}\|-\mathcal{B}}{\mathcal{A}} ,
\quad
\nabla\phi = \frac{\mathcal{C}}{\mathcal{A} + \mathcal{B}/\|\nabla\phi\|}.
\tag{6}
$$

化学势 $\mu_\phi$ 中 $\nabla^2\phi$ 部分则用 $\nabla\phi$ 的差分进行计算。

当然，也有完全使用差分的计算公式
$$
\begin{cases}
    \nabla\phi(\boldsymbol{x}) &= \sum_{i \neq 0} {\omega_i \boldsymbol{c}_i \phi(\boldsymbol{x} + \boldsymbol{c}_i \delta_t)} / {c_s^2 \delta_t} ,\\
    \nabla^2\phi(\boldsymbol{x}) &= \sum_{i \neq 0} {2\omega_i [\phi(\boldsymbol{x} + \boldsymbol{c}_i \delta_t) - \phi(\boldsymbol{x})]} / {c_s^2 \delta_t^2} 
\end{cases}
$$


# 绕过化学势的相场模型

值得一提的是，Tim Reis [[3]](https://timreis.org/wp-content/uploads/2022/01/preprint.pdf) 提出了一种不需要计算 $\nabla^2\phi$ 的相场LBM模型。该方法主要是修改了流场分布函数 $g_i$ 的平衡态和源项，所以这里只讲与流场相关的部分（相场部分可以看他的文章）。

定义表面张力为 $\boldsymbol{F}_s = \nabla\cdot\mathrm{T}$，其中毛细张量（capillary tensor） $\mathrm{T}$ 写作：

$$
\mathrm{T}=\sigma (\boldsymbol{I} - \boldsymbol{n}\boldsymbol{n}) \delta_{\mathrm{S}}
$$

而 $\delta_{\mathrm{S}}\simeq\|\nabla\phi\|$ 是界面delta函数（surface delta function）。

Tim Reis 基于多尺度展开结果，把 $\mathrm{T}$ 的作用写进修正的平衡态分布函数里，

$$
g_{i}^{eq}=\omega_i \left[ \frac{p}{c_s^2} + \rho \frac{\boldsymbol{c}_i \cdot \boldsymbol{u}}{c_s^2} + \frac{1}{2c_s^4}(\rho\boldsymbol{u}\boldsymbol{u}+\mathrm{T}):(\boldsymbol{c}_i \boldsymbol{c}_i - c_s^2 \boldsymbol{I}) \right],
$$

这使得化学势 $\mu_\phi$ 的计算被约去了，且外力项也仅剩一个重力 $\boldsymbol{G}$ 。在这种构造下，源项 $G_i$ 也要被分为两部分

$$
G_i = \left(1 - \frac{1}{2 \tau_g}\right) [R_i(\boldsymbol{G}) + S_i]
$$

其中 $R_i(\boldsymbol{G})$ 为外力项

$$
R_i = \omega_i \left(
    \frac{\boldsymbol{c}_i - \boldsymbol{u}}{c_s^2} + 
    \frac{\boldsymbol{c}_i \cdot \boldsymbol{u}}{c_s^4} \boldsymbol{c}_i
\right) \cdot \boldsymbol{G},
$$

而 $S_i$ 为修正项

$$
\begin{cases}
    S_i = (\Gamma_i - \omega_i) (\boldsymbol{c}_i - \boldsymbol{u}) \cdot \nabla\rho \\
    \Gamma_i =\left(
        1 + \frac{\boldsymbol{c}_i - \boldsymbol{u}}{c_s^2} + 
        \frac{\boldsymbol{u}\boldsymbol{u}:(\boldsymbol{c}_i \boldsymbol{c}_i - c_s^2 \boldsymbol{I})}{2c_s^4}
    \right) \omega_i
\end{cases}
$$

宏观量的计算也要变成

$$
\begin{cases}
    \sum_i g_i = \frac{p}{c_s^2} - \frac{\delta_t}{2} (\boldsymbol{u} \cdot \nabla)\rho \\
    \sum_i g_i \boldsymbol{c}_i = \rho \boldsymbol{u} + \frac{\delta_t}{2} \boldsymbol{G}
\end{cases}.
$$

# 参考文献

<div class="note-block">
<b>
带有 '※' 标记的文章意味着可通过免费渠道获得，如：开放获取 (Open Access)；Research gate；arXiv预印本；Sci hub；作者主页等。
</b>
</div>

[1] ※ WANG H, YUAN X, LIANG H, _et al_. **A brief review of the phase-field-based lattice Boltzmann method for multiphase flows**[J]. Capillarity, 2019, 2(3): 33-52. [DOI:10.26804/capi.2019.03.01](https://www.sciopen.com/article/10.26804/capi.2019.03.01).

[2] ※ LIANG H, XU J, CHEN J, _et al_. **Phase-field-based lattice Boltzmann modeling of large-density-ratio two-phase flows**[J]. Physical Review E, 2018, 97(3): 033309. [DOI:10.1103/PhysRevE.97.033309](https://arxiv.org/abs/1710.09541).

[3] ※ REIS T. **A lattice Boltzmann formulation of the one-fluid model for multiphase flow**[J]. Journal of Computational Physics, 2022, 453: 110962. [DOI:10.1016/j.jcp.2022.110962](https://timreis.org/wp-content/uploads/2022/01/preprint.pdf).


<!-- https://blog.csdn.net/qq_32515081/article/details/125967160 -->
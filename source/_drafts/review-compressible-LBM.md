---
title: 【综述】格子Boltzmann方法对可压缩流体和热流的模拟
date: 2024-10-16 18:22:55
tags: [格子Boltzmann方法]
categories:
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

# 格子 Boltzmann 方程简述

TODO

# 可压缩流体的 LBM 模拟

最基本的 LBM 模型通常适用于**等温不可压缩流体**（或弱可压缩）在低马赫数下的数值模拟，并且其普朗特数 $\mathrm{Pr}=1$。 这些限制使得最简单的 LBM 模型难以满足热流体以及超声速流体模拟的要求。因此，研究者们提出了多种 LBM 的改进方式，以实现更多功能，如：可调节的 $\mathrm{Pr}$、高马赫数、流场温度 $T$ 的计算等。

下文为方便起见，介绍的所有 LBM 模型均为二维形式，且不包含外力项的构建和计算。

## 多速 LBM 模型

多速 LBM （_Multi-Speed (MS) LBM_）是一种完全基于 LBM 经典算法的改进，它通过使用比常规 LBM 更为复杂的离散速度模型还原高阶速度矩，在保持算法不变的基础上实现更精确的模拟。

### MS LBM 的常规形式

常规形式的 MS-LBM 演化方程与常规 LBM 一致，Chen等[^MSLBM_chen_1994]给出平衡态 $f_{i}^{(eq)}$ 的表示为

$$
\begin{aligned}
  f_{i}^{(eq)} &= A + B \boldsymbol{c}_i \cdot \boldsymbol{u} + (\boldsymbol{c}_i \cdot \boldsymbol{u})^2 C 
  \\&\quad + D \boldsymbol{u}^2 + (\boldsymbol{c}_i \cdot \boldsymbol{u})^3 E + F (\boldsymbol{c}_i \cdot \boldsymbol{u}) \boldsymbol{u}^2
\end{aligned}
$$

其中 $A,B,C,D,E,F$ 是取决于 DnQb 离散速度模型的参数，它们是密度 $\rho$ 和内能 $e$ 的函数。 平衡态的速度矩需要满足下列关系[^MSLBM_chen_1994]：

$$
\begin{aligned}
  \sum_{i} f_{i}^{(eq)} &= \rho \\
  \sum_{i} c_{i\alpha} f_{i}^{(eq)} &= \rho {u}_{\alpha} \\
  \sum_{i} c_{i\alpha} c_{i\beta} f_{i}^{(eq)} &= \frac{2}{D} \rho e \delta_{\alpha\beta} + \rho u_{\alpha} u_{\beta} \\
  \sum_{i} c_{i\alpha} c_{i\beta} c_{i\gamma} f_{i}^{(eq)} &= \rho u_{\alpha} u_{\beta} u_{\gamma} + 
  \\&\quad \frac{2}{D} (u_{\alpha} \delta_{\beta\gamma} + u_{\gamma} \delta_{\alpha\gamma} + u_{\gamma} \delta_{\alpha\beta}) \rho e \\
  \sum_{i} c_{i\alpha} c_{i\beta} \boldsymbol{c}_i^2 f_{i}^{(eq)} &= \frac{4(D+2)}{D} \rho e \delta_{\alpha\beta}
  \\&\quad + \frac{2}{D} \rho e \boldsymbol{u}^2 \delta_{\alpha\beta} + \frac{2(D+4)}{D} \rho e u_\alpha u_\beta 
  \\&\quad + \rho \boldsymbol{u}^2 u_\alpha u_\beta 
\end{aligned}
$$

其中 $D$ 为空间维度数目。流体的剪切粘度 $\mu$ 、体积粘度 $\lambda$ 、导热系数 $\kappa$ 、比热比 $\gamma$ 和 LBM 松弛时间 $\tau$ 的关系为：

$$
\begin{aligned}
  \mu &= \frac{2}{D} (\tau - \frac{1}{2}) \rho e \\
\lambda &= -\frac{4}{D^2} (\tau - \frac{1}{2}) \rho e\\
\kappa &= \frac{2+D}{D} (\tau - \frac{1}{2}) \rho e \\
\gamma &= \frac{2+D}{D}.
\end{aligned}
$$

该模型的热力学压强为 $p=\frac{2}{D} \rho e$。由于 Prandtl 数为 $\mathrm{Pr} = \frac{\mu}{\kappa}=\frac{2}{D+2}$。因此这一模型仅能模拟热量变化较小的热流。

### Kataoka-Tsutahara 模型

TODO

### 多速 LBM 模型的熵形式

TODO

## 双分布函数 LBM

TODO


# 参考文献

[^MSLBM_chen_1994]: Chen, Y., Ohashi, H., & Akiyama, M. (1994). Thermal lattice Bhatnagar-Gross-Krook model without nonlinear deviations in macrodynamic equations. Phys. Rev. E, 50, 2776–2783.


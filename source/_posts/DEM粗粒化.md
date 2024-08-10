---
title: 离散单元法的颗粒数据粗粒化
date: 2024-07-30 14:43:10
tags: [离散单元法, 粗粒化]
categories:
- [离散单元法]
---
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex/dist/katex.min.css">

离散单元法（**D**iscrete **E**lement **M**ethod，简称**DEM**）是颗粒流数值模拟的重要数值方法，被用于各种各样的颗粒流定量研究中。颗粒流的粗粒化（**Coarse Graining**），便是将颗粒流流场的大量数据转化为连续场数据的数值方法。

在开始介绍之前，先简要列出主要的参考资料：  
- Breard, E. C. P., Dufek, J., Fullard, L., & Carrara, A. (2020). <strong>The Basal Friction Coefficient of Granular Flows With and Without Excess Pore Pressure: Implications for Pyroclastic Density Currents, Water‐Rich Debris Flows, and Rock and Submarine Avalanches</strong>. _Journal of Geophysical Research: Solid Earth_, 125(12), e2020JB020203. [DOI:10.1029/2020JB020203](https://doi.org/10.1029/2020JB020203).
- Lacaze, L., & Kerswell, R. R. (2009). <strong>Axisymmetric Granular Collapse: A Transient 3D Flow Test of Viscoplasticity</strong>. _Physical Review Letters_, 102(10), 108305. [DOI:10.1103/PhysRevLett.102.108305](https://doi.org/10.1103/PhysRevLett.102.108305).
- Weinhart, T., Labra, C., Luding, S., & Ooi, J. Y. (2016). <strong>Influence of coarse‐graining parameters on the analysis of DEM simulations of silo flow</strong>. _Powder Technology_, 293, 138–148. [DOI:10.1016/j.powtec.2015.11.052](https://doi.org/10.1016/j.powtec.2015.11.052)

# 单种颗粒的粗粒化

为简单起见，先从仅存在单种球形颗粒的颗粒流开始说起，并且只考虑获取空间内一点 $\boldsymbol{r}$ 在 $t$ 时刻的粗粒化数据。

## 粗粒化函数

在这里，引入高斯函数

$$
W(\boldsymbol{r}) = \begin{cases}
\frac{1}{V_w} \exp{\left(- \frac{\vert\boldsymbol{r}\vert^2}{2w} \right)} ,& 
{\vert\boldsymbol{r}\vert} < c \\
0 ,& \text{otherwise}
\end{cases}
$$

其中： $w$ 为粗粒化宽度； $c = 3 w$ 为粗粒化截断长度； $V_w$ 是确保其密度积分等于总质量的系数

$$
V_w = 2 \sqrt{2} \cdot w^3 \pi^{\frac{3}{2}} \mathrm{erf}\left( \frac{c \sqrt{2}}{2 w} \right) - 4 c w^2 \pi \exp{\left( \frac{- c^2}{2 w^2} \right)}
$$

Weinhart 等建议 $w$ 的取值范围应在 $[0.75 d_{\mathrm{mean}}, 1.25 d_{\mathrm{mean}}]$ 之间，其中 $d_{\mathrm{mean}}$ 为平均颗粒直径。

## 密度和动量的计算

在 $t$ 时刻，单种颗粒的集合为 $q$ 。对于其中的第 $i$ 个球，其位置矢量为 $\boldsymbol{r}_i$ ，质量为 $m_i$ ，速度为 $\boldsymbol{u}_i$。

密度的粗粒化表达式为：
$$
\rho(\boldsymbol{r}, t) = \sum\limits_{i \in q} m_i W(\boldsymbol{r} - \boldsymbol{r}_{i}(t))
$$

动量的粗粒化表达式为：
$$
\boldsymbol{j}(\boldsymbol{r}, t) = \sum\limits_{i \in q} m_i \boldsymbol{u}_i W(\boldsymbol{r} - \boldsymbol{r}_{i}(t))
$$

所以，速度的计算需要通过 $\boldsymbol{j}$ 和 $\rho$ 完成：

$$
\boldsymbol{u}(\boldsymbol{r}, t) = \boldsymbol{j}(\boldsymbol{r}, t) / \rho(\boldsymbol{r}, t)
$$

{% note %}  
根据 $W(\boldsymbol{r})$ 的定义，当 $\vert \boldsymbol{r} - \boldsymbol{r}_{i}(t) \vert \ge c$ 时， $W(\boldsymbol{r} - \boldsymbol{r}_{i}(t)) = 0$。  
因此在实际计算中可根据这一点跳过某些球的计算。  
{% endnote %}

## 应力张量的计算

应力张量 $[\boldsymbol{\sigma}] = \sigma_{\alpha\beta}$ 被表示为由接触作用产生的部分 $\sigma_{\alpha\beta}^{(c)}$ ，以及由动能作用产生的部分 $\sigma_{\alpha\beta}^{(k)}$ 。即：

$$
\sigma_{\alpha\beta} = \sigma_{\alpha\beta}^{(k)} + \sigma_{\alpha\beta}^{(c)}
$$

对于 $\sigma_{\alpha\beta}^{(k)}$ ，有：

$$
\sigma_{\alpha\beta}^{(k)}(\boldsymbol{r}, t) = \sum\limits_{i \in q} m_i v_{i,\alpha}^{'} v_{i,\beta}^{'} W(\boldsymbol{r} - \boldsymbol{r}_{i}(t))
$$

其中 $v_{i,\alpha}^{'} = u_{i,\alpha} - u_{\alpha}(\boldsymbol{r}, t)$ 为第 $i$ 个球的速度波动值。

对于 $\sigma_{\alpha\beta}^{(c)}$ ，考虑第 $i$ 个球的所有接触：

$$
\sigma_{\alpha\beta}^{(c)}(\boldsymbol{r}, t) = \sum\limits_{i \in q} \sum\limits_{j \in \bar{Q}}^{j \neq i} F_{ij,\alpha} l_{ij,\beta} \int_{0}^{1} W(\boldsymbol{r} - \boldsymbol{r}_{i}(t) + s \boldsymbol{l}_{ij}) \mathrm{d}s
$$

式中 
- $j$ 是整个颗粒系统 $\bar{Q}$ 中的一个球（$j \neq i$），并且球 $i$ 与球 $j$ 发生接触；
- $\boldsymbol{F}_{ij}$ 是球 $i$ 所受到的来自球 $j$ 接触的力， $F_{ij,\alpha}$ 为其 $\alpha$ 方向分量；
- $\boldsymbol{l}_{ij}$ 是球 $i$ 与球 $j$ 的接触矢量（从 $i$ 向 $j$）， $l_{ij,\alpha}$ 为其 $\alpha$ 方向分量。

{% note %}  
在实际操作中，对于
$$\int_{0}^{1} W(\boldsymbol{r} - \boldsymbol{r}_{i}(t) + s \boldsymbol{l}_{ij}) \mathrm{d}s$$
我选择使用**复化 Simpson 积分**计算其数值，以简化操作。  
{% endnote %}


# 多种颗粒的粗粒化

假设一个颗粒流系统 $\bar{Q}$ 中有 $Q$ 类颗粒，对于第 $q$ 类颗粒 $q \in Q$， 可以使用单类颗粒粗粒化的公式计算出：第 $q$ 类颗粒的粗粒化密度 $\rho^{q}(\boldsymbol{r}, t)$ 、速度 $\boldsymbol{u}^{q}(\boldsymbol{r}, t)$ 和应力张量 $\sigma_{\alpha\beta}^{q}(\boldsymbol{r}, t)$。

所以整个颗粒流系统 $\bar{Q}$ 的粗粒化结果为：

$$
\begin{aligned}
  \rho(\boldsymbol{r}, t) &= \sum_{q \in Q} \rho^{q}(\boldsymbol{r}, t) \\
  \boldsymbol{u}(\boldsymbol{r}, t) &= \sum_{q \in Q} \boldsymbol{u}^{q}(\boldsymbol{r}, t) \\
  \sigma_{\alpha\beta}(\boldsymbol{r}, t) &= \sum_{q \in Q} \sigma_{\alpha\beta}^{q}(\boldsymbol{r}, t) \\
\end{aligned}
$$

# 使用粗粒化数据提取颗粒流流变参数

这里以 $\mu(I)$ 流变关系为例子，以上文的文献作为参照。

记速度剪切率张量 $\gamma_{ij} = \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}$ ，其中 $u_i$ 源自粗粒化速度场。令：

$$
\gamma_{ij}^{d} = \gamma_{ij} - \frac{\delta_{ij}}{3} \mathrm{tr}([\gamma])
$$

其中 $\mathrm{tr}([\gamma]) = \gamma_{xx} + \gamma_{yy} + \gamma_{zz}$ 为 $[\gamma]$ 的迹。

所以剪切应变率定义为：

$$
|\dot{\gamma}| = \sqrt{\frac{1}{2} {\gamma}_{ij}^{d} {\gamma}_{ij}^{d} }
$$

定义压强为 $p = \frac{1}{3} \mathrm{tr}([\sigma]) = (\sigma_{xx} + \sigma_{yy} + \sigma_{zz}) / 3$，偏应力张量为： $[\tau] = \tau_{ij} = \sigma_{ij} - p \cdot \delta_{ij}$， 所以切应力为：
$$
|\tau| = \sqrt{\frac{1}{2} \tau_{ij} \tau_{ij} }
$$

# 简单的粗粒化代码示意

我用 Python 写了一份对**单种颗粒系统**进行粗粒化的代码，详情见[此网页](/2024/08/10/DEM粗粒化-CODE)。
---
title: "格子Boltzmann方法的DnQb模型"
date: 2026-02-18 22:26:51
tags: [格子Boltzmann方法,DnQb]
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
</style>

这是自己的一份笔记，主要是和格子Boltzmann方法的DnQb离散速度模型相关的。还请读者确定自己已经知道常见的DnQb模型长什么样再来读这篇文章，这样会好理解一些。

# Gauss-Hermite积分与DnQb模型

一般来说，密度分布函数的Maxwell-Boltzmann平衡态分布为：
$$
f^{eq} = \frac{\rho}{ (2 \pi R T)^{D/2} } \exp \left[ - \frac{(\boldsymbol{\xi} - \boldsymbol{u})^2}{2 R T} \right],
\tag{1-1}
$$

其中 $\boldsymbol{\xi}$ 为分子速度；$R$ 为气体常数；$D$为问题的维度数；$\rho$ 为流场密度；$\boldsymbol{u}$ 为流场速度。而 LBM 里常用的离散平衡态分布一般写作（沿着 $\boldsymbol{c}_\alpha$ 方向，权重为 $W_\alpha$，$c_s^2 = RT$）：

$$
f^{eq}_{\alpha} = \rho W_{\alpha} \times 
\left[
    1 + \frac{\boldsymbol{\xi}_{\alpha} \cdot \boldsymbol{u}}{c_s^2} + 
    \frac{(\boldsymbol{\xi}_{\alpha} \cdot \boldsymbol{u})^2}{2 c_s^4} - 
    \frac{|\boldsymbol{u}|^2}{2 c_s^2}
\right].
\tag{1-2}
$$

那么这两者是怎么联系起来的呢？参考文献 [1] 讲的就是这个事情。最近自己总算是略懂一点了，这里把笔记写下来，作为备忘录。

## 平衡态分布的展开

当马赫数远小于1时，我们可以对式(1)的平衡态做如下离散：
$$
\begin{aligned}
    f^{eq} (\boldsymbol{\xi}) &= \frac{\rho}{ (2 \pi R T)^{D/2} } \exp \left( - \frac{|\boldsymbol{\xi}|^2}{2 R T} \right)\\
    & \quad \times 
    \left[
        1 + \frac{\boldsymbol{\xi} \cdot \boldsymbol{u}}{R T} + 
        \frac{(\boldsymbol{\xi} \cdot \boldsymbol{u})^2}{2 (R T)^2} - 
        \frac{|\boldsymbol{u}|^2}{2 R T}
    \right]
\end{aligned}
\tag{1-3}
$$

式(1-3)中，我们对 $\boldsymbol{u}$ 进行展开，并且只截断到2阶。而格子Boltzmann方法（LBM）要做的，便是把下面的这个积分
$$
\int_{-\infty}^{+\infty} \psi(\boldsymbol{\xi}) f^{eq} (\boldsymbol{\xi}) \mathrm{d}\boldsymbol{\xi}
$$

用离散的方式近似计算，其中 $\psi(\boldsymbol{\xi})$ 是关于 $\boldsymbol{\xi}$ 的多项式。下面我们将先以一维的情况为例，不去考虑多维的积分，简单讲讲数值离散的具体操作。

## 一维情形的离散 —— 以 D1Q3 为例

由于这一节只讨论一维问题，分子速度 $\xi$ 和宏观速度 $u$ 即是标量值（不需要向量形式的加粗），并且维度 $D=1$。所以，式(3)写作
$$
\begin{aligned}
    f^{eq} (\xi) &= \frac{\rho}{ (2 \pi R T)^{1/2} } \exp \left( - \frac{\xi^2}{2 R T} \right)\\
    & \quad \times 
    \left[
        1 + \frac{\xi u}{R T} + 
        \frac{(\xi u)^2}{2 (R T)^2} - 
        \frac{u^2}{2 R T}
    \right].
\end{aligned}
\tag{1-4}
$$

### 积分形式的抽象

式(1-4)的中括号内的四项可以用积分的加法拆成四个积分，并且每一项都可以看成如下形式：
$$
\text{Constant} \cdot \int_{-\infty}^{+\infty} \exp \left( - \frac{\xi^2}{2 R T} \right) \psi(\xi) \mathrm{d}\xi,
$$

其中 $\psi(\xi)$ 为关于 $\xi$ 的多项式。那么，我们令 $\psi_{m}(\xi) = \xi^{m}$ （$m$ 为自然数），将研究目标变为式(1-5)的通式，即：
$$
\frac{\rho}{ (2 \pi R T)^{1/2} } \cdot \int_{-\infty}^{+\infty} \exp \left( - \frac{\xi^2}{2 R T} \right) \psi_{m}(\xi) \mathrm{d}\xi.
\tag{1-5}
$$

为简化这个积分的计算，这里需要引入变量代换
$$
\zeta = \xi / \sqrt{2 R T},
$$

则
$$
\mathrm{d}\xi = \sqrt{2 R T} \mathrm{d}\zeta.
$$

代换后，有：
$$
\begin{aligned}
&
\frac{\rho}{ (2 \pi R T)^{1/2} } \cdot  \int_{-\infty}^{+\infty} \exp \left( - \zeta^2 \right) \psi_{m}(\sqrt{2 R T} \zeta) \sqrt{2 R T} \mathrm{d}\boldsymbol{\zeta}
\\=& \frac{\rho}{ (2 \pi R T)^{1/2} } \cdot (2 R T)^{(m+1)/2} \cdot I_m
\\=& \frac{\rho}{\sqrt{\pi}} \cdot (2 R T)^{m/2} \cdot I_m
\end{aligned}
$$

其中 $\displaystyle I_m =  \int_{-\infty}^{+\infty} \exp \left( - \zeta^2 \right) \zeta^m \mathrm{d}\boldsymbol{\zeta}$，且 $m$ 为自然数。

所以，下面的积分会被展开为
$$
\begin{aligned}
    I=& \int_{-\infty}^{+\infty} \psi_m(\xi) f^{eq} (\xi) \mathrm{d}\xi
    \\=& \frac{\rho}{\sqrt{\pi}} \left[
        \left(1 - \frac{u^2}{2 R T}\right) (2 R T)^{m/2} I_m + 
        \frac{u}{RT} (2 R T)^{(m+1)/2}  I_{m+1} \right. 
    \\&\qquad \left. + 
        \frac{u^2}{2(RT)^2} (2 R T)^{(m+2)/2} I_{m+2}
    \right].
\end{aligned}
\tag{1-6}
$$

### 数值积分的构造

假设我们想用一个[三点格式的Gauss-Hermite数值积分](https://mathworld.wolfram.com/Hermite-GaussQuadrature.html)进行近似，则：
$$
\int_{-\infty}^{+\infty} \exp(- \zeta^2) \phi(\zeta) \mathrm{d} \zeta = \sum_{i \in [-1, 0, 1]} w_{i} \phi(\zeta_{i}).
$$

其中：
$$
\begin{aligned}
\zeta_0 = 0 ,\quad& \zeta_{1} = \sqrt{3/2} ,\quad& \zeta_{-1} = -\sqrt{3/2} ,\\
w_0 = 2\sqrt{\pi} / 3  ,\quad& w_{1} = \sqrt{\pi} / 6 ,\quad& w_{-1} = \sqrt{\pi} / 6.
\end{aligned}
$$

如果用变换关系代换回 $\xi$，则 $\xi_0 = 0$, $\xi_{-1} = -\sqrt{3RT}$ 和 $\xi_{1} = \sqrt{3RT}$.

所以，积分 $I_m$ 可以被近似成：
$$
\begin{aligned}
I_m &= \int_{-\infty}^{+\infty} \exp \left( - \zeta^2 \right) \zeta^m \mathrm{d}\boldsymbol{\zeta} 
\\&\approx \sum_{i \in [-1, 0, 1]} w_{i} \zeta_{i}^m 
\\&= (2RT)^{-m/2} \sum_{i \in [-1, 0, 1]} w_{i} \xi_{i}^m
\end{aligned}
$$

也就是说，
$$
(2 R T)^{m/2} I_m \approx \sum_{i \in [-1, 0, 1]} w_{i} \xi_{i}^m
$$

综上所述，式(1-6)就会被改写为：
$$
\begin{aligned}
    I=& \int_{-\infty}^{+\infty} \psi_m(\xi) f^{eq} (\xi) \mathrm{d}\xi
    \\ \approx& \sum_{i \in [-1, 0, 1]} \xi_{i}^{m} \cdot \left\{
    \rho \frac{w_{i}}{\sqrt{\pi}} \cdot 
    \left[
        1 + \frac{\xi_{i} u}{R T} + \frac{(\xi_{i} u)^2}{2 (R T)^2} - \frac{u^2}{2 R T}
    \right]
    \right\}.
\end{aligned}
$$

到这一步，我们就可以定义 $\xi_{i}$ 方向的离散平衡态 $f^{eq}_{i}$ 为：
$$
f^{eq}_{i} = \rho W_{i} \cdot 
\left[
    1 + \frac{\xi_{i}  u}{c_s^2} + \frac{(\xi_{i} u)^2}{2 c_s^4} - \frac{u^2}{2 c_s^2}
\right].
\tag{1-7}
$$

其中的权重 $W_{i} = w_i / \sqrt{\pi}$。并且，如果令 $c = \sqrt{3RT}$， 则 $c_s = \sqrt{RT} = c / \sqrt{3}$ 。

至此，我们就整理出了 D1Q3 模型：

| ${i}$ | 离散速度 $\xi_{i}$ | 权重系数 $W_i$ |
|:-----:|:----------------:|:-------------:|
|  0 |  0 | 2/3 |
| -c | -c | 1/6 |
|  c |  c | 1/6 |

注意：在多数LBM实践中，我们通常将 $c$ 无量纲化，即令 $c=1$。

## 更高维的扩展

由于各个坐标轴之间是正交的，这种做法同样可以拓展到高维空间。

打个比方，对于三维问题， $\boldsymbol{\xi} = [\xi_x, \xi_y, \xi_z]^{\mathrm{T}}$， 则 $\mathrm{d}\boldsymbol{\xi} = \mathrm{d}\xi_x \mathrm{d}\xi_y \mathrm{d}\xi_z$。积分时就只需要将不同方向的分量分离出来就可以了。

对于二维问题，记 $\psi_{m,n}(\boldsymbol{\xi})=\xi_x^{m} \xi_y^{n}$，则

$$
\begin{aligned}
    I =& \int\psi_{m , n} (\boldsymbol{\xi}) f^{( \mathrm{s t} )} \mathrm{d}\boldsymbol{\xi} 
    \\=& \frac{\rho}{\pi} ( \sqrt{2 R T} )^{m+n} \Biggl\{\left( 1-\frac{u^{2}} {2 R T} \right) I_{m} I_{n} 
    \\&+ \frac{1} {\sqrt{2 R T}}  2 ( u_{x} I_{m+1} I_{n}+u_{y} I_{m} I_{n+1} )
    \\&+ \frac{1} {R T} (u_{x}^{2} I_{m+2} I_{n}+2 u_{x} u_{y} I_{m+1} I_{n+1}+u_{y}^{2} I_{m} I_{n+2}) \Biggr\}
\end{aligned}
\tag{1-8}
$$

式(1-8)同样可以整理成：
$$
I = \sum_{i,j\in[-1,0,1]} \psi(\boldsymbol{\xi}_{i,j}) \cdot \left\{ \rho \frac{w_i w_j}{\pi} 
\left[
    1 + \frac{\boldsymbol{\xi}_{i,j} \cdot \boldsymbol{u}}{c_s^2} + \frac{(\boldsymbol{\xi}_{i,j} \cdot \boldsymbol{u})^2}{2 c_s^4} - \frac{\|\boldsymbol{u}\|^2}{2 c_s^2}
\right]
\right\}
\tag{1-9}
$$

其中 $\boldsymbol{\xi}_{i,j}=\sqrt{2RT} (\zeta_i,\zeta_j)$。令 $W_{i,j} = w_i w_j / \pi$ ，便得出了 D2Q9 模型的权重。

类似的，对于三维问题，则可以整理出
$$
I = \sum_{i,j,k\in[-1,0,1]} \psi(\boldsymbol{\xi}_{i,j,k}) \cdot \left\{ \rho \frac{w_i w_j w_k}{\pi^{3/2}} 
\left[
    1 + \frac{\boldsymbol{\xi}_{i,j,k} \cdot \boldsymbol{u}}{c_s^2} + \frac{(\boldsymbol{\xi}_{i,j,k} \cdot \boldsymbol{u})^2}{2 c_s^4} - \frac{\|\boldsymbol{u}\|^2}{2 c_s^2}
\right]
\right\}
\tag{1-10}
$$

其中 $\boldsymbol{\xi}_{i,j,k}=\sqrt{2RT} (\zeta_i,\zeta_j,\zeta_k)$。令 $W_{i,j,k} = w_i w_j w_k / \pi^{3/2}$ ，便得出了 D3Q27 模型的权重。



# 待定系数法与DnQb模型

上面的内容可以说从数学上给出了明确的证明。然而实际使用中，其实还是用待定系数法来敲定 DnQb 模型更多一些。下面的内容里，为了保持和何雅玲老师的书标号一致，这里将第 $\alpha$ 个格子速度向量记作 $\boldsymbol{e}_{\alpha}$，其权重也记作 $\omega_{\alpha}$。

一般来说，d 维空间下的格子速度的 n 阶张量为：

$$
E_{i_1 i_2 ... i_n}^{(n)} = \sum_{\alpha} (\boldsymbol{e}_{\alpha})_{i_1} (\boldsymbol{e}_{\alpha})_{i_2} ... (\boldsymbol{e}_{\alpha})_{i_n}
\tag{2-1}
$$

假如你选取的离散速度集（共有 $M$ 个速度）是各向同性的，则应该满足

$$
E^{(2n+1)}_{i_1 i_2 ... i_{2n+1}} = 0,
$$

和

$$
E^{(2n)}_{i_1 i_2 ... i_{2n}} = \frac{M}{d (d+2) ... (d+2n-2)} \Delta^{(2n)}_{i_1 i_2 ... i_{2n}}
$$

下面我们首先从张量算子 $\Delta^{(2n)}$ 开始介绍。详情可见文献 [2] 的第3章，这篇文章在scihub上可以直接下载。

## 格子张量

张量算子 $\Delta^{(2n)}$ 实际上计算的是将 $2n$ 个指标全部两两配对的所有方式之和。它的具体形式为：

$$
\begin{aligned}
    \Delta^{(2)}_{ij}&=\delta_{ij} \\
    \Delta^{(4)}_{ijkl}&=\delta_{ij}\delta_{kl} + \delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk} \\
    \Delta^{(6)}_{ijklpq}&=\delta_{ij}\Delta^{(4)}_{ijkl} + \delta_{ik}\Delta^{(4)}_{jlpq} + \delta_{il}\Delta^{(4)}_{jkpq} + \delta_{ip}\Delta^{(4)}_{jklq} + \delta_{iq}\Delta^{(4)}_{jklp} 
\end{aligned}
$$  

并存在如下递归关系：

$$
\Delta^{(2n)}_{i_1 i_2 \dots i_{2n}} = \sum_{j=2}^{2n} \delta_{i_1, i_{j}} \cdot \Delta^{(2n-2)}_{i_2 \dots i_{j-1} i_{j+1} \dots i_{2n}}
$$

其中 $\delta_{ij}$ 为克罗内克积。

1. $2n$ 个指标分成 $n$ 对再求和，共有 $(2n−1)!!$ 项。
**计算思路为** ：第一个元素有 $(2n−1)$ 个选择配对；剩下 $(2n−2)$ 个元素之中的第一个有 $(2n−3)$ 个选择；后续依此类推。
2. 有的书会从排列组合视角出发，写成 $\frac{(2n)!}{2^n n!} = \frac{\text{总积}}{\text{偶数部分的积}}$，结果是等价的。
**计算思路为** ：$2n$ 个元素的全排列，除以每对内部的顺序($2n$)，再除以 $n$ 对之间的顺序($n!$)。


> 为了简洁表示 $\Delta^{(2n)}$ 的分量，文献 [2] 引入了 upper simplicial components，其索引是非递增的多重指标（如$(i_1,i_2,\dots,i_k)$，满足$i_1 \geq i_2 \geq \dots \geq i_k$）。 假设我们讨论一下二维空间的 $\Delta^{(4)}_{ijkl}$ ，其中的 $i,j,k,l$ 取 1 或 2。 那么下标的所有排列便是：
> 
> - 1111（全1，非递增）→ 值3
> - 2111（1个2、3个1，非递增）→ 值0
> - 2211（2个2、2个1，非递增）→ 值1
> - 2221（3个2、1个1，非递增）→ 值0
> - 2222（全2，非递增）→ 值3
> 
> 即 $\Delta^{(4)}=[3,0,1,0,3]$.

在后续的描述中，为了避免写一大堆下标，我们采用如下简记

$$
E^{(n)} = E^{(n)}_{i_1 i_2 \dots i_{n}}, \quad
\Delta^{(2n)} = \Delta^{(2n)}_{i_1 i_2 \dots i_{2n}}, \quad
\delta^{(2n)} = \delta^{(2n)}_{i_1 i_2 \dots i_{2n}},
$$

这里的 $\Delta^{(2n)}$ 就是 $\{i_1, i_2, \dots i_{2n}\}$ 一系列下标组成的 $\Delta^{(2n)}_{i_1 i_2 \dots i_{2n}}$。而 $\delta^{(2n)}$ 则是2n个下标的克罗内克积，当所有下标相等时才等于1，否则为0。

## 速度矩与权重系数

从 $E^{(2n)}$ 的定义中，若假设所有 $\boldsymbol{e}_{\alpha}$ 一样长，可得：

$$
\frac{1}{M} \sum_{\alpha} (\boldsymbol{e}_{\alpha} \cdot \boldsymbol{v})^{2n} = 
\underbrace{\frac{(2n-1)!!}{d (d+2) ... (d+2n-2)}}_{Q_{2n}} \|\boldsymbol{v}\|^{2n}
$$

对于二维情况 ($d=2$)，则 $Q_2 = 1/2$, $Q_4 = 3/8$, $Q_6 = 5/16$, $Q_8 = 35/128$。对于三维情况 ($d=3$)，则 $Q_{2n} = 1 / (2n+1)$.

同理，拓展到奇数阶，则有：

$$
\frac{1}{M} \sum_{\alpha} (\boldsymbol{e}_{\alpha} \cdot \boldsymbol{v})^{2n} (\boldsymbol{e}_{\alpha} \cdot \boldsymbol{v}) = 
Q_{2n} \|\boldsymbol{v}\|^{2n} \boldsymbol{v}
$$

因此，对于一个所有 $\boldsymbol{e}_{\alpha}$ 一样长的离散速度集来说，它的离散速度矩的权重是确定的一个数 —— $Q_{2n}$。

假如 ${\boldsymbol{e}_{\alpha}}$ 中的所有速度并不等长，我们则可以认为：

$$
E_{i_1 i_2 ... i_n}^{(n)} = \sum_{\alpha} 
\omega(\|\boldsymbol{e}_{\alpha}\|) \cdot (\boldsymbol{e}_{\alpha})_{i_1} (\boldsymbol{e}_{\alpha})_{i_2} ... (\boldsymbol{e}_{\alpha})_{i_n}
$$

其中 $\omega(\|\boldsymbol{e}_{\alpha}\|)$ 指代由特征速度的模 $\|\boldsymbol{e}_{\alpha}\|$ 决定的权重系数。

## DnQb 模型的确定

Qian等的 DnQb 模型采用如下的平衡态分布：

$$
f^{eq}_{\alpha} (\rho, \boldsymbol{u}; \boldsymbol{e}_{\alpha}, \omega_{\alpha}) = \omega_{\alpha} \rho \cdot 
\left[
    1 + \frac{\boldsymbol{e}_{\alpha} \cdot \boldsymbol{u}}{c_s^2} + \frac{(\boldsymbol{e}_{\alpha} \cdot \boldsymbol{u})^2}{2 c_s^4} - \frac{\|\boldsymbol{u}\|^2}{2 c_s^2}
\right].
$$

其中 $\rho$ 为流场密度；$\boldsymbol{u}$ 为流场速度，$c_s^2 = RT$。若是为了模拟Navier-Stokes方程，平衡态分布函数需要满足下面的约束：

$$
\begin{aligned}
    \sum_{\alpha} f^{eq}_{\alpha} &= \rho ,\\
    \sum_{\alpha} f^{eq}_{\alpha} \boldsymbol{e}_{\alpha} &= \rho \boldsymbol{u} ,\\
    \sum_{\alpha} f^{eq}_{\alpha} e_{\alpha i} e_{\alpha j} &= \rho u_i u_j + p \delta_{ij}.
\end{aligned}
\tag{2-2}
$$

在待定系数法的框架下，就是把已经设计好的离散速度集 $\{\boldsymbol{e}_{\alpha}\}$ 代入上面的约束中，然后求解对应的权重 $\omega_{\alpha}$。

在这种视角下，我们要怎么确定就是要用 D2Q9、D3Q15、D3Q9 这些形式呢？

首先，正四边形最多能保证三阶格子张量是各向同性的，而到了四阶就变成了 $E_{ijkl}^{(4)}|_{M=4} = 2 \delta_{ijkl}$。为了解决这个问题，我们可以在在原有最近邻链接基础上，增加更远邻的链接，扩大速度方向集合 $\{\boldsymbol{e}_{\alpha}\}$。 而正方晶格加对角线的方式，就构成了典型的 D2Q9 模型：此时速度方向等效于正八边形顶点，当 $M=8>4$ 时，$E^{(4)}$ 满足各向同性条件（$M$ 不整除4）。 

同样的道理，在三维空间中若是使用立方晶格，此时速度方向为立方体顶点（$M=8$，最近邻链接），对称性较低。$E^{(4)}$ 也是各向异性的，无法导出标准流体方程。此时也需要把一些对角线方向加进来，让 $E^{(4)}$ 满足各向同性条件。

### Example: 以 D2Q9 为例

首先，我们把 D2Q9 分成三组

| 组别 | 速度配置 | 权重 |
|-|-|-|
| 1 | $(0,0)$ | $\omega_0$ |
| 2 | cyc:$(\pm c,0)$ | $\omega_1=\omega_2=\omega_3=\omega_4$ |
| 3 | $(\pm c,\pm c)$ | $\omega_5=\omega_6=\omega_7=\omega_8$ |

对于组别 2，根据格子张量的表格，我们有
$$
\begin{aligned}
    \sum_{\alpha\in[1-4]} e_{\alpha i} &= 0 ,& \sum_{\alpha\in[1-4]} e_{\alpha i} e_{\alpha j} &= 2\delta_{ij} \cdot c^2\\
    \sum_{\alpha\in[1-4]} e_{\alpha i} e_{\alpha j} e_{\alpha k} &= 0 ,& \sum_{\alpha\in[1-4]} e_{\alpha i} e_{\alpha j} e_{\alpha k} e_{\alpha l} &= 2\delta_{ijkl} \cdot c^4\\
\end{aligned}
$$

对于组别 3，同理可得
$$
\begin{aligned}
    \sum_{\alpha\in[5-8]} e_{\alpha i} &= 0 ,& \sum_{\alpha\in[5-8]} e_{\alpha i} e_{\alpha j} &= 2\delta_{ij} \cdot (\sqrt{2} c)^2\\
    \sum_{\alpha\in[5-8]} e_{\alpha i} e_{\alpha j} e_{\alpha k} &= 0 ,& \sum_{\alpha\in[5-8]} e_{\alpha i} e_{\alpha j} e_{\alpha k} e_{\alpha l} &= (4 \Delta^{(4)}_{ijkl} - 8 \delta^{(4)}_{ijkl}) \cdot (\sqrt{2} c)^4\\
\end{aligned}
$$

代入式(2-2)，依次得：

$$
\rho (\omega_0 + 4 \omega_1 + 4 \omega_5) + 
\rho \|\boldsymbol{u}\|^2 \left[ \omega_1 \left(\frac{c^2}{c_s^4} - \frac{2}{c_s^2}\right) \right.\\
\left. + \omega_5 \left(\frac{2 c^2}{c_s^4} - \frac{2}{c_s^2}\right)  - \frac{\omega_0}{2 c_s^2} \right] = \rho
\tag{2-3}
$$

$$
\rho \boldsymbol{u} \left(
    \omega_1 \frac{2 c^2}{c_s^2} + \omega_5 \frac{4 c^2}{c_s^2}
\right) = \rho \boldsymbol{u}
\tag{2-4}
$$

$$
\begin{aligned}
\rho \omega_1 \sum_{\alpha=1}^{4} \left[
    e_{\alpha i} e_{\alpha j} \left(1 - \frac{u^2}{2 c_s^2}\right) + 
    \frac{e_{\alpha i} e_{\alpha j} e_{\alpha k} e_{\alpha l} u_k u_l}{2 c_s^4}
\right] &+ \\
\rho \omega_5 \sum_{\alpha=5}^{8} \left[
    e_{\alpha i} e_{\alpha j} \left(1 - \frac{u^2}{2 c_s^2}\right) + 
    \frac{e_{\alpha i} e_{\alpha j} e_{\alpha k} e_{\alpha l} u_k u_l}{2 c_s^4}
\right] &= \rho u_i u_j + p \delta_{ij}
\end{aligned}
\tag{2-5}
$$

由式(2-3)，得：
$$
\omega_0 + 4 \omega_1 + 4 \omega_5 = 1 \\
\omega_1 \left(\frac{c^2}{c_s^4} - \frac{2}{c_s^2}\right) + \omega_2 \left(\frac{2 c^2}{c_s^4} - \frac{2}{c_s^2}\right)  - \frac{\omega_0}{2 c_s^2} = 0
$$

由式(2-4)，得：
$$
(\omega_1 + 2 \omega_5) \frac{2 c^2}{c_s^2} =1
$$

> 式(2-5)需要用到如下化简：
> $$
> \begin{aligned}
> \Delta^{(4)}_{ijkl} u_{k} u_{l} 
> &= \delta_{ij} (\delta_{kl} u_{k} u_{l}) + (\delta_{ik}\delta_{jl}u_{k} u_{l}) + (\delta_{il}\delta_{jk} u_{k} u_{l}) \\
> &=\delta_{ij} \|\boldsymbol{u}\|^2 + 2 u_i u_j
> \end{aligned}
> $$
> 其中：
> - 第一项仅当 $k=l$ 时不为0。因此 $\delta_{kl} u_k u_l = \sum_{k=l} u_k u_l = \sum_{k} u_k u_k = \|\boldsymbol{u}\|^2$。
> - 第二项 $\delta_{ik}$ 要求 $k=i$ 时才非零，$\delta_{jl}$ 要求 $l=j$ 时才非零。因此对 $k,l$ 求和时，只有 $k=i$ 且 $l=j$ 的项保留，即 $u_i u_j$。
> - 第三项 $\delta_{il}$ 要求 $l=i$，$\delta_{jk}$ 要求 $k=j$。因此只有 $l=i$ 且 $k=j$ 的项保留，即 $u_i u_j$。

由于式(2-5)的左侧可被简化为

$$
\begin{aligned}
(\omega_1 - 4 \omega_5) \frac{c^4}{c_s^4} \rho \delta_{ijkl} u_k u_l &+
\rho c^2 \delta_{ij} \left[ (1 - \frac{\|\boldsymbol{u}\|^2}{2 c_s^2}) (2\omega_1 + 4\omega_5) \right. \\
&\left. + 2 \omega_5 \frac{\|\boldsymbol{u}\|^2 c^2}{c_s^4} \right] + 
\rho u_i u_j w_5 \frac{4 c^2}{c_s^4}
\tag{2-6}
\end{aligned}
$$

所以
$$
\begin{aligned}
    w_5 \frac{4 c^2}{c_s^4} &= 1 \\
    \omega_1 - 4 \omega_5 &= 0 \\
     \left(1 - \frac{\|\boldsymbol{u}\|^2}{2 c_s^2}\right) (2\omega_1 + 4\omega_5) + 2 \omega_5 \frac{\|\boldsymbol{u}\|^2 c^2}{c_s^4} &= \frac{p}{\rho c^2}
\end{aligned}
$$

联立这些方程就可以解得： $\omega_0 = 4/9$, $\omega_1 = 1/9$, $\omega_5 = 1/36$, $c_s = c / \sqrt{3}$。

## 附录

（1）对于二维正M边形网格，$\boldsymbol{e}_{\alpha} = e \left( \cos(\frac{2\pi\alpha}{M}), \sin(\frac{2\pi\alpha}{M})\right)$。此时阶数 $n$ 与 $M$ 的对应关系为：

|各向同性张量 $E^{(n)}$| 正多边形边数 $M$ |
|--|--|
| $E^{(2)}$ | $M \gt 2$ |
| $E^{(3)}$ | $M \ge 2, M \neq 3$ |
| $E^{(4)}$ | $M \gt 2, M \neq 4$ |
| $E^{(5)}$ | $M \ge 2, M \neq [3,5]$ |
| $E^{(6)}$ | $M \gt 4, M \neq 6$ |
| $E^{(7)}$ | $M \ge 2, M \neq [3,5,7]$ |

（2）文献[2]给出了三维情况下 $E^{(n)}$ 的各向同性情况，下表用 Y 和 N 分别表示是和否。其中 cyc 表示循环排列, $\phi=(1+\sqrt{5})/2$。

| 几何形状 | 速度配置 | M | $E^{(2)}$ | $E^{(3)}$ | $E^{(4)}$ | $E^{(5)}$ | $E^{(6)}$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 四面体   | $(1,1,1)$, cyc:$(1,-1,-1)$ | 4  | Y | N | N | N | N |
| 立方体   | cyc:$(\pm 1,\pm 1,\pm 1)$ | 8 | Y | Y | N | Y | N |
| 八面体   | cyc:$(\pm 1,0,0)$ | 6 | Y | Y | N | Y | N |
| 十二面体 | $(\pm 1,\pm 1,\pm 1)$, cyc: $(0, \pm \phi-1, \pm \phi)$ | 20 | Y | Y | Y | Y | N |
| 二十面体 | cyc: $(0, \pm \phi, \pm 1)$ | 12 | Y | Y | Y | Y | N |

（3）两种基本正方形格子模型的偶数阶张量

| 速度配置 | M | $E^{(2)}$ | $E^{(4)}$ | $E^{(6)}$ |
| --- | --- | --- | --- | --- |
| cyc:$(\pm 1,0)$ | 4 | $2 \delta^{(2)}$ | $2 \delta^{(4)}$ | $2 \delta^{(6)}$ |
| $(\pm 1,\pm 1)$ | 8 | $4 \delta^{(2)}$ | $4 \Delta^{(4)} - 8 \delta^{(4)}$ | $\frac{4}{3} \Delta^{(6)} - 16 \delta^{(6)}$ |

（4）几种最对称的三维Bravais晶格的偶数阶张量

|形状| 速度配置 | M | $E^{(2)}$ | $E^{(4)}$ | $E^{(6)}$ |
|---| --- | --- | --- | --- | --- |
|Primitive cubic| cyc:$(\pm 1,0,0)$ | 6 | $2 \delta^{(2)}$ | $2 \delta^{(4)}$ | $2 \delta^{(6)}$ |
|Body-centered cubic| $(\pm 1,\pm 1,\pm 1)$ | 8  | $8 \delta^{(2)}$ | $8 (\Delta^{(4)} - 2\delta^{(4)})$ | $8 (\Delta^{(6)} - 2 \Delta^{(4,2)} + 16\delta^{(6)})$ |
|Face-centered cubic| cyc:$(\pm 1,\pm 1,0)$ | 12 | $8 \delta^{(2)}$ | $4 (\Delta^{(4)} - \delta^{(4)})$ | $4 (\Delta^{(4,2)} + 13\delta^{(6)})$ |

在文献[2]的原话中： 
> $\Delta^{(n,m)}$ is the sum of all possible products of pairs of Kronecker delta symbols with $n$ and $m$ indices。

具体到 $\Delta^{(4,2)} = \Delta^{(4,2)}_{ijklmn}$ ，它应该是 $\mathrm{C}_{6}^{2} = 15$ 个项的和。打个比方：第1个项是 $\delta_{ij}\Delta^{(4)}_{klmn}$；第2个项是 $\delta_{ik}\Delta^{(4)}_{jlmn}$；后续依次类推，直到把15个组合都轮换出来，然后求和。

# 参考

[1] HE X, LUO L S. **Theory of the lattice Boltzmann method: From the Boltzmann equation to the lattice Boltzmann equation**[J]. _Physical Review E_, 1997, 56(6): 6811-6817. [DOI:10.1103/PhysRevE.56.6811](https://link.aps.org/doi/10.1103/PhysRevE.56.6811).  
[2] WOLFRAM S. **Cellular automaton fluids 1: Basic theory**[J]. _Journal of Statistical Physics_, 1986, 45(3-4): 471-526. [DOI:10.1007/BF01021083](http://link.springer.com/10.1007/BF01021083).   
[3] 何雅玲, 王勇, 李庆, 童自翔. **格子Boltzmann方法的理论及应用**[M]. 1 版. 北京: 高等教育出版社, 2023.

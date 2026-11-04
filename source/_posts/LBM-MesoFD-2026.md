---
title: 【笔记】MesoFD - 基于平衡态的有限差分LBM
date: 2026-11-04 21:54:27
tags: [格子Boltzmann方法]
categories:
- [格子Boltzmann方法]
katex: true
---
<script defer src="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.17.0/contrib/copy-tex.min.js"></script>


最近看到了 Shi et al. [1] 的新文章很有意思: MesoFD 格式 —— 从平衡态就能去构建各种偏微分方程的有限差分LBM格式 (并包括其对应的源项和平衡态). 

我这里主要是把自己感兴趣的 “Navier-Stokes 方程的构造” 部分内容摘录了出来, 并做分享 (没想到用Taylor展开来分析真的不难啊). 

至于文章里的稳定性分析、DnQb构造等等, 还请各位阅读施保昌老师的原文.

> [1] Shi B, Du R, Chai Z. **Equilibrium-distribution-function-based mesoscopic finite-difference methods for partial differential equations: modeling and analysis**. Physical Review E, 2026, 113(4)：045305. DOI:10.1103/61zg-dmhh. [Also available at [arXiv:2505.17057](https://arxiv.org/abs/2505.17057)]  
> [2] Xinmeng Chen, Zhenhua Chai, Jinlong Shang, Baochang Shi. **Multiple-relaxation-time finite-difference lattice Boltzmann model for the nonlinear convection-diffusion equation**. Physical Review E, 2021, 104:035308. DOI:10.1103/physreve.104.035308


# 基本构型

Shi et al. [1] 提出了一个 MesoFD 格式, 核心在于将 $\boldsymbol{c}_j$ 方向的离散分布函数 $f_j$ 直接表示为邻近时空点处平衡分布函数 $f_j^{eq}$ 的线性组合 ($j=0,1, \dots b-1$, 共b个离散速度矢量). 

$$
f_j (\boldsymbol{x},t)=\sum_{q} \sum_{k} \left( a_{kq} f_j^{eq}(\boldsymbol{x}+k\boldsymbol{c}_j \delta_t, t+q\delta_t) + b_{kq} \delta_t F_j (\boldsymbol{x}+k\boldsymbol{c}_j \delta_t, t+q\delta_t) \right)
\tag{1}
$$

其中 $a_{kq}, b_{kq}$ 均为参数, $F_j$ 为该方向的源项. $-m_q \le k \le m_q$ 以及 $-n \le q \le 0$ 均为整数.  $n$ 和 $m_q$ 均不小于1, 表明这一组合至少使用单层空间数据和 $t-\delta_t$ 的时间数据. 

为便于书写, Shi et al. 用了下面的简记符号

$$
f_j^{eq}(\boldsymbol{x}+k\boldsymbol{c}_j \delta_t, t+q\delta_t) = f_j^{eq,[k,q]}
,\quad
F_j (\boldsymbol{x}+k\boldsymbol{c}_j \delta_t, t+q\delta_t) = F_j^{[k,q]}
$$

所以式(1)写为： $f_j (\boldsymbol{x},t)=\sum_{k,q} [ a_{kq} f_j^{eq,[k,q]} +  b_{kq} \delta_t F_j^{[k,q]}]$ . 对于经典的Navier-Stokes Equations (NSEs), 依旧有：

$$
\rho(\boldsymbol{x},t) = \sum_{j} \sum_{k,q} \left[ a_{kq} f_j^{eq,[k,q]} +  b_{kq} \delta_t F_j^{[k,q]} \right]
\tag{2a}
$$

以及

$$
(\rho\boldsymbol{u})(\boldsymbol{x},t) = \sum_{j}\sum_{k,q}  \boldsymbol{c}_j \left[ \bar{a}_{kq} f_j^{eq,[k,q]} +  \bar{b}_{kq} \delta_t F_j^{[k,q]} \right]
\tag{2b}
$$

需要强调, 计算动量 $\rho\boldsymbol{u}$ 所用的系数 $\{\bar{a}_{kq}, \bar{b}_{kq}\}$ 和式(1)里所用的那套不一样. 

## 式(1)和中心差分方法的联系

若对一个空间内的函数 $g(\boldsymbol{x})$ 沿着 $\boldsymbol{c}_j$ 方向做差分 (记 $\boldsymbol{c}_{\bar{j}} = - \boldsymbol{c}_j$), 则记：
   * 一阶空间中心差分算符： $\Delta_{j,k}^{(1)} g(\boldsymbol{x}) = g(\boldsymbol{x} + k \boldsymbol{c}_j \delta_t) - g(\boldsymbol{x} - k \boldsymbol{c}_j \delta_t)$； 且 $\Delta_{j}^{(1)} g(\boldsymbol{x}) = \Delta_{j,1}^{(1)} g(\boldsymbol{x})$. 
   * 二阶空间中心差分算符： $\Delta_{j,k}^{(2)} g(\boldsymbol{x}) = g(\boldsymbol{x} + k \boldsymbol{c}_j \delta_t) - g(\boldsymbol{x}) +  g(\boldsymbol{x} - k \boldsymbol{c}_j \delta_t)$； 且 $\Delta_{j}^{(2)} g(\boldsymbol{x}) = \Delta_{j,1}^{(2)} g(\boldsymbol{x})$. 

我们便可发现：  $g(\boldsymbol{x} \pm k \boldsymbol{c}_j \delta_t) = (\Delta_{j,k}^{(2)} \pm \Delta_{j,k}^{(1)} + g(\boldsymbol{x})) / 2$ 

那么式(1)要进一步改写为：

$$
\begin{aligned}
    f_j (\boldsymbol{x}, t) =& \sum_{k,q} a_{kq} [1 + (\Delta_{j,k}^{(2)} + \Delta_{j,k}^{(1)})/2] f_j^{eq}(\boldsymbol{x}, t+q\delta_t) + \\&
    \delta_t \sum_{k,q} b_{kq} [1 + (\Delta_{j,k}^{(2)} + \Delta_{j,k}^{(1)})/2] F_j (\boldsymbol{x}, t+q\delta_t) 
\end{aligned}
$$

# 算法框架

这个算法其实也比较容易. 假定已有 $t$ 时刻的数据, 需要前推到 $t+\delta_t$ 时刻. 下面也就只简单介绍这个单步过程.

(1) 对网格点 $\boldsymbol{x}$, 通过它的宏观物理量去计算平衡态 $f_j^{eq}(\boldsymbol{x}, t)$ 和 外力项 $F_j(\boldsymbol{x}, t)$. [具体公式详见原文的第III.C节, 这里不做介绍]

(2) 用式(1)计算所有点的 $f_j$.

(3) 用式(2)更新宏观物理量.

(4) 施加边界条件.



# MesoFD 的直接 Taylor 展开

## 单个分量的分析

下面讨论对式(1)直接 Taylor 展开的分析, 并最高展开至 $N \geq 1$ 阶后做截断. 

记 Taylor 展开算子

$$
D_{k,q} = q \partial_{t} + k \nabla \cdot \boldsymbol{c}_{j}
$$

那么 $f_j^{eq,[k,q]}$ 和 $F_j^{[k,q]}$ 分别展开成

$$
f_j^{eq}(\boldsymbol{x}+k\boldsymbol{c}_j \delta_t, t+q\delta_t) = \sum_{l=0}^{N} \frac{(\delta_{t})^l}{l!} D_{k,q}^l f_j^{eq}(\boldsymbol{x}, t) + O((\delta_{t})^{N+1})
$$

和

$$
F_j(\boldsymbol{x}+k\boldsymbol{c}_j \delta_t, t+q\delta_t) = \sum_{l=0}^{N} \frac{(\delta_{t})^l}{l!} D_{k,q}^l F_j(\boldsymbol{x}, t) + O((\delta_{t})^{N+1})
$$

考虑到算子 $D_{k,q}$ 高阶导表示为

$$
\begin{aligned}
    D_{k,q}^{l} &= (q \partial_{t} + k \nabla \cdot \boldsymbol{c}_{j})^l = \sum_{m=0}^{l} C_l^m \, (q \partial_{t})^{l-m} \, (k\nabla \cdot \boldsymbol{c}_j)^m ,\\
    C_l^m &= l! / (m! (l-m)!)
\end{aligned}
$$

那么, 我们把和空间、时间有关的系数部分提取出来

$$
A_{lm} = \frac{C_l^m}{l!} \sum_{k,q} a_{kq} q^{l-m} k^{m} ,\quad
B_{lm} = \frac{C_l^m}{l!} \sum_{k,q} b_{kq} q^{l-m} k^{m}
$$

则式(1)的展开要写作：

$$
f_j (\boldsymbol{x},t) = \sum_{l=0}^{N} (\delta_t)^l \sum_{m=0}^{l} \left[
    A_{lm} \partial_{t}^{l-m} (\nabla \cdot \boldsymbol{c}_j)^m f_j^{eq} + 
    \delta_t B_{lm} \partial_{t}^{l-m} (\nabla \cdot \boldsymbol{c}_j)^m F_j
\right] + O((\delta_{t})^{N+1})
\tag{3}
$$

其中 $f_j^{eq} = f_j^{eq}(\boldsymbol{x},t)$, $F_j = F_j(\boldsymbol{x},t)$.

取 $N=2$ , 由于 $F_j$ 部分在 $l=2$ 时直接到了 $O(\delta_t^3)$ 量级, 我们得到

$$
\begin{aligned}
f_j (\boldsymbol{x},t) =& \underbrace{A_{00} f_j^{eq} + \delta_t B_{00} F_j}_{l=0}  
\\& + \underbrace{\delta_t (A_{10} \partial_t + A_{11} \nabla \cdot \boldsymbol{c}_j) f_j^{eq} + \delta_t^2 (B_{10} \partial_t + B_{11} \nabla \cdot \boldsymbol{c}_j) F_j}_{l=1}  
\\& +\underbrace{\delta_t^2 [A_{20} \partial_t^2 + A_{21} \partial_t \, (\nabla \cdot \boldsymbol{c}_j) + A_{22} (\nabla \cdot \boldsymbol{c}_j)^2] f_j^{eq} }_{l=2 \,(f^{eq} \text{ part})} 
\\& + O(\delta_t^3)
\end{aligned}
$$


## 分布函数的速度矩

这里我们写列出 k 阶速度矩的表达式：

$$
\begin{aligned}
    M_0 &= \sum_j f_j^{eq}  = \sum_j f_j ,&   M_{0F} &= \sum_j F_j, \\
    M_1 &= \sum_j \boldsymbol{c}_j f_j^{eq},& M_{1F} &= \sum_j \boldsymbol{c}_j F_j, \\
    M_k &= \sum_j \underbrace{\boldsymbol{c}_j \dots \boldsymbol{c}_j}_{k} f_j^{eq},& M_{kF} &= \sum_j \underbrace{\boldsymbol{c}_j \dots \boldsymbol{c}_j}_{k} F_j, (k>1)\\
\end{aligned}
$$

按理来说, 只要把这些速度矩也匹配上就能去还原宏观方程. 

毕竟我们知道, 对于常用的正方形+对称 DnQb 模型, 如果是还原Navier-Stokes方程, 就有: 
  * $M_0 = \rho$.
  * $M_1 = \rho \boldsymbol{u}$.
  * $M_2 = \rho \boldsymbol{u} \boldsymbol{u} + p \boldsymbol{I}$.

不过, 先别急. 

让我们注意一下式(3)中的 $(\nabla \cdot \boldsymbol{c}_j)^m f_j^{eq}$ 部分 (暂不考虑前面的时间偏导), 如果我们对这部分对所有 $j$ 进行求和 (零阶速度矩), 就会发现

$$
\begin{aligned}
    \sum_{j} (\nabla \cdot \boldsymbol{c}_j)^m f_j^{eq} &= (\nabla \dots \nabla) \otimes \sum_{j} (\underbrace{\boldsymbol{c}_j \dots \boldsymbol{c}_j}_{k} f_j^{eq}) \\
    &= \nabla^{m} \otimes M_{m}
\end{aligned}
$$

也就是说, $M_0$ 的内部就自然包含了更高的 m 阶矩 $M_{m}$ ($m=[0,1,..,l]$, $l=[0,1,...N]$, 且 $N \geq 1$). 这个情况对 $M_k$ 以及源项部分的矩也是一样的. 


将式(3)代入 $M_0$, 可得

$$
\begin{aligned}
    M_0 =& A_{00} M_0 + \sum_{l=1}^{N} \delta_{t}^l \sum_{m=0}^{l} [A_{lm} \, \partial_{t}^{l-m} \, \nabla^{m} \otimes M_{m}] \\
    & + \sum_{l=0}^{N-1} \delta_{t}^{l+1} \sum_{m=0}^{l} [B_{lm} \, \partial_{t}^{l-m} \, \nabla^{m} \otimes M_{mF}] + O(\delta_{t}^{l+1})
\end{aligned}
$$

令 $A_{00}=1$ 以消去左右两侧的 $M_0$, 整理得到

$$
\sum_{l=1}^{N} \delta_{t}^{l-1} \sum_{m=0}^{l} [A_{lm} \, \partial_{t}^{l-m} \, \nabla^{m} \otimes M_{m}] + \sum_{l=0}^{N-1} \delta_{t}^l \sum_{m=0}^{l} [B_{lm} \, \partial_{t}^{l-m} \, \nabla^{m} \otimes M_{mF}] = O(\delta_{t}^l)
\tag{4}
$$

同理, 将式(3)代入 $M_1$, 可得

$$
\begin{aligned}
    M_1 =& \bar{A}_{00} M_1 + \sum_{l=1}^{N} \delta_{t}^l \sum_{m=0}^{l} [\bar{A}_{lm} \, \partial_{t}^{l-m} \, \nabla^{m} \otimes M_{m}] \\
    & + \sum_{l=0}^{N-1} \delta_{t}^{l+1} \sum_{m=0}^{l} [\bar{B}_{lm} \, \partial_{t}^{l-m} \, \nabla^{m} \otimes M_{mF}] + O(\delta_{t}^{l+1})
\end{aligned}
$$

其中 $\bar{A}_{lm}, \bar{B}_{lm}$ 的定义仿照之前的写法:

$$
\bar{A}_{lm} = \frac{C_l^m}{l!} \sum_{k,q} \bar{a}_{kq} q^{l-m} k^{m} ,\quad
\bar{B}_{lm} = \frac{C_l^m}{l!} \sum_{k,q} \bar{b}_{kq} q^{l-m} k^{m}
$$

令 $\bar{A}_{00}=1$ 以消去左右两侧的 $M_1$, 整理得到

$$
\sum_{l=1}^{N} \delta_{t}^{l-1} \sum_{m=0}^{l} [\bar{A}_{lm} \, \partial_{t}^{l-m} \, \nabla^{m} \otimes M_{m+1}] + \sum_{l=0}^{N-1} \delta_{t}^l \sum_{m=0}^{l} [\bar{B}_{lm} \, \partial_{t}^{l-m} \, \nabla^{m} \otimes M_{(m+1)F}] = O(\delta_{t}^l)
\tag{5}
$$



【1】对于式(4), 当 $N=1$ 时, 则:

$$
A_{10} \partial_{t} M_0 + A_{11} \nabla \cdot M_1 = - B_{00} M_{0F} + O(\delta_t)
\tag{6a}
$$

当 $N=2$ 时, 则:

$$
\begin{aligned}
    ( A_{10} \partial_{t} M_0 + A_{11} \nabla \cdot M_1 ) &+ \delta_t [
      A_{20} \partial_{t}^2 M_0 + A_{21} \partial_{t} \nabla \cdot M_1 + A_{22} \nabla\nabla : M_2
    ] \\
    & = -B_{00} M_{0F} - \delta_t (B_{10} \partial_{t} M_{0F} + B_{11} \nabla \cdot M_{1F}) + O(\delta_t^2)
\end{aligned}
\tag{6b}
$$

【2】对于式(5), 当 $N=1$ 时, 则:

$$
\bar{A}_{10} \partial_{t} M_1 + \bar{A}_{11} \nabla \cdot M_2 = -\bar{B}_{00} M_{1F} + O(\delta_t)
\tag{7a}
$$

当 $N=2$ 时, 则:

$$
\begin{aligned}
    ( \bar{A}_{10} \partial_{t} M_1 + \bar{A}_{11} \nabla \cdot M_2 ) &+ \delta_t [
      \bar{A}_{20} \partial_{t}^2 M_1 + \bar{A}_{21} \partial_{t} \nabla \cdot M_2 + \bar{A}_{22} \nabla\nabla : M_3
    ] \\
    & = -\bar{B}_{00} M_{1F} - \delta_t (\bar{B}_{10} \partial_{t} M_{1F} + \bar{B}_{11} \nabla \cdot M_{2F}) + O(\delta_t^2)
\end{aligned}
\tag{7b}
$$

Shi et al. 做出下面代换来简化式(6-7), 这些代换在下文会用到: 
* $A_{10} = A_{11} = - B_{00}$,  
* $\bar{A}_{10} = \bar{A}_{11} = - \bar{B}_{00}$,  
* $\displaystyle S = -\frac{B_{00}}{A_{10}} M_{0F} = M_{0F}$,  
* $\displaystyle \boldsymbol{F} = -\frac{\bar{B}_{00}}{\bar{A}_{10}} M_{1F} = M_{1F}$.


# 与非线性对流扩散方程的匹配 (式(6))

这里我们先把非线性对流扩散方程写出来, 作为参照:

$$
\begin{aligned}
    \partial_{t} \phi + \nabla \cdot \bold{B} &= \boldsymbol{S} + \alpha \nabla\nabla:\bold{D} + O(\delta_t^2)\\
    &= \boldsymbol{S} + \nabla \cdot [\alpha \nabla \cdot \bold{D}] + O(\delta_t^2)
\end{aligned}
\tag{8}
$$

若我们只考虑 _"各向同性"+"不可压缩流"_ 来降低心智负担, 则
- $\phi$ 为输运的守恒量;
- $\boldsymbol{u}$ 为 $\phi$ 的输运速度;
- $\alpha$ 为扩散系数 (常数);
- $\bold{B}=\phi \, \boldsymbol{u}$ (这里对不可压缩流有 $\nabla\cdot\boldsymbol{u}=0$ , 所以 $\nabla\cdot\bold{B}=\boldsymbol{u}\cdot\nabla\phi$);
- $\bold{D}=\phi \, [\bold{I}]$.
- $\boldsymbol{S}$ 为源项向量.

式(6)可化简为

$$
\partial_{t} M_0 + \nabla \cdot M_1 = S + O(\delta_t)
\tag{9a}
$$

$$
\begin{aligned}
    ( \partial_{t} M_0 + \nabla \cdot M_1 ) &+ \frac{\delta_t}{A_{10}} [
      A_{20} \partial_{t}^2 M_0 + A_{21} \partial_{t} \nabla \cdot M_1 + A_{22} \nabla\nabla : M_2
    ] \\
    & = S - \frac{\delta_t}{A_{10}} (B_{10} \partial_{t} M_{0F} + B_{11} \nabla \cdot M_{1F}) + O(\delta_t^2)
\end{aligned}
\tag{9b}
$$

而 $\partial_{t}^2 M_0$ 可以由式(9a)的时间偏导表示, 这就和 Chapman-Enskog 展开那套东西比较相似: 高阶矩的计算需要借用低阶矩的方程来化简.

$$
\partial_{t}^2 M_0 = - \partial_{t} \nabla \cdot M_1 + \partial_{t} S + O(\delta_t)
$$

则式(9b)可以进一步化简为

$$
\begin{aligned}
    ( \partial_{t} M_0 + \nabla \cdot M_1 ) &+ \frac{\delta_t}{A_{10}} (B_{10} + A_{20}) \partial_{t}S \\
    & = S - \frac{\delta_t}{A_{10}} [(A_{21} - A_{20}) \partial_{t} \nabla \cdot M_1 \\
    &\quad + B_{11} \nabla \cdot M_{1F} + A_{22} \nabla\nabla : M_2] \\
    &\quad + O(\delta_t^2)
\end{aligned}
\tag{9c}
$$

令 $(B_{10} + A_{20}) \partial_{t}S = 0$ , 对比式(9c)和非线性对流扩散方程(式(8)), 我们现在就可以开始做参数匹配了:

$$
\begin{aligned}
    M_0 = \phi, \,
    M_1 &= \bold{B}, \,
    M_2 = \beta c_s^2 \bold{D} + \bold{C}, \,
    \alpha = -\frac{A_{22}}{A_{10}} \delta_t \beta c_s^2 ,\\
    M_{1F} &= -[(A_{21} - A_{20}) \partial_{t}\bold{B} + A_{22} \nabla\cdot\bold{C}]
\end{aligned}
$$

这里和 $\bold{C}$ 有关的部分用于消除式(8c)中的 $A_{22} \nabla\nabla : M_2$ 项.

# 与Navier-Stokes方程的匹配 (式(6) + 式(7))

首先, 式(6)的二阶结果即为式(9c). 这里只需要结合式(7)的二阶结果, 再和Navier-Stokes方程组进行参数匹配即可. 所以下面主要讨论式(7).

式(7)可化简为

$$
\partial_{t} M_1 +\nabla \cdot M_2 = \boldsymbol{F} + O(\delta_t)
\tag{10a}
$$


$$
\begin{aligned}
    ( \partial_{t} M_1 + \nabla \cdot M_2 ) &+ \frac{\delta_t}{\bar{A}_{10}} [
      \bar{A}_{20} \partial_{t}^2 M_1 + \bar{A}_{21} \partial_{t} \nabla \cdot M_2 + \bar{A}_{22} \nabla\nabla : M_3
    ] \\
    & = \boldsymbol{F} - \frac{\delta_t}{\bar{A}_{10}} (\bar{B}_{10} \partial_{t} \boldsymbol{F} + \bar{B}_{11} \nabla \cdot M_{2F}) + O(\delta_t^2)
\end{aligned}
\tag{10b}
$$

同理, 用式(10a)表示 $\partial_{t}^2 M_1$, 即:

$$
\partial_{t}^2 M_1 = -\partial_{t} \nabla \cdot M_2 + \partial_{t}\boldsymbol{F} + O(\delta_t)
$$

再次代入式(10b)进行化简, 得到:

$$
\begin{aligned}
    ( \partial_{t} M_1 + \nabla \cdot M_2 ) &+ \frac{\delta_t}{\bar{A}_{10}} (\bar{B}_{10} \partial_{t}\boldsymbol{F} + \bar{A}_{20} \partial_{t}\boldsymbol{F})  \\
    & = \boldsymbol{F} - \frac{\delta_t}{\bar{A}_{10}} [(\bar{A}_{21} - \bar{A}_{20}) \partial_{t} \nabla \cdot M_2 \\
    &\quad + \bar{B}_{11} \nabla \cdot M_{2F} + \bar{A}_{22} \nabla\nabla : M_3] \\
    &\quad + O(\delta_t^2)
\end{aligned}
\tag{10c}
$$

【1】对于式(6)的二阶结果(即式(9c)), 再次利用式(10a), 可得

$$
\partial_{t} \nabla \cdot M_1 = -\nabla \nabla : M_2 + \nabla \cdot \boldsymbol{F} + O(\delta_t)
$$

再令 $-B_{11} \boldsymbol{F} = (A_{21} - A_{20}) \boldsymbol{F}$ 和 $(B_{10} + A_{20}) \partial_{t}S = 0$, 全部代入到式(9c)中, 得到

$$
\partial_{t} M_0 + \nabla\cdot M_1 = S - \frac{\delta_t}{A_{10}} [(A_{22} - (A_{21} - A_{20}))\nabla \nabla : M_2] + O(\delta_t^2)
\tag{11a}
$$

式(11a)将用于匹配连续性方程.

【2】对于式(7)的二阶结果(即式(10c))

令 $(\bar{B}_{10} \partial_{t}\boldsymbol{F} + \bar{A}_{20} \partial_{t}\boldsymbol{F}) = 0$, 得到

$$
\begin{aligned}
    ( \partial_{t} M_1 + \nabla \cdot M_2 ) & = \boldsymbol{F} - \frac{\delta_t}{\bar{A}_{10}} [(\bar{A}_{21} - \bar{A}_{20}) \partial_{t} \nabla \cdot M_2 \\
    &\quad + \bar{B}_{11} \nabla \cdot M_{2F} + \bar{A}_{22} \nabla\nabla : M_3] \\
    &\quad + O(\delta_t^2)
\end{aligned}
\tag{11b}
$$

式(11a)将用于匹配Navier-Stokes方程.

【3】 针对式(11a)和式(11b), 我们需要给出如下参数匹配 


$$
\begin{aligned}
    M_0 &= \rho , \\
    M_1 &= \rho \bold{u} , \\
    M_2 &= \rho \bold{u} \bold{u} + p \bold{I} , \\
    M_3 &= c_s^2 \rho \bold{u} \bold{I} + M_{30} , \\
    M_{2F} &= \frac{-1}{\bar{B}_{11}} [(\bar{A}_{21} - \bar{A}_{20}) \partial_{t} M_2 + \bar{A}_{22} ( \nabla \cdot M_{30} +
    \\&\quad (\nabla \cdot \rho \bold{u} \bold{I} + \bold{u} \nabla\rho + (\bold{u} \nabla\rho)^{T}) c_s^2 \, )] \\
\end{aligned}
$$

且运动粘度 $\nu = \frac{\bar{A}_{22}}{\bar{A}_{10}} \delta_t c_s^2$. 则式(11a)和式(11b)被转化为Navier-Stokes方程组

$$
\begin{aligned}
    \partial_t \rho + \nabla \cdot (\rho \boldsymbol{u}) &= 
    S - \frac{\delta_t}{A_{10}} [(A_{22} - (A_{21} - A_{20}))\nabla \nabla : M_2] + O(\delta_t^2) \\
    \partial_t (\rho \boldsymbol{u}) + \nabla \cdot (\rho \boldsymbol{u}\boldsymbol{u} + p \bold{I}) &= \boldsymbol{F} + \nu \nabla \cdot \rho [\nabla\boldsymbol{u} + (\nabla\boldsymbol{u})^T] + O(\delta_t^2) \\
\end{aligned}
$$


> **注:**
> 1. Shi et al. 原文式(61a)写的 $M_3 = c_s^2 \Delta \cdot \rho \bold{u} + M_{30}$ 表达式有误 —— 拉普拉斯算子 $\Delta$ 没法点乘在一个矢量 $\rho \boldsymbol{u}$ 上.
> 2. 这里的 $\bold{u} \bold{I}$ 展开后应该是 $u_{\alpha} \delta_{\beta\gamma} + u_{\beta} \delta_{\alpha\gamma} + u_{\gamma} \delta_{\alpha\beta}$.
> 3. $p = c_s^2 \rho$.

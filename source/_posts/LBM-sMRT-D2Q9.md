---
title: 【笔记】Single MRT-featured LBM
date: 2026-08-07 09:38:10
tags:
  - 格子Boltzmann方法
categories:
  - 格子Boltzmann方法
---

<script defer src="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.17.0/katex.min.js" integrity="sha512-LQNxIMR5rXv7o+b1l8+N1EZMfhG7iFZ9HhnbJkTp4zjNr5Wvst75AqUeFDxeRUa7l5vEDyUiAip//r+EFLLCyA==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>
<script defer src="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.17.0/contrib/copy-tex.min.js" integrity="sha512-cQxSkTM4RvFAjdBeBDkrllhYfERwZWjM/oijKfPsmhR0JI2U3fXSlbUaJp5SlgBh/FzYmMyWudRTXLnNj3MbEA==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex/dist/katex.css">
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


最近在网上刷到了这篇 arxiv 预印本，

> Zhou J G. **Single MRT-featured Lattice Boltzmann Method (SmrtLBM)**[PP/OL]. arXiv（2024-09-02）. http://arxiv.org/abs/2409.01076.

加上休息时间有空，就简单看了这篇文章并做了笔记。

# 基本原理

这篇文章以 D2Q9 离散速度模型为例，对其多松弛形式进行了简化。

## 多松弛模型 (MRT) 的简单介绍

定义：
* $\boldsymbol{f} = [f_0(\boldsymbol{x}, t), ..., f_8(\boldsymbol{x}, t)]^\mathrm{T}$ 分布函数向量；
* $\boldsymbol{f}^{(\mathrm{pc})} = [f_0^{(\mathrm{pc})}(\boldsymbol{x}, t), ..., f_8^{(\mathrm{pc})}(\boldsymbol{x}, t)]^\mathrm{T}$ 碰撞后的分布函数向量；  
* $\boldsymbol{f}^{(\mathrm{eq})} = [f_0^{(\mathrm{eq})}(\boldsymbol{x}, t), ..., f_8^{(\mathrm{eq})}(\boldsymbol{x}, t)]^\mathrm{T}$ 平衡态分布函数向量。

> 这里的平衡态用的仍旧是 LBGK 模型的形式
> $$f^{(\mathrm{eq})}_i = w_i \rho \left[ 1 + \frac{\boldsymbol{c}_i \cdot \boldsymbol{u}}{c_s^2} + \frac{(\boldsymbol{c}_i \cdot \boldsymbol{u})^2}{2 c_s^4} - \frac{u^2}{2 c_s^2} \right]$$ 
> 其中 $\boldsymbol{c}_i$ 为离散速度方向， $w_i$ 为 $\boldsymbol{c}_i$ 方向的权重系数， $c_s$ 为格子声速。 $\rho$ 和 $\boldsymbol{u}$ 分别为流场密度和速度。
> 
> 在 D2Q9 模型中，$c_s=1/\sqrt{3}$，
> $$\boldsymbol{c}_{i} = \begin{cases} (0,0) & i=0 \\ (\pm 1, 0), (0, \pm 1) c & i=1-4 \\ (\pm 1, \pm 1) c & i=5-8 \end{cases}$$
> 权重系数 $w_0 = 4/9$ , $w_{1-4} = 1/9$ , $w_{5-8} = 1/36$ .

MRT 方程的碰撞步在D2Q9模型下的矩阵形式为：
$$
\boldsymbol{f}^{(\mathrm{pc})} = [\bold{M}^{-1}][\bold{S}][\bold{M}] (\boldsymbol{f} - \boldsymbol{f}^{(\mathrm{eq})})
\tag{1}
$$
松弛步则写作：
$$
f_i (\boldsymbol{x}+\boldsymbol{c}_i \delta_t, t + \delta_t) = f_i^{(\mathrm{pc})} (\boldsymbol{x}, t + \delta_t)
\tag{2}
$$
其中 $\delta_t$ 为时间步长， $[\bold{S}]$ 为松弛系数的对角矩阵。变换矩阵 $[\bold{M}]$ 为：
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

流场密度 $\rho$ 和速度 $\boldsymbol{u}$ 分别为
$$
\rho = \sum_i f_i , \quad \rho\boldsymbol{u} = \sum_i \boldsymbol{c}_i f_i
$$

郭照立等[1]记松弛时间的对角矩阵为：
$$[\bold{S}] = \text{diag}(0, s_e, s_{\varepsilon}, 0, s_q, 0, s_q, s_{\nu}, s_{\nu})$$
则运动粘度 $\nu = c_s^2 (\frac{1}{s_{\nu}} - \frac{1}{2}) \delta_t$ ，以及体黏性系数 $\zeta = c_s^2 (\frac{1}{s_e} - \frac{1}{2}) \delta_t$ 。

## Zhou 的简化

Zhou[2]将矩阵 $[\bold{S}]$ 改写成
$$
\begin{aligned}
    [\bold{S}] &= [\bold{S}_0] + [\bold{S}_t] \\
    &= \text{diag}(
        1/\tau_s, 1/\tau_s, 1/\tau_s, 1/\tau_s,
        1/\tau_s, 1/\tau_s, 1/\tau_s,
        1/\tau, 1/\tau
    ),
\end{aligned}
\tag{3}
$$
其中：
$$
\begin{aligned}
[\bold{S}_0] &= \left( \frac{1}{\tau} - \frac{1}{\tau_s} \right)
\text{diag}(0, 0, 0, 0, 0, 0, 0, 1, 1), \\
[\bold{S}_t] &= \frac{1}{\tau_s} \, \text{diag}(1, 1, 1, 1, 1, 1, 1, 1, 1).
\end{aligned}
$$
可见在这种构造形式下，$\tau= 1 / s_{\nu}$，与LBGK模型一致。

为了实现化简，Zhou[2]令 $\tau_s=1$，并把式(3)代入式(1)中展开，整理得：
$$
\begin{aligned}
f_i^{(\mathrm{pc})}(\boldsymbol{x},t) &= f_{i}^{(\mathrm{eq})}(\boldsymbol{x},t) + (-1)^i \left( 1 - \frac{1}{\tau} \right) \times \\
& \sum_{k=1}^{8} \frac{(-1)^k}{4}
\left[ \boldsymbol{1}_{[1,2,3,4]}(i,k) + \boldsymbol{1}_{[5,6,7,8]}(i,k) \right]
\left[ f_{k}(\boldsymbol{x},t) - f_{k}^{(\mathrm{eq})}(\boldsymbol{x},t) \right]
\end{aligned}
\tag{4}
$$
式(4)中的函数 $\boldsymbol{1}_{[a,b,c,d]}(i,k)$ 定义为
$$
\boldsymbol{1}_{[a,b,c,d]}(i,k) = \begin{cases}
    1, & \text{if i in [a,b,c,d], and k in [a,b,c,d]}\\
    0, & \text{otherwise}
\end{cases}
$$


# 实例代码

这里顺手做了个双周期剪切流的简单算例（拿上次双松弛LBM的代码改的），设置如下。

- 计算域： $[0,L] \times [0,L]$，其中 $L=2\pi$。
- 雷诺数： $Re = 1 / \nu = 1000$。 其中 $\nu$ 为流体运动粘度 (Physical)。
- 初始速度场
$$
\begin{cases}
    u_x =& U_0 \tanh (\sigma\, (\frac{1}{4} - |\frac{y}{L} - \frac{1}{2}| )) \\
    u_y =& U_0 \delta \sin (2 \pi \, (\frac{x}{L} + \frac{1}{2}) )
\end{cases}
$$
- 压强场 $p$ ： 用伪谱法求解压力泊松方程。
- LBM松弛时间 $\tau$: 0.515
- 网格大小： $\delta_x = L/N$ ($N=200$)。
- 时间步 $\delta_t$ ： 用 LBM 中 $\tau$ 的关系反推。
- 初始速度场参数： $U_0 = 1$, $\sigma=30$, $\delta=0.05$

代码我放在了[Github Gist](https://gist.github.com/Ives0721/bc73d950cfb4adb05900b27ef9c3b496)上，计算结果如下图所示

<img src="DoubleShear-sMRT.apng">

# Reference

[1] 郭照立，郑楚光. **格子Boltzmann方法的原理及应用**[M]. 科学出版社，2009.

[2] Zhou J G. **Single MRT-featured Lattice Boltzmann Method (SmrtLBM)**[PP/OL]. arXiv（2024-09-02）. http://arxiv.org/abs/2409.01076.
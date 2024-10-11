---
title: "【例题】球坐标系下的流体力学问题"
date: 2024-10-10 18:21:51
tags: [流体力学]
categories:
- [流体力学]
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

## 题 1

> 一无限流体被单位质量力 $\mu r^{-3/2}$ 的力作用，方向指向原点。若开始时流体静止且其中有一个空腔 $r=c$，证明：在过 $(\frac{2}{5 \mu})^{1/2} c^{5/4}$ 时间之后，空腔被流体填满。

这里假设流体不可压缩。取球坐标系 $(r, \theta, \phi)$ 进行计算，并以空腔中心为原点。

### 思路1：功能关系

假设当 $t=t_0$ 时，空腔边界的球壳形薄液面位于 $r=r_0$ 处。该薄液面厚度为 $\mathrm{d}r_0$, 液面速度 $v_r (r_0, t_0) = 0$。  

若 $t_0$ 时刻，位于 $r=r_0$ 处且厚度为 $\mathrm{d}r_0$ 的液面向内推动了 $\mathrm{d} r$， 由于该液面的质量积为 $\mathrm{d}m = \rho (4 \pi r_0^2) \mathrm{d}r_0$，则这一过程中体积力（$F_{br} = -\mu r^{-3/2}$）所作的功为

$$
\mathrm{d} W = \mathrm{d}m \cdot \int_{r_0}^{+\infty} F_{br} \mathrm{d}r = \mathrm{d}m \cdot (-2 \mu r_0^{-1/2}) = -8 \mu \rho \pi r_0^{3/2} \mathrm{d}r_0 .
$$

当 $t_1 = t_0 + \delta_t$ 时，该液面前进到 $r=r_1$ 位置（$r_1 \le r_0$），此时液面速度为 $v_r(r_1, t_1) = \dfrac{\mathrm{d} r_1}{\mathrm{d} t}$ 。

所以薄液面从 $r_0$ 运动到 $r_1$ 过程中， $F_{br}$ 做的总功为：

$$
W = \int_{r_0}^{r_1} \mathrm{d} W = \frac{16}{5} \mu \rho \pi (r_0^{5/2} - r_1^{5/2})
$$

记速度 $\bold{v} = [v_r, 0, 0]^T$。代入球坐标系下的连续性方程，化简得：

$$
\frac{\partial (r^2 v_r)}{\partial r} = 0
$$

对 $r$ 积分，可得 $r^2 v_r = f(t)$ ，其中 $f(t)$ 为待定函数。根据这一方程，可得在 $t_1$ 时刻满足 $r_0^2 v_r(r_0, t_1) = r_1^2 v_r(r_1, t_1)$ 。

在 $t_1$ 时刻， $r=r_0$ 层液体的动能表示为

$$
\mathrm{d} E = \frac{1}{2} \cdot \mathrm{d}m \cdot v_r^2(r_0, t_1) = 2 \rho \pi \cdot \frac{r_1^4}{r_0^2} \cdot v_r^2(r_1, t_1) \mathrm{d}r_0
$$

这一时刻下整个系统的动能为：
$$
E = \int_{r_1}^{+\infty} \mathrm{d} E = 2 \rho \pi r_1^3 \cdot v_r^2(r_1, t_1)
$$

据此，我们令 $t_0$ 时刻即为初始时刻（$t_0=0$），满足 $r_0 = c$ 和 $v_r(t=0)=0$。所以 $t_0$ 时刻的动能为 0。对 $t_0$ 到 $t_1$ 时刻的过程， $F_{br}$ 所做的功 $W$ 全部转化为流场动能 $E$，即：$W=E$。

$$
\frac{16}{5} \mu \rho \pi (c^{5/2} - r_1^{5/2}) = 2 \rho \pi r_1^3 \cdot v_r^2(r_1, t_1) - 0
$$

也就有：

$$
v_r(r_1, t_1) = \frac{\mathrm{d} r_1}{\mathrm{d} t} = \sqrt{\frac{8 \mu}{5}} \cdot \sqrt{r_1^{-3/2} (c^{5/2} - r_1^{5/2})}
$$

分离变量后解得：

$$
\frac{16}{25} (c^{5/2} - r_1^{5/2}) = \left( \sqrt{\frac{8 \mu}{5}} t + c_1 \right)^2
$$

其中 $c_1$ 为待定常数。当 $t=0$ 时， $r_1 = c$，所以 $c_1 = 0$。

假设当 $t=T_{m}$ 时 $r_1=0$， 则可解得： $T_m = (\frac{2}{5 \mu})^{1/2} c^{5/4}$。


### 思路2：联立动量方程和连续性方程

记速度 $\bold{v} = [v_r, 0, 0]^T$。

代入连续性方程，化简得：

$$
\frac{\partial (r^2 v_r)}{\partial r} = 0
$$

对 $r$ 积分，可得 $r^2 v_r = f(t)$ ，其中 $f(t)$ 为待定函数。

该场景下的运动方程为

$$
\rho \frac{\mathrm{D} \bold{v}}{\mathrm{D} t} = \rho \bold{F}_b
$$

所以，动量方程写作：

$$
\frac{\partial v_r}{\partial t} + v_r \frac{\partial v_r}{\partial r} =  F_{br}
$$

由题目可知 $F_{br} = - \mu r^{-3/2}$。将 $r^2 v_r = f(t)$ 代入到运动方程中，得：

$$
\frac{f'(t)}{r^2} + v_r \frac{\partial v_r}{\partial r} = -\mu r^{-3/2}
$$

记 $R(t)$ 为真空球面的半径变化，则 $v_r (r=R(t)) = R'(t)$。对 $r$ 积分，得到：
$$
\begin{aligned}
    \int_{\infty}^{R(t)} \frac{f'(t)}{r^2} + v_r \frac{\partial v_r}{\partial r} \mathrm{d}r &= -\mu \int_{\infty}^{R(t)} r^{-3/2} \mathrm{d}r \\
    \left.\left[ -\frac{f'(t)}{r} + \frac{(v_r)^2}{2} \right]\right|_{\infty}^{R(t)} &=
    \left. \frac{2 \mu}{\sqrt{r}} \right|_{\infty}^{R(t)}
\end{aligned}
$$

代入 $v_r (r=\infty)=0$，整理得到： 
$$
-\frac{f'(t)}{R(t)} + \frac{(R'(t))^2}{2} = \frac{2 \mu}{\sqrt{R(t)}}
$$

在边界上有 $f(t) = R^2(t) \cdot v_r (r=R(t)) = R^2 \cdot R'$ ， 所以 $f'(t) = 2 R R'^{2} + R^2 \cdot R{''}$ 。上式化简为：

$$
-\frac{3}{2} R'^{2} -R \cdot R{''} = \frac{2 \mu}{\sqrt{R}}
$$


由于 $\frac{\mathrm{d} (R^{1/2} (R')^{2})}{\mathrm{d} R} =\frac{R^{-1/2}}{2} (R')^{2} + 2 R^{1/2} R{''}$，上式整理得： 
$$
-\frac{\mathrm{d}(R^{1/2} (R')^{2})}{8\mu + 5 (R^{1/2} (R')^{2})} = \frac{\mathrm{d}R}{2R} 
$$

分离变量后并整理出 $R'$ ，可表示为：

$$
R' = \sqrt{\frac{R^{-1/2}}{5} (c_1 R^{-5/2} - 8 \mu) }
$$

将 $R(t=0)=c$ 和 $R'(t=0) = 0$ 代入，解得： $c_1 = 8 \mu c^{5/2}$ 。即：

$$
R' = \sqrt{\frac{R^{-1/2}}{5} 8 \mu (c^{5/2} R^{-5/2} - 1)} 
\\= \sqrt{\frac{8 \mu}{5}} \cdot \sqrt{R^{-1/2} (c^{5/2} R^{-5/2} - 1)}
$$

综上，真空球面的总运动时长表示为： 

$$
T = \int_{0}^{c} \frac{1}{R'} \mathrm{d}R = \sqrt{\frac{5}{8 \mu}} \int_{0}^{c} \left[ R^{-1/2} \left( \left(\frac{c}{R}\right)^{5/2} - 1 \right) \right]^{-1/2} \mathrm{d}R
$$

因为：

$$
\begin{aligned}
    & \int_{0}^{c} \left[ R^{-1/2} \left( \left(\frac{c}{R}\right)^{5/2} - 1 \right) \right]^{-1/2} \mathrm{d}R \\
    =& \int_{0}^{c} R^{1/4} \left[ \left( \left(\frac{c}{R}\right)^{5/2} - 1 \right) \right]^{-1/2} \mathrm{d}R \\
    =& \frac{4}{5} \int_{0}^{c}  \left[ \left( \left(\frac{c}{R}\right)^{5/2} - 1 \right) \right]^{-1/2} \mathrm{d}(R^{5/4}) \\
    \xlongequal{t=R^{5/4}} & \frac{4}{5} \int_{0}^{c^{5/4}}  \left[ \left( c^{5/2} t^{-2} - 1 \right) \right]^{-1/2} \mathrm{d}t 
    \quad= \frac{4}{5} \int_{0}^{c^{5/4}} t \left[ \left( c^{5/2} - t^2 \right) \right]^{-1/2} \mathrm{d}t \\
    =& \frac{4}{5} \cdot \left(-\frac{1}{2}\right) \int_{0}^{c^{5/4}} \left[ \left( c^{5/2} - t^2 \right) \right]^{-1/2} \mathrm{d}(c^{5/2} - t^2) \\
    =& \frac{4}{5} \cdot \left(-\frac{1}{2}\right) \left.2 (c^{5/2} - t^2)^{1/2}\right|_{0}^{c^{5/4}} = \frac{4}{5} c^{5/4}
\end{aligned}
$$

所以 $T = \sqrt{\frac{5}{8 \mu}} \cdot \frac{4}{5} c^{5/4} = \sqrt{\frac{2}{5 \mu}} c^{5/4}$ .

证毕。

<!--
https://zhuanlan.zhihu.com/p/692839599
-->

---
title: "流体在球坐标系下的控制方程"
date: 2024-10-10 08:58:51
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

<details><summary>参考资料：</summary>
  
1. https://pages.mtu.edu/~fmorriso/cm310/Navier.pdf  
2. https://www.bilibili.com/read/cv16068248/  
3. https://zhuanlan.zhihu.com/p/550128553  
4. https://zhuanlan.zhihu.com/p/692839599  
5. https://www.tessshebaylo.com/navier-stokes-equation-derivation-in-spherical-coordinates/  
6. https://texmex.mit.edu/pub/emanuel/CLASS/12.340/navier-stokes%282%29.pdf  
7. https://ocw.mit.edu/courses/8-01sc-classical-mechanics-fall-2016/mit8_01scs22_chapter30.pdf  

</details>


## 球坐标系的控制体

![球坐标系的微元示意图](https://www.researchgate.net/profile/Gustav-Rasmussen/publication/314913282/figure/fig27/AS:471590883336194@1489447436698/Spherical-coordinates-r-th-and-ph-The-3-boldtype-coordinate-axes-represent-the-Cartesian.png)


记上面控制体的八个点分别为  
$$
\begin{aligned}
    A=(R, \theta, \phi) ,\quad &
    B=(R, \theta + \mathrm{d}\theta, \phi) \\
    C=(R + \mathrm{d}r, \theta + \mathrm{d}\theta, \phi),\quad &
    D=(R + \mathrm{d}r, \theta, \phi) \\
    E=(R, \theta, \phi+\mathrm{d}\phi),\quad &
    F=(R, \theta + \mathrm{d}\theta, \phi+\mathrm{d}\phi) \\
    G=(R + \mathrm{d}r, \theta + \mathrm{d}\theta, \phi+\mathrm{d}\phi),\quad &
    H=(R + \mathrm{d}r, \theta, \phi+\mathrm{d}\phi)
\end{aligned}
$$

若点 A 的球坐标 $(R, \theta, \phi)$ 对应笛卡尔坐标 $(x,y,z)$，则：  
$$
\begin{cases}
x =& R \sin\theta \cos\phi \\
y =& R \sin\theta \sin\phi \\
z =& R \cos\theta \\
\end{cases} ,\quad
\begin{bmatrix}
    \bold{e}_r \\ \bold{e}_\theta \\ \bold{e}_\phi
\end{bmatrix} =
\begin{bmatrix}
    \sin\theta \cos\phi & \sin\theta \sin\phi & \cos\theta \\
    \cos\theta \cos\phi & \cos\theta \sin\phi & -\sin\theta \\
    -\sin\theta & \cos\theta & 0\\
\end{bmatrix}
\begin{bmatrix}
    \bold{e}_x \\ \bold{e}_y \\ \bold{e}_z
\end{bmatrix}
$$

其中 $[\bold{e}_r, \bold{e}_\theta, \bold{e}_\phi]^T$ 和 $[\bold{e}_x, \bold{e}_y, \bold{e}_z]^T$ 分别为球坐标系和笛卡尔坐标系各坐标方向的单位法向量。 

* 微元的体积为 $\mathrm{dV} = R^2 \sin\theta \mathrm{d}r \mathrm{d}\phi \mathrm{d}\theta$ 。  
* 球坐标系中的流场速度矢量表示为 $\bold{v} = v_R \bold{e}_r + v_{\theta} \bold{e}_\theta + v_{\phi} \bold{e}_\phi$。  
* 球坐标可使用笛卡尔坐标系来表示： $R = \sqrt{x^2 + y^2 + z^2}$； $\theta = \arctan \frac{\sqrt{x^2 + y^2}}{z}$ ； $\phi = \arctan \frac{y}{x}$。  
* 球坐标中的矢量偏导： $\displaystyle \mathrm{d}\bold{x} = \bold{e}_r \cdot \mathrm{d}R + \frac{\mathrm{d}\theta}{R} \bold{e}_\theta + \frac{\mathrm{d} \phi}{R \sin\theta} \bold{e}_\phi$  
* 球坐标系下的散度为
$$
\begin{aligned}
\nabla\cdot\bold{v} =& \mathrm{div}\,\bold{v} \\
&= \frac{1}{\mathrm{dV}} \left\{ \frac{\partial}{\partial R}[((R \sin\theta \mathrm{d}\phi)\cdot(R \mathrm{d}\theta)) v_R] \mathrm{d}R + \right.
\\&\quad \left. \frac{\partial}{\partial \theta}[((R \sin\theta \mathrm{d}\phi)\cdot\mathrm{d} R) v_\theta] \mathrm{d}\theta + 
\frac{\partial}{\partial \phi}[(R \mathrm{d}\theta \cdot\mathrm{d} R) v_\phi] \mathrm{d}\phi \right\}
\\&= 
\frac{1}{R^2} \frac{\partial (R^2 v_R)}{\partial R} + 
\frac{1}{R \sin\theta} \frac{\partial (v_\theta \sin\theta)}{\partial \theta} +
\frac{1}{R \sin\theta} \frac{\partial (v_\phi)}{\partial \phi} 
\end{aligned}
$$

## 连续性方程

控制体内部在 $\Delta t$ 时间段的质量变化表示为 $\displaystyle \frac{\partial \rho}{\partial t} \mathrm{dV} \cdot \Delta t$ ，下文对所有面进行分析：

（1）面 ABCD 的流入为 $\rho v_{\phi} ((R \mathrm{d}\theta) \cdot (\mathrm{d}R)) \cdot \Delta t$ ，面 EFGH 的流出为 $\displaystyle \rho v_{\phi} R \mathrm{d}\theta \mathrm{d}R \cdot \Delta t + \frac{\partial}{\partial \phi}\left( \rho v_{\phi} R \mathrm{d}\theta \mathrm{d}R \right) \mathrm{d}\phi \Delta t$。所以对于 ABCD 和 EFGH 两个面，其净流出为：  
$$
\begin{aligned}
& \frac{\partial}{\partial \phi}\left( \rho v_{\phi} R \mathrm{d}\theta \mathrm{d}R \right) \mathrm{d}\phi \Delta t \\=&
\frac{\partial}{\partial \phi}\left( \rho v_{\phi} R \right) \mathrm{d}R \mathrm{d}\phi \mathrm{d}\theta \cdot \Delta t .
\end{aligned}
$$

（2）面 ABFE 的流入为 $\rho v_{R} ((R \mathrm{d}\theta) \cdot (R \sin\theta \mathrm{d}\phi)) \cdot \Delta t$ ，面 DCGH 的流出为 $\displaystyle \rho R^2 v_{R} \sin\theta \mathrm{d}\theta \mathrm{d}\phi \cdot \Delta t + \frac{\partial}{\partial R}\left( \rho R^2 v_{R} \sin\theta \mathrm{d}\theta \mathrm{d}\phi \right) \mathrm{d}R \Delta t$ 。所以对于 ABFE 和 DCGH 两个面，其净流出为：  
$$
\begin{aligned}
& \frac{\partial}{\partial R}\left( \rho R^2 v_{R} \sin\theta \mathrm{d}\theta \mathrm{d}\phi \right) \mathrm{d}R \Delta t \\=&
\frac{\partial}{\partial R}\left( \rho R^2 v_{R} \sin\theta  \right) \mathrm{d}R \mathrm{d}\phi \mathrm{d}\theta \cdot \Delta t .
\end{aligned}
$$

（3）面 AEHD 的流入为 $\rho v_{\theta} (\mathrm{d}R \cdot (R \sin\theta \mathrm{d}\phi)) \cdot \Delta t$ ，面 BFGC 的流出为 $\displaystyle \rho R v_{\theta} \sin\theta \mathrm{d}R \mathrm{d}\phi \cdot \Delta t + \frac{\partial}{\partial \theta}\left( \rho R v_{\theta} \sin\theta \mathrm{d}R \mathrm{d}\phi \right) \mathrm{d}\theta \cdot \Delta t$ 。所以对于 AEHD 和 BFGC 两个面，其净流出为：  
$$
\begin{aligned}
    & \frac{\partial}{\partial \theta}\left( \rho R v_{\theta} \sin\theta \mathrm{d}R \mathrm{d}\phi \right) \mathrm{d}\theta \cdot \Delta t \\=&
    \frac{\partial}{\partial \theta} (\rho R v_{\theta} \sin\theta) \mathrm{d}R \mathrm{d}\phi \mathrm{d}\theta \cdot \Delta t .
\end{aligned}
$$

因此，所有面的净流出量为
$$
\left[\frac{\partial}{\partial \phi}\left( \rho v_{\phi} R \right) + 
\frac{\partial}{\partial R}\left( \rho R^2 v_{R} \sin\theta  \right) +
\frac{\partial}{\partial \theta} (\rho R v_{\theta} \sin\theta)
\right] \mathrm{d}R \mathrm{d}\phi \mathrm{d}\theta \cdot \Delta t .
$$

**控制体内部在 $\Delta t$ 时间段的质量变化，等于所有面的净流入量，此即连续性方程**。整理可得：

$$
(R^2 \sin\theta) \frac{\partial \rho}{\partial t} + 
\left[\frac{\partial}{\partial \phi}\left( \rho v_{\phi} R \right) +
\frac{\partial}{\partial R}\left( \rho R^2 v_{R} \sin\theta  \right) + 
\frac{\partial}{\partial \theta} (\rho R v_{\theta} \sin\theta)
\right] = 0
$$

对于不可压缩流体，满足 $\frac{\partial \rho}{\partial t} = 0$ 和 $\frac{\partial \rho}{\partial R}=\frac{\partial \rho}{\partial \theta}=\frac{\partial \rho}{\partial \phi}=0$，则上式化简为：

$$
\frac{\partial (R^2 v_{R})}{\partial R} \sin\theta + R \frac{\partial v_\phi}{\partial \phi} +  R \frac{\partial (v_{\theta} \sin\theta)}{\partial \theta} = 0
$$


## 运动方程

### 控制体内的动量变化

由上一节的推导，我们可以类比得到：控制体内部的动量变化为 $\displaystyle \frac{\partial (\rho \bold{v})}{\partial t} \mathrm{dV} \cdot \Delta t$，以及所有面动量净流出量的总和  
$$
\left[\frac{\partial}{\partial \phi}\left( \rho v_{\phi} R \bold{v} \right) +
\frac{\partial}{\partial R}\left( \rho R^2 v_{R} \bold{v} \sin\theta \right) + 
\frac{\partial}{\partial \theta} (\rho R v_{\theta} \bold{v} \sin\theta )
\right] \mathrm{d}R \mathrm{d}\phi \mathrm{d}\theta \cdot \Delta t
$$

所以，控制体内的动量变化为

$$
\left[\frac{\partial (\rho \bold{v})}{\partial t} R^2 \sin\theta + \frac{\partial}{\partial \phi}\left( \rho v_{\phi} R \bold{v} \right) +
\frac{\partial}{\partial R}\left( \rho R^2 v_{R} \bold{v} \sin\theta \right) + 
\frac{\partial}{\partial \theta} (\rho R v_{\theta} \bold{v} \sin\theta)
\right] \mathrm{d}R \mathrm{d}\phi \mathrm{d}\theta \cdot \Delta t
$$

由于球坐标系内 $\bold{x}=[R, \theta, \phi]^T$， $\frac{\partial \bold{x}}{\partial t} = [\frac{\partial R}{\partial t}, \frac{\partial \theta}{\partial t}, \frac{\partial \phi}{\partial t}]^T = [v_R,  v_\theta / R, v_\phi / (R \sin\theta)]^T$， 则球坐标系内的随体导数为 

$$
\frac{\mathrm{D}}{\mathrm{D} t} = \frac{\partial}{\partial t} + \frac{\partial}{\partial \bold{x}} \cdot \frac{\partial \bold{x}}{\partial t} =
\frac{\partial}{\partial t} + (v_R \frac{\partial}{\partial R} + \frac{v_\theta}{R} \frac{\partial}{\partial \theta} + \frac{v_\phi}{R \sin\theta} \frac{\partial}{\partial \phi})
$$

那么对于任意量 $\xi$， $(R^2 \sin\theta) \frac{\mathrm{D \xi}}{\mathrm{D} t}$ 则写作：

$$
(R^2 \sin\theta) \frac{\mathrm{D \xi}}{\mathrm{D} t} = 
\frac{\partial \xi}{\partial t} R^2 \sin\theta + (\frac{\partial \xi}{\partial R} v_R R^2 \sin\theta + \frac{\partial \xi}{\partial \theta} v_\theta R \sin\theta + v_\phi R \frac{\partial \xi}{\partial \phi})
$$


所以，**控制体内的动量变化**被进一步改写为：

$$
(R^2 \sin\theta)\left[
    \frac{\mathrm{D (\rho \bold{v})}}{\mathrm{D} t} + \rho \bold{v} \cdot (\nabla\cdot\bold{v})
\right] \mathrm{d}R \mathrm{d}\phi \mathrm{d}\theta \cdot \Delta t
$$

根据连续性方程，有 $\frac{\mathrm{D \rho}}{\mathrm{D} t} + \rho\, (\nabla\cdot\bold{v})=0$，则**控制体内的动量变化**为： $\displaystyle \rho \frac{\mathrm{D \bold{v}}}{\mathrm{D} t} \mathrm{dV} \cdot \Delta t$


### 控制体各个面所受到的外力作用

流体控制体所受到的外力分为两种作用：体力和面力。

（1）对于体力，假设体力的大小为 $\rho \bold{F}_b$。则体力在 $\Delta t$ 时间段内的作用为： $\displaystyle \rho \bold{F}_b \mathrm{dV} \Delta t$ 。

其中 $\bold{F}_b = [F_{br}, F_{b\theta}, F_{b\phi}]^T$。

（2）对于面力，记 $r$ , $\theta$, $\phi$ 方向的面力向量分别为 $\boldsymbol{p}_r$, $\boldsymbol{p}_\theta$, $\boldsymbol{p}_\phi$。下面对体积为 $\mathrm{dV}$ 的控制体进行分析。


👉 ABCD 和 EFGH 面的受力：  
$$
\frac{\partial (\boldsymbol{p}_\phi \cdot R \mathrm{d}\theta \cdot \mathrm{d}R)}{\partial \phi}  \cdot (R \sin\theta \mathrm{d}\phi)
$$

👉 AEHD 和 BFGC 面的受力：  
$$
\frac{\partial (\boldsymbol{p}_\theta \cdot  R \sin\theta \mathrm{d}\phi \cdot \mathrm{d}R)}{\partial \theta} \cdot (R \mathrm{d}\theta)
$$

👉 ABFE 和 CGHD 面的受力：  
$$
\frac{\partial (\boldsymbol{p}_r \cdot R \sin\theta \mathrm{d}\phi \cdot R \mathrm{d}\theta)}{\partial r}  \cdot \mathrm{d}R
$$

所以面力的和为：
$$
\begin{aligned}
    & \left[
        R \frac{\partial \boldsymbol{p}_\phi}{\partial \phi} +
        R \frac{\partial (\boldsymbol{p}_\theta \sin\theta)}{\partial \theta} +
        \sin\theta \frac{\partial (R^2 \boldsymbol{p}_r)}{\partial r}
    \right]
    \mathrm{d}R \mathrm{d}\phi \mathrm{d}\theta \\
    =& \left[
        \frac{1}{R \sin\theta} \frac{\partial \boldsymbol{p}_\phi}{\partial \phi} +
        \frac{1}{R \sin\theta} \frac{\partial (\boldsymbol{p}_\theta \sin\theta)}{\partial \theta} +
        \frac{1}{R^2} \frac{\partial (R^2 \boldsymbol{p}_r)}{\partial r}
    \right] \mathrm{dV} \\
    =& \mathrm{div}\,\boldsymbol{P} \cdot \mathrm{dV}
\end{aligned}
$$

其中 $\boldsymbol{P}=[\boldsymbol{p}_r, \boldsymbol{p}_\theta, \boldsymbol{p}_\phi]^T$ 。综上所述，运动方程为  
$$
\rho \frac{\mathrm{D \bold{v}}}{\mathrm{D} t} = \rho \bold{F}_{b} + \mathrm{div}\,\boldsymbol{P}
$$

## 球坐标系下的 Navier-Stokes 方程


运动方程的**等号左侧**在球坐标系展开为：

$$
\rho \left[ 
\frac{\partial \bold{v}}{\partial t} + (v_R \frac{\partial \bold{v}}{\partial R} + \frac{v_\theta}{R} \frac{\partial \bold{v}}{\partial \theta} + \frac{v_\phi}{R \sin\theta} \frac{\partial \bold{v}}{\partial \phi}) \right] 
$$


运动方程的**等号右侧**在球坐标系展开为：
$$
\rho \bold{F}_{b} + 
\left[
    \frac{1}{R \sin\theta} \frac{\partial \boldsymbol{p}_\phi}{\partial \phi} +
    \frac{1}{R \sin\theta} \frac{\partial (\boldsymbol{p}_\theta \sin\theta)}{\partial \theta} +
    \frac{1}{R^2} \frac{\partial (R^2 \boldsymbol{p}_r)}{\partial r}
\right]
$$

若需要获得 Navier-Stokes 方程，则需要获得 $\boldsymbol{p}_r, \boldsymbol{p}_\theta, \boldsymbol{p}_\phi$ 在球坐标系下的本构方程。

$$
\boldsymbol{p}_r = \begin{bmatrix}
    p_{rr} \\ p_{r \theta} \\ p_{r \phi}
\end{bmatrix} ,
\boldsymbol{p}_\theta = \begin{bmatrix}
    p_{\theta r} \\ p_{\theta \theta} \\ p_{\theta \phi}
\end{bmatrix} ,
\boldsymbol{p}_\phi = \begin{bmatrix}
    p_{\phi r} \\ p_{\phi \theta} \\ p_{\phi \phi}
\end{bmatrix}
$$
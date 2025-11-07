---
title: 【笔记】格子Bolztmann方法中单点局部边界条件的实现
date: 2025-11-07 09:28:57
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
.container {
    display: flex;
}
.col, p {
    flex: 1;
    font-family: "Arial";
    text-align: left;
}
</style>


## Lattice Boltzmann method

The (force-free) LBM equation is

$$
f_{i}(\boldsymbol{x} + \boldsymbol{c}_{i}, t+\delta_t) - f_{i}(\boldsymbol{x}, t) = \Omega_{i}(\boldsymbol{x}, t).
$$

In LBGK model, $\Omega_{i}(\boldsymbol{x}, t) = -\frac{1}{\tau} \left[ f_{i}(\boldsymbol{x}, t) - f_{i}^{(eq)} \right]$. 

The $f_{i}^{(eq)}$ is the equilibrium state,

$$
f_{i}^{(eq)} = w_{i} \rho \left[
    1 + \frac{\boldsymbol{c}_{i} \cdot \boldsymbol{u}}{c_s^2}
    + \frac{(\boldsymbol{c}_{i} \cdot \boldsymbol{u})^2}{2 c_s^4}
    - \frac{|\boldsymbol{u}|^2}{2 c_s^2}
\right].
$$

$w_i$ is the weight of the $\boldsymbol{c}_{i}$. $c_s$ is the speed of sound. $\rho, \boldsymbol{u}$ are the fluid density and velocity, respectively.

For simplicity, we define $f^{*}_{i} = f_{i} - \Omega_{i}$ as the post-collision population.


## Curved boundary in LBM

<div class="container">
<div class="col">
    <img src="Fig01.png" width="50%">
</div>

<div class="col">

The purpose of boundary algorithms is to reconstruct missing $f_{\bar{i}}(\boldsymbol{x}_{F}, t+\delta_t)$ after the streaming step. 

Interpolated bounce-back method (IBB) ia a general ideas to implement the boundary condition. We will introduce it at below.

</div>
</div>

> NOTE:  ① In the left figure, $\boldsymbol{x}_{w} = \boldsymbol{x}_{F} + q \boldsymbol{c}_{\bar{i}}$, and $\boldsymbol{c}_{\bar{i}} = -\boldsymbol{c}_{i}$. ② Let the yellow node is $\vec{x}_b$, then $q = |\vec{x}_w - \vec{x}_F| / |\vec{x}_b - \vec{x}_F|$ 



### Interpolated bounce-back method
It treats the unknown $f_{i}(\boldsymbol{x}_{F}, t+\delta_t)$ as a polynomial of the known populations.

$$
f_{i}(\boldsymbol{x}_{F}, t+\delta_t) =
a_1 f_{\bar{i}}^{*}(\boldsymbol{x}_{FF}, t) + 
a_2 f_{\bar{i}}^{*}(\boldsymbol{x}_{F}, t) + a_3 f_{i}^{*}(\boldsymbol{x}_{F}, t) + K
$$

where $K$ is a correction term.

<br>

> **|👉 How to find the coefficients $a_i$?**
> (1) Use Chapman-Enskog expansion, Taylor expansion, or Grad's method to find the coefficients.
> (2) BFL: a mesoscopic, geometrical approach proposed by Bouzidi et al.


### BFL scheme

<img src="BFL-sketch-2.png">


> It determinates the boundary scheme based on the value of $q$. The length of the orange line is 1.


When $q \ge \frac{1}{2}$ (Figure (a)), we have

$$
f_{i}(\boldsymbol{x}_{+}, t+\delta_t) = f_{\bar{i}}^{*} (\boldsymbol{x}_{F}, t)
,\quad
f_{i}(\boldsymbol{x}_{FF}, t+\delta_t) = f_{i}^{*} (\boldsymbol{x}_{F}, t).
$$

Then we can do linear interpolation to calculate $f_{i}(\boldsymbol{x}_{F}, t+\delta_t)$:

$$
\begin{aligned}
f_{i}(\boldsymbol{x}_{F}, t+\delta_t) =& f_{i}(\boldsymbol{x}_{+}, t+\delta_t) + \\
&\frac{2q - 1}{2q} [f_{i}(\boldsymbol{x}_{FF}, t+\delta_t)- f_{i}(\boldsymbol{x}_{+}, t+\delta_t)] \\
=& \frac{1}{2q} f_{\bar{i}}^{*} (\boldsymbol{x}_{F}, t) + (1-\frac{1}{2q}) f_{i}^{*} (\boldsymbol{x}_{F}, t)
\end{aligned}
$$

When $q < \frac{1}{2}$ (Figure (b)), we can do linear interpolation, i.e.,

$$
\begin{aligned}
f_{i}(\boldsymbol{x}_{F}, t+\delta_t) =& f_{\bar{i}}^{*}(\boldsymbol{x}_{+}, t)\\
=& f_{\bar{i}}^{*}(\boldsymbol{x}_{F}, t) \\
& + (1-2q) [f_{\bar{i}}^{*}(\boldsymbol{x}_{FF}, t) - f_{\bar{i}}^{*}(\boldsymbol{x}_{F}, t)]\\
=& 2q f_{\bar{i}}^{*}(\boldsymbol{x}_{F}, t) + (1-2q) f_{\bar{i}}^{*}(\boldsymbol{x}_{FF}, t)
\end{aligned}
$$



## One-node local boundary


We can notice that: the general scheme usually employs the data of $\boldsymbol{x}_{FF}$, which is the second layer of fluid. 


> To design local linkwise boundary condition, one must discard the nonlocal contribution $f^{*}(\boldsymbol{x}_{FF}, t)$ that appears in the interpolation scheme.
> 
> Because the existence of $f^{*}(\boldsymbol{x}_{FF}, t)$ would make the scheme hard to describe the narrow region or corner.


### 1st order time approximation


<img src="Zhao等-sketch.png">


> The boundary node $\vec{x}_w$ is at the middle of $\vec{x}_1$ and $\vec{x}_2$, which is easier to use the half-way bounce-back [(*J. Fluid Mech.*, 271 (1994), pp. 285–309)](https://doi.org/10.1017/S0022112094001771).

The first one is based on the following **first-order in time approximation:**

$$
\underbrace{f_{\bar{i}}^{*}(\boldsymbol{x}_{FF}, t) = f_{\bar{i}}(\boldsymbol{x}_{F}, t+\delta_t) }_{\text{LBM streaming step}}
\approx
\underbrace{f_{\bar{i}}(\boldsymbol{x}_{F}, t)}_{\text{Approximation}}
$$

Zhao et al.[3] point out that: their scheme can obtain second-order accuracy under the diffusive scaling ($\delta_t = O(\delta_x^2)$) [3].

Here: 

$\vec{x}_w = \vec{x}_{F} + q \vec{c}_{\bar{i}}$, $\vec{x}_1 = \vec{x}_{F} + l \vec{c}_{\bar{i}}$, $\vec{x}_2 = 2 \vec{x}_{w} - \vec{x}_1$. 

Stability condition: $l \in [\max\{0, 2q-1\}, 2q]$



With $\vec{x}_1$ and $\vec{x}_{FF}$, the interpolation of $\vec{x}_{F}$ is:

$$
\begin{aligned}
f_i (\vec{x}_{F}, t + \delta_t) =& \frac{l}{1+l} f_i (\vec{x}_{FF}, t + \delta_t) + \frac{1}{1+l} f_i (\vec{x}_{1}, t + \delta_t) \\
=& \underbrace{\frac{l}{1+l} f_i^{*} (\vec{x}_{F}, t)}_{\text{LBM Streaming}} + \frac{1}{1+l} f_i (\vec{x}_{1}, t + \delta_t)
\end{aligned}
$$

For the unknown distribution function $f_i (\vec{x}_{1}, t + \delta_t)$, Zhao et al.[3] use the half-way bounce-back purposed by Ladd [(*J. Fluid Mech.*, 271 (1994), pp. 285–309)](https://doi.org/10.1017/S0022112094001771).

$$
f_i (\vec{x}_{1}, t + \delta_t) =
f_{\bar{i}} (\vec{x}_{2}, t + \delta_t) + 2 w_i \rho_0 \frac{\boldsymbol{c}_i \cdot \boldsymbol{u}(\vec{x}_b, t)}{c_s^2}
$$

Here the wall velocity $\boldsymbol{u}(\vec{x}_b, t)$ is a known quantity.

Next, we interpolate $f_{\bar{i}} (\vec{x}_{2}, t + \delta_t)$ with the distribution functions at $\vec{x}_F$ and $\vec{x}_b$:

$$
f_{\bar{i}} (\vec{x}_{2}, t + \delta_t) = 
(1+l-2q) f_{\bar{i}} (\vec{x}_{F}, t + \delta_t) +
(2q-l) f_{\bar{i}} (\vec{x}_{b}, t + \delta_t)
$$

Consider the streaming step from $\vec{x}_F$ and $\vec{x}_b$, we have: $f_{\bar{i}}^{*} (\vec{x}_{F}, t) = f_{\bar{i}} (\vec{x}_{b}, t + \delta_t)$.

Finally, with the approximation

$$
\underbrace{f_{\bar{i}}^{*}(\boldsymbol{x}_{FF}, t) = f_{\bar{i}}(\boldsymbol{x}_{F}, t+\delta_t) }_{\text{LBM streaming step}}
\approx
\underbrace{f_{\bar{i}}(\boldsymbol{x}_{F}, t)}_{\text{Approximation}}, 
$$

we can obtain the local boundary scheme:

$$
f_i (\vec{x}_{F}, t + \delta_t) = 
\underbrace{
    \frac{1+l-2q}{1+l} f_{\bar{i}} (\vec{x}_F, t)
}_{\text{before collision}} +
\underbrace{
    \frac{l}{1+l} f_{i}^{*}(\boldsymbol{x}_{F}, t) + 
    \frac{2q-l}{1+l} f_{\bar{i}}^{*}(\boldsymbol{x}_{F}, t)
}_{\text{post collision terms}} + 
\underbrace{\frac{2 w_i \rho_0}{(1+l)c_s^2} \boldsymbol{c}_i \cdot \boldsymbol{u}(\vec{x}_b, t)}_{\text{Ladd's bounce-back}}
$$


> **[NOTE]**
> (1) This scheme doesn't rely on the collision model.
> (2) If $\boldsymbol{x}_{FF}$ is always available, we can remove the first-order time approximation, and convert the scheme into a two-node scheme.
> (3) When $l = q$, the scheme is equivalent to that purposed by Geier et al. [(*Comput. Math. Appl.*, 70 (2015), pp. 507–547)](https://doi.org/10.1016/j.camwa.2015.05.001)

$$
f_i (\vec{x}_{F}, t + \delta_t) = 
\frac{1-q}{1+q} f_{\bar{i}} (\vec{x}_F, t) +
\frac{q}{1+q} [f_{i}^{*}(\boldsymbol{x}_{F}, t) + 
f_{\bar{i}}^{*}(\boldsymbol{x}_{F}, t)] + 
\frac{2 w_i \rho_0}{(1+q)c_s^2} 
\boldsymbol{c}_i \cdot \boldsymbol{u}(\vec{x}_b, t)
$$

If we still use $f_{\bar{i}}^{*} (\vec{x}_{FF}, t)$ rather than replacing it with $f_{\bar{i}} (\vec{x}_{F}, t)$, the scheme is equivalent to that purposed by Yu et al. [(*AIAA Paper*, 2003-0953, 2003.)](https://doi.org/10.2514/6.2003-953).

### Non-equilibrium state reconstruction

Tao et al.[4] have proposed a scheme to reconstruct the equilibrium and non-equilibrium state of the fluid node near the wall.

For the boundary node $\vec{x}_F$, its unknown distribution function can be calculated by the following interpolation:

$$
f_{i}(\vec{x}_F, t + \delta_t) = \frac{1}{1+q} f_{i}(\vec{x}_w, t + \delta_t) +
\frac{q}{1+q} f_{i}(\vec{x}_{FF}, t + \delta_t).
$$

Then we begin to discuss the reconstruction process.


First, consider the streaming step from $\vec{x}_F$ to $\vec{x}_{FF}$, we have:

$$
f_{i}^{*}(\vec{x}_{F}, t) = f_{i}(\vec{x}_{FF}, t + \delta_t)
$$


Second, the term $f_{i}(\vec{x}_w, t + \delta_t)$ can be divided into two parts of equilibrium ($f_{i}^{(eq)}$) and non-equilibrium states ($f_{i}^{(neq)}$):

$$
\begin{aligned}
f_{i}(\vec{x}_w, t + \delta_t) &= f_{i}^{(eq)}(\vec{x}_w, t + \delta_t) + f_{i}^{(neq)}(\vec{x}_w, t + \delta_t) \\
&= f_{i}^{(eq)}(\vec{x}_w, t + \delta_t) + f_{i}^{(neq)}(\vec{x}_F, t + \delta_t) \\
&= f_{i}^{(eq)}(\vec{x}_w, t + \delta_t) + f_{\bar{i}}^{(neq)}(\vec{x}_F, t)
\end{aligned}
$$

> Consider the substitution from $f_{i}^{(neq)}(\vec{x}_F, t + \delta_t)$ to $f_{\bar{i}}^{(neq)}(\vec{x}_F, t)$ as a non-equilibrium bounce-back, **as mentioned in Tao et al.[4].**

The $f_{i}^{(eq)}(\vec{x}_w, t + \delta_t)$ can be determined by the known velocity and the approximate fluid density 

$$
\begin{aligned}
f_{i}^{(eq)}(\vec{x}_w, t + \delta_t) &= f_{i}^{(eq)} \left( \vec{u}(\vec{x}_w, t + \delta_t), \rho(\vec{x}_w, t + \delta_t) \right) \\
&\approx f_{i}^{(eq)} \left( \vec{u}(\vec{x}_w, t + \delta_t), \rho(\vec{x}_F, t + \delta_t) \right) \\
&\approx f_{i}^{(eq)} \left( \vec{u}(\vec{x}_w, t + \delta_t), \rho(\vec{x}_F, t) \right) \\
\end{aligned}
$$


Here we cite the original sentences from the paper [4]:

> It has been demonstrated that for low speed flow, using $\rho_A (t)$ and $\rho_A (t+\delta_t)$ to approximate $\rho_A (t+\delta_t)$ and $\rho_W (t+\delta_t)$ have second- and third-order accuracies, respectively [23,39].
>
> **== Reference ==**
> [23] Z. Guo, C. Zheng, B. Shi, An extrapolation method for boundary conditions in lattice Boltzmann method, Phys. Fluids 14 (6) (2002) 2007–2010.
> [39] Z. Guo, C. Shu, Lattice Boltzmann Method and Its Applications in Engineering, World Scientific, 2013.

<br>

So the boundary scheme purposed by Tao et al. [4] is:

$$
f_{i}(\vec{x}_F, t + \delta_t) = 
\frac{q}{1+q} f_{i}^{*}(\vec{x}_{F}, t) + 
\frac{1}{1+q} \left\{ f_{i}^{(eq)} \left( \vec{u}(\vec{x}_w, t + \delta_t), \rho(\vec{x}_F, t) \right) + f_{\bar{i}}^{(neq)}(\vec{x}_F, t) \right\}
$$


> **[NOTE]**
> This scheme sounds easy and familar, right? Let's make a summary.
> (1) It comes from the interpolation scheme.
> (2) Employ streaming step to eliminate the unknown $f_{i}(\vec{x}_{FF}, t + \delta_t)$.
> (3) The distribution function at the wall node $\vec{x}_w$ is divieded into two parts: $f_{i}^{(eq)}$ and $f_{i}^{(neq)}$.
> (4) $f_{i}^{(eq)}$ is approximated by known macroscopic variables.
> (5) $f_{i}^{(neq)}$ is approximated by the non-equilibrium bounce-back.


### Enhanced single-node boundary condition


Here, we introduce the enhanced single-node boundary condition purposed by Marson et al.[1]. 

$$
\begin{aligned}
f_{i}(\vec{x}_F, t + \delta_t) =& \cancel{a_1 f_{\bar{i}}^{*}(\vec{x}_{FF})} + 
a_2 \underbrace{(\vec{x}_{F} + 2 f^{(eq)-}(\vec{x}_{w}))}_{\text{1}} + a_3 f_{i}^{*}(\vec{x}_{F})\\
& + a_4 \underbrace{ \tilde{f}_{i}(\vec{x}_w, t + \delta_t) }_{\text{2}} + a_5 \underbrace{\tilde{f}_{i}^{*}(\vec{x}_w)}_{\text{3}} -
\underbrace{K^{-} \frac{f_{i}^{(neq)-}(\vec{x}_F)}{\tau^{-}}}_{\text{4}}
\end{aligned}
$$

> **Term ①**: The haly-way bounce-back from Ladd [(*J. Fluid Mech.*, 271 (1994), pp. 285–309)](https://doi.org/10.1017/S0022112094001771). 
Note that: $f^{(eq)-}(\vec{x}_{w}) = w_i \rho \vec{c}_i \cdot \vec{u}_w / c_s^2$.
> **Term ②**: The term from [4].
> **Term ③**: The new term designed by Marson et al.[1].
> **Term ④**: The new (optional) term mentioned in the Sec. ⅣB of Marson et al.[1].


Marson et al.[1] use TRT collision model. So here we briefly introduce the two-relaxation-time (TRT) collision model.

$$
\begin{cases}
  f_{i}^{+} (\vec{x} + \vec{c}_i \delta_t, t + \delta_t) - f_{i}^{+} (\vec{x}, t)
  = -\frac{1}{\tau^{+}} (f_{i}^{+} (\vec{x}, t) - f_{i}^{(eq)+} (\vec{x}, t)) \\
  f_{i}^{-} (\vec{x} + \vec{c}_i \delta_t, t + \delta_t) - f_{i}^{-} (\vec{x}, t)
  = -\frac{1}{\tau^{-}} (f_{i}^{-} (\vec{x}, t) - f_{i}^{(eq)-} (\vec{x}, t)) \\
\end{cases}
$$

where

$$
\begin{aligned}
  f_{i}^{+} = \frac{1}{2} (f_i + f_{-i}) ,\quad
  f_{i}^{-} = \frac{1}{2} (f_i - f_{-i}) ,\\
  f_{i}^{(eq)+} = \frac{1}{2} (f_i^{(eq)} + f_{-i}^{(eq)}) 
  = w_i \rho (1 + \frac{(\vec{c}_i \cdot \vec{u})^2}{2 c_s^4} - \frac{|\vec{u}|^2}{2 c_s^2})
  ,\\
  f_{i}^{(eq)-} = \frac{1}{2} (f_i^{(eq)} - f_{-i}^{(eq)})
  = w_i \rho \frac{\vec{c}_i \cdot \vec{u}}{c_s^2}
  ,\\
\end{aligned}
$$

and the TRT magic number is: $\Lambda = (\tau^{+} - \frac{1}{2}) (\tau^{-} - \frac{1}{2})$.


Marson et al.[1] design the wall populations $\tilde{f}_{i}^{*}(\vec{x}_w)$ and $\tilde{f}_{i}(\vec{x}_w, t+\delta_t)$.【$\Omega$ is the collision operator.】

$$
\begin{aligned}
  \tilde{f}_{i}^{*}(\vec{x}_w, t) \overset{\text{def}}{=}& \tilde{f}_{i}^{(eq)}(\vec{x}_w, t) + (1 - \Omega) \tilde{f}_{i}^{(neq)}(\vec{x}_w, t) \\
  \tilde{f}_{i}(\vec{x}_w, t+\delta_t) \overset{\text{def}}{=}& \tilde{f}_{i}^{(eq)}(\vec{x}_w, t+\delta_t) + \tilde{f}_{i}^{(neq)}(\vec{x}_w, t+\delta_t)
\end{aligned}
$$

Consider $\tilde{f}_{i}(\vec{x}_w, t+\delta_t)$, its equilibrium part is approximated by

$$
  \tilde{f}_{i}^{(eq)} (\vec{x}_w, t+\delta_t) \approx
  f_{i}^{(eq)+} (\rho(\vec{x}_F, t), \vec{u}(\vec{x}_{\{w, F\}}, t+\delta_t)) + 
  f_{i}^{(eq)-} (\rho(\vec{x}_F, t), \vec{u}(\vec{x}_{w}, t+\delta_t))
$$

Here, $\vec{x}_{\{w, F\}}$ means use $\vec{x}_{w}$ or $\vec{x}_F$ as the position.

The non-equilibrium part is approximated by an nonequilibrium bounce-back. This is a first-order approximation of the nonequilibrium bounce-back method of Zou and He, and it leads to the following second-order accurate approximation:
$$
\tilde{f}_{i}^{(neq)}(\vec{x}_w, t+\delta_t) = \tilde{f}_{i}^{(neq)}(\vec{x}_F, t)
$$


For the $\tilde{f}_{i}^{*}(\vec{x}_w, t)$ , we have

$$
  \tilde{f}_{i}^{(eq)} (\vec{x}_w, t) \approx
  f_{i}^{(eq)+} (\rho(\vec{x}_F, t), \vec{u}(\vec{x}_{\{w, F\}}, t)) + 
  f_{i}^{(eq)-} (\rho(\vec{x}_F, t), \vec{u}(\vec{x}_{w}, t))
$$

The non-equilibrium term is also split into symmetric and anti-symmetric parts: $\tilde{f}_{i}^{(neq)}(\vec{x}_w, t) = f_{i}^{(neq)+}(\vec{x}_w, t) + f_{i}^{(neq)-}(\vec{x}_w, t)$. Marson et al.[1] gives three ways to approximate the non-equilibrium term.


(1) **Non-symmetric enhanced local interpolation** (Denoted as N): 
$$
f_{i}^{(neq)+}(\vec{x}_w, t) \approx f_{i}^{(neq)+}(\vec{x}_F, t) ,\quad 
f_{i}^{(neq)-}(\vec{x}_w, t) \approx -f_{i}^{(neq)-}(\vec{x}_F, t)
$$

(2) **Symmetric enhanced local interpolation** (Denoted as S): 
$$
f_{i}^{(neq)+}(\vec{x}_w, t) \approx f_{i}^{(neq)+}(\vec{x}_F, t) ,\quad 
f_{i}^{(neq)-}(\vec{x}_w, t) \approx f_{i}^{(neq)-}(\vec{x}_F, t)
$$


(3) **Central enhanced local interpolation**  (Denoted as C): 
$$
f_{i}^{(neq)+}(\vec{x}_w, t) \approx f_{i}^{(neq)+}(\vec{x}_F, t) ,\quad 
f_{i}^{(neq)-}(\vec{x}_w, t) = 0
$$

The following table shows the interpolation coefficients provided by Marson et al. [1]. It should be noted that $K^{-}$ is the term used to eliminate viscosity errors in steady-state conditions. In practice, $\tau^{+} \rightarrow \frac{1}{2}$ at high Reynolds numbers, making the contribution of $K^{-} \frac{f_{i}^{(neq)-}(\vec{x}_F)}{\tau^{-}}$ less significant.

<img src="Marson等-插值系数表.png" width=100%>

> **[NOTE]**: The calculation of $K^{-}$ could be found in Marson et al. [1].


## References

<style scoped>
p {
  font-size: 12pt;
}
</style>

[1] Marson, F., Thorimbert, Y., Chopard, B., Ginzburg, I. & Latt, J. **Enhanced single-node lattice Boltzmann boundary condition for fluid flows**. *Phys. Rev. E* 103, 053308 (2021).
[2] 郭照立 & 郑楚光. **格子Boltzmann方法的原理及应用**. (科学出版社, 2009).
[3] Zhao, W., Huang, J. & Yong, W.-A. **Boundary conditions for kinetic theory based models I: Lattice boltzmann models**. *Multiscale Modeling & Simulation* 17, 854–872 (2019).
[4] Tao, S., He, Q., Chen, B., Yang, X. & Huang, S. **One-point second-order curved boundary condition for lattice Boltzmann simulation of suspended particles**. *Computers & Mathematics with Applications* 76, 1593–1607 (2018). 



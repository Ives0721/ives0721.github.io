---
title: "Literature Review: 简化的格子 Boltzmann 方法"
date: 2024-05-02 10:30:33
tags: [格子Boltzmann方法]
categories:
- [格子Boltzmann方法]
---
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex/dist/katex.min.css">


> **注：该文章同步发在[知乎专栏](https://zhuanlan.zhihu.com/p/695671429?)和[CSDN](https://blog.csdn.net/weixin_43890806/article/details/138388248?spm=1001.2014.3001.5502)**  
> 这里所说的 “简化的格子 Boltzmann 方法” 其实是 ***Simplified and Highly Stable Lattice Boltzmann Method (SHSLBM)***([^c_1], [^2], [^3], [^4])。该方法的特色是仅使用平衡态分布函数、宏观密度、宏观速度进行 LBM 计算，并保证无条件稳定。  
> **下面的内容只是我对相关文献进行简单阅读后的部分整理**


## 单松弛格子 Boltzmann 方法

### 基本方程

针对流场的 LBGK 演化方程为

$$
f_{\alpha} (\boldsymbol{r} + \boldsymbol{c}_\alpha + \Delta t, t + \Delta t) - f_{\alpha} (\boldsymbol{r}, t) = 
-\frac{1}{\tau} \left[ f_{\alpha} (\boldsymbol{r}, t) - f_{\alpha}^{(eq)} \right]
\tag{1}
$$

其中 $f_{\alpha} (\boldsymbol{r}, t)$ 为 $t$ 时刻在 $\boldsymbol{r}$ 处沿着 $\boldsymbol{c}_\alpha$ 方向的密度分布函数， $\tau$ 为 BGK 碰撞项的松弛时间。 $f_{\alpha}^{(eq)}$ 为 $\boldsymbol{c}_\alpha$ 方向的平衡态分布

$$
f_{\alpha}^{(eq)} = \rho w_\alpha \left[ 1 + \frac{\boldsymbol{c}_\alpha \cdot \boldsymbol{u}}{c_s^2} + \frac{(\boldsymbol{c}_\alpha \cdot \boldsymbol{u})^2}{2 c_s^4} - \frac{\boldsymbol{u} \cdot \boldsymbol{u}}{2 c_s^2} \right] \tag{2}
$$

式 (2) 中  $\rho$ 为流体密度， $\boldsymbol{u} = (u_x, u_y, u_z)^T$ 为流体宏观速度， $w_\alpha$ 为 $\boldsymbol{c}_\alpha$ 的权重系数， $c_s$ 为格子声速。 $w_\alpha$ 、 $\boldsymbol{c}_\alpha$ 和 $c_s$ 的数值由选取的离散速度模型决定。

对于 LBM 而言，宏观量表示为

$$
\rho = \sum_{\alpha} f_{\alpha} =  \sum_{\alpha} f_{\alpha}^{(eq)}, \quad \rho\boldsymbol{u} = \sum_{\alpha} \boldsymbol{c}_\alpha f_{\alpha} =  \sum_{\alpha} \boldsymbol{c}_\alpha f_{\alpha}^{(eq)}
\tag{3}
$$

### Chapman-Enskog 展开

Boltzmann 方程和 LBGK 方程都可以通过多尺度展开还原为 Navier-Stokes 方程。这里以 Chapman-Enskog 为例：
$$
f_{\alpha} = f_{\alpha}^{(eq)} + \epsilon f_{\alpha}^{(1)} + \epsilon^2 f_{\alpha}^{(2)} + ... , \quad \frac{\partial}{\partial t} = \epsilon \frac{\partial}{\partial t_0} + \epsilon^2 \frac{\partial}{\partial t_1}, \quad \frac{\partial}{\partial \boldsymbol{r}} = \epsilon \frac{\partial}{\partial \boldsymbol{r}_0}
$$

上述展开保留至 $\epsilon^{1}$ ，即可还原为[^nabla注]
$$
\frac{\partial \rho}{\partial t} + \nabla \cdot \left(\sum_{\alpha} \boldsymbol{c}_\alpha f_{\alpha}^{(eq)} \right) = 0 
\tag{4}
$$

和
$$
\frac{\partial \rho \boldsymbol{u}}{\partial t} + \nabla \cdot \Pi = \boldsymbol{0}  \tag{5}
$$

式 (5) 中张量 $\Pi$ 表示为 
$$
\Pi_{\beta \gamma} = \sum_{\alpha} (\boldsymbol{c}_\alpha)_{\beta} (\boldsymbol{c}_\alpha)_{\gamma} \left[ f_{\alpha}^{(eq)} + \left(1 - \frac{1}{2 \tau}\right) \epsilon f_{\alpha}^{(1)} \right]
$$

由于只保留到 $\epsilon^{1}$ ，所以 $\epsilon f_{\alpha}^{(1)} \approx f_{\alpha} - f_{\alpha}^{(eq)} = f_{\alpha}^{(neq)}$ 近似为非平衡态分布函数。此时
$$
\Pi_{\beta \gamma} = \sum_{\alpha} (\boldsymbol{c}_\alpha)_{\beta} (\boldsymbol{c}_\alpha)_{\gamma} \left[ f_{\alpha}^{(eq)} + \left(1 - \frac{1}{2 \tau}\right)  f_{\alpha}^{(neq)} \right]
\tag{6}
$$

并且存在如下关系
$$
f_{\alpha}^{(neq)} \approx \epsilon f_{\alpha}^{(1)} = -\tau \left[ \frac{\partial}{\partial t} + \boldsymbol{c}_{\alpha} \cdot \nabla \right] f_{\alpha}^{(eq)} \Delta t
\tag{7}
$$

令 $\mathrm{D} = \frac{\partial}{\partial t} + \boldsymbol{c}_{\alpha} \cdot \nabla$ ，则 $f_{\alpha}^{(neq)} \approx -\tau \cdot (\Delta t) \cdot \mathrm{D}f_{\alpha}^{(eq)}$。

流体运动粘度 $\nu$ 和松弛时间 $\tau$ 的关系为
$$
\nu = \left( \tau - \frac{1}{2} \right) c_s^2 \Delta t
$$

## SHSLBM

式 (7) 说明：在保留至 $\epsilon^{1}$ 的情况下， $f_{\alpha}^{(neq)}$ 和 $f_{\alpha}^{(eq)}$ 存在联系。考虑到式 (1) 中由 $\frac{-1}{\tau}$ 松弛的部分即为 $f_{\alpha}^{(neq)}$ ，因此式 (1) 理应可以被写作完全由 $f_{\alpha}^{(eq)}$ 进行表示的形式。这便是 SHSLBM 的核心思想。

### 基本算法

SHSLBM 的单个时间步迭代分为预报和校正两部分，具体写作：

1. 预报步：计算 $(\boldsymbol{r},t)$ 的预估密度 $\rho^{*}$ 和预估速度 $\boldsymbol{u}^{*}$

$$
\rho^{*} = \sum_{\alpha} f_{\alpha}^{(eq)} (\boldsymbol{r} - \boldsymbol{c}_{\alpha} \Delta t, t -\Delta t), \quad
\rho^{*} \boldsymbol{u}^{*} = \sum_{\alpha} \boldsymbol{c}_{\alpha} f_{\alpha}^{(eq)} (\boldsymbol{r} - \boldsymbol{c}_{\alpha} \Delta t, t -\Delta t)
\tag{8}
$$

2. 校正步：根据 $\rho^{*}$ 和 $\boldsymbol{u}^{*}$ 得出校正值

$$
\begin{cases}
\rho = \rho^{*} \\
\rho \boldsymbol{u} = \rho^{*} \boldsymbol{u}^{*} - \left(1 - \frac{1}{\tau}\right) \sum\limits_{\alpha} \boldsymbol{c}_{\alpha} [f_{\alpha}^{(neq)}(\boldsymbol{r}+\frac{\Delta t}{2}\boldsymbol{c}_{\alpha},t-\frac{\Delta t}{2}) - f_{\alpha}^{(neq)}(\boldsymbol{r}-\frac{\Delta t}{2}\boldsymbol{c}_{\alpha},t-\frac{\Delta t}{2})]
\end{cases}
\tag{9-1}
$$

式 (9-1) 中非平衡态的计算通过对式 (7) 的差分实现，即：

$$
f_{\alpha}^{(neq)}(\boldsymbol{r} - \frac{\Delta t}{2} \boldsymbol{c}_{\alpha}, t - \frac{\Delta t}{2}) = -\tau \left[ f_{\alpha}^{(eq)}(\boldsymbol{r}, t) - f_{\alpha}^{(eq)}(\boldsymbol{r} - \boldsymbol{c}_{\alpha} \Delta t, t - \Delta t) \right] + O((\Delta t)^3)
$$

同理：

$$
f_{\alpha}^{(neq)}(\boldsymbol{r} + \frac{\Delta t}{2} \boldsymbol{c}_{\alpha}, t - \frac{\Delta t}{2}) = -\tau \left[ f_{\alpha}^{(eq)}(\boldsymbol{r} + \boldsymbol{c}_{\alpha} \Delta t, t) - f_{\alpha}^{(eq)}(\boldsymbol{r}, t - \Delta t) \right] + O((\Delta t)^3)
$$

在这两条式子中， $f_{\alpha}^{(eq)}(\boldsymbol{r}, t)$ 和 $f_{\alpha}^{(eq)}(\boldsymbol{r} + \boldsymbol{c}_{\alpha} \Delta t, t)$ 都是未知的。参照Shu等[^5]的做法，对上面两条公式中的这两个量使用预报步数据进行近似，即

$$
f_{\alpha}^{(eq)}(\boldsymbol{r}, t) = f_{\alpha}^{(eq)}(\rho^{*}(\boldsymbol{r}, t), \boldsymbol{u}^{*}(\boldsymbol{r}, t))
$$

式 (9-1) 中 $\rho \boldsymbol{u}$ 的计算还可以进一步化简。注意到，对于 $f_{\alpha}^{(neq)}(\boldsymbol{r} - \frac{\Delta t}{2} \boldsymbol{c}_{\alpha}, t - \frac{\Delta t}{2})$ ，有

$$
\begin{aligned}
\sum\limits_{\alpha} \boldsymbol{c}_{\alpha} f_{\alpha}^{(neq)}(\boldsymbol{r} - \frac{\Delta t}{2} \boldsymbol{c}_{\alpha}, t - \frac{\Delta t}{2}) &= -\tau \left[ \sum\limits_{\alpha} \boldsymbol{c}_{\alpha} f_{\alpha}^{(eq)}(\boldsymbol{r}, t) \right.\\
&\quad\left. -\sum\limits_{\alpha} \boldsymbol{c}_{\alpha}f_{\alpha}^{(eq)}(\boldsymbol{r} - \boldsymbol{c}_{\alpha} \Delta t, t - \Delta t)\right] + O((\Delta t)^3) \\
&= -\tau (\rho^{*} \boldsymbol{u}^{*} - \rho^{*} \boldsymbol{u}^{*}) + O((\Delta t)^3) \\
&= O((\Delta t)^3) 
\end{aligned}
$$

所以式 (9-1) 可进一步化简为式 (9-2) 

$$
\begin{cases}
\rho = \rho^{*} \\
\rho \boldsymbol{u} = \rho^{*} \boldsymbol{u}^{*} - \left(1 - \frac{1}{\tau}\right) \sum\limits_{\alpha} \boldsymbol{c}_{\alpha} f_{\alpha}^{(neq)}(\boldsymbol{r}+\frac{\Delta t}{2}\boldsymbol{c}_{\alpha},t-\frac{\Delta t}{2})
\end{cases}
\tag{9-2}
$$

### 对应的宏观流程

事实上，式 (8) 和式 (9-1), (9-2) 的预报和校正存在着对应的宏观偏微分方程。下面简要说明一下。

#### 预报步

式 (8) 对应的宏观方程为

$$
\begin{cases}
\frac{\partial \rho}{\partial t} + \nabla \cdot \left(\sum\limits_{\alpha} \boldsymbol{c}_\alpha f_{\alpha}^{(eq)} \right) = 0  \\
\frac{\partial \rho \boldsymbol{u}}{\partial t} + \nabla \cdot \left[ \sum\limits_{\alpha}(\boldsymbol{c}_\alpha)_{\beta} (\boldsymbol{c}_\alpha)_{\gamma} f_{\alpha}^{(eq)} + \frac{1}{2 \tau} \sum\limits_{\alpha}(\boldsymbol{c}_\alpha)_{\beta} (\boldsymbol{c}_\alpha)_{\gamma} f_{\alpha}^{(neq)} \right] = \boldsymbol{0}
\end{cases}
\tag{10}
$$

首先，对 $f_{\alpha}^{(eq)} (\boldsymbol{r} - \boldsymbol{c}_{\alpha} \Delta t, t -\Delta t)$ 在 $(\boldsymbol{r},t)$ 进行展开，并代入式 (7) ，得到：

$$
\begin{aligned}
f_{\alpha}^{(eq)} (\boldsymbol{r} - \boldsymbol{c}_{\alpha} \Delta t, t -\Delta t) &= f_{\alpha}^{(eq)} (\boldsymbol{r}, t) - (\Delta t) \mathrm{D} f_{\alpha}^{(eq)} (\boldsymbol{r}, t)  \\
&\quad+ \frac{(\Delta t)^2}{2} \mathrm{D}^{2} f_{\alpha}^{(eq)} (\boldsymbol{r}, t) + O((\Delta t)^3) \\
&= f_{\alpha}^{(eq)} (\boldsymbol{r}, t) - (\Delta t) \mathrm{D} f_{\alpha}^{(eq)} (\boldsymbol{r}, t)  \\
&\quad - \frac{\Delta t}{2 \tau} \mathrm{D} f_{\alpha}^{(neq)} (\boldsymbol{r}, t) + O((\Delta t)^3) \\
\end{aligned}
\tag{11}
$$

将式 (11) 代入式 (8) 的第一条公式，有：
$$
\begin{aligned}
\rho^{*} &= \sum_{\alpha}f_{\alpha}^{(eq)}(\boldsymbol{r}, t) - (\Delta t) \sum_{\alpha} \mathrm{D} f_{\alpha}^{(eq)}(\boldsymbol{r}, t) \\
&\quad - \frac{\Delta t}{2 \tau} \sum_{\alpha} \mathrm{D} f_{\alpha}^{(neq)}(\boldsymbol{r}, t)  + O((\Delta t)^3)
\end{aligned}
\tag{12}
$$

式 (12) 右侧第一项与 $\rho^{*}$ 相等，右侧第三项根据非平衡态的守恒律可知其为 0 。所以在消去式 (12) 中的 $\Delta t$ 后， $\rho^{*}$ 的计算与式 (10) 中第一条公式对应，并具有 $O((\Delta t)^2)$ 精度。


将式 (11) 代入式 (8) 的第二条公式，整理得：
$$
\begin{aligned}
\rho^{*} \boldsymbol{u}^{*} &= \sum_{\alpha} \boldsymbol{c}_{\alpha} f_{\alpha}^{(eq)}(\boldsymbol{r}, t) - (\Delta t) \sum_{\alpha} \boldsymbol{c}_{\alpha} \mathrm{D} f_{\alpha}^{(eq)}(\boldsymbol{r}, t) \\
&\quad - \frac{\Delta t}{2 \tau} \sum_{\alpha} \boldsymbol{c}_{\alpha} \mathrm{D} f_{\alpha}^{(neq)}(\boldsymbol{r}, t)  + O((\Delta t)^3) \\
&= \sum_{\alpha} \boldsymbol{c}_{\alpha} f_{\alpha}^{(eq)}(\boldsymbol{r}, t) - (\Delta t) \left[ \frac{\partial}{\partial t} \sum_{\alpha} \boldsymbol{c}_{\alpha} f_{\alpha}^{(eq)}(\boldsymbol{r}, t) + \right. \\
&\quad \left. \nabla \cdot \sum_{\alpha} (\boldsymbol{c}_{\alpha})_\beta (\boldsymbol{c}_{\alpha})_\gamma \left( f_{\alpha}^{(eq)}(\boldsymbol{r}, t) + \frac{1}{2 \tau} f_{\alpha}^{(neq)}(\boldsymbol{r}, t) \right) \right] \\
&\quad - \frac{\Delta t}{2 \tau} \sum_{\alpha} \frac{\partial}{\partial t} \boldsymbol{c}_{\alpha} f_{\alpha}^{(neq)}(\boldsymbol{r}, t)  + O((\Delta t)^3)
\end{aligned}
\tag{13}
$$

式 (13) 右侧第一项与 $\rho^{*} \boldsymbol{u}^{*}$ 相等，右侧 $\sum\limits_{\alpha} \frac{\partial}{\partial t} \boldsymbol{c}_{\alpha} f_{\alpha}^{(neq)}(\boldsymbol{r}, t)$ 项根据非平衡态的守恒律可知其为 0 。所以在消去式 (13) 中的 $\Delta t$ 后， $\rho^{*} \boldsymbol{u}^{*}$ 的计算与式 (10) 中第二条公式对应，并具有 $O((\Delta t)^2)$ 精度。

#### 校正步

式 (9) 对应的宏观方程为
$$
\begin{cases}
\frac{\partial \rho}{\partial t} = 0 \\
\frac{\partial \rho \boldsymbol{u}}{\partial t} + \nabla \cdot \left[ \left( 1 - \frac{1}{2 \tau} \right) \sum\limits_{\alpha}(\boldsymbol{c}_\alpha)_{\beta} (\boldsymbol{c}_\alpha)_{\gamma} f_{\alpha}^{(neq)} \right] = \boldsymbol{0}
\end{cases}
\tag{14}
$$

这里主要介绍一下动量方程的还原。对 $f_{\alpha}^{(neq)}(\boldsymbol{r} \pm \frac{\Delta t}{2} \boldsymbol{c}_{\alpha}, t - \frac{\Delta t}{2})$ 进行 Taylor 展开，得到


$$
\begin{aligned}
f_{\alpha}^{(neq)}(\boldsymbol{r} \pm \frac{\Delta t}{2} \boldsymbol{c}_{\alpha}, t - \frac{\Delta t}{2}) &= f_{\alpha}^{(neq)}(\boldsymbol{r}, t - \frac{\Delta t}{2}) \\
&\quad \pm \frac{\Delta t}{2} (\boldsymbol{c}_{\alpha} \cdot \nabla) f_{\alpha}^{(neq)}(\boldsymbol{r}, t - \frac{\Delta t}{2}) \\
&\quad + \frac{(\Delta t)^2}{8} (\boldsymbol{c}_{\alpha} \cdot \nabla)^2 f_{\alpha}^{(neq)}(\boldsymbol{r}, t - \frac{\Delta t}{2}) + O((\Delta t)^3)
\end{aligned}
$$


这里我们以未化简的式 (9-1) 为例。在式 (9-1) 的动量校正步中，有：

$$
\begin{aligned}
&\quad \left(1 - \frac{1}{\tau}\right) \frac{1}{\Delta t} \sum\limits_{\alpha} \boldsymbol{c}_{\alpha} [f_{\alpha}^{(neq)}(\boldsymbol{r}+\frac{\Delta t}{2}\boldsymbol{c}_{\alpha},t-\frac{\Delta t}{2}) - f_{\alpha}^{(neq)}(\boldsymbol{r}-\frac{\Delta t}{2}\boldsymbol{c}_{\alpha},t-\frac{\Delta t}{2})] \\
&= \left(1 - \frac{1}{\tau}\right) \sum\limits_{\alpha} \boldsymbol{c}_{\alpha}  (\boldsymbol{c}_{\alpha} \cdot \nabla) f_{\alpha}^{(neq)}(\boldsymbol{r}, t - \frac{\Delta t}{2}) + O((\Delta t)^2) \\
&= \nabla \cdot \left[ \left(1 - \frac{1}{\tau}\right) \sum\limits_{\alpha}(\boldsymbol{c}_\alpha)_{\beta} (\boldsymbol{c}_\alpha)_{\gamma} f_{\alpha}^{(neq)}(\boldsymbol{r}, t - \frac{\Delta t}{2}) \right]  + O((\Delta t)^2) \\
\end{aligned}
$$

此时联立式 (13) ，即可还原式 (14) 中的动量方程。
## SHSLBM 的稳定性分析

这里令 $\mathbf{Y} = (y_1, y_2, y_3, y_4)^T$。其中

$$
y_1 = \rho,\quad y_2 = \rho u_x ,\quad y_3 = \rho u_y ,\quad y_4 = \rho u_z
$$

> 将这些符号代换到式 (2) 的 $f_{\alpha}^{(eq)}$ 中，得：
> $$
> \begin{aligned}
> f_{\alpha}^{(eq)} &= y_{1 \alpha} w_{\alpha} + \frac{w_\alpha}{c_s^2} ( c_{\alpha x} y_{2 \alpha}+ c_{\alpha y} y_{3 \alpha} + c_{\alpha z} y_{4 \alpha} ) \\
> &\quad + \frac{w_\alpha}{2 y_{1 \alpha} c_s^4} ( c_{\alpha x} y_{2 \alpha}+ c_{\alpha y} y_{3 \alpha} + c_{\alpha z} y_{4 \alpha} )^2 \\
> &\quad - \frac{w_\alpha}{2 y_{1 \alpha} c_s^2} (y_{2 \alpha}^2 + y_{3 \alpha}^2 + y_{4 \alpha}^2)
> \end{aligned}
> \tag{15}
> $$
> where the subscript $\alpha$ denotes the quantities at a space level of $(\boldsymbol{r} - \boldsymbol{c}_{\alpha} \Delta t)$.[^2]


### 预报步

对于第 $n$ 步的结果 $\mathbf{Y}^{n}$，式 (8) 的预报步写作 $\mathbf{Y}^{*} = \mathbf{G}(\mathbf{Y}^{n})$。

为便于进行诺依曼稳定性分析，对上式进行求导，有
$$
\mathrm{d \mathbf{Y}^{*}} = \left[ \frac{\partial \mathbf{Y}^{*}}{\partial y_1} \right]^n \mathrm{d}y_1^n + \left[ \frac{\partial \mathbf{Y}^{*}}{\partial y_2} \right]^n \mathrm{d}y_2^n + \left[ \frac{\partial \mathbf{Y}^{*}}{\partial y_3} \right]^n \mathrm{d}y_3^n + \left[ \frac{\partial \mathbf{Y}^{*}}{\partial y_4} \right]^n \mathrm{d}y_4^n
$$

上标 $n$ 代表第 $n$ 个时间步。由于 $\mathrm{d \mathbf{Y}^{*}} = (\mathrm{d} y_1^{*}, \mathrm{d} y_2^{*}, \mathrm{d} y_3^{*}, \mathrm{d} y_4^{*})^T$ ，所以，对于预报步结果的每个分量而言有： $\mathrm{d} y_{\lambda}^{*} = \sum_{\kappa = 1}^{4}\left[ \frac{\partial y_{\lambda}^{*}}{\partial y_{\kappa}} \right]^n \mathrm{d}y_{\kappa}^n$ （$\lambda, \kappa \in [1,2,3,4]$）。


对于本文列出的方程组， $\left[ \frac{\partial y_{\lambda}^{*}}{\partial y_{\kappa}} \right]^n$ 多达 16 个。因此下文仅拿 $\left[ \frac{\partial y_{1}^{*}}{\partial y_{1}} \right]^n$ 举例进行计算。

$$
\begin{aligned}
\left[ \frac{\partial y_{1}^{*}}{\partial y_{1}} \right]^n
&= \frac{\partial}{\partial y_{1}} \sum_{\alpha} f_{\alpha}^{(eq)} \\
&= \sum_{\alpha} \left[ w_{\alpha} - \frac{w_\alpha}{2 (y_{1\alpha}^n)^2 c_s^4} ( c_{\alpha x} y_{2 \alpha}^n + c_{\alpha y} y_{3 \alpha}^n + c_{\alpha z} y_{4 \alpha}^n )^2 \right. \\
&\quad\left. + \frac{w_\alpha}{2 (y_{1\alpha}^n)^2 c_s^2} ((y_{2 \alpha}^n)^2 + (y_{3 \alpha}^n)^2 + (y_{4 \alpha}^n)^2) \right]
\end{aligned}
\tag{16}
$$

假设式 (16) 中的 $\mathrm{d} y^{*}_{i}$, $\mathrm{d} y^{n}_{i}$, $y^{n}_{i \alpha}$ 均可被展开为如下形式（$i = [1,2,3,4]$，$j=\sqrt{-1}$）

$$
\mathrm{d} y^{*}_{i} = \overline{\mathrm{d} y_{i}}^{*} \exp(j \boldsymbol{k \cdot r}) ,\quad 
\mathrm{d} y^{n}_{i} = \overline{\mathrm{d} y_{i}}^{n} \exp(j \boldsymbol{k \cdot r}) ,\quad 
y^{n}_{i \alpha} = y_{i}^{n} \exp(j \boldsymbol{k \cdot r})
$$

代入到式 (16) 和 $\mathrm{d \mathbf{Y}^{*}}$ 表达式，得

$$
\begin{aligned}
\left[ \frac{\partial y_{1}^{*}}{\partial y_{1}} \right]^n
&= \frac{\partial}{\partial y_{1}} \sum_{\alpha} f_{\alpha}^{(eq)} \\
&= \sum_{\alpha} \left[ w_{\alpha} - \frac{w_\alpha}{2 (y_{1}^n)^2 c_s^4} ( c_{\alpha x} y_{2}^n + c_{\alpha y} y_{3}^n + c_{\alpha z} y_{4}^n )^2 \right. \\
&\quad\left. + \frac{w_\alpha}{2 (y_{1}^n)^2 c_s^2} ((y_{2}^n)^2 + (y_{3}^n)^2 + (y_{4}^n)^2) \right]
\end{aligned}
\tag{17}
$$

> **格子张量的一些性质**[^2] 
> $$
> \begin{aligned}
> & \sum_{\alpha} w_{\alpha} = 1 .&\qquad \alpha \in [0,1, ..., q-1]\\
> & \sum_{\alpha} w_{\alpha} c_{\alpha \beta} =  0 .&\qquad \beta \in [x,y,z] \\
> & \sum_{\alpha} w_{\alpha} c_{\alpha \beta} c_{\alpha \gamma} = c_s^2 \delta_{\beta \gamma} .&\qquad \beta, \gamma \in [x,y,z] 
> \end{aligned}
> $$
> 注: $q$ 为离散速度数目

基于格子张量的性质，式 (17) 右侧仅剩下 $\sum\limits_{\alpha} w_{\alpha}$ ，即：

$$
\left[ \frac{\partial y_{1}^{*}}{\partial y_{1}} \right]^n = 1
$$

同理，其他项也可进行类似分析，并得出

$$
\left[ \frac{\partial y_{\lambda}^{*}}{\partial y_{\kappa}} \right]^n = \delta_{\lambda \kappa}
$$

其中 $\delta_{\lambda \kappa}$ 为狄拉克函数，且 $\lambda,\kappa \in [1,2,3,4]$ 。所以我们可以得到 $\overline{\mathrm{d} \mathbf{Y}}^{*} = \overline{\mathrm{d} \mathbf{Y}}^{n}$ 。

### 校正步

令 ${\mathrm{d} y_i}^{n+1} = \overline{\mathrm{d} y_i}^{n+1} \exp(j \boldsymbol{k \cdot r})$（$i \in [1,2,3,4]$）。由于校正步（式 (9-2) ）中并未对 $\rho$ 做修改，易得：

$$
\overline{\mathrm{d} y_1}^{n+1} = \overline{\mathrm{d} y_1}^{*} = \overline{\mathrm{d} y_1}^{n}
$$

根据式 (9-2) 的校正步，以及 $f_{\alpha}^{(neq)}(\boldsymbol{r} + \frac{\Delta t}{2} \boldsymbol{c}_{\alpha}, t - \frac{\Delta t}{2})$ 的差分表示，有 

$$
\begin{aligned}
y_2^{n+1} = y_2^{*} - (\tau - 1) y_2^{n} + (\tau - 1) \sum\limits_{\alpha} c_{\alpha x} f_{\alpha}^{(eq)}(\boldsymbol{r} + \boldsymbol{c}_{\alpha} \Delta t, t) \\
y_3^{n+1} = y_3^{*} - (\tau - 1) y_3^{n} + (\tau - 1) \sum\limits_{\alpha} c_{\alpha y} f_{\alpha}^{(eq)}(\boldsymbol{r} + \boldsymbol{c}_{\alpha} \Delta t, t) \\
y_4^{n+1} = y_4^{*} - (\tau - 1) y_4^{n} + (\tau - 1) \sum\limits_{\alpha} c_{\alpha z} f_{\alpha}^{(eq)}(\boldsymbol{r} + \boldsymbol{c}_{\alpha} \Delta t, t)
\end{aligned}
\tag{18}
$$

这里同样以 $y_2^{n+1}$ 的计算为例，记

$$
\begin{aligned}
S &= \sum\limits_{\alpha} c_{\alpha x} f_{\alpha}^{(eq)}(\boldsymbol{r} + \boldsymbol{c}_{\alpha} \Delta t, t) \\
&= \sum\limits_{\alpha} c_{\alpha x} w_{\alpha} y_{1(-\alpha)}^{*} \\
&\quad + \sum\limits_{\alpha} \frac{c_{\alpha x} w_{\alpha}}{c_s^2} (c_{\alpha x} y_{2(-\alpha)}^{*} + c_{\alpha y} y_{3(-\alpha)}^{*} + c_{\alpha z} y_{4(-\alpha)}^{*}) \\
&\quad + \sum\limits_{\alpha} \frac{c_{\alpha x} w_{\alpha}}{2 c_s^4 y_{1(-\alpha)}^{*}} (c_{\alpha x} y_{2(-\alpha)}^{*} + c_{\alpha y} y_{3(-\alpha)}^{*} + c_{\alpha z} y_{4(-\alpha)}^{*})^2 \\
&\quad - \sum\limits_{\alpha} \frac{c_{\alpha x} w_{\alpha}}{2 c_s^2 y_{1(-\alpha)}^{*}} ( (y_{2(-\alpha)}^{*})^2 + (y_{3(-\alpha)}^{*})^2 +  (y_{4(-\alpha)}^{*})^2 )^2 \\
\end{aligned}
$$

其中下标 $(-\alpha)$ 代表来自 $\boldsymbol{r} + \boldsymbol{c}_{\alpha} \Delta t$ 的值。则

$$
y_2^{n+1} = y_2^{*} + (\tau - 1) (S - y_2^{n}) \tag{19}
$$

所以其差分表示为

$$
\mathrm{d} y_2^{n+1} = \mathrm{d} y_2^{*} + (\tau - 1) (\mathrm{d} S - \mathrm{d} y_2^{n})   \tag{20}
$$

对于 $\mathrm{d} S$，有

$$
\mathrm{d} S = \left( \frac{\partial S}{\partial y_1} \right)^{*} \mathrm{d}y_1^{*} + \left( \frac{\partial S}{\partial y_2} \right)^{*} \mathrm{d}y_2^{*} + \left( \frac{\partial S}{\partial y_3} \right)^{*} \mathrm{d}y_3^{*} + \left( \frac{\partial S}{\partial y_4} \right)^{*} \mathrm{d}y_4^{*}
$$

以 $\left( \frac{\partial S}{\partial y_1} \right)^{*}$ 为例，有

$$
\begin{aligned}
\left( \frac{\partial S}{\partial y_1} \right)^{*} &=
\sum\limits_{\alpha} c_{\alpha x} w_{\alpha} \\
&\quad - \sum\limits_{\alpha} \frac{c_{\alpha x} w_{\alpha}}{2 c_s^4 (y_{1(-\alpha)}^{*})^2} (c_{\alpha x} y_{2(-\alpha)}^{*} + c_{\alpha y} y_{3(-\alpha)}^{*} + c_{\alpha z} y_{4(-\alpha)}^{*})^2  \\
&\quad + \sum\limits_{\alpha} \frac{c_{\alpha x} w_{\alpha}}{2 c_s^2 (y_{1(-\alpha)}^{*})^2} ( (y_{2(-\alpha)}^{*})^2 + (y_{3(-\alpha)}^{*})^2 +  (y_{4(-\alpha)}^{*})^2 )^2 
\end{aligned}
$$

与预报步的做法类似，将 $y_{i(-\alpha)}^{*}$ 近似为 $y_{i}^{*}$，得到：

$$
\begin{aligned}
\left( \frac{\partial S}{\partial y_1} \right)^{*} &=
\sum\limits_{\alpha} c_{\alpha x} w_{\alpha} \\
&\quad - \sum\limits_{\alpha} \frac{c_{\alpha x} w_{\alpha}}{2 c_s^4 (y_{1}^{*})^2} (c_{\alpha x} y_{2}^{*} + c_{\alpha y} y_{3}^{*} + c_{\alpha z} y_{4}^{*})^2  \\
&\quad + \sum\limits_{\alpha} \frac{c_{\alpha x} w_{\alpha}}{2 c_s^2 (y_{1}^{*})^2} ( (y_{2}^{*})^2 + (y_{3}^{*})^2 +  (y_{4}^{*})^2 )^2 
\end{aligned}
$$

此时所有 $y_{i}^{*}$ 均可被视为常数提取出来。根据格子张量的性质，可得： $\left( \frac{\partial S}{\partial y_1} \right)^{*} = 0$。

同理 
$$
\left( \frac{\partial S}{\partial y_2} \right)^{*}  = 1,\quad
\left( \frac{\partial S}{\partial y_3} \right)^{*} = \left( \frac{\partial S}{\partial y_4} \right)^{*} = 0
$$

此时 $\overline{\mathrm{d} S } = \overline{\mathrm{d} y_2}^{*}$ ，则式 (20) 可化简为：

$$
\overline{\mathrm{d} y_2}^{n+1} = \overline{\mathrm{d} y_2}^{*} = \overline{\mathrm{d} y_2}^{n}
$$

同理，对 $y_3^{n+1}$ 和 $y_4^{n+1}$ 计算可得： $\overline{\mathrm{d} y_3}^{n+1} = \overline{\mathrm{d} y_3}^{n}$ 和 $\overline{\mathrm{d} y_4}^{n+1} = \overline{\mathrm{d} y_4}^{n}$ 。

记 $\mathrm{d \mathbf{Y}^{n}} = (\mathrm{d} y_1^{n}, \mathrm{d} y_2^{n}, \mathrm{d} y_3^{n}, \mathrm{d} y_4^{n})^T$， $\mathrm{d \mathbf{Y}^{n+1}} = (\mathrm{d} y_1^{n+1}, \mathrm{d} y_2^{n+1}, \mathrm{d} y_3^{n+1}, \mathrm{d} y_4^{n+1})^T$ ，那么：

$$
\overline{\mathrm{d \mathbf{Y}}}^{n+1} = \mathbf{I} \, \overline{\mathrm{d \mathbf{Y}}}^{n}
$$

其中 $\mathbf{I}$ 为单位矩阵。显然，特征矩阵的特征值全是 1 ，代表着该格式无条件稳定。



[^c_1]: Chen, Z., Shu, C., Wang, Y., Yang, L. M., & Tan, D. (2017). A Simplified Lattice Boltzmann Method without Evolution of Distribution Function. Advances in Applied Mathematics and Mechanics, 9(1), 1–22. DOI:10.4208/aamm.OA-2016-0029

[^2]: Z. Chen, C. Shu, D. Tan; Three-dimensional simplified and unconditionally stable lattice Boltzmann method for incompressible isothermal and thermal flows. Physics of Fluids. 1 May 2017; 29 (5): 053601. DOI:10.1063/1.4983339 

[^3]: Chen Z, Shu C, Tan D, Wu C. On improvements of simplified and highly stable lattice Boltzmann method: Formulations, boundary treatment, and stability analysis. Int J Numer Meth Fluids. 2018; 87: 161–179. DOI:10.1002/fld.4485

[^4]: Chen, Z., & Shu, C. (2020). Simplified and Highly Stable Lattice Boltzmann Method. WORLD SCIENTIFIC. DOI:10.1142/12047

[^5]: C. Shu, Y. Wang, C. Teo, and J. Wu. Development of lattice Boltzmann flux solver for simulation of incompressible flows. Advances In Applied Mathematics And Mechanics. 6, 436 (2014).

[^nabla注]: 本文中所有 $\nabla$ 均表示空间偏导 $\frac{\partial}{\partial \boldsymbol{r}}$。

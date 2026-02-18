
<style>
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


# 一些基本概念

(1) 速度分布函数 $f(\boldsymbol{x}, \boldsymbol{\xi}, t)$ 定义为：在 $t$ 时刻，以位置 $\boldsymbol{x}$ 为中心的空间微元 $\mathrm{d}\boldsymbol{x}$ 内，粒子速度在 $\boldsymbol{\xi}$ 和 $\boldsymbol{\xi} + \mathrm{d}\boldsymbol{\xi}$ 之间的分子数量为 $f(\boldsymbol{x}, \boldsymbol{\xi}, t) \mathrm{d}\boldsymbol{x} \mathrm{d}\boldsymbol{\xi}$。

类似地，也可由此衍生出诸如密度分布函数、总能分布函数的概念。【当然，下文的全部讨论都是基于 **密度分布函数** 进行的】

(2) Boltzmann 方程表述为

$$
\frac{\partial f}{\partial t} + \boldsymbol{\xi} \cdot \frac{\partial f}{\partial \boldsymbol{x}} + \boldsymbol{a} \cdot \frac{\partial f}{\partial \boldsymbol{\xi}} = \Omega(f)
\tag{1}
$$

其中 $\boldsymbol{a}$ 为加速度， $\Omega(f)$ 为碰撞项。

(3) 速度分布函数的Maxwell-Boltzmann平衡态分布为：

$$
f^{eq} = \frac{n}{ (2 \pi R T)^{D/2} } \exp \left[ - \frac{(\boldsymbol{\xi} - \boldsymbol{u})^2}{2 R T} \right]
$$

其中:
* $n = \int f \mathrm{d}\boldsymbol{\xi}$ 为分子数密度。
* $\boldsymbol{u} = \frac{1}{n} \int \boldsymbol{\xi} f \mathrm{d}\boldsymbol{\xi}$ 为流场速度。
* $R$ 为气体常数； $T$ 为流场温度 (等温系统的恒温)； $D$ 为问题的维度数量。

而密度分布函数的Maxwell-Boltzmann平衡态分布为：

$$
f^{eq} = \frac{\rho}{ (2 \pi R T)^{D/2} } \exp \left[ - \frac{(\boldsymbol{\xi} - \boldsymbol{u})^2}{2 R T} \right]
$$

其中: $\rho$ 为流场密度。

(4) BGK 碰撞项表示为

$$
\Omega(f) = -\frac{1}{\tau} \left[ f - f^{eq} \right]
$$

(5) 碰撞不变量 $\mathcal{\psi} = [1, \boldsymbol{\xi}, |\boldsymbol{\xi}|^2]$

# Part A: 不含外力项的 DUGKS 

若暂时不考虑外力作用的部分（即 $\boldsymbol{a} \cdot \frac{\partial f}{\partial \boldsymbol{\xi}}$），则 Boltzmann 方程表述为

$$
\frac{\partial f}{\partial t} + \boldsymbol{\xi} \cdot \frac{\partial f}{\partial \boldsymbol{x}} = \Omega(f)
\tag{2}
$$


这里，我们关注 $t_n$ 时刻某个网格点 $\boldsymbol{x}_i$ 上特定离散速度方向 $\boldsymbol{\xi}_k$ 的分布函数 $f(\boldsymbol{x}_i, \boldsymbol{\xi}_k, t_n)$，下面简写为 $f_{i,k}^{n}$。并定义以网格点 $\boldsymbol{x}_i$ 为中心的控制体的体积为 $V_i$。下文中的所有 $\Omega$ 都使用 **BGK 碰撞项**。

此时，对于从 $t_n$ 时刻到 $t_{n+1} = t_n + \Delta t$ 时刻的过程，我们就可以对这个控制体内的 Boltzmann 方程进行差分，得到：
$$
\underbrace{f_{i,k}^{n+1} - f_{i,k}^{n}}_{(a)} + 
\underbrace{
\frac{\Delta t}{V_i} \sum_{j \in N(i)} \left[
    (\boldsymbol{\xi}_k \cdot \boldsymbol{n}_{(ij)}) f_{(ij),k}^{n+1/2} S_{(ij)}
\right]    
}_{(b)} \\ = 
\underbrace{\frac{\Delta t}{2} \left[ \Omega(f_{i,k}^{n+1}) + \Omega(f_{i,k}^{n}) \right]}_{(c)} 
\tag{3}
$$

<div class="note-block">

【1】对于 (a) 项，我们很容易能看出来是对 $\frac{\partial f}{\partial t}$ 的差分。

【2】对于 (b) 项，其含义是计算控制体的界面通量。毕竟通量也是造成控制体内部分布函数变化的一个原因。

记 $N(i)$ 是与控制体 $i$ 相邻的所有其他网格， $j \in N(i)$ 是其中一个与之相邻的控制体。则：

* $\boldsymbol{n}_{(ij)} = (\boldsymbol{x}_j - \boldsymbol{x}_i) / \left \| \boldsymbol{x}_j - \boldsymbol{x}_i \right \|_2$ 是控制体 $i$ 与控制体 $j$ 的交界面的单位法向矢量, 方向由控制体 $i$ 指向控制体 $j$。
* $\boldsymbol{x}_{ij}$ 是控制体 $i$ 与控制体 $j$ 的交界面与 $\boldsymbol{x}_j - \boldsymbol{x}_i$ 的交点。也就是说存在着： $(\boldsymbol{x}_{ij} - \boldsymbol{x}_i) \parallel \boldsymbol{n}_{ij}$ 且 $(\boldsymbol{x}_{ij} - \boldsymbol{x}_i) \cdot \boldsymbol{n}_{ij} > 0$。
* $f_{(ij),k}^{n+1/2} = f(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n+1/2})$, 其中 $t_{n+1/2} = t_{n} + \Delta t / 2$。
* $S_{(ij)}$ 是控制体 $i$ 与控制体 $j$ 的交界面的面积。

【3】对于 (c) 项，选择用梯形积分格式描述式(2)等号右侧碰撞项 $\Omega$ 的作用。

</div>

在式(3)中，我们明显可见仅有的未知量是 $f_{(ij),k}^{n+1/2} = f(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n+1/2})$ 和 $f_{i,k}^{n+1} = f(\boldsymbol{x}_{i}, \boldsymbol{\xi}_k, t_{n+1})$，而后者是我们待求的量。

## 引入辅助分布函数

定义

$$
\begin{aligned}
    \tilde{f} &= f - \frac{\Delta t}{2} \Omega(f) \\
    \tilde{f}^{+} &= f + \frac{\Delta t}{2} \Omega(f)
\end{aligned}
$$

则式(3)可以这样改写：

$$
\begin{aligned}
[ f_{i,k}^{n+1} &- \frac{\Delta t}{2} \Omega(f_{i,k}^{n+1}) ] - [f_{i,k}^{n} + \frac{\Delta t}{2} \Omega(f_{i,k}^{n})]\\
&= - \frac{\Delta t}{V_i} \sum_{j \in N(i)} \left[
    (\boldsymbol{\xi}_k \cdot \boldsymbol{n}_{(ij)}) f_{(ij),k}^{n+1/2} S_{(ij)}
\right] 
\end{aligned}
$$

把中括号内的表达式用定义的辅助函数代换掉，即：

$$
\tilde{f}_{i,k}^{n+1} - \tilde{f}_{i,k}^{+,n} =- 
\frac{\Delta t}{V_i} \sum_{j \in N(i)} \left[
    (\boldsymbol{\xi}_k \cdot \boldsymbol{n}_{(ij)}) f_{(ij),k}^{n+1/2} S_{(ij)}
\right]
\tag{4}
$$

<div class="note-block">

当使用BGK碰撞项时，两个辅助函数分别展开为：

$$
\begin{aligned}
\tilde{f} &= f - \frac{\Delta t}{2} \Omega(f) \\
&= f - \frac{\Delta t}{2} \left\{ -\frac{1}{\tau} \left[ f - f^{eq} \right] \right\} \\
&= \frac{2 \tau + \Delta t}{2 \tau} f - \frac{\Delta t}{2 \tau} f^{eq}
\end{aligned}
$$

和

$$
\begin{aligned}
\tilde{f}^{+} &= f + \frac{\Delta t}{2} \Omega(f) \\
&= f + \frac{\Delta t}{2} \left\{ -\frac{1}{\tau} \left[ f - f^{eq} \right] \right\} \\
&= [1 - \frac{\Delta t}{2 \tau}] f + \frac{\Delta t}{2 \tau} f^{eq} \\
&= \frac{2 \tau - \Delta t}{2 \tau + \Delta t} \tilde{f} + \frac{2 \Delta t}{2 \tau + \Delta t} f^{eq}
\end{aligned}
$$

</div>

这种构造方式有一个特点：**不影响守恒量的计算**。由于 BGK 碰撞项满足如下性质

$$
\int \Omega \mathrm{d}\boldsymbol{\xi} = 0 ,\quad
\int \Omega \boldsymbol{\xi} \mathrm{d}\boldsymbol{\xi} = 0 ,\quad
\int \Omega |\boldsymbol{\xi}|^2  \mathrm{d}\boldsymbol{\xi} = 0,
$$

在 $\tilde{f}$ (或 $\tilde{f}^{+}$) 的速度矩积分中，和碰撞项有关的部分总能被单独分离出来 [积分的加减法] ，所以也就有：

$$
\rho = \int \tilde{f} \mathrm{d}\boldsymbol{\xi} ,\quad
\boldsymbol{u} = \frac{1}{\rho} \int \boldsymbol{\xi} \tilde{f} \mathrm{d}\boldsymbol{\xi}.
$$


## 通量的简化计算

对于 $f_{(ij),k}^{n+1/2} = f(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n+1/2})$ ，我们考虑 $t \in [t_{n}, t_{n+1/2}]$ 这段总时长为 $h = \Delta t / 2$ 的时间段。那么我们就可以认为 $f_{(ij),k}^{n+1/2}$ 它是在 $t_n$ 时刻，从点 $\boldsymbol{x}'_{ij} = \boldsymbol{x}_{ij} - h \boldsymbol{\xi}_k$ 沿着 $\boldsymbol{\xi}_k$ 方向迁移过来的。将式(1)从 $t_{n}$ 到 $t_{n+1/2}$ 沿特征线 $\boldsymbol{\xi}_k$ 积分，得：

$$
\begin{aligned}
f(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n+1/2}) &- f(\boldsymbol{x}'_{ij}, \boldsymbol{\xi}_k, t_{n}) \\
&= \frac{h}{2} [\Omega(f(\boldsymbol{x}'_{ij}, \boldsymbol{\xi}_k, t_{n})) \\
&\quad + \Omega(f(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n+1/2}))]
\end{aligned}
$$

再次引入下面的辅助函数协助化简：

$$
\begin{aligned}
\bar{f}     = f - \frac{h}{2} \Omega &= \frac{2 \tau + h}{2 \tau} f - \frac{h}{2 \tau} f^{eq}, \\
\bar{f}^{+} = f + \frac{h}{2} \Omega &= \frac{2 \tau - h}{2 \tau + h} f + \frac{2h}{2 \tau + h} f^{eq}.
\end{aligned}
$$

则

$$
\bar{f}(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n+1/2}) = \bar{f}^{+}(\boldsymbol{x}'_{ij}, \boldsymbol{\xi}_k, t_{n}).
\tag{5}
$$

由于式(5)中的 $\bar{f}^{+}(\boldsymbol{x}'_{ij}, \boldsymbol{\xi}_k, t_{n})$ 是 $t_{n}$ 时刻的值，所以可以使用该时刻上网格点的信息插值得到。也就是说，

$$
\bar{f}^{+}(\boldsymbol{x}'_{ij}, \boldsymbol{\xi}_k, t_{n}) = \bar{f}^{+}(\boldsymbol{x}_{c}, \boldsymbol{\xi}_k, t_{n}) + (\boldsymbol{x}'_{ij} - \boldsymbol{x}_{c}) \cdot \nabla \bar{f}^{+}(\boldsymbol{x}_{c}, \boldsymbol{\xi}_k, t_{n})
\tag{6}
$$ 

上式中的 $\boldsymbol{x}_{c}$ 是进行插值的中心点，可以根据迎风特性选取 $\boldsymbol{x}_{i}$ (或 $\boldsymbol{x}_{j}$)，或者采用中心差分从而选择 $\boldsymbol{x}_{ij}$ 。 下面仅讲解使用 $\boldsymbol{x}_{ij}$ 时的做法。

当 $\boldsymbol{x}_{c}=\boldsymbol{x}_{ij}$ 时，也就意味着 $\nabla \bar{f}^{+}(\boldsymbol{x}_{c}, \boldsymbol{\xi}_k, t_{n})$ 的计算也要跟着用中心差分。所以：

$$
\begin{aligned}
    \bar{f}^{+}(\boldsymbol{x}'_{ij}, \boldsymbol{\xi}_k, t_{n}) &= \bar{f}^{+}(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n}) + (\boldsymbol{x}'_{ij} - \boldsymbol{x}_{ij}) \cdot \nabla \bar{f}^{+}(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n}), \\
    \nabla \bar{f}^{+}(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n}) &=
    \frac{\bar{f}^{+}(\boldsymbol{x}_{j}, \boldsymbol{\xi}_k, t_{n}) - \bar{f}^{+}(\boldsymbol{x}_{i}, \boldsymbol{\xi}_k, t_{n})}{\boldsymbol{x}_{j} - \boldsymbol{x}_{i}} .
\end{aligned}
\tag{7}
$$

将式(7)回代到式(5)中，即可先把 $\bar{f}(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n+1/2})$ 求解出来。而式(4)所需的 $f_{(ij),k}^{n+1/2} = f(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n+1/2})$ 则表示为：

$$
f_{(ij),k}^{n+1/2} = \frac{2 \tau}{2 \tau + h}
\left[
    \bar{f}(\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n+1/2}) + \frac{h}{2 \tau} f^{eq} (\boldsymbol{x}_{ij}, \boldsymbol{\xi}_k, t_{n+1/2})
\right] 
\tag{8}
$$

我们前面已经提到了: 碰撞项的速度矩积分为0 ，所以：计算 $f^{eq}$ 所需要的宏观量是从 $\bar{f}$ 得到的。对 $\bar{f}$ 的速度积分即可得到 $t_{n+1/2}$ 时刻在 $\boldsymbol{x}_{ij}$ 位置的宏观量。

## 分子速度空间的离散化 - 以低马赫数情况为例

### 基本原理

> 这一小节的内容看着很熟悉对吗？格子Boltzmann方法 (LBM) 在处理低马赫数场景时也是这么做的！

<div class="note-block">

前面两个小结已经介绍完了在某个特定分子速度 $\boldsymbol{\xi}_{k}$ 下的离散，这一节则是讨论构建离散速度集及其对应的权重集。

它的最终目标是让分布函数的某些速度矩积分能够被近似，从而便于数值计算。比如：

$$
\rho = \int f \mathrm{d}\boldsymbol{\xi} = \sum f(\boldsymbol{\xi}_{\alpha}), \quad
\rho \boldsymbol{u} = \int \boldsymbol{\xi} f \mathrm{d}\boldsymbol{\xi} = \sum \boldsymbol{\xi}_{\alpha} (w_i f(\boldsymbol{\xi}_{\alpha}))
$$

这里的 $w_\alpha$ 就是 $\boldsymbol{\xi}_{\alpha}$ 的权重系数。上面的式子只写了 $\rho$ 和 $\rho \boldsymbol{u}$ 的原因是：BGK 碰撞项不一定完全守恒，简单来说就是“用分布函数算的速度矩”不一定就全都会等于“用平衡态算的速度矩”。

$$
W \equiv \sum_{\alpha} \mathcal{\psi}(\boldsymbol{\xi}_{\alpha}) \cdot (w_i f(\boldsymbol{\xi}_{\alpha}))
\neq 
\sum_{\alpha} \mathcal{\psi}(\boldsymbol{\xi}_{\alpha}) \cdot (w_i f^{eq}(\boldsymbol{\xi}_{\alpha}, W))
$$

上面的 $\mathcal{\psi}$ 既可以指代碰撞不变量，也可通指由 $\boldsymbol{\xi}$ 构成的所有多项式。这里我们简单介绍如何使用 Gauss-Hermite 积分做到让最简单的零阶矩和一阶矩保持一致的。

</div>

<!-- 假设我们有一个 $D$ 维空间下的离散速度集合 $\zeta = \{\boldsymbol{\xi}_0, \boldsymbol{\xi}_1, ..., \boldsymbol{\xi}_{b-1}\}$ ，对应的权重系数的集合为 $\mathbb{W} = \{w_0, w_1, ..., w_{b-1}\}$ 。 -->

当马赫数远小于1时，我们可以对平衡态做如下离散：

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
\tag{9}
$$

此时我们的目标变成了：

$$
\int \psi(\xi) f^{eq} \mathrm{d}\boldsymbol{\xi} = 
\sum_{\boldsymbol{\alpha}} \mathcal{\psi}(\boldsymbol{\xi}_{\alpha}) \cdot (w_{\boldsymbol{\alpha}} f^{eq}(\boldsymbol{\xi}_{\boldsymbol{\alpha}})), \quad
\psi(\boldsymbol{\xi}) = 1, \boldsymbol{\xi}.
\tag{10}
$$

根据上面的离散平衡态（式(9)），我们可以将权函数定义为 : $\displaystyle\exp\left( - \frac{|\boldsymbol{\xi}|^2}{2 R T} \right) (2 \pi R T)^{-D/2}$ 。假设所选求积公式的积分点及对应权重分别为 $\boldsymbol{\xi}_{\boldsymbol{\alpha}}$ 和 $\boldsymbol{W}_{\boldsymbol{\alpha}}$， 则式(10)的权重 $w_{\bold{i}}$ 表示为

$$
w_{\boldsymbol{\alpha}} = \boldsymbol{W}_{\boldsymbol{\alpha}} (2 \pi R T)^{D/2} \exp\left( \frac{|\boldsymbol{\xi}_{\boldsymbol{\alpha}}|^2}{2 R T} \right)
$$

这里我们可以用 Gauss-Hermite 积分来解决数值积分问题。假设 $\xi_{\mathrm{\alpha}}$ 和 $\mathrm{W}_{\mathrm{\alpha}}$ 是一维Gauss-Hermite积分的积分点及对应权重，那么对于我们要计算的 $D$ 维问题，需要的积分点和积分权重变成了：

$$
\boldsymbol{\xi}_{\boldsymbol{\alpha}} = (\xi_{\mathrm{\alpha}_1}, \xi_{\mathrm{\alpha}_2}, ..., \xi_{\mathrm{\alpha}_D}), \quad
\boldsymbol{W}_{\boldsymbol{\alpha}} = \mathrm{W}_{\mathrm{\alpha}_1} \mathrm{W}_{\mathrm{\alpha}_2} ... \mathrm{W}_{\mathrm{\alpha}_D}
$$

这里这个加粗的 $\boldsymbol{\alpha}=\mathrm{\alpha}_1 \mathrm{\alpha}_2 ... \mathrm{\alpha}_D$ 。

在计算中，权重系数 $w_{\boldsymbol{\alpha}}$ 可以被吸收到离散分布函数中，即:

$$
f_{\boldsymbol{\alpha}} = w_{\boldsymbol{\alpha}} f(\boldsymbol{\xi}_{\boldsymbol{\alpha}}), \quad
f^{eq}_{\boldsymbol{\alpha}} = w_{\boldsymbol{\alpha}} f^{eq}(\boldsymbol{\xi}_{\boldsymbol{\alpha}})
$$

此时，离散平衡态 $f^{eq}_{\boldsymbol{\alpha}}$ 就变成了

$$
f^{eq}_{\boldsymbol{\alpha}} (\boldsymbol{\xi}) = \rho W_{\boldsymbol{\alpha}} \times 
\left[
    1 + \frac{\boldsymbol{\xi} \cdot \boldsymbol{u}}{R T} + 
    \frac{(\boldsymbol{\xi} \cdot \boldsymbol{u})^2}{2 (R T)^2} - 
    \frac{|\boldsymbol{u}|^2}{2 R T}
\right]
$$

## 算法的全貌 [不包含外力]

【TODO】


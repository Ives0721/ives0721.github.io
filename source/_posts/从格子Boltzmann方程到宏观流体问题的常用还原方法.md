---
title: 从格子Boltzmann方程到宏观流体问题的常用还原方法
date: 2023-04-08 17:51
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

![](https://img.shields.io/badge/格子Boltzmann方法-LBE2NS-green)

{% note default %}
这是我之前的一次期末大报告作业，里面的内容也算是从其他文章里搬运的，在此重新誊抄出来。内容如有错误，欢迎各位指正。

[注] 1. 下文提到的 $x$ 、 $u$ 和 $\xi$ 无论有没有加粗均为向量。  
2. 本文已经同步发表在本人的[知乎](https://zhuanlan.zhihu.com/p/476380996)和[CSDN](https://blog.csdn.net/weixin_43890806/article/details/130032260?spm=1001.2014.3001.5501)上，可在其他网站获得同样的阅读效果。
{% endnote %}

## 一、 单松弛格子Boltzmann方程
### 1.1		Boltzmann方程
格子Boltzmann方法基于Boltzmann方程（即(1)式），该方程引入速度分布函数，描述了非热力学平衡状态的热力学系统统计行为：

$$
\frac{\partial f}{\partial t} + \mathbf{\xi} \cdot \frac{\partial f}{\partial \boldsymbol{x}} + \frac{F}{m} \cdot \frac{\partial f}{\partial \mathbf{\xi}}=\Omega(f)  \tag{1}
$$

(1)式中 $\Omega(f)$ 为碰撞项， $\mathbf{\xi}=\frac{\partial x}{\partial t}$ 。速度分布函数的定义是：在时刻$t$，以$x$为中心的空间微元 $\mathrm{d}\mathit{x}$ 内，速度在 $\xi$ 和 $(\xi+\mathrm{d}\it{\xi})$ 之间的分子的数量为 $f \mathrm{d}\it{x} \mathrm{d}\it{\xi}$ 。记分子质量为 $m$ ，流体宏观速度为$u$，单位体积内能为$e$，分子数密度为$n$，则：

$$\rho = mn = m \int{f} \mathrm{d}\xi$$

$$\rho u = m \int{\xi f} \mathrm{d}\xi$$

$$\rho E = \rho e + \frac{1}{2} \rho u^2 = m \int{ f \frac{\xi^2}{2} } \mathrm{d} \xi$$
 
碰撞项  描述了分子碰撞对整个系统运动产生的影响，它具有两个特点：①满足质量、动量和能量守恒；②能够反映系统趋于平衡态的趋势。Bhatnagur、Gross和Krook三人在此基础上提出了最为简单且著名的BGK模型，即(2)式。

$$\Omega(f) = -\frac{1}{\tau_c} \left ( f - f^{(eq)} \right )  \tag{2}$$ 

其中$\tau_c$为松弛时间， $f^{(eq)}$ 为平衡态分布函数。 可根据Boltzmann H定理计算得：

$$f^{(eq)} = \frac{n}{(2 \pi RT)^{3/2}}   \mathrm{exp}\left( -\frac{(\mathbf{\xi-u})^2}{2RT} \right)  \tag{3}$$

### 1.2		格子Boltzmann方程
为了将Boltzmann方程应用于实际的数值模拟，需要将流体视为大量的离散粒子。离散粒子于流体分子不同，它们驻留在网格上，并按一定的规则在网格上发生碰撞和迁移。因此，1.2节将 $mf$ 记为 $f$，此时 $f$ 为密度分布函数。对密度分布函数 $f$ ，不考虑 $F$ 的Boltzmann-BGK方程为式。

$$\frac{\partial f}{\partial t} + \mathbf{\xi} \cdot \frac{\partial f}{\partial \mathit{x}} = - \frac{1}{\tau_c} \left[ f - f^{(eq)} \right] \tag{4}$$

其宏观量按如下方式确定

$$\rho = \int { f } \mathrm{d} \xi ，\rho u = \int { \xi f } \mathrm{d} \xi ，\rho E = \int { f \frac{\xi^2}{2} } \mathrm{d} \xi$$ 

对 $f$ 和其对应的 $f^{(eq)}$ 可进行如下离散：$f_i = w_i f$ ， $f_i^{(eq)} = w_i f^{(eq)}$。其中下标 $i$ 表示特征方向 $\boldsymbol{c}_i$ ， $f_i$ 和$f_i^{(eq)}$均沿着$\boldsymbol{c}_i$方向运动。据此，沿着特征方向离散(4)式得到的差分格式为 

$$f_i (x+ξ_i Δt,t+Δt)-f_i (x,t)=-\frac{1}{τ} \left[f_i (x,t)-f_i^{(eq)} (x,t) \right] \tag{5}$$

其中 $\tau = \frac{\tau_c}{\Delta t}$ 为无量纲松弛时间。由于对应的 $f^{(eq)}$ 为：

$$f^{(eq)} = \rho \cdot \mathrm{exp} ( -  \frac{(\xi - \it{u})^2}{2\it{RT}} ) \tag{6-1}$$

所以将 $f^{(eq)}$ 展开至2阶并离散，可得：

$$f^{(eq)}_i = w_i \rho \left[ 1 + \frac{\xi_i \cdot u}{c_s^2} + \frac{(\xi_i \cdot u)^2}{2 c_s^4} -  \frac{u \cdot u}{2 c_s^2} \right] \tag{6-2}$$ 

其中 $c_s = \sqrt{RT}$ 在LBM模拟中一般被称作格子声速。

## 二、Chapman-Enskog展开分析
将基本碰撞不变量 $φ_i=1,ξ,\frac{ξ^2}{2}$ 乘以(1)式两端，并对 $\xi$ 积分有：

$$\frac{∂ρ}{\partial t}+∇⋅(ρ \mathbf{u})=0 \tag{7-1}$$

$$\frac{∂(ρu)}{\partial t}+∇⋅(ρ \mathbf{uu}+P) = ρ \mathit{\bf{a}} \tag{7-2}$$

$$\frac{∂(\rho E)}{\partial t} + ∇⋅(\rho \mathbf{u}E + P⋅\mathbf{u} + \mathbf{Q})=\rho \mathbf{a⋅u} \tag{7-3}$$

Chapman-Enskog展开（下文简称C-E展开）是处理Boltzmann方程的一种分析方法，该方法引入小量$\varepsilon$作为展开因子（一般与Knudsen数同阶），将BGK方程在不同的尺度上展开，即：

$$f = \sum_{n=0}^{\infty} \varepsilon^n f^{(n)}  \tag{8}$$

并将分布函数、导数、物理量等都按照$\varepsilon$的不同阶次展开。该方法对 $f^{(n)}$ 进行了限制，使其满足：$\int f^{(n)} φ_i \mathrm{d}ξ=0$， $φ_i (i=1,..,5)$ 为基本碰撞不变量。将 $f$ 视为 $\rho$ 、$\mathbf{u}$ 和 $T$ 及其各阶导数的泛函，而时间偏导项则展开为

$$\frac{\partial }{\partial t} = \sum_{n=0}^{\infty }  \frac{\partial _n}{\partial t}  \tag{9}$$

将(8)式和(9)式代入到(7)式中。得 $O({\varepsilon}^0)$ 方程组为：

$$\frac{∂_0 ρ}{\partial t}+∇⋅(ρ \mathbf{u})=0 \tag{10-1}$$


$$\frac{∂_0 (ρu)}{\partial t}+∇⋅(ρ \mathbf{uu}+P^{(0)}) = ρ \mathbf{a}  \tag{10-2}$$


$$\frac{∂_0(\rho E)}{\partial t} + ∇⋅(\rho \mathbf{u}E + P^{(0)} ⋅ \mathbf{u} + \mathbf{Q}^{(0)})=\rho \mathbf{a⋅u}  \tag{10-3}$$

以及 $O({\varepsilon}^n)$ 方程组为 $(n>0)$ ：

$$\frac{∂_n ρ}{\partial t}=0  \tag{11-1}$$


$$\frac{∂_n (ρu)}{\partial t}+∇ ⋅ P^{(n)} = 0  \tag{11-2}$$


$$\frac{∂_n (\rho E)}{\partial t} + ∇⋅( P^{(n)} ⋅ \mathbf{u} + \mathbf{Q}^{(n)} ) = 0  \tag{11-3}$$

其中  
$$C=\boldsymbol{ξ-u}, P^{(n)}=m∫CCf^{(n)} \mathrm{d}ξ , Q^{(n)}=m∫C \frac{C^2}{2} f^{(n)} \mathrm{d}ξ$$
 
只需取精确到$O({\varepsilon}^1)$的宏观方程组进行分析，即可还原出Naiver-Stokes方程。且具有如下关系：

$$P=P^{(0)} + \varepsilon P^{(1)} = p \mathbf{I} - 2 \varepsilon \mu (\mathbf{S} - \frac{\mathrm{Tr}(S)}{3} \mathbf{I})$$
 
$$Q = Q^{(0)} + \varepsilon Q^{(1)} = - \varepsilon \kappa \nabla{T}$$
 
其中：$\kappa$为热传导系数， $\mu$ 为粘性系数。

离散的LBGK方程也可采用类似的方法展开，还原出宏观流体运动方程（见附录A）。展开方式如下：

$$f_i = f_i^{(0)} + \varepsilon f_i^{(1)} + \varepsilon^2 f_i^{(2)} + ...$$
 

$$\frac{\partial}{\partial x} = \varepsilon \frac{\partial}{\partial x_1} , x_1 = \varepsilon x$$
 
略去$O({\varepsilon}^3)$项，并对$O({\varepsilon}^0)$项、$O({\varepsilon}^1)$项和$O({\varepsilon}^2)$项求0阶和1阶速度矩，即可得到不可压缩流体的Navier-Stokes方程。由于在分析过程中略去了$O({\varepsilon}^3)$项和速度的3阶项，导致结果缺失了部分物理信息。

## 三、 Grad-13矩方法及其发展
### 3.1		Grad的矩分析方法
Grad -13矩是一种由数学家Harold Grad提出的分析方法，其基于Hermite多项式进行分析。对于一个长度为 $N$ 的向量 $\boldsymbol{x}$ ，若定义权函数

$$\omega( \boldsymbol{x}) = \frac{1}{(2 \pi)^{N/2}} \mathrm{exp} \left(- \frac{ \boldsymbol{x} \cdot  \boldsymbol{x}}{2} \right)      \tag{12}$$

则其对应的$n$阶Hermite多项式写作

$$\mathcal{H} ^{(n)} ( \boldsymbol{x}) = \frac{(-1)^n}{\omega ( \boldsymbol{x})} \nabla^{n} \omega ( \boldsymbol{x}) \tag{13}$$
 由Hermite多项式的正交性，可证明得到：

$$∫_{-∞}^{+∞} \omega ( \boldsymbol{x}) H_i^{(m)}( \boldsymbol{x}) H_j^{(n)}( \boldsymbol{x}) \mathrm{d} \boldsymbol{x} = \delta_{mn} \delta_{ij} \tag{14}$$ 

据此，可将任意一个函数 $f(x)$ 使用Hermite多项式进行展开[^GRAD_1949]。即：

$$f( \boldsymbol{x})=ω^{\frac{1}{2}} \sum_{n=0}^∞ \sum_{i=1}^N a_i^{(n)} H_i^{(n)}( \boldsymbol{x})  \tag{15}$$
 所以

$$\int_{-\infty}^{-\infty} \omega^{1/2} f(x) \mathcal{H}_j^{(m)} (x) \mathrm{d} \it{x} = \sum_{i=\mathrm{1}}^{N} a_i^{\mathrm{(}m\mathrm{)}} \delta_{ij}^m = m\mathrm{!} a_i^{m}   \mathrm{\tag{16}}$$

其中 $a_i^{(m)}$ 可借由Hermite多项式的正交性进行求解。 $\delta_{ij}^{m}$ 是一个记号，表示 $m$个 $\delta_{ij}$的乘积之和， $i$ 和 $j$ 分别来自下标集 $(i_1 , i_2 , ... , i_m)$ 和 $(j_1 , j_2 , ... , j_m)$ 。
在此数学基础之上，Harold Grad将分布函数 $f$ 展开为[^GRAD_1949_LBM]：

$$f(x, \xi, t) = \omega \sum_{n=0}^{\infty} \frac{1}{n!} a^{(n)} H^{(n)}(\xi) \tag{17}$$

根据式(16)， $a^{(n)}$ 的表达式为：

$$a^{(n)} = \int_{-\infty}^{+\infty} f(x, \xi, t) H^{(n)}(\xi) d\xi \tag{18}$$

可得 $a^{(n)}$ 的前4项为：

$$a^{(0)} = \int_{-\infty}^{+\infty} f(x, \xi, t) d\xi= \rho   \tag{19-1}$$


$$a^{(1)} =  \int_{-\infty}^{+\infty} f(x, \xi, t) \xi d\xi = \rho u  \tag{19-2}$$


$$a^{(2)} =  \int_{-\infty}^{+\infty} f(x, \xi, t) (\xi^2 - \delta) d\xi = \rho (\mathbf{u}^2 - \mathbf{\delta}) +  \mathbf{P}  \tag{19-3}$$


$$a^{(3)} =  \int_{-\infty}^{+\infty} f(x, \xi, t) (\xi^3 - \xi \delta) d\xi = (D-1) \rho \bf{u}^3 + \bf{u} \it{a}^{\mathrm{(2)}} \mathrm{+}  \mathbf{Q}   \tag{19-4}$$

注意到从 $a^{(0)}$ 到 $a^{(3)}$ 都描述了流场的宏观量，所以通过上述表达式的值可以解出流场宏观量。与C-E分析方法不同的是分析过程中没有舍去高阶量。因此这种展开方法得到的式子合理地引入了更高阶的矩，物理意义更加明确。
### 3.2		Gauss-Hermite积分的引入和离散化BGK方程
Grad的思路开创了描述的新方式，而单肖文、袁学锋、陈沪东三位学者在此基础上进一步发展，于2006年在 Journal of Fluid Mechanics 上提出了将原思路拓展到离散的Boltzmann-BGK方程中的方法[^SHAN等_2006]，简化了求解过程。
首先，基于式(16)将 $f$ 截断至前 $N$ 阶，得：

$$f \approx f^N = \omega \sum_{n=0}^{N} \frac{1}{n!} a^{(n)} H^{(n)}(\xi)  \tag{20}$$

$f^N$ 和 $f$ 有着完全相同的前 $N$ 阶矩。通过截断操作将BGK方程对宏观流动方程的近似保持在动力学层面而不是热力学层面。此时$f^N$可表示为：

$$f^{N}(x, \xi , t) H^{(n)}(\xi) = \omega(\xi) p(x, \xi, t)  \tag{21}$$

其中 $p(x, \xi, t)$ 为阶数不超过 $2N$ 的多项式。使用Gauss-Hermite积分可将 $a^{(n)}$ 精确地展开为 $p(x, \xi, t)$ 的加权和。

$$a^{(n)} =\int \omega(\xi) p(x, \xi, t) d \xi = \sum_{a=1}^{d} w_a p(x, \xi_a, t) = \sum_{a=1}^{d} \frac{w_a}{\omega(\xi_a)} f^{N}(x, \xi_a, t) H^{(n)}(\xi_a)   \tag{22}$$

式(22)中 $w_a$ 和 $ \boldsymbol{\xi}_a$ 是Gauss-Hermite积分的权重和特征向量（ $a=1,2,...,d$ ）。所以使用 $f^N (x,  \boldsymbol{\xi}_a, t)$ 表示 $f^N$ 的前 $N$ 阶矩在数学上可行。
据此可定义离散化的分布函数 $f_a$ （ $a=1,2,...,d$ ），其对应的权重和特征方向为 $w_a$ 和 $\xi_a$ ：

$$f_a (x,t) = \frac{w_a}{\omega( \boldsymbol{\xi}_a)} f(x,  \boldsymbol{\xi}_a, t) \tag{23}$$

则宏观量可表达为：

$$\rho = \sum_{a=1}^{d} f_a , \rho \mathbf{u} =  \sum_{a=1}^{d} f_a \mathbf{\xi}_a ,  \rho \mathbf{uu} + \mathbf{P} = \sum_{a=1}^{d} f_a \mathbf{\xi}_a \mathbf{\xi}_a    \tag{24}$$

下面将(1)式的碰撞项换为(2)式，并把 $\frac{\mathbf{F}}{m} \frac{\partial f}{\partial \xi}$ 直接记作外力项 $-\mathbf{F} (\xi)$ ，得

$$\frac{\partial f}{\partial t} + \mathbf{\xi} \cdot \frac{\partial f}{\partial \mathit{x}} = \frac{-1}{\tau} ( f - f^{(eq)} ) + \mathbf{F} (\xi) \tag{25}$$

对式(25)，令 $\boldsymbol{\xi} = \boldsymbol{\xi}_a$ ，并在两端乘上 $\frac{w_a}{\omega( \boldsymbol{\xi}_a)}$ ，即可得到使用 $f_a$ 离散的方程（ $a=1,2,...,d$ ）：

$$\frac{\partial f_a}{\partial t} + \xi_a \cdot \nabla{f_a} = -\frac{1}{\tau} ( f_a - f_a^{(eq)} ) + F_a \tag{26-1}$$


$$f_a^{(eq)} = \frac{w_a}{\omega (\xi_a)} f^{(eq)} (\xi_a)  \tag{26-2}$$


$$F_a = \frac{w_a}{\omega (\xi_a)} \mathbf{F} (\xi_a)  \tag{26-3}$$

此时只需将 $f^{(eq)}$ 和 $F$ 进行类似的展开，即可完成对整个系统的离散。

## 附录A． 使用C-E展开方法将离散的LBGK方程还原为Navier-Stokes方程
LBE方程为(A-1.1)式。源项 $\boldsymbol{F}$ 和 $f^{(\mathrm{eq})}$ 为由郭照立、郑楚光和施保昌提出的格式，即(A-1.2)式和(A-1.3)式。[注意，这里的 $f$ 是密度分布函数]

$$f_i (x+ξ_i Δt,t+Δt)-f_i (x,t)=-\frac{1}{τ} \left[f_i (x,t)-f_i^{(eq)} (x,t) \right]+Δt⋅F_i (x,t)   \tag{A-1.1}$$


$$F_i=(1 - \frac{Δt}{2τ}) w_i \left[ \frac{\mathbf{\xi_i-u} }{c_s^2 } + \frac{\mathbf{(\xi_i \cdot u) \xi_i} }{c_s^4} \right] ⋅ F  \tag{A-1.2}$$

$$f_i^{(eq)} (x,t) = \rho w_i \left[ 1 + \frac{\mathbf{ξ_i⋅u}}{c_s^2} +\frac{ \mathbf{uu} : ( \mathbf{ξ_i ξ_i} - c_s^2 \mathbf{I}) }{2c_s^4} \right]  \tag{A-1.3}$$ 

并且对于 $f_i$ 和 $f_i^{(eq)}$ 有：

$$
\begin{cases} \rho = \sum_i {f_i}  \\ \rho \mathbf{u} = \sum_i {\mathbf{\xi_i} f_i} + \frac{\Delta t}{2} F \end{cases}$$ 

和

$$\begin{cases} \rho = \sum_i {f_i^{(eq)}}  \\ \rho \mathbf{u} = \sum_i {\mathbf{\xi_i} f_i^{(eq)}} \\ \rho (\mathbf{uu} + c_s^2 \mathbf{I}) = \sum_i {\mathbf{\xi_i \xi_i} f_i^{(eq)}} \end{cases}
$$ 

将(A-1.1)式在 $(\mathbf{x} , t)$处展开，取至 2 阶项得：

$$
\begin{aligned}
\Delta t \frac{\partial f_i}{\partial t} + ξ_i \Delta t ⋅ \frac{\partial f_i}{∂x} + &\\
\frac{1}{2} \left[\Delta t^2 \frac{\partial^2 f_i}{\partial t^2} + 2\Delta t(\xi_i \Delta t) ⋅ \frac{\partial^2 f_i}{\partial t∂x} + \Delta t^2 (\xi_i \xi_i): \frac{\partial^2 f_i}{∂x^2}  \right] + &\\
\frac{1}{\tau } \left[f_i - f_i^{(eq)} \right] - \Delta t F_i &= 0  \tag{A-2}
\end{aligned}
$$ 

代入多尺度展开等式，并舍去 $\varepsilon^3$ 项和更高阶的项，可得：

$$
\begin{aligned}
    (\varepsilon \frac{\partial f_i^{(0)}}{\partial t_1} + \varepsilon^2 \frac{\partial f_i^{(1)}}{\partial t_1} + \varepsilon^2  \frac{\partial f_i^{(0)}}{\partial t_2} ) + \xi_i ⋅ (\varepsilon \frac{\partial f_i^{(0)}}{\partial x_1} + \varepsilon^2  \frac{\partial f_i^{(0)}}{\partial x_1}) + &\\
    \frac{Δt}{2} \left[ \varepsilon^2 \frac{\partial^2 f_i^{(0)}}{\partial t_1^2} + 2 \xi_i \varepsilon^2 ⋅ \frac{\partial^2 f_i^{(0)}}{\partial t_1 \partial x_1} + (\xi_i \xi_i) :  (\varepsilon^2 \frac{\partial^2 f_i^{(0)}}{\partial x^2} ) \right] &\\
    +  \frac{1}{τΔt} \left[ f_i^{(0)} + \varepsilon f_i^{(1)} + \varepsilon^2 f_i^{(2)}  -f_i^{(eq)} \right] - \varepsilon F_i^{(1)} + O(\varepsilon^3 ) &=0   \tag{A-3}
\end{aligned}
$$ 

从(A-3)式中提取$O({\varepsilon}^0)$项、$O({\varepsilon}^1)$项和$O({\varepsilon}^2)$项，并令 $D_k = \frac{\partial}{\partial t_k} + \xi_k \cdot \frac{\partial}{\partial \mathbf{x}_k}$  。(A-3)式对任意 $\varepsilon$ 都成立，则关于各阶的系数均为0，即：

$$O({\varepsilon}^0) = \frac{1}{\tau \Delta t} (f_i^{(0)} - f_i^{(eq)}) = 0  \tag{A-4.1}$$


$$O(\varepsilon^1) = D_1 f_i^{(0)} + \frac{ f_i^{(1)} }{\tau \Delta t} - F_i^{(1)} = 0  \tag{A-4.2}$$


$$O(\varepsilon^2) = \frac{\partial f_i^{(0)}}{\partial t_2} -  (1 - \frac{1}{2\tau}) D_1 f_i^{(1)} + \frac{ f_i^{(2)} }{\tau \Delta t} + \frac{\Delta t}{2} D_1 F_i^{(1)} = 0  \tag{A-4.3}$$
 根据多尺度展开的变量关系，可得：
$\sum_i {f_i^{(1)}} = \sum_i {f_i^{(2)}} = 0$ ， $\sum_i { \xi_i f_i^{(2)}} = 0$ ， $\sum_i { \xi_i f_i^{(1)}} = - \frac{\Delta t}{2} F_i^{(1)}$ 

$$\left\{\begin{matrix} \sum_{i} {F_i} = 0 \\ \sum_{i} {\xi_i F_i} = (1 - \frac{1}{2 \tau}) F \\ \sum_{i} {\xi_i \xi_i F_i} = (1 - \frac{1}{2 \tau}) 2 u F \end{matrix}\right.  $$
 ，

$$\left\{\begin{matrix} \sum_{i} {F_i^{(1)}} = 0 \\ \sum_{i} {\xi_i F_i^{(1)}} = (1 - \frac{1}{2 \tau}) F^{(1)} \\ \sum_{i} {\xi_i \xi_i F_i^{(1)}} = (1 - \frac{1}{2 \tau}) 2 u F^{(1)} \end{matrix}\right.$$
 
为方便下文分析，记： $\Pi^{(k)} = \sum_i { \xi_i \xi_i f_i^{(k)} }$ ， $\Gamma^{(k)} = \sum_i { \xi_i \xi_i \xi_i f_i^{(k)} }$ 。
对 $O({\varepsilon}^1)$ 项求1阶和2阶速度矩，得：

$$\frac{\partial \rho}{\partial t_1} + \frac{\partial (\rho u)}{\partial x_1} = 0  \tag{A-5.1}$$


$$\frac{\partial (\rho u)}{\partial t_1} + \frac{\partial \Pi^{(0)}}{\partial x_1} = 0  \tag{A-5.2}$$


$$\Pi^{(0)} = \rho (\it{\bf{uu}} \mathrm{+} c_s^2 \bf{I})   \tag{A-5.3}$$

对 $O({\varepsilon}^2)$ 项求1阶和2阶速度矩，得：

$$\frac{\partial {\rho}}{\partial t_{2}} = 0 \tag{A-6.1}$$

$$\frac{\partial (\rho \boldsymbol{u})}{\partial t_{2}} + (1 - \frac{1}{2 \tau}) \frac{\partial \Pi^{(1)}}{\partial \boldsymbol{x}_{1}} + \Delta t (1 - \frac{1}{2 \tau}) \frac{\partial \boldsymbol{F}^{(1)}}{\partial \boldsymbol{x}_{1}} = 0 \tag{A-6.2}$$

由(A-4.2)式，得： $f_i^{(1)} = \tau \Delta t \{ D_1 f_i^{(0)}  - F_i^{(1)}\}$ ，所以

$$\Pi^{(1)} = \sum_{i} \boldsymbol{\xi}_{i} \boldsymbol{\xi}_{i} f_{i}^{(1)} = -\tau \Delta t [\frac{\partial (\rho \boldsymbol{u} \boldsymbol{u})}{\partial t_{1}} + c_s^2 \frac{\partial (\rho \boldsymbol{I})}{\partial t_{1}} + \frac{\partial \Gamma^{(0)}}{\partial \boldsymbol{x}_{1}} - (1 - \frac{1}{2 \tau}) 2\boldsymbol{uF^{(1)}} ] \tag{A-7}$$

$\boldsymbol{uF^{(1)}}$是Grad的论文中的记号， $(\bf{\it{uF}}^{(1)} )_{\it{lm}} = \mathrm{\frac{1}{2}} (\it{ u_l F_m^{(1)} + u_m F_l^{(1)} } \mathrm{)}$ 。下文为方便分析，取 $\Pi_{lm}^{(1)}$ 分析。

对第1项，有： $\frac{\partial (\rho u_l u_m)}{\partial t_1} = \rho u_l \frac{\partial u_m}{\partial t_1} +  u_m \frac{\partial (\rho u_l)}{\partial t_1}$ 。由(A-5.2)式，可得[^CAO_2019]：

$$\rho \frac{\partial u_m}{\partial t_1} + u_m \frac{\partial \rho}{\partial t_1} + \frac{\partial}{\partial x_{1n}} \left( \rho u_m u_n + \rho c_s^2 \delta_{mn} \right) = F_m^{(1)}  \tag{A-8.1}$$


$$\frac{\partial (\rho u_l)}{\partial t_1}+ \frac{\partial}{\partial x_{1n}} \left( \rho u_l u_n + \rho c_s^2 \delta_{ln} \right) = F_l^{(1)}  \tag{A-8.2}$$
 
并且

$$\frac{\partial (\rho u_l u_m u_n)}{\partial x_{1n}} = - u_l u_m  \frac{\partial (\rho u_n)}{\partial x_{1n}} + u_l \frac{\partial (\rho u_m u_n)}{\partial x_{1n}} + u_m \frac{\partial (\rho u_l  u_n)}{\partial x_{1n}} \tag{A-8.3}$$

所以
$$\frac{\partial (\rho u_l u_m)}{\partial x_{1n}} = (u_l F_m^{(1)} + u_m F_l^{(1)}) - \frac{\partial (\rho u_l u_m u_n)}{\partial x_{1n}} - c_s^2 \left( u_l \frac{\partial \rho}{\partial x_{1m}} + u_m \frac{\partial \rho}{\partial x_{1l}} \right) \tag{A-8.4}$$ 

所以第2项的分量为： $c_s^2 \delta_{lm} \frac{\partial \rho}{\partial t_1}$ 。

对第3项的分量，由于 
$$\Gamma_{lmn}^{(0)} = \rho c_s^2 ( u_l \delta_{mn} + u_n \delta_{lm} + u_m \delta_{ln} )$$
 ，所以：

$$\frac{\partial \Gamma_{lmn}^{(0)}}{\partial x_{1n}} = c_s^2 (  \delta_{mn} \frac{\partial \rho u_l}{\partial x_{1n}} + \delta_{lm} \frac{\partial \rho u_n}{\partial x_{1n}} + \delta_{ln} \frac{\partial \rho u_m}{\partial x_{1n}}  ) \tag{A-8.5}$$

所以

$$\Pi_{lm}^{(1)} = - \tau \Delta t   \left[ - \frac{\partial (\rho u_l u_m u_n)}{\partial x_{1n}} + \frac{1}{2 \tau} (u_l F_m^{(1)} + u_m F_l^{(1)})  + \rho c_s^2 ( \frac{\partial u_l}{\partial x_{1m}} +  \frac{\partial u_m}{\partial x_{1l}}  ) \right]  \tag{A-8.6}$$ 

略去$\Pi_{lm}^{(1)}$ 中的速度的3阶项（即 $\frac{\partial}{\partial x_{1n}} \{ \rho u_l u_m u_n \}$ ），有

$$\Pi^{(1)} = - \Delta t \cdot \boldsymbol{u F}^{(1)} - \tau \rho c_s^2 \Delta t \left[ \frac{\partial \boldsymbol{u}}{\partial \boldsymbol{x_1}}  + \left(  \frac{\partial \boldsymbol{u}}{\partial \boldsymbol{x_1}} \right)^{\mathrm{T}}  \right]  \tag{A-9}$$
 回代（A-6.2）式，得等价形式为

$$\frac{\partial (\rho \boldsymbol{u})}{\partial t_2} + c_s^2 \Delta t (\tau - \frac{1}{2}) \frac{\partial}{\partial \boldsymbol{x_1}} \left\{  \rho \left[  \frac{\partial \boldsymbol{u}}{\partial \boldsymbol{x_1}}  + \left(  \frac{\partial \boldsymbol{u}}{\partial \boldsymbol{x_1}} \right)^{\mathrm{T}}  \right] \right\} = \boldsymbol{0} \tag{A-10}$$

 将(A-10)式和(A-5.2)式联立，得

$$\frac{\partial (\rho u)}{\partial t} + \frac{\partial (\rho u u)}{\partial x} = - \frac{\partial (\rho c_s^2)}{\partial x} + c_s^2 \Delta t (\tau - \frac{1}{2}) \frac{\partial}{\partial x}\{\rho [\frac{\partial u}{\partial x} + (\frac{\partial u}{\partial x})^T ]\} + F   \tag{A-11}$$ 

此时令 $p = c_s^2 \rho$ ， $\nu = c_s^2 \Delta t (\tau - \frac{1}{2})$ ，即可还原回Navier-Stokes方程

## 附录B．Gauss-Hermite积分
对于任意一维函数 $f(\xi)$ ，高斯积分通过寻找积分点 $\xi_k$ 和其对应权重 $w_k$ 实现对 $\int_{a}^{b} \omega(\xi) f(\xi) d \xi$ 的近似，即

$$\int_{a}^{b} \omega(\xi) f(\xi) d \xi \approx \sum_{k=1}^{n} w_k \xi_k  \tag{B-1}$$ 

其中 $\omega (\xi)$ 是任意的权重函数。n 点高斯求积的积分点正是对应的第n个正交多项式的根，此时 $\xi_k$ 对应的权重 $w_k$ 为

$$w_k = \frac{ <P_{n-1}, P_{n-1}> }{ P_{n-1}(\xi_k) P_{n}^{'}(\xi_k) }   \tag{B-2}$$ 

其中 $P_{n}^{'} = \frac{d P_n}{d x}$ 。(B-2)式的精度为 $2n-1$ 。
所以在一维Gauss-Hermite积分中，取 $P_n = \mathcal{H}^{(n)}$ ， $\mathcal{H}^{(n)}$ 的权函数和表达式见(12)式和(13)式（下文同理）。 n点高斯求积的积分点即为 $\mathcal{H}^{(n)}$ 的根。因为 $\frac{d  \mathcal{H}^{(n)}}{d \xi} = \xi  \mathcal{H}^{(n)} -  \mathcal{H}^{(n+1)} =  n \mathcal{H}^{(n-1)}$ ，所以有 

$$w_k = n! / [n \mathcal{H}^{(n-1)} (\xi_k)]^2 \tag{B-3}$$ 

对于更高维的情况，可做类似构造[^SHAN等_2006]。分析如下积分

$$\frac{1}{(2 \pi)^{D/2}} \int \mathrm{exp}(- \frac{\xi \cdot \xi}{2}) p(\xi) \mathrm{d} \xi  \tag{B-4}$$

其中 $\xi$ 是长度为D的向量。 $p(\xi)$ 是D维n次多项式，可写为(B-5)式

$$p(\xi) = \sum\limits_{n_1 + n_2 + ... +n_D \le n}  {a_{n_1 n_2 ... n_D}} \prod\limits_{j=1}\limits^{D}{\xi_{j}^{n_j}}  \tag{B-5}$$

其中 $a_{n_1 n_2 ... n_D}$为常数。

考虑(B-5)式中的任意一项  $\prod\limits_{j=1}\limits^{D}{\xi_{j}^{n_j}}$ （由于$a_{n_1 n_2 ... n_D}$ 为常数所以可暂不处理，在最后一步被归并到常数项中），回代到式(B-4)中，得

$$\frac{1}{(2 \pi)^{D/2}} \int \mathrm{exp}(- \frac{\xi \cdot \xi}{2}) p(\xi) \mathrm{d} \xi = \prod\limits_{j=1}\limits^{D} \left( \frac{1}{\sqrt{2 \pi}} \int \mathrm{exp}(- \frac{\xi_{j}^2}{2}) \xi_j^{n_j} \mathrm{d}\xi_j \right) = \prod\limits_{j=1}\limits^{D} ( \sum\limits_{k=1}\limits^{n} w_k \xi_k )   \tag{B-6.1}$$

其中 $w_k$和 $\xi_k$是一维n次多项式高斯积分的权重和积分点。积分结果被转化为多个一维Gauss-Hermite积分的连乘。注意到(B-6.1)式等价于

$$\frac{1}{(2 \pi)^{D/2}} \int \mathrm{exp}(- \frac{\xi \cdot \xi}{2}) p(\xi) \mathrm{d} \xi = \sum\limits_{k_1=1}\limits^{n} ... \sum\limits_{k_D=1}\limits^{n} [(w_{k_1} w_{k_2} ... w_{k_D}) (\xi_{k_1}^{n_1} \xi_{k_2}^{n_2} ... \xi_{k_D}^{n_D})]  \tag{B-6.2}$$

如果定义 $\xi_{k_1 ... k_D} = (\xi_{k_1}, \xi_{k_2}, ... , \xi_{k_D})$ 和  $w_{k_1 ... k_D} = w_{k_1} w_{k_2} ... w_{k_D}$ ，那么(B-4)式等价为

$$\frac{1}{(2 \pi)^{D/2}} \int \mathrm{exp}(- \frac{\xi \cdot \xi}{2}) p(\xi) \mathrm{d} \xi =  \sum w_{k_1 ... k_D} p(\xi_{k_1 ... k_D})  \tag{B-7}$$

此时则将(B-4)式转换为与一维Gauss-Hermite积分类似的形式，即(B-7)式。


## 参考
 
[^SHAN等_2006]:SHAN X, YUAN X-F, CHEN H. Kinetic Theory Representation of Hydrodynamics: A Way beyond the Navier–Stokes Equation[J]. Journal of Fluid Mechanics, 2006, 550(1): 413.  http://doi.org/10.1017/S0022112005008153    
[^GRAD_1949]:GRAD H. Note on N-Dimensional Hermite Polynomials[J]. Communications on Pure and Applied Mathematics, 1949, 2(4): 325–330. http://doi.org/10.1002/cpa.3160020402  
[^GRAD_1949_LBM]:GRAD H. On the Kinetic Theory of Rarefied Gases[J]. Communications on Pure and Applied Mathematics, 1949, 2(4): 331–407.  http://doi.org/10.1002/cpa.3160020403  
[^CAO_2019]:CAO W. Investigation of the Applicability of the Lattice Boltzmann Method to Free-Surface Hydrodynamic Problems in Marine Engineering[D/OL]. Laboratoire de recherche en Hydrodynamique, Énergétique et Environnement Atmosphérique (LHEEA): École centrale de Nantes, 2019.  https://tel.archives-ouvertes.fr/tel-02383174
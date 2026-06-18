---
title: 【笔记】LBM相场模型中的序参数梯度差分
tags:
  - 格子Boltzmann方法
  - 相场模型
categories:
  - 格子Boltzmann方法
  - 相场模型
date: 2026-06-13 20:39:20
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

在 ["【笔记】格子Boltzmann方法的相场模型"](https://ives0721.github.io/2026/04/03/LBM-phaseField/) 中，我把LBM里的相场模型简单记了一下。笔记里有提到，序参数 $\phi$ 的梯度可以通过下面的差分进行计算。

$$
\displaystyle
\begin{cases}
    \nabla\phi(\boldsymbol{x}) =& \sum\limits_{i \neq 0} \dfrac{w_i \boldsymbol{c}_i \phi(\boldsymbol{x} + \boldsymbol{c}_i \delta_t)}{c_s^2 \delta_t} ,\\
    \nabla^2\phi(\boldsymbol{x}) =& \sum\limits_{i \neq 0} \dfrac{2w_i}{c_s^2 \delta_t^2}  [\phi(\boldsymbol{x} + \boldsymbol{c}_i \delta_t) - \phi(\boldsymbol{x})]
\end{cases}
$$

其中 $i$ 为 DnQb 离散速度模型中的第 $i$ 个离散速度 $\boldsymbol{c}_i$ 的下标， $w_i$ 是 $\boldsymbol{c}_i$ 的权重系数。

## 序参数 Φ 的一阶导数

对于点 $\boldsymbol{x}$ 及其邻近格子 $\boldsymbol{x}+\boldsymbol{c}_i \delta_t$，考虑序参数 $\phi$ 在点 $\boldsymbol{x}$ 的Taylor展开，得：

$$
\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) = \phi(\boldsymbol{x}) + (\boldsymbol{c}_i \delta_t) \cdot \nabla\phi(\boldsymbol{x}) + \frac{\delta_t^2}{2} (\boldsymbol{c}_i \cdot \nabla)^2 \phi(\boldsymbol{x}) + O(\delta_t^3)
\tag{1}
$$

式(1)两边同时乘上 $w_i \boldsymbol{c}_i$，并对所有方向求和，得到：

$$
\begin{aligned}
    \sum_i w_i \boldsymbol{c}_i \phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) &=
    \sum_i w_i \boldsymbol{c}_i \phi(\boldsymbol{x})  \\
    &\quad + \sum_i w_i \boldsymbol{c}_i [(\boldsymbol{c}_i \delta_t) \cdot \nabla\phi(\boldsymbol{x})]  \\
    &\quad + \sum_i w_i \boldsymbol{c}_i \left[
        \frac{\delta_t^2}{2} (\boldsymbol{c}_i \cdot \nabla)^2 \phi(\boldsymbol{x})
    \right] + O(\delta_t^3)
\end{aligned}
\tag{2}
$$

下面分析式(2)右侧的三个项。

（1）对第1项，$\phi(\boldsymbol{x})$ 是定值，所以

$$
\sum_i w_i \boldsymbol{c}_i \phi(\boldsymbol{x}) = \phi(\boldsymbol{x}) \cdot \sum_i (w_i \boldsymbol{c}_i )
$$

由于 $\sum_i w_i \boldsymbol{c}_i = \vec{0}$，所以第1项为 0 。

（2）对第2项，

$$
\sum_i w_i \boldsymbol{c}_i [(\boldsymbol{c}_i \delta_t) \cdot \nabla\phi(\boldsymbol{x})] = \delta_t \nabla\phi(\boldsymbol{x}) \cdot \sum_i (w_i \boldsymbol{c}_i \boldsymbol{c}_i) = \delta_t c_s^2 \nabla\phi(\boldsymbol{x})
$$

其中 $\sum_i (w_i \boldsymbol{c}_i \boldsymbol{c}_i) = c_s^2 [\boldsymbol{I}]$ 。

这个关系可以由常规 DnQb 模型中的 $\sum_i \boldsymbol{c}_i \boldsymbol{c}_i f_i^{eq} = \rho (\boldsymbol{u}\boldsymbol{u} + c_s^2 [\boldsymbol{I}])$ 推出。（取 $\boldsymbol{u}=\vec{0}$）

（3）对第3项，为简化考虑，这里推导其中一个分量，其形式为：

$$
\frac{\delta_t^2}{2} \sum_i w_i c_{i\gamma} \left[
    c_{i\alpha} c_{i\beta} \frac{\partial^2 \phi(\boldsymbol{x})}{\partial x_{\alpha} \partial x_{\beta}}
\right] = 
\frac{\delta_t^2}{2} \frac{\partial^2 \phi(\boldsymbol{x})}{\partial x_{\alpha} \partial x_{\beta}} \sum_i w_i c_{i\alpha} c_{i\beta} c_{i\gamma}
$$

由于在常规 DnQb 模型中，

$$
\sum_i c_{i\alpha} c_{i\beta} c_{i\gamma} f_i^{eq} = c_s^2 \rho [u_{\alpha} \delta_{\beta\gamma}+ u_{\beta} \delta_{\alpha\gamma} + u_{\gamma} \delta_{\alpha\beta}]
$$

这里的 $\delta$ 为克罗内克函数。同样取 $\boldsymbol{u}=\vec{0}$ ，得：

$$
\sum_i w_i c_{i\alpha} c_{i\beta} c_{i\gamma} = c_s^2 [u_{\alpha} \delta_{\beta\gamma}+ u_{\beta} \delta_{\alpha\gamma} + u_{\gamma} \delta_{\alpha\beta}] = 0
$$

所以式(2)被化简为

$$
\sum_i w_i \boldsymbol{c}_i \phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) =
\delta_t c_s^2 \nabla\phi(\boldsymbol{x}) + O(\delta_t^3)
$$

即：

$$
\nabla\phi(\boldsymbol{x}) = \frac{1}{\delta_t c_s^2} \sum_i w_i \boldsymbol{c}_i \phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) + O(\delta_t^2)
\tag{3}
$$

## 序参数 Φ 的二阶导数

将式(1)改写成

$$
\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) - \phi(\boldsymbol{x}) = (\boldsymbol{c}_i \delta_t) \cdot \nabla\phi(\boldsymbol{x}) + \frac{\delta_t^2}{2} (\boldsymbol{c}_i \cdot \nabla)^2 \phi(\boldsymbol{x}) + O(\delta_t^3)
\tag{4}
$$


式(4)两边同时乘上 $w_i$，并对所有方向求和，得到：

$$
\begin{aligned}
    \sum_i w_i [\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) - \phi(\boldsymbol{x})] &= \sum_i w_i [(\boldsymbol{c}_i \delta_t) \cdot \nabla\phi(\boldsymbol{x})]  \\
    &\quad + \sum_i w_i \left[
        \frac{\delta_t^2}{2} (\boldsymbol{c}_i \cdot \nabla)^2 \phi(\boldsymbol{x})
    \right] + O(\delta_t^3)
\end{aligned}
\tag{5}
$$

参照上一节的推导，易得式(5)等式右侧的两项分别满足

$$
\sum_i w_i c_{i\alpha} \delta_t \nabla_{\alpha}\phi(\boldsymbol{x}) =
\delta_t \nabla_{\alpha}\phi(\boldsymbol{x}) \sum_i w_i c_{i\alpha} = 0
$$

以及

$$
\sum_i w_i \frac{\delta_t^2}{2} c_{i\alpha} c_{i\beta} \cdot \nabla_{\alpha\beta}\phi(\boldsymbol{x}) = \frac{\delta_t^2}{2} \nabla_{\alpha\beta}\phi(\boldsymbol{x}) \sum_i w_i c_{i\alpha} c_{i\beta} = \frac{\delta_t^2 c_s^2 \delta_{\alpha\beta}}{2} \nabla_{\alpha\beta}\phi(\boldsymbol{x})
$$

其中 $\nabla_{\alpha} = \frac{\partial}{\partial x_{\alpha}}$ 和 $\nabla_{\alpha\beta} = \frac{\partial^2}{\partial x_{\alpha}\partial x_{\beta}}$ 。而 $\delta_{\alpha\beta}$ 同样是克罗内克函数。

所以式(5)表示为

$$
\sum_i w_i [\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) - \phi(\boldsymbol{x})] = \frac{\delta_t^2 c_s^2}{2} \nabla^2\phi(\boldsymbol{x}) + O(\delta_t^3)
$$

即：

$$
\nabla^2\phi(\boldsymbol{x}) = \frac{2}{\delta_t^2 c_s^2} \sum_i w_i [\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) - \phi(\boldsymbol{x})] + O(\delta_t^2)
\tag{6}
$$


## 关于中心差分

假设我们采用中心差分形式，需要对 $\boldsymbol{x}-\boldsymbol{c}_i \delta_t$ 进行 Taylor 展开：

$$
\phi(\boldsymbol{x}-\boldsymbol{c}_i \delta_t) = \phi(\boldsymbol{x}) - (\boldsymbol{c}_i \delta_t) \cdot \nabla\phi(\boldsymbol{x}) + \frac{\delta_t^2}{2} (\boldsymbol{c}_i \cdot \nabla)^2 \phi(\boldsymbol{x}) + O(\delta_t^3)
\tag{7}
$$

（A）将式(1)与式(7)相减，偶数阶项（零阶项 $\phi(\boldsymbol{x})$ 和二阶项）会被消去，得到：
$$
\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) - \phi(\boldsymbol{x}-\boldsymbol{c}_i \delta_t) = 2(\boldsymbol{c}_i \delta_t) \cdot \nabla\phi(\boldsymbol{x}) + O(\delta_t^3)
$$

两边同时乘上 $w_i \boldsymbol{c}_i$ 并对所有方向求和：
$$
\sum_i w_i \boldsymbol{c}_i [\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) - \phi(\boldsymbol{x}-\boldsymbol{c}_i \delta_t)] = 2\delta_t \nabla\phi(\boldsymbol{x}) \cdot \sum_i (w_i \boldsymbol{c}_i \boldsymbol{c}_i) + O(\delta_t^3)
$$

利用前面提到的二阶矩性质 $\sum_i w_i \boldsymbol{c}_i \boldsymbol{c}_i = c_s^2 \boldsymbol{I}$，等式右侧化简为 $2\delta_t c_s^2 \nabla\phi(\boldsymbol{x})$。因此：

$$
\nabla\phi(\boldsymbol{x}) = \frac{1}{2\delta_t c_s^2} \sum_i w_i \boldsymbol{c}_i [\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) - \phi(\boldsymbol{x}-\boldsymbol{c}_i \delta_t)] + O(\delta_t^2)
\tag{8}
$$

式(8)就是 Lee 等在文献[1]的式(64)的格式。

（B）将式(1)与式(7)相加，奇数阶项（一阶项和三阶项）会被消去，得到：
$$
\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) - 2 \phi(\boldsymbol{x}) + \phi(\boldsymbol{x}-\boldsymbol{c}_i \delta_t) = \delta_t^2 (\boldsymbol{c}_i \cdot \nabla)^2 \phi(\boldsymbol{x}) + O(\delta_t^3)
$$

两边同时乘上 $w_i $ 并对所有方向求和：
$$
\sum_i w_i [\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) - 2 \phi(\boldsymbol{x}) + \phi(\boldsymbol{x}-\boldsymbol{c}_i \delta_t)] = \sum_i \delta_t^2 (\boldsymbol{c}_i \cdot \nabla)^2 \phi(\boldsymbol{x}) + O(\delta_t^3)
$$

同理可得：
$$
\nabla^2 \phi(\boldsymbol{x}) = \frac{1}{c_s^2 \delta_t^2} \sum_i w_i [\phi(\boldsymbol{x}+\boldsymbol{c}_i \delta_t) - 2 \phi(\boldsymbol{x}) + \phi(\boldsymbol{x}-\boldsymbol{c}_i \delta_t)] + O(\delta_t^2)
\tag{9}
$$

式(9)就是 Lee 等在文献[1]的式(65)的格式。

[1] Lee T，Lin C-L. **A stable discretization of the lattice Boltzmann equation for simulation of incompressible two-phase flows at high density ratio**[J]. _Journal of Computational Physics_，2005，206（1）：16-47. DOI:10.1016/j.jcp.2004.12.001.


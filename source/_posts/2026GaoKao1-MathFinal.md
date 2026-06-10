---
title: 2026GaoKao1-MathFinal
date: 2026-06-10 20:47:24
tags:
categories:
- [杂谈]
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



> 最近看到2026年的全国Ⅰ卷高考数学压轴题，就做了一下。几乎全程用反证法做推导的题目实在是太罕见了，特此做一下记录。

# 题目

已知定义在 $\mathbb{R}$ 上的函数 $f(x)$ 满足：当 $x < 0$ 时， $f(x) = 2^x$。 结合 $D(x_0) =\{ d \in \mathbb{R} \vert f(x_0 + d) > f(x_0) \}$。

（1） 若当 $x \geq 0$， $f(x)=1-x$，求 $D(-1)$；  

（2） 若 $f(x)$ 为奇函数， $x_1 x_2 \neq 0$， $f(x_1) < f(x_2)$，证明： $D(x_1) \supseteq D(x_2)$；  

（3） 若 $f(x)$ 满足： ①当 $f(x_1) < f(x_2)$ 时， $D(x_1) \supseteq D(x_2)$； ②当 $0 < x < 1$ 时， $f(x) < f(0)$。证明：
- (i) $f(0) \geq 1$；  
- (ii) $f(x)$ 在区间 $(0, +\infty)$ 上单调递增。  


# 求解

## 小题 (1)

因为 $f(-1) = 1/2$，所以这里先计算 $f(x) > 1/2$ 在两边的解。

- 当 $x < 0$ 时，解得： $x \in (-1, 0)$。  
- 当 $x \geq 0$ 时，解得： $x \in [0, 1/2)$。

$f(x) > 1/2$ 的解便是 $x \in (-1, 1/2)$。

根据集合 $D(x_0)$ 的定义，我们就有：

$$
\begin{aligned}
D(-1) &= \{ d \in \mathbb{R} \vert f(-1 + d) > f(-1) \} \\
&= \{ d \in \mathbb{R} \vert -1 < -1 + d < 1/2 \} \\
&= \{ d \in \mathbb{R} \vert d \in (0, 3/2) \} \\
\end{aligned}
$$

## 小题 (2)

根据题意，此时的 $f(x)$ 写作：

$$
f(x) = \begin{cases}
2^x & x < 0\\
0 & x=0 \\
-2^{-x} & x>0
\end{cases}
$$

这里得强调一下 $f(x)$ 在 $x<0$ 和 $x>0$ 时都是增函数（等会要用！）。

由 $x_1 x_2 \neq 0$ ，我们知道这俩数肯定不是 0 。那么接下来的讨论就只有两个数同号或异号的情况。

【情况1】 $x_1$ 和 $x_2$ 都小于 0 。

由于 $f(x_1) < f(x_2)$， 自然也就有 $x_1 < x_2 < 0$ 。根据定义，我们可以写出
$$
\begin{cases}
D(x_1) &= \{ d \in \mathbb{R} \vert d \in (0, -x_1) \} \\
D(x_2) &= \{ d \in \mathbb{R} \vert d \in (0, -x_2) \} \\
\end{cases}
$$
明显有 $D(x_1) \supseteq D(x_2)$ 。

【情况2】 $x_1$ 和 $x_2$ 都大于 0 。

根据情况1的推理，同理有 $x_2 > x_1 > 0$，并且此时 $f(x_1) < f(x_2) < 0$ 。

打个比方，对于 $D(x_1)$：
- $(x_1 + d) \in (x_1, +\infty)$ 是 $D(x_1)$ 的一部分；
- 同时, $(x_1 + d) \in (-\infty,0]$ 是 $D(x_1)$ 的另一部分。

因此 $D(x_1) = \{ d \in \mathbb{R} \vert d \in (-\infty, -x_1] \cup (0, +\infty) \}$。

同理 $D(x_2) = \{ d \in \mathbb{R} \vert d \in (-\infty, -x_2] \cup (0, +\infty) \}$。

所以 $D(x_1) \supseteq D(x_2)$ 也是成立的。

【情况3】 $x_1 > 0$ 和 $x_2 < 0$ 。

由于 $f(x_1) < f(x_2)$， 所以 $|x_1| > |x_2|$。

根据前面两部分的推理，有

$$
\begin{cases}
D(x_1) &= \{ d \in \mathbb{R} \vert d \in (-\infty, -x_1] \cup (0, +\infty) \} \\
D(x_2) &= \{ d \in \mathbb{R} \vert d \in (0, -x_2) \} \\
\end{cases}
$$

所以 $D(x_1) \supseteq D(x_2)$ 也是成立的。

证明完毕。

## 小题 (3)

> 到这里就很有意思了。
> 两个小问要直接证明的话，相当不好走通。只能上**反证法**了。
> 这也是我对它感兴趣的原因。

### 小题 (3)(i)

首先，我们假设 $f(0) < 1$ 成立。

由于 $f(x)$ 在 $x>0$ 的形式不知道，这里我们就关注 $x<0$ 区间。该区间必然存在一个点 $x_0 < 0$ ，使得 $f(x_0) \geq f(0)$ 成立。


根据题目给的条件①，当 $f(0) < f(x_0)$ ，则 $D(0) \supseteq D(x_0)$ 。对于任意 $d \in D(x_0)$，同样应有 $d \in D(0)$

所以，这里我们取 $\epsilon \in (0, \min\{1, |x_0|\})$，此时 $\epsilon$ 同时满足 $(x_0 + \epsilon) < 0$ 和 $\epsilon < 1$。

因为 $f(x_0 + \epsilon) > f(x_0)$，所以 $\epsilon \in D(x_0)$。但是，根据题目给的条件②， $f(\epsilon) < f(0)$，所以 $\epsilon \notin D(0)$。

此时出现了与 $D(0) \supseteq D(x_0)$ 相矛盾的结论。因此，假设不成立。

所以，也就只剩下 $f(0) \geq 1$ 的情况。

### 小题 (3)(ii)

> 这里是最离谱的。  
> 我一开始习惯性在 $x\in(0,1), y>0$ 的地方画了两条线，一条单调增，一条单调减。结果两条线都推出了相互矛盾的东西。结果往 $y<0$ 的地方画就没什么问题。  
> 合着居然还要先推断 $f(x)$ 在 $x>0$ 时究竟是正还是负，甚至还是继续用**反证法**。  
> 做到这里的同学有苦头了。  

接下来，当 $x_0 > 0$ 时，由定义便有 $D(x_0) =\{ d \in \mathbb{R} \vert f(x_0 + d) > f(x_0) \}$ 。

【Part A】

假设存在点 $x_0 > 0$ 能够满足 $f(x_0) > 0$，如果这个都满足不了那就说明 $f(x) \leq 0$ 在 $x>0$ 的情况下恒成立。

由于当 $x<0$ 时， $f(x) \in (0, 1)$，所以必然存在 $t < 0$，使得 $f(t) < f(x_0)$。那么 $D(x_0) \subseteq D(t)$。

考虑 $D(t) =\{ d \in \mathbb{R} \vert f(t + d) > f(t) \}$ 。
- 当 $(t+d) \leq 0$时，解得 $d \in (0, -t]$。  
- 而当 $(t+d) > 0$ 时，由于 $f(0) > 1$，解得 $(t+d) \in D(0)$。  

记 $\mathbb{X}(t)=\{d \in \mathbb{R} \vert \exist \xi \in D(0), d = \xi - t \}$。所以， $D(t) = (0, -t] \cup \mathbb{X}(t)$ 。  
【注：$\mathbb{X}(t)$ 其实就是把 $D(0)$ 向右平移了 $|t|$ 个单位，后面我还用了一次这个表示来简化描述】

这里需要强调：由于题(3)(i)里面证明了 $f(0) \geq 1$，所以 $D(0)$ 一定是 $(0, +\infty)$ 的子集。也就是只能在 $d>0$ 的区间内才有可能找到 $f(d) > f(0)$ 的情况。
因此 $\mathbb{X}(t) \subseteq (0, +\infty)$ ，得出 $D(x_0) \subseteq D(t) \subseteq (0, +\infty)$。

【Part B-1】

回过头看 $D(x_0) = \{ d \in \mathbb{R} \vert f(x_0+d) > f(x_0)\}$ 。

考虑一个特殊点 $d_1=-x_0 < 0$，此时 $f(x_0+d_1) = f(0)$。  
如果 $x_0 \in (0, 1)$，则根据题(3)的条件②应该有 $f(x_0) < f(0) = f(x_0+d_1)$，也就是说 $d_1=-x_0 < 0$ 应该属于 $D(x_0)$。  
但我们前面说了 $D(x_0) \subseteq (0, +\infty)$，产生了矛盾。  
所以， $f(x_0)$ 在 $x_0 \in (0, 1)$ 不应该是正数。

【Part B-2】

接下来考虑 $x_0 \geq 1$ 的情况，只剩这个区间可以来假设 $f(x_0) > 0$ 的可能。  
还是考虑 $d_1=-x_0 < 0$ 这个一定不属于 $D(x_0)$ 的特殊点。把 $d_1$ 代入，应该有 $f(x_0 + d_1) = f(0) \leq f(x_0)$ ，即 $f(x_0) \geq 1$ 。

任意从负数中选取 $\gamma_1 < 0$，则必然有 $f(\gamma_1) < f(x_0)$ 和 $D(x_0) \subseteq D(\gamma_1)$。  
【注：其实是只能从这个已知信息多的地方入手了，并且还要仿照题(3)(i)的解法，把一个点从负数加到 $(0,1)$ 区间】

基于前面的内容，很容易写出 $D(\gamma_1) = (0, -\gamma_1] \cup \mathbb{X}(\gamma_1)$ 。  
令 $d_2 = x_0 - \gamma_1$，显然 $d_2 > 0$，并且 $f(\gamma_1 + d_2) = f(x_2) \geq 1 > f(\gamma_1)$。 所以 $d_2$ 属于 $D(\gamma_1)$。

构造点 $\gamma_2 \in (-d_2, 1 - d_2)$，此时的 $\gamma_2$ 同时满足 $\gamma_2 < \gamma_1$ 和 $(\gamma_2 + d_2) \in (0, 1)$。  
那么，根据题(3)的条件①应该有 $D(\gamma_1) \subseteq D(\gamma_2)$。 也就有 $d_2$ 也属于 $D(\gamma_2)$。  
但是，根据题(3)的条件②， $f(\gamma_2 + d_2) < f(0)$，对应 $d_2$ 不属于 $D(\gamma_2)$，产生了矛盾。

所以 $f(x_0)$ 在 $x_0 \geq 1$ 时同样不是正数。

【Part C】

经过了 Part B 的逻辑推演， $f(x_0) \leq 0$ 在 $x_0 \in (0, +\infty)$ 的区间内恒成立。相信大伙做到这里完全力竭了吧。但是！接下来就不需要反证法啦！

还是用这个 $x_0 \in (0, +\infty)$ 。此时我们有 $f(x_0) < 0 < f(t)$，其中 $t$ 为任意负数。那么 $D(t) \subseteq D(x_0)$。

$D(t)$ 的形式我们在 Part A 部分写过了，即： $D(t) = (0, -t] \cup \mathbb{X}(t)$ 。 
所以 $(0, -t] \subseteq D(x_0)$。
因此，在 $d \in (0, -t]$ 时， $f(x_0 + d) > f(x_0)$。

由于 $t$ 是任意负数，上面的 $f(x_0 + d) > f(x_0)$ 则可推导出 $f(x)$ 在区间 $x \in (x_0, +\infty)$ 上单调递增。

再结合 $x_0 \in (0, +\infty)$ 的任意性，可知 $f(x)$ 在区间 $x \in (0, +\infty)$ 上单调递增。

证明完毕。

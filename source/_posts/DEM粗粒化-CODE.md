---
title: 【代码】DEM 的颗粒流数据粗粒化
date: 2024-08-10 11:01:48
tags: [离散单元法, 粗粒化]
categories:
- [离散单元法]
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

基于[这篇](/2024/07/30/DEM粗粒化/)的内容，我用 Python 写了一份对**单种颗粒系统**进行粗粒化的代码。


其具体使用流程为：
1. `obj = CG_mono_3D(...)`: 初始化。
2. `obj.set_domain(...)`: 设置矩形的计算域。
3. `obj.get_CG_data(R, CG_width)`: 对 `R` 点进行粗粒化，其粗粒化的长度为 `CG_width`。


<details>
<summary>DEM_CG_mono.py</summary>
这里使用<a href="https://gist.github.com/Ives0721/6796724c0de90bb03d5236af27d4eeb9">Github Gist</a>展示代码。如果无法查看下面内容，青科学上网。
<script src="https://gist.github.com/Ives0721/6796724c0de90bb03d5236af27d4eeb9.js"></script>
</details>

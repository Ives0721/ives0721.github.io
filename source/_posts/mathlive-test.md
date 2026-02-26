---
title: 在 Hexo 中试用 Mathlive 库
date: 2026-02-23 14:36:35
tags: [Mathlive]
categories:
- [杂谈]
---
<script defer src="https://cdn.jsdelivr.net/npm/mathlive"></script>
<script>
window.addEventListener("DOMContentLoaded", () =>
    MathLive.renderMathInDocument(),
);
</script>


最近刚刚刷到一个叫 [Mathlive](https://mathlive.io/) 的库，突然一想 Hexo 也是把 Markdown 转成 html 的，那就尝试了一下 Mathlive 的效果。

# 导入

我只是在 Markdown 开头补了一个从 CDN 加载的 `<script>`，然后又从官方的[script tag](https://mathlive.io/mathfield/guides/integration/#using-a-script-tag)方式里面抄了一小段：

```html
<script defer src="https://cdn.jsdelivr.net/npm/mathlive"></script>
<script>
window.addEventListener("DOMContentLoaded", () =>
    MathLive.renderMathInDocument(),
);
</script>
```

# 尝试

从它的[文档](https://mathlive.io/mathfield/reference/commands/)来看，支持的LaTeX命令好像挺多的。


## math-span

输入：`<math-span>e^{i\pi} + 1 = 0</math-span>`。   
效果为：<math-span>e^{i\pi} + 1 = 0</math-span>

我也试了使用 `$$e^{i\pi} + 1 = 0$$`，  
结果得到 $$e^{i\pi} + 1 = 0$$

【注，用 `$$` 得出的默认渲染结果是会自动跨行的，用 `$` 无法渲染。】

## math-field

输入（关于 styling 可以查看[对应文档](https://mathlive.io/mathfield/guides/customizing/)）：

```html
<math-field style="font-size:1rem; display: block">
  x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}
</math-field>
```

结果：
<math-field style="font-size:1rem; display: block">
  x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}
</math-field>

# 一个简单的可视化 LaTeX 公式编辑器

详情参见 https://mathlive.io/mathfield/guides/interacting/

目前它比较简陋，仅在 PC 浏览器下的渲染效果较好。

<div>
    <math-field id="formula" style="font-size:1.5rem; display: block; max-width:100%"> </math-field>
    <textarea id="latex" autocapitalize="off" autocomplete="off" autocorrect="off" spellcheck="false" style="display: flex;width:100%"></textarea>
    <script>
        const mf = document.getElementById("formula");
        const latex = document.getElementById("latex");
        mf.addEventListener("input",(ev) => latex.value = mf.value);
        latex.value = mf.value;
        // Listen for changes in the "latex" text field,
        // and reflect its value in the mathfield.
        latex.addEventListener("input", (ev) =>
            mf.setValue( ev.target.value, {silenceNotifications: true} )
        );
    </script>
</div>


<details>

<summary>html 实现</summary>

```html
<div>
    <math-field id="formula" style="font-size:1.5rem; display: block; max-width:100%"> </math-field>
    <textarea id="latex" autocapitalize="off" autocomplete="off" autocorrect="off" spellcheck="false" style="display: flex;width:100%"></textarea>
    <script>
        const mf = document.getElementById("formula");
        const latex = document.getElementById("latex");
        mf.addEventListener("input",(ev) => latex.value = mf.value);
        latex.value = mf.value;
        // Listen for changes in the "latex" text field,
        // and reflect its value in the mathfield.
        latex.addEventListener("input", (ev) =>
            mf.setValue( ev.target.value, {silenceNotifications: true} )
        );
    </script>
</div>
```

</details>

---
title: RDPWrap 安装记录
date: 2024-10-28 09:06:30
tags: [RDPWrap]
categories:
- [杂谈]
---

工位上换了新电脑，于是需要重新安装 RDPWrap 让机子上的 Windows 11 家庭版能够打开 RDP 远程连接。

1. 安装 RDPWrap 

在 [RDPWrap 的 Github 网址](https://github.com/stascorp/rdpwrap)上下载它的最新 Release 【其实已经很久没更新了】。 解压后运行 `install.bat` 进行安装。

2. 出现 `Listener state: Not-listening` 的解决方案。

1）首先在 [sebaxakerhtc/rdpwrap.ini](https://github.com/sebaxakerhtc/rdpwrap.ini/) 这个仓库里下一份最新的 `rdpwrap.ini`，一般来说大多数Windows OS 版本号它都有配置。  
2）然后将自行下载的 `rdpwrap.ini` 替换掉原始文件[^1]，原始文件位于 `C:\Program Files\RDP Wrapper` 文件夹。   
3）替换完成后，回到解压 RDPWrap 的文件夹并在此处打开终端，运行 `.\RDPWInst.exe -r` [^2]。这之后我的电脑就设置完成了。


[^1]: 执行此操作前先在**具有管理员权限的终端**里把对应远程服务关了：`net stop TermService`。
[^2]: 这里参考了 Github 上的 [issue 999](https://github.com/stascorp/rdpwrap/issues/999)。


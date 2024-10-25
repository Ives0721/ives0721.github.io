---
title: 我的 WSL 2 安装记录
date: 2024-10-23 09:27:20
tags: [WSL]
categories:
- [杂谈]
---

工位上换了新电脑，于是需要给新电脑安装 Linux 子系统（WSL 2）用来跑代码。  
这里简单记录一下我的安装过程。

1. **安装 WSL 2**

这里主要是按照微软官方的教程[^MS_WSL_1][^MS_WSL_2]来做。

先在 powershell 里启用 WSL 和虚拟机功能



```powershell
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart  
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart  
```

然后安装[适用于 x64 计算机的 WSL2 Linux 内核更新包](https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi)，并在 powershell 里将 WSL 2 设置为默认。

```powershell
wsl --set-default-version 2
```

WSL 2 的自定义安装位置参考的是[CSDN](https://blog.csdn.net/aaada123/article/details/142643762)上给的做法：我下载了[Ubuntu 24.04 的 AppxBundle](https://wslstorestorage.blob.core.windows.net/wslblob/Ubuntu2404-240425.AppxBundle)，从包里提取 `Ubuntu_2404.0.5.0_x64.Appx` ，然后在我想要的位置将它当成 zip 解压。最后运行一下 `ubuntu2404.exe` 即可。

2. **安装CUDA**

[Nvidia 的安装教程](https://docs.nvidia.com/cuda/wsl-user-guide/index.html#getting-started-with-cuda-on-wsl-2) 和它给的[下载链接](https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=WSL-Ubuntu&target_version=2.0&target_type=deb_local)主要是针对 WSL-Ubuntu 的。

我用的网络下载，在shell里就应该是：
```shell
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get -y install cuda-toolkit-12-6
sudo apt-get -y install cuda
```

即使是这么安装完，还是会有非常幽默的 `nvcc` 找不到、 `libcuda.so` 找不到等问题。我就参照网上的做法[^WSL_cuda1][^WSL_cuda2]在 `~/.bashrc` 里补充了这么一段：

```bash
# >>> Nvidia CUDA initialize >>>
export CUDA_HOME=/usr/local/cuda
export PATH=${CUDA_HOME}/bin:${PATH}
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
export LD_LIBRARY_PATH=/usr/lib/wsl/lib/:${LD_LIBRARY_PATH}
# <<< Nvidia CUDA initialize <<<
```

搞定。

[^MS_WSL_1]: https://learn.microsoft.com/zh-cn/windows/wsl/install  
[^MS_WSL_2]: https://learn.microsoft.com/zh-cn/windows/wsl/install-manual  
[^WSL_cuda1]: https://stackoverflow.com/questions/68221962/nvcc-not-found-but-cuda-runs-fine  
[^WSL_cuda2]: https://github.com/microsoft/WSL/issues/8587  
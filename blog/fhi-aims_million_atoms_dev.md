---
title: "FHI-aims百万原子开发协作日志"
date: 2021-06-07 12:09:26 +0800
---

# 开发约定

## 环境平台

- **超算平台**：神威太湖之光

- **程序路径**：`/home/export/online3/para030/wuyangjun/07_aims_million_dev/fhi-aims_MPE_O3_local-index`

  - git clone到自己目录

  ```shell
  git clone /home/export/online3/para030/wuyangjun/07_aims_million_dev/fhi-aims_MPE_O3_local-index
  ```

  - git pull 、git push进行同步，合并，推送等
  - git push 报错可以在主仓库执行配置以下内容

  ```shell
  git config --global receive.denyCurrentBranch updateInstead
  ```

- **测试路径**：`/home/export/online3/para030/wuyangjun/07_aims_million_dev/sw5_test`
  
  - `H2_mini`：氢分子
  - `H2O`：水分子
  - `RBDmini`：RBD蛋白

## 编译环境

- **必须加载710编译器**：`module load sw/compiler/gcc710`

## 精度要求

FHI-aims现有两个版本的程序，大内存版本（global index）和小内存版本（local index），它们的切换可在`control.in`进行配置：
- **local index**：在control.in加上以下配置

```shell
use_local_index .true.
load_balancing .true.
```

- **global index**：在control.in注释以下配置

```shell
#use_local_index .true.
#load_balancing .true.
```

**而local index会引入精度误差，为保证精度一致性，测试分两步骤：**

- 步骤1：先在**global index**下开发测试，测试都确保**无精度损失**后
- 步骤2：再切换为**local index**测试，只要精度误差（和原版global index比）不超过**原始程序**的global index和local index二者的误差即可。
  - 开发后有4个程序版本，原始程序的global index(gidx)和local index(lidx)，优化后的global index(gidx')和local index(lidx')
  - 步骤1保证gidx和gidx'完全一致
  - 步骤2保证lidx‘和gidx/gidx'的误差不超过lidx和gidx的误差


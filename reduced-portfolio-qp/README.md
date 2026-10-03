# reduced-portfolio-qp

针对 8 个策略的投资组合凸二次规划求解器，采用 C++17 和原对偶内点法。源自 V5 降维实现，不依赖第三方数值求解库，目前仅支持 Windows。

16 个变量包含 8 个组合变量和 8 个换手辅助变量。利用辅助变量块的对角加秩一结构，通过 Sherman–Morrison 消元和 Schur 补，将 Newton 矩阵的 Cholesky 分解降到 8 维。遇到消减误差或分解不稳定时回退到 16 维。

## 最小文件分类

```text
reduced-portfolio-qp/
├── README.md                 使用说明
├── LICENSE                   项目 MIT 许可证
├── THIRD_PARTY_NOTICES.md     原源码附带的编译运行库许可声明
├── .gitignore                排除构建、结果和私有输入
├── src/portfolio.cpp         求解器、回测和 XLS 读取源码
├── scripts/run.ps1           编译、运行、基准和自检入口
├── scripts/make_example.py   无第三方依赖的合成数据生成器
├── examples/synthetic.xls    合成示例，可公开复现
└── docs/algorithm.md         算法与输入边界
```

## 快速开始

需要 Windows、PowerShell 和 Zig C++ 编译器；原项目使用 Zig 0.13.0。在仓库根目录执行：

```powershell
.\scripts\run.ps1 -Action build -Compiler 'C:\tools\zig\zig.exe'
.\scripts\run.ps1 -Action check
.\scripts\run.ps1 -Action run
.\scripts\run.ps1 -Action bench -Count 50
```

可把 Zig 加入 PATH 后省略 `-Compiler`。生成程序位于 `build/portfolio.exe`，结果位于 `results/results.json`。
使用自己的输入：`.\scripts\run.ps1 -Action run -InputPath 'D:\data\data.xls'`。
可选 Python 3 用于重新生成示例：`python scripts/make_example.py`。

## 功能和边界

支持最大分散、最大收益和均值方差目标，各有等权和指数加权统计；回看 60 行，持有 15 行。约束包括非负权重、权重和为 1、L1 换手预算 0.1，以及四组权重范围 0.15–0.35。
维度和约束固定，不是通用 QP 库。自检覆盖未来数据隔离和最大回撤穷举对照，不代表对所有输入的收敛保证。整理时其他合成输入曾触发 `Newton matrix is not positive definite`；当前提供的是通过自检的固定示例，原 V5 数值算法保留，任意输入的稳健性仍需完善。

示例由固定随机种子合成，不包含原笔试数据、参考 algo.py、个人报告或历史结果；收益和计时仅用于演示。原文件仍保留在原项目目录。

项目源码采用 MIT 许可证；编译运行库的原附带声明见 THIRD_PARTY_NOTICES.md。

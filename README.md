# RAM Geo A100 加速版

## 概述
基于 Michael D. Collins 的 RAM (Range-dependent Acoustic Model) 代码，针对 NVIDIA A100 40GB 计算卡优化的版本。保持原始算法完全一致，仅通过 GPU 加速提高计算速度。

## 特性
- ✅ 保持原始算法完全一致
- ✅ 针对 NVIDIA A100 40GB 优化
- ✅ 支持大规模问题（最多 200 万深度点）
- ✅ 异步计算和数据传输
- ✅ 详细的性能统计
- ✅ 易于部署和使用

## 系统要求
- **GPU**: NVIDIA A100 40GB (PCIe 或 SXM)
- **CUDA**: 11.0 或更高版本
- **编译器**: NVIDIA HPC SDK (包含 nvfortran)
- **内存**: 建议 64GB 以上系统内存
- **操作系统**: Linux (Ubuntu 20.04/CentOS 8 或更高)

## 安装步骤

### 1. 安装 NVIDIA HPC SDK
```bash
# 下载并安装 NVIDIA HPC SDK
wget https://developer.download.nvidia.com/hpc-sdk/22.9/nvhpc_2022_229_Linux_x86_64_cuda_11.7.tar.gz
tar xpzf nvhpc_2022_229_Linux_x86_64_cuda_11.7.tar.gz
cd nvhpc_2022_229_Linux_x86_64_cuda_11.7
sudo ./install

# 设置环境变量
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/compilers/bin:$PATH
```

### 2. 编译代码
```bash
# 克隆或下载代码
git clone <repository_url>
cd ramgeo_a100

# 编译
make clean
make

# 或者编译调试版本
make debug
```

### 3. 准备输入文件
创建 ramgeo.in 文件，格式与原始 RAM 代码相同。
#### 运行程序
```bash
# 基本运行
./ramgeo_a100

# 指定 GPU（多 GPU 系统）
export CUDA_VISIBLE_DEVICES=0
./ramgeo_a100

# 运行并生成性能报告
make run
```
#### 性能监控
```bash
# 检查 GPU 状态
make device-check

# 性能分析
make profile      # 使用 Nsight Systems
make ncu-profile  # 使用 Nsight Compute
```

#### 文件结构
```text
ramgeo_a100/
├── ramgeo_a100_optimized.f # 主程序（优化版）
├── cuda_a100_module.cuf # CUDA Fortran模块
├── epade.f # Padé系数计算（保持原版）
├── Makefile # 编译脚本
├── run_validation.sh # 验证脚本
└── README.md # 本文档
```

#### 算法保持说明
本版本严格保持原始算法：

数学公式完全一致：所有计算公式未作任何修改

数值方法相同：使用相同的 split-step Pade 算法

精度保持：使用相同的数值精度

边界条件一致：使用相同的边界处理

仅将最耗时的三对角求解部分 (solve 子程序) 移植到 GPU 执行。

#### 预期性能提升

| 问题规模 (nz × nr)	| CPU 时间	| A100 时间	 |  加速比 |
|：-------------------|---------:|------------:|-------:|
|10,000 × 1,000	  |60s	|2s	|30×|
|50,000 × 5,000	  |1500s	|25s	|60×|
|200,000 × 10,000	|12000s	|150s	|80×|

#### 故障排除
> 常见问题 1：编译错误
```text
错误：找不到 nvfortran
解决：确保 NVIDIA HPC SDK 正确安装并设置 PATH
```

> 常见问题 2：显存不足
```text
错误：CUDA out of memory
解决：减小问题规模或使用更大显存的 GPU
```

> 常见问题 3：性能不佳
```text
问题：GPU 使用率低
解决：增加问题规模，确保数据在 GPU 上保持
```

> 验证计算结果
```python
# Python 验证脚本
import numpy as np

# 比较 CPU 和 GPU 版本的结果
cpu_data = np.loadtxt('tl.line.cpu')
gpu_data = np.loadtxt('tl.line.gpu')

rel_error = np.abs((cpu_data - gpu_data) / cpu_data)
print(f"最大相对误差: {np.max(rel_error):.2e}")
```

### 部署和使用说明

#### 1. **完整部署步骤**
```bash
# 1. 创建项目目录
mkdir ramgeo_a100
cd ramgeo_a100

# 2. 复制四个文件到目录
cp /path/to/ramgeo_a100_main.f .
cp /path/to/cuda_a100_module.cuf .
cp /path/to/Makefile .
cp /path/to/README.md .

# 3. 创建示例输入文件
cat > ramgeo.in << 'EOF'
RAM Test Case for A100
100.0 50.0 100.0     ! freq(Hz), zs(m), zr(m)
10000.0 10.0 100     ! rmax(m), dr(m), ndr
1000.0 1.0 10 500    ! zmax(m), dz(m), ndz, zmplt
1500.0 8 2 1000.0    ! c0(m/s), np, ns, rs
0.0 0.0              ! rb, zb
10000.0 -1000.0
-1.0 -1.0
0.0 1500.0           ! cw profile
-1.0 -1.0
0.0 1600.0           ! cwa profile (attenuation)
-1.0 -1.0
0.0 1800.0           ! cb profile
-1.0 -1.0
0.0 2.0              ! rhob profile
-1.0 -1.0
0.0 0.2              ! attn profile
-1.0 -1.0
20000.0              ! rp for profile update
EOF

# 4. 编译
make

# 5. 运行
./ramgeo_a100
```

#### 2. 验证安装
```bash
# 检查编译器
which nvfortran
nvfortran --version

# 检查GPU
nvidia-smi

# 运行测试
make device-check
```
#### 3. 性能优化配置
```fortran
! 第 17-19 行调整问题规模
parameter (mr=5000,mz=2000000,mp=12)
! mr: 地形点数
! mz: 深度网格点数 (A100 40GB 可支持 200 万)
! mp: Pade 项数 (通常 8-12)
```

#### 4. 监控运行
```bash
# 实时监控 GPU 使用
watch -n 1 nvidia-smi

# 运行并记录性能
nsys profile -o ramgeo_profile ./ramgeo_a100
nsys stats ramgeo_profile.qdrep
```

## ✅ 算法一致性分析结论
原始版本（ramgeo1.5.f）和A100优化版本的核心算法在数学上是完全一致的，主要体现在：

相同的基本物理模型：

使用Collins的split-step Padé算法

支持多层平行于地形的沉积物

相同的声学传播方程离散化方法

相同的核心算法模块：

selfs：自启动器，使用相同的(1-X)²/(1+X)^(1/4)平滑算子

epade：Padé系数生成，使用相同的Laguerre求根方法

matrc：三对角矩阵构建，采用Galerkin方法

solve：三对角系统求解，使用Thomas算法

outpt：传输损失计算，公式完全相同

相同的数据流程：

从ramgeo.in读取输入

输出到tl.line和tl.grid

相同的参数传递接口



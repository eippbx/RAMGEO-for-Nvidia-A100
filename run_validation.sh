#!/bin/bash

# RAMGEO A100 优化版本验证脚本

echo "=========================================="
echo "RAMGEO A100 GPU优化版本验证"
echo "=========================================="

# 1. 编译原版
echo "步骤1: 编译原版ramgeo..."
gfortran -o ramgeo_orig ramgeo1.5.f -ffixed-line-length-none
if [ $? -ne 0 ]; then
    echo "错误: 原版编译失败"
    exit 1
fi

# 2. 编译A100优化版
echo "步骤2: 编译A100优化版..."
make clean
make
if [ $? -ne 0 ]; then
    echo "错误: A100优化版编译失败"
    exit 1
fi

# 3. 运行小规模测试
echo "步骤3: 运行小规模测试..."
cat > test_small.in << EOF
* Test case
50.0 100.0 200.0
10.0 0.1 10
500.0 5.0 10 500.0
1500.0 5 2 5.0
0.0 100.0
5.0 200.0
10.0 300.0
-1.0 -1.0
* Sound speed profile
0.0 1500.0
100.0 1480.0
200.0 1470.0
300.0 1490.0
-1.0 -1.0
* Bottom sound speed
0.0 1600.0
100.0 1650.0
200.0 1700.0
300.0 1750.0
-1.0 -1.0
* Density
0.0 1.0
100.0 1.5
200.0 1.8
300.0 2.0
-1.0 -1.0
* Attenuation
0.0 0.1
100.0 0.2
200.0 0.3
300.0 0.5
-1.0 -1.0
EOF

echo "运行原版..."
./ramgeo_orig < test_small.in
mv tl.line tl.line.orig
mv tl.grid tl.grid.orig

echo "运行A100优化版..."
./ramgeo_a100 < test_small.in
mv tl.line tl.line.a100
mv tl.grid tl.grid.a100

# 4. 验证结果
echo "步骤4: 验证结果一致性..."
python3 -c "
import numpy as np

# 读取原版结果
orig = np.loadtxt('tl.line.orig')
a100 = np.loadtxt('tl.line.a100')

# 检查长度
if len(orig) != len(a100):
    print('错误: 结果长度不一致')
    print(f'原版: {len(orig)} 行')
    print(f'A100: {len(a100)} 行')
    exit(1)

# 计算相对误差
diff = np.abs(orig[:,1] - a100[:,1])
rel_diff = diff / np.abs(orig[:,1])
max_rel_diff = np.max(rel_diff[np.abs(orig[:,1]) > 1e-10])

print(f'最大相对误差: {max_rel_diff:.2e}')

if max_rel_diff < 1e-6:
    print('✓ 结果一致性验证通过')
else:
    print('✗ 结果一致性验证失败')
    print('差异较大的行:')
    for i in range(len(diff)):
        if rel_diff[i] > 1e-3 and np.abs(orig[i,1]) > 1e-10:
            print(f'  行 {i+1}: 原版={orig[i,1]:.6f}, A100={a100[i,1]:.6f}, 误差={rel_diff[i]:.2e}')
"

# 5. 性能测试（可选）
if [ "$1" == "perf" ]; then
    echo "步骤5: 性能测试..."
    cat > test_large.in << EOF
* Large test case
100.0 100.0 200.0
50.0 0.05 20
1000.0 2.0 20 1000.0
1500.0 10 2 10.0
0.0 100.0
10.0 150.0
20.0 200.0
30.0 250.0
40.0 300.0
50.0 350.0
-1.0 -1.0
* Sound speed profile
0.0 1500.0
200.0 1480.0
400.0 1470.0
600.0 1490.0
800.0 1510.0
1000.0 1520.0
-1.0 -1.0
* Bottom sound speed (same as above)
0.0 1600.0
200.0 1650.0
400.0 1700.0
600.0 1750.0
800.0 1800.0
1000.0 1850.0
-1.0 -1.0
* Density
0.0 1.0
200.0 1.5
400.0 1.8
600.0 2.0
800.0 2.2
1000.0 2.5
-1.0 -1.0
* Attenuation
0.0 0.1
200.0 0.2
400.0 0.3
600.0 0.5
800.0 0.7
1000.0 1.0
-1.0 -1.0
EOF

    echo "运行性能测试..."
    echo "原版性能:"
    time ./ramgeo_orig < test_large.in > /dev/null 2>&1
    
    echo "A100优化版性能:"
    time ./ramgeo_a100 < test_large.in > /dev/null 2>&1
fi

# 6. 清理
echo "步骤6: 清理临时文件..."
rm -f test_small.in test_large.in

echo "=========================================="
echo "验证完成"
echo "=========================================="
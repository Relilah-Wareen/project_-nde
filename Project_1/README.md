# Project1: 泊松方程数值求解

## 目录结构
- `src/`: 源代码目录，包含 `BVPSolver.cpp` 和 `ProblemDomain.cpp`。
- `test/`: 测试目录，包含测试代码和 `Makefile`。
- `report/`: 报告目录，包含 LaTeX 报告和 `Makefile`。
- `data/`: 数据目录，存储生成的 CSV 文件。
- `plot/`: 绘图目录，存储生成的热力图。

## 命令说明

### 总目录下的 `Makefile`
- `make build`: 编译所有测试代码。
- `make run`: 运行所有测试代码。
- `make plot`: 运行所有 Python 绘图脚本，生成热力图。
- `make clean`: 清理生成的可执行文件和报告文件。
- `make report`: 生成 LaTeX 报告。

### `test/` 目录下的 `Makefile`
- `make build`: 编译所有测试代码。
- `make run`: 运行所有测试代码。
- `make plot`: 运行所有 Python 绘图脚本，生成热力图。
- `make clean`: 清理生成的可执行文件。

### `report/` 目录下的 `Makefile`
- `make report`: 编译 LaTeX 报告，生成 PDF 文件。
- `make clean`: 清理生成的 LaTeX 临时文件。

## 使用步骤
1. **编译代码**: 在总目录下运行 `make build`。
2. **运行测试**: 在总目录下运行 `make run`。
3. **生成热力图**: 在总目录下运行 `make plot`。
4. **生成报告**: 在总目录下运行 `make report`。
5. **清理文件**: 在总目录下运行 `make clean`。

## 依赖
- C++ 编译器（支持 C++20）
- LAPACK 和 LAPACKE 库
- Python 3（用于绘图）
- XeLaTeX（用于生成报告）

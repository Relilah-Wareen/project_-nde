# 编译器
CXX = g++
# 编译选项
CXXFLAGS = -std=c++20 -I../src
# 链接库
LDLIBS = -llapacke -llapack

# 可执行文件列表
EXECUTABLES = neumann_test mixed_ndnd_test dirichlet_test circular_dirichlet_dirichlet circular_neumann_dirichlet circular_dirichlet_neumann func2test func3test

# Python 脚本文件列表（自动查找当前目录下的所有 .py 文件）
PYTHON_SCRIPTS = $(wildcard *.py)

# 默认目标：构建所有可执行文件
build: $(EXECUTABLES)

# 构建规则：使用通配符匹配所有 .cpp 文件并生成对应的可执行文件
%: %.cpp ../src/BVPSolver.cpp ../src/ProblemDomain.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< ../src/BVPSolver.cpp ../src/ProblemDomain.cpp $(LDLIBS)

# 运行所有测试
run: build
	@for test in $(EXECUTABLES); do \
		echo "Running $$test..."; \
		./$$test; \
	done

# 运行所有 Python 绘图脚本
plot:
	@for script in $(PYTHON_SCRIPTS); do \
		echo "Running $$script..."; \
		python3 $$script; \
	done

# 清理生成的可执行文件
clean:
	rm -f $(EXECUTABLES)

.PHONY: build run plot clean
SRC_DIR = ../src
SRC_FILES = $(SRC_DIR)/BoundCon.cpp $(SRC_DIR)/BVPSolver.cpp $(SRC_DIR)/Function.cpp $(SRC_DIR)/Operator.cpp $(SRC_DIR)/SparseMatrix.cpp
TEST_FILE = test_batch.cpp
OUT_FILE = test_batch
PLOT_SCRIPTS = plot_1d_solutions.py plot_2d_solutions.py plot_2d_irregular.py

all: $(OUT_FILE)

build: $(SRC_FILES) $(TEST_FILE)
	g++ $(SRC_FILES) $(TEST_FILE) -o $(OUT_FILE)

clean:
	rm -f $(OUT_FILE)

run: $(OUT_FILE)
	./$(OUT_FILE)

plot: $(PLOT_SCRIPTS)
	python3 plot_1d_solutions.py
	python3 plot_2d_solutions.py
	python3 plot_2d_irregular.py

.PHONY: all build clean run plot

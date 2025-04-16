#include "../src/BVPSolver.hpp"
#include "../src/Operator.hpp"
#include "../src/BoundCon.hpp"
#include "../src/Function.hpp"
#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <map>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
using namespace std;

// 测试函数定义
// 测试函数1D real u = exp(sin(x)+x)
// 测试函数2D real u = exp(sin(x)+y)
double test_func(double x) {
    double exp_part = exp(x + sin(x));
    return -exp_part * (pow(1 + cos(x), 2) - sin(x));
}

double f2D(double x, double y) {
    return -exp(y+sin(x)) + sin(x)*exp(y+sin(x)) - cos(x)*cos(x)*exp(y+sin(x));
}

double Sq_Dirichlet_bound(double x, double y) {
    return exp(y+sin(x));
}

double Sq_Neumann_bound1(double x, double y) { return -exp(y); }
double Sq_Neumann_bound2(double x, double y) { return exp(y+sin(1))*cos(1); }
double Sq_Neumann_bound3(double x, double y) { return -exp(sin(x)); }
double Sq_Neumann_bound4(double x, double y) { return exp(1+sin(x)); }

// 函数映射表
map<string, function<double(double)>> func1d_map = {
    {"test_func", test_func}
};

map<string, function<double(double, double)>> func2d_map = {
    {"f2D", f2D},
    {"Sq_Dirichlet_bound", Sq_Dirichlet_bound},
    {"Sq_Neumann_bound1", Sq_Neumann_bound1},
    {"Sq_Neumann_bound2", Sq_Neumann_bound2},
    {"Sq_Neumann_bound3", Sq_Neumann_bound3},
    {"Sq_Neumann_bound4", Sq_Neumann_bound4}
};

// 创建边界条件
BoundCon createBoundary(const json& config) {
    int dim = config["dim"];
    string boundary_type = config["boundary_type"];
    bool irregular = false;  // 默认不设置为irregular
    
    if (dim == 2) {
        if (config["irregular"] == true) {
            irregular = true;  // 如果配置了irregular，并且为true，则设置为irregular
        }
    }
    
    if (dim == 1) {
        if (boundary_type == "Dirichlet") {
            return BoundCon(1, BCType::Dirichlet, 
                          config["boundary_values"]["left"], 
                          config["boundary_values"]["right"]);
        }
        else if (boundary_type == "Neumann") {
            return BoundCon(1, BCType::Neumann,
                          config["boundary_values"]["left"],
                          config["boundary_values"]["right"]);
        }
        else if (boundary_type == "Mixed") {
            string mixed_type = config["mixed_type"];
            Mixed1DType m1d;
            if (mixed_type == "DD") m1d = Mixed1DType::DD;
            else if (mixed_type == "DN") m1d = Mixed1DType::DN;
            else if (mixed_type == "ND") m1d = Mixed1DType::ND;
            else if (mixed_type == "NN") m1d = Mixed1DType::NN;
            else throw runtime_error("Invalid 1D mixed type");
            
            return BoundCon(1, BCType::Mixed,
                          config["boundary_values"]["left"],
                          config["boundary_values"]["right"],
                          m1d);
        }
    }
    else if (dim == 2) {
        if (boundary_type == "Dirichlet") {
            Function uniform_func(func2d_map["Sq_Dirichlet_bound"]);
            BoundCon bc(2, BCType::Dirichlet, uniform_func);
            
            // 只有在需要irregular的情况下，才设置为不规则
            if (irregular) {
                bc.set_Irregular();  // 设置为不规则
            }
            
            return bc;
        }
        else if (boundary_type == "Neumann") {
            Function left(func2d_map[config["boundary_functions"]["left"]]);
            Function right(func2d_map[config["boundary_functions"]["right"]]);
            Function bottom(func2d_map[config["boundary_functions"]["bottom"]]);
            Function top(func2d_map[config["boundary_functions"]["top"]]);
            return BoundCon(2, BCType::Neumann, left, right, bottom, top);
        }
        else if (boundary_type == "Mixed") {
            string mixed_type = config["mixed_type"];
            Mixed2DType m2d;
            
            if (mixed_type == "DDDD") m2d = Mixed2DType::DDDD;
            else if (mixed_type == "DDDN") m2d = Mixed2DType::DDDN;
            else if (mixed_type == "DDND") m2d = Mixed2DType::DDND;
            else if (mixed_type == "DDNN") m2d = Mixed2DType::DDNN;
            else if (mixed_type == "DNDD") m2d = Mixed2DType::DNDD;
            else if (mixed_type == "DNDN") m2d = Mixed2DType::DNDN;
            else if (mixed_type == "DNND") m2d = Mixed2DType::DNND;
            else if (mixed_type == "DNNN") m2d = Mixed2DType::DNNN;
            else if (mixed_type == "NDDD") m2d = Mixed2DType::NDDD;
            else if (mixed_type == "NDDN") m2d = Mixed2DType::NDDN;
            else if (mixed_type == "NDND") m2d = Mixed2DType::NDND;
            else if (mixed_type == "NDNN") m2d = Mixed2DType::NDNN;
            else if (mixed_type == "NNDD") m2d = Mixed2DType::NNDD;
            else if (mixed_type == "NNDN") m2d = Mixed2DType::NNDN;
            else if (mixed_type == "NNND") m2d = Mixed2DType::NNND;
            else if (mixed_type == "NNNN") m2d = Mixed2DType::NNNN;
            else throw runtime_error("Invalid 2D mixed type");
            
            Function left(func2d_map[config["boundary_functions"]["left"]]);
            Function right(func2d_map[config["boundary_functions"]["right"]]);
            Function bottom(func2d_map[config["boundary_functions"]["bottom"]]);
            Function top(func2d_map[config["boundary_functions"]["top"]]);
            
            return BoundCon(2, BCType::Mixed, m2d, left, right, bottom, top);
        }
    }
    
    throw runtime_error("Unsupported boundary condition configuration");
}


// 创建算子
template <int Dim>
Operator<Dim>* createOperator(const string& op_type) {
    if (op_type == "Injection") return new Injection<Dim>();
    if (op_type == "LinearInterpolation") return new LinearInterpolation<Dim>();
    if (op_type == "FullWeighting") return new FullWeighting<Dim>();
    if (op_type == "QuadraticInterpolation") return new QuadraticInterpolation<Dim>();
    throw runtime_error("Unsupported operator type: " + op_type);
}

// 运行单个测试用例
void runTestCase(const json& config, ofstream& logFile) {
    int dim = config["dim"];
    double epsilon = config["epsilon"];
    int nu1 = config["nu1"];
    int nu2 = config["nu2"];
    string solver_type = config["solver_type"];
    string restriction = config["restriction"];
    string prolongation = config["prolongation"];
    
    // 创建函数对象
    Function f;
    if (dim == 1) {
        f = Function(test_func);
    } else {
        f = Function(f2D);
    }
    
    // 创建边界条件
    BoundCon bc = createBoundary(config);

    for (int grid_num : config["grid_nums"]) {
        logFile << "Running test case: dim=" << dim 
             << ", grid_num=" << grid_num 
             << ", boundary=" << config["boundary_type"]
             << ", solver=" << solver_type << endl;
        
        if (dim == 1) {
            Operator<1>* restrict_op = createOperator<1>(restriction);
            Operator<1>* prolong_op = createOperator<1>(prolongation);
            
            BVPSolver<1> solver(*restrict_op, *prolong_op, grid_num, bc, f, nu1, nu2, epsilon);
            solver.BuildLinearSystem();
            
            auto start = chrono::high_resolution_clock::now();
            
            if (solver_type == "FMG") {
                solver.SolveByFMG();
            } else if (solver_type == "VC") {
                solver.SolveByVC();
            } else if (solver_type == "LU") {
                solver.SolveByLU();
            } else if (solver_type == "Relax") {
                solver.SolveByRelaxOnly(2.0/3);
            }
            
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
            
            string filename = "../data/solution_dim" + to_string(dim) + 
                            "_n" + to_string(grid_num) + 
                            "_" + config["boundary_type"].get<string>();
            if (bc.getBCType() == BCType::Mixed) {
                filename += "_" + string(config["mixed_type"]);
            }
            filename += ".csv";
            
            if(bc.getBCType()==BCType::Neumann) solver.normalizeNeumannSolution();
            solver.printResultToCSV(filename);
            
            logFile << "Completed in " << duration.count() << " μs" << endl;
            logFile << "Results saved to " << filename << endl << endl;
            
            delete restrict_op;
            delete prolong_op;
        }
        else if (dim == 2) {
            Operator<2>* restrict_op = createOperator<2>(restriction);
            Operator<2>* prolong_op = createOperator<2>(prolongation);
            
            BVPSolver<2> solver(*restrict_op, *prolong_op, grid_num, bc, f, nu1, nu2, epsilon);
            solver.BuildLinearSystem();
            
            auto start = chrono::high_resolution_clock::now();
            
            if (solver_type == "FMG") {
                solver.SolveByFMG();
            } else if (solver_type == "VC") {
                solver.SolveByVC();
            } else if (solver_type == "LU") {
                solver.SolveByLU();
            } else if (solver_type == "Relax") {
                solver.SolveByRelaxOnly(4.0/5);
            }
            
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
            
            string filename = "../data/solution_dim" + to_string(dim) + 
                            "_n" + to_string(grid_num) + 
                            "_" + config["boundary_type"].get<string>();
            if (bc.getBCType() == BCType::Mixed) {
                filename += "_" + string(config["mixed_type"]);
            }
            if (bc.IsIrregular()) {
                filename += "_irregular";
            }
            filename += ".csv";
            if(bc.getBCType()==BCType::Neumann) solver.normalizeNeumannSolution();
            solver.printResultToCSV(filename);
            
            logFile << "Completed in " << duration.count() << " μs" << endl;
            logFile << "Results saved to " << filename << endl << endl;
            
            delete restrict_op;
            delete prolong_op;
        }
    }
}

int main() {
    try {
        ofstream logFile("../run_log/test_log.txt");
        if (!logFile.is_open()) {
            throw runtime_error("Could not open test_log.txt for writing");
        }

        ifstream config_file("../json/config.json");
        if (!config_file.is_open()) {
            throw runtime_error("Could not open config.json");
        }
        
        json config = json::parse(config_file);
        
        for (const auto& test_case : config["test_cases"]) {
            runTestCase(test_case, logFile);
        }
        
        logFile.close(); // Ensure the file is properly closed when done
        
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}

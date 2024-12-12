#pragma once;
#include<optional>
#include<vector>

enum MODE
{
    TEST, MAIN
};


class Solver
{
private:

    std::vector<std::vector<double>> table;
    std::vector<double> num_solution;
    std::vector<double> double_solution;
    std::vector<double> real_solution;
    std::vector<double> differencies;

    double max_x_dif_pos{};
    double max_defference{};
    
    const double break_point = 0.4;
    const double m1 = 0.0;
    const double m2 = 1.0;
    int nodes_number{};
    double a{ 1.0 };
    double b{ 0.0 };

    double k1(double x, MODE mode);
    double k2(double x, MODE mode);
    double q1(double x, MODE mode);
    double q2(double x, MODE mode);
    double f1(double x, MODE mode);
    double f2(double x, MODE mode);


    double calc_a(double x, double step, MODE mode);
    double calc_d(double x, double step, MODE mode);
    double calc_phi(double x, double step, MODE mode);


    void Calc_double_solution(std::vector<double>& double_solution, int n);
    void Calc_real_solution(std::vector<double>& true_solution, int n);
    double Calc_differencies(const std::vector<double>& one, const std::vector<double>& other);

    void create_table(int n, MODE mode);


public:

    Solver();
    ~Solver();
    
    void Solve(int nodes_num, MODE mode);

    std::vector<std::vector<double>>& get_table();
    std::vector<double>& get_num_solution();
    std::vector<double>& get_real_solution();
    std::vector<double>& get_differencies();
    std::vector<double>& get_double_solution();
    double get_max_diff();
    double get_x_diff();
    double get_step();
    double get_break_point();
    std::vector<double> get_diaposon();
};


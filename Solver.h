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
    std::vector<double> real_solution;
    std::vector<double> differencies;

    double test_max_defference{};
    
    const double break_point = 0.4;
    const double m1 = 0.0;
    const double m2 = 1.0;
    int n{};


    double integral(const double high_border, const double low_border, double (*func)(std::optional<double>), std::optional<double> x = std::nullopt);
    static double k1(std::optional<double> x = std::nullopt);
    static double k2(std::optional<double> x = std::nullopt);
    static double q1(std::optional<double> x = std::nullopt);
    static double q2(std::optional<double> x = std::nullopt);
    static double f1(std::optional<double> x = std::nullopt);
    static double f2(std::optional<double> x = std::nullopt);


    void calc_a(std::vector<double>& a_coefficients,int n, MODE mode);
    void calc_d(std::vector<double>& d_coefficients, int n, MODE mode);
    void calc_phi(std::vector<double>& phi_coefficients, int n, MODE mode);

    void calc_test_a(std::vector<double>& a_coefficients, int n);
    void calc_main_a(std::vector<double>& a_coefficients, int n);
    void calc_test_d(std::vector<double>& d_coefficients, int n);
    void calc_main_d(std::vector<double>& d_coefficients, int n);
    void calc_test_phi(std::vector<double>& phi_coefficients, int n);
    void calc_main_phi(std::vector<double>& phi_coefficients, int n);


    void Calc_real_solution(std::vector<double>& true_solution, int n);
    double Calc_differencies(const std::vector<double>& one, const std::vector<double>& other);

    void create_table(int n);
    

public:

    Solver();
    ~Solver();

    void Solve(int nodes_num, MODE mode);

    std::vector<std::vector<double>>& get_table();
    std::vector<double>& get_num_solution();
    std::vector<double>& get_real_solution();
    std::vector<double>& get_differencies();
    double get_test_max_diff();
};


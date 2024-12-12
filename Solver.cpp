#include "Solver.h"
#include "SolveMatrix.h"
#include<optional>
#include<cmath>
#include<functional>
#include<vector>
#include<algorithm>

Solver::Solver() {};
Solver::~Solver()
{
};


double Solver::k1(double x, MODE mode)
{
	if (mode == MAIN) {
		return (x + 1);
	}
	else {
		return 1.4;
	}
}

double Solver::k2(double x, MODE mode)
{
	if (mode == MAIN) {
		return x;
	}
	else {
		return 0.4;
	}
}

double Solver::q1(double x, MODE mode)
{
	if (mode == MAIN) {
		return x;
	}
	else {
		return 0.4;
	}
}

double Solver::q2(double x, MODE mode)
{
	if (mode == MAIN) {
		return (x*x);
	}
	else {
		return 0.16;
	}
}

double Solver::f1(double x, MODE mode)
{
	if (mode == MAIN) {
		return x;
	}
	else {
		return 0.4;
	}
}

double Solver::f2(double x, MODE mode)
{
	if (mode == MAIN) {
		return exp(-x);
	}
	else {
		return exp(-0.4);
	}
}



double Solver::calc_a(double x, double step, MODE mode)
{
	if (mode == TEST) {
		double x_cur = x;
		double x_prev = x - step;
		double bp = break_point;

		if (bp >= x_cur) {
			return k1(x_cur - step / 2.0, mode);
		}
		else if (bp <= x_prev) {
			return k2(x_prev - step / 2.0, mode);
		}
		else {
			double tmp1 = (bp - x_prev) / k1((x_prev + bp) / 2.0, mode);
			double tmp2 = (x_cur - bp) / k2((bp + x_cur) / 2.0, mode);
			return 1.0 / ((1.0 / step) * (tmp1 + tmp2));
		}
	}
	else {
		double xi = x;
		double xi_1 = x - step;
		double KSI = break_point;

		if (KSI >= xi)
			return k1(xi - step / 2.0, mode);
		else if (KSI <= xi_1)
			return k2(xi - step / 2.0, mode);
		else
			return 1.0/((1 / step) * (((KSI - xi_1) / (k1((xi_1 + KSI) / 2.0, mode))) + ((xi - KSI) / (k2((KSI + xi) / 2.0, mode)))));
	}
}

double Solver::calc_d(double x, double step,  MODE mode)
{
	double xi_next = x + step / 2.0;
	double xi_prev = x - step / 2.0;
	double bp = break_point;

	if (bp >= xi_next) {
		return q1(x, mode);
	}
	else if (bp <= xi_prev) {
		return q2(x, mode);
	}
	else {
		double term1 = q1((xi_prev + bp) / 2.0, mode) * (bp - xi_prev);
		double term2 = q2((bp + xi_next) / 2.0, mode) * (xi_next - bp);
		return (1.0 / step) * (term1 + term2);
	}
}

double Solver::calc_phi(double x, double step,  MODE mode)
{
	double xi_next = x + step / 2.0;
	double xi_prev= x - step / 2.0;
	double bp = break_point;

	if (bp >= xi_next) {
		return f1(x, mode);
	}
	else if (bp <= xi_prev) {
		return f2(x, mode);
	}
	else {
		double term1 = f1((xi_prev + bp) / 2.0, mode) * (bp - xi_prev);
		double term2 = f2((bp + xi_next) / 2.0, mode) * (xi_next - bp);
		return (1.0 / step) * (term1 + term2);
	}
}

void Solver::Calc_double_solution(std::vector<double>& double_solution, int nodes_num)
{
	int n = 2*nodes_num;
	double step = fabs(a - b) / static_cast<double>(n);
	++n;

	std::vector<std::vector<double>> matrix_A;
	std::vector<double> vector_b;

	matrix_A.push_back({ 1.0, 0.0, 0.0 });
	vector_b.push_back(m1);

	for (int i = 1; i < n - 1; ++i) {
		double xi = i * step;
		double ai_prev = calc_a(xi, step, MAIN);
		double ai_next = calc_a(xi + step, step, MAIN);
		double di = calc_d(xi, step, MAIN);

		matrix_A.push_back({
			ai_prev / (step * step),
			-(ai_prev + ai_next) / (step * step) - di,
			ai_next / (step * step)
			});

		vector_b.push_back(-calc_phi(xi, step, MAIN));
	}

	matrix_A.push_back({ 0.0, 0.0, 1.0 });
	vector_b.push_back(m2);

	std::vector<double> tmp = SolveMatrix(matrix_A, vector_b);

	for (int i = 0; i < tmp.size(); ++i) {
		if (i % 2 == 0)
			double_solution.push_back(tmp[i]);
	}
}

void Solver::Calc_real_solution(std::vector<double>& true_solution, int n)
{
	double c1 = 0.060557222866650585;
	double c2 = -1.0605572228666509;
	double c3 = -0.4720245507344367;
	double c4 = -4.331084823580059;

	double step = 1.0 / static_cast<double>(n);
	double x = step;

	true_solution.emplace_back(m1);
	while (x <= break_point) {
		true_solution.emplace_back(c1 * exp(x * sqrt(14) / 7) + c2 * exp(-x * sqrt(14) / 7) + 1.0);
		x += step;
	}

	while (x <= (m2 + 0.00000000000001) ) {
		true_solution.emplace_back(c3 * exp(x * sqrt(10) / 5) + c4 * exp(-x * sqrt(10) / 5) + 6.25 * exp(-0.4));
		x += step;
	}
}

double Solver::Calc_differencies(const std::vector<double>& one, const std::vector<double>& other)
{
	if (one.size() != other.size()) {
		return -1.0;
	}
	else {
		double maximum = -1.0;

		differencies.resize(one.size());

		for (int i = 0; i < differencies.size(); ++i) {
			differencies[i] = abs(one[i] - other[i]);

			if (differencies[i] > maximum) {
				maximum = differencies[i];
				max_x_dif_pos = i;
			}
		}

		return maximum;
	}
}

void Solver::create_table(int n, MODE mode)
{
	std::vector<double> sol;
	if (mode == TEST) {
		sol = real_solution;
	}
	if (mode == MAIN) {
		sol = double_solution;
	}
	table.resize(num_solution.size());

	double x = m1;
	double step = fabs(a - b) / static_cast<double>(n);

	for (int i = 0; i < table.size(); ++i) {
		table[i] = { static_cast<double>(i), x, sol[i], num_solution[i], differencies[i] };
		x += step;
	}
}

void Solver::Solve(int nodes_num, MODE mode )
{
	nodes_number = nodes_num;
	int n = nodes_num;
	double step = fabs(a - b) / static_cast<double>(n);
	++n;

	std::vector<std::vector<double>> matrix_A;
	std::vector<double> vector_b;

	matrix_A.push_back({ 1.0, 0.0, 0.0 });
	vector_b.push_back(m1);

	for (int i = 1; i < n - 1; ++i) {
		double xi = i * step;
		double ai_prev = calc_a(xi, step, mode);
		double ai_next = calc_a(xi + step, step, mode);
		double di = calc_d(xi, step, mode);

		matrix_A.push_back({
			ai_prev / (step * step),
			-(ai_prev + ai_next) / (step * step) - di,
			ai_next / (step * step)
			});

		vector_b.push_back(-calc_phi(xi, step, mode));
	}

	matrix_A.push_back({ 0.0, 0.0, 1.0 });
	vector_b.push_back(m2);

	
	if (mode == TEST) {
		num_solution = SolveMatrix(matrix_A, vector_b);
		Calc_real_solution(real_solution, n - 1);
		max_defference = Calc_differencies(num_solution, real_solution);
	}
	if (mode == MAIN) {
		num_solution = SolveMatrix(matrix_A, vector_b);
		Calc_double_solution(double_solution, n - 1);
		max_defference = Calc_differencies(num_solution,double_solution);
	}

	create_table(n - 1,mode);
}



std::vector<std::vector<double>>& Solver::get_table()
{
	return table;
}

std::vector<double>& Solver::get_num_solution()
{
	return num_solution;
}

std::vector<double>& Solver::get_real_solution()
{
	return real_solution;
}

std::vector<double>& Solver::get_differencies()
{
	return differencies;
}

std::vector<double>& Solver::get_double_solution()
{
	return double_solution;
}

double Solver::get_max_diff()
{
	return max_defference;
}

double Solver::get_x_diff()
{
	return table[max_x_dif_pos][1];;
}

double Solver::get_step()
{
	return fabs(a- b)/static_cast<double>(nodes_number);
}

double Solver::get_break_point()
{
	return break_point;
}

std::vector<double> Solver::get_diaposon()
{
	std::vector<double> d{ b , a };
	return d;
}

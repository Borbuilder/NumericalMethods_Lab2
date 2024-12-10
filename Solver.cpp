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


double Solver::integral(const double high_border, const double low_border, double (*func)(std::optional<double>), std::optional<double> x)
{
	if (!x.has_value()) {
		if (func == k1 || func == k2) {
			return 1 / func(x)* (high_border - low_border);
		}
		else {
			return func(x) * (high_border - low_border);
		}
	}
	else {
		if (func == k1 || func == k2) {
			return 1/func((high_border + low_border) / 2) * (high_border - low_border);
		}
		else{
			return func((high_border + low_border) / 2) * (high_border - low_border);
		}
	}
}

double Solver::k1(std::optional<double> x)
{
	if (x.has_value()) {
		return *x + 1;
	}
	else {
		return 1.4;
	}
}

double Solver::k2(std::optional<double> x)
{
	if (x.has_value()) {
		return *x;
	}
	else {
		return 0.4;
	}
}

double Solver::q1(std::optional<double> x)
{
	if (x.has_value()) {
		return *x;
	}
	else {
		return 0.4;
	}
}

double Solver::q2(std::optional<double> x)
{
	if (x.has_value()) {
		return (*x)*(*x);
	}
	else {
		return 0.16;
	}
}

double Solver::f1(std::optional<double> x)
{
	if (x.has_value()) {
		return *x;
	}
	else {
		return 0.4;
	}
}

double Solver::f2(std::optional<double> x)
{
	if (x.has_value()) {
		return exp(-(*x));
	}
	else {
		return exp(-0.4);
	}
}



void Solver::calc_a(std::vector<double>& a_coefficients, int n, MODE mode)
{
	if (mode == TEST) 
		calc_test_a(a_coefficients,n);
	else
		calc_main_a(a_coefficients,n);
}

void Solver::calc_d(std::vector<double>& d_coefficients, int n, MODE mode)
{
	if (mode == TEST)
		calc_test_d(d_coefficients, n);
	else
		calc_main_d(d_coefficients, n);
}

void Solver::calc_phi(std::vector<double>& phi_coefficients, int n, MODE mode)
{
	if (mode == TEST)
		calc_test_phi(phi_coefficients, n);
	else
		calc_main_phi(phi_coefficients, n);
}

void Solver::calc_test_a(std::vector<double>& a_coefficients, int n)
{
	double step = 1.0 / static_cast<double>(n);
	double x = step;

	while (x <= break_point) {
		a_coefficients.emplace_back( 1/((1/step)*integral(x, x - step, k1)) );
		x += step;
	}

	if (x - step < break_point && break_point < x) 
		a_coefficients.emplace_back( 1/((1/step)*(integral(break_point, x - step, k1) + integral(x, break_point, k2))) );

	while (x <= m2) {
		a_coefficients.emplace_back(1 / ((1 / step) * integral(x, x - step, k2)));
		x += step;
	}
}

void Solver::calc_main_a(std::vector<double>& a, int n)
{
	const double step = 1.0 / static_cast<double>(n);
	double x = step;
	double x_half = step / 2.0;

	while (x <= break_point)
	{
		a.push_back(k1(x_half));
		x += step;
		x_half += step;
	}

	a.push_back(step / ((break_point - x + step) / k1((break_point + x - step) / 2.0) + (x - break_point) / k2((x + break_point) / 2.0)));
	x += step;
	x_half += step;

	while (x < 1.0)
	{
		a.push_back(k2(x_half));
		x += step;
		x_half += step;
	}
	if (std::abs(x - 1.0) < 1e-8)
	{
		a.push_back(k2(x_half));
	}
}

void Solver::calc_test_d(std::vector<double>& d_coefficients, int n)
{
	double step = 1.0 / static_cast<double>(n);
	double x = step;

	while (x + 0.5 * step <= break_point) {
		d_coefficients.emplace_back( (1 / step) * integral(x + 0.5 * step, x - 0.5 * step, q1) );
		x += step;
	}

	if (x - 0.5 * step < break_point && break_point < x + 0.5 * step)
		d_coefficients.emplace_back( (1 / step) * (integral(break_point, x - 0.5 * step, q1) + integral(x + 0.5 * step, break_point, q2)) );
	x += step;

	while (x <= m2 - step) {
		d_coefficients.emplace_back( (1 / step) * integral(x + 0.5 * step, x - 0.5 * step, q2) );
		x += step;
	}
}

void Solver::calc_main_d(std::vector<double>& d, int n)
{
	const double step = 1.0 / static_cast<double>(n);
	double x_half = step * 3.0 / 2.0;
	double x = step;
	while (x_half <= break_point)
	{
		d.push_back(q1(x));
		x += step;
		x_half += step;
	}
	d.push_back((q1((x_half - step + break_point) / 2.0) * (break_point - x_half + step) + q2((x_half + break_point) / 2.0) * (x_half - break_point)) / step);
	x += step;
	x_half += step;
	while (x_half < 1.0)
	{
		d.push_back(q2(x));
		x += step;
		x_half += step;
	}
}

void Solver::calc_test_phi(std::vector<double>& phi_coefficients, int n)
{
	double step = 1.0 / static_cast<double>(n);
	double x = step;

	while (x + 0.5 * step <= break_point) {
		phi_coefficients.emplace_back((1 / step) * integral(x + 0.5 * step, x - 0.5 * step, f1));
		x += step;
	}

	if (x - 0.5 * step < break_point && break_point < x + 0.5 * step)
		phi_coefficients.emplace_back((1 / step) * (integral(break_point, x - 0.5 * step, f1) + integral(x + 0.5 * step, break_point, f2)));
	x += step;

	while (x <= m2-step) {
		phi_coefficients.emplace_back((1 / step) * integral(x + 0.5 * step, x - 0.5 * step, f2));
		x += step;
	}
}

void Solver::calc_main_phi(std::vector<double>& phi, int n)
{
	const double step = 1.0 / static_cast<double>(n);
	double x_half = step * 3.0 / 2.0;
	double x = step;

	while (x_half <= break_point)
	{
		phi.push_back(q1(x));
		x += step;
		x_half += step;
	}
	phi.push_back((f1((x_half - step + break_point) / 2.0) * (break_point - x_half + step) + f2((x_half + break_point) / 2.0) * (x_half - break_point)) / step);
	x += step;
	x_half += step;
	while (x_half < 1.0)
	{
		phi.push_back(f2(x));
		x += step;
		x_half += step;
	}
}



void Solver::Calc_double_solution(std::vector<double>& double_solution, int n)
{
	std::vector<double> a_coefficients;
	calc_a(a_coefficients, n, MAIN);
	std::vector<double> d_coefficients;
	calc_d(d_coefficients, n, MAIN);
	std::vector<double> phi_coefficients;
	calc_phi(phi_coefficients, n, MAIN);

	const unsigned long long size = d_coefficients.size();
	double step = 1.0 / static_cast<double>(2 * n);

	//const unsigned long long size = d_coefficients.size();

	std::vector<double> A(size);
	std::vector<double> B(size);
	std::vector<double> C(size);

	for (unsigned long long i = 0; i < size; i++)
	{
		A[i] = a_coefficients[i] / (step * step);
		B[i] = a_coefficients[i + 1] / (step * step);
		C[i] = (a_coefficients[i] + a_coefficients[i + 1]) / (step * step) + d_coefficients[i];
	}

	double_solution = SolveMatrix(A, C, B, phi_coefficients, { 0.0, 0.0 }, { m1, m2 });

	/*double_solution = std::vector<double>(n + 1);
	for (unsigned long long i = 0; i < n + 1; i++)
	{
		double_solution[i] = tmp_sol[i * 2];
	}*/
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

	while (x < m2) {
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
	double step = 1.0 / static_cast<double>(n);

	for (int i = 0; i < table.size(); ++i) {
		table[i] = { static_cast<double>(i), x, sol[i], num_solution[i], differencies[i] };
		x += step;
	}
}

void Solver::Solve(int nodes_num, MODE mode)
{
	int n = nodes_num;

	std::vector<double> a_coefficients;
	calc_a(a_coefficients, n, mode);
	std::vector<double> d_coefficients;
	calc_d(d_coefficients, n, mode);
	std::vector<double> phi_coefficients;
	calc_phi(phi_coefficients, n, mode);

	double step = 1 / static_cast<double>(n);
	const unsigned long long size = d_coefficients.size();

	std::vector<double> A(size);
	std::vector<double> B(size);
	std::vector<double> C(size);

	for (unsigned long long i = 0; i < size; i++)
	{
		A[i] = a_coefficients[i] / (step * step);
		B[i] = a_coefficients[i + 1] / (step * step);
		C[i] = (a_coefficients[i] + a_coefficients[i + 1]) / (step * step) + d_coefficients[i];
	}

	if (mode == TEST) {
		num_solution = SolveMatrix(A, C, B, phi_coefficients, { 0.0, 0.0 }, { m1, m2 });
		Calc_real_solution(real_solution, n);
		max_defference = Calc_differencies(num_solution, real_solution);
	}
	if (mode == MAIN) {
		num_solution = SolveMatrix(A, C, B, phi_coefficients, { 0.0, 0.0 }, { m1, m2 });
		Calc_double_solution(double_solution, n);
		max_defference = Calc_differencies(num_solution,double_solution);
	}

	create_table(n,mode);
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

double Solver::get_max_diff()
{
	return max_defference;
}

double Solver::get_x_diff()
{
	return table[max_x_dif_pos][1];;
}

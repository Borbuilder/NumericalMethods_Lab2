#include "Solver.h"
#include<iostream>
#include<iomanip>

int main()
{
	Solver solution1;
/*
	solution.Solve(15, TEST);
	std::vector<std::vector<double>> table = solution.get_table();

	for (auto& vec : table) {
		for (auto& el : vec) {
			std::cout << std::setprecision(15)<< el << "  ";
		}
		std::cout << std::endl;
	}

	std::cout << solution.get_test_max_diff();

	Solver solution1;*/

	solution1.Solve(1000, MAIN);
	std::vector<std::vector<double>> table1 = solution1.get_table();

	for (auto& vec : table1) {
		for (auto& el : vec) {
			std::cout << std::setprecision(15) << el << "  ";
		}
		std::cout << std::endl;
	}

	std::cout << solution1.get_max_diff();
	return 0;
}
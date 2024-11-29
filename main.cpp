#include "Solver.h"
#include<iostream>
#include<iomanip>

int main()
{
	Solver solution;

	solution.Solve(19, TEST);
	std::vector<std::vector<double>> table = solution.get_table();

	for (auto& vec : table) {
		for (auto& el : vec) {
			std::cout << std::setprecision(15)<< el << "  ";
		}
		std::cout << std::endl;
	}

	std::cout << solution.get_test_max_diff();
	return 0;
}
#include "Solver.h"
#include<iostream>
#include<iomanip>

int main()
{
	
	Solver solution1;

	solution1.Solve(1000, TEST);
	auto table = solution1.get_table();
	for (auto& vec : table) {
		for (auto& el : vec) {
			std::cout << std::setprecision(15) << el << "  ";
		}
		std::cout << std::endl;
	}

	std::cout << solution1.get_max_diff() << std::endl;
	std::cout << solution1.get_x_diff() << std::endl;
	std::cout << solution1.get_step() << std::endl;
	std::cout << solution1.get_break_point() << std::endl;
	std::cout << solution1.get_diaposon()[0] << " , " << solution1.get_diaposon()[1];
	return 0;
}

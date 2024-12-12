#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "SolveMatrix.h"
#include "Solver.h"

namespace py = pybind11;

PYBIND11_MODULE(Lab_2_Module, m) {
    py::class_<Solver>(m, "Solver")
        .def(py::init<>())
        .def("Solve", &Solver::Solve)
        .def("get_table", &Solver::get_table)
        .def("get_num_solution", &Solver::get_num_solution)
        .def("get_real_solution", &Solver::get_real_solution)
        .def("get_differencies", &Solver::get_differencies)
        .def("get_double_solution", &Solver::get_double_solution)
        .def("get_max_diff", &Solver::get_max_diff)
        .def("get_x_diff", &Solver::get_x_diff)
        .def("get_step", &Solver::get_step)
        .def("get_break_point", &Solver::get_break_point)
        .def("get_diaposon", &Solver::get_diaposon);

    py::enum_<MODE>(m, "MODE")
        .value("TEST", MODE::TEST)
        .value("MAIN", MODE::MAIN)
        .export_values();
}

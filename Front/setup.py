from setuptools import setup, Extension
import pybind11

pybind11_include = pybind11.get_include()

ext_modules = [
    Extension(
        'Lab_2_Module',
        ['TaskModule.cpp', 'SolveMatrix.cpp', 'Solver.cpp'],
        include_dirs=[pybind11_include],
        extra_compile_args=['/std:c++17'],
        extra_link_args=[],
    ),
]

setup(
    name='Lab_2_Module',
    version='0.1',
    ext_modules=ext_modules,
    zip_safe=False,
)

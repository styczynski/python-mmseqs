#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "unafold/hybrid-min.hpp"

namespace py = pybind11;

/*
struct hybrid_min_args {
    std::string input1;
    std::string input2;
    std::optional<double> tMin;
    std::optional<double> tMax;
    std::optional<double> tInc;
    bool version;
    bool help;
    bool debug;
    std::optional<std::string> na;
    std::optional<std::string> output;
    std::optional<std::string> suffix;
    std::optional<double> sodium;
    std::optional<double> magnesium;
    bool polymer;
    bool quiet;
    std::optional<std::string> mfold;
    std::optional<std::string> constraints;
};
*/

PYBIND11_PLUGIN(unafold_python_native) {
  py::module m("unafold_python_native", R"doc(
        Python module
        -----------------------
        .. currentmodule:: unafold_python
        .. autosummary::
           :toctree: _generate
           
           add
           subtract
    )doc");

  pybind11::class_<hybrid_min_args>(m, "HybridMinArgs")
    .def(pybind11::init<>())
    .def_readwrite("input1", &hybrid_min_args::input1)
    .def_readwrite("input2", &hybrid_min_args::input2)
    .def_readwrite("tmin", &hybrid_min_args::tMin)
    .def_readwrite("tmax", &hybrid_min_args::tMax)
    .def_readwrite("tinc", &hybrid_min_args::tInc)
    .def_readwrite("version", &hybrid_min_args::version)
    .def_readwrite("help", &hybrid_min_args::help)
    .def_readwrite("debug", &hybrid_min_args::debug)
    .def_readwrite("na", &hybrid_min_args::na)
    .def_readwrite("output", &hybrid_min_args::output)
    .def_readwrite("suffix", &hybrid_min_args::suffix)
    .def_readwrite("sodium", &hybrid_min_args::sodium)
    .def_readwrite("magnesium", &hybrid_min_args::magnesium)
    .def_readwrite("polymer", &hybrid_min_args::polymer)
    .def_readwrite("quiet", &hybrid_min_args::quiet)
    .def_readwrite("mfold", &hybrid_min_args::mfold)
    .def_readwrite("constraints", &hybrid_min_args::constraints);

  m.def("_hybrid_min", &unafold_hybrid_min, R"doc(
        Run hybrid-min utility
    )doc");

  return m.ptr();
}
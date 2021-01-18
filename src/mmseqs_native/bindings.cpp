#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "src/src/commons/Application.h"

namespace py = pybind11;

PYBIND11_PLUGIN(unafold_python_native) {
  py::module m("mmseqs_native", R"doc(
        Python module
        -----------------------
        .. currentmodule:: unafold_python
        .. autosummary::
           :toctree: _generate
           
           add
           subtract
    )doc");

    pybind11::class_<mmseqs_call_args>(m, "MMSeqsCallArgs")
    .def(pybind11::init<>())
    .def_readwrite("cli_args", &mmseqs_call_args::cli_args);

  m.def("_call_mmseqs", &call_mmseqs, R"doc(
        Run mmseqs2
    )doc");

  return m.ptr();
}
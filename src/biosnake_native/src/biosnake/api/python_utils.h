#ifndef API_PYTHON_UTILS_H
#define API_PYTHON_UTILS_H

#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

static py::array sequence_to_numpy(const char* data, const size_t len) {
    auto dtype = pybind11::dtype(pybind11::format_descriptor<char>::format());
    auto base = pybind11::array(dtype, {len}, {sizeof(char)});
    return pybind11::array(dtype, {len}, {sizeof(int)}, data, base);
}

template<typename T>
py::array vector_to_numpy(const std::vector<T>& data) {
    auto v = new std::vector<T>(data);
    auto capsule = py::capsule(v, [](void *v) { delete reinterpret_cast<std::vector<T>*>(v); });
    return py::array(v->size(), v->data(), capsule);
}

#endif
#include <pybind11/pybind11.h>
namespace py = pybind11;

void _native_summary(py::module_);

PYBIND11_MODULE(_native, m) {
    _native_summary(m.def_submodule("summary"));
}

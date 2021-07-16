#include <cstdef>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

namespace py = pybind11;
using std::uintptr_t;

class CWrap {
    void *m_ptr;
    bool m_is_reference;

public:
    CWrap(uintptr_t ptr, py::object parent, bool is_reference)
        : m_ptr(reinterpret_cast<void *>(ptr)), m_is_reference(is_reference) {}

    uintptr_t address() const { return reinterpret_cast<uintptr_t>(m_ptr); }

    operator bool() const { return m_ptr; }
};

void cwrap_support(py::class_ c) {
    c.def("from_param", [](py::object cls, py::object other));
}

#include <pybind11/detail/common.h>
#include <string>
#include <vector>
#include <optional>
#include <chrono>
#include <iostream>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/chrono.h>
#include <pybind11/stl.h>

#include "ert/ecl/ecl_sum.hpp"
#include "ert/ecl/ecl_sum_data.hpp"
#include "ert/ecl/smspec_node.hpp"

#include "detail/ecl/ecl_sum.hpp"
#include "detail/ecl/ecl_sum_data.hpp"
#include "detail/ecl/ecl_sum_file_data.hpp"
#include "detail/ecl/ecl_smspec.hpp"

extern "C" double ecl_sum_data_vector_iget(const ecl_sum_data_type *data,
                                           time_t sim_time, int params_index,
                                           bool is_rate, int time_index1,
                                           int time_index2, double weight1,
                                           double weight2);

void ecl_sum_data_init_interp_from_sim_time(const ecl_sum_data_type *data,
                                            time_t sim_time, int *index1,
                                            int *index2, double *weight1,
                                            double *weight2);

namespace py = pybind11;
using std::chrono::system_clock;

using namespace pybind11::literals;

using std::optional;
using std::string;
using std::vector;

namespace ecl {
namespace data {

void assert_ecl_type(const void *obj, int type_id) {
    if (!(*reinterpret_cast<const int *>(obj) == type_id)) {
        throw py::type_error("Unexpected ecl type");
    }
}

} // namespace data
} // namespace ecl

namespace {
template <typename T> const T *_get_object_from_cwrap(py::object &obj) {
    static auto basecclass = py::module_::import("cwrap").attr("BaseCClass");
    if (!py::isinstance(obj, basecclass)) {
        throw py::type_error("Expected a BaseCClass type");
    }

    auto address = py::cast<uint64_t>(obj.attr("_address")());
    auto object = reinterpret_cast<const T *>(address);
    ecl::data::assert_ecl_type(object, ECL_SUM_ID);
    return object;
}

const ecl::smspec_node &_get_smspec_by_name(const ecl_sum_type *sum,
                                            const std::string &key) {
    const auto &index = sum->smspec->gen_var_index;
    auto it = index.find(key.c_str());
    if (it == index.end())
        throw py::key_error("No such key");
    return *it->second;
}

const ecl::smspec_node &_get_node_by_params_index(const ecl_smspec_type *smspec,
                                                  int params_index) {
    auto node_index = smspec->inv_index_map.at(params_index);
    return *smspec->smspec_nodes[node_index];
}

vector<const ecl::smspec_node *>
_get_nodes_by_keys(const ecl_smspec_type *smspec,
                   const optional<vector<string>> &keys) {
    std::vector<const ecl::smspec_node *> nodes{};
    if (keys.has_value()) {
        for (const auto &key : *keys) {
            try {
                nodes.push_back(smspec->gen_var_index.at(key));
            } catch (std::out_of_range &e) {
                throw py::value_error(e.what());
            }
        }
    } else {
        for (const auto &node : smspec->smspec_nodes) {
            nodes.push_back(node.get());
        }
    }
    return nodes;
}

} // namespace

/**
 * @brief Generate a numpy 1d-array at the specified time points
 */
py::array to_numpy(py::object self, const std::string &key,
                   const vector<system_clock::time_point> &time_points) {
    const auto ecl_sum = _get_object_from_cwrap<ecl_sum_type>(self);
    const auto data = ecl_sum->data;
    const auto &smspec = _get_smspec_by_name(ecl_sum, key);

    time_t start_time = ecl_sum_data_get_data_start(data);
    time_t end_time = ecl_sum_data_get_sim_end(data);
    float start_value = 0;
    float end_value = 0;

    if (!smspec.is_rate()) {
        start_value =
            ecl_sum_data_iget_first_value(data, smspec.get_params_index());
        end_value =
            ecl_sum_data_iget_last_value(data, smspec.get_params_index());
    }

    py::array_t<float, py::array::c_style> array(time_points.size());
    auto output_data = array.mutable_data();
    for (size_t i{}; i < time_points.size(); ++i) {
        auto sim_time = system_clock::to_time_t(time_points[i]);

        float value{};
        if (sim_time < start_time) {
            value = start_value;
        } else if (sim_time > end_time) {
            value = end_value;
        } else {
            int time_index1, time_index2;
            double weight1, weight2;
            ecl_sum_data_init_interp_from_sim_time(
                data, sim_time, &time_index1, &time_index2, &weight1, &weight2);
            value = ecl_sum_data_vector_iget(
                data, sim_time, smspec.get_params_index(), smspec.is_rate(),
                time_index1, time_index2, weight1, weight2);
        }

        output_data[i] = value;
    }
    return std::move(array);
}

/**
 * @brief Generate a numpy array report_only
 */
py::array to_numpy_report_only(py::object self, const std::string &key) {
    const auto ecl_sum = _get_object_from_cwrap<ecl_sum_type>(self);
    const auto data = ecl_sum->data;
    const auto &smspec = _get_smspec_by_name(ecl_sum, key);

    size_t size = 1 + ecl_sum_data_get_last_report_step(data) -
                  ecl_sum_data_get_first_report_step(data);

    py::array_t<float, py::array::c_style> array(size);
    int offset = 0;
    auto output_data = array.mutable_data();
    for (const auto &index_node : data->index) {
        const auto &data_file = data->data_files[index_node.data_index];
        const auto &params_map = index_node.params_map;
        int params_index = params_map[smspec.get_params_index()];

        const ecl::smspec_node &smspec_node =
            _get_node_by_params_index(data->smspec, smspec.get_params_index());
        float default_value = smspec_node.get_default();
        offset +=
            data_file->get_data_report(params_index, index_node.length,
                                       &output_data[offset], default_value);
    }

    return std::move(array);
}

/**
 * @brief Generate a numpy array
 */
py::array to_numpy(py::object self, const std::string &key) {
    const auto ecl_sum = _get_object_from_cwrap<ecl_sum_type>(self);
    const auto data = ecl_sum->data;
    const auto &smspec = _get_smspec_by_name(ecl_sum, key);

    py::array_t<float, py::array::c_style> array(data->index.length());
    int offset = 0;
    auto output_data = array.mutable_data();
    for (const auto &index_node : data->index) {
        const auto &data_file = data->data_files[index_node.data_index];
        const auto &params_map = index_node.params_map;
        int params_index = params_map[smspec.get_params_index()];

        if (params_index >= 0)
            data_file->get_data(params_index, index_node.length,
                                &output_data[offset]);
        else {
            const ecl::smspec_node &smspec_node = _get_node_by_params_index(
                data->smspec, smspec.get_params_index());
            for (int i = 0; i < index_node.length; i++)
                output_data[offset + i] = smspec_node.get_default();
        }
        offset += index_node.length;
    }

    return std::move(array);
}

py::array to_dataframe_values(py::object self,
                              const optional<vector<string>> &keys) {
    const auto ecl_sum = _get_object_from_cwrap<ecl_sum_type>(self);
    const auto data = ecl_sum->data;
    const auto smspec = ecl_sum->smspec;

    // Find all relevant smspec_nodes
    const auto nodes = _get_nodes_by_keys(smspec, keys);

    const ssize_t num_rows = ecl_sum_data_get_length(data);
    const ssize_t num_cols = nodes.size();

    py::array_t<float, py::array::c_style> array{{num_rows, num_cols}};
    for (ssize_t time_index{}; time_index < num_rows; ++time_index) {
        for (ssize_t key_index{}; key_index < num_cols; ++key_index) {
            const auto node = nodes[key_index];
            const auto param_index = node->get_params_index();

            array.mutable_at(time_index, key_index) =
                ecl_sum_data_iget(data, time_index, param_index);
        }
    }
    return std::move(array);
}

py::array to_dataframe_values_interp(
    py::object self, const optional<vector<string>> &keys,
    const vector<system_clock::time_point> &time_points) {
    const auto ecl_sum = _get_object_from_cwrap<ecl_sum_type>(self);
    const auto data = ecl_sum->data;
    const auto smspec = ecl_sum->smspec;

    // Find all relevant smspec_nodes
    const auto nodes = _get_nodes_by_keys(smspec, keys);

    const auto start_time = ecl_sum_data_get_data_start(data);
    const auto end_time = ecl_sum_data_get_sim_end(data);

    const ssize_t num_rows = time_points.size();
    const ssize_t num_cols = nodes.size();

    py::array_t<float, py::array::c_style> array{{num_rows, num_cols}};

    auto inner_loop_clamp = [&](auto time_index, auto fn) {
        for (ssize_t key_index{}; key_index < num_cols; ++key_index) {
            const auto node = nodes[key_index];
            const auto param_index = node->get_params_index();

            float value = 0.f;
            if (!node->is_rate())
                value = fn(data, param_index);
            array.mutable_at(time_index, key_index) = value;
        }
    };

    auto inner_loop = [&](auto time_index, auto sim_time) {
        double weight1, weight2;
        int time_index1, time_index2;

        ecl_sum_data_init_interp_from_sim_time(
            data, sim_time, &time_index1, &time_index2, &weight1, &weight2);

        for (ssize_t key_index{}; key_index < num_cols; ++key_index) {
            const auto node = nodes[key_index];
            const auto param_index = node->get_params_index();

            float value = ecl_sum_data_vector_iget(
                data, sim_time, param_index, node->is_rate(), time_index1,
                time_index2, weight1, weight2);
            array.mutable_at(time_index, key_index) = value;
        }
    };

    for (ssize_t time_index{}; time_index < num_rows; ++time_index) {
        auto sim_time = system_clock::to_time_t(time_points[time_index]);
        if (sim_time < start_time) {
            inner_loop_clamp(time_index, &ecl_sum_data_iget_first_value);
        } else if (sim_time > end_time) {
            inner_loop_clamp(time_index, &ecl_sum_data_iget_last_value);
        } else {
            inner_loop(time_index, sim_time);
        }
    }
    return std::move(array);
}

void _native_summary(py::module_ m) {
    using time_index_t = optional<vector<system_clock::time_point>>;

    m.def(
        "to_numpy",
        [](py::object self, std::string &key, const time_index_t &time_index,
           bool report_only) -> py::array {
            if (time_index && report_only)
                throw std::invalid_argument(
                    "Either time_index=None or report_only=False");
            if (time_index)
                return to_numpy(self, key, *time_index);
            if (report_only)
                return to_numpy_report_only(self, key);
            else
                return to_numpy(self, key);
        },
        py::arg{"self"}, py::arg{"key"}, py::arg{"time_index"},
        py::arg{"report_only"});

    m.def("to_dataframe_values", &to_dataframe_values, py::arg{"self"},
          py::arg{"keys"});
    m.def("to_dataframe_values_interp", &to_dataframe_values_interp,
          py::arg{"self"}, py::arg{"keys"}, py::arg{"time_points"});
}

#pragma once

#include <vector>

typedef struct ecl_smspec_struct ecl_smspec_type;

namespace detail {

struct IndexNode {
    IndexNode(int d, int o, int l) {
        this->data_index = d;
        this->offset = o;
        this->length = l;
    }

    int end() const { return this->offset + this->length; }

    int data_index;
    int offset;
    int length;
    int report1;
    int report2;
    time_t time1;
    time_t time2;
    double days1;
    double days2;
    std::vector<int> params_map;
};

class CaseIndex {
public:
    IndexNode &add(int length) {
        int offset = 0;
        int data_index = this->index.size();

        if (!this->index.empty())
            offset = this->index.back().end();

        this->index.emplace_back(data_index, offset, length);
        return this->index.back();
    }

    /*
  The lookup_time() and lookup_report() methods will lookup which file_data
  instance corresponds to the time/report argument. The methods will return two
  pointers to file_data instances, if the argument is inside one file_data
  instance the pointers will be equal - otherwise they will point to the
  file_data instance before and after the argument:

  File 1                     File 2
  |------|-----|------|      |----|----------|---|
      /|\                /|\
       |                  |
       |                  |
       A                  B

  For time A the lookup_time function will return <file1,file1> whereas for time
  B the function will return <file1,file2>.
 */

    std::pair<const IndexNode *, const IndexNode *>
    lookup_time(time_t sim_time) const {
        auto iter = this->index.begin();
        auto next = this->index.begin();
        if (sim_time < iter->time1)
            throw std::invalid_argument("Simulation time out of range");

        ++next;
        while (true) {
            double t1 = iter->time1;
            double t2 = iter->time2;

            if (sim_time >= t1) {
                if (sim_time <= t2)
                    return std::make_pair<const IndexNode *, const IndexNode *>(
                        &(*iter), &(*iter));

                if (next == this->index.end())
                    throw std::invalid_argument("Simulation days out of range");

                if (sim_time < next->time1)
                    return std::make_pair<const IndexNode *, const IndexNode *>(
                        &(*iter), &(*next));
            }
            ++next;
            ++iter;
        }
    }

    std::pair<const IndexNode *, const IndexNode *>
    lookup_days(double days) const {
        auto iter = this->index.begin();
        auto next = this->index.begin();
        if (days < iter->days1)
            throw std::invalid_argument("Simulation days out of range");

        ++next;
        while (true) {
            double d1 = iter->days1;
            double d2 = iter->days2;

            if (days >= d1) {
                if (days <= d2)
                    return std::make_pair<const IndexNode *, const IndexNode *>(
                        &(*iter), &(*iter));

                if (next == this->index.end())
                    throw std::invalid_argument("Simulation days out of range");

                if (days < next->days1)
                    return std::make_pair<const IndexNode *, const IndexNode *>(
                        &(*iter), &(*next));
            }
            ++next;
            ++iter;
        }
    }

    const IndexNode &lookup(int internal_index) const {
        for (const auto &node : this->index)
            if (internal_index >= node.offset && internal_index < node.end())
                return node;

        throw std::invalid_argument("Internal error when looking up index: " +
                                    std::to_string(internal_index));
    }

    const IndexNode &lookup_report(int report) const {
        for (const auto &node : this->index)
            if (node.report1 <= report && node.report2 >= report)
                return node;

        throw std::invalid_argument("Internal error when looking up report: " +
                                    std::to_string(report));
    }

    /*
      This will check that we have a datafile which report range covers the
      report argument, in adition there can be 'holes' in the series - that must
      be checked by actually querying the data_file object.
    */

    bool has_report(int report) const {
        for (const auto &node : this->index)
            if (node.report1 <= report && node.report2 >= report)
                return true;

        return false;
    }

    IndexNode &back() { return this->index.back(); }

    void clear() { this->index.clear(); }

    size_t length() const { return this->index.back().end(); }

    std::vector<IndexNode>::const_iterator begin() const {
        return this->index.begin();
    }

    std::vector<IndexNode>::const_iterator end() const {
        return this->index.end();
    }

private:
    std::vector<IndexNode> index;
};

} // namespace detail

namespace ecl {
class ecl_sum_file_data;
}

typedef struct ecl_sum_data_struct {
    const ecl_smspec_type *smspec;
    std::vector<ecl::ecl_sum_file_data *> data_files;
    detail::CaseIndex index;
} ecl_sum_data_type;

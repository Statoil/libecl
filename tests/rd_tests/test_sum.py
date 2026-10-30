import csv
import datetime
import os
import os.path
import shutil
import stat

import cwrap
import pytest


def assert_frame_equal(a, b):
    if not a.equals(b):
        raise AssertionError("Expected dataframes to be equal")


try:
    from pandas.testing import assert_frame_equal  # type: ignore
except ImportError:
    pass

from contextlib import contextmanager

from resdata import ResDataType, UnitSystem
from resdata.resfile import FortIO, ResdataFile, ResdataKW, openFortIO
from resdata.summary import Summary, SummaryKeyWordVector, SummaryVarType
from resdata.util.util import CTime, TimeVector

from tests import ResdataTest
from tests.util import TestAreaContext
from tests.util.mock import createSummary


@contextmanager
def pushd(path):
    if not os.path.isdir(path):
        os.makedirs(path)
    cwd = os.getcwd()
    os.chdir(path)

    yield

    os.chdir(cwd)


def create_prediction(history, pred_path):
    restart_case = os.path.join(os.getcwd(), history.base)
    restart_step = history.last_report
    with pushd(pred_path):
        prediction = create_case(
            case="PREDICTION",
            restart_case=restart_case,
            restart_step=restart_step,
            data_start=history.end_date,
        )
        prediction.fwrite()


def fopr(days):
    return days


def fopt(days):
    return days


def fgpt(days):
    if days < 50:
        return days
    else:
        return 100 - days


def create_case(case="CSV", restart_case=None, restart_step=-1, data_start=None):
    length = 100
    return createSummary(
        case,
        [
            ("FOPT", None, 0, "SM3"),
            ("FOPR", None, 0, "SM3/DAY"),
            ("FGPT", None, 0, "SM3"),
        ],
        sim_length_days=length,
        num_report_step=10,
        num_mini_step=10,
        data_start=data_start,
        func_table={"FOPT": fopt, "FOPR": fopr, "FGPT": fgpt},
        restart_case=restart_case,
        restart_step=restart_step,
    )


def create_case2(case="CSV", restart_case=None, restart_step=-1, data_start=None):
    length = 100
    return createSummary(
        case,
        [
            ("WOPT", "OPX", 0, "SM3"),
            ("FOPR", None, 0, "SM3/DAY"),
            ("BPR", None, 10, "SM3"),
            ("RPR", None, 3, "BARS"),
            ("COPR", "OPX", 421, "BARS"),
        ],
        sim_length_days=length,
        num_report_step=10,
        num_mini_step=10,
        data_start=data_start,
        func_table={"FOPT": fopt, "FOPR": fopr, "FGPT": fgpt},
        restart_case=restart_case,
        restart_step=restart_step,
    )


class SumTest(ResdataTest):
    def test_mock(self):
        case = createSummary(
            "CSV", [("FOPT", None, 0, "SM3"), ("FOPR", None, 0, "SM3/DAY")]
        )
        self.assertTrue("FOPT" in case)
        self.assertFalse("WWCT:OPX" in case)

    def test_deprecated_mock_module(self):
        # Exercise the deprecated resdata.util.test.mock.rd_sum_mock helper
        # so that its mainline lines (variable creation, addTStep, key access)
        # are covered. We cannot import the parent package normally because
        # resdata.util.test.__init__ registers C types that conflict with
        # tests.util, so we install stub parent packages first and load the
        # target module directly from its source path under its real name.
        import importlib.util
        import os
        import sys
        import types
        import warnings

        import resdata

        for parent in ("resdata.util.test", "resdata.util.test.mock"):
            if parent not in sys.modules:
                stub = types.ModuleType(parent)
                stub.__path__ = []
                sys.modules[parent] = stub

        path = os.path.join(
            os.path.dirname(resdata.__file__),
            "util",
            "test",
            "mock",
            "rd_sum_mock.py",
        )
        module_name = "resdata.util.test.mock.rd_sum_mock"
        spec = importlib.util.spec_from_file_location(module_name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            case = module.createSummary(
                "CSV",
                [("FOPT", None, 0, "SM3"), ("WOPR", "OP_1", 0, "SM3/DAY")],
            )
        self.assertIn("FOPT", case)
        self.assertIn("WOPR:OP_1", case)

    def test_TIME_special_case(self):
        case = createSummary(
            "CSV", [("FOPT", None, 0, "SM3"), ("FOPR", None, 0, "SM3/DAY")]
        )
        keys = case.keys()
        self.assertEqual(len(keys), 2)
        self.assertIn("FOPT", keys)
        self.assertIn("FOPR", keys)

        keys = case.keys(pattern="*")
        self.assertEqual(len(keys), 2)
        self.assertIn("FOPT", keys)
        self.assertIn("FOPR", keys)

    def test_identify_var_type(self):
        self.assertEnumIsFullyDefined(
            SummaryVarType, "rd_smspec_var_type", "lib/include/resdata/smspec_node.hpp"
        )
        self.assertEqual(
            Summary.var_type("WWCT:OP_X"), SummaryVarType.RD_SMSPEC_WELL_VAR
        )
        self.assertEqual(Summary.var_type("RPR"), SummaryVarType.RD_SMSPEC_REGION_VAR)
        self.assertEqual(Summary.var_type("WNEWTON"), SummaryVarType.RD_SMSPEC_MISC_VAR)
        self.assertEqual(Summary.var_type("YEAR"), SummaryVarType.RD_SMSPEC_MISC_VAR)
        self.assertEqual(Summary.var_type("MONTH"), SummaryVarType.RD_SMSPEC_MISC_VAR)
        self.assertEqual(
            Summary.var_type("AARQ:4"), SummaryVarType.RD_SMSPEC_AQUIFER_VAR
        )

        self.assertEqual(
            Summary.var_type("RXFT"), SummaryVarType.RD_SMSPEC_REGION_2_REGION_VAR
        )
        self.assertEqual(
            Summary.var_type("RxxFT"), SummaryVarType.RD_SMSPEC_REGION_2_REGION_VAR
        )
        self.assertEqual(
            Summary.var_type("RXFR"), SummaryVarType.RD_SMSPEC_REGION_2_REGION_VAR
        )
        self.assertEqual(
            Summary.var_type("RxxFR"), SummaryVarType.RD_SMSPEC_REGION_2_REGION_VAR
        )
        self.assertEqual(Summary.var_type("RORFR"), SummaryVarType.RD_SMSPEC_REGION_VAR)

        case = createSummary(
            "CSV",
            [
                ("FOPT", None, 0, "SM3"),
                ("FOPR", None, 0, "SM3/DAY"),
                ("AARQ", None, 10, "???"),
                ("RGPT", None, 1, "SM3"),
            ],
        )

        node1 = case.smspec_node("FOPT")
        self.assertEqual(node1.var_type(), SummaryVarType.RD_SMSPEC_FIELD_VAR)
        self.assertIsNone(node1.wgname)

        node2 = case.smspec_node("AARQ:10")
        self.assertEqual(node2.var_type(), SummaryVarType.RD_SMSPEC_AQUIFER_VAR)
        self.assertEqual(node2.get_num(), 10)
        self.assertFalse(node2.is_rate())

        node3 = case.smspec_node("RGPT:1")
        self.assertEqual(node3.var_type(), SummaryVarType.RD_SMSPEC_REGION_VAR)
        self.assertEqual(node3.get_num(), 1)
        self.assertEqual(node3.num, 1)
        self.assertEqual(node3.default, 0)
        self.assertTrue(node3.is_total())
        self.assertFalse(node3.is_historical())

        self.assertEqual(node3.keyword, "RGPT")

        self.assertLess(node1, node3)
        self.assertGreater(node2, node3)
        self.assertEqual(node1, node1)
        self.assertNotEqual(node1, node2)

        with self.assertRaises(TypeError):
            a = node1 < 1

    def test_csv_export(self):
        case = createSummary(
            "CSV", [("FOPT", None, 0, "SM3"), ("FOPR", None, 0, "SM3/DAY")]
        )
        sep = ";"
        with TestAreaContext("resdata/csv"):
            case.export_csv("file.csv", sep=sep)
            self.assertTrue(os.path.isfile("file.csv"))
            input_file = csv.DictReader(open("file.csv"), delimiter=sep)
            for row in input_file:
                self.assertIn("DAYS", row)
                self.assertIn("DATE", row)
                self.assertIn("FOPT", row)
                self.assertIn("FOPR", row)
                self.assertEqual(len(row), 4)
                break

            self.assertEqual(case.unit("FOPT"), "SM3")

        with TestAreaContext("resdata/csv"):
            case.export_csv("file.csv", keys=["FOPT"], sep=sep)
            self.assertTrue(os.path.isfile("file.csv"))
            input_file = csv.DictReader(open("file.csv"), delimiter=sep)
            for row in input_file:
                self.assertIn("DAYS", row)
                self.assertIn("DATE", row)
                self.assertIn("FOPT", row)
                self.assertEqual(len(row), 3)
                break

        with TestAreaContext("resdata/csv"):
            date_format = "%y-%m-%d"
            sep = ","
            case.export_csv("file.csv", keys=["F*"], sep=sep, date_format=date_format)
            self.assertTrue(os.path.isfile("file.csv"))
            with open("file.csv") as f:
                time_index = -1
                for line in f:
                    tmp = line.split(sep)
                    self.assertEqual(len(tmp), 4)

                    if time_index >= 0:
                        d = datetime.datetime.strptime(tmp[1], date_format)
                        self.assertEqual(case.iget_date(time_index), d)

                    time_index += 1

    def test_solve(self):
        length = 100
        case = create_case()
        self.assert_solve(case)

    def assert_solve(self, case):
        with self.assertRaises(KeyError):
            case.solve_days("MISSING:KEY", 0.56)

        sol = case.solve_days("FOPT", 150)
        self.assertEqual(len(sol), 0)

        sol = case.solve_days("FOPT", -10)
        self.assertEqual(len(sol), 0)

        sol = case.solve_days("FOPT", 50)
        self.assertEqual(len(sol), 1)
        self.assertFloatEqual(sol[0], 50)

        sol = case.solve_days("FOPT", 50.50)
        self.assertEqual(len(sol), 1)
        self.assertFloatEqual(sol[0], 50.50)

        sol = case.solve_days("FOPR", 50.90)
        self.assertEqual(len(sol), 1)
        self.assertFloatEqual(sol[0], 50.00 + 1.0 / 86400)

        sol = case.solve_dates("FOPR", 50.90)
        t = (
            case.get_data_start_time()
            + datetime.timedelta(days=50)
            + datetime.timedelta(seconds=1)
        )
        self.assertEqual(sol[0], t)

        sol = case.solve_days("FOPR", 50.90, rates_clamp_lower=False)
        self.assertEqual(len(sol), 1)
        self.assertFloatEqual(sol[0], 51.00)

        sol = case.solve_days("FGPT", 25.0)
        self.assertEqual(len(sol), 2)
        self.assertFloatEqual(sol[0], 25.00)
        self.assertFloatEqual(sol[1], 75.00)

        sol = case.solve_dates("FGPT", 25)
        self.assertEqual(len(sol), 2)
        t0 = case.get_data_start_time()
        t1 = t0 + datetime.timedelta(days=25)
        t2 = t0 + datetime.timedelta(days=75)
        self.assertEqual(sol[0], t1)
        self.assertEqual(sol[1], t2)

    def test_different_names(self):
        case = create_case()
        with TestAreaContext("sum_different"):
            case.fwrite()
            shutil.move("CSV.SMSPEC", "CSVX.SMSPEC")
            with self.assertRaises(IOError):
                case2 = Summary.load("Does/not/exist", "CSV.UNSMRY")

            with self.assertRaises(IOError):
                case2 = Summary.load("CSVX.SMSPEC", "CSVX.UNSMRY")

            case2 = Summary.load("CSVX.SMSPEC", "CSV.UNSMRY")
            self.assert_solve(case2)
            self.assertEqual(case.unit("FOPR"), "SM3/DAY")

    def test_invalid(self):
        case = create_case()
        with TestAreaContext("sum_invalid"):
            case.fwrite()
            with open("CASE.txt", "w") as f:
                f.write("No - this is not ResdataKW file ....")

            with self.assertRaises(IOError):
                case2 = Summary.load("CSV.SMSPEC", "CASE.txt")

            with self.assertRaises(IOError):
                case2 = Summary.load("CASE.txt", "CSV.UNSMRY")

            kw1 = ResdataKW("TEST1", 30, ResDataType.RD_INT)
            kw2 = ResdataKW("TEST2", 30, ResDataType.RD_INT)

            with openFortIO("CASE.KW", FortIO.WRITE_MODE) as f:
                kw1.fwrite(f)
                kw2.fwrite(f)

            with self.assertRaises(IOError):
                case2 = Summary.load("CSV.SMSPEC", "CASE.KW")

            with self.assertRaises(IOError):
                case2 = Summary.load("CASE.KW", "CSV.UNSMRY")

    def test_kw_vector(self):
        case1 = create_case()
        case2 = createSummary(
            "CSV",
            [
                ("FOPR", None, 0, "SM3/DAY"),
                ("FOPT", None, 0, "SM3"),
                ("FWPT", None, 0, "SM3"),
            ],
            sim_length_days=100,
            num_report_step=10,
            num_mini_step=10,
            func_table={"FOPT": fopt, "FOPR": fopr, "FWPT": fgpt},
        )
        self.assertEqual(case2.get_key_index("FOPR"), 1)

        kw_list = SummaryKeyWordVector(case1)
        kw_list.add_keyword("FOPT")
        kw_list.add_keyword("FGPT")
        kw_list.add_keyword("FOPR")

        t = case1.get_data_start_time() + datetime.timedelta(days=43)
        data = case1.get_interp_row(kw_list, t)
        for d1, d2 in zip(
            data,
            [
                case1.get_interp("FOPT", date=t),
                case1.get_interp("FOPT", date=t),
                case1.get_interp("FOPT", date=t),
            ],
        ):
            self.assertFloatEqual(d1, d2)

        tmp = []
        for key in kw_list:
            tmp.append(key)

        for k1, k2 in zip(kw_list, tmp):
            self.assertEqual(k1, k2)

        kw_list2 = kw_list.copy(case2)
        self.assertIn("FOPT", kw_list2)
        self.assertIn("FOPR", kw_list2)
        self.assertIn("FGPT", kw_list2)
        data2 = case2.get_interp_row(kw_list2, t)

        self.assertEqual(len(data2), 3)
        self.assertEqual(data[0], data2[0])
        self.assertEqual(data[2], data2[2])

        with TestAreaContext("sum_vector"):
            with cwrap.open("f1.txt", "w") as f:
                case1.dump_csv_line(t, kw_list, f)

            with cwrap.open("f2.txt", "w") as f:
                case2.dump_csv_line(t, kw_list2, f)

            with open("f1.txt") as f:
                d1 = f.readline().split(",")

            with open("f2.txt") as f:
                d2 = f.readline().split(",")

            self.assertEqual(d1[0], d2[0])
            self.assertEqual(d1[2], d2[2])
            self.assertEqual(d2[1], "")

    def test_vector_select_all(self):
        case = create_case()
        rd_sum_vector = SummaryKeyWordVector(case, True)
        keys = case.keys()
        self.assertEqual(len(keys), len(rd_sum_vector))
        for key in keys:
            self.assertIn(key, rd_sum_vector)

    def test_first_last(self):
        case = create_case()
        with self.assertRaises(KeyError):
            case.last_value("NO_SUCH_KEY")
        last_fopt = case.last_value("FOPT")
        values = case.get_values("FOPT")
        self.assertEqual(last_fopt, values[-1])

        with self.assertRaises(KeyError):
            case.first_value("NO_SUCH_KEY")

        first_fopt = case.first_value("FOPT")
        self.assertEqual(first_fopt, values[0])

    def test_time_range(self):
        case = create_case()
        with self.assertRaises(ValueError):
            case.time_range(num_timestep=1)

        time_range = case.time_range(num_timestep=10)
        self.assertEqual(len(time_range), 10)
        self.assertEqual(time_range[0], case.get_data_start_time())
        self.assertEqual(time_range[-1], case.get_end_time())

        # Cover clamping when start is before data_start_time and end is after end_time.
        data_start = case.get_data_start_time()
        data_end = case.get_end_time()
        early_start = data_start - datetime.timedelta(days=10)
        late_end = data_end + datetime.timedelta(days=10)
        clamped = case.time_range(start=early_start, end=late_end, num_timestep=5)
        self.assertEqual(clamped[0], data_start)
        self.assertEqual(clamped[-1], data_end)

    def test_resample(self):
        case = create_case()
        time_vector = case.alloc_time_vector(False)
        case2 = case.resample("RS", time_vector)
        time_vector_resample = case2.alloc_time_vector(False)
        first_diff = time_vector_resample.first_neq(time_vector)
        self.assertEqual(time_vector_resample, time_vector)

    # The purpose of this test is to reproduce a slightly contrived error situation.
    #
    # 1. A history simulation is created and stored somewhere in the
    #    filesystem.
    #
    # 2. We create a prediction, which has 'RESTART' reference to
    #    the history case.
    #
    # 3. The prediction case is loaded from disk, with a cwd different from the
    #    location of the predition case.
    #
    # This configuration would previously lead to a bug in the path used to
    # resolve the history case, and the history would silently be ignored.

    def test_restart_abs_path(self):
        with TestAreaContext("restart_test"):
            history = create_case(case="HISTORY")
            history.fwrite()

            pred_path = "prediction"
            create_prediction(history, pred_path)

            pred = Summary(os.path.join(pred_path, "PREDICTION"))
            # The restart case has a maximum length of 72 characters, depending
            # on the path used for $TMP and so on we do not really know here if
            # the restart_case has been set or not.
            if pred.restart_case:
                self.assertTrue(isinstance(pred.restart_case, Summary))
                self.assertEqual(
                    pred.restart_case.case, os.path.join(os.getcwd(), history.case)
                )
                self.assertEqual(pred.restart_step, history.last_report)

                length = pred.sim_length
                pred_times = pred.alloc_time_vector(False)
                hist_times = history.alloc_time_vector(False)

                for index in range(len(hist_times)):
                    self.assertEqual(hist_times[index], pred_times[index])

    def test_restart_too_long_history_path(self):
        with TestAreaContext(
            "restart_test_too_long_path_for_the_eclipse_restart_keyword_1234567890123456789012345678901234567890"
        ):
            history = create_case(case="HISTORY")
            history.fwrite()

            pred_path = "prediction"
            create_prediction(history, pred_path)

            pred = Summary(os.path.join(pred_path, "PREDICTION"))
            self.assertIsNone(pred.restart_case)

    def test_restart_perm_denied(self):
        with TestAreaContext("restart_test"):
            with pushd("history/case1"):
                history = create_case(case="HISTORY")
                history.fwrite()

            prediction = create_case(
                case="PREDICTION",
                restart_case="history/case1/HISTORY",
                data_start=history.end_date,
            )
            prediction.fwrite()

            os.chmod("history", 0)

            # This just tests that we can create a summary instance even if we do not
            # have access to load the history case.
            pred = Summary("PREDICTION")

            os.chmod("history", stat.S_IRWXU)

    def test_units(self):
        case = create_case()
        self.assertEqual(case.unit_system, UnitSystem.METRIC)

        # We do not really have support for writing anything else than the
        # default MERIC unit system. To be able to test the read functionality
        # we therefor monkey-patch the summary files in place.
        with TestAreaContext("unit_test"):
            case = create_case("UNITS")
            case.fwrite()
            case2 = Summary("UNITS")

            kw_list = []
            f = ResdataFile("UNITS.SMSPEC")
            for kw in f:
                if kw.name == "INTEHEAD":
                    kw[0] = 3
                kw_list.append(kw.copy())

            f.close()
            with openFortIO("UNITS.SMSPEC", mode=FortIO.WRITE_MODE) as f:
                for kw in kw_list:
                    kw.fwrite(f)

            case = Summary("UNITS")
            self.assertEqual(case.unit_system, UnitSystem.LAB)

    def test_numpy_vector(self):
        case = create_case()

        with self.assertRaises(KeyError):
            case.numpy_vector("NO_SUCH_KEY")

        numpy_vector = case.numpy_vector("FOPT")
        self.assertEqual(len(numpy_vector), len(case))
        numpy_dates = case.numpy_dates
        self.assertEqual(numpy_dates[0].tolist(), case.get_data_start_time())
        self.assertEqual(numpy_dates[-1].tolist(), case.get_end_time())

        dates = case.dates
        self.assertEqual(dates[0], case.get_data_start_time())
        self.assertEqual(dates[-1], case.get_end_time())

        dates = (
            [datetime.datetime(2000, 1, 1)]
            + case.dates
            + [datetime.datetime(2020, 1, 1)]
        )
        fopr = case.numpy_vector("FOPR", time_index=dates)
        fopt = case.numpy_vector("FOPT", time_index=dates)
        self.assertEqual(len(fopt), len(dates))

        self.assertEqual(fopr[0], 0)
        self.assertEqual(fopr[-1], 0)

        self.assertEqual(fopt[0], 0)
        self.assertEqual(fopt[0], case.first_value("FOPT"))
        self.assertEqual(fopt[-1], case.last_value("FOPT"))

        with self.assertRaises(ValueError):
            v = case.numpy_vector("FOPR", time_index=dates, report_only=True)

        v = case.numpy_vector("FOPR", report_only=True)
        self.assertEqual(len(v), len(case.report_dates))

    def test_vector(self):
        case = create_case()

        # The get_vector method is extremely deprecated.
        v1 = case.get_vector("FOPT")
        v2 = case.get_vector("FOPT", report_only=True)
        s1 = sum([x.value for x in v1])
        s2 = sum([x.value for x in v2])

    def test_wells_and_groups(self):
        case = create_case()
        self.assertEqual(case.wells(), [])
        self.assertEqual(case.groups(), [])

    def test_csv_load(self):
        case = create_case2()
        frame = case.pandas_frame()
        rd_sum = Summary.from_pandas("PANDAS", frame, dims=[20, 10, 5])

        for key in frame.columns:
            self.assertTrue(key in rd_sum)

        df = rd_sum.pandas_frame()
        assert_frame_equal(frame, df)

        rd_sum_less = Summary.from_pandas(
            "PANDAS",
            frame,
            dims=[20, 10, 5],
            headers=["BPR:10", "RPR:3", "COPR:OPX:1,2,3"],
        )
        del frame["WOPT:OPX"]
        del frame["FOPR"]
        df_less = rd_sum_less.pandas_frame()
        assert_frame_equal(frame, df_less)

    def test_total_and_rate(self):
        self.assertTrue(Summary.is_total("FOPT"))
        self.assertTrue(Summary.is_total("WWPT:OP_3"))
        self.assertFalse(Summary.is_total("RPR:2"))

        self.assertTrue(Summary.is_rate("WOPR:OP_4"))
        self.assertFalse(Summary.is_rate("BPR:123"))
        self.assertTrue(Summary.is_rate("FWIR"))

    def test_load_case(self):
        path = os.path.join(
            self.TESTDATA_ROOT, "local/ECLIPSE/cp_simple3/SIMPLE_SUMMARY3"
        )
        case = Summary(path)
        self.assertFloatEqual(case.sim_length, 545.0)

        fopr = case.numpy_vector("FOPR")
        for time_index, value in enumerate(fopr):
            self.assertEqual(fopr[time_index], value)

    def test_load_case_lazy_and_eager(self):
        path = os.path.join(self.TESTDATA_ROOT, "local/ECLIPSE/cp_simple3/SHORT.UNSMRY")
        lazy_dates = Summary(path, lazy_load=True).numpy_dates
        eager_dates = Summary(path, lazy_load=False).numpy_dates
        self.assertEqual(len(lazy_dates), 107)
        self.assertEqual(len(eager_dates), 107)
        for l, e in zip(lazy_dates, eager_dates):
            self.assertEqual(l, e)

    def test_write_not_implemented(self):
        path = os.path.join(
            self.TESTDATA_ROOT, "local/ECLIPSE/cp_simple3/SIMPLE_SUMMARY3"
        )
        case = Summary(path, lazy_load=True)
        self.assertFalse(case.can_write())
        with self.assertRaises(NotImplementedError):
            case.fwrite()

    def test_directory_conflict(self):
        with TestAreaContext("dir_conflict"):
            case = create_case("UNITS")
            case.fwrite()
            os.mkdir("UNITS")
            case2 = Summary("./UNITS")

    def test_resample_extrapolate(self):
        """
        Test resampling of summary with extrapolate option of lower and upper boundaries enabled
        """

        time_points = TimeVector()

        path = os.path.join(
            self.TESTDATA_ROOT, "local/ECLIPSE/cp_simple3/SIMPLE_SUMMARY3"
        )
        rd_sum = Summary(path, lazy_load=True)

        start_time = rd_sum.get_data_start_time() - datetime.timedelta(seconds=86400)
        end_time = rd_sum.get_end_time() + datetime.timedelta(seconds=86400)
        delta = end_time - start_time

        N = 25
        time_points.initRange(
            CTime(start_time),
            CTime(end_time),
            CTime(int(delta.total_seconds() / (N - 1))),
        )
        time_points.append(CTime(end_time))
        resampled = rd_sum.resample(
            "OUTPUT_CASE",
            time_points,
            lower_extrapolation=True,
            upper_extrapolation=True,
        )

        for key in rd_sum.keys():
            self.assertIn(key, resampled)

        self.assertEqual(
            rd_sum.get_data_start_time() - datetime.timedelta(seconds=86400),
            resampled.get_data_start_time(),
        )

        key_not_rate = "FOPT"
        for time_index, t in enumerate(time_points):
            if t < rd_sum.get_data_start_time():
                self.assertFloatEqual(
                    resampled.iget(key_not_rate, time_index),
                    rd_sum.first_value(key_not_rate),
                )
            elif t > rd_sum.get_end_time():
                self.assertFloatEqual(
                    resampled.iget(key_not_rate, time_index),
                    rd_sum.get_last_value(key_not_rate),
                )
            else:
                self.assertFloatEqual(
                    resampled.iget(key_not_rate, time_index),
                    rd_sum.get_interp_direct(key_not_rate, t),
                )

        key_rate = "FOPR"
        for time_index, t in enumerate(time_points):
            if t < rd_sum.get_data_start_time():
                self.assertFloatEqual(resampled.iget(key_rate, time_index), 0)
            elif t > rd_sum.get_end_time():
                self.assertFloatEqual(resampled.iget(key_rate, time_index), 0)
            else:
                self.assertFloatEqual(
                    resampled.iget(key_rate, time_index),
                    rd_sum.get_interp_direct(key_rate, t),
                )


def create_time_vector(lst):
    vec = TimeVector()
    for l in lst:
        vec.append(l)
    return vec


@pytest.mark.parametrize("time_index_type", [list, create_time_vector, tuple])
def test_summary_to_pandas_frame(time_index_type):
    case = create_case()
    dates = time_index_type(
        [datetime.datetime(2000, 1, 1)] + case.dates + [datetime.datetime(2020, 1, 1)]
    )
    frame = case.pandas_frame(column_keys=["FOPT", "FOPR"], time_index=dates)

    fopr = frame["FOPR"]
    fopt = frame["FOPT"]

    assert fopr.iloc[0] == 0
    assert fopr.iloc[-1] == 0

    assert fopt.iloc[0] == 0
    assert fopt.iloc[0] == case.first_value("FOPT")
    assert fopt.iloc[-1] == case.last_value("FOPT")

    frame = case.pandas_frame()
    rows, columns = frame.shape
    assert len(case.keys()) == columns
    assert len(case) == rows


def test_t_step():
    summary = createSummary(
        "CASE",
        [
            ("FOPT", None, 0, "SM3"),
            ("FOPR", None, 0, "SM3/DAY"),
            ("FGPT", None, 0, "SM3"),
        ],
        sim_length_days=100,
        num_report_step=10,
        num_mini_step=10,
        sim_start=datetime.date(2010, 1, 1),
        func_table={"FOPT": fopt, "FOPR": fopr, "FGPT": fgpt},
    )
    assert summary.path is None
    assert summary.end_time.date() == summary.end_date

    t_step = summary.add_t_step(11, 101)
    assert t_step.get_report() == 11
    assert t_step.get_sim_days() == 101
    assert t_step.get_mini_step() == 10 * 10
    assert summary.iget_days(100) == 101
    assert "FOPT" in t_step
    assert t_step["FOPT"] == 0.0
    t_step["FOPT"] = 100.0
    assert t_step["FOPT"] == 100.0
    assert t_step.get_sim_time().datetime() == datetime.datetime(2010, 4, 12, 0, 0)

    node = summary.get_last("FOPT")
    assert node.days == 101
    assert summary.get_report(days=101) == 11

    assert list(summary.get_interp_vector("FOPT", days_list=[1, 10, 20, 30, 101])) == [
        1.0,
        10.0,
        20.0,
        30.0,
        100.0,
    ]

    assert summary.get_from_report("FOPT", 11) == 100.0
    assert summary.first_gt("FOPT", 10).days == 11
    assert summary.first_lt("FOPT", 10).days == 0
    assert summary.first_gt_index("FOPT", 10) == 11
    assert summary.first_lt_index("FOPT", 10) == 0

    assert summary.get_last("FOPT").days == 101.0
    assert summary["FOPT"].first.days == 0


def test_that_end_time_is_sim_start_when_there_are_no_steps():
    summary = createSummary(
        "CASE",
        [],
        sim_length_days=100,
        num_report_step=0,
        num_mini_step=0,
        sim_start=datetime.date(2010, 1, 1),
        func_table={},
    )
    assert summary.end_date == summary.start_date
    # Note that in this situation
    # sum.data_start will raise


@pytest.mark.usefixtures("use_tmpdir")
def test_that_numpy_vector_of_a_writer_without_time_steps_is_empty():
    writer = Summary.writer("WRITER", datetime.datetime(2010, 1, 1), 5, 5, 5)
    writer.add_variable("FOPT", unit="SM3")

    assert len(writer.numpy_vector("FOPT")) == 0


@pytest.mark.usefixtures("use_tmpdir")
def test_that_sim_length_of_a_writer_without_time_steps_is_zero():
    writer = Summary.writer("WRITER", datetime.datetime(2010, 1, 1), 5, 5, 5)

    assert writer.sim_length == 0.0


@pytest.mark.usefixtures("use_tmpdir")
def test_that_first_day_of_a_writer_without_time_steps_raises_index_error():
    writer = Summary.writer("WRITER", datetime.datetime(2010, 1, 1), 5, 5, 5)

    with pytest.raises(IndexError):
        _ = writer.first_day


@pytest.mark.usefixtures("use_tmpdir")
def test_that_first_day_is_the_time_of_the_first_step_in_days():
    writer = Summary.writer("WRITER", datetime.datetime(2010, 1, 1), 5, 5, 5)
    writer.add_variable("FOPT", unit="SM3")
    writer.add_t_step(1, 3.0)
    writer.add_t_step(2, 4.0)

    assert writer.first_day == 3.0

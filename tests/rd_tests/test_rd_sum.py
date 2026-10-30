import datetime
import os.path
from textwrap import dedent

import numpy as np
import pandas as pd
import pytest
import resfo
from cwrap import open as copen
from hypothesis import assume, given
from pandas.testing import assert_frame_equal
from resdata.resfile import FortIO, ResdataKW, openFortIO, openResdataFile
from resdata.summary import (
    Summary,
    SummaryKeyWordVector,
    SummaryVarType,
)
from resdata.summary.rd_sum import date2num
from resdata.util.util import DoubleVector, StringList, TimeVector
from resfo_utilities.testing import (
    Date,
    Simulator,
    Smspec,
    SmspecIntehead,
    SummaryMiniStep,
    SummaryStep,
    UnitSystem,
    Unsmry,
    summaries,
)

from tests import ResdataTest, equinor_test
from tests.util import TestAreaContext
from tests.util.mock import createSummary

from ._grid_fixtures import write_char_kw, write_float_kw, write_int_kw


def test_write_repr():
    """repr(Summary) used to segfault when there is only a startdate"""
    writer = Summary.writer("TEST", datetime.date(2000, 2, 3), 10, 10, 10)
    assert repr(writer).startswith(
        'Summary(name="writer", time=[2000-02-03 00:00:00, 2000-02-03 00:00:00], keys=0) at 0x'
    )


def test_that_year2263_through_pandas_works():
    """Date 2262-04-11 is the upper limit for nanosecond 64-bit
    timestamps, which is what pandas<2 prefers.
    """
    eclsum = Summary.from_pandas(
        "TESTCASE",
        pd.DataFrame(
            [
                {"DATE": datetime.date(2000, 1, 1), "FPR": 200},
                {"DATE": datetime.date(2263, 1, 1), "FPR": 1},
            ]
        ).set_index("DATE"),
    )
    assert eclsum.numpy_dates[0] == np.datetime64("2000-01-01T00:00:00.000")

    # (This is hit by a round-off error due to 32-bit floats):
    assert eclsum.numpy_dates[1] == np.datetime64("2262-12-31T23:57:52.000")

    # Roundtrip:
    assert eclsum.pandas_frame().index[1] == datetime.datetime(2262, 12, 31, 23, 57, 52)


@pytest.mark.usefixtures("use_tmpdir")
def test_dump_csv_line():
    case = createSummary(
        "CASE1",
        [("FOPT", None, 0, "SM3"), ("WOPT", "OP1", 0, "SM3")],
    )
    case.fwrite()

    rd_sum = Summary("CASE1")
    kw_vector = SummaryKeyWordVector(rd_sum)
    kw_vector.add_keywords("F*")

    dtime = datetime.datetime(2010, 6, 1, 0, 0, 0)
    out_path = "dump.csv"
    with copen(out_path, "w") as out_handle:
        rd_sum.dump_csv_line(dtime, kw_vector, out_handle)

    assert os.path.isfile(out_path)
    assert os.path.getsize(out_path) > 0


@equinor_test()
class SummaryTest(ResdataTest):
    def setUp(self):
        self.test_file = self.createTestPath("Equinor/ECLIPSE/Gurbat/ECLIPSE.SMSPEC")
        self.rd_sum = Summary(self.test_file)

    def test_time_range_year(self):
        real_range = self.rd_sum.time_range(interval="1y", extend_end=False)
        extended_range = self.rd_sum.time_range(interval="1y", extend_end=True)
        assert real_range[-1] < extended_range[-1]

    def test_time_range_day(self):
        real_range = self.rd_sum.time_range(interval="1d", extend_end=False)
        extended_range = self.rd_sum.time_range(interval="1d", extend_end=True)
        assert real_range[-1] == extended_range[-1]

    def test_time_range_month(self):
        real_range = self.rd_sum.time_range(interval="1m", extend_end=False)
        extended_range = self.rd_sum.time_range(interval="1m", extend_end=True)
        assert real_range[-1] < extended_range[-1]

    def test_dump_csv_line(self):
        rd_sum_vector = SummaryKeyWordVector(self.rd_sum)
        rd_sum_vector.add_keywords("F*")

        with self.assertRaises(KeyError):
            rd_sum_vector.add_keyword("MISSING")

        dtime = datetime.datetime(2002, 1, 1, 0, 0, 0)
        with TestAreaContext("Summary/csv_dump"):
            test_file_name = self.createTestPath("dump.csv")
            outputH = copen(test_file_name, "w")
            self.rd_sum.dump_csv_line(dtime, rd_sum_vector, outputH)
            assert os.path.isfile(test_file_name)

    def test_truncated_smspec(self):
        with TestAreaContext("Summary/truncated_smspec") as ta:
            ta.copy_file(self.test_file)
            ta.copy_file(self.createTestPath("Equinor/ECLIPSE/Gurbat/ECLIPSE.UNSMRY"))

            file_size = os.path.getsize("ECLIPSE.SMSPEC")
            with open("ECLIPSE.SMSPEC", "r+") as f:
                f.truncate(file_size // 2)

            with self.assertRaises(IOError):
                Summary("ECLIPSE")

    def test_truncated_data(self):
        with TestAreaContext("Summary/truncated_data") as ta:
            ta.copy_file(self.test_file)
            ta.copy_file(self.createTestPath("Equinor/ECLIPSE/Gurbat/ECLIPSE.UNSMRY"))

            file_size = os.path.getsize("ECLIPSE.UNSMRY")
            with open("ECLIPSE.UNSMRY", "r+") as f:
                f.truncate(file_size // 2)

            with self.assertRaises(IOError):
                Summary("ECLIPSE")

    def test_missing_smspec_keyword(self):
        with TestAreaContext("Summary/truncated_data") as ta:
            ta.copy_file(self.test_file)
            ta.copy_file(self.createTestPath("Equinor/ECLIPSE/Gurbat/ECLIPSE.UNSMRY"))

            with openResdataFile("ECLIPSE.SMSPEC") as f:
                kw_list = []
                for kw in f:
                    kw_list.append(ResdataKW.copy(kw))

            with openFortIO("ECLIPSE.SMSPEC", mode=FortIO.WRITE_MODE) as f:
                for kw in kw_list:
                    if kw.get_name() == "KEYWORDS":
                        continue
                    kw.fwrite(f)

            with self.assertRaises(IOError):
                Summary("ECLIPSE")

    def test_missing_unsmry_keyword(self):
        with TestAreaContext("Summary/truncated_data") as ta:
            ta.copy_file(self.test_file)
            ta.copy_file(self.createTestPath("Equinor/ECLIPSE/Gurbat/ECLIPSE.UNSMRY"))

            with openResdataFile("ECLIPSE.UNSMRY") as f:
                kw_list = []
                for kw in f:
                    kw_list.append(ResdataKW.copy(kw))

            with openFortIO("ECLIPSE.UNSMRY", mode=FortIO.WRITE_MODE) as f:
                c = 0
                for kw in kw_list:
                    if kw.get_name() == "PARAMS":
                        if c % 5 == 0:
                            continue
                    c += 1
                    kw.fwrite(f)

            with self.assertRaises(IOError):
                Summary("ECLIPSE")

    def test_labscale(self):
        case = self.createTestPath("Equinor/ECLIPSE/LabScale/HDMODEL")
        summary = Summary(case, lazy_load=True)
        self.assertEqual(
            summary.get_start_time(), datetime.datetime(2013, 1, 1, 0, 0, 0)
        )
        self.assertEqual(
            summary.get_end_time(), datetime.datetime(2013, 1, 1, 19, 30, 0)
        )
        self.assertFloatEqual(summary.getSimulationLength(), 19.50)

        summary = Summary(case, lazy_load=False)
        self.assertEqual(
            summary.get_start_time(), datetime.datetime(2013, 1, 1, 0, 0, 0)
        )
        self.assertEqual(
            summary.get_end_time(), datetime.datetime(2013, 1, 1, 19, 30, 0)
        )
        self.assertFloatEqual(summary.getSimulationLength(), 19.50)


@given(summaries())
@pytest.mark.usefixtures("use_tmpdir")
def test_to_from_pandas(summary):
    smspec, unsmry = summary
    assume(len(smspec.keywords) == len(set(smspec.keywords)))
    smspec.to_file("TEST.SMSPEC")
    unsmry.to_file("TEST.UNSMRY")
    summary = Summary("TEST", lazy_load=False)

    baseline = summary.pandas_frame()
    roundtrip = Summary.from_pandas(
        "TEST", baseline, dims=(smspec.nx, smspec.ny, smspec.nz)
    ).pandas_frame()

    roundtrip.index = roundtrip.index.values.astype(np.int64)
    baseline.index = baseline.index.values.astype(np.int64)

    assert_frame_equal(roundtrip, baseline, check_exact=False, atol=0, rtol=1e-6)


@given(summaries())
@pytest.mark.usefixtures("use_tmpdir")
def test_that_non_matching_dataframe_gives_empty_columns(summary):
    smspec, unsmry = summary
    assume("BOGUS" not in (smspec.keywords))
    smspec.to_file("TEST.SMSPEC")
    unsmry.to_file("TEST.UNSMRY")

    assert (
        Summary("TEST", lazy_load=False)
        .pandas_frame(column_keys=["BOGUS"])
        .columns.empty
    )


def create_summary(
    summary_keys=("FOPR",),
    time_units="DAYS",
    case="TEST",
    formatted=False,
    start_date=None,
    times=(0.0, 1.0, 2.0),
    values=None,
):
    """Helper function to create test summary files

    Args:
        summary_keys: List of summary keywords
        time_units: Time units (DAYS, HOURS, etc)
        case: Case name
        formatted: Whether to create formatted files
        start_date: Simulation start date
        times: List of time values
        values: List of value lists for each time
    """
    start_date = start_date or Date(
        day=1, month=1, year=2000, hour=0, minutes=0, micro_seconds=0
    )
    num_values = len(summary_keys)

    if values is None:
        values = [[100.0 + i * 10.0] * num_values for i in range(len(times))]

    unsmry = Unsmry(
        steps=[
            SummaryStep(
                seqnum=i,
                ministeps=[
                    SummaryMiniStep(mini_step=0, params=[time] + values[i]),
                ],
            )
            for i, time in enumerate(times)
        ]
    )
    smspec = Smspec(
        nx=10,
        ny=10,
        nz=10,
        restarted_from_step=0,
        num_keywords=1 + len(summary_keys),
        restart="        ",
        keywords=["TIME    "] + list(summary_keys),
        well_names=[":+:+:+:+"] + [":+:+:+:+"] * len(summary_keys),
        region_numbers=[-32676] + [0] * len(summary_keys),
        units=[time_units.ljust(8)] + ["SM3     "] * len(summary_keys),
        start_date=start_date,
        intehead=SmspecIntehead(
            unit=UnitSystem.METRIC,
            simulator=Simulator.ECLIPSE_100,
        ),
    )
    _format = resfo.Format.FORMATTED if formatted else resfo.Format.UNFORMATTED
    smspec.to_file(f"{case}.{'F' if formatted else ''}SMSPEC", file_format=_format)
    unsmry.to_file(f"{case}.{'F' if formatted else ''}UNSMRY", file_format=_format)


@pytest.mark.usefixtures("use_tmpdir")
def test_has_key():
    create_summary(summary_keys=("FOPR", "FOPT", "FWCT"))

    summary = Summary("TEST")
    assert summary.has_key("FOPR")
    assert summary.has_key("FOPT")
    assert summary.has_key("FWCT")
    assert not summary.has_key("NONEXISTENT")


@pytest.mark.usefixtures("use_tmpdir")
def test_that_the_start_date_property_is_the_start_date_from_the_file():
    create_summary(
        start_date=Date(day=1, month=1, year=2000, hour=0, minutes=0, micro_seconds=0)
    )
    summary = Summary("TEST")
    assert summary.start_date == datetime.date(2000, 1, 1)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_numpy_vector():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 1.0, 2.0, 3.0),
        values=[[100.0], [110.0], [120.0], [130.0]],
    )

    summary = Summary("TEST")
    data = summary.numpy_vector("FOPR")

    assert data.tolist() == pytest.approx(
        [
            100.0,
            110.0,
            120.0,
            130.0,
        ]
    )


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_first_and_last_value():
    create_summary(
        summary_keys=("FOPT",),
        times=(0.0, 1.0, 2.0),
        values=[[100.0], [200.0], [300.0]],
    )

    summary = Summary("TEST")
    assert summary.first_value("FOPT") == pytest.approx(100.0)
    assert summary.last_value("FOPT") == pytest.approx(300.0)
    assert summary.get_last("FOPT").value == pytest.approx(300.0)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_dates_and_times():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0, 3.0))

    summary = Summary("TEST")
    dates = summary.dates

    assert len(dates) == 4
    assert dates[0] == datetime.datetime(2000, 1, 1, 0, 0, 0)
    assert dates[1] == datetime.datetime(2000, 1, 2, 0, 0, 0)

    days = summary.days
    assert len(days) == 4
    assert days[0] == pytest.approx(0.0)
    assert days[1] == pytest.approx(1.0)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_keys_filtering():
    create_summary(summary_keys=("FOPR", "FOPT", "FWCT", "FGPR", "FGPT"))

    summary = Summary("TEST")
    assert list(summary.keys()) == ["FGPR", "FGPT", "FOPR", "FOPT", "FWCT"]

    f_keys = summary.keys(pattern="F*")
    assert len(f_keys) >= 5

    fop_keys = summary.keys(pattern="FOP*")
    assert "FOPR" in fop_keys
    assert "FOPT" in fop_keys


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_pandas_frame():
    create_summary(
        summary_keys=("FOPR", "FOPT"),
        times=(0.0, 1.0, 2.0),
        values=[[100.0, 1000.0], [110.0, 1100.0], [120.0, 1200.0]],
    )

    summary = Summary("TEST")
    df = summary.pandas_frame(column_keys=["FOPR", "FOPT"])

    assert set(df.columns) == {"FOPR", "FOPT"}
    assert df["FOPR"].to_numpy().tolist() == pytest.approx([100.0, 110.0, 120.0])
    assert df["FOPT"].to_numpy().tolist() == pytest.approx([1000.0, 1100.0, 1200.0])


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_iget_access():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 1.0, 2.0),
        values=[[100.0], [200.0], [300.0]],
    )

    summary = Summary("TEST")
    assert summary.iget("FOPR", 0) == pytest.approx(100.0)
    assert summary.iget("FOPR", 1) == pytest.approx(200.0)
    assert summary.iget("FOPR", 2) == pytest.approx(300.0)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_lazy_vs_eager_load():
    create_summary(summary_keys=("FOPR", "FOPT"), times=(0.0, 1.0, 2.0))

    sum_lazy = Summary("TEST", lazy_load=True)
    sum_eager = Summary("TEST", lazy_load=False)

    assert len(sum_lazy) == len(sum_eager)

    lazy_data = sum_lazy.numpy_vector("FOPR")
    eager_data = sum_eager.numpy_vector("FOPR")

    assert np.allclose(lazy_data, eager_data)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_contains_operator():
    create_summary(summary_keys=("FOPR", "FOPT"))

    summary = Summary("TEST")
    assert "FOPR" in summary
    assert "FOPT" in summary
    assert "NONEXISTENT" not in summary


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_iteration():
    keys = ["FOPR", "FOPT", "FWCT"]
    create_summary(summary_keys=keys)

    assert list(Summary("TEST")) == keys


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_start_and_end_time():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 5.0, 10.0),
        start_date=Date(day=15, month=3, year=2010, hour=0, minutes=0, micro_seconds=0),
    )

    summary = Summary("TEST")
    assert summary.start_time == datetime.datetime(2010, 3, 15, 0, 0, 0)
    assert summary.end_time == datetime.datetime(2010, 3, 25, 0, 0, 0)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_sim_length():
    create_summary(summary_keys=("FOPR",), times=(0.0, 10.0, 20.0, 30.0))

    summary = Summary("TEST")
    assert summary.sim_length == pytest.approx(30.0)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_iget_date():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))

    summary = Summary("TEST")
    date0 = summary.iget_date(0)
    date1 = summary.iget_date(1)

    assert date0 == datetime.datetime(2000, 1, 1, 0, 0, 0)
    assert date1 == datetime.datetime(2000, 1, 2, 0, 0, 0)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_iget_days():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.5, 3.0))

    summary = Summary("TEST")
    assert summary.iget_days(0) == pytest.approx(0.0)
    assert summary.iget_days(1) == pytest.approx(1.5)
    assert summary.iget_days(2) == pytest.approx(3.0)


def test_summary_var_type():
    assert Summary.var_type("AAQR") == SummaryVarType.RD_SMSPEC_AQUIFER_VAR
    assert Summary.var_type("BPR") == SummaryVarType.RD_SMSPEC_BLOCK_VAR
    assert Summary.var_type("FOPT") == SummaryVarType.RD_SMSPEC_FIELD_VAR
    assert Summary.var_type("GWPR") == SummaryVarType.RD_SMSPEC_GROUP_VAR
    assert Summary.var_type("RPR") == SummaryVarType.RD_SMSPEC_REGION_VAR
    assert Summary.var_type("WOPR") == SummaryVarType.RD_SMSPEC_WELL_VAR
    assert Summary.var_type("CGORL") == SummaryVarType.RD_SMSPEC_COMPLETION_VAR
    assert Summary.var_type("SOPR") == SummaryVarType.RD_SMSPEC_SEGMENT_VAR
    assert Summary.var_type("LBPR") == SummaryVarType.RD_SMSPEC_LOCAL_BLOCK_VAR
    assert Summary.var_type("LCGOR") == SummaryVarType.RD_SMSPEC_LOCAL_COMPLETION_VAR
    assert Summary.var_type("LWOPR") == SummaryVarType.RD_SMSPEC_LOCAL_WELL_VAR


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_export_csv():
    create_summary(summary_keys=("FOPR", "FOPT"), times=(0.0, 1.0, 2.0))

    summary = Summary("TEST")
    summary.export_csv("output.csv", keys=["FOPR", "FOPT"])

    assert os.path.exists("output.csv")
    with open("output.csv") as f:
        content = f.read()
        assert "FOPR" in content
        assert "FOPT" in content


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_check_sim_time():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 10.0, 20.0),
        start_date=Date(day=1, month=1, year=2000, hour=0, minutes=0, micro_seconds=0),
    )

    summary = Summary("TEST")
    assert summary.check_sim_time(datetime.datetime(2000, 1, 1, 0, 0, 0))
    assert summary.check_sim_time(datetime.datetime(2000, 1, 15, 0, 0, 0))
    assert summary.check_sim_time(datetime.datetime(2000, 1, 21, 0, 0, 0))
    assert not summary.check_sim_time(datetime.datetime(2000, 1, 31, 0, 0, 0))


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_numpy_dates():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))

    summary = Summary("TEST")
    np_dates = summary.numpy_dates

    assert np_dates.tolist() == [
        datetime.datetime(2000, 1, 1, 0, 0),
        datetime.datetime(2000, 1, 2, 0, 0),
        datetime.datetime(2000, 1, 3, 0, 0),
    ]
    assert np_dates.dtype == np.dtype("datetime64[ms]")


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_key_error_on_missing_key():
    create_summary(summary_keys=("FOPR",))

    summary = Summary("TEST")
    with pytest.raises(KeyError):
        summary.numpy_vector("NONEXISTENT")

    with pytest.raises(KeyError):
        summary.last_value("NONEXISTENT")

    with pytest.raises(KeyError):
        summary.first_value("NONEXISTENT")


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_writer_basic():
    start_date = datetime.datetime(2000, 1, 1)
    writer = Summary.writer("WRITER_TEST", start_date, 10, 10, 10)

    var = writer.add_variable("FOPT", unit="SM3")
    t_step = writer.add_t_step(1, 1.0)
    t_step[var.get_key1()] = 1000.0

    t_step = writer.add_t_step(2, 2.0)
    t_step[var.get_key1()] = 2000.0

    writer.fwrite()

    loaded = Summary("WRITER_TEST")
    assert "FOPT" in loaded
    assert loaded.first_value("FOPT") == pytest.approx(1000.0)
    assert loaded.last_value("FOPT") == pytest.approx(2000.0)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_writer_multiple_variables():
    start_date = datetime.datetime(2005, 5, 10)
    writer = Summary.writer("MULTI_VAR", start_date, 5, 5, 5)
    assert writer.base == "MULTI_VAR"
    assert writer.path is None

    fopt = writer.add_variable("FOPT", unit="SM3")
    fopr = writer.add_variable("FOPR", unit="SM3/DAY")
    fwpt = writer.add_variable("FWPT", unit="SM3")
    wopt = writer.add_variable("WOPT", unit="SM3", wgname="OP1")
    rpr = writer.add_variable("RPR", unit="SM3/DAY", num=1)
    aaqr = writer.add_variable("AAQR", unit="SM3/DAY", num=1)
    boft = writer.add_variable("BOFT", unit="SM3/DAY", num=1)
    gvpr = writer.add_variable("GVPR", unit="SM3/DAY", wgname="GROUP1")
    sofr = writer.add_variable("SOFR", unit="SM3/DAY", wgname="OP1", num=2)
    cgor = writer.add_variable("CGOR", unit="SM3/SM3", wgname="OP1", num=2)
    rsft = writer.add_variable("RSFT", unit="SM3", num=2)

    for i in range(1, 5):
        t_step = writer.add_t_step(i, float(i))
        t_step[fopt.get_key1()] = 100.0 * i
        t_step[fopr.get_key1()] = 10.0 * i
        t_step[fwpt.get_key1()] = 50.0 * i
        t_step[wopt.get_key1()] = 30.0 * i
        t_step[rpr.get_key1()] = 1000.0 * i
        t_step[aaqr.get_key1()] = 2000.0 * i
        t_step[boft.get_key1()] = 0.1 * i
        t_step[gvpr.get_key1()] = 0.3 * i
        t_step[sofr.get_key1()] = 1.0 * i
        t_step[cgor.get_key1()] = 0.9 * i
        t_step[rsft.get_key1()] = 3000.0 * i

    writer.fwrite()

    loaded = Summary("MULTI_VAR")
    assert len(loaded) == 4

    assert loaded.get_values("FOPT") == pytest.approx([100.0 * i for i in range(1, 5)])
    assert loaded.groups() == ["GROUP1"]
    assert loaded.groups("*") == ["GROUP1"]
    assert loaded.groups("NO") == []

    assert loaded.pandas_frame().to_dict() == {
        "FOPT": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(100.0),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(200.0),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(300.0),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(400.0),
        },
        "FOPR": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(10.0),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(20.0),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(30.0),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(40.0),
        },
        "FWPT": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(50.0),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(100.0),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(150.0),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(200.0),
        },
        "WOPT:OP1": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(30.0),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(60.0),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(90.0),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(120.0),
        },
        "RPR:1": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(1000.0),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(2000.0),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(3000.0),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(4000.0),
        },
        "AAQR:1": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(2000.0),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(4000.0),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(6000.0),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(8000.0),
        },
        "BOFT:1,1,1": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(0.1),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(0.2),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(0.3),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(0.4),
        },
        "GVPR:GROUP1": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(0.3),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(0.6),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(0.9),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(1.2),
        },
        "SOFR:OP1:2": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(1.0),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(2.0),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(3.0),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(4.0),
        },
        "CGOR:OP1:2,1,1": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(0.9),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(1.8),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(2.7),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(3.6),
        },
        "RSFT:2--10": {
            pd.Timestamp("2005-05-11 00:00:00"): pytest.approx(3000.0),
            pd.Timestamp("2005-05-12 00:00:00"): pytest.approx(6000.0),
            pd.Timestamp("2005-05-13 00:00:00"): pytest.approx(9000.0),
            pd.Timestamp("2005-05-14 00:00:00"): pytest.approx(12000.0),
        },
    }
    np.testing.assert_allclose(
        loaded.numpy_vector("FOPT"), [100.0 * i for i in range(1, 5)]
    )
    np.testing.assert_allclose(
        loaded.numpy_vector("FOPR"), [10.0 * i for i in range(1, 5)]
    )
    np.testing.assert_allclose(
        loaded.numpy_vector("FWPT"), [50.0 * i for i in range(1, 5)]
    )
    np.testing.assert_allclose(
        loaded.numpy_vector(wopt.get_key1()), [30.0 * i for i in range(1, 5)]
    )
    np.testing.assert_allclose(
        loaded.numpy_vector(rpr.get_key1()), [1000.0 * i for i in range(1, 5)]
    )
    np.testing.assert_allclose(
        loaded.numpy_vector(aaqr.get_key1()), [2000.0 * i for i in range(1, 5)]
    )
    np.testing.assert_allclose(
        loaded.numpy_vector(boft.get_key1()), [0.1 * i for i in range(1, 5)]
    )
    np.testing.assert_allclose(
        loaded.numpy_vector(gvpr.get_key1()), [0.3 * i for i in range(1, 5)]
    )
    np.testing.assert_allclose(
        loaded.numpy_vector(sofr.get_key1()), [1.0 * i for i in range(1, 5)]
    )
    np.testing.assert_allclose(
        loaded.numpy_vector(cgor.get_key1()), [0.9 * i for i in range(1, 5)]
    )
    np.testing.assert_allclose(
        loaded.numpy_vector(rsft.get_key1()), [3000.0 * i for i in range(1, 5)]
    )


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_writer_local_variables():
    start_date = datetime.datetime(2005, 5, 10)
    # writer does not actually support writing local
    # variables but will add them
    writer = Summary.writer("MULTI_VAR", start_date, 5, 5, 5)
    lwopt = writer.add_variable(
        "LWOPT", unit="SM3", wgname="OP1", lgr="LGR1", lgr_ijk=(0, 0, 0)
    )
    assert lwopt.get_key1() == "LWOPT:LGR1:OP1"
    assert writer.has_key(lwopt.get_key1())

    for i in range(1, 5):
        t_step = writer.add_t_step(i, float(i))
        t_step[lwopt.get_key1()] = 100.0 * i


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_from_pandas_roundtrip():
    dates = pd.date_range("2000-01-01", periods=5, freq="D")
    df = pd.DataFrame(
        {
            "FOPR": [100.0, 110.0, 120.0, 130.0, 140.0],
            "FOPT": [1000.0, 1100.0, 1200.0, 1300.0, 1400.0],
        },
        index=dates,
    )

    sum_from_df = Summary.from_pandas("FROM_PANDAS", df)

    assert "FOPR" in sum_from_df
    assert "FOPT" in sum_from_df
    assert sum_from_df.first_value("FOPR") == pytest.approx(100.0)
    assert sum_from_df.last_value("FOPT") == pytest.approx(1400.0)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_time_range():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 30.0, 60.0, 90.0, 120.0),
        start_date=Date(day=1, month=1, year=2000, hour=0, minutes=0, micro_seconds=0),
    )

    summary = Summary("TEST")
    time_range = summary.time_range(interval="10d")
    assert isinstance(time_range, TimeVector)
    assert [t.datetime() for t in time_range] == [
        datetime.datetime(2000, 1, 1, 0, 0),
        datetime.datetime(2000, 1, 11, 0, 0),
        datetime.datetime(2000, 1, 21, 0, 0),
        datetime.datetime(2000, 1, 31, 0, 0),
        datetime.datetime(2000, 2, 10, 0, 0),
        datetime.datetime(2000, 2, 20, 0, 0),
        datetime.datetime(2000, 3, 1, 0, 0),
        datetime.datetime(2000, 3, 11, 0, 0),
        datetime.datetime(2000, 3, 21, 0, 0),
        datetime.datetime(2000, 3, 31, 0, 0),
        datetime.datetime(2000, 4, 10, 0, 0),
        datetime.datetime(2000, 4, 20, 0, 0),
        datetime.datetime(2000, 4, 30, 0, 0),
    ]

    assert len(time_range) > 0


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_unit():
    create_summary(summary_keys=("FOPR",))

    summary = Summary("TEST")
    unit = summary.unit("FOPR")
    assert unit == "SM3"


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_case_property():
    create_summary(summary_keys=("FOPR",), case="MYCASE")

    summary = Summary("MYCASE")
    assert summary.case == "MYCASE"


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_base_and_path(tmpdir):
    create_summary(summary_keys=("FOPR",))

    summary = Summary("TEST")
    assert summary.case == "TEST"
    assert summary.base == "TEST"
    assert summary.path is None
    assert summary.abs_path == tmpdir


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_length_property():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0, 3.0, 4.0))

    summary = Summary("TEST")
    assert summary.length == 5
    assert len(summary) == 5


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_smspec_node():
    create_summary(summary_keys=("FOPR",))

    summary = Summary("TEST")
    node = summary.smspec_node("FOPR")

    assert node is not None
    assert node.keyword == "FOPR"


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_get_interp():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 2.0, 4.0),
        values=[[100.0], [200.0], [300.0]],
    )

    summary = Summary("TEST")

    # Test at exact time points
    value_at_0 = summary.get_interp("FOPR", days=0.0)
    assert value_at_0 == pytest.approx(100.0)

    value_at_2 = summary.get_interp("FOPR", days=2.0)
    assert value_at_2 == pytest.approx(200.0)

    value_at_4 = summary.get_interp("FOPR", days=4.0)
    assert value_at_4 == pytest.approx(300.0)

    assert summary.get_interp_direct("FOPR", date=datetime.date(2000, 1, 1)) == 100.0
    keyvec = SummaryKeyWordVector(summary)
    keyvec.add_keyword("FOPR")
    assert list(
        summary.get_interp_row(keyvec, datetime.date(2000, 1, 1))
    ) == pytest.approx([100.0])


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_get_interp_vector():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 2.0, 4.0),
        values=[[100.0], [200.0], [300.0]],
    )

    summary = Summary("TEST")

    interp_values = summary.get_interp_vector(
        "FOPR", days_list=[0.0, 1.0, 2.0, 3.0, 4.0]
    )
    assert list(interp_values) == pytest.approx(
        [
            100.0,
            200.0,
            200.0,
            300.0,
            300.0,
        ]
    )


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_pandas_frame_with_time_index():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 5.0, 10.0, 15.0),
        values=[[100.0], [150.0], [200.0], [250.0]],
    )

    summary = Summary("TEST")

    df = summary.pandas_frame(
        time_index=[
            datetime.datetime(2000, 1, 1),
            datetime.datetime(2000, 1, 6),
            datetime.datetime(2000, 1, 11),
        ],
        column_keys=["FOPR"],
    )

    assert list(df.columns) == ["FOPR"]
    assert df["FOPR"].to_numpy().tolist() == pytest.approx([100.0, 150.0, 200.0])


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_static_methods():
    """Test static methods is_rate and is_total"""
    assert Summary.is_rate("FOPR")
    assert not Summary.is_rate("FOPT")
    assert Summary.is_total("FOPT")
    assert not Summary.is_total("FOPR")


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_get_key_index():
    create_summary(summary_keys=("FOPR", "FOPT"))

    summary = Summary("TEST")

    assert summary.get_key_index("FOPR") == 1
    assert summary.get_key_index("FOPT") == 2


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_iiget():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 1.0, 2.0),
        values=[[100.0], [200.0], [300.0]],
    )

    summary = Summary("TEST")
    key_index = summary.get_key_index("FOPR")

    assert summary.iiget(0, key_index) == pytest.approx(100.0)
    assert summary.iiget(1, key_index) == pytest.approx(200.0)
    assert summary.iiget(2, key_index) == pytest.approx(300.0)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_load_method():
    create_summary(summary_keys=("FOPR",), case="LOADTEST")

    summary = Summary.load("LOADTEST.SMSPEC", "LOADTEST.UNSMRY")

    assert "FOPR" in summary
    assert len(summary) == 3


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_data_start():
    create_summary(
        summary_keys=("FOPR",),
        start_date=Date(day=20, month=6, year=2015, hour=0, minutes=0, micro_seconds=0),
    )

    summary = Summary("TEST")
    assert summary.data_start == datetime.datetime(2015, 6, 20, 0, 0, 0)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_first_day():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))

    summary = Summary("TEST")
    assert summary.first_day == pytest.approx(0.0)


def test_date2num():
    dt = datetime.datetime(2000, 1, 1, 12, 30, 0)
    num = date2num(dt)

    assert isinstance(num, float)
    assert num == pytest.approx(730120.520833333)


def test_that_summary_raises_value_error_when_load_case_is_falsy():
    with pytest.raises(ValueError, match="load_case must be the basename"):
        Summary(None)

    with pytest.raises(ValueError, match="load_case must be the basename"):
        Summary("")


@pytest.mark.usefixtures("use_tmpdir")
def test_that_summary_load_raises_ioerror_when_smspec_is_not_a_file():
    with pytest.raises(IOError, match="No such file"):
        Summary.load("NONEXISTENT.SMSPEC", "TEST.UNSMRY")


@pytest.mark.usefixtures("use_tmpdir")
def test_that_summary_load_raises_ioerror_when_unsmry_is_not_a_file():
    create_summary(summary_keys=("FOPR",))

    with pytest.raises(IOError, match="No such file"):
        Summary.load("TEST.SMSPEC", "NONEXISTENT.UNSMRY")


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_wells():
    smspec = Smspec(
        nx=10,
        ny=10,
        nz=10,
        restarted_from_step=0,
        num_keywords=3,
        restart="        ",
        keywords=["TIME    ", "WOPR    ", "WWCT    "],
        well_names=[":+:+:+:+", "WELL1   ", "WELL2   "],
        region_numbers=[-32676, 0, 0],
        units=["DAYS    ", "SM3/DAY ", "        "],
        start_date=Date(day=1, month=1, year=2000, hour=0, minutes=0, micro_seconds=0),
        intehead=SmspecIntehead(
            unit=UnitSystem.METRIC, simulator=Simulator.ECLIPSE_100
        ),
    )
    unsmry = Unsmry(
        steps=[
            SummaryStep(
                seqnum=0,
                ministeps=[SummaryMiniStep(mini_step=0, params=[0.0, 100.0, 0.5])],
            )
        ]
    )
    smspec.to_file("TEST.SMSPEC")
    unsmry.to_file("TEST.UNSMRY")

    summary = Summary("TEST")
    wells = summary.wells()
    assert isinstance(wells, StringList)
    assert list(wells) == ["WELL1", "WELL2"]


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_groups():
    create_summary(summary_keys=("FOPR",))

    summary = Summary("TEST")
    groups = summary.groups()
    assert isinstance(groups, StringList)
    assert list(groups) == []


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_numpy_vector_with_time_index():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 5.0, 10.0),
        values=[[100.0], [150.0], [200.0]],
    )

    summary = Summary("TEST")

    time_index = [
        datetime.datetime(2000, 1, 1),
        datetime.datetime(2000, 1, 6),
    ]

    data = summary.numpy_vector("FOPR", time_index=time_index)
    assert data.tolist() == pytest.approx([100.0, 150.0])


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_time_range_with_start():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 30.0, 60.0, 90.0),
        start_date=Date(day=1, month=1, year=2000, hour=0, minutes=0, micro_seconds=0),
    )

    summary = Summary("TEST")

    custom_start = datetime.datetime(2000, 1, 15)
    time_range = summary.time_range(start=custom_start, interval="10d")

    assert isinstance(time_range, TimeVector)
    assert [t.datetime() for t in time_range] == [
        datetime.datetime(2000, 1, 15, 0, 0),
        datetime.datetime(2000, 1, 25, 0, 0),
        datetime.datetime(2000, 2, 4, 0, 0),
        datetime.datetime(2000, 2, 14, 0, 0),
        datetime.datetime(2000, 2, 24, 0, 0),
        datetime.datetime(2000, 3, 5, 0, 0),
        datetime.datetime(2000, 3, 15, 0, 0),
        datetime.datetime(2000, 3, 25, 0, 0),
        datetime.datetime(2000, 4, 4, 0, 0),
    ]


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_time_range_with_end():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 30.0, 60.0, 90.0),
        start_date=Date(day=1, month=1, year=2000, hour=0, minutes=0, micro_seconds=0),
    )

    summary = Summary("TEST")

    custom_end = datetime.datetime(2000, 2, 15)
    time_range = summary.time_range(end=custom_end, interval="10d")

    assert isinstance(time_range, TimeVector)
    assert [t.datetime() for t in time_range] == [
        datetime.datetime(2000, 1, 1, 0, 0),
        datetime.datetime(2000, 1, 11, 0, 0),
        datetime.datetime(2000, 1, 21, 0, 0),
        datetime.datetime(2000, 1, 31, 0, 0),
        datetime.datetime(2000, 2, 10, 0, 0),
        datetime.datetime(2000, 2, 20, 0, 0),
    ]


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_blocked_production():
    """Test blocked_production method"""
    create_summary(
        summary_keys=("FOPT",),
        times=(0.0, 10.0, 20.0, 30.0),
        values=[[0.0], [1000.0], [2500.0], [4000.0]],
    )

    summary = Summary("TEST")

    time_range = [
        datetime.datetime(2000, 1, 1),
        datetime.datetime(2000, 1, 11),
        datetime.datetime(2000, 1, 21),
    ]

    blocked = summary.blocked_production("FOPT", time_range)
    assert isinstance(blocked, DoubleVector)
    assert list(blocked) == pytest.approx([1000.0, 1500.0])


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_get_report_time():
    create_summary(
        start_date=Date(day=1, month=1, year=2000, hour=0, minutes=0, micro_seconds=0)
    )

    summary = Summary("TEST")

    assert summary.get_report_time(0) == datetime.date(day=1, month=1, year=2000)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_get_interp_vector_with_date_list():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 2.0, 4.0),
        values=[[100.0], [200.0], [300.0]],
    )

    summary = Summary("TEST")

    date_list = [
        datetime.datetime(2000, 1, 1),
        datetime.datetime(2000, 1, 3),
        datetime.datetime(2000, 1, 5),
    ]

    interp_values = summary.get_interp_vector("FOPR", date_list=date_list)
    assert list(interp_values) == pytest.approx([100.0, 200.0, 300.0])


@pytest.mark.usefixtures("use_tmpdir")
def test_that_restart_case_for_cases_without_restart_is_none():
    create_summary(summary_keys=("FOPR",))

    summary = Summary("TEST")
    assert summary.restart_case is None


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_get_days():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))

    summary = Summary("TEST")
    assert summary.get_days() == pytest.approx([0.0, 1.0, 2.0])


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_get_dates():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))

    summary = Summary("TEST")
    assert summary.get_dates() == [
        datetime.datetime(2000, 1, 1, 0, 0),
        datetime.datetime(2000, 1, 2, 0, 0),
        datetime.datetime(2000, 1, 3, 0, 0),
    ]


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_report_step_property():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))

    summary = Summary("TEST")
    assert summary.report_step == [1, 2, 3]
    assert summary.get_report_step() == [1, 2, 3]
    assert [summary.iget_report(i) for i in range(3)] == [1, 2, 3]
    assert summary.first_report == 1
    assert summary.last_report == 3


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_end_date():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 1.0, 2.0),
        start_date=Date(day=10, month=5, year=2010, hour=0, minutes=0, micro_seconds=0),
    )

    summary = Summary("TEST")
    assert summary.end_date == datetime.date(2010, 5, 12)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_first_gt_index():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 1.0, 2.0),
        values=[[100.0], [200.0], [300.0]],
    )

    summary = Summary("TEST")
    assert summary.first_gt_index("FOPR", 150.0) == 1


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_first_lt_index():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 1.0, 2.0),
        values=[[300.0], [200.0], [100.0]],
    )

    summary = Summary("TEST")
    assert summary.first_lt_index("FOPR", 250.0) == 1


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_solve_dates():
    """Test solve_dates method"""
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 1.0, 2.0, 3.0),
        values=[[100.0], [150.0], [200.0], [250.0]],
    )

    summary = Summary("TEST")
    assert summary.solve_dates("FOPR", 175.0) == [
        datetime.datetime(2000, 1, 2, 0, 0, 1)
    ]


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_solve_days():
    """Test solve_days method"""
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 1.0, 2.0, 3.0),
        values=[[100.0], [150.0], [200.0], [250.0]],
    )

    summary = Summary("TEST")
    days = summary.solve_days("FOPR", 175.0)
    assert isinstance(days, DoubleVector)
    assert list(days) == pytest.approx([1.000011])


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_fwrite():
    start_date = datetime.datetime(2000, 1, 1)
    writer = Summary.writer("FWRITE_TEST", start_date, 5, 5, 5)

    var = writer.add_variable("FOPT")
    t_step = writer.add_t_step(1, 1.0)
    t_step[var.get_key1()] = 1000.0

    assert writer.can_write()
    writer.fwrite()

    summary = Summary("FWRITE_TEST")
    assert len(summary) == 1
    assert list(summary.keys()) == ["FOPT"]
    assert summary.days == [1.0]
    kw = summary["FOPT"]
    assert kw.values.tolist() == [1000.0]


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_alloc_time_vector():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))

    summary = Summary("TEST")
    time_vec = summary.alloc_time_vector(False)
    assert isinstance(time_vec, TimeVector)
    assert [t.datetime() for t in time_vec] == [
        datetime.datetime(2000, 1, 1, 0, 0),
        datetime.datetime(2000, 1, 2, 0, 0),
        datetime.datetime(2000, 1, 3, 0, 0),
    ]


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_alloc_data_vector():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))

    summary = Summary("TEST")
    key_index = summary.get_key_index("FOPR")
    data_vec = summary.alloc_data_vector(key_index, False)
    assert isinstance(data_vec, DoubleVector)
    assert list(data_vec) == [100.0, 110.0, 120.0]


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_get_general_var_index():
    create_summary(summary_keys=("FOPR",))

    summary = Summary("TEST")
    assert summary.get_general_var_index("FOPR") == 1


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_export_csv_with_keys_none():
    create_summary(summary_keys=("FOPR", "FOPT", "FWCT"), times=(0.0, 1.0, 2.0))

    summary = Summary("TEST")
    summary.export_csv("all_keys.csv", keys=None)

    assert os.path.exists("all_keys.csv")

    with open("all_keys.csv") as f:
        content = f.read()
    assert content == dedent("""\
        DAYS;DATE;FOPR;FOPT;FWCT
           0.00;2000-01-01;100;100;100
           1.00;2000-01-02;110;110;110
           2.00;2000-01-03;120;120;120
    """)


@pytest.mark.usefixtures("use_tmpdir")
def test_summary_resample():
    create_summary(
        summary_keys=("FOPR",),
        times=(0.0, 5.0, 10.0, 15.0),
        values=[[100.0], [150.0], [200.0], [250.0]],
    )

    summary = Summary("TEST")

    time_points = TimeVector()
    time_points.append(datetime.datetime(2000, 1, 1))
    time_points.append(datetime.datetime(2000, 1, 6))
    time_points.append(datetime.datetime(2000, 1, 11))

    resampled = summary.resample("RESAMPLED", time_points)

    assert resampled.alloc_time_vector(False) == time_points
    assert resampled["FOPR"].values.tolist() == pytest.approx([100.0, 150.0, 200.0])


def test_that_reading_non_existent_summary_raises_oserror(tmp_path):
    with pytest.raises(OSError, match="Failed to create summary"):
        Summary(str(tmp_path / "does_not_exist"))

    with pytest.raises(OSError, match="No such file"):
        Summary.load(
            str(tmp_path / "does_not_exist.SMSPEC"),
            str(tmp_path / "does_not_exist.UNSMRY"),
        )


def test_that_missing_permissions_raises_oserror(tmp_path):
    (smspec_path := tmp_path / "CASE.SMSPEC").write_text("", encoding="utf-8")
    (unsmry_path := tmp_path / "CASE.UNSMRY").write_text("", encoding="utf-8")
    mode = smspec_path.stat().st_mode
    smspec_path.chmod(0o000)  # no permissions
    try:
        with pytest.raises(OSError, match="Failed to create summary"):
            Summary(str(smspec_path))
        with pytest.raises(OSError, match="Failed to create summary"):
            Summary.load(str(smspec_path), str(unsmry_path))
    finally:
        smspec_path.chmod(mode)


@pytest.mark.usefixtures("use_tmpdir")
def test_add_variable_after_tstep_raises():
    writer = Summary.writer("AFTER_TSTEP", datetime.datetime(2000, 1, 1), 5, 5, 5)
    var = writer.add_variable("FOPT")
    t_step = writer.add_t_step(1, 1.0)
    t_step[var.get_key1()] = 1.0

    with pytest.raises(ValueError, match="interchange variable adding and timesteps"):
        writer.add_variable("FWPT")


@pytest.mark.usefixtures("use_tmpdir")
def test_add_local_variable_after_tstep_raises():
    writer = Summary.writer("AFTER_TSTEP_LGR", datetime.datetime(2000, 1, 1), 5, 5, 5)
    var = writer.add_variable("FOPT")
    t_step = writer.add_t_step(1, 1.0)
    t_step[var.get_key1()] = 1.0

    with pytest.raises(ValueError, match="interchange variable adding and timesteps"):
        writer.add_variable("BLOCK", lgr="LGR1", lgr_ijk=(1, 1, 1))


@pytest.mark.usefixtures("use_tmpdir")
def test_iget_date_with_invalid_index_raises():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))
    summary = Summary("TEST")

    with pytest.raises(ValueError, match="Internal error when looking up index"):
        summary.iget_date(99999)


@pytest.mark.usefixtures("use_tmpdir")
def test_get_from_report_with_invalid_step_raises():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))
    summary = Summary("TEST")

    with pytest.raises(ValueError, match="Internal error when looking up report"):
        summary.get_from_report("FOPR", 99999)


@pytest.mark.usefixtures("use_tmpdir")
def test_get_interp_direct_with_out_of_range_time_raises():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))
    summary = Summary("TEST")

    with pytest.raises(IndexError, match="Invalid time_t instance"):
        summary.get_interp_direct("FOPR", datetime.datetime(2099, 1, 1))


@pytest.mark.usefixtures("use_tmpdir")
def test_alloc_data_vector_with_out_of_range_index_raises():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))
    summary = Summary("TEST")

    with pytest.raises(IndexError, match="Out of range"):
        summary.alloc_data_vector(99999, False)


@pytest.mark.usefixtures("use_tmpdir")
def test_get_general_var_index_with_unknown_key_raises():
    create_summary(summary_keys=("FOPR",), times=(0.0, 1.0, 2.0))
    summary = Summary("TEST")

    with pytest.raises(IndexError, match="Invalid lookup summary object"):
        summary.get_general_var_index("NO_SUCH_KEY")


@pytest.mark.parametrize("startdat", [[1], [1, 1]])
def test_that_smspec_with_too_short_startdat_is_rejected(tmp_path, startdat):
    with openFortIO(str(tmp_path / "CASE.SMSPEC"), mode=FortIO.WRITE_MODE) as f:
        write_int_kw(f, "DIMENS", [1, 1, 1, 1, 0, -1])
        write_char_kw(f, "KEYWORDS", ["TIME"])
        write_char_kw(f, "WGNAMES", [":+:+:+:+"])
        write_char_kw(f, "UNITS", ["DAYS"])
        write_int_kw(f, "STARTDAT", startdat)
    with openFortIO(str(tmp_path / "CASE.UNSMRY"), mode=FortIO.WRITE_MODE) as f:
        write_int_kw(f, "SEQHDR", [0])
        write_int_kw(f, "MINISTEP", [0])
        write_float_kw(f, "PARAMS", [0.0])

    with pytest.raises(ValueError, match="Invalid size of STARTDAT"):
        Summary.load(str(tmp_path / "CASE.SMSPEC"), str(tmp_path / "CASE.UNSMRY"))

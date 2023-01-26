#!/usr/bin/env python3

import numpy as np
import os
import pathlib
import shutil
import subprocess
import xarray as xr


def get_test_directory() -> pathlib.Path:
    """Get the directory of this test file"""
    return pathlib.Path(__file__).parent


def get_stella_path() -> pathlib.Path:
    """Returns the absolute path to the stella executable

    Can be controlled by setting the STELLA_EXE_PATH environment variable
    """
    default_path = get_test_directory() / "../../../stella"
    stella_path = pathlib.Path(os.environ.get("STELLA_EXE_PATH", default_path))
    return stella_path.absolute()


def run_stella(stella_path: str, input_filename: str) -> int:
    """Run stella with a given input file"""
    subprocess.run([stella_path, input_filename], check=True)


def copy_input_file(input_filename: str, destination):
    """Copy input_filename to destination directory"""
    shutil.copyfile(get_test_directory() / input_filename, destination / input_filename)


def test_simple_run(tmp_path):
    """Run a short test case to check that the output is generated with
    the expected fields, and the expected number of timesteps

    """

    input_filename = "simple.in"
    copy_input_file(input_filename, tmp_path)
    os.chdir(tmp_path)
    run_stella(get_stella_path(), input_filename)

    with xr.open_dataset(tmp_path / input_filename.replace("in", "out.nc")) as df:
        expected_keys = [
            "phi2",
            "phi_vs_t",
        ]

        for field in expected_keys:
            assert field in df

        assert len(df["t"]) == 11
        assert np.allclose(df["t"][-1], 10.0)


def convert_byte_array(array):
    return "\n".join((str(line, encoding="utf-8").strip() for line in array.data))


def test_golden_answer(tmp_path):
    """Check results are identical to a previous run.

    PLEASE NOTE: IF RESULTS ARE EXPECTED TO CHANGE,
    GOLDEN_ANSWER.out.nc MUST BE UPDATED

    """

    input_filename = "simple.in"
    copy_input_file(input_filename, tmp_path)
    os.chdir(tmp_path)
    run_stella(get_stella_path(), input_filename)

    test_file = tmp_path / input_filename.replace("in", "out.nc")
    golden_file = get_test_directory() / "GOLDEN_ANSWER.out.nc"

    with xr.open_dataset(test_file) as df, xr.open_dataset(golden_file) as golden:
        assert set(df.keys()) == set(golden.keys())

        for key in golden:
            if golden[key].shape == ():
                continue
            if golden[key].dtype.kind == "S":
                golden_str = convert_byte_array(golden[key])
                df_str = convert_byte_array(df[key])
                assert df_str == golden_str, key
            else:
                assert np.allclose(df[key], golden[key], equal_nan=True), key

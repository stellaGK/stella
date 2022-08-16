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


def convert_byte_array(array):
    return "\n".join((str(line, encoding="utf-8").strip() for line in array.data))


def compare_two_output_files(file1, file2, tolerance=1e-6):
    """Compares all the variables in file1 and file2 using `numpy.allclose`"""

    VARIABLES_TO_SKIP = [
        "input_file",
        "code_info",
    ]

    with xr.open_dataset(file1) as df, xr.open_dataset(file2) as golden:
        assert set(df.keys()) == set(golden.keys())

        for key in golden:
            if key in VARIABLES_TO_SKIP:
                continue

            if golden[key].shape == ():
                continue
            if golden[key].dtype.kind == "S":
                golden_str = convert_byte_array(golden[key])
                df_str = convert_byte_array(df[key])
                assert df_str == golden_str, key
            else:
                assert np.allclose(df[key], golden[key], equal_nan=True), key


def test_restart(tmp_path):
    full_input_file = "full_run.in"
    start_half_input_file = "start_half_run.in"
    end_half_input_file = "end_half_run.in"

    for input_file in [full_input_file, start_half_input_file, end_half_input_file]:
        copy_input_file(input_file, tmp_path)

    os.chdir(tmp_path)

    run_stella(get_stella_path(), full_input_file)
    run_stella(get_stella_path(), start_half_input_file)

    shutil.copy(tmp_path / "start_half_run.out.nc", tmp_path / "end_half_run.out.nc")

    run_stella(get_stella_path(), end_half_input_file)

    compare_two_output_files("full_run.out.nc", "end_half_run.out.nc")

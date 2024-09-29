
# Update input file

This is an executable, as well as part of the stella code, which will read the user's `input.in` file. It will first remove comments and write the results to `.input.in`. Next, all default input variables are added, and backwards compatibility is ensured with old stella input files. The results are written to `input_with_defaults.in` which is the input file stella will actually read.

<br>

## Create update-input-file executable

The executable is compiled using 12 threads through:
```
make update-input-file -j 12
```

<br>

## Default stella input file

In order to obtain all input variables used in stella, as well as their default values, the following command can be used:
```
./update-input-file --default-input
```

<br>

## Update the input file

The input file `input.in` can be manually updated through:
```
./update-input-file input.in
```
The results are found inside `input_with_defaults.in`.

<br>

## Useful python scripts

Two python scripts are located in `STELLA_CODE/parameters/input_file_tool`:

- make_default_input_pretty.py
- namelist_to_code.py

The first one is for users, in order to make any input file produced by `./update-input-file` pretty. Which means that variables are written in lower case, trailing digits are removed from numbers, and extra spacing is removed. The second script is useful for developers to write the pieces of code for backwards compatibility.
```
./update-input-file --default-input
python3 make_default_input_pretty.py
```

<br>

## Temporary files

For the cleaning of stella I have included `default_stella_input.in` which is the current state of the input variables. As well as `default_stella_input_after_cleaning.in` which is the structure of the namelist we want to have in a cleaned version of stella.

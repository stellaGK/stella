# stella compilation

The code stella can be compiled with `make` or `CMake`. 
The compilation steps are detailed below.

[TOC]

## Dependencies

stella requires MPI, and has several optional dependencies:

- netCDF Fortran
- FFTW3
- LAPACK

## Compilation

There are two ways to build stella: with CMake (experimental) or with
plain `make`.

### CMake

**Note**: If you have previously built stella with plain `make` you
_must_ run `make clean` before attempting to build with CMake, or the
existing built objects will interfere with the CMake build.

Building stella with CMake requires CMake >= 3.16. You can download
the latest version from the [CMake
website](https://cmake.org/download/), but it is often easier to
install with `pip`:

```
pip install cmake
```

Building stella is then a matter of first configuring the build:

```
cmake . -B COMPILATION/build_cmake
```

and then building proper:

```
cmake --build COMPILATION/build_cmake
```

You may need to pass a few flags to the first `cmake` command to tell
it where to find some dependencies:

```
cmake . -B build \
  -DnetCDFFortran_ROOT=/path/to/netcdf/fortran
  -DFFTW_ROOT=/path/to/fftw
```

There are a few build options:

- `STELLA_ENABLE_LAPACK`: Enable LAPACK (default: on)
- `STELLA_ENABLE_FFT`: Enable FFTs (default: on)
- `STELLA_ENABLE_NETCDF`: Enable NetCDF (default: on)
- `STELLA_ENABLE_DOUBLE`: Promotes precisions of real and complex to double
  (default: on)
- `STELLA_ENABLE_LOCAL_SPFUNC`: Enable local special functions" (default: off)
- `STELLA_ENABLE_NAGLIB`: Use the NAG library (default: off)
- `STELLA_ENABLE_POSIX`: Enable POSIX functions for command line functionality
  (default: off)
- `STELLA_ENABLE_F200X`: Enable use of F2003/F2008 functionality (default: on)

You can turn these on or off with `-D<option name>=ON/OFF`. You can
get a complete list of options by running the following in a build
directory:

```
cmake -LH
```

### Make

The other build system uses plain `make`:

1. Set `STELLA_SYSTEM='system'`, with `system` replaced by the appropriate system on
   which you are running. See the `Makefiles` directory for a list of supported
   systems.
2. Optionally, set the following environment variables to override the locations
   in the `STELLLA_SYSTEM` Makefile:
   - `FFTW_LIB_DIR`: directory containing libfftw3
   - `FFTW_INC_DIR`: directory including fftw3.f
   - `NETCDF_LIB_DIR`: directory containing libnetcdff
   - `NETCDF_INC_DIR`: directory including netcdf.inc
4. Set the environment variable `MAKEFLAGS=-IMakefiles`, or set `-IMakefiles`
   when you run `make`
5. Run `make`

For example, to compile on Ubuntu:

```bash
# using bash:
export STELLA_SYSTEM=gnu_ubuntu
export MAKEFLAGS=-IMakefiles
make

# or in one line:
make -IMakefiles STELLA_SYSTEM=gnu_ubuntu
```


## Guide to the MakeFiles

### Overview of Makefile's

To compile stella with `make` the Makefile, which contains the compilation commands, has been split into three main files:

- COMPILATION/Makefile
- COMPILATION/Makefile.extermals
- COMPILATION/Makefile.stella 

Some of the external libraries and automated tests are compiled using their own Makefile's:

- EXTERNALS/mini_libstell/Makefile
- EXTERNALS/git_version/Makefile
- AUTOMATIC_TESTS/test_fortran_routines/Makefile
- AUTOMATIC_TESTS/test_fortran_routines/Makefile

For each operating system, a specific Makefile is written, e.g.:

- COMPILATION/Makefiles/Makefile.gnu_ubuntu
- COMPILATION/Makefiles/Makefile.debian
- COMPILATION/Makefiles/Makefile.macosx

Which will include one of the compiler Makefile's, e.g.: 

- COMPILATION/Makefiles/Compilers/Makefile.intel
- COMPILATION/Makefiles/Compilers/Makefile.gnu-gfortran 

To include your own Makefile commands, you can add:

- COMPILATION/Makefiles/Makefile.local


### Dependencies of Fortran scripts - fortdep

The dependencies of the Fortran scripts are written automatically by the program `fortdep` using the following command
```
make depend
```

Which generates the following file:

- COMPILATION/Makefile.depend

### Makefile

### Makefile.externals

### Makefile.stella

## Introduction to make

Here we give a [quick introduction](https://makefiletutorial.com/) into `make`. 

### static libraries

In a static library, the modules are bound into the executable file before execution. 
Static libraries are commonly named libname.a. The .a suffix refers to archive.  

A Makefile consists of a set of rules. A rule generally looks like this:

```Makefile
targets: prerequisites
	command
	command
	command
```

- The targets are file names, separated by spaces. Typically, there is only one per rule.
- The commands are a series of steps typically used to make the target(s). 
  These need to start with a tab character, not spaces.
- The prerequisites are also file names, separated by spaces. These files need to 
  exist before the commands for the target are run. These are also called dependencies.
- When we run a target, 'make' will only run it, if the target doesn't exist, or if the
  prerequisites are newer than the target
  
### PHONY targets

When we define a target as `target: ...`, the command `make` will 
build a file based on other files. If our target does not build a 
file we label it as a PHONY target. Therefore, a phony target is 
not really the name of a file; rather it is just a name for a recipe 
to be executed when you make an explicit request. 

### Automatic variables

- $@ is the name of the target being generated 
- $^ are all the prerequisites  
- $< is the first prerequisite 
- $? are all the prerequisites newer than the target

### Compiler options
- -o is a compiler option that specifies the name of the output file
- -c is a compiler option that tells the compiler to generate an object file




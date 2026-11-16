# stella compilation


[TOC]

The dependencies and the compilation with `make` or `CMake` are described in the [main README](../README.md#installation-and-compilation). This file describes how the Makefiles are organised.

## Guide to the Makefiles

### Overview of Makefiles

To compile stella with `make` the Makefile, which contains the compilation commands, has been split into three main files:

- COMPILATION/Makefile
- COMPILATION/Makefile.externals
- COMPILATION/Makefile.stella 

Some of the external libraries and automated tests are compiled using their own Makefiles:

- EXTERNALS/mini_libstell/Makefile
- EXTERNALS/git_version/Makefile
- AUTOMATIC_TESTS/test_fortran_routines/Makefile
- AUTOMATIC_TESTS/test_fortran_routines/Makefile

For each operating system, a specific Makefile is written, e.g.:

- COMPILATION/Makefiles/Makefile.gnu_ubuntu
- COMPILATION/Makefiles/Makefile.debian
- COMPILATION/Makefiles/Makefile.macosx

Which will include one of the compiler Makefiles, e.g.: 

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




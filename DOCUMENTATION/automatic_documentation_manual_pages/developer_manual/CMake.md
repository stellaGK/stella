---
title: CMake in stella
subtitle: Some notes on using and developing CMake for stella
---

Stella now has (experimental) support for building with CMake.


# Developing the stella CMake build system

One important consideration when developing stella is that any new
files _must_ be listed in the `STELLA_SOURCES_*` variables: either
`STELLA_SOURCES_f90` _or_ `STELLA_SOURCES_fpp` as appropriate. If you
add a new file and do not add it to exactly one of these variables,
you will get a build error.

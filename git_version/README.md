Fortran-git-version
===================

A Fortran module plus helper Makefile and CMake module for capturing
the git version of your project in your application or library.

The implementation is split out into a Fortran submodule in order to
avoid re-compilation cascades from `use` of the parent module.

Uses the popular `GetGitRevisionDescription.cmake` from
[rpavlik's](https://github.com/rpavlik/cmake-modules)
[cmake-modules](https://github.com/rpavlik/cmake-modules) collection,
used under the [Boost Software
Licence](https://www.boost.org/LICENSE_1_0.txt).

Requirements
------------

Fortran-git-version requires `git` and a Fortran compiler that supports
submodules -- this should be any recent-ish compiler, for example
gfortran has supported submodules since version 6.

Usage
-----

You should include a copy of Fortran-git-version in your project. In
particular, you should **NOT** do any of the following:

- include Fortran-git-version as a git submodule
- compile or install Fortran-git-version separately from your project
- use Fortran-git-version via CMake's `FetchContent`

Doing any of the above will cause Fortran-git-version to capture _its
own_ version, and not your software's.

You probably want to delete the `.git` directory from your bundled
copy, otherwise git will print a warning like:

```
warning: adding embedded git repository: fortran-git-version
hint: You've added another git repository inside your current repository.
hint: Clones of the outer repository will not contain the contents of
hint: the embedded repository and will not know how to obtain it.
hint: If you meant to add a submodule, use:
hint:
hint:   git submodule add <url> fortran-git-version
hint:
hint: If you added this path by mistake, you can remove it from the
hint: index with:
hint:
hint:   git rm --cached fortran-git-version
hint:
hint: See "git help submodule" for more information.
```

If you see this warning, run:

```bash
$ git rm --cached --force fortran-git-version
$ rm -rf fortran-git-version/.git
$ git add fortran-git-version
```

CMake
-----

You can use Fortran-git-version in your CMake project like so:

```cmake
add_subdirectory(fortran-git-version)

target_link_libraries(<your target> PRIVATE fortran_git::fortran_git)
```

---
title: Input files
subtitle: An overview on stella's input files and some tips on how to use them.
---

# Input files used by stella

Stella uses three different input files, denoted by the following extensions:
<ul><li>`.in` These are the standard input files that stella uses </li>
<li>`.list` A list of input files that stella can run in parallel, provided enough cores. </li>
<li>`.multi` A list of three input files for use with the multiple-flux-tube radial boundary condition. </li></ul>


## The `.in` input file

### `!include`  directive in `.in` files 

When stella processes `.in` files, it recognizes an include directive that takes the form

`!include NEW_FILE.in`

Using this directive, 

Note that the included files can have their own include directive, and so stella handles these directives recursively. 


In practice, the `!include` directive is useful for parameters scans where only one parameter

## The `.list` input file


## The `.multi` input file


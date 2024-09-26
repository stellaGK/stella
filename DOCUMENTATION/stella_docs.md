---
project: stella
project_github: https://github.com/stellaGK/stella/
summary: stella is a flux tube gyrokinetic code for micro-stability and turbulence simulations of strongly magnetised plasma
author: The stella team
author_description:
    stella has been developed by many developers
    See the [citation
    file](https://github.com/stellaGK/stella/blob/master/CITATION.cff)
    for a complete list
src_dir: ..
page_dir: automatic_documentation_manual_pages
exclude_dir: ../EXTERNALS
             ../tests
             ./automatic_documentation_FORD
output_dir: ./automatic_documentation_FORD
predocmark: >
docmark: <
fpp_extensions: fpp
                F90
display: public
         protected
         private
print_creation_date: true
md_extensions: markdown.extensions.toc
               ford.md_striped_table
---

stella solves the gyrokinetic-Poisson system of equations in the local limit
using an operator-split, implicit-explicit numerical scheme. It is capable of
evolving electrostatic fluctuations with fully kinetic electrons and an
arbitrary number of ion species in general magnetic geometry, including
stellarators.

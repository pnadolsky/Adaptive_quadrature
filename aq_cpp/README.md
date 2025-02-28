# Adaptive Quadrature C++ library
This README serves as a directory for this library.  Currently, descriptions of the functionality is split between the files README_QUADRATURE.md, README_MISC.md, and README_ADAPTIVE.md.

*  README_QUADRATURE.md:  This library contains classes the run single quadratures on finctions of one variable.
*  README_ADAPTIVE.md:  This library contains classes that implement the adaptive quadrature tree and the code to read and write the jsons 
*  README_MISC.md:  This library contains miscellaneous functions (like the ploylog integrand) that have various uses

Note that json.hpp is required for weights_loader to load the quadrature roots.  If you have this with your C++ installation then everyhting should compile normally.  if you choose download a copy and use -Iinclude as in the makefile, then put the header here: include/nlohmann and it will work as normal.

The makefile should work on Windows, Cygwin, Linux and MacOS if you have make/mingw32-make installed.  It will generate a 
"build" directory for the object files and a "bin" directory for the executables.  "make clean" will remove these directories.  

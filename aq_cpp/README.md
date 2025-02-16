# Adaptive Quadrature C++ library
This README serves as a directory for this library.  Currently, descriptuions of the functionality is split between the files, README_QUADRATURE.md and README_MISC.md.

Note that json.hpp is required for weights_loader to load the quadrature roots.  Currently, the code use a working directory copy (not included in the repo).  If you have this with your C++ installation then you can change the headers to use that (but please do not merge your local headers back into main or you will break mine).  If I ever get around to loading json.hpp correctly I will update the headers accordingly and make annoucement here.

*  README_QUADRATURE.md:  This library contains classes the run single quadratures on finctions of one variable.
*  README_ADAPTIVE.md:  This library contains classes that implement the adaptive quadrature tree and the code to read and write the jsons 
*  README_MISC.md:  This library contains miscellaneous functions (like the ploylog integrand) that have various uses

# DMphyClus
A R package for phylogenetic clustering

To install DMphyClus, simply download DMphyClus_1.0.tar.gz and run "R CMD INSTALL DMphyClus_1.0.tar.gz". Note that, for now, the GNU Scientific Library (GSL) is required for the package to compile.

The GSL can be found at http://mirror.csclub.uwaterloo.ca/gnu/gsl/ . Installation instructions can be found in the INSTALL file within the downloaded archive. 

IMPORTANT: If you have trouble installing the gsl library in R, it might be that the LD_LIBRARY_PATH environment variable in R is not properly set. 
You can check this with,

Sys.getenv()[["LD_LIBRARY_PATH"]]

Find where  libgsl.so.23.0.0 is located and add the directory with a command like,

Sys.setenv(LD_LIBRARY_PATH = paste(Sys.getenv()[["LD_LIBRARY_PATH"]], "/usr/local/lib", sep = "::"))

You can find the location by copying in your OS console:

find / -name libgsl.so.23.0.0 -printf '%h\n' 2>/dev/null

If this command returns nothing, please make sure the GSL is properly installed.

I strongly suggest reading the vignette before starting to use the package. To do so, start an R session and type "vignette("DMphyClus_Example", package = "DMphyClus")". This should open a browser window and display it. Else, you can simply go to the directory where the package vignette file is found, e.g. ~/R/i686-pc-linux-gnu-library/3.2/DMphyClus/doc/DMphyClusExample, and open DMphyClus_Example.html directly. Be aware that running the example in the vignette can take several minutes!

Note that you will need to have the OpenMP libraries available for the package to compile properly.

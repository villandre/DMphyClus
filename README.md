# DMphyClus
A R package for phylogenetic clustering

To install DMphyClus, simply download DMphyClus_1.0.tar.gz and run "R CMD INSTALL DMphyClus_1.0.tar.gz". Note that, for now, the GNU Scientific Library (GSL) and the Boost C++ libraries are required for the package to compile. 

The GSL can be found at http://mirror.csclub.uwaterloo.ca/gnu/gsl/ . Installation instructions can be found in the INSTALL file within the downloaded archive. The Boost C++ libraries can be found at https://dl.bintray.com/boostorg/release/1.67.0/source/ . Once again, installation instructions can be found in the README file. 

IMPORTANT: You might get error message:
TreeNode.h:7:25: fatal error: gsl/gsl_rng.h: No such file or directory

This means the GSL header files were not found. You can find their location by copying in the console:

find / -name gsl_rng.h -printf '%h\n' 2>/dev/null

If this command returns nothing, please make sure the GSL is properly installed. Otherwise, on line 2 of DMphyClus/src/Makevars, replace

-I/usr/local/include/gsl with -I/path/returned/by/the/find/command

where /path/returned/by/the/find/command is replaced by the output of the find command. 

I strongly suggest reading the vignette before starting to use the package. To do so, start an R session and type "vignette("DMphyClus_Example", package = "DMphyClus")". This should open a browser window and display it. Else, you can simply go to the directory where the package vignette file is found, e.g. ~/R/i686-pc-linux-gnu-library/3.2/DMphyClus/doc/DMphyClusExample, and open DMphyClus_Example.html directly. Be aware that running the example in the vignette can take several minutes!

Note that you will need to have the OpenMP libraries available for the package to compile properly.

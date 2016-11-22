# VIPR
Verifying Integer Programming Results

This repository contains software to check certificates for integer programming results and provides technical specifications for the .vipr certificate file format.

Supplementary information is provided for the article [Verifying Integer Programming Results](http://nbn-resolving.de/urn:nbn:de:0297-zib-61044), ZIB-Report 16-58, November 2016.

Detailed instance-wise results for the computational experiments in the article are provided [here](experiments/README.md).  Links to all associated .vipr certificate files are also provided.

## Software tools

This repository contains three C++ programs:
- [viprchk](code/viprchk.cpp): A program that verifies mixed-integer linear programming certificate files specified in the .vipr file format, described below.
- [viprttn](code/viprchk.cpp): A program that tightens and improves .vipr files, potentially reducing their size and allowing for easier checking.
- [vipr2html](code/viprchk.cpp): A program that converts .vipr certificate files to a human readable HTML format (not recommended for large files.

A [makefile](code/makefile) is provided to build files, they rely on the [GNU GMP library](https://gmplib.org/) for exact rational arithmetic.


## File format specification (.vipr)

A conceptual description of the verified integer programming result (.vipr) file format is given in [Verifying Integer Programming Results](VerifyingIPResults.pdf).  A more detailed technical description is provided [here](http://rawgit.com/ambros-gleixner/VIPR/master/cert_spec_v1_0.html).

A small example is given as [paper_eg3.vipr](code/paper_eg3.vipr), and certificates for large instances can be found [here](experiments/README.md).

## Software for generating (.vipr) certificates

A branch of the [exact rational version](http://scip.zib.de/#exact) of the [SCIP optimization software](http://scip.zib.de) has been developed to produce certificates of MIP results.  It is currently under active development and is available by request from the authors of this project.  It will be included in future releases of the SCIP software.

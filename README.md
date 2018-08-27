# VIPR. Verifying Integer Programming Results

## About

*VIPR* is new a software project to verify, in exact rational arithmetic, the correctness of results computed by mixed-integer linear programming solvers.  It is based on an elementary file format for LP-based branch-and-cut certificates proposed in the article

> Kevin K.H. Cheung, Ambros Gleixner, and Daniel E. Steffy: [Verifying Integer Programming Results](http://dx.doi.org/10.1007/978-3-319-59250-3_13). In: F. Eisenbrand and J. Koenemann, eds., Integer Programming and Combinatorial Optimization: 19th International Conference, IPCO 2017, pp. 148-160, 2017, [`doi:10.1007/978-3-319-59250-3_13`](http://dx.doi.org/10.1007/978-3-319-59250-3_13).

This repository contains a detailed technical [specification of the certificate file format](http://rawgit.com/ambros-gleixner/VIPR/master/cert_spec_v1_0.html), [software](code/) to check, display, and compress certificate files, and [supplementary information](experiments/) on the computational experiments conducted for the article above.


## Software

*VIPR* currently comes with three C++ programs:
- [`viprchk`](code/viprchk.cpp): A program that verifies mixed-integer linear programming certificate files specified in the `.vipr` file format.
- [`viprttn`](code/viprchk.cpp): A program that tightens and improves .vipr files, potentially reducing their size and allowing for easier checking.
- [`vipr2html`](code/viprchk.cpp): A program that converts `.vipr` certificate files to a human readable HTML format (not recommended for large files).

The code should compile with most modern compilers supporting the C++11 standard.  Directory [`code`](code/) provides a basic [`makefile`](code/makefile).  We have successfully built the tools with GNU g++ version 5.4.0 under Ubuntu 16.04.4.  It depends on the [GNU Multiple Precision library](https://gmplib.org/).


## File format specification `.vipr`

A conceptual description of the verified integer programming result (`.vipr`) file format is given in the article [above](http://nbn-resolving.de/urn:nbn:de:0297-zib-61044).  A more detailed technical specification is provided [here](http://rawgit.com/ambros-gleixner/VIPR/master/cert_spec_v1_0.html).

A small example is given as [paper_eg3.vipr](code/paper_eg3.vipr).  Certificates for large MIP instances from the literature can be found as part of the [supplementary information](experiments/) of the article.


## Developers

- [Kevin K.H. Cheung](https://carleton.ca/math/people/kevin-cheung/), School of Mathematics and Statistics, Carleton University
- [Ambros Gleixner](http://www.zib.de/gleixner), Department Optimization, Zuse Institute Berlin
- [Daniel E. Steffy](https://files.oakland.edu/users/steffy/web/), Department of Mathematics and Statistics, Oakland University


## Software for generating `.vipr` certificates

We have created an extension of the [exact rational version](http://scip.zib.de/#exact) of the [SCIP optimization software](http://scip.zib.de) that is able to produce certificates of MIP results.  It is currently under active development and is available by request from the authors of this project.  It will be included in future releases of the SCIP software.

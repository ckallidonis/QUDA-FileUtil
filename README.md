# QUDA File Utilities

### Package for manipulating QUDA output files

This package contains a set of programs and scripts,
particularly useful to manipulate ASCII and HDF5 files.
The repository supports ASCII/HDF5 files containing disconnected
loops, as well as two- and three-point correlation functions.
These ASCII/HDF5 files are expected to be output of
[quda-QKXTM-Multigrid](https://github.com/ETMC-QUDA/quda-QKXTM-Multigrid).

## Package contents

* **extract:**
This folder contains C/C++ programs to extract HDF5 files
into their corresponding ASCII format.

* **convert:**
This folder contains C/C++ programs to convert ASCII files
into the HDF5 format.

* **scripts:**
This folder contains scripts which:
  * run the programs that extract HDF5 files.
  * reshape files already in ASCII format into a different,
    predefined ASCII form to ensure compatibility with the
    extracted HDF5 files.
  * compare two- and three-point functions and disconnected
   quark loops.

## Dependencies
The package requires linking to the HDF5 and FFTW libraries.

## Compilation/Installation
As a rule of thumb, the user should copy [make.inc.SAMPLE](make.inc.SAMPLE) and 
[extract/make_extract.inc.SAMPLE](./extract/make_extract.inc.SAMPLE) to
a name corresponding to the local workstation they are compiling the package.
Then modify the new files to contain the paths to the HDF5 and FFTW libraries
and the MPI on the local workstation and create soft links

`ln -s make.inc.<local> make.inc`

`cd extract`

`ln -s make_extract.inc.<local> make_extract.inc`.

Finally, type `make` on the main directory of the package.

## Author
**Christos Kallidonis**

E-mail: c.kallidonis@cyi.ac.cy

Web: [christoskallidonis.com](http://christoskallidonis.com)

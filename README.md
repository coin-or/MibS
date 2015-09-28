MibS (Mixed Integer Bilevel Solver) 0.5
=======================================

SUPPORTED PLATFORMS

MiBS should work on OS X, Linux, and Windows. The project files for building
on Windows have not been tested in some time, however. the recommended way of
building is to use the Unix-style Makefile. On Windows, this can be done
with the MSys shell. On other platforms, the make command comes pre-installed.

DEPENDENCIES

MibS depends on the CHiPPS (https://projects.coin-or.org/CHiPPS) and SYMPHONY
(https://projects.coin-or.org/SYMPHONY) projects of COIN-OR. First install the
libraries of these two projects and their dependencies. This can be done by
building from source or downloading a binary distribution.

BUILDING

With CHiPPS and SYMPHONY installed, modify the Makefile to reflect the installation 
directory of the above COIN-OR projects (edit the value of the variable '''COININSTDIR'''). 
Then execute the command

'''
make
'''

RUNNING

To solve a bilevel program, you must provide both an MPS file and an auxiliary
information file that specifies which variables and constraints are associated
with the each level. Then call mibs like this:

mibs -Alps_instance file.mps -MibS_auxiliaryInfoFile aux_file.txt

It is also possible to specify additional settings in a parameter file with, e.g.,

mibs -param mibs.par

MibS has many parameters. See the example parameter file '''mibs.par''' and
the header file '''MibParams.h''' for explanations.

HELP

Please post questions and issues to the github project page for MibS.

http://github.com/tkralphs/MibS

Enjoy!
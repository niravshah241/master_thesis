# RBmatlab

Release version 1.16.09, Last change in this file 26. Sept. 2016.

RBmatlab is a MATLAB library for model order reduction with Reduced
Basis Methods for various discretizations types and application
settings.

The package comprises for example:
-Grid management (2D rectangular, triangular, multi-D adaptive 
  nonconforming cubegrid)
-Discretization techniques: Finite Elements, Finite Volumes, 
  Finite Differences, Local Discontinuos Galerkin
  including discrete function types, visualization, etc.
-Problem types: stationary and instationary advection-diffusion, 
  two-phase flow, Stokes, Navier-Stokes
-Demos and example scripts
-Interfaces to external packages (dune, ALBERTA-Grid, Comsol)


## Installation:

After unpacking (or git checkout) of the package, you need to adjust
MATLAB startup script to load the RBmatlab package. Under UNIX-like
environments, the MATLAB startup file is located in
$HOME/matlab/startup.m. This file is run on every start of MATLAB and
can be used to load specific packages such as RBmatlab.

If your MATLAB installation supports the "setenv" command, the installation
is very simple:

Make sure, that the following code is executed during MATLAB startup,
e.g. by putting it into the $HOME/matlab/startup.m file or by
executing it every time you start MATLAB by cut'n paste.  Of course
replace the pathname by your local version. In particular you must set
the path to the RBmatlab directory and choose the temporary-data
directory RBMATLABTEMP, where some GB of space are available. 

    %-------- start of init code ---------
    setenv('RBMATLABHOME','/path/to/your/rbmatlab/directory/');
    setenv('RBMATLABTEMP','/tmp/matlab');
    addpath(getenv('RBMATLABHOME'));
    startup_rbmatlab
    %-------- end of init code ---------

## Basic Overview

The overall idea of the RBmatlab software is to define global interfaces
that can be used for applying the Reduced Basis method.
As a simple example, consider the following code for a simple thermal-block
example:

Define the model. In this example, we use a simple structure-based
definition of the model:
    model = thermalblock_model_struct;

The idea of RBmatlab is to use so-called `model_data` structures, that
contain all "expensive" components, such as the grid or the information
about the discrete functions:
    model_data = gen_model_data(model);

Perform a detailed (expensive, high dimensional) simulation:
    dsim = detailed_simulation(model, model_data);

Construction of the reduced basis and operator components:
    detailed_data = gen_detailed_data(model, model_data);

Construction of the reduced data (for error estimation, ...)
    reduced_data = gen_reduced_data(model, detailed_data);

Performing a reduced basis simulation:
    rbsim = rb_simulation(model, reduced_data);
       
Reconstructing the full solution (if needed)
    rbrec = rb_reconstruction(model, detailed_data, rbsim);

In general, all models in RBmatlab should be callable with the above interface
functions. Please make sure that your model follows this principle by either implementing
the AbstractModel interfaces, or by providing function handels in the structure returned
by your model constructor. The function `gen_model_data(model)` for example is just a wrapper
that calls `model.gen_model_data()`. All other interfaces work in the same manner.

There are two different types of approaches, both are supported and
serve different goals:

-struct-based models: the models and data-instances are realized by
 MATLAB structures. This allows simple "1-file-for-one-model"
 implementation, and is very useful for instructional purposes, and
 suitable for users not familiar with object oriented programming.

-object oriented programming based models: This allows modular
 exchange of subparts of models (e.g. basis generation, etc.) and
 should be used by experienced users.


## RB-Tutorial:

The package is used for RB summerschools and the RB-Tutorial book chapter

  B. Haasdonk: Reduced Basis Methods for Parametrized {PDE}s --  
  A Tutorial Introduction for Stationary and Instationary Problems.
  Chapter in P. Benner, A. Cohen, M. Ohlberger and K. Willcox (eds.): 
  "Model Reduction and Approximation: Theory and Algorithms", SIAM, Philadelphia, 2016

A Preprint version of this chapter is available at:

  http://www.simtech.uni-stuttgart.de/publikationen/prints.php?ID=938

We suggest to read that chapter and "play along" the presented
experiments with the script `rb_tutorial.m`. In particular the commands
`rb_tutorial(1)` until `rb_tutorial(17)` reproduce the experiments/plots
of that tutorial comprising thermal diffusion and advection-diffusion
problems.

As special variant of that MATLAB script we offer
`scripts/rb_tutorial_standalone.m` which is decoupled from RBmatlab in
the sense that only that file together with the file 
`data_rb_tutorial_standalone.mat` is sufficient to reproduce the
experiments for the thermal diffusion in the tutorial.

## Demos:

After starting MATLAB, have a look at the "rbmatlabdemos" (simply call
`rbmatlabdemos` from the MATLAB command line) and start them one by one to
get an impression of the extent of the package. In particular the following
demos are suggested. View their source to see the use of the package's
routines:

1. Example of a detailed simulation (finite volume time evolution)
   `demo_explicit_FV`

2. Steps of model reduction (inspect source-code while executing!)
   `demo_rb_steps_struct`

3. Interactive Gui
   `demo_rb_gui`

Many further demos are available, see the corresponding directory.

## Documentation:

The documentation of the package is realized via Doxygen and a tool
mtoc++,

The documentation of the package is generated by Doxygen - an
automatic documentation tool. Although designed for C-like languages,
Doxygen allows to pre-process the source code by a filter program
that makes MATLAB-code parseable by Doxygen as well.  This generates a
browsable html documentation of the files, classes and
dependencies. The corresponding starting file is
`doxygen/html/index.html`.

In order to extend and re-build this documentation by your own, you
therefore need to download and install a recent version of doxygen
available from its website. and the filter program mtoc++ currently 
available from 

http://www.ians.uni-stuttgart.de/MoRePaS/software/mtocpp/docs/index.html

In order to build the documentation, you need to run the make_docu.sh 
script in the base directory of the RBmatlab package. (Please ignore the
errors and warnings)

As preliminary optional step, in order to obtain a clean configuration file,
you can delete the file `/doxygen/configuration`


## Contact:
The package is jointly developed by the University of Stuttgart, the 
University of Muenster, the University of Ulm.

The package is available via the website http://www.morepas.org

Please send remarks/comments/questions to rbmatlab developers mailing
list: rbmatlab@listserv.uni-stuttgart.de
# master_thesis

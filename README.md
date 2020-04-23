# multilayer
Code for multilayer solver.
## Working example
1. Compile the required fortran library by running
  ```
  cd $BASILISK/ppr/
  make
  ````
2. Under ./working_example, run `sh produce.sh NAME` to generate a folder fNAME that contains the executable. NAME can be either test or test_adaptive (without extension). test is quadtree without adapt_wavelet.
3. Copy ./pre_10layer or ./pre_5layer into the fNAME folder and name it ./pre. This folder contains the synthesized initial wave field and velocity field. Run `mkdir surface` to creat a folder that the program writes output to.
4. Run `./test_adaptive NLAYER MAXLEVEL MINLEVEL ETAE TEND`. The input parameters are number of layers, max level, min level, max error for eta and the ending time respectively.

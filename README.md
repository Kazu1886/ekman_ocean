# ekman_ocean
The model runs on C and requires GSL packages to be installed, which can be found at https://www.gnu.org/software/gsl/ and proper instillation can be done by following https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903#osx
Once the path has been set up, to run the model with the commands:
gcc -Wall -I/PATH_TO_GSL/2.7.1/include -c Ocean_model.c
gcc -L/PATH_TO_GSL/lib Ocean_model.o -lgsl -lgslcblas -lm
./a.out -v 2

where PATH_TO_GSL for me is usr/local/Cellar/gsl/2.7.1
-v is the option for what version of the model to run:
  0: no horizontal transport
  1: difusion only
  2: full model

# CFTR
State estimation for CFTR using the sum-product algorithm

Use Python 3 with the Anaconda Python distribution. (Minimum requirements: numpy, matplotlib)

By default, the simulation lasts 200 time instants. Two plots are generated: the a posteriori probabilities of each state at each time instant, and a comparison of the estimated (max a posteriori) state and the actual kinetic state. Two values are printed: the number of time instants in which the state estimate was correct, and the number in which the state estimate was within one of the correct state.

To run:

from main import main; main()

To change the length of the simulation, for example to 1000 time instants, add "numTimeInstants=1000" as a parameter to main().

To suppress plots, add "plots=False" as a parameter to main(). (In this case, only the two values mentioned above are printed.)

To change the receptor parameters, edit the CFTR.py file.

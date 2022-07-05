# ResearchCode

These projects are the result of work done as an undergrad research student under Professor Richard Scalettar at the University of California Davis.
Both projects are based in the field of condensed matter physics and take two fairly common models, the Ising Model and the Hubbard Model, and work to use
Monte Carlo algorithms as a means of optimization. 

The Ising Model at this point is written in both C++ and python and both codes have slight bugs and different physical setups. One project in particular
takes the standard Ising model and adjusts the coupling parameter linearly with respect to the position within the lattice. i.e. stronger coupling in the
middle than on the ends. The results of this change in coupling strength is then evaluated by looking at the changes in magnetization, binder ration, and 
specific heat of the lattice. 

The Hubbard Model code is used to investigate quantum state transfer in the a 1 dimensional spin chain of length N with m up and l down electrons. 
State transfer is then investigated using time evolution of this quantum system and overall Fidelity, ie. transfer to the final state, is optimized using 
a Monte Carlo algorithm which determines the best hopping values between neighboring sites. Future investigations into extended hubbard are planned for the 
foreseeable future. 

It should also be noted at this time that there is a bug in the Monte Carlo Algorithm for the hubbard model code. Mainly there is a tendency for the system to 
get stuck in a local minumum within the space of possible coupling coefficients instead of finding the true minimum. This might have something to do with 
the step size that I've chosen or with the temperature dependence of the probability distribution that I used. 


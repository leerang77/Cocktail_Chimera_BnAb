This folder contains the code and data used to simulate affinity maturation.

runMultipleGC.m is the execution file, and should be used by passing appropriate arguments.

SharedConstants.m, UseSharedConstants.m, and InitializeSharedConstants.m allows management of the simulation parameters.

Bcell.m and Tcells.m are MATLAB class definitions, used to represent B cells and T cells.

E2AgTable.mat and E2Ag.m are used to determine the amount of antigen captured by B cells based on their affinities.

getAffinityChange.m is used to sample mutations. 

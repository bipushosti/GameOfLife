Basic implementation of Conway's Game of Life. 

Initializes a random grid of desired size and number of iterations with 1, or 0.1 = Occupied 
0 = Unoccupied

Contains the files:
gameOfLife.cu
plot_grid.m

gameOfLife.cu contains the CUDA code to generate a result grid whose size is GridSize x Number of Iterations.

plot_grid.m contains octave code to iterate over and create animation of the result file.



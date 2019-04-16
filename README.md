Basic implementation of Conway's [Game of Life](https://en.wikipedia.org/wiki/Conway's_Game_of_Life) using CUDA.

Initializes a random grid of desired size and number of iterations with 1, or 0, where: 

1 = Occupied, 0 = Unoccupied

The repository contains the following files:

1. gameOfLife.cu  
2. plot_grid.m

gameOfLife.cu contains the CUDA code to generate a result grid whose size is GridSize x Number of Iterations.
   
plot_grid.m contains octave code to iterate over and create animation of the result file.







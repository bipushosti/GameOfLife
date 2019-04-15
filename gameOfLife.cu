

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

static void HandleError( cudaError_t err,const char *file, int line );

__global__ void get_new_status(uint8_t* gridStatus,uint8_t* newStatus, uint32_t gridSize_x, uint32_t gridSize_y, uint32_t iteration_number);
__device__ uint8_t get_neighbour_status(uint32_t xsize, uint32_t ysize,uint32_t xcoord, uint32_t ycoord, uint8_t* gridStatus, uint32_t iteration_number);
__global__ void update_cell_grid(uint8_t* newStatus, uint8_t* gridStatus, uint32_t gridSize_x, uint32_t gridSize_y, uint32_t iteration_number);


//Function that returns an integer array of the status of neighbours when the 
//x,y dimensions of the grid and the current location of cell underconsideration is provided

__device__ uint8_t get_neighbour_status(uint32_t xsize, uint32_t ysize,uint32_t xcoord, uint32_t ycoord, uint8_t* gridStatus, uint32_t iteration_number)
{
	//Sum of the neighbour_status array; Denotes the total alive cells in the neighbours
	uint8_t neighbour_sum = 0;

	//Total size of the cell grid
	uint32_t totalSize = xsize * ysize;


	//Y axis goes from top to bottom
	//Coordinates of neighbours relative to current position
	//The order is right to left and top to bottom excluding the center cell

	int8_t neighbour_x_positions[8] = {-1, 0, 1, -1, 1,-1, 0, 1};
	int8_t neighbour_y_positions[8] = {-1, -1, -1, 0, 0, 1, 1, 1} ;
	
	int neighbour_xcoord, neighbour_ycoord;

	for(int i=0;i<8;i++) {

		neighbour_xcoord = xcoord + (int32_t)neighbour_x_positions[i];
		neighbour_ycoord = ycoord + (int32_t)neighbour_y_positions[i];
	
	//	printf("(nieghbour_xcoord, neighbout_ycoord) :: (%d,%d)\n",neighbour_xcoord, neighbour_ycoord);
	
		//If the neighbour cell is out of bounds assume the cell is dead
		if( (neighbour_xcoord < xsize) && (neighbour_xcoord >= 0) && (neighbour_ycoord < ysize) && (neighbour_ycoord >= 0) ){
			//Getting the status of the neighbour that is still inside the grid
			neighbour_sum = neighbour_sum + gridStatus[(iteration_number-1) * totalSize + neighbour_ycoord * xsize + neighbour_xcoord];
		}
		else {
			neighbour_sum += 0;
		}
	}

	return neighbour_sum;

}

__global__ void get_new_status(uint8_t* gridStatus,uint8_t* newStatus, uint32_t gridSize_x, uint32_t gridSize_y, uint32_t iteration_number)
{
	int x_id = blockIdx.x * blockDim.x + threadIdx.x; 
	int y_id = blockIdx.y * blockDim.y + threadIdx.y;	

	//Total Grid Size
	uint32_t totalSize = gridSize_x * gridSize_y;

	//For each thread/grid box use rules to decide what to do
	//Rules are:
		//Any live cell with fewer than two live neighbors dies, as if by underpopulation.
		//Any live cell with two or three live neighbors lives on to the next generation.
		//Any live cell with more than three live neighbors dies, as if by overpopulation.
		//Any dead cell with exactly three live neighbors becomes a live cell, as if by reproduction.
	
	
	if( (x_id < gridSize_x) && (y_id < gridSize_y) ) {
	
		//Status of the cell in question
		uint8_t cell_status = gridStatus[(iteration_number-1) * totalSize + y_id * gridSize_x + x_id];

		//Getting the total number of alive neighbours
		uint8_t neighbour_status = get_neighbour_status(gridSize_x, gridSize_y, x_id, y_id, gridStatus, iteration_number);

		if((cell_status == 0) && (neighbour_status==3)) {			

		/*	if(neighbour_status == 3){
				newStatus[y_id * gridSize_x + x_id] = 1;
			}
			else{
				newStatus[y_id * gridSize_x + x_id] = 0;
			}
		*/
			newStatus[y_id * gridSize_x + x_id] = 1;

		}
		else if((cell_status == 1) && ((neighbour_status == 2) || (neighbour_status == 3) )) {
		/*	if(neighbour_status  < 2) {				
				newStatus[y_id * gridSize_x + x_id] = 0;
			}
			else if((neighbour_status  == 2) || (neighbour_status == 3)){
				newStatus[y_id * gridSize_x + x_id] = 1;
			}
			if(neighbour_status  > 3){
				newStatus[y_id * gridSize_x + x_id] = 0;
			}
		*/

			newStatus[y_id * gridSize_x + x_id] = 1;			
		}
					
	}
}

//Function that copies oves the new grid values into the original one
__global__ void update_cell_grid(uint8_t* newStatus, uint8_t* gridStatus, uint32_t gridSize_x, uint32_t gridSize_y, uint32_t iteration_number)
{

	int x_id = blockIdx.x * blockDim.x + threadIdx.x; 
	int y_id = blockIdx.y * blockDim.y + threadIdx.y;	
	//Total Grid Size
	uint32_t totalSize = gridSize_x * gridSize_y;
	
	if( (x_id < gridSize_x) && (y_id < gridSize_y) ) {
		gridStatus[iteration_number * totalSize + y_id * gridSize_x + x_id] = newStatus[ y_id * gridSize_x + x_id];
		newStatus[y_id * gridSize_x + x_id] = 0;
	}
}



static void HandleError( cudaError_t err,const char *file, int line ) {
    if (err != cudaSuccess) {
  		printf( "%s in %s at line %d\n", cudaGetErrorString( err ),file, line );
        exit( EXIT_FAILURE );
    }
}


int main()
{

	int number_of_iterations = 1000;
	uint32_t cellGridSize_x = 100;
	uint32_t cellGridSize_y = 100;
	uint32_t totalCellGridSize = cellGridSize_x * cellGridSize_y;

	uint8_t* cell_grids = (uint8_t*)malloc(totalCellGridSize * number_of_iterations * sizeof(uint8_t));
	
	//-----------------------------------------//

	uint8_t* d_cell_grids;
	uint8_t* d_new_positions;
	HANDLE_ERROR(cudaMalloc((void**)&d_cell_grids, totalCellGridSize * number_of_iterations * sizeof(uint8_t)));
	HANDLE_ERROR(cudaMalloc((void**)&d_new_positions, totalCellGridSize * sizeof(uint8_t)));
	//-----------------------------------------//

	srand(time(0));
	
	//Generating random bits for the first cell grid
	for(uint32_t i = 0; i<totalCellGridSize; i++) {
		cell_grids[i] = rand() % 2;
	}

	//-----------------------------------------//
	//Copying the first grid from host to device
	HANDLE_ERROR(cudaMemcpy(d_cell_grids, cell_grids, totalCellGridSize * sizeof(uint8_t),cudaMemcpyHostToDevice));
	//-----------------------------------------//
	printf("xsize: %d, ysize: %d\n",(cellGridSize_x + 31)/32,(cellGridSize_y + 31)/32);

	//Setting the block and grid size; Setting the block size to have 32 threads
	dim3 gridSize((cellGridSize_x + 31)/32,cellGridSize_y,1);
	dim3 blockSize(32,1,1);

	for(uint32_t i=1; i<number_of_iterations; i++) {
		get_new_status<<<gridSize,blockSize>>>(d_cell_grids,d_new_positions,cellGridSize_x, cellGridSize_y,i);
		HANDLE_ERROR(cudaDeviceSynchronize());	
		update_cell_grid<<<gridSize,blockSize>>>(d_new_positions, d_cell_grids, cellGridSize_x, cellGridSize_y, i);
		HANDLE_ERROR(cudaDeviceSynchronize());	
	}
	//-----------------------------------------//
	HANDLE_ERROR(cudaMemcpy(cell_grids, d_cell_grids, number_of_iterations * totalCellGridSize * sizeof(uint8_t), cudaMemcpyDeviceToHost));
	//-----------------------------------------//
	//Putting the data copied from the Device into a file
	FILE * fp;
	fp = fopen("result.txt","w+");

	for(int h=0; h < number_of_iterations; h++) {
		for(int i=0; i<cellGridSize_y; i++) {
			for(int j=0; j<cellGridSize_x; j++) {
				fprintf(fp,"%d ",cell_grids[h * totalCellGridSize + i * cellGridSize_x + j]);	
			}
			fprintf(fp,"\n");	
		}
	}
	//-----------------------------------------//
	fclose(fp);
	HANDLE_ERROR(cudaFree(d_cell_grids));
	HANDLE_ERROR(cudaFree(d_new_positions));
	free(cell_grids);
	return 0;
}

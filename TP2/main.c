#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define TRUE 1
#define FALSE 0

int main( int argc, char **argv ) {
  int rank, grid_rank, col_rank, size, grid_size, col_size;
  MPI_Comm grid_2D, column;
  int dims[2], periods[2], remain_dims[2], grid_coords[2], col_coords[1], reorder;
   
  MPI_Status status; 
  MPI_Init( NULL, NULL ); 
  MPI_Comm_rank( MPI_COMM_WORLD, &rank); 
  MPI_Comm_size( MPI_COMM_WORLD, &size);

  dims[0] = 3; 
  dims[1] = 4; 
  periods[0] = FALSE; 
  periods[1] = TRUE; 
  reorder = 1;

  remain_dims[0] = 1;
  remain_dims[1] = 0;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &grid_2D);
  MPI_Comm_rank(grid_2D, &grid_rank);
  MPI_Comm_size(grid_2D, &grid_size);

  MPI_Cart_sub(grid_2D, remain_dims, &column); 
  MPI_Comm_rank(column, &col_rank);
  MPI_Comm_size(column, &col_size);

  MPI_Cart_coords(grid_2D, grid_rank, 2, grid_coords);
  MPI_Cart_coords(column, col_rank, 1, col_coords);

  int W = 0;
  int V[3];

  if (col_rank == 1) {
    W = rank;
    V[0] = rank;
    V[1] = rank;
    V[2] = rank;
  }
  
  // Envoi de V aux voisins avec le rang de la deuxième cellule 
  //MPI_Scatter(V, 1, MPI_INT, &W, 1, MPI_INT, 1, column);

  // Communication avec les voisins (dessus et dessous)
  if (col_rank == 1) {
    int rank_top, rank_bottom;
    MPI_Cart_shift(column, 0, 1, &rank_bottom, &rank_top);
    MPI_Send(&W, 1, MPI_INT, rank_bottom, 99, column);
    MPI_Send(&W, 1, MPI_INT, rank_top, 99, column); 
  } else {
    MPI_Recv(&W, 1, MPI_INT, 1, 99, column, &status);
  }

  printf("Proc %d is %d/%d (%d;%d) in grid and %d/%d (%d) in col - W : %d\n", rank, grid_rank, grid_size-1, grid_coords[0], grid_coords[1], col_rank, col_size-1, col_coords[0], W);  

  // MPI_Send(&value, 1, MPI_INT, 0, 99, MPI_COMM_WORLD);
  // MPI_Recv(&value, 1, MPI_INT, 8, 99, MPI_COMM_WORLD, &status);


  MPI_Finalize();
} 

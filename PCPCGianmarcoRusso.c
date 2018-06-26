/*
 ============================================================================
 Name        : PCPCGianmarcoRusso.c
 Author      : Gianmarco Russo
 Version     :
 Copyright   : copyright GianmarcoRusso
 Description : Jacobi iteration
 ============================================================================
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

#define MAXSTEPS 100
#define DIFFNORMLIMIT 0.01f

int main(int argc, char* argv[]){
	int  my_rank; /* rank of process */
	int  p;       /* number of processes */

	struct timeval stop, start;

	int steps = 0;
	int n=50000,m=1000;
	float *x, *xnew;

	int subMsize = 0;
	float *xFinal;

	double diffnorm = 1.0f, diffnormTot = 1.0f;

	/* start up MPI */
	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &p); 

	//dimensione da linea di comando per facilitare test weak scalability
	if(argc > 1) n=atoi(argv[1]);

	int r = n%p;
	int size = my_rank<r?n/p+1:n/p;

	x = malloc((size+2)*m*sizeof(float));
	xnew = malloc((size+2)*m*sizeof(float));

	/*
	 * offset to initialize sub-matrices with
	 * different random values f the same sequence (like a only big matrix)
	 */
	int offset = my_rank<r?my_rank*size:my_rank*size+r;

	int randVal = 0;
	srand(1);
	for(int i = 1; i <= n; i++){
		for(int j = 0; j < m; j++){
			randVal = rand() % 10;
			if(i>offset && i<=offset+size){
				x[(i-offset)*m+j] = randVal;
			}else{
				//do nothing
			}
		}
	}

	if(my_rank==0){
		/*
		 * master allocates a big matrix to
		 * collect results of sub-matrices
		 * received from slaves (plus his) to give
		 * the complete solution at the end
		 */

		xFinal = malloc(n*m*sizeof(float));
	}

	gettimeofday(&start, NULL);

	while(diffnormTot > DIFFNORMLIMIT && steps < MAXSTEPS){
		if (my_rank == 0){
			/* send last row */
			MPI_Send(&x[size*m], m, MPI_FLOAT, my_rank+1, 0, MPI_COMM_WORLD);

			/* receive first row of successor */
			MPI_Recv(&x[(size+1)*m], m, MPI_FLOAT, my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		} else if(my_rank == p-1){
			/* receive last row of predecessor */
			MPI_Recv(&x[0], m, MPI_FLOAT, my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			/* send first row */
			MPI_Send(&x[1*m], m, MPI_FLOAT, my_rank-1, 0, MPI_COMM_WORLD);

		}
		else{
			/* receive last row of predecessor */
			MPI_Recv(&x[0], m, MPI_FLOAT, my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			/* send first row */
			MPI_Send(&x[1*m], m, MPI_FLOAT, my_rank-1, 0, MPI_COMM_WORLD);

			/* send last row */
			MPI_Send(&x[size*m], m, MPI_FLOAT, my_rank+1, 0, MPI_COMM_WORLD);

			/* receive first row of successor */
			MPI_Recv(&x[(size+1)*m], m, MPI_FLOAT, my_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		}

		for(int i = 1; i <= size; i++){
			if(  (my_rank == 0 && i == 1) || (my_rank == p-1 && i == size) ){
				// 0 doesn't change first row an p-1 doesn't change last row
			}else{
				for(int j = 1; j < m-1; j++){
					xnew[i*m+j] = (x[(i+1)*m+j] + x[(i-1)*m+j] + x[i*m+j+1] + x[i*m+j-1])/4;
				}
			}
		}

		diffnorm = 0;
		for(int i = 1; i <= size; i++){
			if(  (my_rank == 0 && i == 1) || (my_rank == p-1 && i == size) ){
				// 0 doesn't calculate diffnorm on first row an p-1 doesn't c. d. on last row
			}else{
				for(int j = 1; j < m-1; j++){
					diffnorm += (xnew[i*m+j] - x[i*m+j]) * (xnew[i*m+j] - x[i*m+j]);
				}
			}
		}

		for(int i = 1; i <= size; i++){
			if(  (my_rank == 0 && i == 1) || (my_rank == p-1 && i == size) ){
				// 0 doesn't copy first row an p-1 doesn't copy last row
			}else{
				for(int j = 1; j < m-1; j++){
					x[i*m+j] = xnew[i*m+j];
				}
			}
		}

		MPI_Reduce(&diffnorm, &diffnormTot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if(my_rank == 0) diffnormTot = sqrt(diffnormTot);

		MPI_Bcast(&diffnormTot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		steps++;
	}

	if(my_rank != 0){
		//send submatrix
		MPI_Send(&x[m], size*m, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
	}else{
		//master collect submatrices in a big matrix
		//his submatrix
		for(int i = 1; i <= size; i++){
			for(int j = 0; j < m; j++){
				xFinal[(i-1)*m+j] = x[i*m+j];
			}
		}
		//other's submatrices
		int index = size*m;

		for(int i = 1; i < p; i++){
			//receive submatrix and store in the right position of big matrix

			if(i<r || r == 0) subMsize = size*m;
			else subMsize = (size -1)*m;
			//printf("%d: %d\n",i,subMsize);
			MPI_Recv(&xFinal[index], subMsize, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			index += subMsize;
		}

		gettimeofday(&stop, NULL);

		/* print xFinal matrix (commented because take too long to print)
			
 		for(int i = 0; i < n; i++){
			for(int j = 0; j < m; j++){
				printf("%.2f ",xFinal[i*m+j]);
			}
			printf("\n");
		}
		printf("\n");
		
		*/
	}

	if(my_rank == 0){
		double msec = ((stop.tv_sec - start.tv_sec) * 1000.0)
		            				+ ((stop.tv_usec - start.tv_usec) / 1000.0);
		printf("diffnormTot: %.4f, steps: %d, time: %f ms\n",diffnormTot,steps,msec);
	}

	/* shut down MPI */
	MPI_Finalize();


	return 0;
}


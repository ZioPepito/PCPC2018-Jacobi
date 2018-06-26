#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define MAXSTEPS 100
#define DIFFNORMLIMIT 0.01f

int main(int argc, char* argv[]){
	struct timeval stop, start;
	int steps = 0;
	int n=2000,m=1000;
	float *x, *xnew;

	x = malloc(n*m*sizeof(float));
	xnew = malloc(n*m*sizeof(float));

	srand(1);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			x[i*m+j] = rand() % 10;
		}
	}

	/* print initial matrix
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			printf("%.2f ",x[i*m+j]);
		}
		printf("\n");
	}
	printf("\n");*/

	double diffnorm = 1;

	gettimeofday(&start, NULL);

	while(diffnorm > DIFFNORMLIMIT && steps < MAXSTEPS){
		for(int i = 1; i < n-1; i++){
			for(int j = 1; j < m-1; j++){
				xnew[i*m+j] = (x[(i+1)*m+j] + x[(i-1)*m+j] + x[i*m+j+1] + x[i*m+j-1])/4;
			}
		}

		diffnorm = 0;

		for(int i = 1; i < n-1; i++){
			for(int j = 1; j < m-1; j++){
				diffnorm += (xnew[i*m+j] - x[i*m+j]) * (xnew[i*m+j] - x[i*m+j]);
			}
		}
		
		diffnorm = sqrt(diffnorm);

		for(int i = 1; i < n-1; i++){
			for(int j = 1; j < m-1; j++){
				x[i*m+j] = xnew[i*m+j];
			}
		}

		steps++;
	}

	gettimeofday(&stop, NULL);

	/* print final matrix
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			printf("%.2f ",x[i*m+j]);
		}
		printf("\n");
	}*/

	double msec = ((stop.tv_sec - start.tv_sec) * 1000.0)
            + ((stop.tv_usec - start.tv_usec) / 1000.0);

	printf("\nPassi: %d, diffnorm: %.4f, tempo: %f ms \n",steps,diffnorm,msec);
}

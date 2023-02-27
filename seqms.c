#include <stdlib.h> 
#include <stdio.h>
#include <mpi.h>
#include <complex.h>

 
#define WIDTH 800
#define HEIGHT 800
#define MAX_ITER 1000

 
int main(int argc, char **argv) {

    int rank, size, i, j, k, count;
    double start, end, total_time, comm_time, comp_time, ratio;
    double x_min = -2.0, x_max = 2.0, y_min = -2.0, y_max = 2.0;
    double x, y;
    int *local_counts, *global_counts;
    complex double c, z;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    local_counts = (int *)malloc(sizeof(int) * HEIGHT / size);
    global_counts = (int *)malloc(sizeof(int) * WIDTH * HEIGHT);
    start = MPI_Wtime();
    
    for (i = rank * HEIGHT / size; i < (rank + 1) * HEIGHT / size; i++) {
        for (j = 0; j < WIDTH; j++) {
            c = ((x_max - x_min) * j / (double)WIDTH + x_min) + ((y_max - y_min) * i / (double)HEIGHT + y_min) * I;
            z = 0;
            count = 0;
            while (cabs(z) < 2 && count < MAX_ITER) {
                z = z * z + c;
                count++;
            }
            local_counts[i - rank * HEIGHT / size] = count;
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);
    comm_time = MPI_Wtime() - start;
    MPI_Gather(local_counts, HEIGHT / size, MPI_INT, global_counts, HEIGHT / size, MPI_INT, 0, MPI_COMM_WORLD);
    end = MPI_Wtime();
    total_time = end - start;
    comp_time = total_time - comm_time;
    ratio = comp_time / comm_time;
    if (rank == 0) {
        printf("Communication time: %f seconds\n", comm_time);
        printf("Computation time: %f seconds\n", comp_time);
       printf("Communication-to-computation ratio: %f\n", ratio);
        FILE *fp;
        fp = fopen("mandelbrot2.ppm", "wb");
        fprintf(fp, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
        for (i = 0; i < HEIGHT; i++) {
            for (j = 0; j < WIDTH; j++) {
                k = global_counts[i * WIDTH + j];
                if (k == MAX_ITER) {
                    fputc(0, fp);
                    fputc(0, fp);
                    fputc(0, fp);
                } else {
                    fputc((k * 23) % 256, fp);
                    fputc((k * 41) % 256, fp);
                    fputc((k * 67) % 256, fp);
                }
            }
        }
        fclose(fp);
    }

    MPI_Finalize();
    free(local_counts);
    free(global_counts);
    
    return 0;

}

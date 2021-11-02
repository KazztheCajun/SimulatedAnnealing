#include "BenchmarkPSA.h"
#include "BenchmarkSSA.h"

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <errno.h>

#define MAX_LOOPS 100
#define MAX_ANNEALING_STEPS 200
#define MAX_HEAT 100
#define MAX_REPEATED_VALUE 200
#define MAX_ALPHA 0.99
#define MAX_THREADS 12

void main()
{
    int i, j, thread, annealing = 0;
    double alpha, heat = 0.0;
/*    for(thread = 2; thread <= 12; thread += 2)
    {
        printf("Beginning Thread Test #%d!\n", thread);
        //heat = i*25;
        char *s = (char*)malloc(34 * sizeof(char));
        sprintf(s, "./results/Thread_Benchmark#%d.txt", thread);
        //printf("%s\n", s);
        FILE * f = fopen(s, "w+");
        if (!f)
        {
            printf("Error opening/creating \"%s\": %d\n Exiting Benchmark!\n", s, errno);
            return;
        }
        else
        {
            printf("File \"%s\" opened sucessfully\n", s);
        }
        benchmarkPSA(100, 200, .9, thread, f); 
        free(s);
        fclose(f);
    }
*/    
    printf("%s\n", "Starting Serial Simulated Annealing Benchmark!");
    char *s = "./results/SearialSA_Benchmark.txt";
    FILE * f = fopen(s, "w+");
    if (!f)
    {
        printf("Error opening/creating \"%s\": %d\n Exiting Benchmark!\n", s, errno);
        return;
    }
    else
    {
        printf("File \"%s\" opened sucessfully\n", s);
    }
    benchmarkSSA(MAX_HEAT, MAX_ANNEALING_STEPS, MAX_ALPHA, f); 
    fclose(f);
    
}
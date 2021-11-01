#include "BenchmarkPSA.h"

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <errno.h>

#define MAX_LOOPS 100
#define STRING_LENGTH 200
//#define ANNEALING_STEPS 100
//#define STARTING_HEAT 100
#define MAX_REPEATED_VALUE 200
//#define ALPHA 0.9
//#define MAX_THREADS 10
//#define SA_CHUNK_SIZE 10
#define _USE_MATH_DEFINES 




// generate a random floating point number from min to max
double rand_double(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

// Rastrigin function as fitness function
double rastrigin(double *s)
{
    int i;
    double sum = 0.0;
    for(i = 0; i < STRING_LENGTH; i++)
    {
        sum += 10 + ((s[i]*s[i]) - 10*cos(M_TWOPI*s[i]));
    }
    return sum;
}


// Generate a new, unique string of numbers
void new_string(double *c)
{
    int i;
    for (i = 0; i < STRING_LENGTH; i++)
    {
        c[i] = rand_double(-5.12, 5.12); // generate random value
        //c[i] = 5.0; // generate debug values
        //c[i] = 0.0;
    }
}

// Generate a neighborhood for a given string
void get_neighbors(double *n, double *c)
{
    int i, j, index = 0;
    for (i = 0; i < STRING_LENGTH; i++) // for each column in neighborhood
    {
        for (j = 0; j < STRING_LENGTH; j++) // for each row in neighbor string
        {
            if (index == j) // if change index == row in neighbor string, generate new number
            {
                n[i*STRING_LENGTH + j] = rand_double(-5.12, 5.12);
            }
            else // otherwise copy current value
            {
                n[i*STRING_LENGTH + j] = c[j];
            }
        }
        index++;
    }
}

void copyString(double *this, double *that) // copies the info from this into that
{
    int i;
    for(i = 0; i < STRING_LENGTH; i++)
    {
        that[i] = this[i];
    }
}

int isSameString(double *a, double *b) // returns 1 if strings contain the same values, returns 0 if not
{
    int i = 0;
    while(i < STRING_LENGTH)
    {
        if(a[i] != b[i])
        {
            return 0;
        }
        i++;
    }
    return 1;
}

void printString(double *s)
{
    int i;
    for(i=0;i<STRING_LENGTH-1;i++)
    {
        printf("%.10f, ", s[i]); 
    }
    printf("%.10f\n", s[STRING_LENGTH - 1]);
}


void simulated_annealing(double *min_strings, double *origin_strings, double *local_minima, double *timeData, int *loopNum, int loop, int STARTING_HEAT, int ANNEALING_STEPS, int ALPHA)
{
    
    // initialize time, temp, repeat counter, and loop counter
    int time, repeat, step, i , j, index = 0;
    // initialize the current value and new value
    double vc, vn, tdata, temp = 0.0;

    double *neighborhood = (double *)malloc(sizeof(double)*STRING_LENGTH*STRING_LENGTH); // allocate memory for the neighborhood strings array: [number of loops][length of string of doubles]
    double *current = (double *)malloc(sizeof(double)*STRING_LENGTH); // allocate memory for the best strings array: [length of string of doubles]
    double *last = (double *)malloc(sizeof(double)*STRING_LENGTH); // allocate memory for the current array: [length of string of doubles]
    double *neighbor = (double *)malloc(sizeof(double)*STRING_LENGTH); // allocate memory for the best strings array: [length of string of doubles]
    double *best = (double *)malloc(sizeof(double)*STRING_LENGTH); // allocate memory for the current array: [length of string of doubles]

    // select a starting point
    new_string(current);
    for(i=0;i<STRING_LENGTH;i++) // save origin string for each loop
    {
        origin_strings[loop*MAX_LOOPS + i] = current[i];
    }
    tdata = omp_get_wtime(); // save start time
    time = 1; // initialize the time counter
    temp = STARTING_HEAT; // initialize algorithm temp
    repeat = 0; // initialize repeated value counter
    // first loop  that loops until the algorithm has suffeciently cooled
    while(temp > 1 || repeat < MAX_REPEATED_VALUE)
    {
        //printf("\nRestart # %d ", time);
        step = 0; // reinitialize the step counter
        // second loop that tracks the loops during an annealing phase
        while(step < ANNEALING_STEPS)
        {
            //printf("Loop %d\n", step + 1);
            // generate a neighborhood of the current value
            get_neighbors(neighborhood, current);
            // save the last current value for comparison
            copyString(current, last);
            //printf("Current: ");
            //printString(current);
            // select a neighboor at random
            index = (rand() % STRING_LENGTH); // get random index number
            for(i = 0; i < STRING_LENGTH; i++)
            {
                neighbor[i] = neighborhood[index*MAX_LOOPS + i];
            }
            //printf("Neighbor: ");
            //printString(neighbor);
            // calculate the value of the current (vc) and neighbor (vn)
            vc, vn = 0;
            vc = rastrigin(current);
            vn = rastrigin(neighbor);
            //printf("Current best: %.10f\n", vc);
            //printf("Current new: %.10f\n", vn);
            // if neighboor is better than current, set current = neighbor
            if (vc > vn) // rastrigin is a minimization function, so check if vc is greater than vn
            {
                copyString(neighbor, current);
                //vc = vn;
                repeat = 0;
                //printf("New best found: %.10f\n", vc);
            }
            // else if random(0-.99) < e^(vc-vn/temp), set vc = vn
            else if (rand_double(0.0, 0.999999999999999) < exp(-(vn-vc)/temp))
            {
                copyString(neighbor, current);
                //vc = vn;
                repeat = 0;
                //printf("Annealling jump to: %.10f\n", vc);
            }
            else
            {
                //printf("No change in current\n");
            }
            // otherwise do nothing and try a new neighboor until annealing phase is done
            step++;
        }
        if (isSameString(current, last))
        {
            repeat++;
            //printf("Repeated value %.10f %d times\n", vc, repeat);
        }
        //printf("End of restart #%d\n", time);
        temp = (ALPHA * temp) / log(1 + time); // modified logarithmic decrease
        time++; // incriment time
        
    }
    for(i=0;i<STRING_LENGTH;i++) // copy local min string
    {
        min_strings[loop*MAX_LOOPS + i] = current[i];
    }
    local_minima[loop] = rastrigin(current); // save evaluation
    timeData[loop] = omp_get_wtime() - tdata; // save SA thread run time
    loopNum[loop] = time; // save total iterations of SA for thread
    free(neighborhood);
    free(current);
    free(last);
    free(neighbor);
    free(best);
}

// add test suite

// once algorithm is done, print simple results to console, and detailed stats to a document

void benchmarkPSA(int STARTING_HEAT, int ANNEALING_STEPS, double ALPHA, int MAX_THREADS, FILE *output)
{
    // initialize random
    time_t ts; 
    srand(time(&ts));
    //FILE *output;
    int i, j, l, threads = 0;
    double sum, total = 0.0;
    double min = 200000.0;
    double *ms = (double *)malloc(sizeof(double)*MAX_LOOPS*STRING_LENGTH); // allocate memory for the minimum strings array: [number of loops][length of string of doubles]
    double *os = (double *)malloc(sizeof(double)*MAX_LOOPS*STRING_LENGTH); // allocate memory for the origin strings array: [number of loops][length of string of doubles]
    double *lm = (double *)malloc(sizeof(double)*MAX_LOOPS); // allocate memory for the results of the simulated annealing run: [number of loops]
    double *rt = (double *)malloc(sizeof(double)*MAX_LOOPS); // allocate memory for the time results of each run of the sa algorithm: [number of loops]
    double *tt = (double *)malloc(sizeof(double)*MAX_LOOPS); // allocate memeory for the total time of each simulated annealing: [number of loops]
    int *tl = (int *)malloc(sizeof(int)*MAX_LOOPS); // allocate memeory for the total iterations in a simulated annealing run: [number of loops]
    omp_set_num_threads(MAX_THREADS);
    total = omp_get_wtime();
    printf("Beginning Simulated Annealing with %d heat, %d steps, %.2f alpha, and %d threads!\n", STARTING_HEAT, ANNEALING_STEPS, ALPHA, MAX_THREADS);
    #pragma omp parallel for schedule(dynamic, (MAX_LOOPS / MAX_THREADS))
    for(i = 0; i < MAX_LOOPS; i++)
    {
        if (i == 0)
        {
            threads = omp_get_num_threads();
        }
        
        simulated_annealing(ms, os, lm, rt, tl, i, 
                            STARTING_HEAT, 
                            ANNEALING_STEPS, 
                            ALPHA);
    }

    total = omp_get_wtime() - total;
    for(i=0;i<MAX_LOOPS;i++)
    {
        sum += rt[i];
        if(min > lm[i])
        {
            min = lm[i];
        }
        //printf("Local Min: %.5f\n", lm[i]);
    }
    /*
    printf("Loop %d:\nOrigin String: ", l + 1);
    for(i=0;i<STRING_LENGTH-1;i++)
    {
        printf("%.5f, ", os[l*MAX_LOOPS + i]); 
    }
    printf("%.5f\n", os[l*MAX_LOOPS + (STRING_LENGTH - 1)]);
    printf("Local Min String: ");
    for(i=0;i<STRING_LENGTH-1;i++)
    {
        printf("%.5f, ", ms[l*MAX_LOOPS + i]);
    }
    printf("%.5f\n", ms[l*MAX_LOOPS + (STRING_LENGTH - 1)]); */
    fprintf(output, "Local Mins:\n");
    for(i=0;i<MAX_LOOPS;i++)
    {
        fprintf(output, "%.5f\n", lm[i]);
    }
    fprintf(output, "\nTotal Execution time: %.5f seconds\n\nIndividual PSA execution times:\n", total);
    for(l = 0; l < MAX_LOOPS; l++)
    {
        fprintf(output, "%.5f\n", rt[l]);
    }
    fprintf(output, "\n# of Threads: %d\nStarting Heat: %.2f\nAlpha: %.2f", threads, STARTING_HEAT, ALPHA);
    fprintf(output, "\n# of Iterations: \n");
    for(l = 0; l < MAX_LOOPS; l++)
    {
        fprintf(output, "%d\n", tl[l]);
    }
   
    sum /= MAX_LOOPS;
    printf("\n\nThe global min is: %.5f.\nTotal execution time was %.5f seconds, with %.5f seconds on average per Simulated Annealing run.\nThe algorithm was run on %d threads\n", min, total, sum, threads);
    free(ms);
    free(os);
    free(lm);
    fclose(output);
}
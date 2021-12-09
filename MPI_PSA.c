#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <errno.h>


#define MAX_REPEATED_VALUE 200
#define _USE_MATH_DEFINES 




// generate a random floating point number from min to max
double rand_double(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

// Rastrigin function as fitness function
double rastrigin(double *s, int length)
{
    int i;
    double sum = 0.0;
    for(i = 0; i < length; i++)
    {
        sum += 10 + ((s[i]*s[i]) - 10*cos(2*M_PI*s[i]));
    }
    return sum;
}


// Generate a new, unique string of numbers
void new_string(double *c, int length)
{
    int i;
    for (i = 0; i < length; i++)
    {
        c[i] = rand_double(-5.12, 5.12); // generate random value
        //c[i] = 5.0; // generate debug values
        //c[i] = 0.0;
    }
}

// Generate a neighborhood for a given string
void get_neighbors(double *n, double *c, int length)
{
    int i, j, index = 0;
    for (i = 0; i < length; i++) // for each column in neighborhood
    {
        for (j = 0; j < length; j++) // for each row in neighbor string
        {
            if (index == j) // if change index == row in neighbor string, generate new number
            {
                n[i*length + j] = rand_double(-5.12, 5.12);
            }
            else // otherwise copy current value
            {
                n[i*length + j] = c[j];
            }
        }
        index++;
    }
}

void copyString(double *this, double *that, int length) // copies the info from this into that
{
    int i;
    for(i = 0; i < length; i++)
    {
        that[i] = this[i];
    }
}

int isSameString(double *a, double *b, int length) // returns 1 if strings contain the same values, returns 0 if not
{
    int i = 0;
    while(i < length)
    {
        if(a[i] != b[i])
        {
            return 0;
        }
        i++;
    }
    return 1;
}

void printString(double *s, int length)
{
    int i;
    for(i=0;i<length-1;i++)
    {
        printf("%.10f, ", s[i]); 
    }
    printf("%.10f\n", s[length - 1]);
}


void simulated_annealing(double *min_strings, double *origin_strings, double *local_minima, int *loopNum, int loop, int max_loops, int string_length, int core, double STARTING_HEAT, int ANNEALING_STEPS, int ALPHA)
{
    
    // initialize time, temp, repeat counter, and loop counter
    int time, repeat, step, i , j, index = 0;
    // initialize the current value and new value
    double vc, vn, temp = 0.0;

    double *neighborhood = (double *)malloc(sizeof(double)*string_length*string_length); // allocate memory for the neighborhood strings array: [number of loops][length of string of doubles]
    double *current = (double *)malloc(sizeof(double)*string_length); // allocate memory for the best strings array: [length of string of doubles]
    double *last = (double *)malloc(sizeof(double)*string_length); // allocate memory for the current array: [length of string of doubles]
    double *neighbor = (double *)malloc(sizeof(double)*string_length); // allocate memory for the best strings array: [length of string of doubles]
    double *best = (double *)malloc(sizeof(double)*string_length); // allocate memory for the current array: [length of string of doubles]
    struct timeval start_run_time, stop_run_time;
    new_string(current, string_length); // select a starting point
    for(i=0;i<string_length;i++) // save origin string for each loop
    {
        origin_strings[loop*max_loops + i] = current[i];
    }
    // save start time
    //gettimeofday(&start_run_time, NULL);
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
            get_neighbors(neighborhood, current, string_length);
            // save the last current value for comparison
            copyString(current, last, string_length);
            //printf("Current: ");
            //printString(current);
            // select a neighboor at random
            index = (rand() % string_length); // get random index number
            for(i = 0; i < string_length; i++)
            {
                neighbor[i] = neighborhood[index*max_loops + i];
            }
            //printf("Neighbor: ");
            //printString(neighbor);
            // calculate the value of the current (vc) and neighbor (vn)
            vc, vn = 0;
            vc = rastrigin(current, string_length);
            vn = rastrigin(neighbor, string_length);
            //printf("Current best: %.10f\n", vc);
            //printf("Current new: %.10f\n", vn);
            // if neighboor is better than current, set current = neighbor
            if (vc > vn) // rastrigin is a minimization function, so check if vc is greater than vn
            {
                copyString(neighbor, current, string_length);
                vc = vn;
                repeat = 0;
                //printf("New best found: %.10f\n", vc);
            }
            // else if random(0-.99) < e^(vc-vn/temp), set vc = vn
            else if (rand_double(0.0, 0.999999999999999) < exp(-(vn-vc)/temp))
            {
                copyString(neighbor, current, string_length);
                vc = vn;
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
        if (isSameString(current, last, string_length))
        {
            repeat++;
            //printf("Core %d: Repeated value %.10f %d times\n", core, vc, repeat);
        }
        //printf("End of restart #%d\n", time);
        temp = (ALPHA * temp) / log(1 + time); // modified logarithmic decrease
        time++; // incriment time
        
    }
    for(i=0;i<string_length;i++) // copy local min string
    {
        min_strings[loop*max_loops + i] = current[i];
    }
    local_minima[loop] = vc; // save evaluation
    //gettimeofday(&stop_run_time, NULL);
    //timeData[loop] = (stop_run_time.tv_sec + stop_run_time.tv_usec/1000000) - (start_run_time.tv_sec + start_run_time.tv_usec/1000000);
    printf("Local min for core %d: %.5f | Found in %d loops\n", core, vc, time);
    loopNum[loop] = time; // save total iterations of SA for thread
    free(neighborhood);
    free(current);
    free(last);
    free(neighbor);
    free(best);
}

// add test suite

// once algorithm is done, print simple results to console, and detailed stats to a document

int main (int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    // initialize random
    time_t ts; 
    srand(time(&ts));

    int size, rank, i, r, max_loops, length;
    double min = 200000.0;
    int* local_loops;
    struct timeval start_total_time, stop_total_time;
    long int tt; // total time of SA run

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(argc != 3) 
    {
        if(rank == 0)
        {
            printf("%s\n", "Usage: \"Max # of Loops\" \"Size of number String\"");
            MPI_Finalize();
            exit(-1);
        }
    }
    else
    {
        max_loops = atoi(argv[1]);
        if(rank == 0)
        {
            //printf("Allocating %d loops to %d cores\n", max_loops, size);
            local_loops = (int*)malloc(sizeof(int*) * size);
            if(rank == 0) // setup simulated annealing run
            {
                r = max_loops % size; // calculate the remainder of loops after int division
                for(i = 0; i < size; i++)
                {
                    local_loops[i] =  max_loops / size; // evenly distribute the divisible number of loops to this core
                    if(r > 0) // if there are still loops in the remainder
                    {
                        local_loops[i]++; // add another one to this core
                        r--; // decriment remainder
                    }
                }
            }
            //printf("%s\n", "Broadcasting to cores");
        }
    }
    //printf("Rank: %d\n", rank);

    MPI_Bcast(local_loops, size, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the loop numbers for each core
    //printf("%d loops for rank %d\n", local_loops[rank], rank);

    char *s = (char*)malloc(29 * sizeof(char));
    sprintf(s, "./results/Core%d Results.txt", rank);
    FILE* output = fopen(s, "w+");
    //printf("Opening file %s\n", s);
    if(!output)
    {
        printf("%s\n", "Error opening file.");
        return -1;
    }

    length = atoi(argv[2]);

    double *ms = (double *)malloc(sizeof(double)*local_loops[rank]*length); // allocate memory for the minimum strings array: [number of loops][length of string of doubles]
    double *os = (double *)malloc(sizeof(double)*local_loops[rank]*length); // allocate memory for the origin strings array: [number of loops][length of string of doubles]
    double *lm = (double *)malloc(sizeof(double)*local_loops[rank]); // allocate memory for the results of the simulated annealing run: [number of loops]
    //long int *rt = (long int*)malloc(sizeof(long int*)*local_loops[rank]); // allocate memory for the time results of each run of the sa algorithm: [number of loops]
    int *tl = (int *)malloc(sizeof(int)*local_loops[rank]); // allocate memeory for the total iterations in a simulated annealing run: [number of loops]
    gettimeofday(&start_total_time, NULL);
    for(i = 0; i < local_loops[rank]; i++) // for the # of times on each core
    {
        printf("Beginning run #%d on core %d\n", i+1, rank);
        simulated_annealing(ms, os, lm, tl, i, local_loops[rank], length, rank, 100, 100, .8 );
    }
    gettimeofday(&stop_total_time, NULL);
    tt = (stop_total_time.tv_sec + stop_total_time.tv_usec/1000000) - (start_total_time.tv_sec + start_total_time.tv_usec/1000000);


    
    // print results to files
    for(i=0;i<local_loops[rank];i++)
    {
        if(min > lm[i])
        {
            min = lm[i];
        }
        printf("Core %d | Run %d | Local Min: %.10f\n", rank, i + 1, lm[i]);
    }
/*   
    printf("Loop %d:\nOrigin String: ", l + 1);
    for(i=0;i<STRING_LENGTH-1;i++)
    {
        printf("%.5f, ", os[l*local_loops[rank] + i]); 
    }
    printf("%.5f\n", os[l*local_loops[rank] + (STRING_LENGTH - 1)]);
    printf("Local Min String: ");
    for(i=0;i<STRING_LENGTH-1;i++)
    {
        printf("%.5f, ", ms[l*local_loops[rank] + i]);
    }
    printf("%.5f\n", ms[l*local_loops[rank] + (STRING_LENGTH - 1)]);
    fprintf(output, "Exection time per run:\n") 
    for(i = 0; i < local_loops[rank]; i++)
    {
        fprintf(output, "%.10f\n", rt[i]);
    }
    */
    fprintf(output, "Local Mins:\n");
    for(i=0;i<local_loops[rank];i++)
    {
        fprintf(output, "%.10f\n", lm[i]);
    }
    fprintf(output, "\nTotal Execution time: %ld seconds\n\n", tt);
    
    fprintf(output, "\n# of Threads: %d\nStarting Heat: %.2f\nAlpha: %.2f", size, 100.0, 0.8);
    fprintf(output, "\n# of Iterations: \n");
    for(i = 0; i < local_loops[rank]; i++)
    {
        fprintf(output, "%d\n", tl[i]);
    }

    free(ms);
    free(os);
    free(lm);
    //free(rt);
    free(tl);
    fclose(output); 
    MPI_Finalize();
    return 0;
}
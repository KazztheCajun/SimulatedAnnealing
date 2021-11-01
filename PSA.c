#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define MAX_LOOPS 1000
#define STRING_LENGTH 10
#define MAX_REPEAT_VALUE 100
#define _USE_MATH_DEFINES 
#define ALPHA 0.8



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
        //c[i] = 10.0; // generate debug values
        //c[i] = 1.0;
    }
}

// Generate a neighborhood for a given string
void get_neighbors(double n[][STRING_LENGTH], double *c)
{
    int i, j, index = 0;
    for (i = 0; i < STRING_LENGTH; i++) // for each column in neighborhood
    {
        for (j = 0; j < STRING_LENGTH; j++) // for each row in neighbor string
        {
            if (index == j) // if change index == row in neighbor string, generate new number
            {
                n[i][j] = rand_double(-5.12, 5.12);
            }
            else // otherwise copy current value
            {
                n[i][j] = c[j];
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


void simulated_annealing(double min_strings[][STRING_LENGTH], double origin_strings[][STRING_LENGTH], double local_minima[MAX_LOOPS])
{
    
    // initialize time, temp, repeat counter, and loop counter
    int time, temp, repeat, loop, i , index = 0;
    // arrays to track found local maxima, current string, last used string, & previously used strings
    double vc, vn, best[STRING_LENGTH], current[STRING_LENGTH], last[STRING_LENGTH], neighbor[STRING_LENGTH], neighborhood[STRING_LENGTH][STRING_LENGTH];
    // select a starting point
    new_string(current);
    // first loop that tracks total iterations of the algorithm
    while(loop < MAX_LOOPS)
    {
        //printf("\nRestart # %d ", loop + 1);
        // (re)initialize time, current point, temp, and repeated maxima values
        time = 1;
        temp = 100;
        repeat = 0;
        
        for(i=0;i<STRING_LENGTH;i++) // save origin string for each loop
        {
            origin_strings[loop][i] = current[i];
        }
        // second loop that loops until the current value is reapeated enough
        while(repeat < MAX_REPEAT_VALUE)
        {
            //printf("Loop %d\n", time + 1);
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
                neighbor[i] = neighborhood[index][i];
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
            // otherwise do nothing and try a new neighboor until max repeated value reached
            if (isSameString(current, last))
            {
                repeat++;
                //printf("Repeated value %.10f %d times\n", vc, repeat);
            }
            temp = (ALPHA * temp) / log(1 + time); // logarithmic decrease
            time++; // incriment time
        }
        for(i=0;i<STRING_LENGTH;i++) // copy local min string
        {
            min_strings[loop][i] = current[i];
        }
        local_minima[loop] = rastrigin(current); // save evaluation
        loop++;
        //printf("End of restart #%d\n", loop -1);
    }
}





// once algorithm is done, print simple results to console, and detailed stats to a document

void main()
{
    // initialize random
    srand(time(NULL));
    int i, j, l;
    double min = 2000.0;
    double ms[MAX_LOOPS][STRING_LENGTH];
    double os[MAX_LOOPS][STRING_LENGTH];
    double lm[MAX_LOOPS];
    simulated_annealing(ms, os, lm);
    
    for(l = 0; l < MAX_LOOPS; l++)
    {
        double temp[STRING_LENGTH];
        printf("Original String:\n");
        for(i=0;i<STRING_LENGTH-1;i++)
        {
            temp[i] = os[l][i];
            printf("%.10f, ", os[l][i]); 
        }
        printf("%.10f (%.10f)\n", os[l][STRING_LENGTH - 1], rastrigin(temp));
        printf("Local minima of %f found by:\n", lm[l]);
        for(i=0;i<STRING_LENGTH-1;i++)
        {
            printf("%.10f, ", ms[l][i]);
        }
        printf("%.10f\n\n", ms[l][STRING_LENGTH - 1]);
    }

    for(i=0;i<MAX_LOOPS;i++)
    {
        if(min > lm[i])
        {
            min = lm[i];
        }
    }
    printf("The global min is: %.5f", min);
    
}
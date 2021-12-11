# Simulated Annealing Agorithms

### SerialSA.c 
- An implimentation of my simulated annealing algorithm designed to execute total iterations in serial (loop 1 then loop 2 then loop 3)

### ParallelSA.c 
- Initial Open MP implimentation of my SA algorithm, designed to execute the total iterations in parallel in a Shared Memory environment.

### BenchmarkPSA_OMP.c 
- Modified implimentation of ParallelSA, using Open MP programming practices and is designed to accept input for the hyperparameters for tuning and experimentation purposes.

### BenchmarkSSA.c
- Modified implimentation of SerialSA, designed to execute serially and accept input for the hyperparameters for tuning and experimentation purposes.

### MPI_PSA.c 
- Initial MPI implimentation of my SA algorithm, designed to execute the total iterations in parallel in a Distrubuted Memory environment.

# How to run the code

I have created a convenient benchmarking script to test the BenchmarkPSA_OMP and BenchmarkSSA

Just extract the object files for Benchmark.c, BenchmarkPSA_OMP.c, and BenchmarkSSA.c and then compile them into an executable

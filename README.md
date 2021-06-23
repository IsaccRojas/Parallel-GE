# Parallel GE

Simple implementation of a parallel Gauss-Jordan elimination algorithm with block-cyclic distribution in C with OpenMPI, for demonstration purposes. Includes a helper program to generate valid random systems of equations as matrices with a column vector.

## Usage

### parallel_ge.c

This program is to be compiled and run in an MPI environment, as it will query the MPI system for the number of processors. It has the following parameter syntax:

>./parallel_ge [inputfile] [n] [r] [bs] [genUX]

where

- `inputfile` is the input file containing an nxn matrix with an additional vector column,
- `n` is the matrix size,
- `r` is the block-cyclic parameter,
- `bs` is the back substitution method to use (0 for sequential, 1 for parallel),
- `genUX` is whether to generate a file containing the final UX matrix (reduced row-echelon form matrix U and solution vector X),

e.g.

>./parallel_ge matrix_1024.txt 1024 32 0 1

which will take a 1024x1024 matrix input file with an additional vector column, use 32 blocks per cycle of distribution, use the sequential back substition method, and will generate the UX file.

### random_gen.c

This program is to be linked with the C math library. It has the following parameter syntax:

>./random_gen [n] [outputfile]

where

- `n` is the matrix size,
- `outputfile` is the output file name,

e.g.

>./random_gen 1024 matrix_1024.txt

which will generate a 1024x1024 random system of equation as a matrix with an additional column vector, with the filename matrix_1024.txt.

## License

This software is licensed under the MIT License.

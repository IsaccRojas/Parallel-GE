# Parallel GE

Simple SPMD implementation of a parallel Gauss-Jordan elimination algorithm with block-cyclic distribution in C with OpenMPI, for demonstration purposes. Includes a helper program to generate valid random systems of equations as matrices with a column vector.

## Algorithm

### 1. Data Distribution

For **N** rows, **P** processors, and **R** cycles, data is loaded into a root processor and blocks of **N** / (**R** * **P**) rows are distributed to each processor, cycling back to the root processor every **P** send operations, and completing after **R** * **P** total send operations:

![Data distribution figure 1.](https://i.imgur.com/g5qbnXh.png)
![Data distribution figure 2.](https://i.imgur.com/0t4E8G8.png)

### 2. Elimination

Each row is iterated on in the complete matrix, determining which processor owns the pivot row and that processor distributing it to all other processors, to eliminate the leading coefficient of all rows "below" the pivot row, eventually yielding the upper triangular matrix U:

![Elimination figure 1.](https://i.imgur.com/KhRLdZ6.png)

### 3. Back Substitution

The following equation is then solved for each element of the solution vector **x**:

![Back substitution figure 1.](https://i.imgur.com/YmNYgnZ.gif)

#### Sequential Back Substitution

Each row is iterated on in the complex matrix from the final row, with a processor sending the partial solution vector **x** to the next processor in the "sequence":

![Back substitution figure 2.](https://i.imgur.com/miuILeK.png)

#### Parallel Back Substitution

Each row is iterated on in the complex matrix from the final row, with a processor sending the partial solution vector **x** *within a cycle* to the next processor in the "sequence"; the last processor to be computing at the end of a cycle distributes the **N** / **R** rows of **x** computed in that cycle to all other processors, to allow forward computation of the equation:

![Back substitution figure 3.](https://i.imgur.com/lyNt0RM.png)

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

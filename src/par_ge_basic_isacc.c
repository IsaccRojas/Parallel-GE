#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <time.h>

//compute log2 of an int (truncated)
int log2_int(int n) {
	int i;
	for (i = 31; i >= 0; i--)
		if (n >> i == 1)
			return i;
	return 0;
}

//compute b^exp for ints (truncated)
int pow_int(int b, int exp) {
	int b_start = b;
	int i;
	for (i = 1; i < exp; i++)
		b = b * b_start;
	return b;
}

int max(int a, int b) {
	if (a > b)
		return a;
	else
		return b;
}

int min(int a, int b) {
	if (a < b)
		return a;
	else
		return b;
}

//read string at offset up to end of word and write into buffer; fd must be an open file descriptor
//returned pointer must be freed later
char* read_word(int fd, off_t offset, off_t *end) {
	lseek(fd, offset, SEEK_SET);
	
	char *c = malloc(sizeof(char));
	int pos = offset;
	char *buf;
	while (1) {
		read(fd, c, 1);
		if (*c == ' ' || *c == '\n' || *c == '\0') {
			buf = malloc(sizeof(char) * (pos - offset + 1));
			
			//set file offset back to given offset and read bytes
			lseek(fd, offset, SEEK_SET);
			read(fd, buf, pos - offset + 1);
			buf[pos - offset] = '\0';
			*end = pos;
			
			free(c);	
			return buf;
		}
		pos++;
	}
	return NULL;
}

int main(int argc, char *argv[]) {
	//parameters/variables
	int n, id, bs, r, config, genUX;
	int nodes, myid, numprocs, i, j, k, offset, rows, cols, fd, starting_i, data_size, block_size, block_rows, cycle, cycle_rows, pivotproc, local_i, global_i, global_k, partial_sum;
	double PI25DT = 3.141592653589793238462643;
	double val, pivotval, sum, runtime_dd, runtime_ge, runtime_bs, x_checksum_even, x_checksum_odd, u_checksum_even, u_checksum_odd;
	double *mat;
	double *global_mat;
	double *pivotrow;
	double *data;
	double *x;
	double *partial_sums;
	off_t word_start, word_end;
	char *c;
	char *word;

	//initialize OpenMPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if (myid == 0) {
		if (argc < 9) {
			printf("error: invalid number of arguments\n");
			return -1;
		}

		//get n and r
		n = atoi(argv[1]);
		id = atoi(argv[2]);
		fd = open(argv[3], O_RDONLY);
		bs = atoi(argv[4]);
		r = atoi(argv[5]);
		if (strcmp(argv[6], "SM") == 0)
			nodes = max(1, numprocs / 16);
		else
			nodes = min(8, numprocs);
		genUX = atoi(argv[7]);
		
		//send n and r to processors
		for (k = 1; k < numprocs; k++) {
			MPI_Send(&n, 1, MPI_INT, k, 0, MPI_COMM_WORLD);	
			MPI_Send(&r, 1, MPI_INT, k, 0, MPI_COMM_WORLD);	
			MPI_Send(&bs, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
		}

		global_mat = malloc(sizeof(double) * (n * (n + 1)));
		c = malloc(sizeof(char));
		offset = 0; i = 0; j = 0; k = 0; cols = 0;
		while (1) {
			//reset and get value
			lseek(fd, 0, SEEK_SET);
			word = read_word(fd, offset, &word_end);
			val = atof(word);
			free(word);
			global_mat[k] = val;
				
			j++;
			if (j == n + 1) {
				i++;
				j = 0;
			}
			if (i == n)
				break;
			offset = word_end + 1;
			
			//get character at end of word
			lseek(fd, word_end, SEEK_SET);
			read(fd, c, 1);
			if (*c == '\0')
				break;

			k++;	
		}
		close(fd);
	
		rows = n;
		cols = n + 1;
	
		if (i != n) {
			printf("error: did not get n rows with n + 1 columns (got %d rows)\n", rows);
			return -1;
		}

		//distribute rows to all processors
		runtime_dd = MPI_Wtime();

		cycle_rows = n / r;
		block_rows = cycle_rows / numprocs;
		block_size = (n + 1) * block_rows;
		data_size = r * block_size;
		mat = malloc(sizeof(double) * data_size);
		for (k = 0; k < numprocs * r; k++) {
			if (k % numprocs != 0)
				MPI_Send(global_mat + (k * block_size), block_size, MPI_DOUBLE, k % numprocs, k / numprocs, MPI_COMM_WORLD);
			else
				memcpy(mat + ((k / numprocs) * block_size), global_mat + (k * block_size), block_size * sizeof(double));
		}
	} else {
		//receive n, r, and bs from "root" processor
		MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&r, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&bs, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
		//receive r rows from "root" processor
		cycle_rows = n / r;
		block_rows = cycle_rows / numprocs;
		block_size = (n + 1) * block_rows;
		data_size = (n * (n + 1)) / numprocs;
		mat = malloc(sizeof(double) * data_size);
		for (k = 0; k < r; k++)
			MPI_Recv(mat + (k * block_size), block_size, MPI_DOUBLE, 0, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	//perform elimination for each column
	runtime_dd = MPI_Wtime() - runtime_dd;
	runtime_ge = MPI_Wtime();
		
	//starting iteration (disposition) within a cycle for processor myid
	starting_i = (int)(cycle_rows * ((float)myid / (float)numprocs));
	pivotrow = malloc(sizeof(double) * (n + 1));
	for (i = 0; i < n; i++) {
		cycle = i / cycle_rows;

		//processor id with pivot row
		pivotproc = (i / block_rows) % numprocs;		

		//relative row position; if my id == pivotproc, then must be in [0, n/rp)
		//shift by a cycle disposition and number of blocks in a cycle minus one
		local_i = i - starting_i - (cycle * block_rows * (numprocs - 1));
				
		//processor with pivot row distributes to all other procs
		if (myid == pivotproc) {
			for (k = 0; k < numprocs; k++)
				if (k != myid && (cycle + 1 < r || k > myid))
					MPI_Send(mat + ((n + 1) * local_i), n + 1, MPI_DOUBLE, k, i + numprocs, MPI_COMM_WORLD);

			memcpy(pivotrow, mat + ((n + 1) * local_i), sizeof(double) * (n + 1));
		} else if (cycle + 1 < r || myid > pivotproc)
			MPI_Recv(pivotrow, n + 1, MPI_DOUBLE, pivotproc, i + numprocs, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		//apply pivot row to local rows "below" pivot row
		for (k = 0; k < block_rows * r; k++) {
			//global row position
			//shift by a cycle disposition and number of blocks in cycle of k minus one
			global_k = k + starting_i + ((k / block_rows) * block_rows * (numprocs - 1)); 
			
			if (global_k > i) {
				pivotval = mat[i + (k * (n + 1))];
				
				//apply to each column of row k
				for (j = i; j < n + 1; j++) {
					//mat[i][j] -= pivotrow[j] * (mat[i][i] / pivotrow[i])
					//a_ij = a_ij - R_j * (a_ii / R_i)				
					mat[j + (k * (n + 1))] -= pivotrow[j] * (pivotval / pivotrow[i]);
				}
			}
		}
	}

	//perform back substitution
	runtime_ge = MPI_Wtime() - runtime_ge;
	runtime_bs = MPI_Wtime();

	x = malloc(sizeof(double) * n);
	if (bs == 0) {	
		for (k = 0; k < r; k++) {
			//receive x values (if not at start of "chain")
			if (!((k == 0) && (myid == numprocs - 1)))
				MPI_Recv(
					x + (block_rows * (myid + 1)) + ((r - k - 1) * cycle_rows),
					(block_rows * (numprocs - myid - 1)) + (k * cycle_rows),
					MPI_DOUBLE,
					myid == numprocs - 1 ? 0 : myid + 1,
					k + n + numprocs,
					MPI_COMM_WORLD,
					MPI_STATUS_IGNORE
				);
			
			//get x for all local rows
			for (i = ((r - k) * block_rows) - 1; i >= (r - k - 1) * block_rows; i--) {
				sum = 0.0;
				global_i = i + starting_i + ((i / block_rows) * block_rows * (numprocs - 1));
				
				//substitute previous x's into current row
				for (j = global_i + 1; j < n; j++)
					sum += mat[j + (i * (n + 1))] * x[j];
	
				//solve for x
				x[global_i] = (1 / mat[global_i + (i * (n + 1))]) * (mat[n + (i * (n + 1))] - sum);
			}
			
			//send x values (if not at end of "chain"; non-blocking)
			if (!((k == r - 1) && (myid == 0)))
				MPI_Send(
					x + (block_rows * myid) + ((r - k - 1) * cycle_rows),
					(block_rows * (numprocs - myid)) + (k * cycle_rows),
					MPI_DOUBLE, 
					myid == 0 ? numprocs - 1 : myid - 1, 
					n + numprocs + (myid == 0 ? k + 1 : k), 
					MPI_COMM_WORLD
				);		
		}
	} else {
		partial_sums = malloc(sizeof(double) * n);
		for (k = 0; k < n; k++)
			partial_sums[k] = 0.0;

		for (k = 0; k < r; k++) {
			//receive n/r "broadcasted" x values (if not at start of "chain" and not proc 0)
			if ((k != 0) && (myid != 0)) {	
				MPI_Recv(
					x + (cycle_rows * (r - k)),
					cycle_rows,
					MPI_DOUBLE,
					0,
					(2 * k) + n + numprocs,
					MPI_COMM_WORLD,
					MPI_STATUS_IGNORE
				);
				
				//use received x's to get partition of sum for computation of local x's
				for (i = ((r - k) * block_rows) - 1; i >= 0; i--) {
					global_i = i + starting_i + ((i / block_rows) * block_rows * (numprocs - 1));
					
					//substitute previous "broadcasted" received x's into current row
					for (j = n - (k * cycle_rows); j < n - ((k - 1) * cycle_rows); j++)
						partial_sums[global_i] += mat[j + (i * (n + 1))] * x[j];
				}		
			}
			
			//receive 1 to n/r "sequential" x values (if not proc numprocs - 1)
			if (!(myid == numprocs - 1))
				MPI_Recv(
					x + (block_rows * (myid + 1)) + ((r - k - 1) * cycle_rows),
					(block_rows * (numprocs - myid - 1)),
					MPI_DOUBLE,
					myid + 1,
					((2 * k) + 1) + n + numprocs,
					MPI_COMM_WORLD,
					MPI_STATUS_IGNORE
				);

			//get x for all local rows
			for (i = ((r - k) * block_rows) - 1; i >= (r - k - 1) * block_rows; i--) {
				global_i = i + starting_i + ((i / block_rows) * block_rows * (numprocs - 1));
				sum = partial_sums[global_i];
				
				//substitute previous "sequentially" received x's into current row
				for (j = global_i + 1; j < n - (k * cycle_rows); j++)
					sum += mat[j + (i * (n + 1))] * x[j];
				
				//solve for x
				x[global_i] = (1 / mat[global_i + (i * (n + 1))]) * (mat[n + (i * (n + 1))] - sum);
			}
			
			//"broadcast" n/rp x values (if at end of cycle, else send "sequential" x values)
			if (myid == 0 && !(k + 1 >= r)) {
				for (i = numprocs - 1; i >= 0; i--)
					if (i != myid)
						MPI_Send(
							x + ((r - k - 1) * cycle_rows),
							cycle_rows,
							MPI_DOUBLE,
							i,
							(2 * (k + 1)) + n + numprocs,
							MPI_COMM_WORLD
						);

				//use local x's to get partition of sum for computation of later local x's
				for (i = ((r - (k + 1)) * block_rows) - 1; i >= 0; i--) {
					global_i = i + starting_i + ((i / block_rows) * block_rows * (numprocs - 1));
					
					//substitute previous "broadcasted" received x's into current row
					for (j = n - ((k + 1) * cycle_rows); j < n - (k * cycle_rows); j++)
						partial_sums[global_i] += mat[j + (i * (n + 1))] * x[j];
				}		
	

			} else
				MPI_Send(
					x + (block_rows * myid) + ((r - k - 1) * cycle_rows),
					(block_rows * (numprocs - myid)),
					MPI_DOUBLE, 
					myid - 1, 
					((2 * k) + 1) + n + numprocs, 
					MPI_COMM_WORLD
				);
		}
	}

	//get matrix U and get checksums
	runtime_bs = MPI_Wtime() - runtime_bs;
	
	//send each block to 0
	for (k = 0; k < r * numprocs; k++) {
		if (k % numprocs == myid) {
			if (myid == 0) {
				memcpy(
					global_mat + (block_size * k), 
					mat + (block_size * (k / numprocs)), 
					block_size * sizeof(double)
				);
			} else {
				MPI_Send(
					mat + (block_size * (k / numprocs)),
					block_size,
					MPI_DOUBLE,
					0,
					n + numprocs + r + k,
					MPI_COMM_WORLD
				);
			}
		} else if (myid == 0)
			MPI_Recv(
				global_mat + (block_size * k),
				block_size,
				MPI_DOUBLE,
				k % numprocs,
				n + numprocs + r + k,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE
			);
	}
	
	if (myid == 0) {
		//get X checksums
		x_checksum_even = 0.0;
		x_checksum_odd = 0.0;
		for (k = 0; k < n; k++) {
			if (k % 2 == 0)
				x_checksum_even += x[k];
			else
				x_checksum_odd += x[k];
		}

		//get U checksums
		u_checksum_even = 0.0;
		u_checksum_odd = 0.0;
		for (k = 0; k < n * (n + 1); k++) {
			if (k % 2 == 0)
				u_checksum_even += global_mat[k];
			else
				u_checksum_odd += global_mat[k];
		}	
	}
	
	if (myid == 0) {
		//get and open stats filename
		char filename[256];
		sprintf(filename, "op_P%d_%s_bs%d_r%d_n%d_%d_stats_%s.txt", numprocs, argv[6], bs, r, n, id, argv[8]);
		fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0666);

		//print metrics
		dprintf(fd, "P: %d\n", numprocs);
		dprintf(fd, "n: %d\n", n);
		dprintf(fd, "r: %d\n", r);
		dprintf(fd, "bs: %d\n", bs);
		dprintf(fd, "combination: nodes = %d, ppn = %d\n", nodes, numprocs / nodes);
		dprintf(fd, "total runtime: \t%.16fs\n", runtime_dd + runtime_ge + runtime_bs);
		dprintf(fd, "ge runtime: \t%.16fs\n", runtime_ge);
		dprintf(fd, "ge+bs runtime: \t%.16fs\n", runtime_ge + runtime_bs);
		dprintf(fd, "X checksum (even): %.16f\n", x_checksum_even);
		dprintf(fd, "X checksum (odd): %.16f\n", x_checksum_odd);
		dprintf(fd, "U checksum (even): %.16f\n", u_checksum_even);
		dprintf(fd, "U checksum (odd): %.16f\n", u_checksum_odd);
		close(fd);

		if (genUX) {
			//get and open UX filename
			sprintf(filename, "op_P%d_%s_bs%d_r%d_n%d_%d_UX_%s.txt", numprocs, argv[6], bs, r, n, id, argv[8]);
			fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0666);
			
			//print metrics
			dprintf(fd, "P: %d\n", numprocs);
			dprintf(fd, "n: %d\n", n);
			dprintf(fd, "r: %d\n", r);
			dprintf(fd, "bs: %d\n", bs);
			dprintf(fd, "combination: nodes = %d, ppn = %d\n", nodes, numprocs / nodes);
			for (i = 0; i < n; i++) {
				dprintf(fd, "%.16f | ", x[i]);
				for (j = 0; j < n; j++)
					dprintf(fd, "%.16f ", global_mat[j + (i * (n + 1))]);
				dprintf(fd, "\n");
			}
			dprintf(fd, "X checksum (even): %.16f\n", x_checksum_even);
			dprintf(fd, "X checksum (odd): %.16f\n", x_checksum_odd);
			dprintf(fd, "U checksum (even): %.16f\n", u_checksum_even);
			dprintf(fd, "U checksum (odd): %.16f\n", u_checksum_odd);		
		}
	}	
	
	MPI_Finalize();	
}

const char* dgemm_desc = "Recursive dgemm.";

#define likely(x)	__builtin_expect((x),1)
#define unlikely(x)	__builtin_expect((x),0)

#if !defined(BLOCK_SIZE2)
#define BLOCK_SIZE2 160
#endif

#if !defined(BLOCK_SIZE1)
#define BLOCK_SIZE1 40
#endif

#define min(a,b) (((a)<(b))?(a):(b))

#include <stdlib.h>

/* This auxiliary subroutine performs a smaller dgemm operation
 *	C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, double* restrict A, double* restrict B, double* restrict C) {
	//__builtin_prefetch(A + k*lda, 0, 1);
	/* For each column j of B */ 
	for (int j = 0; j < N; ++j) {
		/* For each row i of A */
		for (int i = 0; i < M; ++i) {
			/* Compute C(i,j) */
			double cij = C[i+j*lda];
			for (int k = 0; k < K; ++k)
				cij += A[i+k*lda] * B[k+j*lda];
				C[i+j*lda] = cij;
		}
	}
}

//transposed version
static void do_block_T (int lda, int M, int N, int K, double* restrict A_T, double* restrict B, double* restrict C) {
	//__builtin_prefetch(A + k*lda, 0, 1);
	/* For each column j of B */ 
	for (int j = 0; j < N; ++j) {
		/* For each row i of A */
		for (int i = 0; i < M; ++i) {
			/* Compute C(i,j) */
			double cij = C[i+j*lda];
			for (int k = 0; k < K; ++k)
				cij += A_T[k+i*lda] * B[k+j*lda];
				C[i+j*lda] = cij;
		}
	}
}

/* This auxiliary subroutine performs a smaller dgemm operation
 *	C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N.
 * This function assumes that the submatrix is size BLOCK_SIZE
 */
static void do_exact_block (int lda, double* restrict A, double* restrict B, double* restrict C) {
	//	__builtin_prefetch(A + k*lda, 0, 1);
	/* For each column j of B */ 
	for (int j = 0; j < BLOCK_SIZE1; ++j) {
		/* For each row i of A */
		for (int i = 0; i < BLOCK_SIZE1; ++i) {
			/* Compute C(i,j) */
			double cij = C[i+j*lda];
			for (int k = 0; k < BLOCK_SIZE1; ++k) {
				cij += A[i+k*lda] * B[k+j*lda];
			}
			C[i+j*lda] = cij;
		}
	}
}

//transposed version
static void do_exact_block_T (int lda, double* restrict A_T, double* restrict B, double* restrict C) {
	//	__builtin_prefetch(A + k*lda, 0, 1);
	/* For each column j of B */ 
	for (int j = 0; j < BLOCK_SIZE1; ++j) {
		/* For each row i of A */
		for (int i = 0; i < BLOCK_SIZE1; ++i) {
			/* Compute C(i,j) */
			double cij = C[i+j*lda];
			for (int k = 0; k < BLOCK_SIZE1; ++k) {
				cij += A_T[k+i*lda] * B[k+j*lda];
			}
			C[i+j*lda] = cij;
		}
	}
}

static void L1_dgemm (int lda, int I, int J, int K, double* restrict A_T, double* restrict B, double* restrict C) {
	/* For each block-column of B */
	for (int j = 0; j < J; j += BLOCK_SIZE1)
		/* For each block-row of A */ 
		for (int i = 0; i < I; i += BLOCK_SIZE1)
			/* Accumulate block dgemms into block of C */
			for (int k = 0; k < K; k += BLOCK_SIZE1) {
				if ( (i<=(I-BLOCK_SIZE1)) && (j<=(J-BLOCK_SIZE1)) && (k<=(K-BLOCK_SIZE1)) ){
					//do_exact_block(lda, A + i + k*lda, B + k + j*lda, C + i + j*lda);
					do_exact_block_T(lda, A_T + k + i*lda, B + k + j*lda, C + i + j*lda);
				} else {
					/* Correct block dimensions if block "goes off edge of" the matrix */
					int I_ = min (BLOCK_SIZE1, I-i);
					int J_ = min (BLOCK_SIZE1, J-j);
					int K_ = min (BLOCK_SIZE1, K-k);
					/* Perform individual block dgemm */
					do_block_T(lda, I_, J_, K_, A_T + k + i*lda, B + k + j*lda, C + i + j*lda);
				}
			}
}

inline static void L2_dgemm (int lda, double* restrict A_T, double* restrict B, double* restrict C) {
	/* For each block-column of B */
	for (int j = 0; j < lda; j += BLOCK_SIZE2)
		/* For each block-row of A */ 
		for (int i = 0; i < lda; i += BLOCK_SIZE2)
			/* Accumulate block dgemms into block of C */
			for (int k = 0; k < lda; k += BLOCK_SIZE2) {
				/* Correct block dimensions if block "goes off edge of" the matrix */
				int I_ = min (BLOCK_SIZE2, lda-i);
				int J_ = min (BLOCK_SIZE2, lda-j);
				int K_ = min (BLOCK_SIZE2, lda-k);
				/* Perform individual block dgemm */
				//do_block_T(lda, M, N, K, A_T + k + i*lda, B + k + j*lda, C + i + j*lda);

				L1_dgemm(lda, I_, J_, K_, A_T + k + i*lda, B + k + j*lda, C + i + j*lda);
			}
}

// Resursive matrix multiplication
// lda: size of original square matrix
// I, J, K: size of the current blocks
// is, is, is: starting indixes
static void rmm(int lda, int I, int J, int K, int is, int js, int ks, double* restrict A_T, double* restrict B, double* restrict C) {
	int i, j, k;
	int Ih1, Jh1, Kh1; //first half sizes
	int Ih2, Jh2, Kh2; //second half sizes

	//can we solve directly?
	if ((I==1)||(J==1)||(K==1)) {
		//TODO!
		return;
	}

	//calculate half sizes
	Ih2=I>>1; Jh2=J>>1; Kh2=K>>1;
	Ih1=I-Ih2; Jh1=J-Jh2; Kh1=K-Kh2;
	//C11
	rmm(lda, Ih1, Jh1, Kh1, is, js, ks, A_T, B, C);
	rmm(lda, Ih1, Jh1, Kh2, is, js, ks + Kh1, A_T+Kh1, B+Kh1, C);
	//C22
	


	//if 

}

/* This routine performs a dgemm operation
 *	C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format. 
 * On exit, A and B maintain their input values. */	
void square_dgemm (int lda, double* restrict A, double* restrict B, double* restrict C) {
	double* restrict A_T;

	//transpose A
	//TODO: blocked version?
	A_T = (double*) malloc(sizeof(double)*lda*lda);
	for (int i=0; i<lda; i++)
		for (int j=0; j<lda; j++)
			A_T[i + lda*j] = A[i*lda + j];

	rmm(lda, lda, lda, lda, 0, 0, 0, A_T, B, C);

	free(A_T);
}

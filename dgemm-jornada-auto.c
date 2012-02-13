const char* dgemm_desc = "jornada1";

#define likely(x)	__builtin_expect((x),1)
#define unlikely(x)	__builtin_expect((x),0)

#define ALIGNED_BLOCKS
#define TRANSPOSE

#if !defined(BLOCK_SIZE2)
#define BLOCK_SIZE2 160
#endif

#if !defined(BLOCK_SIZE1)
#define BLOCK_SIZE1 40
#endif

#define min(a,b) (((a)<(b))?(a):(b))

#include <stdlib.h>

//array of pointers to functions
//typedef static void (func)(int, double* restrict, double* restrict, double* restrict);
//typedef func* func_ptr;
//func_ptr [BLOCK_SIZE1];
//static void (*exact_blocks[BLOCK_SIZE1])(int, double* restrict, double* restrict, double* restrict);

#include "auto_blocks_inc.c"

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

	/* For each column j of B */
	int j_lda = 0;
	for (int j = 0; j < N; ++j, j_lda+=lda) {
		//__builtin_prefetch(B + j_lda, 0, 2);
		//__builtin_prefetch(C + j_lda, 1, 2);
		/* For each row i of A */
		int i_lda = 0;
		int ij_lda = j_lda;
		for (int i = 0; i < M; ++i, i_lda+=lda, ij_lda++) {
			//__builtin_prefetch(A_T + i_lda, 0, 2);
			/* Compute C(i,j) */
			double cij = C[ij_lda];
			int ki_lda = i_lda;
			int kj_lda = j_lda;
			for (int k = 0; k < K; ++k) {
				cij += A_T[ki_lda++] * B[kj_lda++];
			}
			C[ij_lda] = cij;
		}
	}
}

/* This auxiliary subroutine performs a smaller dgemm operation
 *	C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N.
 * This function assumes that the submatrix is size BLOCK_SIZE
 */
static void do_exact_block (int lda, double* restrict A, double* restrict B, double* restrict C) {
	/* For each column j of B */ 
	for (int j = 0; j < BLOCK_SIZE1; ++j) {
		int j_lda = j * lda;
		/* For each row i of A */
		for (int i = 0; i < BLOCK_SIZE1; ++i) {
			int ij_lda = i + j_lda;
			/* Compute C(i,j) */
			double cij = C[ij_lda];
			for (int k = 0; k < BLOCK_SIZE1; ++k) {
				cij += A[i+k*lda] * B[k+j_lda];
			}
			C[ij_lda] = cij;
		}
	}
}

//transposed version
static void do_exact_block_T (int lda, double* restrict A_T, double* restrict B, double* restrict C) {
	/* For each column j of B */ 
	int j_lda = 0;
	for (int j = 0; j < BLOCK_SIZE1; ++j, j_lda+=lda) {
		/* For each row i of A */
		int i_lda = 0;
		int ij_lda = j_lda;
		for (int i = 0; i < BLOCK_SIZE1; ++i, i_lda+=lda, ij_lda++) {
			/* Compute C(i,j) */
			double cij = C[ij_lda];
			int ki_lda = i_lda;
			int kj_lda = j_lda;
			for (int k = 0; k < BLOCK_SIZE1; ++k) {
				cij += A_T[ki_lda++] * B[kj_lda++];
			}
			C[ij_lda] = cij;
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
				//__builtin_prefetch(A_T + k + i*lda, 0, 2);
				//__builtin_prefetch(B   + k + j*lda, 0, 2);
				//__builtin_prefetch(C   + i + j*lda, 1, 2);
				//__builtin_prefetch(C   + i + j*lda, 0, 2);
				if ( (i<=(I-BLOCK_SIZE1)) && (j<=(J-BLOCK_SIZE1)) && (k<=(K-BLOCK_SIZE1))){
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

//size of submatrix == BLOCK_SIZE2
static void L1_dgemm_exact (int lda, double* restrict A_T, double* restrict B, double* restrict C) {
	/* For each block-column of B */
	for (int j = 0; j < BLOCK_SIZE2; j += BLOCK_SIZE1)
		/* For each block-row of A */ 
		for (int i = 0; i < BLOCK_SIZE2; i += BLOCK_SIZE1)
			/* Accumulate block dgemms into block of C */
			for (int k = 0; k < BLOCK_SIZE2; k += BLOCK_SIZE1) {
#ifdef ALIGNED_BLOCKS
				do_exact_block_T(lda, A_T + k + i*lda, B + k + j*lda, C + i + j*lda);
#else
				if ( (i<=(BLOCK_SIZE2-BLOCK_SIZE1)) && (j<=(BLOCK_SIZE2-BLOCK_SIZE1)) && (k<=(BLOCK_SIZE2-BLOCK_SIZE1)) ){
					do_exact_block_T(lda, A_T + k + i*lda, B + k + j*lda, C + i + j*lda);
				} else {
					/* Correct block dimensions if block "goes off edge of" the matrix */
					int I_ = min (BLOCK_SIZE1, BLOCK_SIZE2-i);
					int J_ = min (BLOCK_SIZE1, BLOCK_SIZE2-j);
					int K_ = min (BLOCK_SIZE1, BLOCK_SIZE2-k);
					/* Perform individual block dgemm */
					do_block_T(lda, I_, J_, K_, A_T + k + i*lda, B + k + j*lda, C + i + j*lda);
				}
#endif
			}
}

inline static void L2_dgemm (int lda, double* restrict A_T, double* restrict B, double* restrict C) {
	/* For each block-column of B */
	for (int j = 0; j < lda; j += BLOCK_SIZE2)
		/* For each block-row of A */ 
		for (int i = 0; i < lda; i += BLOCK_SIZE2)
			/* Accumulate block dgemms into block of C */
			for (int k = 0; k < lda; k += BLOCK_SIZE2) {
				if ( (i<=(lda-BLOCK_SIZE2)) && (j<=(lda-BLOCK_SIZE2)) && (k<=(lda-BLOCK_SIZE2)) ){
					//do_exact_block(lda, A + i + k*lda, B + k + j*lda, C + i + j*lda);
					L1_dgemm_exact(lda, A_T + k + i*lda, B + k + j*lda, C + i + j*lda);
				} else {
					/* Correct block dimensions if block "goes off edge of" the matrix */
					int I_ = min (BLOCK_SIZE2, lda-i);
					int J_ = min (BLOCK_SIZE2, lda-j);
					int K_ = min (BLOCK_SIZE2, lda-k);
					/* Perform individual block dgemm */
					//do_block_T(lda, M, N, K, A_T + k + i*lda, B + k + j*lda, C + i + j*lda);

					L1_dgemm(lda, I_, J_, K_, A_T + k + i*lda, B + k + j*lda, C + i + j*lda);
				}
			}
}

/* This routine performs a dgemm operation
 *	C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format. 
 * On exit, A and B maintain their input values. */	
void square_dgemm (int lda, double* restrict A, double* restrict B, double* restrict C) {
	double* restrict A_T;

	//transpose A
	//TODO: blocked version?
#ifdef TRANSPOSE
	A_T = (double*) malloc(sizeof(double)*lda*lda);
	for (int i=0; i<lda; i++)
		for (int j=0; j<lda; j++)
			A_T[i + lda*j] = A[i*lda + j];
#endif

	if (lda<=BLOCK_SIZE1) {
		//do_block_T(lda, lda, lda, lda, A_T, B, C);
		exact_blocks[lda](lda, A_T, B, C);
	} else if (lda<=BLOCK_SIZE2){
		L1_dgemm(lda, lda, lda, lda, A_T, B, C);
	} else {
		L2_dgemm(lda, A_T, B, C);
	}

	free(A_T);
}

const char* dgemm_desc = "jornada1";

#define likely(x)	__builtin_expect((x),1)
#define unlikely(x)	__builtin_expect((x),0)

#define ALIGNED_BLOCKS
#define TRANSPOSE

#if !defined(BLOCK_SIZE2)
#define BLOCK_SIZE2 144
#endif

#if !defined(BLOCK_SIZE1)
#define BLOCK_SIZE1 24
#endif

#if !defined(BLOCK_DIRECT)
#define BLOCK_DIRECT 48
#endif

#define ROUND(a,b) ((a/b)*b)
#define min(a,b) (((a)<(b))?(a):(b))

#define EXACT do_exact_block_4x2

#include <stdlib.h>

static void do_block_2x2 (const int lda, const int M, const int N, const int K, const double* A_T, const double* B, double* restrict C);
static void do_exact_block_2x2 (const int lda, const double* A, const double* B, double* restrict C);
static void do_exact_block_3x2 (const int lda, const double* A, const double* B, double* restrict C);
static void do_exact_block_4x2 (const int lda, const double* A, const double* B, double* restrict C);
static void do_exact_block_8x1 (const int lda, const double* A, const double* B, double* restrict C);

#include "auto_blocks_inc.c"
//#include "auto_L1_inc.c"

static int lda_A;

/* This auxiliary subroutine performs a smaller dgemm operation
 *	C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block_2x2 (const int lda, const int M, const int N, const int K, const double* A_T, const double* B, double* restrict C) {
  int i, j, k;
  register double acc_00, acc_01, acc_10, acc_11;
  int mb = ROUND(M, 2);
  int nb = ROUND(N, 2);
  
  for (i = 0; i < mb; i+=2)
  {
    for (j = 0; j < nb; j+=2)
    {
      acc_00 = *(C + (i+0) + (j+0)*lda);
      acc_01 = *(C + (i+0) + (j+1)*lda);
      acc_10 = *(C + (i+1) + (j+0)*lda);
      acc_11 = *(C + (i+1) + (j+1)*lda);
      for (k = 0; k < K; k++)
      {
        acc_00 += (*(A_T + (i+0)*lda_A + k)) * (*(B + k + (j+0)*lda));
        acc_01 += (*(A_T + (i+0)*lda_A + k)) * (*(B + k + (j+1)*lda));
        acc_10 += (*(A_T + (i+1)*lda_A + k)) * (*(B + k + (j+0)*lda));
        acc_11 += (*(A_T + (i+1)*lda_A + k)) * (*(B + k + (j+1)*lda));
      }
      *(C + (i+0) + (j+0)*lda) = acc_00;
      *(C + (i+0) + (j+1)*lda) = acc_01;
      *(C + (i+1) + (j+0)*lda) = acc_10;
      *(C + (i+1) + (j+1)*lda) = acc_11;
    }
  }

// If any dim is odd, work on borders
  if ( (M!=mb) || (N!=nb) )
  {
    // note: LDA_ = LDA-1
  
    if (M!=mb)
    {
      // last row
      for (j = 0; j < nb; j+=2)
      {
        acc_00 = *(C + (mb) + (j+0)*lda);
        acc_01 = *(C + (mb) + (j+1)*lda);
        for (k = 0; k < K; k++)
        {
          acc_00 += (*(A_T + (mb)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_01 += (*(A_T + (mb)*lda_A + k)) * (*(B + k + (j+1)*lda));
        }
        *(C + (mb) + (j+0)*lda) = acc_00;
        *(C + (mb) + (j+1)*lda) = acc_01;
      }
    }
  
    // last column
    if (N!=nb)
    {
      for (i = 0; i < mb; i+=2)
      {
        acc_00 = *(C + (i+0) + (nb)*lda);
        acc_10 = *(C + (i+1) + (nb)*lda);
        for (k = 0; k < K; k++)
        {
          acc_00 += (*(A_T + (i+0)*lda_A + k)) * (*(B + k + (nb)*lda));
          acc_10 += (*(A_T + (i+1)*lda_A + k)) * (*(B + k + (nb)*lda));
        }
        *(C + (i+0) + (nb)*lda) = acc_00;
        *(C + (i+1) + (nb)*lda) = acc_10;
      }
    }
  
    // last element
    if ( (M!=mb) && (N!=nb) )
    {
      acc_00 = *(C + (mb) + (nb)*lda);
      for (k = 0; k < K; k++)
      {
        acc_00 += (*(A_T + (mb)*lda_A + k)) * (*(B + k + (nb)*lda));
      }
      *(C + (mb) + (nb)*lda) = acc_00;
    }
  }
}

/* This auxiliary subroutine performs a smaller dgemm operation
 *	C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N.
 * This function assumes that the submatrix is size BLOCK_SIZE
 */
static void do_exact_block_2x2 (const int lda, const double* A, const double* B, double* restrict C) {
  int i, j, k;
  register double acc_00, acc_01, acc_10, acc_11;
  
  for (i = 0; i < BLOCK_SIZE1; i+=2)
  {
    for (j = 0; j < BLOCK_SIZE1; j+=2)
    {
      acc_00 = *(C + (i+0) + (j+0)*lda);
      acc_01 = *(C + (i+0) + (j+1)*lda);
      acc_10 = *(C + (i+1) + (j+0)*lda);
      acc_11 = *(C + (i+1) + (j+1)*lda);
      for (k = 0; k < BLOCK_SIZE1; k++)
      {
        acc_00 += (*(A + (i+0)*lda_A + k)) * (*(B + k + (j+0)*lda));
        acc_01 += (*(A + (i+0)*lda_A + k)) * (*(B + k + (j+1)*lda));
        acc_10 += (*(A + (i+1)*lda_A + k)) * (*(B + k + (j+0)*lda));
        acc_11 += (*(A + (i+1)*lda_A + k)) * (*(B + k + (j+1)*lda));
      }
      *(C + (i+0) + (j+0)*lda) = acc_00;
      *(C + (i+0) + (j+1)*lda) = acc_01;
      *(C + (i+1) + (j+0)*lda) = acc_10;
      *(C + (i+1) + (j+1)*lda) = acc_11;
    }
  }

}

static void do_exact_block_3x2 (const int lda, const double* A, const double* B, double* restrict C)
{
  register double acc_00,acc_01,acc_10,acc_11,acc_20,acc_21;

  int i, j, k;
  
  for (i = 0; i < BLOCK_SIZE1; i+=3)
  {
    for (j = 0; j < BLOCK_SIZE1; j+=2)
    {
      acc_00 = *(C + (i+0) + (j+0)*lda);
      acc_01 = *(C + (i+0) + (j+1)*lda);
      acc_10 = *(C + (i+1) + (j+0)*lda);
      acc_11 = *(C + (i+1) + (j+1)*lda);
      acc_20 = *(C + (i+2) + (j+0)*lda);
      acc_21 = *(C + (i+2) + (j+1)*lda);
      for (k=0; k<BLOCK_SIZE1; k++)
      {
          acc_00 += (*(A + (i+0)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_01 += (*(A + (i+0)*lda_A + k)) * (*(B + k + (j+1)*lda));
          acc_10 += (*(A + (i+1)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_11 += (*(A + (i+1)*lda_A + k)) * (*(B + k + (j+1)*lda));
          acc_20 += (*(A + (i+2)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_21 += (*(A + (i+2)*lda_A + k)) * (*(B + k + (j+1)*lda));
      }
      *(C + (i+0) + (j+0)*lda) = acc_00;
      *(C + (i+0) + (j+1)*lda) = acc_01;
      *(C + (i+1) + (j+0)*lda) = acc_10;
      *(C + (i+1) + (j+1)*lda) = acc_11;
      *(C + (i+2) + (j+0)*lda) = acc_20;
      *(C + (i+2) + (j+1)*lda) = acc_21;
    }
  }
}

static void do_exact_block_4x2 (const int lda, const double* A, const double* B, double* restrict C)
{
  register double acc_00,acc_01,acc_10,acc_11,acc_20,acc_21,acc_30,acc_31;

  int i, j, k;
  
  for (i = 0; i < BLOCK_SIZE1; i+=4)
  {
    for (j = 0; j < BLOCK_SIZE1; j+=2)
    {
      acc_00 = *(C + (i+0) + (j+0)*lda);
      acc_01 = *(C + (i+0) + (j+1)*lda);
      acc_10 = *(C + (i+1) + (j+0)*lda);
      acc_11 = *(C + (i+1) + (j+1)*lda);
      acc_20 = *(C + (i+2) + (j+0)*lda);
      acc_21 = *(C + (i+2) + (j+1)*lda);
      acc_30 = *(C + (i+3) + (j+0)*lda);
      acc_31 = *(C + (i+3) + (j+1)*lda);
      for (k=0; k<BLOCK_SIZE1; k++)
      {
          acc_00 += (*(A + (i+0)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_01 += (*(A + (i+0)*lda_A + k)) * (*(B + k + (j+1)*lda));
          acc_10 += (*(A + (i+1)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_11 += (*(A + (i+1)*lda_A + k)) * (*(B + k + (j+1)*lda));
          acc_20 += (*(A + (i+2)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_21 += (*(A + (i+2)*lda_A + k)) * (*(B + k + (j+1)*lda));
          acc_30 += (*(A + (i+3)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_31 += (*(A + (i+3)*lda_A + k)) * (*(B + k + (j+1)*lda));
      }
      *(C + (i+0) + (j+0)*lda) = acc_00;
      *(C + (i+0) + (j+1)*lda) = acc_01;
      *(C + (i+1) + (j+0)*lda) = acc_10;
      *(C + (i+1) + (j+1)*lda) = acc_11;
      *(C + (i+2) + (j+0)*lda) = acc_20;
      *(C + (i+2) + (j+1)*lda) = acc_21;
      *(C + (i+3) + (j+0)*lda) = acc_30;
      *(C + (i+3) + (j+1)*lda) = acc_31;
    }
  }
}

static void do_exact_block_8x1 (const int lda, const double* A, const double* B, double* restrict C)
{
  register double acc_00,acc_10,acc_20,acc_30,acc_40,acc_50,acc_60,acc_70;

  int i, j, k;
  
  for (i = 0; i < BLOCK_SIZE1; i+=8)
  {
    for (j = 0; j < BLOCK_SIZE1; j+=1)
    {
      acc_00 = *(C + (i+0) + (j+0)*lda);
      acc_10 = *(C + (i+1) + (j+0)*lda);
      acc_20 = *(C + (i+2) + (j+0)*lda);
      acc_30 = *(C + (i+3) + (j+0)*lda);
      acc_40 = *(C + (i+4) + (j+0)*lda);
      acc_50 = *(C + (i+5) + (j+0)*lda);
      acc_60 = *(C + (i+6) + (j+0)*lda);
      acc_70 = *(C + (i+7) + (j+0)*lda);
      for (k=0; k<BLOCK_SIZE1; k++)
      {
          acc_00 += (*(A + (i+0)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_10 += (*(A + (i+1)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_20 += (*(A + (i+2)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_30 += (*(A + (i+3)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_40 += (*(A + (i+4)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_50 += (*(A + (i+5)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_60 += (*(A + (i+6)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_70 += (*(A + (i+7)*lda_A + k)) * (*(B + k + (j+0)*lda));
      }
      *(C + (i+0) + (j+0)*lda) = acc_00;
      *(C + (i+1) + (j+0)*lda) = acc_10;
      *(C + (i+2) + (j+0)*lda) = acc_20;
      *(C + (i+3) + (j+0)*lda) = acc_30;
      *(C + (i+4) + (j+0)*lda) = acc_40;
      *(C + (i+5) + (j+0)*lda) = acc_50;
      *(C + (i+6) + (j+0)*lda) = acc_60;
      *(C + (i+7) + (j+0)*lda) = acc_70;
    }
  }
}

static void L1_dgemm (const int lda, const int I, const int J, const int K, const double* A_T, const double* B, double* restrict C) {
	/* For each block-column of B */
	for (int j = 0; j < J; j += BLOCK_SIZE1)
		/* For each block-row of A */ 
		for (int i = 0; i < I; i += BLOCK_SIZE1)
			/* Accumulate block dgemms into block of C */
			for (int k = 0; k < K; k += BLOCK_SIZE1) {
				if ( (i<=(I-BLOCK_SIZE1)) && (j<=(J-BLOCK_SIZE1)) && (k<=(K-BLOCK_SIZE1))){
					EXACT (lda, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
				} else {
					/* Correct block dimensions if block "goes off edge of" the matrix */
					int I_ = min (BLOCK_SIZE1, I-i);
					int J_ = min (BLOCK_SIZE1, J-j);
					int K_ = min (BLOCK_SIZE1, K-k);
					/* Perform individual block dgemm */
					do_block_2x2(lda, I_, J_, K_, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
				}
			}
}

//size of submatrix == BLOCK_SIZE2
static void L1_dgemm_exact (const int lda, const double* A_T, const double* B, double* restrict C) {
	/* For each block-column of B */
	for (int j = 0; j < BLOCK_SIZE2; j += BLOCK_SIZE1)
		/* For each block-row of A */ 
		for (int i = 0; i < BLOCK_SIZE2; i += BLOCK_SIZE1)
			/* Accumulate block dgemms into block of C */
			for (int k = 0; k < BLOCK_SIZE2; k += BLOCK_SIZE1) {
#ifdef ALIGNED_BLOCKS
				EXACT (lda, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
#else
				if ( (i<=(BLOCK_SIZE2-BLOCK_SIZE1)) && (j<=(BLOCK_SIZE2-BLOCK_SIZE1)) && (k<=(BLOCK_SIZE2-BLOCK_SIZE1)) ){
					do_exact_block_4x2(lda, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
				} else {
					/* Correct block dimensions if block "goes off edge of" the matrix */
					int I_ = min (BLOCK_SIZE1, BLOCK_SIZE2-i);
					int J_ = min (BLOCK_SIZE1, BLOCK_SIZE2-j);
					int K_ = min (BLOCK_SIZE1, BLOCK_SIZE2-k);
					/* Perform individual block dgemm */
					do_block_2x2(lda, I_, J_, K_, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
				}
#endif
			}
}

static void L2_dgemm (const int lda, const double* A_T, const double* B, double* restrict C) {
	/* For each block-column of B */
	for (int j = 0; j < lda; j += BLOCK_SIZE2)
		/* For each block-row of A */ 
		for (int i = 0; i < lda; i += BLOCK_SIZE2)
			/* Accumulate block dgemms into block of C */
			for (int k = 0; k < lda; k += BLOCK_SIZE2) {
				if ( (i<=(lda-BLOCK_SIZE2)) && (j<=(lda-BLOCK_SIZE2)) && (k<=(lda-BLOCK_SIZE2)) ){
					//do_exact_block(lda, A + i + k*lda, B + k + j*lda, C + i + j*lda);
					L1_dgemm_exact(lda, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
				} else {
					/* Correct block dimensions if block "goes off edge of" the matrix */
					int I_ = min (BLOCK_SIZE2, lda-i);
					int J_ = min (BLOCK_SIZE2, lda-j);
					int K_ = min (BLOCK_SIZE2, lda-k);
					/* Perform individual block dgemm */
					L1_dgemm(lda, I_, J_, K_, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
				}
			}
}

/* This routine performs a dgemm operation
 *	C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format. 
 * On exit, A and B maintain their input values. */	
void square_dgemm (int lda, const double* A, const double* B, double* restrict C) {
	double* restrict A_T;

	//transpose A
#ifdef TRANSPOSE
	lda_A = ((lda+1)>>1)<<1;
	posix_memalign((void**) &A_T, 16, sizeof(double)*lda_A*lda);
	for (int i=0; i<lda; i++) {
		for (int j=i; j<lda; j++) {
			A_T[i + lda_A*j] = A[i*lda + j];
			A_T[j + lda_A*i] = A[j*lda + i];
		}
	}
#endif

	if (lda<=BLOCK_DIRECT) {
		exact_blocks[lda](A_T, B, C);
	} else if (lda<=BLOCK_SIZE2){
		//L1_blocks[lda](A_T, B, C);
		L1_dgemm(lda, lda, lda, lda, A_T, B, C);
	} else {
		L2_dgemm(lda, A_T, B, C);
	}

	free(A_T);
}

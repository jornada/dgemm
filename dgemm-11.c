#include <stdlib.h>
const char* dgemm_desc = "dgemm - group #11";
static int lda_A;

#define ALIGNED_BLOCKS
#define EXACT_SMALL
#define EXACT_L1

#if !defined(BLOCK_SIZE2)
#define BLOCK_SIZE2 416
#endif

#if !defined(BLOCK_SIZE1)
#define BLOCK_SIZE1 208
#endif

#if !defined(BLOCK_DIRECT)
#define BLOCK_DIRECT 64
#endif

#define ROUND(a,b) ((a/b)*b)
#define min(a,b) (((a)<(b))?(a):(b))

static void do_block_2x2 (const int lda, const int M, const int N, const int K, const double* A_T, const double* B, double* restrict C);
static void do_block_4x2 (const int lda, const int M, const int N, const int K, const double* A_T, const double* B, double* restrict C);
static void do_block_8x1 (const int lda, const int M, const int N, const int K, const double* A_T, const double* B, double* restrict C);
static void do_exact_block_2x2 (const int lda, const double* A, const double* B, double* restrict C);
static void do_exact_block_3x2 (const int lda, const double* A, const double* B, double* restrict C);
static void do_exact_block_4x2 (const int lda, const double* A, const double* B, double* restrict C);
static void do_exact_block_8x1 (const int lda, const double* A, const double* B, double* restrict C);

#define EXACT(a,b,c,d)\
	((lda_A%256)==0)? do_exact_block_8x1(a,b,c,d) : do_exact_block_4x2(a,b,c,d)
#define INEXACT(a,b,c,d,e,f,g)\
	((lda_A%256)==0)? do_block_8x1(a,b,c,d,e,f,g) : do_block_4x2(a,b,c,d,e,f,g)

#ifdef EXACT_SMALL
#include "auto_blocks_inc.c"
#endif

#ifdef EXACT_L1
#include "auto_L1_inc.c"
#endif

/* This auxiliary subroutine performs a smaller dgemm operation
 *	C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block_2x2 (const int lda, const int M, const int N, const int K, const double* A_T, const double* B, double* restrict C) {
  register double acc_00, acc_01, acc_10, acc_11;
  int i, j, k;
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

  // last row
  if (M!=mb)
  {
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

static void do_block_8x1 (const int lda, const int M, const int N, const int K, const double* A_T, const double* B, double* restrict C) {
  register double acc_00,acc_10,acc_20,acc_30,acc_40,acc_50,acc_60,acc_70;
  int i, j, k;
  int mb = ROUND(M, 8);
  
  for (i = 0; i < mb; i+=8)
  {
    for (j = 0; j < N; j+=1)
    {
      acc_00 = *(C + (i+0) + (j+0)*lda);
      acc_10 = *(C + (i+1) + (j+0)*lda);
      acc_20 = *(C + (i+2) + (j+0)*lda);
      acc_30 = *(C + (i+3) + (j+0)*lda);
      acc_40 = *(C + (i+4) + (j+0)*lda);
      acc_50 = *(C + (i+5) + (j+0)*lda);
      acc_60 = *(C + (i+6) + (j+0)*lda);
      acc_70 = *(C + (i+7) + (j+0)*lda);
      for (k=0; k<K; k++)
      {
          acc_00 += (*(A_T + (i+0)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_10 += (*(A_T + (i+1)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_20 += (*(A_T + (i+2)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_30 += (*(A_T + (i+3)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_40 += (*(A_T + (i+4)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_50 += (*(A_T + (i+5)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_60 += (*(A_T + (i+6)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_70 += (*(A_T + (i+7)*lda_A + k)) * (*(B + k + (j+0)*lda));
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
  
  // last rows
  if (M!=mb)
  {
    for (i = mb; i < M; i++)
    {
      for (j = 0; j < N; j++)
      {
        acc_00 = *(C + (i) + (j+0)*lda);
        for (k = 0; k < K; k++)
        {
          acc_00 += (*(A_T + (i)*lda_A + k)) * (*(B + k + (j+0)*lda));
        }
        *(C + (i) + (j+0)*lda) = acc_00;
      }
    }
  }
}

static void do_block_4x2 (const int lda, const int M, const int N, const int K, const double* A_T, const double* B, double* restrict C) {
  register double acc_00,acc_01,acc_10,acc_11,acc_20,acc_21,acc_30,acc_31;
  int i, j, k;
  int mb = ROUND(M, 4);
  int nb = ROUND(N, 2);
  
  for (i = 0; i < mb; i+=4)
  {
    for (j = 0; j < nb; j+=2)
    {
      acc_00 = *(C + (i+0) + (j+0)*lda);
      acc_01 = *(C + (i+0) + (j+1)*lda);
      acc_10 = *(C + (i+1) + (j+0)*lda);
      acc_11 = *(C + (i+1) + (j+1)*lda);
      acc_20 = *(C + (i+2) + (j+0)*lda);
      acc_21 = *(C + (i+2) + (j+1)*lda);
      acc_30 = *(C + (i+3) + (j+0)*lda);
      acc_31 = *(C + (i+3) + (j+1)*lda);
      for (k=0; k<K; k++)
      {
          acc_00 += (*(A_T + (i+0)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_01 += (*(A_T + (i+0)*lda_A + k)) * (*(B + k + (j+1)*lda));
          acc_10 += (*(A_T + (i+1)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_11 += (*(A_T + (i+1)*lda_A + k)) * (*(B + k + (j+1)*lda));
          acc_20 += (*(A_T + (i+2)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_21 += (*(A_T + (i+2)*lda_A + k)) * (*(B + k + (j+1)*lda));
          acc_30 += (*(A_T + (i+3)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_31 += (*(A_T + (i+3)*lda_A + k)) * (*(B + k + (j+1)*lda));
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

  // last rows
  if (M!=mb)
  {
    for (i = mb; i < M; i++)
    {
      for (j = 0; j < nb; j+=2)
      {
        acc_00 = *(C + (i) + (j+0)*lda);
        acc_01 = *(C + (i) + (j+1)*lda);
        for (k = 0; k < K; k++)
        {
          acc_00 += (*(A_T + (i)*lda_A + k)) * (*(B + k + (j+0)*lda));
          acc_01 += (*(A_T + (i)*lda_A + k)) * (*(B + k + (j+1)*lda));
        }
        *(C + (i) + (j+0)*lda) = acc_00;
        *(C + (i) + (j+1)*lda) = acc_01;
      }
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

  // last elements
  if ( (M!=mb) && (N!=nb) )
  {
    for (i = mb; i < M; i++)
    {
      acc_00 = *(C + (i) + (nb)*lda);
      for (k = 0; k < K; k++)
      {
        acc_00 += (*(A_T + (i)*lda_A + k)) * (*(B + k + (nb)*lda));
      }
    *(C + (i) + (nb)*lda) = acc_00;
    }
  }
}

/* This auxiliary subroutine performs a smaller dgemm operation
 *	C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N.
 * This function assumes that the submatrix is size BLOCK_SIZE
 */
static void do_exact_block_2x2 (const int lda, const double* A, const double* B, double* restrict C) {
  register double acc_00, acc_01, acc_10, acc_11;
  int i, j, k;
  
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
					//INEXACT (lda, BLOCK_SIZE1, BLOCK_SIZE1, BLOCK_SIZE1, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
					EXACT (lda, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
				} else {
					/* Correct block dimensions if block "goes off edge of" the matrix */
					int I_ = min (BLOCK_SIZE1, I-i);
					int J_ = min (BLOCK_SIZE1, J-j);
					int K_ = min (BLOCK_SIZE1, K-k);
					/* Perform individual block dgemm */
					INEXACT (lda, I_, J_, K_, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
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
				//INEXACT (lda, BLOCK_SIZE1, BLOCK_SIZE1, BLOCK_SIZE1, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
				EXACT (lda, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
#else
				if ( (i<=(BLOCK_SIZE2-BLOCK_SIZE1)) && (j<=(BLOCK_SIZE2-BLOCK_SIZE1)) && (k<=(BLOCK_SIZE2-BLOCK_SIZE1)) ){
					EXACT (lda, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
				} else {
					/* Correct block dimensions if block "goes off edge of" the matrix */
					int I_ = min (BLOCK_SIZE1, BLOCK_SIZE2-i);
					int J_ = min (BLOCK_SIZE1, BLOCK_SIZE2-j);
					int K_ = min (BLOCK_SIZE1, BLOCK_SIZE2-k);
					/* Perform individual block dgemm */
					INEXACT (lda, I_, J_, K_, A_T + k + i*lda_A, B + k + j*lda, C + i + j*lda);
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
	lda_A = ((lda+1)>>1)<<1;
	posix_memalign((void**) &A_T, 16, sizeof(double)*lda_A*lda);
	for (int i=0; i<lda; i++) {
		for (int j=i; j<lda; j++) {
			A_T[i + lda_A*j] = A[i*lda + j];
			A_T[j + lda_A*i] = A[j*lda + i];
		}
	}

	if (lda<=BLOCK_DIRECT) {
#ifdef EXACT_SMALL
		exact_blocks[lda](A_T, B, C);
#else
		do_block_4x2(lda, lda, lda, lda, A_T, B, C);
#endif
	} else if (lda<=BLOCK_SIZE2){
#ifdef EXACT_L1
		L1_blocks[lda](A_T, B, C);
#else
		L1_dgemm(lda, lda, lda, lda, A_T, B, C);
#endif
	} else {
		L2_dgemm(lda, A_T, B, C);
	}

	free(A_T);
}

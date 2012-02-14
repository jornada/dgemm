//#ifdef BORDER
static void LABEL (const double* A_T, const double* B, double* restrict C) {
  int i, j, k;
  register double acc_00, acc_01, acc_10, acc_11;
  
  for (i = 0; i < LDA_; i+=2)
  {
    for (j = 0; j < LDA_; j+=2)
    {
      acc_00 = *(C + (i+0) + (j+0)*LDA);
      acc_01 = *(C + (i+0) + (j+1)*LDA);
      acc_10 = *(C + (i+1) + (j+0)*LDA);
      acc_11 = *(C + (i+1) + (j+1)*LDA);
      for (k = 0; k < LDA; k++)
      {
        acc_00 += (*(A_T + (i+0)*LDA_A + k)) * (*(B + k + (j+0)*LDA));
        acc_01 += (*(A_T + (i+0)*LDA_A + k)) * (*(B + k + (j+1)*LDA));
        acc_10 += (*(A_T + (i+1)*LDA_A + k)) * (*(B + k + (j+0)*LDA));
        acc_11 += (*(A_T + (i+1)*LDA_A + k)) * (*(B + k + (j+1)*LDA));
      }
      *(C + (i+0) + (j+0)*LDA) = acc_00;
      *(C + (i+0) + (j+1)*LDA) = acc_01;
      *(C + (i+1) + (j+0)*LDA) = acc_10;
      *(C + (i+1) + (j+1)*LDA) = acc_11;
    }
  }

// If dim is odd, work on borders
  // note: LDA_ = LDA-1
#ifdef BORDER
  // last row
  for (j = 0; j < LDA_; j+=2)
  {
    acc_00 = *(C + (LDA_) + (j+0)*LDA);
    acc_01 = *(C + (LDA_) + (j+1)*LDA);
    for (k = 0; k < LDA; k++)
    {
      acc_00 += (*(A_T + (LDA_)*LDA_A + k)) * (*(B + k + (j+0)*LDA));
      acc_01 += (*(A_T + (LDA_)*LDA_A + k)) * (*(B + k + (j+1)*LDA));
    }
    *(C + (LDA_) + (j+0)*LDA) = acc_00;
    *(C + (LDA_) + (j+1)*LDA) = acc_01;
  }

  // last column
  for (i = 0; i < LDA_; i+=2)
  {
    acc_00 = *(C + (i+0) + (LDA_)*LDA);
    acc_10 = *(C + (i+1) + (LDA_)*LDA);
    for (k = 0; k < LDA; k++)
    {
      acc_00 += (*(A_T + (i+0)*LDA_A + k)) * (*(B + k + (LDA_)*LDA));
      acc_10 += (*(A_T + (i+1)*LDA_A + k)) * (*(B + k + (LDA_)*LDA));
    }
    *(C + (i+0) + (LDA_)*LDA) = acc_00;
    *(C + (i+1) + (LDA_)*LDA) = acc_10;
  }

  // last element
  acc_00 = *(C + (LDA_) + (LDA_)*LDA);
  for (k = 0; k < LDA; k++)
  {
    acc_00 += (*(A_T + (LDA_)*LDA_A + k)) * (*(B + k + (LDA_)*LDA));
  }
  *(C + (LDA_) + (LDA_)*LDA) = acc_00;
#endif
}
/*
#else
static void LABEL (const double* A, const double* B, double* restrict C)
{
  register double acc_00,acc_01,acc_10,acc_11,acc_20,acc_21,acc_30,acc_31;

  int i, j, k;
  
  for (i = 0; i < LDA_; i+=4)
  {
    for (j = 0; j < LDA_; j+=2)
    {
      acc_00 = *(C + (i+0) + (j+0)*LDA);
      acc_01 = *(C + (i+0) + (j+1)*LDA);
      acc_10 = *(C + (i+1) + (j+0)*LDA);
      acc_11 = *(C + (i+1) + (j+1)*LDA);
      acc_20 = *(C + (i+2) + (j+0)*LDA);
      acc_21 = *(C + (i+2) + (j+1)*LDA);
      acc_30 = *(C + (i+3) + (j+0)*LDA);
      acc_31 = *(C + (i+3) + (j+1)*LDA);
      for (k=0; k<LDA; k++)
      {
          acc_00 += (*(A + (i+0)*LDA_A + k)) * (*(B + k + (j+0)*LDA));
          acc_01 += (*(A + (i+0)*LDA_A + k)) * (*(B + k + (j+1)*LDA));
          acc_10 += (*(A + (i+1)*LDA_A + k)) * (*(B + k + (j+0)*LDA));
          acc_11 += (*(A + (i+1)*LDA_A + k)) * (*(B + k + (j+1)*LDA));
          acc_20 += (*(A + (i+2)*LDA_A + k)) * (*(B + k + (j+0)*LDA));
          acc_21 += (*(A + (i+2)*LDA_A + k)) * (*(B + k + (j+1)*LDA));
          acc_30 += (*(A + (i+3)*LDA_A + k)) * (*(B + k + (j+0)*LDA));
          acc_31 += (*(A + (i+3)*LDA_A + k)) * (*(B + k + (j+1)*LDA));
      }
      *(C + (i+0) + (j+0)*LDA) = acc_00;
      *(C + (i+0) + (j+1)*LDA) = acc_01;
      *(C + (i+1) + (j+0)*LDA) = acc_10;
      *(C + (i+1) + (j+1)*LDA) = acc_11;
      *(C + (i+2) + (j+0)*LDA) = acc_20;
      *(C + (i+2) + (j+1)*LDA) = acc_21;
      *(C + (i+3) + (j+0)*LDA) = acc_30;
      *(C + (i+3) + (j+1)*LDA) = acc_31;
    }
  }
}
#endif
*/

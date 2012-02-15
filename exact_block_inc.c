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

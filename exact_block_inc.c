
//Original version
#if 0
static void LABEL (double* restrict A_T, double* restrict B, double* restrict C) {
	/* For each column j of B */ 
	int j_lda = 0;
	for (int j = 0; j < LDA; ++j, j_lda+=LDA) {
		/* For each row i of A */
		int i_lda = 0;
		int ij_lda = j_lda;
		for (int i = 0; i < LDA; ++i, i_lda+=LDA_A, ij_lda++) {
			/* Compute C(i,j) */
			double cij = C[ij_lda];
			int ki_lda = i_lda;
			int kj_lda = j_lda;
			for (int k = 0; k < LDA; k++) {
				cij += A_T[ki_lda++] * B[kj_lda++];
			}
			C[ij_lda] = cij;
		}
	}
}
#else
static void LABEL (double* restrict A_T, double* restrict B, double* restrict C) {
	/* For each column j of B */ 
	for (int j = 0; j < LDA; j++) {
		/* For each row i of A */
		for (int i = 0; i < LDA; i++) {
			/* Compute C(i,j) */
			double cij = C[i + j*LDA];
			for (int k = 0; k < LDA; k++) {
				cij += A_T[i*LDA_A + k] * B[j*LDA + k];
			}
			C[i + j*LDA] = cij;
		}
	}
}
#endif


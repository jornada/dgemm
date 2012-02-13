
static void LABEL (int lda, double* restrict A_T, double* restrict B, double* restrict C) {
	/* For each column j of B */ 
	int j_lda = 0;
	for (int j = 0; j < FUNC_SIZE; ++j, j_lda+=lda) {
		/* For each row i of A */
		int i_lda = 0;
		int ij_lda = j_lda;
		for (int i = 0; i < FUNC_SIZE; ++i, i_lda+=lda, ij_lda++) {
			/* Compute C(i,j) */
			double cij = C[ij_lda];
			int ki_lda = i_lda;
			int kj_lda = j_lda;
			for (int k = 0; k < FUNC_SIZE; ++k) {
				cij += A_T[ki_lda++] * B[kj_lda++];
			}
			C[ij_lda] = cij;
		}
	}
}



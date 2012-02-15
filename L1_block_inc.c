static void LABEL (const double* A_T, const double* B, double* restrict C) {
	/* For each block-column of B */
	for (int j = 0; j < LDA_; j += BLOCK_SIZE1) {
		/* For each block-row of A */ 
		for (int i = 0; i < LDA_; i += BLOCK_SIZE1) {
			/* Accumulate block dgemms into block of C */
			for (int k = 0; k < LDA_; k += BLOCK_SIZE1) {
				do_exact_block_4x2(LDA, A_T + k + i*LDA_A, B + k + j*LDA, C + i + j*LDA);
			}
		}
	}

#ifdef BORDER
	//fix k
	for (int j = 0; j < LDA_; j += BLOCK_SIZE1) {
		for (int i = 0; i < LDA_; i += BLOCK_SIZE1) {
			do_block_4x2(LDA, BLOCK_SIZE1, BLOCK_SIZE1, LDA_BORDER, \
				A_T + LDA_ + i*LDA_A, B + LDA_ + j*LDA, C + i + j*LDA);
		}
	}
	//fix j
	for (int i = 0; i < LDA_; i += BLOCK_SIZE1) {
		for (int k = 0; k < LDA_; k += BLOCK_SIZE1) {
			do_block_4x2(LDA, BLOCK_SIZE1, LDA_BORDER, BLOCK_SIZE1, \
				A_T + k + i*LDA_A, B + k + LDA_*LDA, C + i + LDA_*LDA);
		}
	}
	//fix i
	for (int j = 0; j < LDA_; j += BLOCK_SIZE1) {
		for (int k = 0; k < LDA_; k += BLOCK_SIZE1) {
			do_block_4x2(LDA, LDA_BORDER, BLOCK_SIZE1, BLOCK_SIZE1, \
				A_T + k + LDA_*LDA_A, B + k + j*LDA, C + LDA_ + j*LDA);
		}
	}

	//fix kj
	for (int i = 0; i < LDA_; i += BLOCK_SIZE1) {
		do_block_4x2(LDA, BLOCK_SIZE1, LDA_BORDER, LDA_BORDER, \
			A_T + LDA_ + i*LDA_A, B + LDA_ + LDA_*LDA, C + i + LDA_*LDA);
	}
	//fix ki
	for (int j = 0; j < LDA_; j += BLOCK_SIZE1) {
		do_block_4x2(LDA, LDA_BORDER, BLOCK_SIZE1, LDA_BORDER, \
			A_T + LDA_ + LDA_*LDA_A, B + LDA_ + j*LDA, C + LDA_ + j*LDA);
	}
	//fix ij
	for (int k = 0; k < LDA_; k += BLOCK_SIZE1) {
		do_block_4x2(LDA, LDA_BORDER, LDA_BORDER, BLOCK_SIZE1, \
			A_T + k + LDA_*LDA_A, B + k + LDA_*LDA, C + LDA_ + LDA_*LDA);
	}

	//fix ijk
	do_block_4x2(LDA, LDA_BORDER, LDA_BORDER, LDA_BORDER, \
		A_T + LDA_ + LDA_*LDA_A, B + LDA_ + LDA_*LDA, C + LDA_ + LDA_*LDA);

#endif
}

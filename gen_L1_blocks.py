#!/usr/bin/env python

import sys, re

f=open('dgemm-11.c')
str = f.read()
f.close()

obj = re.search('#define\s*BLOCK_DIRECT\s*([0-9]+)', str)
if obj is None:
	print 'Could not determine BLOCK_DIRECT'
	sys.exit()
BLOCK_DIRECT = int(obj.group(1))

obj = re.search('#define\s*BLOCK_SIZE1\s*([0-9]+)', str)
if obj is None:
	print 'Could not determine BLOCK_SIZE1'
	sys.exit()
BLOCK_SIZE1 = int(obj.group(1))

obj = re.search('#define\s*BLOCK_SIZE2\s*([0-9]+)', str)
if obj is None:
	print 'Could not determine BLOCK_SIZE2'
	sys.exit()
BLOCK_SIZE2 = int(obj.group(1))

f=open('auto_L1_inc.c', 'w')

for bs in range(BLOCK_DIRECT, BLOCK_SIZE2+1):
	#size of transposed matrix
	bs_A = ((bs+1)//2)*2
	#bs rounded down by BLOCK_SIZE1
	bs_ = ((bs)//BLOCK_SIZE1)*BLOCK_SIZE1
	#border used in L1_*
	border = bs-bs_;
	#bs_A = bs
	if (border): f.write('#define BORDER\n')
	f.write('#define LDA %i\n'%(bs))
	f.write('#define LDA_BORDER %i\n'%(border))
	f.write('#define LDA_ %i\n'%(bs_))
	f.write('#define LDA_A %i\n'%(bs_A))
	f.write('#define LABEL L1_%i\n'%(bs))
	f.write('#include "L1_block_inc.c"\n')
	f.write('#undef LABEL\n')
	f.write('#undef LDA\n')
	f.write('#undef LDA_\n')
	f.write('#undef LDA_BORDER\n')
	f.write('#undef LDA_A\n')
	if (border): f.write('#undef BORDER\n')
	f.write('\n')

f.write('static void (*L1_blocks[%i])(const double*, const double*, double* restrict)={'%(BLOCK_SIZE2+1));
f.write('NULL');
for bs in range(1,BLOCK_DIRECT):
	f.write(', NULL');
for bs in range(BLOCK_DIRECT,BLOCK_SIZE2+1):
	f.write(', L1_%i'%(bs));

f.write('};\n');
f.close()


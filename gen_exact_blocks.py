#!/usr/bin/env python

import sys, re

f=open('dgemm-jornada-auto.c')
obj = re.search('#define\s*BLOCK_SIZE1\s*([0-9]+)', f.read())
f.close()

if obj is None:
	print 'Could not determine BLOCK_SIZE1'
	sys.exit()


BLOCK_SIZE1 = int(obj.group(1))

f=open('auto_blocks_inc.c', 'w')

for bs in range(1,BLOCK_SIZE1):
	bs_A = ((bs+1)//2)*2
	#bs_A = bs
	f.write('#define LDA %i\n'%(bs))
	f.write('#define LDA_A %i\n'%(bs_A))
	f.write('#define LABEL block_%i\n'%(bs))
	f.write('#include "exact_block_inc.c"\n')
	f.write('#undef LABEL\n')
	f.write('#undef LDA\n\n')
	f.write('#undef LDA_A\n\n')

f.write('static void (*exact_blocks[%i])(double* restrict, double* restrict, double* restrict)={'%(BLOCK_SIZE1));
f.write('NULL');
for bs in range(1,BLOCK_SIZE1):
	#f_init.write('exact_blocks[%i] = block_%i;\n'%(bs,bs));
	f.write(', block_%i'%(bs));

f.write('};\n');
f.close()


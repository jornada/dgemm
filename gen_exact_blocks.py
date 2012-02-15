#!/usr/bin/env python

import sys, re

f=open('dgemm-11.c')
obj = re.search('#define\s*BLOCK_DIRECT\s*([0-9]+)', f.read())
f.close()

if obj is None:
	print 'Could not determine BLOCK_DIRECT'
	sys.exit()


BLOCK_DIRECT = int(obj.group(1))

f=open('auto_blocks_inc.c', 'w')

for bs in range(1,BLOCK_DIRECT+1):
	#size of transposed matrix
	bs_A = ((bs+1)//2)*2
	#bs rounded down, used for register blocking
	bs_ = ((bs)//2)*2
	#bs_A = bs
	if (bs%2): f.write('#define BORDER\n')
	f.write('#define LDA %i\n'%(bs))
	f.write('#define LDA_ %i\n'%(bs_))
	f.write('#define LDA_A %i\n'%(bs_A))
	f.write('#define LABEL block_%i\n'%(bs))
	f.write('#include "exact_block_inc.c"\n')
	f.write('#undef LABEL\n')
	f.write('#undef LDA\n')
	f.write('#undef LDA_\n')
	f.write('#undef LDA_A\n')
	if (bs%2): f.write('#undef BORDER\n')
	f.write('\n')

f.write('static void (*exact_blocks[%i])(const double*, const double*, double* restrict)={'%(BLOCK_DIRECT+1));
f.write('NULL');
for bs in range(1,BLOCK_DIRECT+1):
	#f_init.write('exact_blocks[%i] = block_%i;\n'%(bs,bs));
	f.write(', block_%i'%(bs));

f.write('};\n');
f.close()


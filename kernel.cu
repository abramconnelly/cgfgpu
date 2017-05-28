
#include <stdio.h>


extern "C" __global__ void kernel ( int cnt, int* buf  )
{
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	if ( x > cnt ) return;

	buf[x] = 1000-x;

}

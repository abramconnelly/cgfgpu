

#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>

CUdevice		cuDevice;
CUcontext		cuContext;
CUmodule		cuCustom;
CUfunction		cuKernel;

bool cudaCheck ( CUresult status, char* msg )
{
	if ( status != CUDA_SUCCESS ) {
		const char* stat = "";
		cuGetErrorString ( status, &stat );
		printf ( "CUDA ERROR: %s (in %s)\n", stat, msg  );	
		exit(-1);
		return false;
	} 
	return true;
}

void cudaStart ()
{
	int version = 0;
    char name[100];
	int count;
    
	cuInit ( 0 );
	cudaCheck ( cuDeviceGet ( &cuDevice, 0 ), "cuDeviceGet" );
	cudaCheck ( cuDeviceGetName ( name, 100, cuDevice ), "cuDeviceGetName" );
	cudaCheck ( cuCtxCreate ( &cuContext, 0, cuDevice ), "cuCtxCreate" );

	printf ( "CUDA Device: %s\n", name );	
			
}

int main (int argc, char* argv)
{
	cudaStart ();

	// Load CUDA module and kernel
	printf ( "Loading module\n");
	cudaCheck ( cuModuleLoad ( &cuCustom, "kernel.ptx" ), "cuModuleLoad (kernel.ptx)" );
	cudaCheck ( cuModuleGetFunction ( &cuKernel, cuCustom, "kernel" ), "cuModuleGetFunction (kernel)" );	

	int cnt = 100000;	
	int threads = 256;
	int blocks = int(cnt / threads)+1;

	CUdeviceptr gmem;
	cudaCheck ( cuMemAlloc ( &gmem, cnt*sizeof(int) ), "cuMemAlloc" );
	
	void* args[2] = { &cnt, &gmem };
	cudaCheck ( cuLaunchKernel ( cuKernel, threads, 1, 1, blocks, 1, 1, 0, NULL, args, NULL ), "cuLaunch(user_kernel)" );
	
	// retrieve results
	int* dat = (int*) malloc ( cnt * sizeof(int) );
	cudaCheck ( cuMemcpyDtoH ( dat, gmem, cnt*sizeof(int) ), "cuMemcpyDtoH" );
	for (int n=0; n < 100; n++ ) {
		printf ( "%d: %d\n", n, dat[n] );
	}

 	return 1;
}



#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>

#include "cgft.hpp"

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

	cgf_t* cgf1;
	cgf_t* cgf2;
	
	FILE* fp1 = fopen ( "hg19.cgfv3", "rb" );
	if ( fp1 == 0x0 ) printf ( "Cannot find cgf.\n" );
	cgf1 = cgft_read ( fp1 );
	if ( cgf1 == 0x0 ) printf ( "Cannot read cgf.\n" );

	FILE* fp2 = fopen ( "hu34D5B9-GS01670-DNA_E02.cgfv3", "rb" );
	if ( fp2 == 0x0 ) printf ( "Cannot find cgf.\n" );
	cgf2 = cgft_read ( fp2 );
	if ( cgf2 == 0x0 ) printf ( "Cannot read cgf.\n" );

	
	FILE* fout = fopen ( "out.txt", "wt" );
	tilepath_t tp = cgf1->Path[0];
	cgft_output_band_format( cgf1, &tp, fout );
	tp = cgf2->Path[0];
	cgft_output_band_format( cgf2, &tp, fout );
	fclose (fout);

 	return 1;
}

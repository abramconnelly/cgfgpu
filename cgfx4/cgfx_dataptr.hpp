
#ifndef DEF_DATAPTR
	#define DEF_DATAPTR

	#include <stdint.h>

	#ifdef WIN32
		#include <unistd_win.h>
		#include <getopt_win.h>
		#include <time.h>


		typedef uint32_t			uint32;
		typedef uint64_t			uint64;
		typedef int16_t				sint16;
		typedef int32_t				sint32;
		typedef int64_t				sint64;
		typedef unsigned char		byte;
		typedef unsigned char		uchar;
		typedef signed char			schar;
		typedef uint16_t			ushort;
		typedef uint32_t			uint;
		typedef int64_t				slong;		// note: keyword 'ulong' cannot be used with NV_ARM
	#else
		#include <unistd.h>
		#include <sys/time.h>
		#include <getopt.h>
		#include <cinttypes>
	#endif

	#include "cuda.h"		

	struct DataPtr {
		DataPtr();
	
		char		type;				// data type
		uint64_t	max;				// max element count
		uint64_t	num;				// curr element count
		uint64_t	size;				// size of data
		uint64_t	stride;				// stride of data	
		char*		cpu;				// cpu pointer		
		CUdeviceptr	gpu;				// gpu pointer
	};

	class Allocator {
	public:
		void Create		(DataPtr& p, int max, int stride );		
		void Destroy	(DataPtr& p);
		int  Alloc		(DataPtr& p);
		void Commit		(DataPtr& p);
		void Retrieve	(DataPtr& p);
		void Zero		(DataPtr& p);

		char* getCPU (DataPtr& p, int i)		{ return p.cpu + i*p.stride; }
	};

#endif

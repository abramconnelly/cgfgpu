
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>

#include "cgft.hpp"

// 32x   Tiles => TileGroup (64-bits)
// N(t)  Tiles per Path (variable, ave. 10k)
// 863x  Paths 

// TileGroup: 
// |<-- 32 bits -->|<-4bit-|-4bit-|-4bit-|-4bit->|  ..
//     canonical     kndx0   kndx1   ..    kndx7 
//
// TileBlocks: (index into tilebuf by tile path), SIZE= # of tilepaths
// | Block0: P0, cnt | Block1: P1, cnt | Block2: P2, cnt | .. 
//
// TileBuffer: (list of tile groups), SIZE= Total # tiles in genome / 32
// | Group1 | Group2 | Group3 | ..
// P0                         P1
//
// KnotMap:   {maps from a knot index to a knot table position}
// |-- kt0 --|-- kt1 --|-- kt2 --|
//    kndx0     kndx1     kndx2 
//
// KnotTable: {groups of knots|
// |- var0 ------|- var1 -|- var2 -|- var3 -|- var4 ------|- var5 ---| ..
// kt0                                 kt1                    kt2
//
// Variant:
// |<--- 32 bits -->|<--- 32 bits --->|<--- 32 bits --->|<--- 32 bits --->|
//        Avariant        Aspan             Bvariant           Bspan
//
// OverflowMap:  {maps from a tilepath to an overflow position}
// |-- ot0 --|-- ot1 --|-- ot2 --|    
//
// Overflow:
// |- var0 -|- var1 -|- var2 -|- var3 -| ..
// ot0               ot1               ot2


// Struct for memory management
struct data_ptr {
	int			num;
	uint64_t	stride;
	uint64_t	size;
	CUdeviceptr	gpu;
	char*		cpu;
};


// CGFT GPU
//
class cgft_gpu {

	void create ( cgf_t* src );

	void reserve ( data_ptr& ptr, int cnt, int stride );

	std::string	name;
									
	data_ptr	tileblks;		// tile blks  = index into tile buffer by tile path
	data_ptr	tilebuf;		// tile buf   = set of tile groups (TileGroup = Cache + KIndex)
	data_ptr    knotmap;		// knot map   = map from KIndex -> KnotPos {TileMap}
	data_ptr	knottable;		// knot table = list of knots by KnotPos
	data_ptr	overflowmap;	// overflow map
	data_ptr	overflow;		// overflow

};


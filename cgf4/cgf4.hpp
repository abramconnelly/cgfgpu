#ifndef GCGF_H
#define GCGF_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/time.h>

#include <getopt.h>

#include <cstdlib>
#include <map>
#include <vector>
#include <string>
#include <complex>
#include <iterator>

#include <cstdio>
#include <cinttypes>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>

#define CGF_MAGIC "{\"cgf.b\""
#define CGF_VERSION "0.4.0"
#define CGLF_VERSION "0.1.0"

#define OVF16_MAX 0xffff
#define OVF64_MAX 0xffffffffffffffff

#define SPAN_SDSL_ENC_VAL (1<<30)

extern char DEFAULT_TILEMAP[];

typedef struct tilepath4_type {
  std::string Name;

  uint64_t ExtraDataSize;
  std::vector< char > ExtraData;

  // Low quality information.
  // The "Het" and "Hom" portions indicate whether
  // the low quality data is same across both alleles
  // and doesn't indicate anything about whether the
  // variants are heterozygous or homozygous.
  // Hom represents the bulk of the data so
  // we get size savings by splitting it out
  // into Het and Hom parts.
  //
  // - Step is strictly increasing TileStep positions
  // - Variant are allele interleaved tile variants
  // - NocSum is the current inclusive count of the Noc elements
  // - NocStart is the start of the nocall run
  // - NocLen is the run of nocall elements
  //
  // Since the TileStep and NocSum are strictly non-decreasing,
  // they benefit from being an 'enc_vector' whereas the rest
  // benefit from being variable length encoded.
  //

  uint64_t LoqTileStepHomSize;
  uint64_t LoqTileVariantHomSize;
  uint64_t LoqTileNocSumHomSize;
  uint64_t LoqTileNocStartHomSize;
  uint64_t LoqTileNocLenHomSize;

  uint64_t LoqTileStepHetSize;
  uint64_t LoqTileVariantHetSize;
  uint64_t LoqTileNocSumHetSize;
  uint64_t LoqTileNocStartHetSize;
  uint64_t LoqTileNocLenHetSize;

  sdsl::enc_vector<> LoqTileStepHom;
  sdsl::vlc_vector<> LoqTileVariantHom;
  sdsl::enc_vector<> LoqTileNocSumHom;
  sdsl::vlc_vector<> LoqTileNocStartHom;
  sdsl::vlc_vector<> LoqTileNocLenHom;

  sdsl::enc_vector<> LoqTileStepHet;
  sdsl::vlc_vector<> LoqTileVariantHet;
  sdsl::enc_vector<> LoqTileNocSumHet;
  sdsl::vlc_vector<> LoqTileNocStartHet;
  sdsl::vlc_vector<> LoqTileNocLenHet;

} tilepath4_t;
typedef tilepath4_t tilepath_t;

typedef struct cgf4_type {

  unsigned char Magic[8];
  std::string CGFVersion;
  std::string LibraryVersion;
  uint64_t TilePathCount;
  std::string TileMap;

  // Number of elements in the Loq, Span, Canon and
  // CacheOverflow per 'block'.
  // CacheOverflow will be padded with Stride elements
  // per block even if there aren't that many cache overflows.
  // Loq, Span and Canon elements will be padded at the
  // end of tile path boundaries to fill the Stride.
  //
  uint32_t Stride;

  // Number of tiles proper in each tile path
  //
  std::vector< uint64_t > TileStepCount;

  // Offset in each vector taking into account padding at
  // end.
  // For example, Tile Path 0 might hold 7 tile steps
  // but the offset for the next tile path will occur in
  // position 1 for Loq, Span and Canon, and 4 for Cache Overflow.
  //
  std::vector< uint64_t > StrideOffset;

  // Loq quality vector. 1 represents a 'low quality' tile whereas '0'
  // represents a high quality tile.
  // LSB represents tile position p+0
  // MSB represents tile position p+7
  //
  std::vector< unsigned char > Loq;

  // Canon | Span
  // ------------
  //   0      1     -> anchor spanning tile
  //   1      1     -> non-anchor spanning tile
  //
  std::vector< unsigned char > Span;

  // 1 -> canonical tile (tile variant 0)
  // 0 -> non canoncial tile (overflow mustbe consulted for tile variant)
  // with exception as above
  //
  std::vector< unsigned char > Canon;

  // packed vector of hexist (4bit elements)
  // Stride indicates how many chars are associated with each of
  // the other vectors.
  // For now this should be 4 chars (8 hexits max) assocaited with each
  // Loq, Span and Canon entry.
  // The number of valid hexits is informed by the non canonical
  // bits set in Canon.
  // If the number of non-canonical elemetns in this stride (set of 8 tile steps, say)
  // the Overflow structures need to be consulted for tile variants
  //
  // The hexit encoding is as follows:
  // 0x0              complex
  // 0x1 to 0xe       lookup value in tile map (0x1 maps to tile map entry 1 etc.)
  // 0xf              overflow
  //
  std::vector< unsigned char > CacheOverflow;


  std::vector< unsigned uint64_t > OverflowOffset;
  std::vector< unsigned uint16_t > Overflow;

  std::vector< unsigned uint64_t > Overflow64Offset;
  std::vector< unsigned uint64_t > Overflow64;

  std::vector<uint64_t> TilePathStructOffset;
  std::vector<tilepath_t> TilePath;


} cgf4_t;

typedef cgf4_t cgf_t;

//--------------

typedef struct tilepath_vec_type {
  std::string name;
  std::vector<int> allele[2];
  std::vector<int> loq_flag[2];
  std::vector< std::vector<int> > loq_info[2];
} tilepath_vec_t;

typedef struct tilepath_ez_type {

  int tilepath;

  int N;

  // hiq
  //
  std::vector<uint64_t> cache;
  std::vector<unsigned char> span_bv;

  // interleaved overflow
  //
  std::vector<int16_t> ovf_vec;
  std::vector<int32_t> ovf32_vec;
  std::vector<int64_t> ovf64_vec;

  std::vector<char> data_vec;

  // loq
  //
  std::vector<unsigned char> loq_bv;  // floor( (N+7)/8 )

  int n_loq;

  std::vector<int> loq_info_pos;

  // interleaved for multi-allelic
  //
  std::vector<int> loq_info_variant;
  std::vector<int> loq_info_sn;
  std::vector<int> loq_info_noc;

  std::vector<int> loq_info_pos_hom;
  std::vector<int> loq_info_variant_hom;
  std::vector<int> loq_info_sn_hom;
  std::vector<int> loq_info_noc_hom;

  std::vector<int> loq_info_pos_het;
  std::vector<int> loq_info_variant_het;
  std::vector<int> loq_info_sn_het;
  std::vector<int> loq_info_noc_het;

} tilepath_ez_t;

typedef struct cgf_ez_type {

  std::string tilemap_str;

  std::map< std::string, int > tilemap;
  std::vector<tilepath_ez_t> tilepath;

} cgf_ez_t;


const char *read_tilemap_from_file(std::string &, const char *);

// io functions
//
//

void cgf_create_container(FILE *fp, const char *cgf_version, const char *cglf_version, const char *tilemap);
cgf_t *cgf_read(FILE *fp);
cgf_t *cgf_read_hiq(FILE *fp);
int cgf_read_band_tilepath(cgf_t *cgf, int idx, FILE *fp);

uint64_t cgf_write_to_file(cgf_t *cgf, const char *ofn);


void cgf_print(cgf_t *cgf);

void cgf_output_band_format(cgf_t *cgf, int tilepath_idx, FILE *fp, int hiq);

uint64_t cgf_write_to_file(cgf_t *cgf, const char *ofn);

// ez functions
//
//

void print_bgf(tilepath_vec_t &);
void print_tilepath_vec(tilepath_vec_t &);
void ez_create(tilepath_ez_t &, tilepath_vec_t &, std::map< std::string, int > &);
void ez_print(tilepath_ez_t &);

int load_tilemap(std::string &, std::map< std::string, int > &);
void mk_tilemap_key(std::string &key, tilepath_vec_t &tilepath, int tilestep, int n);


void ez_create_enc_vector(sdsl::enc_vector<> &enc_vec, std::vector<int> &v);
void ez_create_vlc_vector(sdsl::vlc_vector<> &vlc_vec, std::vector<int> &v);
void ez_to_tilepath(tilepath_t *tilepath, tilepath_ez_t *ez);



// extra functions
//

int cgf_sanity(cgf_t *cgf);

// concordance
//

int cgf_hiq_concordance(int *r_match, int *r_tot,
                        cgf_t *a, cgf_t *b,
                        int start_tile_path, int start_tile_step,
                        int end_tile_path, int end_tile_step);



// helper functions


inline int NumberOfSetBits32(uint32_t u)
{
  u = u - ((u >> 1) & 0x55555555);
  u = (u & 0x33333333) + ((u >> 2) & 0x33333333);
  return (((u + (u >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

// This is slower than the above but is more explicit
//
inline int NumberOfSetBits8(uint8_t u)
{
  u = (u & 0x55) + ((u>>1) & 0x55);
  u = (u & 0x33) + ((u>>2) & 0x33);
  u = (u & 0x0f) + ((u>>4) & 0x0f);
  return u;
}

#endif

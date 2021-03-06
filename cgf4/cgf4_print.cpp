#include "cgf4.hpp"

#ifdef WIN32
	#define PRIu64	"lld"
#endif

void cgf_print(cgf_t *cgf) {
  int i, j, n, prev_sum=0;
  size_t sz;
  int stride, tilepath;
  int block, prev_block;;
  tilepath_t *tp;

  printf("Magic: %c%c%c%c%c%c%c%c\n",
      cgf->Magic[0], cgf->Magic[1],
      cgf->Magic[2], cgf->Magic[3],
      cgf->Magic[4], cgf->Magic[5],
      cgf->Magic[6], cgf->Magic[7]);
  printf("CGFVersion: %s\n", cgf->CGFVersion.c_str());
  printf("LibVersion: %s\n", cgf->LibraryVersion.c_str());
  printf("Path: %" PRIu64 "\n", cgf->TilePathCount);

  sz = cgf->TileMap.size();
  if (sz<6) {
    printf("TileMap: ");
    for (i=0; i<sz; i++) { printf("%c", cgf->TileMap[i]); }
    printf("\n");
  } else {
    printf("TileMap: %c%c%c ... %c%c%c\n",
        cgf->TileMap[0], cgf->TileMap[1], cgf->TileMap[2],
        cgf->TileMap[sz-3], cgf->TileMap[sz-2], cgf->TileMap[sz-1]);
  }

  printf("Stride: %u\n",
      (unsigned int)cgf->Stride);

  printf("StepCount:");
  for (i=0; i<cgf->TileStepCount.size(); i++) {
    printf(" %" PRIu64, cgf->TileStepCount[i]);
  }
  printf("\n");

  printf("StrideOffset:");
  for (i=0; i<cgf->StrideOffset.size(); i++) {
    printf(" %" PRIu64, cgf->StrideOffset[i]);
  }
  printf("\n");

  printf("Loq | Span | Canon | Ovf\n");
  stride = 4;
  prev_block=0;
  for (tilepath=0; tilepath<cgf->StrideOffset.size(); tilepath++) {

    printf("tilepath[%i]:\n", tilepath);
    for (block=prev_block; block<cgf->StrideOffset[tilepath]; block++) {
      printf("[%08x] %02x %02x %02x %02x | %02x %02x %02x %02x | %02x %02x %02x %02x | %02x %02x %02x %02x\n",
          stride*(block - prev_block),
          cgf->Loq[stride*block+0], cgf->Loq[stride*block+1], cgf->Loq[stride*block+2], cgf->Loq[stride*block+3],
          cgf->Span[stride*block+0], cgf->Span[stride*block+1], cgf->Span[stride*block+2], cgf->Span[stride*block+3],
          cgf->Canon[stride*block+0], cgf->Canon[stride*block+1], cgf->Canon[stride*block+2], cgf->Canon[stride*block+3],
          cgf->CacheOverflow[stride*block+0], cgf->CacheOverflow[stride*block+1], cgf->CacheOverflow[stride*block+2], cgf->CacheOverflow[stride*block+3]);
    }
    prev_block = block;

    uint64_t idx64=0, end64=0;
    end64 = cgf->OverflowOffset[tilepath];
    if (tilepath>0) { idx64 = cgf->OverflowOffset[tilepath-1]; }

    printf("Overflow[%i]:", (int)(end64-idx64));
    for (; idx64<end64; idx64+=3) {
      printf(" {%i,%i,%i}",
          (int)cgf->Overflow[idx64+0],
          (int)cgf->Overflow[idx64+1],
          (int)cgf->Overflow[idx64+2]);
    }
    printf("\n");


    idx64=0;
    end64 = cgf->Overflow64Offset[tilepath];
    if (tilepath>0) { idx64 = cgf->Overflow64Offset[tilepath-1]; }

    printf("Overflow64[%i]:", (int)(end64-idx64));
    for (; idx64<end64; idx64+=3) {
      printf(" {%i,%i,%i}",
          (int)cgf->Overflow64[idx64+0],
          (int)cgf->Overflow64[idx64+1],
          (int)cgf->Overflow64[idx64+2]);
    }
    printf("\n");

    printf("TilePathStructOffset[%i]: %llu\n", tilepath, (unsigned long long int)cgf->TilePathStructOffset[tilepath]);

    tp = &(cgf->TilePath[tilepath]);
    printf("  Name: %s\n", tp->Name.c_str());
    printf("  ExtraDataSize: %llu\n", (unsigned long long int)tp->ExtraDataSize);
    printf("  ExtraData(%i):", (int)tp->ExtraData.size());
    for (j=0; j<tp->ExtraData.size(); j++) {
      if ((j%32)==0) { printf("\n   "); }
      printf(" %02x", (uint8_t)tp->ExtraData[j]);
    }
    printf("\n");

#ifdef USE_SDSL
    printf("  LoqHom_Size: %llu %llu %llu %llu %llu\n",
        (unsigned long long int)tp->LoqTileStepHomSize,
        (unsigned long long int)tp->LoqTileVariantHomSize,
        (unsigned long long int)tp->LoqTileNocSumHomSize,
        (unsigned long long int)tp->LoqTileNocStartHomSize,
        (unsigned long long int)tp->LoqTileNocLenHomSize);

    printf("  LoqHet_Size: %llu %llu %llu %llu %llu\n",
        (unsigned long long int)tp->LoqTileStepHetSize,
        (unsigned long long int)tp->LoqTileVariantHetSize,
        (unsigned long long int)tp->LoqTileNocSumHetSize,
        (unsigned long long int)tp->LoqTileNocStartHetSize,
        (unsigned long long int)tp->LoqTileNocLenHetSize);

    printf("...\n");

    printf("LoqHom:");
    n = (int)tp->LoqTileStepHom.size();
    printf(" ");
    prev_sum = 0;
    for (i=0; i<n; i++) {
      if ((i>0) && ((i%4)==0)) { printf("\n "); }
      printf(" {");
      printf("%i(%i,%i)[%i]:",
          (int)tp->LoqTileStepHom[i],
          (int)((tp->LoqTileVariantHom[2*i] == SPAN_SDSL_ENC_VAL) ? -1 : tp->LoqTileVariantHom[2*i]),
          (int)((tp->LoqTileVariantHom[2*i+1] == SPAN_SDSL_ENC_VAL) ? -1 : tp->LoqTileVariantHom[2*i+1]),
          //(int)tp->LoqTileVariantHom[2*i+1],
          (int)tp->LoqTileNocSumHom[i]);
      for (j=prev_sum; j<tp->LoqTileNocSumHom[i]; j++) {
        if (j>prev_sum) { printf(" "); }
        printf("[%i+%i]",
            (int)tp->LoqTileNocStartHom[j],
            (int)tp->LoqTileNocLenHom[j]);
      }
      prev_sum = (int)tp->LoqTileNocSumHom[i];
      printf("}");
    }
    printf("\n");

    printf("LoqHet:\n");
    n = (int)tp->LoqTileStepHet.size();
    printf(" ");
    prev_sum = 0;
    for (i=0; i<n; i++) {
      if ((i>0) && ((i%4)==0)) { printf("\n "); }
      printf(" {");
      printf("%i(%i,%i)[%i,%i]:",
          (int)tp->LoqTileStepHet[i],
          (int)tp->LoqTileVariantHet[2*i],
          (int)tp->LoqTileVariantHet[2*i+1],
          (int)tp->LoqTileNocSumHet[2*i],
          (int)tp->LoqTileNocSumHet[2*i+1]);

      for (j=prev_sum; j<tp->LoqTileNocSumHet[2*i]; j++) {
        if (j>prev_sum) { printf(" "); }
        printf("[%i+%i]",
            (int)tp->LoqTileNocStartHet[j],
            (int)tp->LoqTileNocLenHet[j]);
      }
      prev_sum = (int)tp->LoqTileNocSumHet[2*i];

      for (j=prev_sum; j<tp->LoqTileNocSumHet[2*i+1]; j++) {
        if (j>prev_sum) { printf(" "); }

        if (j >= (int)tp->LoqTileNocStartHet.size()) {
          printf("SANITY: j (%i) >= LoqTileNocStartHet (%i)\n", j, (int)tp->LoqTileNocStartHet.size());
        }

        if (j >= (int)tp->LoqTileNocLenHet.size()) {
          printf("SANITY: j (%i) >= LoqTileNocLenHet (%i)\n", j, (int)tp->LoqTileNocLenHet.size());
        }

        printf("[%i+%i]",
            (int)tp->LoqTileNocStartHet[j],
            (int)tp->LoqTileNocLenHet[j]);
      }
      prev_sum = (int)tp->LoqTileNocSumHet[2*i+1];

      printf("}");
    }
    printf("\n");
#endif



  }

  printf("\n");
  printf("\n");


}


void mk_vec_tilemap(std::vector< std::vector< std::vector<int> > > &vtm, const char *tm) {
  int i, j, k, ii, jj;
  char *chp;
  std::string s;

  std::vector< std::vector<int> > tm_entry;
  std::vector< std::vector<int> > x;
  std::vector< int > y;

  int entry_count=0;
  int tval;

  int enc_val=0;

  std::string s_tm = tm;

  std::stringstream line_stream(s_tm);
  std::vector<std::string> lines, alleles, variants, vals;
  std::string item, item1, item2, item0;

  while (std::getline(line_stream, item, '\n')) {
    lines.push_back(item);
    entry_count++;
    if (entry_count==16) { break; }
  }

  vtm.clear();
  for (i=0; i<lines.size(); i++) {
    std::stringstream ss(lines[i]);

    tm_entry.clear();

    alleles.clear();
    while (std::getline(ss, item, ':')) {
      alleles.push_back(item);
    }

    for (j=0; j<alleles.size(); j++) {

      y.clear();
      std::stringstream s0(alleles[j]);
      while (std::getline(s0, item0, ';')) {

        int eo=0;
        std::stringstream s1(item0);
        while (std::getline(s1, item1, '+')) {

          tval = (int)strtol(item1.c_str(), NULL, 16);
          if (eo==1) {
            for (ii=1; ii<tval; ii++) { y.push_back(-1); }
          } else {
            y.push_back(tval);
          }

          eo = 1-eo;
        }
      }
      tm_entry.push_back(y);
    }
    vtm.push_back(tm_entry);

    enc_val++;
  }

}



//-----
//
int cgf_output_band_format(cgf_t *cgf, int tilepath_idx, FILE *fp, int hiq) {
  int i, j, k;
  int ii, jj;
  int ntile, n_q, n_r, n_ovf, tilestep=0;

  tilepath_t *tilepath;
  int loc_verbose = 0;
  int ovf_count=0;

  uint64_t s64, i64, n64;

  std::vector<int> variant_v[2];
  std::vector< std::vector<int> > noc_v[2];
  std::vector<int> v;

  int output_noc=1;
  unsigned char *loq, *span, *canon, *cache_ovf;

  int block, cpos;
  unsigned char loq8, hiq8, canon8, span8, anchor8;

  uint32_t loq_mask, span_mask, hiq_mask, cache_mask, xspan_mask, anchor_mask, cache_ovf_mask, canon_mask;
  uint32_t lo_cache;
  int start_offset=0;

  int hexit[8], cur_hexit, n_hexit, n_cache_ovf;
  int stride;

  std::vector< std::vector< std::vector<int> > > tilemap_vec;


  //--

  if (tilepath_idx >= (int)cgf->TileStepCount.size()) { return -1; }

  mk_vec_tilemap(tilemap_vec, cgf->TileMap.c_str());

  for (i=0; i<8; i++) { hexit[8] = 0; }

  ntile = (int)cgf->TileStepCount[tilepath_idx];
  stride = (int)cgf->Stride;

  if (tilepath_idx>0) { start_offset = cgf->StrideOffset[tilepath_idx-1]; }
  loq = &(cgf->Loq[start_offset*stride]);
  span = &(cgf->Span[start_offset*stride]);
  canon = &(cgf->Canon[start_offset*stride]);
  cache_ovf = &(cgf->CacheOverflow[start_offset*stride]);

  tilepath = &(cgf->TilePath[tilepath_idx]);

  if (hiq>0) {
    output_noc=0;
  }

  variant_v[0].resize(ntile);
  variant_v[1].resize(ntile);
  noc_v[0].resize(ntile);
  noc_v[1].resize(ntile);

  for (i=0; i<ntile; i++) {
    variant_v[0][i] = -1;
    variant_v[1][i] = -1;

    noc_v[0][i] = v;
    noc_v[1][i] = v;
  }


  for (i=0; i<ntile; i++) {
    n_q = i/8;
    n_r = i%8;

    hiq8 = ~(loq[n_q]);

    if ( (hiq8&(1<<n_r)) && (canon[n_q]&(1<<n_r)) ) {
      variant_v[0][i] = 0;
      variant_v[1][i] = 0;
      continue;
    }

  }

  n_q = (ntile + ((8*stride)-1))/(stride*8);

  for (ii=0; ii<n_q; ii++) {

    ovf_count=0;
    loq_mask = loq[4*ii] | (loq[4*ii+1]<<8) | (loq[4*ii+2]<<16) | (loq[4*ii+3]<<24);
    span_mask = span[4*ii] | (span[4*ii+1]<<8) | (span[4*ii+2]<<16) | (span[4*ii+3]<<24);
    cache_mask = canon[4*ii] | (canon[4*ii+1]<<8) | (canon[4*ii+2]<<16) | (canon[4*ii+3]<<24);
    lo_cache = cache_ovf[4*ii] | (cache_ovf[4*ii+1]<<8) | (cache_ovf[4*ii+2]<<16) | (cache_ovf[4*ii+3]<<24);

    xspan_mask = ~span_mask;
    hiq_mask = ~loq_mask;
    canon_mask = cache_mask & xspan_mask & hiq_mask;
    anchor_mask = span_mask & hiq_mask & (~cache_mask);
    cache_ovf_mask = (anchor_mask & hiq_mask) | ((~span_mask) & (~canon_mask) & hiq_mask);

    for (i=0; i<8; i++) {
      hexit[i] = (int)((lo_cache & (0xf<<(4*i))) >> (4*i));
    }

    cur_hexit=0;
    n_hexit=0;
    n_cache_ovf = 0;
    for (i=0; i<32; i++) {
      if (cache_ovf_mask & (1<<i)) {

        if ((cur_hexit<8) && (hexit[cur_hexit] != 0xf)) {

          int hexit_val = hexit[cur_hexit];

          for (j=0; j<tilemap_vec[hexit_val][0].size(); j++) {
            int cur_tilestep = 32*ii + i + j;

            if (cur_tilestep >= ntile) { continue; }

            if (loc_verbose) {
              fprintf(fp, "    cur_tilestep %i : hexit_val %i -> %i %i\n",
                  cur_tilestep,
                  hexit_val,
                  tilemap_vec[hexit_val][0][j],
                  tilemap_vec[hexit_val][1][j] );
            }

            variant_v[0][cur_tilestep] = tilemap_vec[hexit_val][0][j];
            variant_v[1][cur_tilestep] = tilemap_vec[hexit_val][1][j];
          }
        }
        cur_hexit++;

        n_cache_ovf++;
        n_hexit++;
      }
    }
    if (n_hexit>8) { n_hexit=8; }


    for (i=0; i<32; i++) {

      if (tilestep >= ntile) { continue; }

      if (canon_mask&(1<<i)) {
        variant_v[0][tilestep] = 0;
        variant_v[1][tilestep] = 0;

        if (loc_verbose) { fprintf(fp, "  tilestep %i -> 0 0\n", tilestep); }

      }
      else if (cache_ovf_mask & (1<<i)) {
        if (cur_hexit<n_hexit) {

          if (loc_verbose) {
            fprintf(fp, ">>> tilestep %i hexit %i\n", tilestep, hexit[cur_hexit]);
          }

        }
        cur_hexit++;
      }
      //else if (loq_mask & (1<<i)) { if (loc_verbose) { fprintf(fp, "  loq %i\n", tilestep); } }

      tilestep++;
    }



  }


  s64 = 0;
  if (tilepath_idx>0) { s64 = cgf->OverflowOffset[tilepath_idx-1]; }
  for (i64=s64; i64<cgf->OverflowOffset[tilepath_idx]; i64+=3) {

    tilestep = (int)cgf->Overflow[i64];
    int vara = (int)cgf->Overflow[i64+1];
    int varb = (int)cgf->Overflow[i64+2];

    if (vara >= OVF16_MAX) { vara = -1; }
    if (varb >= OVF16_MAX) { varb = -1; }

    variant_v[0][tilestep] = vara;
    variant_v[1][tilestep] = varb;
  }

  s64 = 0;
  if (tilepath_idx>0) { s64 = cgf->Overflow64Offset[tilepath_idx-1]; }
  for (i64=s64; i64<cgf->Overflow64Offset[tilepath_idx]; i64+=3) {

    tilestep = (int)cgf->Overflow64[i64];
    int vara = (int)cgf->Overflow64[i64+1];
    int varb = (int)cgf->Overflow64[i64+2];

    if (vara >= OVF64_MAX) { vara = -1; }
    if (varb >= OVF64_MAX) { varb = -1; }

    variant_v[0][tilestep] = vara;
    variant_v[1][tilestep] = varb;
  }

#ifdef USE_SDSL

  if (output_noc) {

    int prev_noc_start = 0;
    for (i=0; i<tilepath->LoqTileStepHom.size(); i++) {
      tilestep = tilepath->LoqTileStepHom[i];
      int vara = tilepath->LoqTileVariantHom[2*i];
      int varb = tilepath->LoqTileVariantHom[2*i+1];

      if (vara==SPAN_SDSL_ENC_VAL) { vara = -1; }
      if (varb==SPAN_SDSL_ENC_VAL) { varb = -1; }

      variant_v[0][tilestep] = vara;
      variant_v[1][tilestep] = varb;

      if (loc_verbose) {
        fprintf(fp, "  loq hom %i -> %i %i\n", tilestep, vara, varb);
      }

      for (j=prev_noc_start; j<tilepath->LoqTileNocSumHom[i]; j++) {
        noc_v[0][tilestep].push_back( (int)tilepath->LoqTileNocStartHom[j] );
        noc_v[0][tilestep].push_back( (int)tilepath->LoqTileNocLenHom[j] );

        noc_v[1][tilestep].push_back( (int)tilepath->LoqTileNocStartHom[j] );
        noc_v[1][tilestep].push_back( (int)tilepath->LoqTileNocLenHom[j] );
      }
      prev_noc_start = (int)tilepath->LoqTileNocSumHom[i];
    }

    prev_noc_start=0;
    for (i=0; i<tilepath->LoqTileStepHet.size(); i++) {
      tilestep = tilepath->LoqTileStepHet[i];
      int vara = tilepath->LoqTileVariantHet[2*i];
      int varb = tilepath->LoqTileVariantHet[2*i+1];

      if (vara==SPAN_SDSL_ENC_VAL) { vara = -1; }
      if (varb==SPAN_SDSL_ENC_VAL) { varb = -1; }

      variant_v[0][tilestep] = vara;
      variant_v[1][tilestep] = varb;

      if (loc_verbose) {
        fprintf(fp, "  loq het %i -> %i %i\n", tilestep, vara, varb);
      }

      for (j=prev_noc_start; j<tilepath->LoqTileNocSumHet[2*i]; j++) {
        noc_v[0][tilestep].push_back( (int)tilepath->LoqTileNocStartHet[j] );
        noc_v[0][tilestep].push_back( (int)tilepath->LoqTileNocLenHet[j] );
      }
      prev_noc_start = (int)tilepath->LoqTileNocSumHet[2*i];

      for (j=prev_noc_start; j<tilepath->LoqTileNocSumHet[2*i+1]; j++) {
        noc_v[1][tilestep].push_back( (int)tilepath->LoqTileNocStartHet[j] );
        noc_v[1][tilestep].push_back( (int)tilepath->LoqTileNocLenHet[j] );
      }
      prev_noc_start = (int)tilepath->LoqTileNocSumHet[2*i+1];

    }
  }
#endif

  // output all vectors
  //

  for (i=0; i<2; i++) {
    fprintf(fp, "[");
    for (j=0; j<variant_v[i].size(); j++) {
      fprintf(fp, " %i", variant_v[i][j]);
    }
    fprintf(fp, "]\n");
  }

  if (output_noc) {
    for (i=0; i<2; i++) {
      fprintf(fp, "[");
      for (j=0; j<noc_v[i].size(); j++) {
        fprintf(fp, "[");
        for (k=0; k<noc_v[i][j].size(); k++) {
          fprintf(fp, " %i", noc_v[i][j][k]);
        }
        fprintf(fp, " ]");
      }
      fprintf(fp, "]\n");
    }
  }

  return 0;
}

//----

#define _ZCHUNK 16384



// takenf rom http://zlib.net/zpipe.c (public domain by Mark Adler)
//
int _infz(std::vector< uint8_t > &src, std::vector< uint8_t > &dst) {

#ifdef USE_ZSTREAM
  int i, ret;
  unsigned int have;
  z_stream strm;
  unsigned char in[_ZCHUNK];
  unsigned char out[_ZCHUNK];

  dst.clear();

  // allocate inflate state
  //
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.next_in = Z_NULL;
  ret = inflateInit(&strm);
  if (ret != Z_OK) { return ret; }

  strm.avail_in = src.size();
  if (strm.avail_in == 0) { return 0; }

  strm.next_in = (unsigned char *)(&(src[0]));


  // decompress until deflate stream ends or end of file
  //
  do {

    if (strm.avail_in==0) { break; }

    // run inflate() on input until output buffer not full
    //
    do {
      strm.avail_out = _ZCHUNK;
      strm.next_out = out;
      ret = inflate(&strm, Z_NO_FLUSH);
      switch (ret) {
        case Z_STREAM_ERROR:
        case Z_NEED_DICT: ret = Z_DATA_ERROR;
        case Z_DATA_ERROR:
        case Z_MEM_ERROR:
          (void)inflateEnd(&strm);
          return ret;
      }
      have = _ZCHUNK - strm.avail_out;
      for (i=0; i<have; i++) { dst.push_back(out[i]); }
    } while (strm.avail_out == 0);

  } while (ret != Z_STREAM_END);

  (void)inflateEnd(&strm);
  return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;

#else
	return 0;	// fail, no zstream
#endif
}

//-----
//
int cgf_output_band_format2(cgf_t *cgf, int tilepath_idx, FILE *fp, int step_start, int step_n, uint32_t fill_level, int default_fill) {
  int i, j, k;
  int ii, jj;
  int ntile, n_q, n_r, n_ovf, tilestep=0;

  int idx, byte_pos, a;

  tilepath_t *tilepath;
  int loc_verbose = 0;
  int ovf_count=0;

  uint32_t n32, i32, u32, step32, pos32;
  uint64_t s64, i64, n64;

  size_t end_noninc;

  std::vector<int> variant_v[2];
  std::vector< std::vector<int> > noc_v[2];
  std::vector<int> v;

  unsigned char *loq, *span, *canon, *cache_ovf;

  int block, cpos;
  unsigned char loq8, hiq8, canon8, span8, anchor8;

  uint32_t loq_mask, span_mask, hiq_mask, cache_mask, xspan_mask, anchor_mask, cache_ovf_mask, canon_mask;
  uint32_t lo_cache;
  int start_offset=0;

  int hexit[8], cur_hexit, n_hexit, n_cache_ovf;
  int stride;

  std::vector< std::vector< std::vector<int> > > tilemap_vec;

  std::vector< uint8_t > _tbuf, _encbuf;
  std::string sbuf;


  //--

  mk_vec_tilemap(tilemap_vec, cgf->TileMap.c_str());

  for (i=0; i<8; i++) { hexit[i] = 0; }

  if (tilepath_idx >= (int)cgf->TileStepCount.size()) { return -1; }

  ntile = (int)cgf->TileStepCount[tilepath_idx];
  stride = (int)cgf->Stride;

  if (tilepath_idx>0) { start_offset = cgf->StrideOffset[tilepath_idx-1]; }
  loq = &(cgf->Loq[start_offset*stride]);
  span = &(cgf->Span[start_offset*stride]);
  canon = &(cgf->Canon[start_offset*stride]);
  cache_ovf = &(cgf->CacheOverflow[start_offset*stride]);

  tilepath = &(cgf->TilePath[tilepath_idx]);

  variant_v[0].resize(ntile);
  variant_v[1].resize(ntile);
  noc_v[0].resize(ntile);
  noc_v[1].resize(ntile);


  for (i=0; i<ntile; i++) {
    variant_v[0][i] = default_fill;
    variant_v[1][i] = default_fill;
    noc_v[0][i] = v;
    noc_v[1][i] = v;
  }

  for (i=0; i<ntile; i++) {
    n_q = i/8;
    n_r = i%8;

    hiq8 = ~(loq[n_q]);

    if ( (hiq8&(1<<n_r)) && (canon[n_q]&(1<<n_r)) ) {

      if (fill_level & 1) {
        variant_v[0][i] = 0;
        variant_v[1][i] = 0;
      }
      continue;
    }

  }

  n_q = (ntile + ((8*stride)-1))/(stride*8);

  for (ii=0; ii<n_q; ii++) {

    ovf_count=0;
    loq_mask = loq[4*ii] | (loq[4*ii+1]<<8) | (loq[4*ii+2]<<16) | (loq[4*ii+3]<<24);
    span_mask = span[4*ii] | (span[4*ii+1]<<8) | (span[4*ii+2]<<16) | (span[4*ii+3]<<24);
    cache_mask = canon[4*ii] | (canon[4*ii+1]<<8) | (canon[4*ii+2]<<16) | (canon[4*ii+3]<<24);
    lo_cache = cache_ovf[4*ii] | (cache_ovf[4*ii+1]<<8) | (cache_ovf[4*ii+2]<<16) | (cache_ovf[4*ii+3]<<24);

     xspan_mask = ~span_mask;
    hiq_mask = ~loq_mask;
    canon_mask = cache_mask & xspan_mask & hiq_mask;
    anchor_mask = span_mask & hiq_mask & (~cache_mask);
    cache_ovf_mask = (anchor_mask & hiq_mask) | ((~span_mask) & (~canon_mask) & hiq_mask);

    for (i=0; i<8; i++) {
      hexit[i] = (int)((lo_cache & (0xf<<(4*i))) >> (4*i));
    }

    cur_hexit=0;
    n_hexit=0;
    n_cache_ovf = 0;
    for (i=0; i<32; i++) {
      if (cache_ovf_mask & (1<<i)) {

        if ((cur_hexit<8) && (hexit[cur_hexit] != 0xf)) {

          int hexit_val = hexit[cur_hexit];

          for (j=0; j<tilemap_vec[hexit_val][0].size(); j++) {
            int cur_tilestep = 32*ii + i + j;

            if (cur_tilestep >= ntile) { continue; }

            if (loc_verbose) {
              fprintf(fp, "    cur_tilestep %i : hexit_val %i -> %i %i\n",
                  cur_tilestep,
                  hexit_val,
                  tilemap_vec[hexit_val][0][j],
                  tilemap_vec[hexit_val][1][j] );
            }

            if (fill_level & 2) {
              variant_v[0][cur_tilestep] = tilemap_vec[hexit_val][0][j];
              variant_v[1][cur_tilestep] = tilemap_vec[hexit_val][1][j];
            }
          }
        }
        cur_hexit++;

        n_cache_ovf++;
        n_hexit++;
      }
    }
    if (n_hexit>8) { n_hexit=8; }


    for (i=0; i<32; i++) {

      if (tilestep >= ntile) { continue; }

      if (canon_mask&(1<<i)) {

        if (fill_level & 1) {
          variant_v[0][tilestep] = 0;
          variant_v[1][tilestep] = 0;
        }

        if (loc_verbose) { fprintf(fp, "  tilestep %i -> 0 0\n", tilestep); }

      }
      else if (cache_ovf_mask & (1<<i)) {
        if (cur_hexit<n_hexit) {

          if (loc_verbose) {
            fprintf(fp, ">>> tilestep %i hexit %i\n", tilestep, hexit[cur_hexit]);
          }

        }
        cur_hexit++;
      }

      tilestep++;
    }


  }


  s64 = 0;
  if (tilepath_idx>0) { s64 = cgf->OverflowOffset[tilepath_idx-1]; }
  for (i64=s64; i64<cgf->OverflowOffset[tilepath_idx]; i64+=3) {

    tilestep = (int)cgf->Overflow[i64];
    int vara = (int)cgf->Overflow[i64+1];
    int varb = (int)cgf->Overflow[i64+2];

    if (vara >= OVF16_MAX) { vara = -1; }
    if (varb >= OVF16_MAX) { varb = -1; }

    if (fill_level & 4) {
      variant_v[0][tilestep] = vara;
      variant_v[1][tilestep] = varb;
    }
  }

  s64 = 0;
  if (tilepath_idx>0) { s64 = cgf->Overflow64Offset[tilepath_idx-1]; }
  for (i64=s64; i64<cgf->Overflow64Offset[tilepath_idx]; i64+=3) {

    tilestep = (int)cgf->Overflow64[i64];
    int vara = (int)cgf->Overflow64[i64+1];
    int varb = (int)cgf->Overflow64[i64+2];

    if (vara >= OVF64_MAX) { vara = -1; }
    if (varb >= OVF64_MAX) { varb = -1; }


    if (fill_level & 4) {
      variant_v[0][tilestep] = vara;
      variant_v[1][tilestep] = varb;
    }
  }

#ifdef USE_SDSL
  if (fill_level & 8) {

    int prev_noc_start = 0;
    for (i=0; i<tilepath->LoqTileStepHom.size(); i++) {
      tilestep = tilepath->LoqTileStepHom[i];
      int vara = tilepath->LoqTileVariantHom[2*i];
      int varb = tilepath->LoqTileVariantHom[2*i+1];

      if (vara==SPAN_SDSL_ENC_VAL) { vara = -1; }
      if (varb==SPAN_SDSL_ENC_VAL) { varb = -1; }

      variant_v[0][tilestep] = vara;
      variant_v[1][tilestep] = varb;

      if (loc_verbose) {
        fprintf(fp, "  loq hom %i -> %i %i\n", tilestep, vara, varb);
      }

      for (j=prev_noc_start; j<tilepath->LoqTileNocSumHom[i]; j++) {
        noc_v[0][tilestep].push_back( (int)tilepath->LoqTileNocStartHom[j] );
        noc_v[0][tilestep].push_back( (int)tilepath->LoqTileNocLenHom[j] );

        noc_v[1][tilestep].push_back( (int)tilepath->LoqTileNocStartHom[j] );
        noc_v[1][tilestep].push_back( (int)tilepath->LoqTileNocLenHom[j] );
      }
      prev_noc_start = (int)tilepath->LoqTileNocSumHom[i];
    }

    prev_noc_start=0;
    for (i=0; i<tilepath->LoqTileStepHet.size(); i++) {
      tilestep = tilepath->LoqTileStepHet[i];
      int vara = tilepath->LoqTileVariantHet[2*i];
      int varb = tilepath->LoqTileVariantHet[2*i+1];

      if (vara==SPAN_SDSL_ENC_VAL) { vara = -1; }
      if (varb==SPAN_SDSL_ENC_VAL) { varb = -1; }

      variant_v[0][tilestep] = vara;
      variant_v[1][tilestep] = varb;

      if (loc_verbose) {
        fprintf(fp, "  loq het %i -> %i %i\n", tilestep, vara, varb);
      }

      for (j=prev_noc_start; j<tilepath->LoqTileNocSumHet[2*i]; j++) {
        noc_v[0][tilestep].push_back( (int)tilepath->LoqTileNocStartHet[j] );
        noc_v[0][tilestep].push_back( (int)tilepath->LoqTileNocLenHet[j] );
      }
      prev_noc_start = (int)tilepath->LoqTileNocSumHet[2*i];

      for (j=prev_noc_start; j<tilepath->LoqTileNocSumHet[2*i+1]; j++) {
        noc_v[1][tilestep].push_back( (int)tilepath->LoqTileNocStartHet[j] );
        noc_v[1][tilestep].push_back( (int)tilepath->LoqTileNocLenHet[j] );
      }
      prev_noc_start = (int)tilepath->LoqTileNocSumHet[2*i+1];

    }
  }
#endif

  // output all vectors
  //

  end_noninc = variant_v[0].size();

  if (step_n > 0) {
    if ( (size_t)(step_start + step_n) < end_noninc ) {
      end_noninc = (size_t)(step_start + step_n);
    }
  }

  for (i=0; i<2; i++) {
    fprintf(fp, "[");
    for (j=step_start; j<end_noninc; j++) {
      fprintf(fp, " %i", variant_v[i][j]);
    }
    fprintf(fp, "]\n");
  }

  // genotype filel
  //
  if (fill_level & 16) {

    byte_pos = 0;

    while (byte_pos < tilepath->ExtraData.size()) {

      if ((byte_pos + 4) > tilepath->ExtraData.size()) { break; }

      sbuf.clear();
      for (i=0; i<4; i++) {
        if (tilepath->ExtraData[byte_pos]>0) {
          //buf += tilepath->ExtraData[byte_pos];
          sbuf.push_back(tilepath->ExtraData[byte_pos]);
        }
        byte_pos++;
      }

      if ((sbuf=="gt0.") || (sbuf=="gt1.") ||
          (sbuf=="gtz0") || (sbuf=="gtz1")) {

        int gtz_flag = 0;
        if ((sbuf=="gtz0") || (sbuf=="gtz1")) { gtz_flag = 1; }

        if      ((sbuf=="gt0.") || (sbuf=="gtz0")) { a = 0; }
        else if ((sbuf=="gt1.") || (sbuf=="gtz1")) { a = 1; }
        else { break; }

        if ((byte_pos + 4) > tilepath->ExtraData.size()) { break; }
        n32 = (uint32_t)( *((uint32_t *)(&(tilepath->ExtraData[byte_pos]))) );
        byte_pos += 4;

        if ((byte_pos + n32) > tilepath->ExtraData.size()) { break; }

        _tbuf.clear();
        _encbuf.clear();
        if (gtz_flag) {

          for (i=0; i<n32; i++) {
            _tbuf.push_back(tilepath->ExtraData[byte_pos+i]);
          }
          _infz(_tbuf, _encbuf);
        }
        else {
          for (i=0; i<n32; i++) {
            _encbuf.push_back(tilepath->ExtraData[byte_pos+i]);
          }
        }
        byte_pos += n32;

        if ((_encbuf.size()%8)!=0) { break; }

        for (i=0; i<_encbuf.size(); i+=8) {

          step32 = (uint32_t)( *((uint32_t *)(&(_encbuf[i]))) );
          pos32 = (uint32_t)( *((uint32_t *)(&(_encbuf[i+4]))) );

          if (step32 < noc_v[a].size()) {
            noc_v[a][step32].push_back(pos32);
            noc_v[a][step32].push_back(1);
          }
          else { break; }

        }
      }

    }

    if (byte_pos != tilepath->ExtraData.size()) {
      //error
    }

  }

  if (fill_level & 8) {
    for (i=0; i<2; i++) {
      fprintf(fp, "[");
      for (j=step_start; j<end_noninc; j++) {
        fprintf(fp, "[");
        for (k=0; k<noc_v[i][j].size(); k++) {
          fprintf(fp, " %i", noc_v[i][j][k]);
        }
        fprintf(fp, " ]");
      }
      fprintf(fp, "]\n");
    }
  }

  return 0;
}

//---


void cgf4_print_tilepath_stats(cgf_t *cgf, cgf_opt_t *cgf_opt) {
  int i, j, k, ii, jj;
  int tilepath;

  int match=0, tot=0;

  int start_pos, end_pos, n_pos,
      n_q, n_q_end;

  unsigned char *loq,
                *span,
                *canon,
                *cache_ovf;
  uint16_t *overflow;
  uint64_t *overflow64;

  int hexit[8],
      hexit_relative_step[8];

  uint32_t loq_mask,
           hiq_mask,
           span_mask,
           xspan_mask,
           cache_mask,
           lo_cache,
           canon_mask,
           anchor_mask,
           nonnchor_span_mask,
           cache_ovf_mask;

  uint32_t env_mask;

  uint32_t u32, t_u32;
  uint16_t prev_tile_step16;
  int prev_tile_step;

  int tile_step_block_start;
  int tile_offset;

  int stride;
  int tilepath_idx;

  int start_block_tile = 0;
  int idx, z;

  int t_match, t_tot;

  int tot_hiq, tot_hiq_knot, tot_hiq_anchor, tot_span;
  int n_tilepath=1;

  uint64_t iistart;
  uint64_t jjstart;

  if (cgf_opt->verbose) {
    printf("#tilepath,n_tot,n_hiq,n_hiq_knot,n_hiq_anchor\n");
  }

  loq         = &(cgf->Loq[0]);
  span         = &(cgf->Span[0]);
  canon       = &(cgf->Canon[0]);
  cache_ovf   = &(cgf->CacheOverflow[0]);
  overflow     = &(cgf->Overflow[0]);
  overflow64   = &(cgf->Overflow64[0]);


  stride = cgf->Stride;

  if (cgf_opt->endtilepath>=0) {
    n_tilepath = cgf_opt->endtilepath - cgf_opt->tilepath + 1;
  }

  for (tilepath = cgf_opt->tilepath; tilepath < (cgf_opt->tilepath + n_tilepath); tilepath++) {

    start_pos = 0;
    if (tilepath>0) {
      start_pos = 8 * stride * (cgf->StrideOffset[tilepath-1]);
    }

    end_pos = 0;
    if (tilepath>0) {
      end_pos = 8 * stride * (cgf->StrideOffset[tilepath-1]);
    }
    end_pos += cgf->TileStepCount[tilepath];

    n_q = start_pos / (8*stride);
    n_q_end = end_pos / (8*stride);

    tot_hiq = 0;
    tot_hiq_knot = 0;
    tot_hiq_anchor = 0;
    tot_span=0;

    for (ii=n_q; ii<=n_q_end; ii++) {

      // collect the uint32_t bit vectors into a convenient form
      //
      loq_mask    = loq[4*ii] | (loq[4*ii+1]<<8) | (loq[4*ii+2]<<16) | (loq[4*ii+3]<<24);
      span_mask   = span[4*ii] | (span[4*ii+1]<<8) | (span[4*ii+2]<<16) | (span[4*ii+3]<<24);
      cache_mask  = canon[4*ii] | (canon[4*ii+1]<<8) | (canon[4*ii+2]<<16) | (canon[4*ii+3]<<24);
      lo_cache    = cache_ovf[4*ii] | (cache_ovf[4*ii+1]<<8) | (cache_ovf[4*ii+2]<<16) | (cache_ovf[4*ii+3]<<24);

      xspan_mask      = ~span_mask;
      hiq_mask        = ~loq_mask;

      // non anchor spanning tiles are indicated with a span bit set and a canon bit set
      // so make sure to account for them to get the actual canononical bits out.
      //
      canon_mask      = cache_mask & xspan_mask & hiq_mask;

      // anchor tile bit vector for convenience.
      //
      anchor_mask     = hiq_mask & (~cache_mask) & span_mask;


      tot_hiq += NumberOfSetBits32(hiq_mask);
      tot_hiq_knot += NumberOfSetBits32(anchor_mask) + NumberOfSetBits32(hiq_mask & xspan_mask);
      tot_hiq_anchor += NumberOfSetBits32(anchor_mask);

      tot_span += NumberOfSetBits32(span_mask);

    }

    printf("%i,%i,%i,%i,%i\n", tilepath, (int)cgf->TileStepCount[tilepath], tot_hiq, tot_hiq_knot, tot_hiq_anchor);



  }

}

//---

void cgf4_print_header_json(cgf_t *cgf, FILE *ofp) {
  int i, j, prv;

  if (cgf->TileMapCacheInit==0) {
    str2tilemap(cgf->TileMap, &(cgf->TileMapCache));
    cgf->TileMapCacheInit=1;
  }

  fprintf(ofp, "{\n");
  fprintf(ofp, "  \"Magic\":\"");
  for (i=0; i<6; i++) {
    if (cgf->Magic[i] == '"') {
      fprintf(ofp, "\\");
    }
    fprintf(ofp, "%c", cgf->Magic[i]);
  }
  fprintf(ofp, "\",\n");

  fprintf(ofp, "  \"Version\":\"%s\",\n", cgf->CGFVersion.c_str());
  fprintf(ofp, "  \"LibraryVersion\":\"%s\",\n", cgf->LibraryVersion.c_str());
  fprintf(ofp, "  \"TilePathCount\":%llu,\n", (long long unsigned int)cgf->TilePathCount);
  fprintf(ofp, "  \"Stride\":%llu,\n", (long long unsigned int)cgf->Stride);
  fprintf(ofp, "  \"TileMap\":[");

  prv=0;
  for (i=0; i<cgf->TileMapCache.offset.size(); i++) {
    if ((i%8)==0) {
      if (i>0) { fprintf(ofp, ","); }
      fprintf(ofp, "\n    ");
    }
    else {
      fprintf(ofp, ", ");
    }

    fprintf(ofp, "[[");
    for (j=prv; j<cgf->TileMapCache.offset[i]; j++) {
      if (j>prv) { fprintf(ofp, ","); }
	  if (j<cgf->TileMapCache.variant[0].size()) fprintf(ofp, "%i", cgf->TileMapCache.variant[0][j]);
    }
    fprintf(ofp, "],[");
    for (j=prv; j<cgf->TileMapCache.offset[i]; j++) {
      if (j>prv) { fprintf(ofp, ","); }
      if (j<cgf->TileMapCache.variant[1].size()) fprintf(ofp, "%i", cgf->TileMapCache.variant[1][j]);
    }
    fprintf(ofp, "]]");

    prv = cgf->TileMapCache.offset[i];

  }
  fprintf(ofp, "  ],\n");

  fprintf(ofp, "  \"TileStepCount\":[");
  for (i=0; i<cgf->TileStepCount.size(); i++) {
    if ((i%16)==0) {
      if (i>0) { fprintf(ofp, ","); }
      fprintf(ofp, "\n    ");
    }
    else {
      fprintf(ofp, ",");
    }

    fprintf(ofp, "%llu", (unsigned long long int)cgf->TileStepCount[i]);

  }
  fprintf(ofp, "]\n");


  fprintf(ofp, "}\n");
}

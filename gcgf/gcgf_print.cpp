#include "gcgf.hpp"

void gcgf_print(cgf_t *cgf) {
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
      printf(" %02x", tp->ExtraData[j]);
    }
    printf("\n");

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



// TODO:
//
// - cache overflow
//
// THIS FUNCTION IS NOT FINISHED UNTIL CACHE OVEFLOW IS POPULATED
//
void cgf_output_band_format(cgf_t *cgf, int tilepath_idx, FILE *fp, int hiq) {
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

  mk_vec_tilemap(tilemap_vec, cgf->TileMap.c_str());

  for (i=0; i<8; i++) { hexit[8] = 0; }

  ntile = (int)cgf->TileStepCount[tilepath_idx];
  stride = (int)cgf->Stride;

  if (tilepath_idx>0) { start_offset = cgf->StrideOffset[tilepath_idx-1]; }
  loq = &(cgf->Loq[start_offset*stride]);
  span = &(cgf->Span[start_offset*stride]);
  canon = &(cgf->Canon[start_offset*stride]);
  cache_ovf = &(cgf->CacheOverflow[start_offset*stride]);

  //DEBUG
  //printf("start_offset: %i\n", start_offset);
  //printf("ntile: %i\n", ntile);
  //printf("stride: %i\n", stride);

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

  //DEBUG
  //printf("ntile: %i, n_q: %i (%i)\n", ntile, n_q, n_q*8*stride);
  //loc_verbose=1;

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

    //DEBUG
    //printf("ii: %i\n", ii);
    //printf("  canon_mask:  %08x\n", canon_mask);
    //printf("  anchor_mask: %08x\n", anchor_mask);
    //printf("  cach_mask:   %08x\n", cache_ovf_mask);

    for (i=0; i<8; i++) {
      hexit[i] = (int)((lo_cache & (0xf<<(4*i))) >> (4*i));
    }

    cur_hexit=0;
    n_hexit=0;
    n_cache_ovf = 0;
    for (i=0; i<32; i++) {
      if (cache_ovf_mask & (1<<i)) {

        //DEBUG
        //printf("    i%i, cur_hexit: %i\n", i, cur_hexit);

        if ((cur_hexit<8) && (hexit[cur_hexit] != 0xf)) {

          //DEBUG
          //printf("    ... cur_hexit %i\n", cur_hexit);fflush(stdout);

          int hexit_val = hexit[cur_hexit];

          //DEBUG
          //printf("    ... hexit_val %i\n", hexit_val); fflush(stdout);

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

}


/*
void cgf_output_band_format(cgf_t *cgf, int tilepath_idx, FILE *fp, int hiq) {
  uint64_t start_pos, pos, block_stride=0, ii, jj;
  int ovf_count=0;
  unsigned char loq_mask, hiq_mask, xspan_mask,
                cache_mask, canon_mask, anchor_mask, cache_ovf_mask;

  std::vector< int > variant_v[2];
  std::vector< std::vector< int > > noc_v[2];
  std::vector< int > v;

  int hexit[8], cur_hexit, n_hexit, n_cache_ovf;
  unsigned int s;
  unsigned char tmask;

  // Initialize
  //
  for (ii=0; ii<cgf->TileStepCount[tilepath_idx]; ii++) {
    variant_v[0].push_back(-1);
    variant_v[1].push_back(-1);
    noc_v[0].push_back(v);
    noc_v[1].push_back(v);
  }

  block_stride = (uint64_t)cgf->Stride;

  start_pos = 0;
  if (tilepath_idx>0) {
    start_pos = cgf->StrideOffset[tilepath_idx-1];
  }
  for (pos = start_pos;
       pos < cgf->StrideOffset[tilepath_idx];
       pos += block_stride) {

    ovf_count=0;
    for (ii=0; ii<block_stride; ii++) {

      hiq_mask = ~cgf->Loq[pos+ii];
      canon_mask = cgf->Canon[pos+ii];


      for (s=0; s<8; s++) {
        tmask = (1<<s);

        if ((hiq_mask & tmask) && (canon_mask & tmask)) {
          variant_v[0][pos+ii] = 0;
          variant_v[1][pos+ii] = 0;
        }


      }

    }

  }

  printf("[");
  for (ii=0; ii<variant_v[0].size(); ii++) {
    printf(" %i", variant_v[0][ii]);
  }
  printf("]\n");

  printf("[");
  for (ii=0; ii<variant_v[1].size(); ii++) {
    printf(" %i", variant_v[1][ii]);
  }
  printf("]\n");

  if (!hiq) {

    printf("[");
    for (ii=0; ii<noc_v[0].size(); ii++) {
      printf("[");
      for (jj=0; jj<noc_v[0].size(); jj++) {
        printf(" %i", noc_v[0][ii][jj]);
      }
      printf("]\n");
    }
    printf("]\n");

    printf("[");
    for (ii=0; ii<noc_v[1].size(); ii++) {
      printf("[");
      for (jj=0; jj<noc_v[1].size(); jj++) {
        printf(" %i", noc_v[1][ii][jj]);
      }
      printf("]\n");
    }
    printf("]\n");

  }

}
*/

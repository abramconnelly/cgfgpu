#include "gcgf.hpp"

static uint64_t _calc_tilepath_size(tilepath_t *tilepath) {
  int i, j, k;
  uint32_t u32;
  uint64_t u64;
  uint64_t byte_count = 0, n_8, n_32;

  // Name
  //
  byte_count += (uint64_t)(sizeof(uint32_t) + tilepath->Name.size());

  // ExtraDataSize
  //
  byte_count += (uint64_t)sizeof(uint64_t);

  // ExtraData
  //
  byte_count += (uint64_t)(sizeof(char)*(tilepath->ExtraDataSize));

  // uint64_t sizes of the SDSL arrays
  //
  byte_count += (uint64_t)(sizeof(uint64_t)*10);

  // The SDSL arrays themselves
  //
  byte_count += (uint64_t)(sdsl::size_in_bytes(tilepath->LoqTileStepHom));
  byte_count += (uint64_t)(sdsl::size_in_bytes(tilepath->LoqTileVariantHom));
  byte_count += (uint64_t)(sdsl::size_in_bytes(tilepath->LoqTileNocSumHom));
  byte_count += (uint64_t)(sdsl::size_in_bytes(tilepath->LoqTileNocStartHom));
  byte_count += (uint64_t)(sdsl::size_in_bytes(tilepath->LoqTileNocLenHom));

  byte_count += (uint64_t)(sdsl::size_in_bytes(tilepath->LoqTileStepHet));
  byte_count += (uint64_t)(sdsl::size_in_bytes(tilepath->LoqTileVariantHet));
  byte_count += (uint64_t)(sdsl::size_in_bytes(tilepath->LoqTileNocSumHet));
  byte_count += (uint64_t)(sdsl::size_in_bytes(tilepath->LoqTileNocStartHet));
  byte_count += (uint64_t)(sdsl::size_in_bytes(tilepath->LoqTileNocLenHet));

  return byte_count;

}


void ez_add_tilepath(cgf_t *cgf, int tilepath_idx, tilepath_ez_t &ez) {
  int i;
  size_t sz;
  unsigned char u8;
  uint32_t u32, i32;
  uint64_t u64, i64, sz64;
  std::vector< uint32_t > data;

  char buf[65];

  tilepath_t empty_tilepath;
  tilepath_t *tp;

  sz = (ez.N + 7)/8;

  //DEBUG
  //printf("ez_add_tilepath: N %i, sz %i\n", (int)ez.N, (int)sz);

  for (i64=0; i64<sz; i64++) {
    cgf->Loq.push_back(ez.loq_bv[i64]);
    cgf->Span.push_back(ez.span_bv[i64]);
  }

  for (i64=0; i64<(4-(sz%4)); i64++) {
    cgf->Loq.push_back(ez.loq_bv[i64]);
    cgf->Span.push_back(ez.span_bv[i64]);
  }

  sz = (ez.N + 31)/32;
  for (i64=0; i64<sz; i64++) {
    u32 = (uint32_t)( (ez.cache[i64] & 0xffffffff00000000) >> 32 );
    for (i32=0; i32<4; i32++) {
      u8 = (unsigned char)((u32>>(i32*8)) & 0xff);
      cgf->Canon.push_back(u8);
    }

    u32 = (uint32_t)( ez.cache[i64] & 0xffffffff );
    for (i32=0; i32<4; i32++) {
      u8 = (unsigned char)((u32>>(i32*8)) & 0xff);
      cgf->CacheOverflow.push_back(u8);
    }

  }

  //DEBUG
//  printf("---------------------------\n\n");
//  printf("%i %i %i %i\n",
//      (int)cgf->Loq.size(),
//      (int)cgf->Span.size(),
//      (int)cgf->Canon.size(),
//      (int)cgf->CacheOverflow.size());
//
//  //printf(" %02x %02x %02x %02x | %02x %02x %02x %02x | %02x %02x %02x %02x | %02x %02x %02x %02x",
//  printf("         .     Loq     |     Span    |    Canon    |    Ovf    \n");
//  sz = cgf->Loq.size();
//  for (i64=0; i64<sz; i64+=4) {
//    printf("[%08llx]", (unsigned long long int)i64);
//    printf(" %02x %02x %02x %02x | %02x %02x %02x %02x | %02x %02x %02x %02x | %02x %02x %02x %02x\n",
//        cgf->Loq[i64+0], cgf->Loq[i64+1], cgf->Loq[i64+2], cgf->Loq[i64+3],
//        cgf->Span[i64+0], cgf->Span[i64+1], cgf->Span[i64+2], cgf->Span[i64+3],
//        cgf->Canon[i64+0], cgf->Canon[i64+1], cgf->Canon[i64+2], cgf->Canon[i64+3],
//        cgf->CacheOverflow[i64+0], cgf->CacheOverflow[i64+1], cgf->CacheOverflow[i64+2], cgf->CacheOverflow[i64+3]);
//  }
//  printf("\n\n");

  cgf->TileStepCount.push_back((uint64_t)ez.N);

  sz64=(ez.N+31)/32;
  //if ((sz%32)!=0) { sz64+=32-(sz%32); }
  //sz64=32*sz;
  if (cgf->StrideOffset.size()>0) {
    sz = cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    sz64 += (uint64_t)sz;
  }
  cgf->StrideOffset.push_back(sz64);


  if (ez.ovf_vec.size()>0) {
    //cgf->Overflow.resize(ez.ovf_vec.size());
    for (i=0; i<ez.ovf_vec.size(); i++) {
      //cgf->Overflow[i] = ez.ovf_vec[i];
      cgf->Overflow.push_back(ez.ovf_vec[i]);
    }

    sz64=0;
    if (cgf->OverflowOffset.size()>0) { sz64 = cgf->OverflowOffset[ cgf->OverflowOffset.size()-1 ]; }
    sz64+=(uint64_t)ez.ovf_vec.size();
    cgf->OverflowOffset.push_back(sz64);
  } else {
    sz64=0;
    if (cgf->OverflowOffset.size()>0) { sz64 = cgf->OverflowOffset[ cgf->OverflowOffset.size()-1 ]; }
    cgf->OverflowOffset.push_back(sz64);
  }

  if (ez.ovf64_vec.size()>0) {
    //cgf->Overflow64.resize(ez.ovf64_vec.size());
    for (i=0; i<ez.ovf64_vec.size(); i++) {
      //cgf->Overflow64[i] = ez.ovf64_vec[i];
      cgf->Overflow64.push_back(ez.ovf64_vec[i]);
    }
    sz64=0;
    if (cgf->Overflow64Offset.size()>0) { sz64 = cgf->Overflow64Offset[ cgf->Overflow64Offset.size()-1 ]; }
    sz64+=(uint64_t)ez.ovf64_vec.size();
    cgf->Overflow64Offset.push_back(sz64);
  } else {
    sz64=0;
    if (cgf->Overflow64Offset.size()>0) { sz64 = cgf->Overflow64Offset[ cgf->Overflow64Offset.size()-1 ]; }
    cgf->Overflow64Offset.push_back(sz64);
  }

  //TODO: add low quality and extra data fields
  //

  cgf->TilePath.push_back(empty_tilepath);
  tp = &(cgf->TilePath[ cgf->TilePath.size()-1 ]);

  sprintf(buf, "%04x", tilepath_idx);
  tp->Name.clear();
  tp->Name = buf;

  tp->ExtraDataSize = 0;
  tp->ExtraData.clear();

  tp->LoqTileStepHomSize = 0;
  tp->LoqTileVariantHomSize = 0;
  tp->LoqTileNocSumHomSize = 0;
  tp->LoqTileNocStartHomSize = 0;
  tp->LoqTileNocLenHomSize = 0;

  tp->LoqTileStepHetSize = 0;
  tp->LoqTileVariantHetSize = 0;
  tp->LoqTileNocSumHetSize = 0;
  tp->LoqTileNocStartHetSize = 0;
  tp->LoqTileNocLenHetSize = 0;

  ez_to_tilepath(tp, &ez);

  sz64 = _calc_tilepath_size(tp);
  if (cgf->TilePathStructOffset.size()>0) {
    sz64 += cgf->TilePathStructOffset[ cgf->TilePathStructOffset.size()-1 ];
  }
  cgf->TilePathStructOffset.push_back(sz64);


  cgf->TilePathCount++;
}

int cgf_read_band_tilepath(cgf_t *cgf, int idx, FILE *fp) {
  int i, j, k, ch=1;
  int read_line = 0;
  int step=0;

  std::vector<std::string> names;
  std::string s;

  std::vector<tilepath_vec_t> ds;
  tilepath_vec_t cur_ds;

  int pcount=0;
  int state_mod = 0;
  int cur_allele = 0;
  int loq_flag = 0;
  int cur_tilestep = 0;

  const char *fn_tilemap = "default_tile_map_v0.1.0.txt";
  std::map< std::string, int > tilemap;
  std::map< std::string, int >::iterator ent;

  std::vector<int> loq_vec;

  load_tilemap(cgf->TileMap, tilemap);

  s.clear();
  while (ch!=EOF) {
    ch = fgetc(fp);

    if (ch==EOF) { break; }
    if (ch=='\n') {
      state_mod = (state_mod+1)%4;
      pcount=0;
      cur_tilestep=0;

      loq_vec.clear();

      if (state_mod==0) {
        cur_tilestep=0;
        ds.push_back(cur_ds);
        cur_ds.allele[0].clear();
        cur_ds.allele[1].clear();
        cur_ds.loq_flag[0].clear();
        cur_ds.loq_flag[1].clear();
        cur_ds.loq_info[0].clear();
        cur_ds.loq_info[1].clear();
        cur_ds.name.clear();
      }

      continue;
    }
    if (ch=='[') {

      loq_flag=0;
      if (state_mod>=2) {

        s.clear();
        loq_vec.clear();

        pcount++;
        while (pcount>1) {
          cur_allele = state_mod%2;
          ch = fgetc(fp);

          if (ch==EOF) {
            printf("ERROR: premature eof\n");
            return -1;
          }

          if ((ch==' ') || (ch==']')) {
            if (s.size() > 0) {
              loq_vec.push_back(atoi(s.c_str()));
            }
            s.clear();
          }

          if (ch==']') {
            cur_ds.loq_flag[cur_allele].push_back(loq_flag);
            pcount--;

            cur_ds.loq_info[cur_allele].push_back(loq_vec);
            loq_vec.clear();
            cur_tilestep++;
            continue;
          }
          if (ch=='[') { pcount++; continue; }
          if (ch==' ') { continue; }

          s += ch;
          loq_flag=1;
        }
      }

      continue;
    }
    if ((ch==' ') || (ch==']')) {

      if (s.size() == 0) { continue; }

      if (state_mod==0) {
        cur_ds.allele[0].push_back(atoi(s.c_str()));
      } else if (state_mod==1) {
        cur_ds.allele[1].push_back(atoi(s.c_str()));
      }
      s.clear();
      continue;
    }
    s += ch;

  }

  tilepath_ez_t ez;
  ez_create(ez, ds[0], tilemap);

  //ez_print(ez);

  ez_add_tilepath(cgf, idx, ez);

  //printf("\n\n");
  //gcgf_print(cgf);


  return 0;

}

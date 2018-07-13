#include "cgf4.hpp"

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
#ifdef USE_SDSL
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
#endif

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

  for (i64=0; i64<sz; i64++) {
    cgf->Loq.push_back(ez.loq_bv[i64]);
    cgf->Span.push_back(ez.span_bv[i64]);
  }

  if ((sz%4)!=0) {
    for (i64=0; i64<(4-(sz%4)); i64++) {
      cgf->Loq.push_back(0xff);
      cgf->Span.push_back(0x00);
    }
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

  cgf->TileStepCount.push_back((uint64_t)ez.N);

  sz64=(ez.N+31)/32;
  if (cgf->StrideOffset.size()>0) {
    sz = cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    sz64 += (uint64_t)sz;
  }
  cgf->StrideOffset.push_back(sz64);


  if (ez.ovf_vec.size()>0) {
    for (i=0; i<ez.ovf_vec.size(); i++) {
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
    for (i=0; i<ez.ovf64_vec.size(); i++) {
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

#ifdef USE_SDSL
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
#endif

  sz64 = _calc_tilepath_size(tp);
  if (cgf->TilePathStructOffset.size()>0) {
    sz64 += cgf->TilePathStructOffset[ cgf->TilePathStructOffset.size()-1 ];
  }
  cgf->TilePathStructOffset.push_back(sz64);

  cgf->TilePathCount++;
}

int cgf_read_band(FILE *fp, tilepath_vec_t &ds) {
  int i, j, k, ch=1;
  std::string s;

  int pcount=0;
  int state_mod = 0;
  int cur_allele = 0;
  int loq_flag = 0;
  int cur_tilestep = 0;

  std::vector<int> loq_vec;

  s.clear();
  while (ch!=EOF) {
    ch = fgetc(fp);

    if (ch==EOF) { break; }
    if (ch=='\n') {
      state_mod = (state_mod+1)%4;
      if (state_mod==0) { return 0; }

      pcount=0;
      cur_tilestep=0;
      loq_vec.clear();

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

          // premature EOF
          //
          if (ch==EOF) { return -1; }

          if ((ch==' ') || (ch==']')) {
            if (s.size() > 0) {
              loq_vec.push_back(atoi(s.c_str()));
            }
            s.clear();
          }

          if (ch==']') {
            ds.loq_flag[cur_allele].push_back(loq_flag);
            pcount--;

            ds.loq_info[cur_allele].push_back(loq_vec);
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
        ds.allele[0].push_back(atoi(s.c_str()));
      } else if (state_mod==1) {
        ds.allele[1].push_back(atoi(s.c_str()));
      }
      s.clear();
      continue;
    }
    s += ch;

  }

  return 0;
}

int cgf_read_band_tilepath(FILE *fp, cgf_t *cgf, int idx) {
  int r;
  tilepath_vec_t ds;
  tilepath_ez_t ez;
  std::map< std::string, int > tilemap;
  std::map< std::string, int >::iterator ent;

  //---
  //

  load_tilemap(cgf->TileMap, tilemap);

  r = cgf_read_band(fp, ds);
  if (r<0) { return r; }

  ez_create(ez, ds, tilemap);
  ez_add_tilepath(cgf, idx, ez);

  return 0;
}

static int _invert_noc(tilepath_vec_t &tpv) {
  int i, j, knot_len, is_hiq=0;
  size_t n;


  n = tpv.allele[0].size();
  for (i=0; i<n; i++) {
    tpv.loq_flag[0][i] = 1;
    tpv.loq_flag[1][i] = 1;
  }

  for (i=0; i<n; i+=knot_len) {
    is_hiq=0;
    knot_len=0;
    do {
      if ((tpv.loq_info[0][i+knot_len].size() > 0) ||
          (tpv.loq_info[1][i+knot_len].size() > 0)) {
        is_hiq=1;
      }
      knot_len++;
    } while ( ((i+knot_len)<n) &&
              ( (tpv.allele[0][i+knot_len] < 0) ||
                (tpv.allele[1][i+knot_len] < 0) ) );

    if (is_hiq) {
      for (j=0; j<knot_len; j++) { tpv.loq_flag[0][i+j] = 0; }
      for (j=0; j<knot_len; j++) { tpv.loq_flag[1][i+j] = 0; }
    }
  }

  for (i=0; i<n; i++) {
    tpv.loq_info[0][i].clear();
    tpv.loq_info[1][i].clear();

    if (tpv.loq_flag[0][i]==1) {
      tpv.loq_info[0][i].push_back(0);
      tpv.loq_info[0][i].push_back(0);
      tpv.loq_info[1][i].push_back(0);
      tpv.loq_info[1][i].push_back(0);
    }
  }

  return 0;
}


int cgf_read_genotype_band_tilepath(FILE *fp, cgf_t *cgf, int idx, int gtz_flag) {
  int i, j, r, a;
  tilepath_vec_t ds, orig_ds;
  tilepath_ez_t ez;
  std::map< std::string, int > tilemap;
  std::map< std::string, int >::iterator ent;
  tilepath_t *tp;

  std::vector< int > empty_vec;

  std::vector< uint32_t > gt_pos_info;
  uint32_t u32;

  uint8_t *u8v;

  //int gtz_flag = 1;
  std::vector< unsigned char > zbuf;

#ifdef USE_ZSTREAM
  z_stream defz;
#endif

  //---
  //

  load_tilemap(cgf->TileMap, tilemap);

  r = cgf_read_band(fp, ds);
  if (r<0) { return r; }

  // 'reverse' band
  //
  orig_ds = ds;
  _invert_noc(ds);

  ez_create(ez, ds, tilemap);
  ez_add_tilepath(cgf, idx, ez);

  // clear all data from the tilepath (low quality)
  // structure
  //
  tp = &(cgf->TilePath[ cgf->TilePath.size()-1 ]);

  //--------
  //--------
  //--------

  // ecnode genotype position information
  // as tile step and offset integers
  //
  tp->ExtraDataSize = 0;
  tp->ExtraData.clear();

  // first allele
  //

  if (gtz_flag) {

    gt_pos_info.clear();

    for (a=0; a<2; a++) {

      tp->ExtraData.push_back('g');
      tp->ExtraData.push_back('t');
      tp->ExtraData.push_back('z');
      tp->ExtraData.push_back( (a==0) ? '0' : '1');

      gt_pos_info.clear();
      for (i=0; i<orig_ds.loq_info[a].size(); i++) {
        for (j=0; j<orig_ds.loq_info[a][i].size(); j+=2) {
          gt_pos_info.push_back((uint32_t)i);
          gt_pos_info.push_back( (orig_ds.loq_info[a][i][j] < 0) ? (uint32_t)SPAN_SDSL_ENC_VAL : (uint32_t)orig_ds.loq_info[a][i][j] );
        }
      }

#ifdef USE_ZSTREAM
      zbuf.clear();
      uLong ulsz = compressBound( gt_pos_info.size()*sizeof(uint32_t) );
      zbuf.resize(ulsz);

      compress(&(zbuf[0]), &ulsz, (const Bytef *)(&(gt_pos_info[0])), gt_pos_info.size()*sizeof(uint32_t));

      u32 = (uint32_t)ulsz;
      u8v = (uint8_t *)(&u32);
      for (j=0; j<4; j++) { tp->ExtraData.push_back(u8v[j]); }
      for (i=0; i<ulsz; i++) { tp->ExtraData.push_back(zbuf[i]); }
#endif

    }

  }
  else {

    for (a=0; a<2; a++) {
      tp->ExtraData.push_back('g');
      tp->ExtraData.push_back('t');
      tp->ExtraData.push_back((a==0) ? '0' : '1');
      tp->ExtraData.push_back('.');

      gt_pos_info.clear();
      for (i=0; i<orig_ds.loq_info[a].size(); i++) {
        for (j=0; j<orig_ds.loq_info[a][i].size(); j+=2) {
          gt_pos_info.push_back((uint32_t)i);
          gt_pos_info.push_back( (orig_ds.loq_info[a][i][j] < 0) ? (uint32_t)SPAN_SDSL_ENC_VAL : (uint32_t)orig_ds.loq_info[a][i][j] );
        }
      }

      u32 = (uint32_t)(gt_pos_info.size()*sizeof(uint32_t));
      u8v = (uint8_t *)(&u32);
      tp->ExtraData.push_back(u8v[0]);
      tp->ExtraData.push_back(u8v[1]);
      tp->ExtraData.push_back(u8v[2]);
      tp->ExtraData.push_back(u8v[3]);

      for (i=0; i<gt_pos_info.size(); i++) {
        u32 = gt_pos_info[i];
        u8v = (uint8_t *)(&u32);
        for (j=0; j<4; j++) { tp->ExtraData.push_back(u8v[j]); }
      }

    }

  }

  tp->ExtraDataSize = tp->ExtraData.size();


  //--------
  //--------
  //--------


  tp->LoqTileStepHomSize      = 0;
  tp->LoqTileVariantHomSize   = 0;
  tp->LoqTileNocSumHomSize    = 0;
  tp->LoqTileNocStartHomSize  = 0;
  tp->LoqTileNocLenHomSize    = 0;

  tp->LoqTileStepHetSize      = 0;
  tp->LoqTileVariantHetSize   = 0;
  tp->LoqTileNocSumHetSize    = 0;
  tp->LoqTileNocStartHetSize  = 0;
  tp->LoqTileNocLenHetSize    = 0;

#ifdef USE_SDSL
  ez_create_enc_vector(tp->LoqTileStepHom,      empty_vec);
  ez_create_vlc_vector(tp->LoqTileVariantHom,   empty_vec);
  ez_create_enc_vector(tp->LoqTileNocSumHom,    empty_vec);
  ez_create_vlc_vector(tp->LoqTileNocStartHom,  empty_vec);
  ez_create_vlc_vector(tp->LoqTileNocLenHom,    empty_vec);

  ez_create_enc_vector(tp->LoqTileStepHet,      empty_vec);
  ez_create_vlc_vector(tp->LoqTileVariantHet,   empty_vec);
  ez_create_enc_vector(tp->LoqTileNocSumHet,    empty_vec);
  ez_create_vlc_vector(tp->LoqTileNocStartHet,  empty_vec);
  ez_create_vlc_vector(tp->LoqTileNocLenHet,    empty_vec);

  tp->LoqTileStepHomSize      = (uint64_t)(sdsl::size_in_bytes(tp->LoqTileStepHom));
  tp->LoqTileVariantHomSize   = (uint64_t)(sdsl::size_in_bytes(tp->LoqTileVariantHom));
  tp->LoqTileNocSumHomSize    = (uint64_t)(sdsl::size_in_bytes(tp->LoqTileNocSumHom));
  tp->LoqTileNocStartHomSize  = (uint64_t)(sdsl::size_in_bytes(tp->LoqTileNocStartHom));
  tp->LoqTileNocLenHomSize    = (uint64_t)(sdsl::size_in_bytes(tp->LoqTileNocLenHom));

  tp->LoqTileStepHetSize      = (uint64_t)(sdsl::size_in_bytes(tp->LoqTileStepHet));
  tp->LoqTileVariantHetSize   = (uint64_t)(sdsl::size_in_bytes(tp->LoqTileVariantHet));
  tp->LoqTileNocSumHetSize    = (uint64_t)(sdsl::size_in_bytes(tp->LoqTileNocSumHet));
  tp->LoqTileNocStartHetSize  = (uint64_t)(sdsl::size_in_bytes(tp->LoqTileNocStartHet));
  tp->LoqTileNocLenHetSize    = (uint64_t)(sdsl::size_in_bytes(tp->LoqTileNocLenHet));
#endif

  return 0;
}

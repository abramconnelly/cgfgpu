#include "cgf4.hpp"

// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------


static int _write_string(std::string &s, FILE *fp) {
  uint32_t u32;
  size_t sz;

  u32 = (uint32_t)s.size();
  sz = fwrite(&(u32), sizeof(uint32_t), 1, fp);
  if (sz!=1) { return -1; }

  sz = fwrite(&(s[0]), sizeof(char), s.size(), fp);
  if (sz!=s.size()) { return -1; }

  return 0;
}

struct membuf : std::streambuf {
  membuf(char *beg, char *end) { this->setg(beg, beg, end); }
};

#ifdef USE_SDSL

static int _read_vlc_vector(sdsl::vlc_vector<> &vlc, size_t sz, FILE *fp) {
  size_t s;
  int ch;
  std::vector<char> b;

  for (s=0; s<sz; s++) {
    ch = fgetc(fp);
    if (ch==EOF) { return -1; }
    b.push_back(ch);
  }

  membuf sbuf(&(b[0]), &(b[0]) + b.size());
  std::istream in(&sbuf);
  vlc.load(in);

  return 0;
}


static int _write_vlc_vector(sdsl::vlc_vector<> &vlc, FILE *fp) {
  std::ostringstream bufstream;
  std::string s;
  size_t sz;

  vlc.serialize(bufstream);
  s = bufstream.str();

  sz = fwrite( s.c_str(), sizeof(char), s.size(), fp);
  if (sz!=s.size()) { return -1; }

  return 0;
}

static int _read_enc_vector(sdsl::enc_vector<> &enc, size_t sz, FILE *fp) {
  size_t s;
  int ch;
  std::vector<char> b;

  for (s=0; s<sz; s++) {
    ch = fgetc(fp);
    if (ch==EOF) { return -1; }
    b.push_back(ch);
  }

  membuf sbuf(&(b[0]), &(b[0]) + b.size());
  std::istream in(&sbuf);
  enc.load(in);

  return 0;
}

static int _write_enc_vector(sdsl::enc_vector<> &enc, FILE *fp) {
  std::ostringstream bufstream;
  std::string s;
  size_t sz;

  enc.serialize(bufstream);
  s = bufstream.str();

  sz = fwrite( s.c_str(), sizeof(char), s.size(), fp);
  if (sz!=s.size()) { return -1; }

  return 0;
}
#endif

static uint64_t _calc_path_size(tilepath_t *tilepath) {
  int i, j, k;
  uint32_t u32;
  uint64_t u64;
  uint64_t byte_count = 0, n_8, n_32;

  // Name
  //
  byte_count += (uint64_t)(sizeof(uint32_t) + tilepath->Name.size());

  // Extra Data
  //
  byte_count += (uint64_t)(sizeof(uint64_t) + tilepath->ExtraData.size());

  // Size structures for the low quality SDSL vectors
  //
  byte_count += (uint64_t)(sizeof(uint64_t)*10);

  // The SDSL structure sizes
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


// -----------------------------------------------------------
// -----------------------------------------------------------
// -----------------------------------------------------------

size_t gcgf_tilepath_read(cgf_t *cgf, int tilepath_idx, FILE *fp) {
  int i, j, k, ch;
  char buf[1024];
  //cgf_t *cgf=NULL;
  tilepath_t *tilepath;

  size_t sz;

  uint64_t u64, i64;
  uint32_t u32;
  unsigned char ub[32];

  uint64_t offset;

  tilepath_t empty_tilepath;

  if (cgf->TilePath.size() <= tilepath_idx) {
    for (i=cgf->TilePath.size(); i<=tilepath_idx; i++) {
      cgf->TilePath.push_back(empty_tilepath);
    }
  }

  tilepath = &(cgf->TilePath[tilepath_idx]);

  // Name of tilepath.
  //
  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->Name.clear();
  tilepath->Name.resize(u32);
  sz = fread(&(tilepath->Name[0]), sizeof(char), u32, fp);
  if (sz!=u32) { return -1; }

  // Extra Data (size and bytes
  //
  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->ExtraDataSize = u64;
  tilepath->ExtraData.resize(u64);
  if ( tilepath->ExtraData.size() > 0 ) {
	sz = fread(&(tilepath->ExtraData[0]), sizeof(char), u64, fp);
	if (sz!=u64) { return -1; }
  }

  // Size of SDSL homozygous arrays
  //

  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->LoqTileStepHomSize = u64;

  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->LoqTileVariantHomSize = u64;

  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->LoqTileNocSumHomSize = u64;

  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->LoqTileNocStartHomSize = u64;

  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->LoqTileNocLenHomSize = u64;

  // Size of SDSL heterozygous arrays
  //

  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->LoqTileStepHetSize = u64;

  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->LoqTileVariantHetSize = u64;

  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->LoqTileNocSumHetSize = u64;

  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->LoqTileNocStartHetSize = u64;

  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { return -1; }
  tilepath->LoqTileNocLenHetSize = u64;


  // Read SDSL structures in
  //
#ifdef USE_SDSL
  if (tilepath->LoqTileStepHomSize>0) {
    k = _read_enc_vector(tilepath->LoqTileStepHom, tilepath->LoqTileStepHomSize, fp);
    if (k<0) { return k; }
  }

  if (tilepath->LoqTileVariantHomSize>0) {
    k = _read_vlc_vector(tilepath->LoqTileVariantHom, tilepath->LoqTileVariantHomSize, fp);
    if (k<0) { return k; }
  }

  if (tilepath->LoqTileNocSumHomSize>0) {
    k = _read_enc_vector(tilepath->LoqTileNocSumHom, tilepath->LoqTileNocSumHomSize, fp);
    if (k<0) { return k; }
  }

  if (tilepath->LoqTileNocStartHomSize>0) {
    k = _read_vlc_vector(tilepath->LoqTileNocStartHom, tilepath->LoqTileNocStartHomSize, fp);
    if (k<0) { return k; }
  }

  if (tilepath->LoqTileNocLenHomSize>0) {
    k = _read_vlc_vector(tilepath->LoqTileNocLenHom, tilepath->LoqTileNocLenHomSize, fp);
    if (k<0) { return k; }
  }

  //--

  if (tilepath->LoqTileStepHetSize>0) {
    k = _read_enc_vector(tilepath->LoqTileStepHet, tilepath->LoqTileStepHetSize, fp);
    if (k<0) { return k; }
  }

  if (tilepath->LoqTileVariantHetSize>0) {
    k = _read_vlc_vector(tilepath->LoqTileVariantHet, tilepath->LoqTileVariantHetSize, fp);
    if (k<0) { return k; }
  }

  if (tilepath->LoqTileNocSumHetSize>0) {
    k = _read_enc_vector(tilepath->LoqTileNocSumHet, tilepath->LoqTileNocSumHetSize, fp);
    if (k<0) { return k; }
  }

  if (tilepath->LoqTileNocStartHetSize>0) {
    k = _read_vlc_vector(tilepath->LoqTileNocStartHet, tilepath->LoqTileNocStartHetSize, fp);
    if (k<0) { return k; }
  }

  if (tilepath->LoqTileNocLenHetSize>0) {
    k = _read_vlc_vector(tilepath->LoqTileNocLenHet, tilepath->LoqTileNocLenHetSize, fp);
    if (k<0) { return k; }
  }
#endif

  return 0;
}


cgf_t *cgf_read_hiq(FILE *fp) {
  int i, j, k, ch, t;
  cgf_t *cgf=NULL;

  size_t sz;

  uint64_t u64, i64, offset;
  uint32_t u32;
  unsigned char ub[32];

  cgf = new cgf_t;

  // Read Magic
  //
  sz = fread(ub, sizeof(char), 8, fp);
  if (sz!=8) { goto cgf_read_hiq_error; }
  for (i=0; i<8; i++) {
    if (ub[i] != CGF_MAGIC[i]) { goto cgf_read_hiq_error; }
    cgf->Magic[i] = ub[i];
  }

  // CGF Version
  //
  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto cgf_read_hiq_error; }

  cgf->CGFVersion.clear();
  cgf->CGFVersion.reserve(u32);
  for (i=0; i<u32; i++) {
    ch = fgetc(fp);
    if (ch==EOF) { goto cgf_read_hiq_error; }
    cgf->CGFVersion += ch;
  }

  // Library Version
  //
  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto cgf_read_hiq_error; }

  cgf->LibraryVersion.clear();
  cgf->LibraryVersion.reserve(u32);
  for (i=0; i<u32; i++) {
    ch = fgetc(fp);
    if (ch==EOF) { goto cgf_read_hiq_error; }
    cgf->LibraryVersion += ch;
  }

  // Tile Path Count
  //
  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { goto cgf_read_hiq_error; }
  cgf->TilePathCount = u64;

  // TileMap
  //
  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto cgf_read_hiq_error; }
  cgf->TileMap.clear();
  cgf->TileMap.reserve(u32);
  for (i=0; i<u32; i++) {
    ch = fgetc(fp);
    if (ch==EOF) { goto cgf_read_hiq_error; }
    cgf->TileMap += ch;
  }

  // Stride
  //

  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto cgf_read_hiq_error; }
  cgf->Stride = u32;

  if (cgf->TilePathCount>0) {

    // Tile Step Count vector
    //
    cgf->TileStepCount.clear();
    for (i=0; i<(int)cgf->TilePathCount; i++) {
      sz = fread(&u64, sizeof(uint64_t), 1, fp);
      if (sz!=1) { goto cgf_read_hiq_error; }
      cgf->TileStepCount.push_back(u64);
    }

    // Stride Offset
    //

    cgf->StrideOffset.clear();
    for (i=0; i<(int)cgf->TilePathCount; i++) {
      sz = fread(&u64, sizeof(uint64_t), 1, fp);
      if (sz!=1) { goto cgf_read_hiq_error; }
      cgf->StrideOffset.push_back(u64);
    }

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->Loq.resize(u64);
    sz = fread(&(cgf->Loq[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto cgf_read_hiq_error; }

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->Span.resize(u64);
    sz = fread(&(cgf->Span[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto cgf_read_hiq_error; }

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->Canon.resize(u64);
    sz = fread(&(cgf->Canon[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto cgf_read_hiq_error; }

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->CacheOverflow.resize(u64);
    sz = fread(&(cgf->CacheOverflow[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto cgf_read_hiq_error; }


    // simple overflow
    //

    cgf->OverflowOffset.resize( cgf->TilePathCount );
    sz = fread(&(cgf->OverflowOffset[0]), sizeof(uint64_t), cgf->TilePathCount, fp);
    if (sz!=cgf->TilePathCount) { goto cgf_read_hiq_error; }

    u64 = cgf->OverflowOffset[ cgf->OverflowOffset.size()-1 ];

    if (u64>0) {

      cgf->Overflow.resize( u64 );
      sz = fread(&(cgf->Overflow[0]), sizeof(uint16_t), u64, fp);
      if (sz!=(u64)) { goto cgf_read_hiq_error; }

    }

    // Overflow that couldn't be stored in the above
    //

    cgf->Overflow64Offset.resize( cgf->TilePathCount );
    sz = fread(&(cgf->Overflow64Offset[0]), sizeof(uint64_t), cgf->TilePathCount, fp);
    if (sz!=cgf->TilePathCount) { goto cgf_read_hiq_error; }

    u64 = cgf->Overflow64Offset[ cgf->Overflow64Offset.size()-1 ];


    if (u64>0) {
      cgf->Overflow64.resize( u64 );
      sz = fread(&(cgf->Overflow64[0]), sizeof(uint64_t), u64, fp);
      if (sz!=(u64)) { goto cgf_read_hiq_error; }
    }

  }

  return cgf;

cgf_read_hiq_error:
  if (feof(fp)) { printf("Unexpected end-of-file.\n"); }
  if (cgf) { delete cgf; }
  return NULL;

}

cgf_t *cgf_read(FILE *fp) {
  int i, j, k, ch;
  cgf_t *cgf=NULL;
  tilepath_t *tilepath;

  int t;

  size_t sz;

  uint64_t u64, i64;
  uint32_t u32;
  unsigned char ub[32];

  uint64_t offset;

  cgf = new cgf_t;

  // Read Magic
  //
  sz = fread(ub, sizeof(char), 8, fp);
  if (sz!=8) { goto gcgf_read_error; }
  for (i=0; i<8; i++) {
    if (ub[i] != CGF_MAGIC[i]) { goto gcgf_read_error; }
    cgf->Magic[i] = ub[i];
  }

  // CGF Version
  //
  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto gcgf_read_error; }

  cgf->CGFVersion.clear();
  cgf->CGFVersion.reserve(u32);
  for (i=0; i<u32; i++) {
    ch = fgetc(fp);
    if (ch==EOF) { goto gcgf_read_error; }
    cgf->CGFVersion += ch;
  }

  // Library Version
  //
  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto gcgf_read_error; }

  cgf->LibraryVersion.clear();
  cgf->LibraryVersion.reserve(u32);
  for (i=0; i<u32; i++) {
    ch = fgetc(fp);
    if (ch==EOF) { goto gcgf_read_error; }
    cgf->LibraryVersion += ch;
  }

  // Tile Path Count
  //
  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { goto gcgf_read_error; }
  cgf->TilePathCount = u64;

  // TileMap
  //
  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto gcgf_read_error; }
  cgf->TileMap.clear();
  cgf->TileMap.reserve(u32);
  for (i=0; i<u32; i++) {
    ch = fgetc(fp);
    if (ch==EOF) { goto gcgf_read_error; }
    cgf->TileMap += ch;
  }

  // Stride
  //

  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto gcgf_read_error; }
  cgf->Stride = u32;

  if (cgf->TilePathCount>0) {

    // Tile Step Count vector
    //
    cgf->TileStepCount.clear();
    for (i=0; i<(int)cgf->TilePathCount; i++) {
      sz = fread(&u64, sizeof(uint64_t), 1, fp);
      if (sz!=1) { goto gcgf_read_error; }
      cgf->TileStepCount.push_back(u64);
    }

    // Stride Offset
    //
    cgf->StrideOffset.clear();
    for (i=0; i<(int)cgf->TilePathCount; i++) {
      sz = fread(&u64, sizeof(uint64_t), 1, fp);
      if (sz!=1) { goto gcgf_read_error; }
      cgf->StrideOffset.push_back(u64);
    }

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->Loq.resize(u64);
    sz = fread(&(cgf->Loq[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto gcgf_read_error; }

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->Span.resize(u64);
    sz = fread(&(cgf->Span[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto gcgf_read_error; }

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->Canon.resize(u64);
    sz = fread(&(cgf->Canon[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto gcgf_read_error; }

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->CacheOverflow.resize(u64);
    sz = fread(&(cgf->CacheOverflow[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto gcgf_read_error; }

    // simple overflow
    //
    cgf->OverflowOffset.resize( cgf->TilePathCount );
    sz = fread(&(cgf->OverflowOffset[0]), sizeof(uint64_t), cgf->TilePathCount, fp);
    if (sz!=cgf->TilePathCount) { goto gcgf_read_error; }

    u64 = cgf->OverflowOffset[ cgf->OverflowOffset.size()-1 ];

    if (u64>0) {
      cgf->Overflow.resize( u64 );
      sz = fread(&(cgf->Overflow[0]), sizeof(uint16_t), u64, fp);
      if (sz!=(u64)) { goto gcgf_read_error; }
    }

    // Overflow that couldn't be stored in the above
    //
    cgf->Overflow64Offset.resize( cgf->TilePathCount );
    sz = fread(&(cgf->Overflow64Offset[0]), sizeof(uint64_t), cgf->TilePathCount, fp);
    if (sz!=cgf->TilePathCount) { goto gcgf_read_error; }

    u64 = cgf->Overflow64Offset[ cgf->Overflow64Offset.size()-1 ];

    if (u64>0) {
      cgf->Overflow64.resize( u64 );
      sz = fread(&(cgf->Overflow64[0]), sizeof(uint64_t), u64, fp);
      if (sz!=(u64)) { goto gcgf_read_error; }
    }

    cgf->TilePathStructOffset.resize(cgf->TilePathCount);
    sz = fread(&(cgf->TilePathStructOffset[0]), sizeof(uint64_t), cgf->TilePathCount, fp);
    if (sz!=cgf->TilePathCount) { goto gcgf_read_error; }

    for (i64=0; i64<cgf->TilePathCount; i64++) {
      sz = gcgf_tilepath_read(cgf, (int)i64, fp);
    }

  }

  return cgf;

gcgf_read_error:
  if (cgf) { delete cgf; }
  return NULL;
}

void cgf_create_container(FILE *fp,
                          const char *cgf_version,
                          const char *cglf_version,
                          const char *tilemap) {
  uint32_t u32;
  uint64_t u64;
  int i, n;

  cgf_version = (cgf_version ? cgf_version : CGF_VERSION);
  cglf_version = (cglf_version ? cglf_version : CGLF_VERSION);

  // CGF magic string
  //
  fwrite(CGF_MAGIC, sizeof(char), 8, fp);

  // CGF version
  //
  u32 = (uint32_t)strlen(cgf_version);
  fwrite(&u32, sizeof(uint32_t), 1, fp);
  fwrite(cgf_version, sizeof(char), u32, fp);

  // Library version
  //
  u32 = (uint32_t)strlen(cglf_version);
  fwrite(&u32, sizeof(uint32_t), 1, fp);
  fwrite(cglf_version, sizeof(char), u32, fp);

  // TilePathCount
  //
  u64 = 0;
  fwrite(&u64, sizeof(uint64_t), 1, fp);

  // TileMap (as string)
  u32 = (uint32_t) strlen(tilemap);
  fwrite(&u32, sizeof(uint32_t), 1, fp);
  fwrite(tilemap, sizeof(char), u32, fp);

  // Stride
  //
  u32=4; fwrite(&u32, sizeof(uint32_t), 1, fp);

}


//---


uint64_t cgf_write_to_file(cgf_t *cgf, const char *ofn) {
  int i, j, k, n, ii;
  size_t sz;
  uint64_t u64;
  uint32_t u32;
  uint16_t u16;
  FILE *ofp;
  char c;
  int n_8, n_32;
  tilepath_t *tilepath;
  uint64_t byte_count=0;

  unsigned char *xx;

  if ((ofp=fopen(ofn, "w"))==NULL) { return -1; }

  for (i=0; i<8; i++) {
    sz = fwrite(&(cgf->Magic[i]), sizeof(char), 1, ofp);
    if (sz!=1) { return -1; }
    byte_count++;
  }

  k = _write_string(cgf->CGFVersion, ofp);
  if (k<0) { return k; }
  byte_count+=(uint64_t)sizeof(uint32_t) + (uint64_t)cgf->CGFVersion.size();

  k = _write_string(cgf->LibraryVersion, ofp);
  if (k<0) { return k; }
  byte_count+=(uint64_t)sizeof(uint32_t) + (uint64_t)cgf->LibraryVersion.size();

  sz = fwrite(&(cgf->TilePathCount), sizeof(uint64_t), 1, ofp);
  if (sz!=1) { return -1; }
  byte_count+=(uint64_t)sizeof(uint64_t);

  k = _write_string(cgf->TileMap, ofp);
  if (k<0) { return k; }
  byte_count+=(uint64_t)sizeof(uint32_t) + (uint64_t)cgf->TileMap.size();

  sz = fwrite(&(cgf->Stride), sizeof(uint32_t), 1, ofp);
  if (sz!=1) { return -1; }
  byte_count+=(uint64_t)sizeof(uint32_t);

  if (cgf->TilePathCount>0) {

    sz = fwrite(&(cgf->TileStepCount[0]), sizeof(uint64_t), cgf->TileStepCount.size(), ofp);
    if (sz!=cgf->TileStepCount.size()) { return -1; }
    byte_count += (uint64_t)(sizeof(uint64_t)*(cgf->TileStepCount.size()));

    sz = fwrite(&(cgf->StrideOffset[0]), sizeof(uint64_t), cgf->StrideOffset.size(), ofp);
    if (sz!=cgf->StrideOffset.size()) { return -1; }
    byte_count += (uint64_t)(sizeof(uint64_t)*(cgf->StrideOffset.size()));

    sz = fwrite(&(cgf->Loq[0]), sizeof(unsigned char), cgf->Loq.size(), ofp);
    if (sz!=cgf->Loq.size()) { return -1; }
    byte_count += (uint64_t)(cgf->Loq.size());

    sz = fwrite(&(cgf->Span[0]), sizeof(unsigned char), cgf->Span.size(), ofp);
    if (sz!=cgf->Span.size()) { return -1; }
    byte_count += (uint64_t)(cgf->Span.size());

    sz = fwrite(&(cgf->Canon[0]), sizeof(unsigned char), cgf->Canon.size(), ofp);
    if (sz!=cgf->Canon.size()) { return -1; }
    byte_count += (uint64_t)(cgf->Canon.size());

    sz = fwrite(&(cgf->CacheOverflow[0]), sizeof(unsigned char), cgf->CacheOverflow.size(), ofp);
    if (sz!=cgf->CacheOverflow.size()) { return -1; }
    byte_count += (uint64_t)(cgf->CacheOverflow.size());

    //--

    sz = fwrite(&(cgf->OverflowOffset[0]), sizeof(uint64_t), cgf->OverflowOffset.size(), ofp);
    if (sz!=cgf->OverflowOffset.size()) { return -1; }
    byte_count += (uint64_t)(cgf->OverflowOffset.size());

    if (cgf->Overflow.size()>0) {
      sz = fwrite(&(cgf->Overflow[0]), sizeof(uint16_t), cgf->Overflow.size(), ofp);
      if (sz!=cgf->Overflow.size()) { return -1; }
      byte_count += (uint64_t)(cgf->Overflow.size());
    }

    //--

    sz = fwrite(&(cgf->Overflow64Offset[0]), sizeof(uint64_t), cgf->Overflow64Offset.size(), ofp);
    if (sz!=cgf->Overflow64Offset.size()) { return -1; }
    byte_count += (uint64_t)(cgf->Overflow64Offset.size());

    if (cgf->Overflow64.size()>0) {
      sz = fwrite(&(cgf->Overflow64[0]), sizeof(uint64_t), cgf->Overflow64.size(), ofp);
      if (sz!=cgf->Overflow64.size()) { return -1; }
      byte_count += (uint64_t)(cgf->Overflow64.size());
    }

    //--

    sz = fwrite(&(cgf->TilePathStructOffset[0]), sizeof(uint64_t), cgf->TilePathStructOffset.size(), ofp);
    if (sz!=cgf->TilePathStructOffset.size()) { return -1; }
    byte_count += (uint64_t)(cgf->TilePathStructOffset.size());

    for (ii=0; ii<cgf->TilePathStructOffset.size(); ii++) {
      tilepath = &(cgf->TilePath[ii]);

      k = _write_string(tilepath->Name, ofp);
      if (k<0) { return k; }
      byte_count+=(uint64_t)sizeof(uint32_t) + (uint64_t)(tilepath->Name.size());

      sz = fwrite(&(tilepath->ExtraDataSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }
      byte_count+=(uint64_t)sizeof(uint64_t);

      sz = fwrite(&(tilepath->ExtraData[0]), sizeof(char), tilepath->ExtraData.size(), ofp);
      if (sz!=cgf->TilePath[ii].ExtraData.size()) { return -1; }
      byte_count+=(uint64_t)(tilepath->ExtraData.size());


      // size information for hom low quality structures
      //

      sz = fwrite(&(tilepath->LoqTileStepHomSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }

      sz = fwrite(&(tilepath->LoqTileVariantHomSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }

      sz = fwrite(&(tilepath->LoqTileNocSumHomSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }

      sz = fwrite(&(tilepath->LoqTileNocStartHomSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }

      sz = fwrite(&(tilepath->LoqTileNocLenHomSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }

      byte_count += (uint64_t)(sizeof(uint64_t)*5);

      // size information for het low quality structures
      //

      sz = fwrite(&(tilepath->LoqTileStepHetSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }

      sz = fwrite(&(tilepath->LoqTileVariantHetSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }

      sz = fwrite(&(tilepath->LoqTileNocSumHetSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }

      sz = fwrite(&(tilepath->LoqTileNocStartHetSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }

      sz = fwrite(&(tilepath->LoqTileNocLenHetSize), sizeof(uint64_t), 1, ofp);
      if (sz!=1) { return -1; }

      byte_count += (uint64_t)(sizeof(uint64_t)*5);

      // low quality hom structures
      //
#ifdef USE_SDSL
	  k = _write_enc_vector(tilepath->LoqTileStepHom, ofp);
      if (k<0) { return -1; }

      k = _write_vlc_vector(tilepath->LoqTileVariantHom, ofp);
      if (k<0) { return -1; }

      k = _write_enc_vector(tilepath->LoqTileNocSumHom, ofp);
      if (k<0) { return -1; }

      k = _write_vlc_vector(tilepath->LoqTileNocStartHom, ofp);
      if (k<0) { return -1; }

      k = _write_vlc_vector(tilepath->LoqTileNocLenHom, ofp);
      if (k<0) { return -1; }

      // low quality het structures
      //

      k = _write_enc_vector(tilepath->LoqTileStepHet, ofp);
      if (k<0) { return -1; }

      k = _write_vlc_vector(tilepath->LoqTileVariantHet, ofp);
      if (k<0) { return -1; }

      k = _write_enc_vector(tilepath->LoqTileNocSumHet, ofp);
      if (k<0) { return -1; }

      k = _write_vlc_vector(tilepath->LoqTileNocStartHet, ofp);
      if (k<0) { return -1; }

      k = _write_vlc_vector(tilepath->LoqTileNocLenHet, ofp);
      if (k<0) { return -1; }

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

    }

  }

  fclose(ofp);

  return byte_count;
}

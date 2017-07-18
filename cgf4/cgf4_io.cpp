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
  sz = fread(&(tilepath->ExtraData[0]), sizeof(char), u64, fp);
  if (sz!=u64) { return -1; }

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


  return 0;

}


cgf_t *cgf_read(FILE *fp) {
  int i, j, k, ch;
  char buf[1024];
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

  //DEBUG
  //printf("magic: %c%c%c%c%c%c%c%c\n",
  //    cgf->Magic[0], cgf->Magic[1], cgf->Magic[2], cgf->Magic[3],
  //    cgf->Magic[4], cgf->Magic[5], cgf->Magic[6], cgf->Magic[7]);

  // CGF Version
  //
  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto gcgf_read_error; }

  //DEBUG
  //printf(">>> %i\n", (int)u32);

  cgf->CGFVersion.clear();
  cgf->CGFVersion.reserve(u32);
  for (i=0; i<u32; i++) {
    ch = fgetc(fp);
    if (ch==EOF) { goto gcgf_read_error; }
    cgf->CGFVersion += ch;

    //DEBUG
    //printf(">> %c %i\n", ch, ch);
  }

  //DEBUG
  //printf("version: %s\n", cgf->CGFVersion.c_str());

  // Library Version
  //
  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto gcgf_read_error; }

  //DEBUG
  //printf(">>> %i\n", (int)u32);

  cgf->LibraryVersion.clear();
  cgf->LibraryVersion.reserve(u32);
  for (i=0; i<u32; i++) {
    ch = fgetc(fp);
    if (ch==EOF) { goto gcgf_read_error; }
    cgf->LibraryVersion += ch;

    //DEBUG
    //printf(">> %c %i\n", ch, ch);
  }

  //DEBUG
  //printf("libver: %s\n", cgf->LibraryVersion.c_str());


  // Tile Path Count
  //
  sz = fread(&u64, sizeof(uint64_t), 1, fp);
  if (sz!=1) { goto gcgf_read_error; }
  cgf->TilePathCount = u64;

  //DEBUG
  //printf("pathcount: %i\n", (int)cgf->TilePathCount);


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

  //DEBUG
  //t = (int)cgf->TileMap.size();
  //printf("tilemap %i '%c%c .. %c%c'\n",
  //    (int)cgf->TileMap.size(),
  //    cgf->TileMap[0], cgf->TileMap[1],
  //    cgf->TileMap[t-2], cgf->TileMap[t-1]);
  //printf("cp1\n");

  // Stride
  //

  sz = fread(&u32, sizeof(uint32_t), 1, fp);
  if (sz!=1) { goto gcgf_read_error; }
  cgf->Stride = u32;

  //DEBUG
  //printf("stride: %i\n", (int)cgf->Stride);
  //printf("cp2\n");

  if (cgf->TilePathCount>0) {

    // Tile Step Count vector
    //
    cgf->TileStepCount.clear();
    for (i=0; i<(int)cgf->TilePathCount; i++) {
      sz = fread(&u64, sizeof(uint64_t), 1, fp);
      if (sz!=1) { goto gcgf_read_error; }
      cgf->TileStepCount.push_back(u64);

      //DEBUG
      //printf("  step[%i] %llu\n", i, (unsigned long long)u64);
    }

    // Stride Offset
    //

    cgf->StrideOffset.clear();
    for (i=0; i<(int)cgf->TilePathCount; i++) {
      sz = fread(&u64, sizeof(uint64_t), 1, fp);
      if (sz!=1) { goto gcgf_read_error; }
      cgf->StrideOffset.push_back(u64);

      //DEBUG
      //printf("  ofst[%i] %llu\n", i, (unsigned long long)u64);
    }

    //DEBUG
    //printf("StrideOffset[%i]:", (int)cgf->StrideOffset.size());
    //for (i=0; i<cgf->TilePathCount; i++) {
    //  printf(" %i", (int)cgf->StrideOffset[i]);
    //}
    //printf("\n");

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->Loq.resize(u64);
    sz = fread(&(cgf->Loq[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto gcgf_read_error; }

    //printf("  loq(%i):\n", (int)cgf->Loq.size());
    //for (i=0; i<cgf->Loq.size(); i++) { printf(" %02x", cgf->Loq[i]); }
    //printf("\n");



    //DEBUIG
    //printf("loq[%i]:", (int)cgf->Loq.size());
    //for (i=0; i<cgf->Loq.size(); i++) { printf(" %02x", cgf->Loq[i]); }
    //printf("\n");

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->Span.resize(u64);
    sz = fread(&(cgf->Span[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto gcgf_read_error; }

    //printf("  span(%i):\n", (int)cgf->Span.size());
    //for (i=0; i<cgf->Span.size(); i++) { printf(" %02x", cgf->Span[i]); }
    //printf("\n");



    //DEBUIG
    //printf("span[%i]:", (int)cgf->Span.size());
    //for (i=0; i<cgf->Span.size(); i++) { printf(" %02x", cgf->Span[i]); }
    //printf("\n");

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->Canon.resize(u64);
    sz = fread(&(cgf->Canon[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto gcgf_read_error; }

    //printf("  strideoffset(%i):\n", (int)cgf->Canon.size());
    //for (i=0; i<cgf->Canon.size(); i++) { printf(" %02x", cgf->Canon[i]); }
    //printf("\n");


    //DEBUIG
    //printf("Canon[%i]:", (int)cgf->Canon.size());
    //for (i=0; i<cgf->Canon.size(); i++) { printf(" %02x", cgf->Canon[i]); }
    //printf("\n");

    u64 = ((uint64_t)cgf->Stride) * cgf->StrideOffset[ cgf->StrideOffset.size()-1 ];
    cgf->CacheOverflow.resize(u64);
    sz = fread(&(cgf->CacheOverflow[0]), sizeof(unsigned char), u64, fp);
    if (sz!=u64) { goto gcgf_read_error; }


    //DEBUIG
    //printf("cp3\n");
    //printf("ovf[%i]:\n", (int)cgf->CacheOverflow.size());
    //for (i=0; i<cgf->CacheOverflow.size(); i++) { printf(" %02x", cgf->CacheOverflow[i]); }
    //printf("\n");

    //printf("  cacheovf(%i):", (int)cgf->CacheOverflow.size());
    //for (i=0; i<(int)cgf->CacheOverflow.size(); i++) {
    //  printf(" %llu(%016llx)", (unsigned long long)cgf->CacheOverflow[i], (unsigned long long)cgf->CacheOverflow[i]);
    //}
    //printf("\n");


    // simple overflow
    //

    cgf->OverflowOffset.resize( cgf->TilePathCount );
    sz = fread(&(cgf->OverflowOffset[0]), sizeof(uint64_t), cgf->TilePathCount, fp);
    if (sz!=cgf->TilePathCount) { goto gcgf_read_error; }

    //DEBUG
    //printf("cp4\n");
    //printf("  tilepathcount %i\n", (int)cgf->TilePathCount);
    //printf("  sz %i\n", (int)cgf->OverflowOffset.size());
    //printf("  ovfoffst(%i):\n", (int)cgf->OverflowOffset.size());
    //for (i=0; i<(int)cgf->OverflowOffset.size(); i++) {
    //  printf(" %llu(%016llx)", (unsigned long long)cgf->OverflowOffset[i], (unsigned long long)cgf->OverflowOffset[i]);
    //}
    //printf("\n");

    u64 = cgf->OverflowOffset[ cgf->OverflowOffset.size()-1 ];
    //cgf->Overflow.resize( u64 * 3 );

    //DEBUG
    //printf("  u64 %llu\n", (unsigned long long)u64);

    if (u64>0) {

      cgf->Overflow.resize( u64 );
      sz = fread(&(cgf->Overflow[0]), sizeof(uint16_t), u64, fp);
      if (sz!=(u64)) { goto gcgf_read_error; }

      //DEBUG
      //printf("Overflow(%i):", (int)cgf->OverflowOffset[ cgf->OverflowOffset.size()-1 ]);
      //for (i=0; i<cgf->OverflowOffset[ cgf->OverflowOffset.size()-1 ]; i++) { printf(" %i", (int)cgf->Overflow[i]); }
      //printf("\n");


    }

    //DEBUG
    //printf("cp5\n");
    //printf("  overflowoffset(%i):", (int)cgf->OverflowOffset.size());
    //for (i=0; i<cgf->OverflowOffset.size(); i++) { printf(" %i", (int)cgf->OverflowOffset[i]); }
    //printf("\n");


    // Overflow that couldn't be stored in the above
    //

    cgf->Overflow64Offset.resize( cgf->TilePathCount );
    sz = fread(&(cgf->Overflow64Offset[0]), sizeof(uint64_t), cgf->TilePathCount, fp);
    if (sz!=cgf->TilePathCount) { goto gcgf_read_error; }

    //printf("cp6\n");
    //printf("ovf64.size: %i\n", (int)cgf->Overflow64Offset.size());


    u64 = cgf->Overflow64Offset[ cgf->Overflow64Offset.size()-1 ];

    //DEBUG
    //for (i=0; i<cgf->Overflow64Offset.size(); i++) {
    //  printf(" (%i) %llu\n", i, (unsigned long long int)cgf->Overflow64Offset[i]);
    //}

    //printf(">>>> %llu\n", (unsigned long long int)u64);

    if (u64>0) {
      //cgf->Overflow64.resize( u64 * 3 );
      cgf->Overflow64.resize( u64 );
      sz = fread(&(cgf->Overflow64[0]), sizeof(uint64_t), u64, fp);
      if (sz!=(u64)) { goto gcgf_read_error; }
    }

    //printf("cp7\n");

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

  // CGF magci string
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
  u32 = (uint32_t)strlen(tilemap);
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
	//tilepath_t *tilepath;

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

  //DEBUG
  //printf("stride: %i\n", (int)cgf->Stride);

  if (cgf->TilePathCount>0) {

    sz = fwrite(&(cgf->TileStepCount[0]), sizeof(uint64_t), cgf->TileStepCount.size(), ofp);
    if (sz!=cgf->TileStepCount.size()) { return -1; }
    byte_count += (uint64_t)(sizeof(uint64_t)*(cgf->TileStepCount.size()));

    //DEBUG
    //printf("tilestepcount(%i):", (int)cgf->TileStepCount.size());
    //for (i=0; i<cgf->TileStepCount.size(); i++) { printf(" %i", (int)cgf->TileStepCount[i]); }
    //printf("\n");

    sz = fwrite(&(cgf->StrideOffset[0]), sizeof(uint64_t), cgf->StrideOffset.size(), ofp);
    if (sz!=cgf->StrideOffset.size()) { return -1; }
    byte_count += (uint64_t)(sizeof(uint64_t)*(cgf->StrideOffset.size()));

    //DEBUG
    //printf("stride(%i):", (int)cgf->StrideOffset.size());
    //for (i=0; i<cgf->StrideOffset.size(); i++) { printf(" %i", (int)cgf->StrideOffset[i]); }
    //printf("\n");


    sz = fwrite(&(cgf->Loq[0]), sizeof(unsigned char), cgf->Loq.size(), ofp);
    if (sz!=cgf->Loq.size()) { return -1; }
    byte_count += (uint64_t)(cgf->Loq.size());

    //DEBUG
    //printf("loq(%i):", (int)cgf->Loq.size());
    //for (i=0; i<cgf->Loq.size(); i++) { printf(" %02x", cgf->Loq[i]); }
    //printf("\n");


    sz = fwrite(&(cgf->Span[0]), sizeof(unsigned char), cgf->Span.size(), ofp);
    if (sz!=cgf->Span.size()) { return -1; }
    byte_count += (uint64_t)(cgf->Span.size());

    //DEBUG
    //printf("span(%i):", (int)cgf->Span.size());
    //for (i=0; i<cgf->Span.size(); i++) { printf(" %02x", cgf->Span[i]); }
    //printf("\n");


    sz = fwrite(&(cgf->Canon[0]), sizeof(unsigned char), cgf->Canon.size(), ofp);
    if (sz!=cgf->Canon.size()) { return -1; }
    byte_count += (uint64_t)(cgf->Canon.size());

    //DEBUG
    //printf("canon(%i):", (int)cgf->Canon.size());
    //for (i=0; i<cgf->Canon.size(); i++) { printf(" %02x", cgf->Canon[i]); }
    //printf("\n");

    sz = fwrite(&(cgf->CacheOverflow[0]), sizeof(unsigned char), cgf->CacheOverflow.size(), ofp);
    if (sz!=cgf->CacheOverflow.size()) { return -1; }
    byte_count += (uint64_t)(cgf->CacheOverflow.size());

    //DEBUG
    //printf("ovf(%i):", (int)cgf->CacheOverflow.size());
    //for (i=0; i<cgf->CacheOverflow.size(); i++) { printf(" %02x", cgf->CacheOverflow[i]); }
    //printf("\n");



    //--

    sz = fwrite(&(cgf->OverflowOffset[0]), sizeof(uint64_t), cgf->OverflowOffset.size(), ofp);
    if (sz!=cgf->OverflowOffset.size()) { return -1; }
    byte_count += (uint64_t)(cgf->OverflowOffset.size());

    //DEBUG
    //printf("overflowoffset(%i):", (int)cgf->OverflowOffset.size());
    //for (i=0; i<cgf->OverflowOffset.size(); i++) { printf(" %i", (int)cgf->OverflowOffset[i]); }
    //printf("\n");

    //DEBUG
    //printf("overflow(%i):", (int)cgf->Overflow.size());

    if (cgf->Overflow.size()>0) {
      sz = fwrite(&(cgf->Overflow[0]), sizeof(uint16_t), cgf->Overflow.size(), ofp);
      if (sz!=cgf->Overflow.size()) { return -1; }
      byte_count += (uint64_t)(cgf->Overflow.size());

      //DEBUG
      //for (i=0; i<cgf->Overflow.size(); i++) { printf(" %i", (int)cgf->Overflow[i]); }

    }

    //DEBUG
    //printf("\n");

    //--

    sz = fwrite(&(cgf->Overflow64Offset[0]), sizeof(uint64_t), cgf->Overflow64Offset.size(), ofp);
    if (sz!=cgf->Overflow64Offset.size()) { return -1; }
    byte_count += (uint64_t)(cgf->Overflow64Offset.size());

    //DEBUG
    //printf("overflow64offset(%i):", (int)cgf->Overflow64Offset.size());
    //for (i=0; i<cgf->Overflow64Offset.size(); i++) { printf(" %i", (int)cgf->Overflow64Offset[i]); }
    //printf("\n");

    if (cgf->Overflow64.size()>0) {
      sz = fwrite(&(cgf->Overflow64[0]), sizeof(uint64_t), cgf->Overflow64.size(), ofp);
      if (sz!=cgf->Overflow64.size()) { return -1; }
      byte_count += (uint64_t)(cgf->Overflow64.size());

      //DEBUG
      //for (i=0; i<cgf->Overflow64.size(); i++) { printf(" %i", (int)cgf->Overflow64[i]); }

    }

    //DEBUG
    //printf("\n");


    //--

    sz = fwrite(&(cgf->TilePathStructOffset[0]), sizeof(uint64_t), cgf->TilePathStructOffset.size(), ofp);
    if (sz!=cgf->TilePathStructOffset.size()) { return -1; }
    byte_count += (uint64_t)(cgf->TilePathStructOffset.size());

    //DEBUG
    //printf("tilepathstruct(%i):", (int)cgf->TilePathStructOffset.size());
    //for (i=0; i<cgf->TilePathStructOffset.size(); i++) { printf(" %i", (int)cgf->TilePathStructOffset[i]); }
    //printf("\n");


    for (ii=0; ii<cgf->TilePathStructOffset.size(); ii++) {
      tilepath = &(cgf->TilePath[ii]);

      //printf("  name(%i): %s\n", (int)tilepath->Name.size(), tilepath->Name.c_str());

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


    }

  }


  fclose(ofp);

	return byte_count;
}


/*
void cgf_output_band_format(cgf_t *cgf, int tilepath_idx) {
  int i, j, k;
  int ntile, n_q, n_r;

  unsigned char loq, hiq;

  std::vector<int> variant_v[2];
  std::vector< std::vector<int> > noc_v[2];

  std::vector<int> v;

  ntile = (int)cgf->TileStepCount[tilepath_idx];

  variant_v[0].resize(ntile);
  variant_v[1].resize(ntile);
  noc_v[0].resize(ntile);
  noc_v[1].resize(ntile);

  for (i=0; i<ntile; i++) {
    variant_v[0] = -1;
    variant_v[1] = -1;
    noc_v[0] = v;
    noc_v[1] = v;
  }

  for (i=0; i<ntile; i++) {
    n_q = i/8;
    n_r = i%8;

    hiq = ~(cgf->Loq[n_q]);

    if ( (hiq&(1<<n_r)) && (cgf->Canon[n_q]&(1<<n_r)) ) {
      variant_v[0][i] = 0;
      variant_v[1][i] = 0;
      continue;
    }

  }




  for (i=0; i<n_ovf; i+=3) {

    tilestep = (int)ovf[i];
    int vara = (int)ovf[i+1];
    int varb = (int)ovf[i+2];

    if (vara >= OVF16_MAX) { vara = -1; }
    if (varb >= OVF16_MAX) { varb = -1; }

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
*/

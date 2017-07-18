#include "cgf4.hpp"

int cgf_sanity(cgf_t *cgf) {
  int i, j, k;
  uint64_t n64, ii, jj;
  int stride = 0;

  for (i=0; i<8; i++) {
    if (cgf->Magic[i] != CGF_MAGIC[i]) { return -1; }
  }

  if (cgf->CGFVersion.size()==0) { return -2; }
  if (cgf->LibraryVersion.size()==0) { return -3; }

  n64 = cgf->TilePathCount;
  stride = (int)cgf->Stride;

  if (stride<=0) { return -4; }

  if (n64==0) { return 0; }

  if (cgf->Loq.size() != cgf->Span.size()) { return -10; }
  if (cgf->Loq.size() != cgf->Canon.size()) { return -11; }
  if (cgf->Loq.size() != cgf->CacheOverflow.size()) { return -12; }

  for (ii=0; ii<n64; ii++) {
    if ((ii>0) && (cgf->StrideOffset[ii-1] >= cgf->StrideOffset[ii])) { return -5; }
    if (cgf->StrideOffset[ii] >= cgf->Loq.size()) { return -6; }
    if (cgf->StrideOffset[ii] >= cgf->Span.size()) { return -7; }
    if (cgf->StrideOffset[ii] >= cgf->Canon.size()) { return -8; }
    if (cgf->StrideOffset[ii] >= (stride*cgf->CacheOverflow.size())) { return -9; }
  }

  return 0;
}

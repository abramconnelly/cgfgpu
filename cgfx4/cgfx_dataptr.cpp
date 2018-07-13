#include <stdio.h>
#include <string.h>

#include "cgfx_dataptr.hpp"


DataPtr::DataPtr()
{
	type = 0;
	max = 0;
	num = 0;
	size = 0;
	stride = 0;
	cpu = 0;
	gpu = 0;
}

void Allocator::Create(DataPtr& p, int max, int stride)
{
	p.max = max;
	p.num = 0;
	p.stride = stride;
	p.size = max*stride;

	if (p.cpu != 0) { free(p.cpu); p.cpu = 0; }
	if (p.gpu != 0) { cuMemFree(p.gpu); p.gpu = 0; }

	p.cpu = (char*) malloc(p.size);
	cuMemAlloc( &p.gpu, p.size );
}

void Allocator::Destroy(DataPtr& p)
{
	if (p.cpu != 0) { free(p.cpu); p.cpu = 0; }
	if (p.gpu != 0) { cuMemFree(p.gpu); p.gpu = 0; }
	p.max = 0;
	p.num = 0;
	p.size = 0;
	p.stride = 0;
}
int Allocator::Alloc(DataPtr& p)
{
	int id = p.num;
	p.num++;
	return id;
}

void Allocator::Zero(DataPtr& p)
{
	memset (p.cpu, 0, p.size);
}

void Allocator::Commit(DataPtr& p)
{
	cuMemcpyHtoD(p.gpu, p.cpu, p.size);
}

void Allocator::Retrieve(DataPtr& p)
{
	cuMemcpyDtoH(p.cpu, p.gpu, p.size);
}


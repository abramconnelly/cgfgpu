
#include <stdio.h>

#include "app_perf.h"
#include "cgfx.hpp"


#include <time.h>

int cgfx_main(int argc, char **argv)
{
	clock_t t;
	int match = 0, tot = 0;

	int start = 0x000;
	int end = 0x35E;

	/*printf("Run cgf4\n");

	cgf_t *cgfa, *cgfb;
	int f;
	FILE* ifp;
	if ((ifp = fopen("D:\\Codes\\cgfgpu_new\\data\\HG00472-GS000016706-ASM.cgfv4", "rb")) == NULL) { printf("ERROR: Cannot read\n"); return 0; }
	cgfa = cgf_read(ifp);
	fclose(ifp);

	if ((ifp = fopen("D:\\Codes\\cgfgpu_new\\data\\HG02147-GS000016341-ASM.cgfv4", "rb")) == NULL) { printf("ERROR: Cannot read\n"); return 0; }
	cgfb = cgf_read(ifp);
	fclose(ifp);

	printf("Concordance (cgf4)\n");
	t = clock();
	cgf_hiq_concordance(&match, &tot, cgfa, cgfb, start, 0, end, 34, NULL);
	t = clock() - t;
	printf("Match: %d of %d, time: %f ms\n\n", match, tot, (float)t * 1000.0f / CLOCKS_PER_SEC);

	int blk_start, blk_end;
	cgf_get_block_start_end(cgfa, blk_start, blk_end, start, 0, end, 34);
	uchar* match_list = (uchar*)malloc((blk_end - blk_start + 1) * sizeof(uchar));
	uchar* tot_list = (uchar*)malloc((blk_end - blk_start + 1) * sizeof(uchar));

	t = clock();
	cgf_hiq_concordance_no_overflow(&match, &tot, match_list, tot_list, cgfa, cgfb, start, 0, end, 34, NULL);
	t = clock() - t;

	printf("Concordance (cgf4), no overflow\n");
	for (int n = 0; n < 25; n++) {
		printf("%02d ", tot_list[n] == 0 ? 0 : ((int)match_list[n] * 99) / (int)tot_list[n]);
	}
	printf("\nMatch: %d of %d, blocks: %d, time: %f ms\n\n", match, tot, blk_end - blk_start, (float)t * 1000.0f / CLOCKS_PER_SEC);
	*/

	printf("Run cgfx\n");

	CGFX_t cgfx;

	int num = 8;

	cgfx.Initialize(201, 10669088);
	cgfx.Reserve( num );

	cgfx.LoadCgf4("data/HG00472-GS000016706-ASM.cgfv4");
	cgfx.LoadCgf4("data/HG02147-GS000016341-ASM.cgfv4");
	cgfx.LoadCgf4("data/hu4BF398-GS01175-DNA_A06.cgfv4");
	cgfx.LoadCgf4("data/hu34D5B9-GS01173-DNA_C07.cgfv4");
	cgfx.LoadCgf4("data/hu34D5B9-GS01670-DNA_E02.cgfv4");
	cgfx.LoadCgf4("data/hu826751-GS03052-DNA_B01.cgfv4");
	cgfx.LoadCgf4("data/huFFB09D-GS01669-DNA_D04.cgfv4");
	cgfx.LoadCgf4("data/NA12776-GS000016396-ASM.cgfv4");

/*	cgfx.LoadCgf4("data/HG00472-GS000016706-ASM.cgfv4");
	cgfx.LoadCgf4("data/HG02147-GS000016341-ASM.cgfv4");
	cgfx.LoadCgf4("data/hu4BF398-GS01175-DNA_A06.cgfv4");
	cgfx.LoadCgf4("data/hu34D5B9-GS01173-DNA_C07.cgfv4");
	cgfx.LoadCgf4("data/hu34D5B9-GS01670-DNA_E02.cgfv4");
	cgfx.LoadCgf4("data/hu826751-GS03052-DNA_B01.cgfv4");
	cgfx.LoadCgf4("data/huFFB09D-GS01669-DNA_D04.cgfv4"); 
	cgfx.LoadCgf4("data/NA12776-GS000016396-ASM.cgfv4"); */


	cgfx.SetConcordanceRange(start, 0, end, 34);
	cgfx.SetConcordanceOutput( num );

	/*PERF_PUSH("CPU");
	for (int i = 0; i < num; i++)
		for (int j = 0; j < num; j++)
		    cgfx.Concordance(i, j);
	PERF_POP(); */

	FILE* fp = fopen("concord.raw", "wb");

	fwrite(&num, 1, sizeof(int), fp);		// number of groups

	int i = 0;

	PERF_PUSH("GPU");
	cgfx.ConcordanceGPU( i );
	cuCtxSynchronize();
	PERF_POP();
		
	PERF_PUSH("CPU");
	for (int j = 0; j < num; j++) {
		cgfx.Concordance(i, j, j);

/*		int sz = cgfx.getMatchSize();
		fwrite(&num, 1, sizeof(int), fp);	// number of concordances in group
		fwrite(&sz, 1, sizeof(int), fp);	// number of entries
		for (int j = 0; j < num; j++) {
			fwrite( cgfx.getMatchList(j), 1, sz, fp);
			fwrite( cgfx.getTotalList(j), 1, sz, fp); 
		} */
	}
	PERF_POP();

	fclose(fp);


	/*int amt;
	PERF_PUSH("Full");
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < i; j++) {
			if (i != j) {
				amt = cgfx.Concordance(i, j);
				printf("%d, %d: %d\n", i, j, amt);
			}
		}
		printf("\n");
	}
	PERF_POP();*/

	printf("Done.\n");

	return 1;
}

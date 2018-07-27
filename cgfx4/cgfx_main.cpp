
#include <stdio.h>

#include "app_perf.h"
#include "cgfx.hpp"

#include <time.h>

int cgfx_main(int argc, char **argv)
{
	clock_t t;
	int match = 0, tot = 0;

	//int start = 0x2F8;
	//int end = 0x2F9;  // 0x35E = 862 last tilepath (inclusive);
	
	int start = 1;
	int end = 862; 

	int startview = 0;

	cgf_t *cgfa, *cgfb;
	int f;
	FILE* ifp;
	if ((ifp = fopen("data/HG00472-GS000016706-ASM.cgfv4", "rb")) == NULL) { printf("ERROR: Cannot read\n"); return 0; }
	cgfa = cgf_read_hiq(ifp);
	fclose(ifp);

	if ((ifp = fopen("data/HG02147-GS000016341-ASM.cgfv4", "rb")) == NULL) { printf("ERROR: Cannot read\n"); return 0; }
	cgfb = cgf_read_hiq(ifp);
	fclose(ifp);

	int blk_start, blk_end;
	cgf_get_block_start_end(cgfa, blk_start, blk_end, start, 0, end, 0);
	int num_block = blk_end-blk_start+1;

	uchar* match_list_a = (uchar*)	malloc( num_block * sizeof(uchar));
	uchar* match_list_b = (uchar*)	malloc( num_block * sizeof(uchar));
	uchar* tot_list = (uchar*)		malloc( num_block * sizeof(uchar));

	// HIQ Concordance, no overflow
	t = clock();
	cgf_hiq_concordance_no_overflow(&match, &tot, match_list_b, tot_list, cgfa, cgfb, start, 0, end, 0, NULL);
	t = clock() - t;
	printf("Match: %3.2f%% (%d / %d), t: %4.1f ms ---- cgf4 cpu, no overflow\n", match*100.0f / tot, match, tot, (float)t * 1000.0f / CLOCKS_PER_SEC);
	for (int n = startview; n < startview + 25; n++) {
		printf("%02d ", match_list_b[n]);
	}
	printf("\n\n");

	// HIQ Concordance, OVERFLOW
	t = clock();
	cgf_hiq_concordance(&match, &tot, cgfa, cgfb, start, 0, end, 0, NULL);
	t = clock() - t;
	printf("Match: %3.2f%% (%d / %d), t: %4.1f ms ---- cgf4 cpu, OVERFLOW\n", match*100.0f / tot, match, tot, (float)t * 1000.0f / CLOCKS_PER_SEC);
	

	int k = 0;
	int num_block_in_path = cgfa->TileStepCount[start] / 32;
	for (int s=0; s < num_block_in_path; s++ ) {
		cgf_hiq_concordance(&match, &tot, cgfa, cgfb, start, s*32, start, s* 32 + 31, NULL);
		match_list_a[k] = match;
		k++;
	}
	for (int n = startview; n < startview + 25; n++) {
		printf("%02d ", match_list_a[n]);
	}
	printf("\n\n");
 
	// Measure overflow
	int over_amt = 0;
	int over_cnt = 0;
	int over_min=-1, over_max=-1;
	int diff;
	for (int n = 0; n <  num_block_in_path; n++) {
		if (match_list_a[n] != match_list_b[n]) {
			diff = fabs(match_list_b[n] - match_list_a[n]);
			if (diff < over_min || over_min==-1) over_min = diff;
			if (diff > over_max || over_max== -1) over_max = diff;
			over_amt += diff;
			over_cnt++;
			//printf("%d %d\n", n, diff);
		}
	}
	printf("num blk/path: %d\n", num_block_in_path);
	printf("overflow cnt: %d (%3.1f%%)\n", over_cnt, float(over_cnt)*100 / num_block_in_path);
	printf("overflow amt: %f (%3.1f%%)\n", float(over_amt) / over_cnt, float(over_amt)*100/(over_cnt*32));
	printf("overflow min: %d\n", over_min);
	printf("overflow max: %d\n", over_max);

	CGFX_t cgfx;


	cgfx.Initialize(201, 10669088);

	cgfx.LoadCgf4("data/HG00472-GS000016706-ASM.cgfv4");
	cgfx.LoadCgf4("data/HG02147-GS000016341-ASM.cgfv4");
	/*cgfx.LoadCgf4("data/hu4BF398-GS01175-DNA_A06.cgfv4");
	cgfx.LoadCgf4("data/hu34D5B9-GS01173-DNA_C07.cgfv4");
	cgfx.LoadCgf4("data/hu34D5B9-GS01670-DNA_E02.cgfv4");
	cgfx.LoadCgf4("data/hu826751-GS03052-DNA_B01.cgfv4");
	cgfx.LoadCgf4("data/huFFB09D-GS01669-DNA_D04.cgfv4");
	cgfx.LoadCgf4("data/NA12776-GS000016396-ASM.cgfv4");*/

	cgfx.SetConcordanceRange(start, 0, end, 0);
	cgfx.SetConcordanceOutput( cgfx.getNumGenome() );

	int i = 0;

	PERF_PUSH("CPU");
	cgfx.SetVerbose(true);
	cgfx.SetStartView(startview);
	cgfx.Concordance(0, 1, 1);
	//for (int j = 0; j < num_genomes; j++)
	//	cgfx.Concordance(i, j, j);
	PERF_POP();

	PERF_PUSH("GPU");
	cgfx.ConcordanceGPU(i);
	PERF_POP();	

	// write single group
	FILE* fp = fopen("concord.raw", "wb");
	int grps = 1;
	fwrite(&grps, 1, sizeof(int), fp);			// number of groups
	int sz = cgfx.getMatchSize();
	int num = cgfx.getNumGenome();
	fwrite(&num, 1, sizeof(int), fp);	// number of concordances in group
	fwrite(&sz, 1, sizeof(int), fp);			// number of entries
	for (int j = 0; j < num; j++) {
		fwrite(cgfx.getMatchList(j), 1, sz, fp);
		fwrite(cgfx.getTotalList(j), 1, sz, fp);
	}

	fclose(fp);


	printf("Done.\n");

	return 1;
}
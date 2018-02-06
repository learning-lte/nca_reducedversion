#include "merge and reduce.h"
//#include <time.h>

#define SIZE 40

int main()
{
	u32 i, j, k, tt,p[4];
	u8 z[SIZE], NFSR[40], LFSR[40], ns[24], prefix[2], pattern[4];
	struct ncastate state2[80];
	
	struct isd2 *isd2_table[4][4];               // This is the pointer array to construct all the isd tables with the 4 ksds being the indices
	struct isd3 *isd2_table2[4][4];
	struct isd4 *isd2_table3[16][16];
	struct initialdata indata;
	struct datalist *list;
	struct datalist2 *list2, *finallist;
	double counter, Pr[4][4],ttc;
	u16  *try1, *try2;
	struct inner_link *node1,*pn;
	struct datalist2_inner *list2a;
	
	struct initialdata2 *initiallist;
	

	srand((unsigned)time(NULL));
	rc4_setup();

	list = NULL;
	try1 = NULL;
	try2 = NULL;
	initiallist = NULL;
	finallist = NULL;
	list2a = NULL;
	list2 = NULL;
	node1 = NULL;
	pn = NULL;
	
	memset(p, 0, 4);
	memset(z, 0, SIZE);
	memset(ns, 0, 24);
	memset(pattern, 0, 4);
	memset(prefix, 0, 2);
	memset(NFSR, 0, 40);
	memset(LFSR, 0, 40);
	memset(&indata, 0, sizeof(struct initialdata));
	memset(state2, 0, sizeof(struct ncastate) * 80);
	

	for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)isd2_table[i][j] = NULL;
	for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)isd2_table2[i][j] = NULL;
	
	
	for (i = 0; i < 16; i++)
	for (j = 0; j < 16; j++)isd2_table3[i][j] = NULL;
	for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)Pr[i][j] = 0.0;


	list2 = (struct datalist2 *)calloc(1, sizeof(struct datalist2));
	if (list2 == NULL){
		printf("error in calloc in main function\n");
	}

	list2a = (struct datalist2_inner *)calloc(1, sizeof(struct datalist2_inner));
	if (list2a == NULL){
		printf("error in calloc in main function\n");
	}

	
	finallist = (struct datalist2 *)calloc(1, sizeof(struct datalist2));
	if (finallist == NULL){
		printf("error in calloc in main function\n");
	}

	node1 = (struct inner_link *)calloc(4, sizeof(struct inner_link));
	if (node1 == NULL){
		printf("error in calloc in node1 in main function\n");
	}

	try1 = (u16 *)calloc(256, sizeof(u16));
	if (try1 == NULL){
		printf("error in calloc in try1 in main function\n");
	}

	
	Precomputation(ns, isd2_table);
	Isd2Pr_diversity(isd2_table, Pr);
	output_Pr(Pr);
	printf("-------------------------------------------------------------------\n");



	initiallist = (struct initialdata2 *)calloc(1, sizeof(struct initialdata2));
	if (initiallist == NULL){
		printf("error in calloc in main function\n");
	}

	Attacktarget_preparation2(initiallist, LFSR, NFSR, list2, 0);

	try1 = FirstSet_preparation(try1, list2, isd2_table, LFSR, NFSR);
//	Multiple_reduction(try1, try2, list2, isd2_table, LFSR, NFSR);
	

	Test_constant2(list2, isd2_table, LFSR, NFSR);
	Test_average2(list2, isd2_table, LFSR, NFSR);

	Merge_normal(initiallist, list2, isd2_table, LFSR, NFSR);
	Merge_normal_combin2(initiallist, list2, isd2_table, LFSR, NFSR);

	
	free(try1);
	free(list);
	free(list2);
	free(list2a);
	
	free(node1);

	return(0);
}



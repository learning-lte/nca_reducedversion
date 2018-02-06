#include "near collision attack.h"

void Precomputation(u8 *ns, struct isd2 *isd2_table[4][4])
{
	u8 i, j, k;
	u8 keystr[2];

	memset(keystr, 0, 2);

	for (i = 0; i < 4; i++)
	for (k = 0; k < 4; k++)
	{
		for (j = 0; j < 2; j++)keystr[j] = (k >> j) & 0x1;
		isd2_table[i][k] = Ksd2isd2_nca_arbitrary(ns, isd2_table[i][k], i, keystr);
		Sort_isd2(isd2_table[i][k], 0, length_isd2(isd2_table[i][k]));
		printf("The KSD is %u and one of the prefix is %u\n", i, k);
		printf("The length is %u\n", length_isd2(isd2_table[i][k]));
//		Output_isd2(isd2_table[i][k]);
		printf("***************************************\n");
//		TotalISD_BSW_arbitrary(LFSR, NFSR, isd2_table, i, keystr);
	}
	printf("Finished pre-computation\n");

	/*for (i = 0; i < 4; i++)
	for (k = 0; k < 4; k++){
	printf("The KSD is %u and one of the prefix is %u: ", i, k);
	printf("The length is %u\n", length_isd2(isd2_table[i][k]));
	printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	}
	printf("-------------------------------------------------------------------\n");*/
}


void Inner_Precomputation(u8 *ns, struct isd3 *isd2_table[4][4], u8 pattern[2])
{
	u8 i, j, k;
	u8 keystr[2];

	memset(keystr, 0, 2);

	for (i = 0; i < 4; i++)
	for (k = 0; k < 4; k++)
	{
		printf("i=%u,k=%u\n", i, k);
		for (j = 0; j < 2; j++)keystr[j] = (k >> j) & 0x1;
		isd2_table[i][k] = Ksd2isd2_innernca_arbitrary(ns, isd2_table[i][k], i, keystr,pattern);
//		Sort_isd3(isd2_table[i][k], 0, length_isd3(isd2_table[i][k]));
//		printf("The KSD is %u and one of the prefix is %u\n", i, k);
//		printf("The length is %u\n", length_isd3(isd2_table[i][k]));
		Output_isd3(isd2_table[i][k]);
		printf("***************************************\n");
	}
//	printf("Finished inner pre-computation\n");

	
}

void Inner_Precomputation6(u8 *ns, struct isd4 *isd2_table[16][16], u8 pattern[4])
{
	u8 i, j, k;
	u8 keystr[4];

	memset(keystr, 0, 4);

	for (i = 0; i < 16; i++)
	for (k = 0; k < 16; k++)
	{
		printf("i=%u,k=%u\n", i, k);
		for (j = 0; j < 4; j++)keystr[j] = (k >> j) & 0x1;
//		isd2_table[i][k] = Ksd2isd6_innernca_arbitrary(ns, isd2_table[i][k], i, keystr, pattern);
		isd2_table[i][k] = Ksd2isd6_innernca_arbitrary4(ns, isd2_table[i][k], i, keystr, pattern);
//		Sort_isd3(isd2_table[i][k], 0, length_isd3(isd2_table[i][k]));
//		printf("The KSD is %u and one of the prefix is %u\n", i, k);
//		printf("The length is %u\n", length_isd3(isd2_table[i][k]));
//		Output_isd3(isd2_table[i][k]);
		Output_isd4(isd2_table[i][k]);
		printf("------------\n");
	}
	//	printf("Finished inner pre-computation\n");


}

void Attacktarget_preparation_inner(struct initialdata2 *initiallist, u8 *LFSR, u8 *NFSR,  u8 begin_point,u8 pattern[4])
{
	initiallist = Targetkeystream_generation_inner(LFSR, NFSR, initiallist,pattern);
	output_initialdata2(initiallist);

//	Prepare_ncadata(initiallist, list, begin_point);
}

void Attacktarget_preparation_pr(struct initialdata2 *initiallist, u8 *LFSR, u8 *NFSR, struct datalist2 *list, u8 begin_point, u8 prefix[2])
{
	initiallist = Targetkeystream_generation2(LFSR, NFSR, initiallist, prefix);
	output_initialdata2(initiallist);

	Prepare_ncadata(initiallist, list, begin_point);
}


void Attacktarget_preparation2(struct initialdata2 *initiallist, u8 *LFSR, u8 *NFSR, struct datalist2 *list, u8 begin_point)
{
	initiallist = Targetkeystream_generation(LFSR, NFSR, initiallist);
	output_initialdata2(initiallist);

	Prepare_ncadata(initiallist, list, begin_point);
}

u16 *FirstSet_preparation(u16 *try1, struct datalist2 *list2, struct isd2 *isd2_table[4][4], u8 *LFSR, u8 *NFSR)
{
	u8 in1;
	double counter;

	in1 = 0;
	counter = 0.0;
	do{
		try1 = (u16 *)calloc(256, sizeof(u16));
		if (try1 == NULL){
			printf("error in calloc in try1 in FirstSet_preparation\n");
		}

		try1 = ISrecovery_nca_Selfcontained(LFSR, NFSR, list2, isd2_table, try1);
		Self_contained2_ncacounting(try1);

		in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
			printf("Corret1\n");
			printf("++++++++++++++++++++++++++++\n");
			counter++;
		}
		else free(try1);
	} while (in1 == 0);

	return(try1);
}

u16 *FirstSet_preparation_inner(u16 *try1, struct datalist2 *list2, struct isd3 *isd2_table[4][4], u8 *LFSR, u8 *NFSR,u8 pattern[2])
{
	u8 in1;
	double counter;

	in1 = 0;
	counter = 0.0;
	do{
		try1 = (u16 *)calloc(4096, sizeof(u16));
		if (try1 == NULL){
			printf("error in calloc in try1 in FirstSet_preparation_inner\n");
		}

		try1 = ISrecovery_innernca_Selfcontained(LFSR, NFSR, list2, isd2_table, try1,pattern);
		Self_contained2_ncacounting(try1);

		in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
			printf("Corret1\n");
			printf("++++++++++++++++++++++++++++\n");
			counter++;
		}
		else free(try1);
	} while (in1 == 0);

	return(try1);
}

u16 *FirstSet_preparation_inner4(u16 *try1, struct datalist2_inner *list2, struct isd4 *isd2_table[16][16], u8 *LFSR, u8 *NFSR, u8 pattern[4])
{
	u8 in1;
	double counter;

	in1 = 0;
	counter = 0.0;
	do{
		try1 = (u16 *)calloc(256, sizeof(u16));
		if (try1 == NULL){
			printf("error in calloc in try1 in FirstSet_preparation_inner\n");
		}

		try1 = ISrecovery_innernca_Selfcontained4(LFSR, NFSR, list2, isd2_table, try1, pattern);
		Self_contained2_ncacounting4(try1);

		in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
			printf("Corret1\n");
			printf("++++++++++++++++++++++++++++\n");
			counter++;
		}
		else free(try1);
	} while (in1 == 0);

	return(try1);
}

void Multiple_reduction(u16 *try1, u16 *try2, struct datalist2 *list2, struct isd2 *isd2_table[4][4], u8 *LFSR, u8 *NFSR)
{
	clock_t start, finish;
	double duration, speed;
	u32 i,j,k;
	u8 in1;
	
	try2 = (u16 *)calloc(256, sizeof(u16));
	if (try2 == NULL){
		printf("error in calloc in try2 in main function\n");
	}

	start = clock();
	for (i = 0; i < 0x1; i++){

		start = clock();
		for (k = 0; k < 120; k++){
			for (j = 0; j < 256; j++)try2[j] = 0;
			try2 = ISrecovery_nca_Selfcontained(LFSR, NFSR, list2, isd2_table, try2);
//			Self_contained2_intersection2(try1, try2);
			printf("k=%u: ",k);
			Self_contained2_intersection_nca(try1, try2);

		}

		Self_contained2_ncacounting(try1);
		in1 = 0;
		in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
			printf("Success!\n");
		}
	}
	finish = clock();

	duration = ((double)(finish - start)) / CLOCKS_PER_SEC;
	printf("The routine takes %f seconds\n", duration);

	free(try2);
}

void Multiple_reduction_inner(u16 *try1, u16 *try2, struct datalist2 *list2, struct isd3 *isd2_table[4][4], u8 *LFSR, u8 *NFSR,u8 pattern[2])
{
	clock_t start, finish;
	double duration, speed;
	u32 i, j, k;
	u8 in1;

	try2 = (u16 *)calloc(4096, sizeof(u16));
	if (try2 == NULL){
		printf("error in calloc in try2 in main function\n");
	}

	start = clock();
	for (i = 0; i < 0x1; i++){

		start = clock();
		for (k = 0; k < 100; k++){
			for (j = 0; j < 4096; j++)try2[j] = 0;
			try2 = ISrecovery_innernca_Selfcontained(LFSR, NFSR, list2, isd2_table, try2, pattern);
//			Self_contained2_intersection2(try1, try2);
			printf("k=%u: ", k);
			Self_contained2_intersection_nca(try1, try2);

		}

		Self_contained2_ncacounting(try1);
		in1 = 0;
		in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
			printf("Success!\n");
		}
	}
	finish = clock();

	duration = ((double)(finish - start)) / CLOCKS_PER_SEC;
	printf("The routine takes %f seconds\n", duration);

	free(try2);
}

void Multiple_reduction_inner4(u16 *try1, u16 *try2, struct datalist2_inner *list2, struct isd4 *isd2_table[16][16], u8 *LFSR, u8 *NFSR, u8 pattern[4])
{
	clock_t start, finish;
	double duration, speed;
	u32 i, j, k;
	u8 in1;

	try2 = (u16 *)calloc(256, sizeof(u16));
	if (try2 == NULL){
		printf("error in calloc in try2 in Multiple_reduction_inner4\n");
	}

	start = clock();
	for (i = 0; i < 0x1; i++){

		start = clock();
		for (k = 0; k < 100; k++){
			for (j = 0; j < 256; j++)try2[j] = 0;
			try2 = ISrecovery_innernca_Selfcontained4(LFSR, NFSR, list2, isd2_table, try2, pattern);
//			Self_contained2_intersection2(try1, try2);
			printf("k=%u: ", k);
			Self_contained2_intersection_nca4(try1, try2);

		}

		Self_contained2_ncacounting4(try1);
		in1 = 0;
		in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
			printf("Success!\n");
		}
	}
	finish = clock();

	duration = ((double)(finish - start)) / CLOCKS_PER_SEC;
	printf("The routine takes %f seconds\n", duration);

	free(try2);
}

void Test_constant(struct datalist2 *list2, struct isd2 *isd2_table[4][4], u8 *LFSR, u8 *NFSR)
{
	u8 in1, *try1;
	double counter;
	u32 j,k;
	
	in1 = 0;
	k = 0;
	try1 = NULL;
	counter = 0.0;

	list2 = SetA_generation_Self_contained2(LFSR, NFSR, list2);
	printf("The searched element is:");
	for (j = 0; j < 23; j++)printf("%x", (list2 + 0)->ISnca[j]);
	printf("\n");
//	output_datalist2(list2);
	try1 = (u8 *)calloc(8388608, sizeof(u8));
	if (try1 == NULL){
		printf("error in calloc in try1 in Test_constant\n");
	}
	do{
		for (j = 0; j < 8388608; j++)try1[j] = 0;
		try1 = SetB_generation_Self_contained4(LFSR, NFSR, list2, isd2_table, try1);
		Self_contained2_counting(try1);
		in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
			printf("Corret1\n");
			counter++;
		}

		k++;
	} while (k < 8000);
	
	printf("The ratio that the correct key will be in the self_contained list is %f\n", counter / 8000);

	free(try1);
}

void Test_constant2(struct datalist2 *list2, struct isd2 *isd2_table[4][4], u8 *LFSR, u8 *NFSR)
{
	u8 in1;
	u16 *try1;
	double counter;
	u32 j, k;

	in1 = 0;
	k = 0;
	try1 = NULL;
	counter = 0.0;

	try1 = (u16 *)calloc(256, sizeof(u16));
	if (try1 == NULL){
		printf("error in calloc in try1 in Test_constant2\n");
	}
	do{
		for (j = 0; j < 256; j++)try1[j] = 0;
		try1 = ISrecovery_nca_Selfcontained(LFSR, NFSR, list2, isd2_table, try1);
		Self_contained2_ncacounting(try1);
		in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
//			printf("Corret1\n");
			counter++;
		}

		k++;
	} while (k < 16777216);

	printf("The ratio that the correct key will be in the self_contained list is %f\n", counter / 16777216);

	free(try1);
}

void Test_constant2_inner(struct datalist2 *list2, struct isd3 *isd2_table[4][4], u8 *LFSR, u8 *NFSR, u8 pattern[2])
{
	u8 in1;
	u16 *try1;
	double counter;
	u32 j, t = 0, k, kk[4];

	in1 = 0;
	k = 0;
	try1 = NULL;
	counter = 0.0;
	for (j = 0; j < 4;j++)kk[j] = 0;

	try1 = (u16 *)calloc(4096, sizeof(u16));
	if (try1 == NULL){
		printf("error in calloc in try1 in Test_constant2\n");
	}
	do{
		for (j = 0; j < 4096; j++)try1[j] = 0;
		try1 = ISrecovery_innernca_Selfcontained(LFSR, NFSR, list2, isd2_table, try1,pattern);
//		Self_contained2_ncacounting(try1);
		t = Self_contained2_ncacounting2(try1);
		kk[t-1]++;
//		printf("t=%u\n", t);
		in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
//			printf("Corret1\n");
			counter++;
		}

		k++;
	} while (k < 4194304);

	printf("The ratio that the correct key will be in the self_contained list is %f\n", counter / 4194304);

	for (j = 0; j < 4; j++)printf("%f\n",((double)kk[j])/4194304);

	free(try1);
}

void Test_constant2_inner4(struct datalist2_inner *list2, struct isd4 *isd2_table[16][16], u8 *LFSR, u8 *NFSR, u8 pattern[4])
{
	u8 in1;
	u16 *try1;
	double counter;
	u32 j, t = 0, k, kk[16];

	in1 = 0;
	k = 0;
	try1 = NULL;
	counter = 0.0;
	for (j = 0; j < 16; j++)kk[j] = 0;

	try1 = (u16 *)calloc(256, sizeof(u16));
	if (try1 == NULL){
		printf("error in calloc in try1 in Test_constant2_inner4\n");
	}
	do{
		for (j = 0; j < 256; j++)try1[j] = 0;
		try1 = ISrecovery_innernca_Selfcontained4(LFSR, NFSR, list2, isd2_table, try1, pattern);
//		Self_contained2_ncacounting(try1);
		t = Self_contained2_ncacounting4a(try1);
		kk[t - 1]++;
//		printf("t=%u\n", t);
		in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
//			printf("Corret1\n");
			counter++;
		}

		k++;
	} while (k < 4194304);

	printf("The ratio that the correct key will be in the self_contained list is %f\n", counter / 4194304);

	for (j = 0; j < 16; j++)printf("%f\n", ((double)kk[j]) / 4194304);

	free(try1);
}

void Test_average2(struct datalist2 *list2, struct isd2 *isd2_table[4][4], u8 *LFSR, u8 *NFSR)
{
	u8 in1;
	u16 *try1;
	double counter;
	u32 j, k,xx=0;

	in1 = 0;
	k = 0;
	try1 = NULL;
	counter = 0.0;

	try1 = (u16 *)calloc(256, sizeof(u16));
	if (try1 == NULL){
		printf("error in calloc in try1 in Test_constant2\n");
	}
	do{
		for (j = 0; j < 256; j++)try1[j] = 0;
		try1 = ISrecovery_nca_Selfcontained(LFSR, NFSR, list2, isd2_table, try1);
		xx += Self_contained2_ncacounting2(try1);
		/*in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
			printf("Corret1\n");
			counter++;
		}*/

		k++;
	} while (k < 16777216);

	printf("The ratio that the correct key will be in the self_contained list is %f\n", (double)xx / 16777216);

	free(try1);
}

void Test_average2_inner(struct datalist2 *list2, struct isd3 *isd2_table[4][4], u8 *LFSR, u8 *NFSR, u8 pattern[2])
{
	u8 in1;
	u16 *try1;
	double counter;
	u32 j, k, xx = 0;

	in1 = 0;
	k = 0;
	try1 = NULL;
	counter = 0.0;

	try1 = (u16 *)calloc(4096, sizeof(u16));
	if (try1 == NULL){
		printf("error in calloc in try1 in Test_constant2\n");
	}
	do{
		for (j = 0; j < 4096; j++)try1[j] = 0;
		try1 = ISrecovery_innernca_Selfcontained(LFSR, NFSR, list2, isd2_table, try1,pattern);
		xx += Self_contained2_ncacounting2(try1);
		/*in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
		printf("Corret1\n");
		counter++;
		}*/

		k++;
	} while (k < 4194304);

	printf("The ratio that the correct key will be in the self_contained list is %f\n", (double)xx / 4194304);

	free(try1);
}

void Test_average2_inner4(struct datalist2_inner *list2, struct isd4 *isd2_table[16][16], u8 *LFSR, u8 *NFSR, u8 pattern[4])
{
	u8 in1;
	u16 *try1;
	double counter;
	u32 j, k, xx = 0;

	in1 = 0;
	k = 0;
	try1 = NULL;
	counter = 0.0;

	try1 = (u16 *)calloc(256, sizeof(u16));
	if (try1 == NULL){
		printf("error in calloc in try1 in Test_average2_inner4\n");
	}
	do{
		for (j = 0; j < 256; j++)try1[j] = 0;
		try1 = ISrecovery_innernca_Selfcontained4(LFSR, NFSR, list2, isd2_table, try1, pattern);
		xx += Self_contained2_ncacounting4a(try1);
		/*in1 = check_candidate1(try1, list2->ISnca_sum);
		if (in1 != 0){
		printf("Corret1\n");
		counter++;
		}*/

		k++;
	} while (k < 4194304);

	printf("The ratio that the correct key will be in the self_contained list is %f\n", (double)xx / 4194304);

	free(try1);
}

void clear_node(struct inner_link *node4)
{
	u8 i;

	for (i = 0; i < 4; i++){
		(node4 + i)->x0 = 0xff;
		(node4 + i)->n63 = 0xff;
		(node4 + i)->x1 = 0xff;
		(node4 + i)->n64 = 0xff;
	}
}


struct inner_link *Twelve2four(struct datalist2 *list2, struct isd3 *isd2_table[4][4], u8 *LFSR, u8 *NFSR, u8 pattern[2], struct inner_link *node4, u8 begin_point)
{
	u8 in1;
	u16 *try1;
	u32 j, i, k, k1, xx = 0, t = 0, position[4],r=0;
	double counter;
	struct inner_link *p1,*p2;
	
	try1 = NULL;
	p1 = NULL;
	p2 = NULL;
	counter = 0.0;
	in1 = 0;
	memset(position,0,4);
	
	try1 = (u16 *)calloc(4096, sizeof(u16));
	if (try1 == NULL){
		printf("error in calloc in try1 in Twelve2four\n");
	}
	for (j = 0; j < 4096; j++)try1[j] = 0;
	try1 = ISrecovery_innernca_Selfcontained(LFSR, NFSR, list2, isd2_table, try1, pattern);

	xx = Self_contained2_ncacounting2(try1);
	printf("xx=%u\n",xx);
	/*in1 = check_candidate1(try1, list2->ISnca_sum);
	if (in1 != 0){
		printf("Corret1\n");
		counter++;
	}*/
	for (j = 0; j < 4; j++)position[j] = 4097;
	k = 0;
	for (j = 0; j < 4096;j++)
	if (try1[j] != 0){
		position[k] = j;
		k++;
	}

	printf("The positions are:");
	/*for (j = 0; j < 4;j++)
	if (position[j] != 4097)printf("j=%u,%u\n",j,position[j]);*/
	for (j = 0; j < 4; j++)
	if (position[j] != 4097)printf("%u,", position[j]);
	printf("\n");
	
	/*r = 0;
	if (begin_point == 0){
		k = 0;
		for (j = 0; j < 4; j++)
		if (position[j] != 4097){
			t = position[j] & 0xc03;
			(node4 + k)->x0 = t & 0x1;
			(node4 + k)->n63 = (t >> 10) & 0x1;
			(node4 + k)->x1 = (t >> 1) & 0x1;
			(node4 + k)->n64 = (t >> 11) & 0x1;
			(node4 + k)->left = NULL;
			(node4 + k)->right = NULL;
			k++;
		}
		r++;
		p2 = node4;
	}
	else{
		p1 = (struct inner_link *)calloc(4, sizeof(struct inner_link));
		if (p1 == NULL){
			printf("error in calloc in p1 in Twelve2four\n");
		}
		k1 = 0;
		for (j = 0; j < 4;j++)
		if (position[j] != 4097){
			
			t = position[j] & 0xc03;
			
			for (k = 0; k < 4;k++)
			if ((((t >> 1) & 0x1) == (p2 + k)->x1) && (((t >> 11) & 0x1) == (p2 + k)->n64)){
				(p1 + k1)->x0 = t & 0x1;
				(p1 + k1)->n63 = (t >> 10) & 0x1;
				(p1 + k1)->x1 = (t >> 1) & 0x1;
				(p1 + k1)->n64 = (t >> 11) & 0x1;
				(p1 + k1)->left = NULL;
				(p1 + k1)->right = NULL;
				(p2 + k)->left = p1;
			}
			k1++;
		}
		p2 = p1;
	}*/
	
	p2 = node4;
	if (begin_point == 0){
		k = 0;
		for (j = 0; j < 4; j++)
		if (position[j] != 4097){
			t = position[j] & 0xc03;
			(node4 + k)->x0 = t & 0x1;
			(node4 + k)->n63 = (t >> 10) & 0x1;
			(node4 + k)->x1 = (t >> 1) & 0x1;
			(node4 + k)->n64 = (t >> 11) & 0x1;
			(node4 + k)->left = NULL;
			(node4 + k)->right = NULL;
			k++;
		}
		
	}
	else{
		p1 = (struct inner_link *)calloc(4, sizeof(struct inner_link));
		if (p1 == NULL){
			printf("error in calloc in p1 in Twelve2four\n");
		}
		clear_node(p1);
		k = 0;
		for (j = 0; j < 4; j++)
		if (position[j] != 4097){
			t = position[j] & 0xc03;
			(p1 + k)->x0 = t & 0x1;
			(p1 + k)->n63 = (t >> 10) & 0x1;
			(p1 + k)->x1 = (t >> 1) & 0x1;
			(p1 + k)->n64 = (t >> 11) & 0x1;
			(p1 + k)->left = NULL;
			(p1 + k)->right = NULL;
			k++;
		}
	}

	// the following routine is used to connect the nodes
	if (begin_point != 0){
		for (i = 0; i < 4; i++){           // the outtest loop for p1
			if (((p1 + i)->x0) != 0xff){
				k = 0;
				while( ((((p1 + i)->x0) != ((p2 + k)->x1)) || (((p1 + i)->n63) != ((p2 + k)->n64))) && (k<4) )k++;

				if (k < 4){
					if ((p2 + k)->left == NULL){
						(p2 + k)->left = p1 + i;
						p2 = p1;
					}
					else{
						(p2 + k)->right = p1 + i;
						p2 = p1;
					}
					
				}
			}
			else break;
			
		}
	}

	free(try1);
	return(node4);
}

struct inner_link *Twelve2four_a(struct inner_link *node4,u32 position)
{
	u32 j,t=0;

	j = position;
	t = j & 0xc03;
	node4->x0 = t & 0x1;
	node4->n63 = (t >> 10) & 0x1;
	node4->x1 = (t >> 1) & 0x1;
	node4->n64 = (t >> 11) & 0x1;
	node4->left = NULL;
	node4->right = NULL;

	return(node4);
}



void Output_node(struct inner_link *node4)
{
	printf("x0=%x,", node4->x0);
	printf("n63=%x,", node4->n63);
	printf("x1=%x,", node4->x1);
	printf("n64=%x,", node4->n64);
	printf("\n");
	
}

void preorder_output(struct inner_link *node4)
{
	if (node4){
		Output_node(node4);
		preorder_output(node4->left);
		preorder_output(node4->right);
	}
	
}

void Merge_normal(struct initialdata2 *initiallist, struct datalist2 *list2, struct isd2 *isd2_table[4][4], u8 *LFSR, u8 *NFSR)
{
	u16 *try1, *try2, *result[50000];
	u8 in1, in11;
	u32 i, j, k, tt,xx=0;
	struct mask_extraction *rx[2], *pp1;
	double counter, c1;
	clock_t start, finish;
	double duration, speed;

	try1 = NULL;
	try2 = NULL;
	for (i = 0; i < 50000; i++)result[i] = NULL;
	for (i = 0; i < 2; i++)rx[i] = NULL;
	pp1 = NULL;
	
	for (i = 0; i < 2; i++)rx[i] = NULL;
	
	c1 = 0.0;
	for (tt = 0; tt < 50000; tt++){
	
//		Prepare_data(initiallist, list2, tt);

		result[tt] = (u16 *)calloc(256, sizeof(u16));
		if (result[tt] == NULL){
			printf("error in calloc in result in Merge_normal\n");
		}
		in11 = 0;
		for (i = 0; i < 4; i++){

			
			try1 = (u16 *)calloc(256, sizeof(u16));                
			if (try1 == NULL){
				printf("error in calloc in try1 in Merge_normal\n");
			}

			try1 = ISrecovery_nca_Selfcontained(LFSR, NFSR, list2, isd2_table, try1);
			Self_contained2_ncacounting(try1);


			try2 = (u16 *)calloc(256, sizeof(u16));
			if (try2 == NULL){
				printf("error in calloc in try2 in Merge_normal\n");
			}

			
			for (k = 0; k < 10; k++){
				for (j = 0; j < 256; j++)try2[j] = 0;
				try2 = ISrecovery_nca_Selfcontained(LFSR, NFSR, list2, isd2_table, try2);
//				printf("k=%u: ", k);
				Self_contained2_intersection_nca(try1, try2);
			}


			Self_contained2_ncacounting(try1);
			in1 = 0;
			in1 = check_candidate1(try1, list2->ISnca_sum);
			if (in1 != 0){
//				printf("Success!\n");
				in11 = 1;
			}
//			printf("----------------------------------------------------------\n");
			free(try2);
			Merge_nca(try1, result[tt]);
			free(try1);
		}

		if (in11 != 0)c1 = c1 + 1.0;
//		printf("tt=%u: ", tt);
		xx += Counting_merge_nca(result[tt]);

		
	}

	printf("The ratio is %f\n",c1/50000.0);
	printf("The average is %f\n", (double)xx/50000.0);
	for (i = 0; i < 50000; i++)
	if (result[i] != NULL)free(result[i]);
	//for (i = 0; i < 2; i++){
	//	rx[i] = (struct mask_extraction *)calloc(1048576, sizeof(struct mask_extraction));
	//	if (rx[i] == NULL){
	//		printf("error in calloc in rx in Merge_normal\n");
	//	}
	//}
	//Extraction1(result[0], 0x000556aa, rx[0]);
	//Extraction2(result[1], 0x0002ab55, rx[1]);

	//Counting_mask(rx[0]);
	//Counting_mask(rx[1]);

	//Sort_mask(rx[0]);
	//Sort_mask(rx[1]);

	//output_mask(rx[0]);
	//output_mask(rx[1]);

	//printf("^^^^^^^^^^^^^^^^^^^^^^^^\n");
	//Counting_mv(rx[0]);
	//Counting_mv(rx[1]);
	//printf("^^^^^^^^^^^^^^^^^^^^^^^^\n");

	//Reduction_union(rx[0], rx[1]);


	//printf("To 1\n");

	//pp1 = (struct mask_extraction *)calloc(8388608, sizeof(struct mask_extraction));                    // there are no more than 2^23 candidates left
	//if (pp1 == NULL){
	//	printf("error in calloc in Merge_normal\n");
	//}
	//pp1 = Reduction_union2(rx[0], rx[1], pp1);
	//Counting_union(pp1);

	//printf("Begins\n");
	//Sort_mask_nc3(pp1);
//	output_mask2(pp1);

//	for (i = 0; i < 2; i++)free(rx[i]);
	
//	free(pp1);
	
}


void Merge_normal_combin(struct initialdata2 *initiallist, struct datalist2 *list2, struct isd2 *isd2_table[4][4], u8 *LFSR, u8 *NFSR)
{
	u16 *try1, *try2, *result[50000], *try[19];
	u8 in1, in11;
	u32 i, j, k, tt, xx = 0;
	struct mask_extraction *rx[2], *pp1;
	double counter, c1;
	clock_t start, finish;
	double duration, speed;

	try1 = NULL;
	try2 = NULL;
	for (i = 0; i < 50000; i++)result[i] = NULL;
	for (i = 0; i < 2; i++)rx[i] = NULL;
	pp1 = NULL;

	for (i = 0; i < 19; i++)try[i] = NULL;

	c1 = 0.0;
	for (tt = 0; tt < 50000; tt++){

//		Prepare_data(initiallist, list2, tt);

		result[tt] = (u16 *)calloc(4096, sizeof(u16));
		if (result[tt] == NULL){
			printf("error in calloc in result in Merge_normal\n");
		}
		in11 = 0;
		for (i = 0; i < 6; i++){
			
			for (k = 0; k < 19; k++){
				
				try[k] = (u16 *)calloc(4096, sizeof(u16));
				if (try[k] == NULL){
					printf("error in calloc in try[%u] in Merge_normal\n", k);
				}
				for (j = 0; j < 4096; j++)try[k][j] = 0;
				try[k] = ISrecovery_nca_Selfcontained(LFSR, NFSR, list2, isd2_table, try[k]);
//				Self_contained2_ncacounting(try[k]);
			}
				
			for (j = 0; j < 18; j++){
				for (k = j+1; k < 19; k++){
					printf("(j=%u,k=%u): ", j, k);
					Self_contained2_intersection_nca(try[j], try[k]);
				}
			}

			for (k = 0; k < 18; k++)Self_contained2_intersection_nca(try[0], try[k]);

			Self_contained2_ncacounting(try[0]);
			in1 = 0;
			in1 = check_candidate1(try[0], list2->ISnca_sum);
			if (in1 != 0){
				printf("Success!\n");
				in11 = 1;
			}
			printf("----------------------------------------------------------\n");
			for (k = 1; k < 19;k++)free(try[k]);
			Merge_nca(try[0], result[tt]);
			free(try[0]);
		}

		if (in11 != 0)c1 = c1 + 1.0;
		printf("tt=%u: ", tt);
		xx += Counting_merge_nca(result[tt]);


	}

	printf("The ratio is %f\n", c1 / 50000.0);
	printf("The average is %f\n", (double)xx / 50000.0);
	for (i = 0; i < 50000; i++)
	if (result[i] != NULL)free(result[i]);
	
}

void Merge_normal_combin2(struct initialdata2 *initiallist, struct datalist2 *list2, struct isd2 *isd2_table[4][4], u8 *LFSR, u8 *NFSR)
{
	u16 *try1, *try2[19], *result[50000];
	u8 in1, in11;
	u32 i, j, k, tt, xx = 0;
	struct mask_extraction *rx[2], *pp1;
	double counter, c1;
	clock_t start, finish;
	double duration, speed;

	try1 = NULL;
	
	for (i = 0; i < 50000; i++)result[i] = NULL;
	for (i = 0; i < 2; i++)rx[i] = NULL;
	pp1 = NULL;

	for (i = 0; i < 19; i++)try2[i] = NULL;

	c1 = 0.0;
	for (tt = 0; tt < 50000; tt++){

//		Prepare_data(initiallist, list2, tt);

		result[tt] = (u16 *)calloc(256, sizeof(u16));
		if (result[tt] == NULL){
			printf("error in calloc in result in Merge_normal_combin2\n");
		}
		in11 = 0;
		for (i = 0; i < 4; i++){


			try1 = (u16 *)calloc(256, sizeof(u16));
			if (try1 == NULL){
				printf("error in calloc in try1 in Merge_normal_combin2\n");
			}

			try1 = ISrecovery_nca_Selfcontained(LFSR, NFSR, list2, isd2_table, try1);
//			Self_contained2_ncacounting(try1);

			for (k = 0; k < 10; k++){
				try2[k] = (u16 *)calloc(256, sizeof(u16));
				if (try2[k] == NULL){
					printf("error in calloc in try2[%u] in Merge_normal_combin2\n",k);
				}
			}


			for (k = 0; k < 10; k++){
				for (j = 0; j < 256; j++)try2[k][j] = 0;
				try2[k] = ISrecovery_nca_Selfcontained(LFSR, NFSR, list2, isd2_table, try2[k]);
//				printf("k=%u: \n", k);
				Self_contained2_intersection_nca(try1, try2[k]);
			}

			for (k = 1; k < 10; k++){
				Self_contained2_intersection_nca(try2[0], try2[k]);
			}
			Self_contained2_ncacounting(try1);
//			printf("start\n");
			Self_contained2_intersection_nca(try1, try2[0]);
//			printf("close\n");
			Self_contained2_ncacounting(try1);
			in1 = 0;
			in1 = check_candidate1(try1, list2->ISnca_sum);
			if (in1 != 0){
//				printf("Success!\n");
				in11 = 1;
			}
//			printf("----------------------------------------------------------\n");
			for (k = 0; k < 10;k++)free(try2[k]);
			Merge_nca(try1, result[tt]);
			free(try1);
		}

		if (in11 != 0)c1 = c1 + 1.0;
//		printf("tt=%u: ", tt);
		xx += Counting_merge_nca(result[tt]);


	}

	printf("The ratio is %f\n", c1 / 50000.0);
	printf("The average is %f\n", (double)xx / 50000.0);
	for (i = 0; i < 50000; i++)
	if (result[i] != NULL)free(result[i]);
}
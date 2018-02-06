#include "cipher.h"

#define SamplingSIZE 4096
#define SamSIZE 4096
#define ListSIZE 61148                       // 2^15.9
#define NCAC 3                               // 2^6/Sum[Binomial[6,i],{i,0,2}]
#define INCAC 2                               // 2^2/Sum[Binomial[2,i],{i,0,1}]
#define VL 4
#define INCAC4 3 

struct ncastate                                 // This is the extracted state to form the nca collision 
{
	u8 state;                               // the position between 0 to 79, 0-39 is the NFSR state, 40-79 is the LFSR state 
	u8 flag;                                // flag =1 means the nca state, otherwise meaning the non-nca state
};



struct isd                                 // to construct the link of ksd2isd
{
	u8 position[3];
	struct isd *next;
};

struct isd2                                 // to construct the link of ksd2isd with the proportion of the isd
{
	u8 position[2];
	double fc;
	struct isd2 *next;
};

struct isd3                                 // to construct the link of ksd2isd with the proportion of the isd
{
	u8 position[1];
	double fc;
	struct isd3 *next;
};

struct isd4                                 // to construct the link of ksd2isd with the proportion of the isd
{
	u8 position[2];
	double fc;
	struct isd4 *next;
};

struct datalist
{
	u8 prefix[2];
	u8 ISnca[8];
	double ISnca_sum;
	u32 position;
};

struct datalist2
{
	u8 prefix[2];
	u8 ISnca[8];
	double ISnca_sum;
	u32 position;
	struct datalist2 *next;
};

struct datalist2_inner
{
	u8 prefix[4];
	u8 ISnca[24];
	double ISnca_sum;
	u32 position;
	struct datalist2_inner *next;
};

struct initialdata
{
	u8 lfsr[40];
	u8 nfsr[40];
//	u8 pattern[0];
};

struct inner_link
{
	u8 x0;
	u8 n63;
	u8 x1;
	u8 n64;
	struct inner_link *left;
	struct inner_link *right;
};

struct initialdata2
{
	u8 lfsr[40];
	u8 nfsr[40];
	u8 keystream[100];
};

struct mask_extraction
{
	u8 bvalue;
	u32 mv;
};

struct mask_extraction2
{
	u8 bvalue;
	u32 mv;
	u32 mv_new;
};

int structcmp(const void *a, const void *b)
{
	return (*(struct datalist *)a).ISnca_sum > (*(struct datalist *)b).ISnca_sum ? 1 : -1;
}


int structcmp2(const void *a, const void *b)
{
	return (*(struct datalist2 *)a).ISnca_sum > (*(struct datalist2 *)b).ISnca_sum ? 1 : -1;
}

int structcmp3(const void *a, const void *b)
{
	return (*(struct mask_extraction *)a).mv > (*(struct mask_extraction *)b).mv ? 1 : -1;
}

int structcmp4(const void *a, const void *b)
{
	return (*(struct mask_extraction2 *)a).mv_new > (*(struct mask_extraction2 *)b).mv_new ? 1 : -1;
}


void Extractstate(struct ncastate state2[160], u8 LFSR[80], u8 *NFSR)       // the extracted state for the first 2 keystream bits of the original grain v1
{
	u32 i;

	memset(state2, 0, sizeof(struct ncastate) * 160);
	/*printf("The extracted state after initialization: \n");
	for (i = 0; i < 80; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 80; i < 160; i++)printf("%x ", state2[i].state);
	printf("\n");   */                                                  // initialize the state2 array

	for (i = 0; i < 80; i++)(state2[i]).state = NFSR[i];
	for (i = 80; i < 160; i++)(state2[i]).state = LFSR[i - 80];

	for (i = 1; i <= 5; i++)(state2[i]).flag = 1;                   // the extracted NFSR part
	for (i = 10; i <= 11; i++)(state2[i]).flag = 1;
	for (i = 31; i <= 32; i++)(state2[i]).flag = 1;
	for (i = 43; i <= 44; i++)(state2[i]).flag = 1;
	for (i = 56; i <= 57; i++)(state2[i]).flag = 1;
	for (i = 63; i <= 64; i++)(state2[i]).flag = 1;

	for (i = 80 + 3; i <= 80 + 4; i++)(state2[i]).flag = 1;            // the extracted LFSR part
	for (i = 80 + 25; i <= 80 + 26; i++)(state2[i]).flag = 1;
	for (i = 80 + 46; i <= 80 + 47; i++)(state2[i]).flag = 1;
	for (i = 80 + 64; i <= 80 + 65; i++)(state2[i]).flag = 1;

	/*printf("The extracted internal state after assigning: \n");
	for (i = 0; i < 80; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 0; i < 80; i++)printf("%x ", state2[i].flag);
	printf("\n");
	for (i = 80; i < 160; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 80; i < 160; i++)printf("%x ", state2[i].flag);
	printf("\n");
	printf("\n");*/
}


void Extractstate_BSW(struct ncastate state2[80], u8 LFSR[40], u8 *NFSR)           // the extracted state for the first 2 keystream bits of the original grain v1
{
	u32 i;

	memset(state2, 0, sizeof(struct ncastate) * 80);
	/*printf("The extracted state after initialization: \n");
	for (i = 0; i < 40; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 40; i < 80; i++)printf("%x ", state2[i].state);
	printf("\n");   */                                                  // initialize the state2 array

	for (i = 0; i < 31; i++)(state2[i]).state = NFSR[i];
	for (i = 33; i < 80; i++)(state2[i]).state = NFSR[i];
	for (i = 80; i < 160; i++)(state2[i]).state = LFSR[i - 80];

	for (i = 31; i <= 32; i++)state2[i].state = NFSR[1 + i - 31] ^ NFSR[2 + i - 31] ^ NFSR[4 + i - 31] ^ NFSR[10 + i - 31] ^ NFSR[43 + i - 31] ^ NFSR[56 + i - 31] ^ LFSR[25 + i - 31] ^ NFSR[63 + i - 31] ^ (LFSR[3 + i - 31] & LFSR[64 + i - 31])
		^ (LFSR[46 + i - 31] & LFSR[64 + i - 31]) ^ (LFSR[64 + i - 31] & NFSR[63 + i - 31]) ^ (LFSR[3 + i - 31] & LFSR[25 + i - 31] & LFSR[46 + i - 31]) ^ (LFSR[3 + i - 31] & LFSR[46 + i - 31] & LFSR[64 + i - 31])
		^ (LFSR[3 + i - 31] & LFSR[46 + i - 31] & NFSR[63 + i - 31]) ^ (LFSR[25 + i - 31] & LFSR[46 + i - 31] & NFSR[63 + i - 31]) ^ (LFSR[46 + i - 31] & LFSR[64 + i - 31] & NFSR[63 + i - 31]);

	for (i = 1; i <= 5; i++)(state2[i]).flag = 1;                   // the extracted NFSR part
	for (i = 10; i <= 11; i++)(state2[i]).flag = 1;
	for (i = 31; i <= 32; i++)(state2[i]).flag = 1;
	for (i = 43; i <= 44; i++)(state2[i]).flag = 1;
	for (i = 56; i <= 57; i++)(state2[i]).flag = 1;
	for (i = 63; i <= 64; i++)(state2[i]).flag = 1;

	for (i = 80 + 3; i <= 80 + 4; i++)(state2[i]).flag = 1;            // the extracted LFSR part
	for (i = 80 + 25; i <= 80 + 26; i++)(state2[i]).flag = 1;
	for (i = 80 + 46; i <= 80 + 47; i++)(state2[i]).flag = 1;
	for (i = 80 + 64; i <= 80 + 65; i++)(state2[i]).flag = 1;

	/*printf("The extracted internal state after assigning: \n");
	for (i = 0; i < 40; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 0; i < 40; i++)printf("%x ", state2[i].flag);
	printf("\n");
	for (i = 40; i < 80; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 40; i < 80; i++)printf("%x ", state2[i].flag);
	printf("\n");
	printf("\n");*/
}


void Extractstate_BSW2(struct ncastate state2[160], u8 LFSR[80], u8 *NFSR)      // the extracted state for the first 2 keystream bits (0x03) of grain v1
{
	u32 i;

	memset(state2, 0, sizeof(struct ncastate) * 160);
	/*printf("The extracted state after initialization: \n");
	for (i = 0; i < 40; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 40; i < 80; i++)printf("%x ", state2[i].state);
	printf("\n");   */                                                  // initialize the state2 array

	for (i = 0; i < 31; i++)(state2[i]).state = NFSR[i];
	for (i = 33; i < 80; i++)(state2[i]).state = NFSR[i];
	for (i = 80; i < 160; i++)(state2[i]).state = LFSR[i - 80];

	for (i = 31; i <= 32; i++)state2[i].state = 0x1 ^ NFSR[1 + i - 31] ^ NFSR[2 + i - 31] ^ NFSR[4 + i - 31] ^ NFSR[10 + i - 31] ^ NFSR[43 + i - 31] ^ NFSR[56 + i - 31] ^ LFSR[25 + i - 31] ^ NFSR[63 + i - 31] ^ (LFSR[3 + i - 31] & LFSR[64 + i - 31])
		^ (LFSR[46 + i - 31] & LFSR[64 + i - 31]) ^ (LFSR[64 + i - 31] & NFSR[63 + i - 31]) ^ (LFSR[3 + i - 31] & LFSR[25 + i - 31] & LFSR[46 + i - 31]) ^ (LFSR[3 + i - 31] & LFSR[46 + i - 31] & LFSR[64 + i - 31])
		^ (LFSR[3 + i - 31] & LFSR[46 + i - 31] & NFSR[63 + i - 31]) ^ (LFSR[25 + i - 31] & LFSR[46 + i - 31] & NFSR[63 + i - 31]) ^ (LFSR[46 + i - 31] & LFSR[64 + i - 31] & NFSR[63 + i - 31]);

	for (i = 1; i <= 5; i++)(state2[i]).flag = 1;                   // the extracted NFSR part
	for (i = 10; i <= 11; i++)(state2[i]).flag = 1;
	for (i = 31; i <= 32; i++)(state2[i]).flag = 1;
	for (i = 43; i <= 44; i++)(state2[i]).flag = 1;
	for (i = 56; i <= 57; i++)(state2[i]).flag = 1;
	for (i = 63; i <= 64; i++)(state2[i]).flag = 1;

	for (i = 80 + 3; i <= 80 + 4; i++)(state2[i]).flag = 1;            // the extracted LFSR part
	for (i = 80 + 25; i <= 80 + 26; i++)(state2[i]).flag = 1;
	for (i = 80 + 46; i <= 80 + 47; i++)(state2[i]).flag = 1;
	for (i = 80 + 64; i <= 80 + 65; i++)(state2[i]).flag = 1;

	/*printf("The extracted internal state after assigning: \n");
	for (i = 0; i < 40; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 0; i < 40; i++)printf("%x ", state2[i].flag);
	printf("\n");
	for (i = 40; i < 80; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 40; i < 80; i++)printf("%x ", state2[i].flag);
	printf("\n");
	printf("\n");*/
}

void Extractstate_BSW_arbitrary(struct ncastate state2[160], u8 LFSR[80], u8 *NFSR, u8 z[2])        // the extracted state for the first 2 keystream bits of the reduced grain v1
{
	u32 i;

	memset(state2, 0, sizeof(struct ncastate) * 160);
	/*printf("The extracted state after initialization: \n");
	for (i = 0; i < 40; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 40; i < 80; i++)printf("%x ", state2[i].state);
	printf("\n");   */                                                  // initialize the state2 array

	for (i = 0; i < 31; i++)(state2[i]).state = NFSR[i];
	for (i = 33; i < 80; i++)(state2[i]).state = NFSR[i];
	for (i = 80; i < 160; i++)(state2[i]).state = LFSR[i - 80];

	for (i = 31; i <= 32; i++)state2[i].state = z[i-31] ^ NFSR[1 + i - 31] ^ NFSR[2 + i - 31] ^ NFSR[4 + i - 31] ^ NFSR[10 + i - 31] ^ NFSR[43 + i - 31] ^ NFSR[56 + i - 31] ^ LFSR[25 + i - 31] ^ NFSR[63 + i - 31] ^ (LFSR[3 + i - 31] & LFSR[64 + i - 31])
		^ (LFSR[46 + i - 31] & LFSR[64 + i - 31]) ^ (LFSR[64 + i - 31] & NFSR[63 + i - 31]) ^ (LFSR[3 + i - 31] & LFSR[25 + i - 31] & LFSR[46 + i - 31]) ^ (LFSR[3 + i - 31] & LFSR[46 + i - 31] & LFSR[64 + i - 31])
		^ (LFSR[3 + i - 31] & LFSR[46 + i - 31] & NFSR[63 + i - 31]) ^ (LFSR[25 + i - 31] & LFSR[46 + i - 31] & NFSR[63 + i - 31]) ^ (LFSR[46 + i - 31] & LFSR[64 + i - 31] & NFSR[63 + i - 31]);

	for (i = 1; i <= 5; i++)(state2[i]).flag = 1;                   // the extracted NFSR part
	for (i = 10; i <= 11; i++)(state2[i]).flag = 1;
	for (i = 31; i <= 32; i++)(state2[i]).flag = 1;
	for (i = 43; i <= 44; i++)(state2[i]).flag = 1;
	for (i = 56; i <= 57; i++)(state2[i]).flag = 1;
	for (i = 63; i <= 64; i++)(state2[i]).flag = 1;

	for (i = 80 + 3; i <= 80 + 4; i++)(state2[i]).flag = 1;            // the extracted LFSR part
	for (i = 80 + 25; i <= 80 + 26; i++)(state2[i]).flag = 1;
	for (i = 80 + 46; i <= 80 + 47; i++)(state2[i]).flag = 1;
	for (i = 80 + 64; i <= 80 + 65; i++)(state2[i]).flag = 1;

	/*printf("The extracted internal state after assigning: \n");
	for (i = 0; i < 40; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 0; i < 40; i++)printf("%x ", state2[i].flag);
	printf("\n");
	for (i = 40; i < 80; i++)printf("%x ", state2[i].state);
	printf("\n");
	for (i = 40; i < 80; i++)printf("%x ", state2[i].flag);
	printf("\n");
	printf("\n");*/
}

void Extract_state(u8 LFSR[80], u8 *NFSR, u8 z[2])
{
	randomIV(LFSR, 80);
	randomIV(NFSR, 80);

	NFSR[31] = z[0] ^ NFSR[1] ^ NFSR[2] ^ NFSR[4] ^ NFSR[10] ^ NFSR[43] ^ NFSR[56] ^ h_function(LFSR[3], LFSR[25], LFSR[46], LFSR[64], NFSR[63]);
	NFSR[32] = z[1] ^ NFSR[2] ^ NFSR[3] ^ NFSR[5] ^ NFSR[11] ^ NFSR[44] ^ NFSR[57] ^ h_function(LFSR[4], LFSR[26], LFSR[47], LFSR[65], NFSR[64]);
}


void Extractnca_arbitrary(u8 ns[8], u8 *LFSR, u8 *NFSR, u8 z[2])
{
	u8 i;

	memset(ns, 0, 8);
	ns[2] = LFSR[1];
	ns[3] = LFSR[2];
	ns[4] = LFSR[21];
	ns[5] = LFSR[22];
	ns[6] = NFSR[23];
	ns[7] = NFSR[24];
	

	ns[0] = z[0] ^ h_function(LFSR[1], LFSR[21], NFSR[23]);
	ns[1] = z[1] ^ h_function(LFSR[2], LFSR[22], NFSR[24]);
}

void Extractnca_arbitrary_inner(u8 ns[12], u8 *LFSR, u8 *NFSR, u8 z[2],u8 pattern[2])
{
	u8 i;

	memset(ns, 0, 12);
	
	LFSR[3] = pattern[0] & 0x1;
	LFSR[4] = pattern[1] & 0x1;
	LFSR[25] = (pattern[0] >> 1) & 0x1;
	LFSR[26] = (pattern[1] >> 1) & 0x1;
	LFSR[46] = (pattern[0] >> 2) & 0x1;
	LFSR[47] = (pattern[1] >> 2) & 0x1;
	LFSR[64] = (pattern[0] >> 3) & 0x1;
	LFSR[65] = (pattern[1] >> 3) & 0x1;
//	printf("to 1a\n");
	ns[2] = LFSR[3];
	ns[3] = LFSR[4];
	ns[4] = LFSR[25];
	ns[5] = LFSR[26];
	ns[6] = LFSR[46];
	ns[7] = LFSR[47];
	ns[8] = LFSR[64];
	ns[9] = LFSR[65];
	ns[10] = NFSR[63];
	ns[11] = NFSR[64];

	ns[0] = z[0] ^ h_function(LFSR[3], LFSR[25], LFSR[46], LFSR[64], NFSR[63]);
	ns[1] = z[1] ^ h_function(LFSR[4], LFSR[26], LFSR[47], LFSR[65], NFSR[64]);
}

void Extractnca_arbitrary_inner4(u8 ns[24], u8 *LFSR, u8 *NFSR, u8 z[4], u8 pattern[4])
{
	u8 i;

	memset(ns, 0, 24);

	LFSR[3] = pattern[0] & 0x1;
	LFSR[4] = pattern[1] & 0x1;
	LFSR[5] = pattern[2] & 0x1;
	LFSR[6] = pattern[3] & 0x1;

	LFSR[25] = (pattern[0] >> 1) & 0x1;
	LFSR[26] = (pattern[1] >> 1) & 0x1;
	LFSR[27] = (pattern[2] >> 1) & 0x1;
	LFSR[28] = (pattern[3] >> 1) & 0x1;
	
	LFSR[46] = (pattern[0] >> 2) & 0x1;
	LFSR[47] = (pattern[1] >> 2) & 0x1;
	LFSR[48] = (pattern[2] >> 2) & 0x1;
	LFSR[49] = (pattern[3] >> 2) & 0x1;
	
	LFSR[64] = (pattern[0] >> 3) & 0x1;
	LFSR[65] = (pattern[1] >> 3) & 0x1;
	LFSR[66] = (pattern[0] >> 3) & 0x1;
	LFSR[67] = (pattern[1] >> 3) & 0x1;
	
	ns[4] = LFSR[3];
	ns[5] = LFSR[4];
	ns[6] = LFSR[5];
	ns[7] = LFSR[6];
	
	ns[8] = LFSR[25];
	ns[9] = LFSR[26];
	ns[10] = LFSR[27];
	ns[11] = LFSR[28];
	
	ns[12] = LFSR[46];
	ns[13] = LFSR[47];
	ns[14] = LFSR[48];
	ns[15] = LFSR[49];

	ns[16] = LFSR[64];
	ns[17] = LFSR[65];
	ns[18] = LFSR[66];
	ns[19] = LFSR[67];

	ns[20] = NFSR[63];
	ns[21] = NFSR[64];
	ns[22] = NFSR[65];
	ns[23] = NFSR[66];

	ns[0] = z[0] ^ h_function(LFSR[3], LFSR[25], LFSR[46], LFSR[64], NFSR[63]);
	ns[1] = z[1] ^ h_function(LFSR[4], LFSR[26], LFSR[47], LFSR[65], NFSR[64]);
	ns[2] = z[2] ^ h_function(LFSR[5], LFSR[27], LFSR[48], LFSR[66], NFSR[65]);
	ns[3] = z[3] ^ h_function(LFSR[6], LFSR[28], LFSR[49], LFSR[67], NFSR[66]);
}

void Extractnca_arbitrary_inner5(u8 ns[2*VL], u8 *LFSR, u8 *NFSR, u8 z[VL], u8 pattern[2])
{
	u8 i;

	memset(ns, 0, 12);

	LFSR[3] = pattern[0] & 0x1;
	LFSR[4] = pattern[1] & 0x1;
	LFSR[25] = (pattern[0] >> 1) & 0x1;
	LFSR[26] = (pattern[1] >> 1) & 0x1;
	LFSR[46] = (pattern[0] >> 2) & 0x1;
	LFSR[47] = (pattern[1] >> 2) & 0x1;
	LFSR[64] = (pattern[0] >> 3) & 0x1;
	LFSR[65] = (pattern[1] >> 3) & 0x1;
	//	printf("to 1a\n");
	ns[2] = LFSR[3];
	ns[3] = LFSR[4];
	ns[4] = LFSR[25];
	ns[5] = LFSR[26];
	ns[6] = LFSR[46];
	ns[7] = LFSR[47];
	ns[8] = LFSR[64];
	ns[9] = LFSR[65];
	ns[10] = NFSR[63];
	ns[11] = NFSR[64];

	ns[0] = z[0] ^ h_function(LFSR[3], LFSR[25], LFSR[46], LFSR[64], NFSR[63]);
	ns[1] = z[1] ^ h_function(LFSR[4], LFSR[26], LFSR[47], LFSR[65], NFSR[64]);
}

void Keystreamgen_nca2(u8 *ks, u8 ns[8])
{
	ks[0] = ns[0] ^ h_function(ns[2], ns[4], ns[6]);
	ks[1] = ns[1] ^ h_function(ns[3], ns[5], ns[7]);
}

void Keystreamgen_nca6(u8 *ks, u8 ns[24])
{
	ks[0] = ns[0] ^ h_function(ns[4], ns[8], ns[12], ns[16], ns[20]);
	ks[1] = ns[1] ^ h_function(ns[5], ns[9], ns[13], ns[17], ns[21]);
	ks[2] = ns[2] ^ h_function(ns[6], ns[10], ns[14], ns[18], ns[22]);
	ks[3] = ns[3] ^ h_function(ns[7], ns[11], ns[15], ns[19], ns[23]);
	
}

void verify_ncabsw(u8 ns[12], u8 *LFSR, u8 *NFSR, u8 keystr[2])
{
	u8 z[2];
	

	randomIV(LFSR, 40);
	randomIV(NFSR, 40);
	Extractnca_arbitrary(ns, LFSR, NFSR, keystr);
	Keystreamgen_nca2(z, ns);
	
	if ( !memcmp(z, keystr, 2) )printf("ok\n");
	
}

void Diff_back(struct ncastate *state2, u8 *LFSR, u8 *NFSR)
{
	u32 i;

	for (i = 0; i < 80; i++)NFSR[i] = (state2 + i)->state;
	for (i = 0; i < 80; i++)LFSR[i] = (state2 + i + 80)->state;
}

u32 length_isd2(struct isd2 *ph)                             // find the length of the link isd2_table[i]
{
	struct isd2 *p;
	u32 n;

	n = 0;
	p = ph;

	if (ph != NULL)
	{
		do
		{
			n++;
			p = p->next;
		} while (p != NULL);
	}

	return(n-1);                                   // return the index of the last element  
}

u32 length_isd3(struct isd3 *ph)                             // find the length of the link isd2_table2[i]
{
	struct isd3 *p;
	u32 n;

	n = 0;
	p = ph;

	if (ph != NULL)
	{
		do
		{
			n++;
			p = p->next;
		} while (p != NULL);
	}

	return(n - 1);                                   // return the index of the last element  
}

void Swap_isd2(struct isd2 *ph, u32 left, u32 right)          // swap the two structures in the link, indexed by left and right position pointers
{
	struct isd2 *p1, *head, *p2,*p3;
	u32 i;

	p1 = NULL;
	p2 = NULL;
	p3 = NULL;
	head = NULL;

	if (left == 0)p1 = ph;
	else{
		for (p1 = ph, i = 0; i < left; i++)p1 = p1->next;
	}

	for (p2 = ph, i = 0; i < right; i++)p2 = p2->next;

	p3 = (struct isd2 *)calloc(1, sizeof(struct isd2));
	if (NULL == p3)
	{
		printf("Error in calloc in Swap_isd2.\n");
		return 0;
	}

	p3->position[0] = p1->position[0];
	p3->position[1] = p1->position[1];
	p3->position[2] = p1->position[2];
//	p3->position[3] = p1->position[3];
//	p3->position[4] = p1->position[4];
	p3->fc = p1->fc;

	p1->position[0] = p2->position[0];
	p1->position[1] = p2->position[1];
	p1->position[2] = p2->position[2];
//	p1->position[3] = p2->position[3];
//	p1->position[4] = p2->position[4];
	p1->fc = p2->fc;

	p2->position[0] = p3->position[0];
	p2->position[1] = p3->position[1];
	p2->position[2] = p3->position[2];
//	p2->position[3] = p3->position[3];
//	p2->position[4] = p3->position[4];
	p2->fc = p3->fc;

	free(p3);                         // new added instruction
}

void Swap_isd3(struct isd3 *ph, u32 left, u32 right)          // swap the two structures in the link, indexed by left and right position pointers
{
	struct isd3 *p1, *head, *p2, *p3;
	u32 i;

	p1 = NULL;
	p2 = NULL;
	p3 = NULL;
	head = NULL;

	if (left == 0)p1 = ph;
	else{
		for (p1 = ph, i = 0; i < left; i++)p1 = p1->next;
	}

	for (p2 = ph, i = 0; i < right; i++)p2 = p2->next;

	p3 = (struct isd3 *)calloc(1, sizeof(struct isd3));
	if (NULL == p3)
	{
		printf("Error in calloc in Swap_isd3.\n");
		return 0;
	}

	p3->position[0] = p1->position[0];
//	p3->position[1] = p1->position[1];
//	p3->position[2] = p1->position[2];
//	p3->position[3] = p1->position[3];
//	p3->position[4] = p1->position[4];
	p3->fc = p1->fc;

	p1->position[0] = p2->position[0];
//	p1->position[1] = p2->position[1];
//	p1->position[2] = p2->position[2];
//	p1->position[3] = p2->position[3];
//	p1->position[4] = p2->position[4];
	p1->fc = p2->fc;

	p2->position[0] = p3->position[0];
//	p2->position[1] = p3->position[1];
//	p2->position[2] = p3->position[2];
//	p2->position[3] = p3->position[3];
//	p2->position[4] = p3->position[4];
	p2->fc = p3->fc;

	free(p3);                         // new added instruction
}

u32 partition(struct isd2 *ph, u32 left, u32 right)
{
	u32 i, pivotIndex, length, storeIndex,r;
	double pivotValue;
	struct isd2 *p1, *p2, *p3;

	p1 = NULL;
	p2 = NULL;
	p3 = NULL;
	length = length_isd2(ph);

		
	pivotIndex = (u32)(floor((left + right) / 2));                 // choose the pivot element
	p1 = ph;
	for (i = 0; i < pivotIndex; i++)p1 = p1->next;
	pivotValue = p1->fc;

	Swap_isd2(ph, pivotIndex, right);
	storeIndex = left;
	p2 = ph;

	for (i = left; i <= right - 1; i++){
		if (i == 0)p2 = ph;
		else for (p2 = ph, r = 0; r < i; r++)p2 = p2->next;

//		if ((p2->fc) < pivotValue){
		if ((p2->fc) > pivotValue){
			Swap_isd2(ph, i, storeIndex);
			storeIndex++;
		}
	}
	Swap_isd2(ph, storeIndex,right);
		
	return(storeIndex);

}

u32 partition3(struct isd3 *ph, u32 left, u32 right)
{
	u32 i, pivotIndex, length, storeIndex, r;
	double pivotValue;
	struct isd3 *p1, *p2, *p3;

	p1 = NULL;
	p2 = NULL;
	p3 = NULL;
	length = length_isd3(ph);


	pivotIndex = (u32)(floor((left + right) / 2));                 // choose the pivot element
	p1 = ph;
	if(p1 != NULL)for (i = 0; i < pivotIndex; i++)p1 = p1->next;
	pivotValue = p1->fc;

	Swap_isd3(ph, pivotIndex, right);
	storeIndex = left;
	p2 = ph;

	for (i = left; i <= right - 1; i++){
		if (i == 0)p2 = ph;
		else for (p2 = ph, r = 0; r < i; r++)p2 = p2->next;

		//		if ((p2->fc) < pivotValue){
		if ((p2->fc) > pivotValue){
			Swap_isd3(ph, i, storeIndex);
			storeIndex++;
		}
	}
	Swap_isd3(ph, storeIndex, right);

	return(storeIndex);

}

void Sort_isd2(struct isd2 *ph, u32 left, u32 right)          // sort the link of isd2_table[i] according to the field of fc
{
	u32 p;

	
	if (left < right){
		p = partition(ph, left, right);
		Sort_isd2(ph, left, p);
		Sort_isd2(ph, p + 1, right);
	}

}

void Sort_isd3(struct isd3 *ph, u32 left, u32 right)          // sort the link of isd2_table[i] according to the field of fc
{
	u32 p;


	if (left < right){
		p = partition3(ph, left, right);
		Sort_isd3(ph, left, p);
		Sort_isd3(ph, p + 1, right);
	}

}


void Restrictedkeystr_generation(u8 *LFSR, u8 *NFSR, struct datalist *list)
{
	u8 z[2],j,k;
	u32 i;
	struct ncastate state2[160];
	double sum;

	memset(z, 0, 2);
	memset(state2, 0, sizeof(struct ncastate) * 160);

	printf("****************************************\n");
	// to ensure the continuous generation of 2-bit keystream prefix
	for (i = 0; i < 2*ListSIZE; i++){
		Extractstate(state2, LFSR, NFSR);
		Diff_back(state2, LFSR, NFSR);
		
		k = 0;
		j = 0;
		while (j < 23){
			if ((state2 + k)->flag == 1){
				(list + i)->ISnca[j] = (state2 + k)->state;
				k++;
				j++;
			}
			else k++;
		}
		
		Keystreamgen(2, z, LFSR, NFSR);
		for (j = 0; j < 2; j++)(list + i)->prefix[j] = z[j];
		for (sum = 0.0, j = 0; j < 23; j++)sum += ((list + i)->ISnca[j]) * pow(2, j);
		(list + i)->ISnca_sum = sum;
		(list + i)->position = i;
	}
}


void RestrictedData_generation(struct datalist *list,struct initialdata *indata)
{
	u32 i;
	u8 LFSR[80], NFSR[80];

//	memset(pattern, 0, 1);
	
//	randomIV(pattern, 1);
	randomIV(LFSR, 80);
	randomIV(NFSR, 80);

//	printf("1--------------------------------------1\n");
	for (i = 0; i < 80;i++)indata->lfsr[i] = LFSR[i];
	for (i = 0; i < 80; i++)indata->nfsr[i] = NFSR[i];
//	for (i = 0; i < 1; i++)indata->pattern[i] = pattern[i];
//	printf("2--------------------------------------2\n");

//	Restrictedkeystr_generation(LFSR, NFSR, pattern, list);
	Restrictedkeystr_generation(LFSR, NFSR, list);
	printf("*--------------------------------------*\n");
	
}


void Restricted_IS_Collision_finding(struct datalist *list)
{
	u8 i1, i2, i3, i4, i5, inter[23];
	u32 i, j;
	double num, sum, *searchp;
	double *copy;

	for (i = 0; i < 23; i++)inter[i] = 0;
	printf("1--------------------------------------1\n");

	copy = NULL;
	copy = (double *)calloc(2 * ListSIZE, sizeof(double));
	if (copy == NULL){
		printf("error in calloc in Restricted_IS_Collision_finding \n");
	}

	qsort(list, 2 * ListSIZE, sizeof(struct datalist), structcmp);

	for (i = 0; i < 2 * ListSIZE; i++)copy[i] = (list + i)->ISnca_sum;

	printf("1*--------------------------------------1*\n");
	for (i = 0; i < 100; i++)printf("%f\n", (list + i)->ISnca_sum);
	printf("2*--------------------------------------2*\n");

	for (num = 0.0, i = 0; i < 2 * ListSIZE; i++){
		for (i1 = 0; i1 < 23; i1++)
		{
			for (j = 0; j < 23; j++)inter[j] = (list + i)->ISnca[j];                    // first flip 1-5 bits in the specified 19 bits nca state
			inter[i1] ^= 0x1;

			for (sum = 0.0, j = 0; j < 23; j++)sum += inter[j] * pow(2, j);
			searchp = (double *)bsearch(&sum, copy, 2 * ListSIZE, sizeof(double), intcmp);

			if (searchp != NULL){
				num++;
				printf("found! position %u and position %u are collisions\n", (list + i)->position, (list + (int)(searchp - copy))->position);
				for (j = 0; j < 23; j++)printf("%x ", (list + i)->ISnca[j]);
				printf("\n");
				for (j = 0; j < 2; j++)printf("%x ", (list + i)->prefix[j]);
				printf("\n");
				for (j = 0; j < 23; j++)printf("%x ", (list + (int)(searchp - copy))->ISnca[j]);
				printf("\n");
				for (j = 0; j < 2; j++)printf("%x ", (list + (int)(searchp - copy))->prefix[j]);
				printf("\n");
				printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
			}
			else continue;
		}
		//		printf("i=%u\n", i);
	}

	printf("1: There are totally %f collisions\n", num);

	for (num = 0.0, i = 0; i < 2 * ListSIZE; i++){
		for (i1 = 0; i1 < 22; i1++)
		for (i2 = i1 + 1; i2 < 23; i2++)
		{
			for (j = 0; j < 23; j++)inter[j] = (list + i)->ISnca[j];                    //  flip 2 bits in the specified 19 bits nca state
			inter[i1] ^= 0x1;
			inter[i2] ^= 0x1;

			for (sum = 0.0, j = 0; j < 23; j++)sum += inter[j] * pow(2, j);
			searchp = (double *)bsearch(&sum, copy, 2 * ListSIZE, sizeof(double), intcmp);

			if (searchp != NULL){
				num++;
				printf("found! position %u and position %u are collisions\n", (list + i)->position, (list + (int)(searchp - copy))->position);
				for (j = 0; j < 23; j++)printf("%x ", (list + i)->ISnca[j]);
				printf("\n");
				for (j = 0; j < 2; j++)printf("%x ", (list + i)->prefix[j]);
				printf("\n");
				for (j = 0; j < 23; j++)printf("%x ", (list + (int)(searchp - copy))->ISnca[j]);
				printf("\n");
				for (j = 0; j < 2; j++)printf("%x ", (list + (int)(searchp - copy))->prefix[j]);
				printf("\n");
				printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
			}
			else continue;
		}
		//		printf("i=%u\n", i);
	}

	printf("2: There are totally %f collisions\n", num);


	for (num = 0.0, i = 0; i < 2 * ListSIZE; i++){
		for (i1 = 0; i1 < 21; i1++)
		for (i2 = i1 + 1; i2 < 22; i2++)
		for (i3 = i2 + 1; i3 < 23; i3++)
		{
			for (j = 0; j < 23; j++)inter[j] = (list + i)->ISnca[j];                    // first flip 3 bits in the specified 19 bits nca state
			inter[i1] ^= 0x1;
			inter[i2] ^= 0x1;
			inter[i3] ^= 0x1;

			for (sum = 0.0, j = 0; j < 23; j++)sum += inter[j] * pow(2, j);
			searchp = (double *)bsearch(&sum, copy, 2 * ListSIZE, sizeof(double), intcmp);

			if (searchp != NULL){
				num++;
				printf("found! position %u and position %u are collisions\n", (list + i)->position, (list + (int)(searchp - copy))->position);
				for (j = 0; j < 23; j++)printf("%x ", (list + i)->ISnca[j]);
				printf("\n");
				for (j = 0; j < 2; j++)printf("%x ", (list + i)->prefix[j]);
				printf("\n");
				for (j = 0; j < 23; j++)printf("%x ", (list + (int)(searchp - copy))->ISnca[j]);
				printf("\n");
				for (j = 0; j < 2; j++)printf("%x ", (list + (int)(searchp - copy))->prefix[j]);
				printf("\n");
				printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
			}
			else continue;
		}
		//		printf("i=%u\n", i);
	}

	for (num = 0.0, i = 0; i < 2 * ListSIZE; i++){
		for (i1 = 0; i1 < 20; i1++)
		for (i2 = i1 + 1; i2 < 21; i2++)
		for (i3 = i2 + 1; i3 < 22; i3++)
		for (i4 = i3 + 1; i4 < 23; i4++)
		{
			for (j = 0; j < 23; j++)inter[j] = (list + i)->ISnca[j];                    // first flip 4 bits in the specified 19 bits nca state
			inter[i1] ^= 0x1;
			inter[i2] ^= 0x1;
			inter[i3] ^= 0x1;
			inter[i4] ^= 0x1;

			for (sum = 0.0, j = 0; j < 23; j++)sum += inter[j] * pow(2, j);
			searchp = (double *)bsearch(&sum, copy, 2 * ListSIZE, sizeof(double), intcmp);

			if (searchp != NULL){
				num++;
				printf("found! position %u and position %u are collisions\n", (list + i)->position, (list + (int)(searchp - copy))->position);
				for (j = 0; j < 23; j++)printf("%x ", (list + i)->ISnca[j]);
				printf("\n");
				for (j = 0; j < 2; j++)printf("%x ", (list + i)->prefix[j]);
				printf("\n");
				for (j = 0; j < 23; j++)printf("%x ", (list + (int)(searchp - copy))->ISnca[j]);
				printf("\n");
				for (j = 0; j < 2; j++)printf("%x ", (list + (int)(searchp - copy))->prefix[j]);
				printf("\n");
				printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
			}
			else continue;
		}
		//		printf("i=%u\n", i);
	}

	for (num = 0.0, i = 0; i < 2 * ListSIZE; i++){
		for (i1 = 0; i1 < 19; i1++)
		for (i2 = i1 + 1; i2 < 20; i2++)
		for (i3 = i2 + 1; i3 < 21; i3++)
		for (i4 = i3 + 1; i4 < 22; i4++)
		for (i5 = i4 + 1; i5 < 23; i5++)
		{
			for (j = 0; j < 23; j++)inter[j] = (list + i)->ISnca[j];           // first flip 5 bits in the specified 19 bits nca state
			inter[i1] ^= 0x1;
			inter[i2] ^= 0x1;
			inter[i3] ^= 0x1;
			inter[i4] ^= 0x1;
			inter[i5] ^= 0x1;

			for (sum = 0.0, j = 0; j < 23; j++)sum += inter[j] * pow(2, j);
			searchp = (double *)bsearch(&sum, copy, 2 * ListSIZE, sizeof(double), intcmp);

			if (searchp != NULL){
				num++;
				printf("found! position %u and position %u are collisions\n", (list + i)->position, (list + (int)(searchp - copy))->position);
				for (j = 0; j < 23; j++)printf("%x ", (list + i)->ISnca[j]);
				printf("\n");
				for (j = 0; j < 2; j++)printf("%x ", (list + i)->prefix[j]);
				printf("\n");
				for (j = 0; j < 23; j++)printf("%x ", (list + (int)(searchp - copy))->ISnca[j]);
				printf("\n");
				for (j = 0; j < 2; j++)printf("%x ", (list + (int)(searchp - copy))->prefix[j]);
				printf("\n");
				printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
			}
			else continue;
		}
		//		printf("i=%u\n", i);
	}


	printf("3: There are totally %f collisions\n", num);

}




u32 Computing_dif(struct ncastate *state2, u8 *LFSR, u8 *NFSR, u8 i1, u8 i2, u8 i3,u8 i4,u8 i5)
{
	u8 z[2], z1[2], lfsr1[40], nfsr1[40];
	u8 j;
	double pattern;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(lfsr1, 0, 40);
	memset(nfsr1, 0, 40);
//	pattern = 0;

	randomIV(NFSR, 40);
	randomIV(LFSR, 40);
	for (j = 0; j < 40; j++)nfsr1[j] = NFSR[j];
	for (j = 0; j < 40; j++)lfsr1[j] = LFSR[j];
	Extractstate(state2, LFSR, NFSR);
	Keystreamgen(2, z, LFSR, NFSR);

	state2[i1].state ^= 0x1;
	state2[i2].state ^= 0x1;
	state2[i3].state ^= 0x1;
	state2[i4].state ^= 0x1;
	state2[i5].state ^= 0x1;

	Diff_back(state2, lfsr1, nfsr1);
	Keystreamgen(2, z1, lfsr1, nfsr1);

	for (pattern = 0.0, j = 0; j < 2; j++)pattern += (z[j] ^ z1[j]) * pow(2, j);
	return((u32)pattern);

}

void verify_BSW(struct ncastate *state2, u8 *LFSR, u8 *NFSR)
{
	u8 z[2];
	u8 i;

	randomIV(LFSR, 40);
	randomIV(NFSR, 40);
	Extractstate_BSW(state2, LFSR, NFSR);
	Diff_back(state2, LFSR, NFSR);
	Keystreamgen(2, z, LFSR, NFSR);
	for (i = 0; i < 2; i++)printf("%x", z[i]);
	printf("\n");
}

void verify_BSW2(struct ncastate *state2, u8 *LFSR, u8 *NFSR)
{
	u8 z[2];
	u8 i;

	randomIV(LFSR, 40);
	randomIV(NFSR, 40);
	Extractstate_BSW2(state2, LFSR, NFSR);
	Diff_back(state2, LFSR, NFSR);
	Keystreamgen(2, z, LFSR, NFSR);
	for (i = 0; i < 2; i++)printf("%x", z[i]);
	printf("\n");
}



u32 Computing_difBSW_arbitray(struct ncastate *state2, u8 *LFSR, u8 *NFSR, u8 i1, u8 i2, u8 i3,u8 i4,u8 i5,u8 keystr[2])
{
	u8 z[2], z1[2], lfsr1[80], nfsr1[80], zero[2];
	u8 j;
	double pattern;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(zero, 0, 2);
	memset(lfsr1, 0, 80);
	memset(nfsr1, 0, 80);
//	pattern = 0;

	randomIV(NFSR, 80);
	randomIV(LFSR, 80);
	/*for (j = 0; j < 80; j++)nfsr1[j] = NFSR[j];
	for (j = 0; j < 80; j++)lfsr1[j] = LFSR[j];*/
	Extractstate_BSW_arbitrary(state2, LFSR, NFSR,keystr);
	for (j = 0; j < 80; j++)nfsr1[j] = NFSR[j];
	for (j = 0; j < 80; j++)lfsr1[j] = LFSR[j];
	Diff_back(state2, LFSR, NFSR);
	Keystreamgen(2, z, LFSR, NFSR);

	state2[i1].state ^= 0x1;
	state2[i2].state ^= 0x1;
	state2[i3].state ^= 0x1;
	state2[i4].state ^= 0x1;
	state2[i5].state ^= 0x1;

	Diff_back(state2, lfsr1, nfsr1);
	Keystreamgen(2, z1, lfsr1, nfsr1);

	for (pattern = 0.0, j = 0; j < 2; j++)pattern += (z[j] ^ z1[j]) * pow(2, j);
	return((u32)pattern);

}

void Computing_difnca_arbitray(u32 ksd[4], u8 ns[8], u8 i1, u8 i2, u8 keystr[2])
{
	u8 z[2];
	u32 k,j;
	double pattern;

	memset(z, 0, 2);
	pattern = 0;

	
	for (k = 0; k < 64; k++){
		for (j = 2; j < 8; j++)ns[j] = (k >> (j - 2)) & 0x1;
		ns[0] = keystr[0] ^ h_function(ns[2], ns[4], ns[6]);
		ns[1] = keystr[1] ^ h_function(ns[3], ns[5], ns[7]);

		if( (i1 != 0) && (i1 != 1) )ns[i1] ^= 0x1;
		else ns[i1] = ns[i1];
		if ((i2 != 0) && (i2 != 1))ns[i2] ^= 0x1;
		else ns[i2] = ns[i2];
		

		Keystreamgen_nca2(z, ns);

		for (j = 0; j < 2; j++)z[j] ^= keystr[j];
		for (pattern = 0.0, j = 0; j < 2; j++)pattern += z[j] * pow(2, j);
		ksd[(u32)pattern] += 1;
	}
}

void Computing_innerdifnca_arbitray(u32 ksd[4], u8 ns[12], u8 i1, u8 keystr[2],u8 pattern[2])
{
	u8 z[2];
	u32 k, j;
	double sum;

	memset(z, 0, 2);
//	pattern = 0;

//	printf("to 2a\n");
	for (k = 0; k < 4; k++){
		ns[10] = (k >> 0) & 0x1;
		ns[11] = (k >> 1) & 0x1;

//		printf("to 2b\n");
		ns[2] = (pattern[0] ) & 0x1;
		ns[3] = (pattern[1] ) & 0x1;
		ns[4] = (pattern[0] >> 1) & 0x1;
		ns[5] = (pattern[1] >> 1) & 0x1;
		ns[6] = (pattern[0] >> 2) & 0x1;
		ns[7] = (pattern[1] >> 2) & 0x1;
		ns[8] = (pattern[0] >> 3) & 0x1;
		ns[9] = (pattern[1] >> 3) & 0x1;
		
//		printf("to 2c\n");
		ns[0] = (keystr[0] ^ h_function(ns[2], ns[4], ns[6], ns[8], ns[10]));
		ns[1] = (keystr[1] ^ h_function(ns[3], ns[5], ns[7], ns[9], ns[11]));

		ns[i1] ^= 0x1;
				
		Keystreamgen_nca2(z, ns);

		for (j = 0; j < 2; j++)z[j] ^= keystr[j];
		for (sum = 0.0, j = 0; j < 2; j++)sum += z[j] * pow(2, j);
		ksd[(u32)sum] += 1;
	}
}

void Computing_innerdifnca_arbitray6(u32 ksd[16], u8 ns[24], u8 i1, u8 keystr[4], u8 pattern[4])
{
	u8 z[4];
	u32 k, j;
	double sum;

	memset(z, 0, 4);
//	pattern = 0;

//	printf("to 2a\n");
	for (k = 0; k < 16; k++){
		ns[20] = (k >> 0) & 0x1;
		ns[21] = (k >> 1) & 0x1;
		ns[22] = (k >> 2) & 0x1;
		ns[23] = (k >> 3) & 0x1;
		
//		printf("to 2b\n");
		ns[4] = (pattern[0] ) & 0x1;
		ns[5] = (pattern[1] ) & 0x1;
		ns[6] = (pattern[2] ) & 0x1;
		ns[7] = (pattern[3] ) & 0x1;
		
		ns[8] = (pattern[0] >> 1) & 0x1;
		ns[9] = (pattern[1] >> 1) & 0x1;
		ns[10] = (pattern[2] >> 1) & 0x1;
		ns[11] = (pattern[3] >> 1) & 0x1;
		
		ns[12] = (pattern[0] >> 2) & 0x1;
		ns[13] = (pattern[1] >> 2) & 0x1;
		ns[14] = (pattern[2] >> 2) & 0x1;
		ns[15] = (pattern[3] >> 2) & 0x1;
		
		ns[16] = (pattern[0] >> 3) & 0x1;
		ns[17] = (pattern[1] >> 3) & 0x1;
		ns[18] = (pattern[2] >> 3) & 0x1;
		ns[19] = (pattern[3] >> 3) & 0x1;
		
//		printf("to 2c\n");
		ns[0] = (keystr[0] ^ h_function(ns[4], ns[8], ns[12], ns[16], ns[20]));
		ns[1] = (keystr[1] ^ h_function(ns[5], ns[9], ns[13], ns[17], ns[21]));
		ns[2] = (keystr[2] ^ h_function(ns[6], ns[10], ns[14], ns[18], ns[22]));
		ns[3] = (keystr[3] ^ h_function(ns[7], ns[11], ns[15], ns[19], ns[23]));
				
		if(i1 != 0)ns[i1] ^= 0x1;
		else ns[i1] = ns[i1];

		Keystreamgen_nca6(z, ns);

		for (j = 0; j < 4; j++)z[j] ^= keystr[j];
		for (sum = 0.0, j = 0; j < 4; j++)sum += z[j] * pow(2, j);
		ksd[(u32)sum] += 1;
	}
}

void Computing_innerdifnca_arbitray6a(u32 ksd[16], u8 ns[24], u8 i1, u8 i2, u8 keystr[4], u8 pattern[4])
{
	u8 z[4];
	u32 k, j;
	double sum;

	memset(z, 0, 4);
	//	pattern = 0;

	//	printf("to 2a\n");
	for (k = 0; k < 16; k++){
		ns[20] = (k >> 0) & 0x1;
		ns[21] = (k >> 1) & 0x1;
		ns[22] = (k >> 2) & 0x1;
		ns[23] = (k >> 3) & 0x1;

		//		printf("to 2b\n");
		ns[4] = (pattern[0]) & 0x1;
		ns[5] = (pattern[1]) & 0x1;
		ns[6] = (pattern[2]) & 0x1;
		ns[7] = (pattern[3]) & 0x1;

		ns[8] = (pattern[0] >> 1) & 0x1;
		ns[9] = (pattern[1] >> 1) & 0x1;
		ns[10] = (pattern[2] >> 1) & 0x1;
		ns[11] = (pattern[3] >> 1) & 0x1;

		ns[12] = (pattern[0] >> 2) & 0x1;
		ns[13] = (pattern[1] >> 2) & 0x1;
		ns[14] = (pattern[2] >> 2) & 0x1;
		ns[15] = (pattern[3] >> 2) & 0x1;

		ns[16] = (pattern[0] >> 3) & 0x1;
		ns[17] = (pattern[1] >> 3) & 0x1;
		ns[18] = (pattern[2] >> 3) & 0x1;
		ns[19] = (pattern[3] >> 3) & 0x1;

//		printf("to 2c\n");
		ns[0] = (keystr[0] ^ h_function(ns[4], ns[8], ns[12], ns[16], ns[20]));
		ns[1] = (keystr[1] ^ h_function(ns[5], ns[9], ns[13], ns[17], ns[21]));
		ns[2] = (keystr[2] ^ h_function(ns[6], ns[10], ns[14], ns[18], ns[22]));
		ns[3] = (keystr[3] ^ h_function(ns[7], ns[11], ns[15], ns[19], ns[23]));

		if (i1 != 0)ns[i1] ^= 0x1;
		else ns[i1] = ns[i1];
		if (i2 != 0)ns[i2] ^= 0x1;
		else ns[i2] = ns[i2];

		Keystreamgen_nca6(z, ns);

		for (j = 0; j < 4; j++)z[j] ^= keystr[j];
		for (sum = 0.0, j = 0; j < 4; j++)sum += z[j] * pow(2, j);
		ksd[(u32)sum] += 1;
	}
}

double Computing_dfrate(struct ncastate *state2, u8 *LFSR, u8 *NFSR, u8 i1, u8 i2, u8 i3, u8 i4,u8 i5, u32 index)
{
	u8 z[2], z1[2], lfsr1[40], nfsr1[40];
	u8 j;
	double pattern, counter, fc;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(lfsr1, 0, 40);
	memset(nfsr1, 0, 40);
	counter = 0.0;
	fc = 0.0;

	while (counter < SamSIZE)
	{
		randomIV(NFSR, 40);
		randomIV(LFSR, 40);
		for (j = 0; j < 40; j++)nfsr1[j] = NFSR[j];
		for (j = 0; j < 40; j++)lfsr1[j] = LFSR[j];
		Extractstate(state2, LFSR, NFSR);
		Diff_back(state2, LFSR, NFSR);
		Keystreamgen(2, z, LFSR, NFSR);

		state2[i1].state ^= 0x1;
		state2[i2].state ^= 0x1;
		state2[i3].state ^= 0x1;
		state2[i4].state ^= 0x1;
		state2[i5].state ^= 0x1;

		Diff_back(state2, lfsr1, nfsr1);
		Keystreamgen(2, z1, lfsr1, nfsr1);

		for (pattern = 0.0, j = 0; j < 2; j++)pattern += (z[j] ^ z1[j]) * pow(2, j);
		if ((u32)pattern == index)fc += 1.0;

		counter += 1;
	}
	return((fc) / (double)SamSIZE);
}

double Computing_dfrateBSW(struct ncastate *state2, u8 *LFSR, u8 *NFSR, u8 i1, u8 i2, u8 i3, u8 i4, u8 i5, u32 index)  // This routine is used to compute the rate of the input difference ISD with respect to the KSD=0x1f   
{
	u8 z[2], z1[2], lfsr1[40], nfsr1[40],zero[2];
	u8 j;
	double pattern, counter, fc;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(zero, 0, 2);
	memset(lfsr1, 0, 40);
	memset(nfsr1, 0, 40);
	counter = 0.0;
	fc = 0.0;
	
	while (counter < SamSIZE)
	{
		randomIV(NFSR, 40);
		randomIV(LFSR, 40);
		for (j = 0; j < 40; j++)nfsr1[j] = NFSR[j];
		for (j = 0; j < 40; j++)lfsr1[j] = LFSR[j];
		Extractstate_BSW2(state2, LFSR, NFSR);
		Diff_back(state2, LFSR, NFSR);
		Keystreamgen(2, z, LFSR, NFSR);

		state2[i1].state ^= 0x1;
		state2[i2].state ^= 0x1;
		state2[i3].state ^= 0x1;
		state2[i4].state ^= 0x1;
		state2[i5].state ^= 0x1;

		Diff_back(state2, lfsr1, nfsr1);
		Keystreamgen(2, z1, lfsr1, nfsr1);

		if (memcmp(z1, zero, 2) == 0){
//				for (pattern = 0.0, j = 0; j < 5; j++)pattern += (z[j] ^ z1[j]) * pow(2, j);
//				if ((u32)pattern == index)fc += 1.0;
			fc += 1.0;
		}

//		counter += 1;
		counter += 1;
	}
	return((fc) / (double)SamSIZE);
}

double Computing_dfrateBSW2(struct ncastate *state2, u8 *LFSR, u8 *NFSR, u8 i1, u8 i2, u8 i3, u8 i4, u8 i5, u32 index)  // This routine is used to compute the rate of the input difference ISD with respect to the KSD=0x1f   
{
	u8 z[2], z1[2], lfsr1[40], nfsr1[40], zero[2];
	u8 j;
	double pattern, counter, fc;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(zero, 0, 2);
	memset(lfsr1, 0, 40);
	memset(nfsr1, 0, 40);
	counter = 0.0;
	fc = 0.0;
	for (j = 0; j < 2; j++)zero[j] = 0x1;

	while (counter < SamSIZE)
	{
		randomIV(NFSR, 40);
		randomIV(LFSR, 40);
		for (j = 0; j < 40; j++)nfsr1[j] = NFSR[j];
		for (j = 0; j < 40; j++)lfsr1[j] = LFSR[j];
		Extractstate_BSW(state2, LFSR, NFSR);
		Diff_back(state2, LFSR, NFSR);
		Keystreamgen(2, z, LFSR, NFSR);

		state2[i1].state ^= 0x1;
		state2[i2].state ^= 0x1;
		state2[i3].state ^= 0x1;
		state2[i4].state ^= 0x1;
		state2[i5].state ^= 0x1;

		Diff_back(state2, lfsr1, nfsr1);
		Keystreamgen(2, z1, lfsr1, nfsr1);

		if (memcmp(z1, zero, 2) == 0){
			fc += 1.0;
		}

		counter += 1;
	}
	return((fc) / (double)SamSIZE);
}

double Computing_dfrateBSW_arbitrary(struct ncastate *state2, u8 *LFSR, u8 *NFSR, u8 i1, u8 i2, u8 i3, u8 i4, u8 i5,u8 index,u8 keystr[2])  // This routine is used to compute the rate of the input difference ISD with respect to the KSD=0x1f   
{
	u8 z[2], z1[2], lfsr1[80], nfsr1[80], zero[2];
	u8 j;
	double pattern, counter, fc;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(zero, 0, 2);
	memset(lfsr1, 0, 80);
	memset(nfsr1, 0, 80);
	counter = 0.0;
	fc = 0.0;
	for (j = 0; j < 2; j++)zero[j] = ( (index >> j) & 0x1 ) ^ keystr[j];

	while (counter < SamSIZE)
	{
		randomIV(NFSR, 80);
		randomIV(LFSR, 80);
		for (j = 0; j < 80; j++)nfsr1[j] = NFSR[j];
		for (j = 0; j < 80; j++)lfsr1[j] = LFSR[j];
		Extractstate_BSW_arbitrary(state2, LFSR, NFSR,keystr);
		Diff_back(state2, LFSR, NFSR);
		Keystreamgen(2, z, LFSR, NFSR);

		state2[i1].state ^= 0x1;
		state2[i2].state ^= 0x1;
		state2[i3].state ^= 0x1;
		state2[i4].state ^= 0x1;
		state2[i5].state ^= 0x1;

		Diff_back(state2, lfsr1, nfsr1);
		Keystreamgen(2, z1, lfsr1, nfsr1);

		if (memcmp(z1, zero, 2) == 0){
			fc += 1.0;
		}

		counter += 1;
	}
	return((fc) / (double)SamSIZE);
}

double Computing_dfratenca_arbitrary(u8 ns[12], u8 i1, u8 i2, u8 i3, u8 index, u8 keystr[2])  // This routine is used to compute the rate of the input difference ISD with respect to the KSD=0x03   
{
	u8 z[2];
	u32 k, j, ksd[4];
	double pattern;

	memset(z, 0, 2);
	memset(ksd, 0, 4);
	pattern = 0;


	for (k = 0; k < 1024; k++){
		for (j = 2; j < 12; j++)ns[j] = (k >> (j - 2)) & 0x1;
		ns[0] = keystr[0] ^ h_function(ns[2], ns[4], ns[6], ns[8], ns[10]);
		ns[1] = keystr[1] ^ h_function(ns[3], ns[5], ns[7], ns[9], ns[11]);

		if ((i1 != 0) && (i1 != 1))ns[i1] ^= 0x1;
		else ns[i1] = ns[i1];
		if ((i2 != 0) && (i2 != 1))ns[i2] ^= 0x1;
		else ns[i2] = ns[i2];
		if ((i3 != 0) && (i3 != 1))ns[i3] ^= 0x1;
		else ns[i3] = ns[i3];

		Keystreamgen_nca2(z, ns);

		for (j = 0; j < 2; j++)z[j] ^= keystr[j];
		for (pattern = 0.0, j = 0; j < 2; j++)pattern += z[j] * pow(2, j);
		ksd[(u32)pattern] += 1;
	}
	return((double)(ksd[index])/1024);
}


struct isd2 *Ksd2isd2_nca_arbitrary(u8 ns[8], struct isd2 *ph, u8 index, u8 z[2])       //index is the KSD and z[2] is one of the 2-bit keystream prefix
{
	u8 i1, i2, i3, i, flag=0;
	u32 ksd[4], r;
	struct isd2 *p1, *p2, *head, *p3;
	
	for (i = 0; i < 4; i++)ksd[i] = 0;
	
	p1 = (struct isd2 *)calloc(1, sizeof(struct isd2));
	if (NULL == p1)
	{
		printf("Error in calloc in Ksd2isd2_nca_arbitrary.\n");
		return 0;
	}
	p2 = NULL;
	p3 = NULL;
	
	r = 0;
	head = NULL;
	
	for (i1 = 2; i1 < 7; i1++)                                     // index starts from 0 to include the less weight isd
	for (i2 = i1 + 1; i2 < 8; i2++)
	{
		flag = 0;
		for (i = 0; i < 4; i++)ksd[i] = 0;
		Computing_difnca_arbitray(ksd, ns, i1, i2, z);

		if ((ksd[index] != 0) && (r == 0) && (flag == 0))
		{
			r++;

			head = p1; 
			p1->position[0] = i1;
			p1->position[1] = i2;
		
			
			p1->fc = (double)(ksd[index])/64;  
//			printf("%u,%u,%u: %f\n", p1->position[0], p1->position[1], p1->position[2],p1->fc);
			p1->next = NULL;
			p2 = p1;
			flag = 1;
			
		}
		else if ((ksd[index] != 0) && (r > 0) && (flag == 0))
		{
			p3 = (struct isd2 *)calloc(1, sizeof(struct isd2));
			if (NULL == p3)
			{
				printf("Error in calloc in Ksd2isd2_nca_arbitrary.\n");
				return 0;
			}

			p3->position[0] = i1;
			p3->position[1] = i2;
			
			p3->fc = (double)(ksd[index]) / 64;
//			printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
			p3->next = NULL;
			p2->next = p3;
			p2 = p3;
			flag = 1;
			r++;

		}
		else continue;
	}

	
	
	for (i1 = 2; i1 < 8; i1++){

		flag = 0;
		for (i = 0; i < 4; i++)ksd[i] = 0;
		Computing_difnca_arbitray(ksd, ns, i1, 0, z);

		if ((ksd[index] != 0) && (flag == 0)){
			p3 = (struct isd2 *)calloc(1, sizeof(struct isd2));
			if (NULL == p3)
			{
				printf("Error in calloc in Ksd2isd2_nca_arbitrary.\n");
				return 0;
			}

			p3->position[0] = i1;
			p3->position[1] = 0;


			p3->fc = (double)(ksd[index]) / 64;
//			printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
			p3->next = NULL;
			p2->next = p3;
			p2 = p3;
			flag = 1;
			r++;
		}
		else continue;
	}

	flag = 0;
	for (i = 0; i < 4; i++)ksd[i] = 0;
	Computing_difnca_arbitray(ksd, ns, 0, 0, z);

	if ((ksd[index] != 0) && (flag == 0)){
		p3 = (struct isd2 *)calloc(1, sizeof(struct isd2));
		if (NULL == p3)
		{
			printf("Error in calloc in Ksd2isd2_nca_arbitrary.\n");
			return 0;
		}

		p3->position[0] = 0;
		p3->position[1] = 0;
	
		p3->fc = (double)(ksd[index]) / 64;
//		printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
		p3->next = NULL;
		p2->next = p3;
		p2 = p3;
		flag = 1;
		r++;
	}
	

	ph = head;
	return(head);
}



struct isd3 *Ksd2isd2_innernca_arbitrary(u8 ns[12], struct isd3 *ph, u8 index, u8 z[2], u8 pattern[2])       //index is the KSD and z[2] is one of the 2-bit keystream prefix, pattern[2] is the 2 4 lfsr bits in h-function
{
	u8 i1, i, flag = 0;
	u32 ksd[4], r;
	struct isd3 *p1, *p2, *head, *p3;

	for (i = 0; i < 4; i++)ksd[i] = 0;

	p1 = (struct isd3 *)calloc(1, sizeof(struct isd3));
	if (NULL == p1)
	{
		printf("Error in calloc in Ksd2isd2_innernca_arbitrary.\n");
		return 0;
	}
	p2 = (struct isd3 *)calloc(1, sizeof(struct isd3));
	if (NULL == p2)
	{
		printf("Error in calloc in Ksd2isd2_innernca_arbitrary.\n");
		return 0;
	}
//	p2 = NULL;
	p3 = NULL;

	r = 0;
	head = NULL;

//	printf("to 1\n");
	for (i1 = 10; i1 < 12; i1++)                                     // index starts from 0 to include the less weight isd
	{
		flag = 0;
		for (i = 0; i < 4; i++)ksd[i] = 0;
//		printf("to 2\n");
		Computing_innerdifnca_arbitray(ksd, ns, i1, z, pattern);
//		printf("to 3\n");
		if ((ksd[index] != 0) && (r == 0) && (flag == 0))
		{
			r++;

			head = p1;
			p1->position[0] = i1;
		
			p1->fc = (double)(ksd[index]) / 4;
//			printf("%u,%u,%u: %f\n", p1->position[0], p1->position[1], p1->position[2],p1->fc);
			p1->next = NULL;
			p2 = p1;
			flag = 1;

		}
		else if ((ksd[index] != 0) && (r > 0) && (flag == 0))
		{
			p3 = (struct isd3 *)calloc(1, sizeof(struct isd3));
			if (NULL == p3)
			{
				printf("Error in calloc in Ksd2isd2_innernca_arbitrary.\n");
				return 0;
			}

			p3->position[0] = i1;
		
			p3->fc = (double)(ksd[index]) / 4;
//			printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
			p3->next = NULL;
			p2->next = p3;
			p2 = p3;
			flag = 1;
			r++;

		}
		else continue;
	}

	flag = 0;
	for (i = 0; i < 4; i++)ksd[i] = 0;
	Computing_innerdifnca_arbitray(ksd, ns, 0, z, pattern);

	if ((ksd[index] != 0) && (flag == 0)){
		p3 = (struct isd3 *)calloc(1, sizeof(struct isd3));
		if (NULL == p3)
		{
			printf("Error in calloc in Ksd2isd2_innernca_arbitrary.\n");
			return 0;
		}

		p3->position[0] = 0;
		
		p3->fc = (double)(ksd[index]) / 4;
//		printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
		p3->next = NULL;
		p2->next = p3;
		p2 = p3;
		flag = 1;
		r++;
	}


	ph = head;
	return(head);
}

struct isd4 *Ksd2isd6_innernca_arbitrary(u8 ns[24], struct isd4 *ph, u8 index, u8 z[4], u8 pattern[4])       //index is the KSD and z[4] is one of the 4-bit keystream prefix, pattern[4] is the 4 4-lfsr-bits in h-function
{
	u8 i1, i2, i, flag = 0;
	u32 ksd[16], r;
	struct isd4 *p1, *p2, *head, *p3;

	for (i = 0; i < 16; i++)ksd[i] = 0;

	p1 = (struct isd4 *)calloc(1, sizeof(struct isd4));
	if (NULL == p1)
	{
		printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
		return 0;
	}
	p2 = (struct isd4 *)calloc(1, sizeof(struct isd4));
	if (NULL == p2)
	{
		printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
		return 0;
	}
//	p2 = NULL;
	p3 = NULL;

	r = 0;
	head = NULL;

//	printf("to 1\n");
	for (i1 = 20; i1 < 24; i1++)                                     // index starts from 0 to include the less weight isd
	{
		flag = 0;
		for (i = 0; i < 16; i++)ksd[i] = 0;
//		printf("to 2\n");
		Computing_innerdifnca_arbitray6a(ksd, ns, i1, 0,z, pattern);
//		printf("to 3\n");
		if ((ksd[index] != 0) && (r == 0) && (flag == 0))
		{
			r++;

			head = p1;
			p1->position[0] = i1;
			p1->position[1] = 0;

			p1->fc = (double)(ksd[index]) / 16;
//			printf("%u,%u,%u: %f\n", p1->position[0], p1->position[1], p1->position[2],p1->fc);
			p1->next = NULL;
			p2 = p1;
			flag = 1;

		}
		else if ((ksd[index] != 0) && (r > 0) && (flag == 0))
		{
			p3 = (struct isd4 *)calloc(1, sizeof(struct isd4));
			if (NULL == p3)
			{
				printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
				return 0;
			}

			p3->position[0] = i1;
			p3->position[0] = 0;

			p3->fc = (double)(ksd[index]) / 16;
//			printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
			p3->next = NULL;
			p2->next = p3;
			p2 = p3;
			flag = 1;
			r++;

		}
		else continue;
	}

	flag = 0;
	for (i = 0; i < 16; i++)ksd[i] = 0;
	Computing_innerdifnca_arbitray6a(ksd, ns, 0,0, z, pattern);

	if ((ksd[index] != 0) && (flag == 0)){
		p3 = (struct isd4 *)calloc(1, sizeof(struct isd4));
		if (NULL == p3)
		{
			printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
			return 0;
		}

		p3->position[0] = 0;
		p3->position[1] = 0;

		p3->fc = (double)(ksd[index]) / 16;
//		printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
		p3->next = NULL;
		p2->next = p3;
		p2 = p3;
		flag = 1;
		r++;
	}

	flag = 0;
	for (i1 = 20; i1 < 23; i1++)
	for (i2 = i1 + 1; i2 < 24;i2++)                             // index starts from 0 to include the less weight isd
	{
		flag = 0;
		for (i = 0; i < 16; i++)ksd[i] = 0;
//		printf("to 2\n");
		Computing_innerdifnca_arbitray6a(ksd, ns, i1, i2, z, pattern);
//		printf("to 3\n");
		if ((ksd[index] != 0) && (r > 0) && (flag == 0))
		{
			r++;

			p3 = (struct isd4 *)calloc(1, sizeof(struct isd4));
			if (NULL == p3)
			{
				printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
				return 0;
			}
			p3->position[0] = i1;
			p3->position[1] = i2;

			p3->fc = (double)(ksd[index]) / 16;
//			printf("%u,%u,%u: %f\n", p1->position[0], p1->position[1], p1->position[2],p1->fc);
			p3->next = NULL;
			p2 = p3;
			flag = 1;

		}
		else if ((ksd[index] != 0) && (r > 0) && (flag == 0))
		{
			p3 = (struct isd4 *)calloc(1, sizeof(struct isd4));
			if (NULL == p3)
			{
				printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
				return 0;
			}

			p3->position[0] = i1;
			p3->position[0] = i2;

			p3->fc = (double)(ksd[index]) / 16;
//			printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
			p3->next = NULL;
			p2->next = p3;
			p2 = p3;
			flag = 1;
			r++;

		}
		else continue;
	}


	ph = head;
	return(head);
}


struct isd4 *Ksd2isd6_innernca_arbitrary4(u8 ns[24], struct isd4 *ph, u8 index, u8 z[4], u8 pattern[4])       //index is the KSD and z[4] is one of the 4-bit keystream prefix, pattern[4] is the 4 4-lfsr-bits in h-function
{
	u8 i1, i2, i, flag = 0;
	u32 ksd[16], r;
	struct isd4 *p1, *p2, *head, *p3;

	for (i = 0; i < 16; i++)ksd[i] = 0;

	p1 = (struct isd4 *)calloc(1, sizeof(struct isd4));
	if (NULL == p1)
	{
		printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
		return 0;
	}
	p2 = (struct isd4 *)calloc(1, sizeof(struct isd4));
	if (NULL == p2)
	{
		printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
		return 0;
	}
	//	p2 = NULL;
	p3 = NULL;

	r = 0;
	head = NULL;

	//	printf("to 1\n");
	for (i1 = 20; i1 < 24; i1++)                                     // index starts from 0 to include the less weight isd
	{
		flag = 0;
		for (i = 0; i < 16; i++)ksd[i] = 0;
//		printf("to 2\n");
		Computing_innerdifnca_arbitray6a(ksd, ns, i1, 0, z, pattern);
//		printf("to 3\n");
		if ((ksd[index] != 0) && (r == 0) && (flag == 0))
		{
			r++;

			head = p1;
			p1->position[0] = i1;
			p1->position[1] = 0;

			p1->fc = (double)(ksd[index]) / 16;
//			printf("%u,%u,%u: %f\n", p1->position[0], p1->position[1], p1->position[2],p1->fc);
			p1->next = NULL;
			p2 = p1;
			flag = 1;

		}
		else if ((ksd[index] != 0) && (r > 0) && (flag == 0))
		{
			p3 = (struct isd4 *)calloc(1, sizeof(struct isd4));
			if (NULL == p3)
			{
				printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
				return 0;
			}

			p3->position[0] = i1;
			p3->position[0] = 0;

			p3->fc = (double)(ksd[index]) / 16;
//			printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
			p3->next = NULL;
			p2->next = p3;
			p2 = p3;
			flag = 1;
			r++;

		}
		else continue;
	}

	flag = 0;
	for (i = 0; i < 16; i++)ksd[i] = 0;
	Computing_innerdifnca_arbitray6a(ksd, ns, 0, 0, z, pattern);

	if ((ksd[index] != 0) && (flag == 0)){
		p3 = (struct isd4 *)calloc(1, sizeof(struct isd4));
		if (NULL == p3)
		{
			printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
			return 0;
		}

		p3->position[0] = 0;
		p3->position[1] = 0;

		p3->fc = (double)(ksd[index]) / 16;
//		printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
		p3->next = NULL;
		p2->next = p3;
		p2 = p3;
		flag = 1;
		r++;
	}

	flag = 0;
	for (i1 = 20; i1 < 23; i1++)
	for (i2 = i1 + 1; i2 < 24; i2++)                             // index starts from 0 to include the less weight isd
	{
		flag = 0;
		for (i = 0; i < 16; i++)ksd[i] = 0;
//		printf("to 2\n");
		Computing_innerdifnca_arbitray6a(ksd, ns, i1, i2, z, pattern);
//		printf("to 3\n");
		if ((ksd[index] != 0) && (r > 0) && (flag == 0))
		{
			r++;

			p3 = (struct isd4 *)calloc(1, sizeof(struct isd4));
			if (NULL == p3)
			{
				printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
				return 0;
			}
			p3->position[0] = i1;
			p3->position[1] = i2;

			p3->fc = (double)(ksd[index]) / 16;
//			printf("%u,%u,%u: %f\n", p1->position[0], p1->position[1], p1->position[2],p1->fc);
			p3->next = NULL;
			p2 = p3;
			flag = 1;

		}
		else if ((ksd[index] != 0) && (r > 0) && (flag == 0))
		{
			p3 = (struct isd4 *)calloc(1, sizeof(struct isd4));
			if (NULL == p3)
			{
				printf("Error in calloc in Ksd2isd6_innernca_arbitrary.\n");
				return 0;
			}

			p3->position[0] = i1;
			p3->position[0] = i2;

			p3->fc = (double)(ksd[index]) / 16;
			//			printf("%u,%u,%u: %f\n", p3->position[0], p3->position[1], p3->position[2], p3->fc);
			p3->next = NULL;
			p2->next = p3;
			p2 = p3;
			flag = 1;
			r++;

		}
		else continue;
	}


	ph = head;
	return(head);
}


struct isd2 *Ksd2isd2_BSW_arbitrary(struct ncastate *state2, struct isd2 *ph, u8 index,u8 z[2])       //index is the KSD and z[2] is one of the 2-bit keystream prefix
{
	u8 i1, i2, i3, i4,i5,flag;
	u8 LFSR[80], NFSR[80];
	u32 r, pattern, counter;
	struct isd2 *p1, *p2, *head, *p3;
	double fc;

	memset(LFSR, 0, 80);
	memset(NFSR, 0, 80);
	counter = 0;
	fc = 0.0;

	p1 = (struct isd2 *)calloc(1, sizeof(struct isd2));
	if (NULL == p1)
	{
		printf("Error in calloc in Ksd2isd2.\n");
		return 0;
	}
	p2 = NULL;
	p3 = NULL;
//	printf("Index=%u\n", index);

	// the following embedded loop is repeated 7140 times, which is the exact value of Binomial(41-5=36,3)       
	r = 0;
	head = NULL;
	for (i1 = 1; i1 < 142; i1++)
	for (i2 = i1 + 1; i2 < 143; i2++)
	for (i3 = i2 + 1; i3 < 144; i3++)
	for (i4 = i3 + 1; i4 < 145; i4++)
	for (i5 = i4 + 1; i5 < 146; i5++)
	{
//		printf("%u,%u,%u,%u,%u\n",i1,i2,i3,i4,i5);
		if ((((i1 >= 1) && (i1 <= 5)) || ((i1 >= 10) && (i1 <= 11)) || ((i1 >= 43) && (i1 <= 44)) || ((i1 >= 56) && (i1 <= 57)) || ((i1 >= 83) && (i1 <= 84)) || ((i1 >= 105) && (i1 <= 106)) || ((i1 >= 126) && (i1 <= 127)) || ((i1 >= 144) && (i1 <= 145)) || ((i1 >= 63) && (i1 <= 64))) &&
			(((i2 >= 1) && (i2 <= 5)) || ((i2 >= 10) && (i2 <= 11)) || ((i2 >= 43) && (i2 <= 44)) || ((i2 >= 56) && (i2 <= 57)) || ((i2 >= 83) && (i2 <= 84)) || ((i2 >= 105) && (i2 <= 106)) || ((i2 >= 126) && (i2 <= 127)) || ((i2 >= 144) && (i2 <= 145)) || ((i2 >= 63) && (i2 <= 64))) &&
			(((i3 >= 1) && (i3 <= 5)) || ((i3 >= 10) && (i3 <= 11)) || ((i3 >= 43) && (i3 <= 44)) || ((i3 >= 56) && (i3 <= 57)) || ((i3 >= 83) && (i3 <= 84)) || ((i3 >= 105) && (i3 <= 106)) || ((i3 >= 126) && (i3 <= 127)) || ((i3 >= 144) && (i3 <= 145)) || ((i3 >= 63) && (i3 <= 64))) &&
			(((i4 >= 1) && (i4 <= 5)) || ((i4 >= 10) && (i4 <= 11)) || ((i4 >= 43) && (i4 <= 44)) || ((i4 >= 56) && (i4 <= 57)) || ((i4 >= 83) && (i4 <= 84)) || ((i4 >= 105) && (i4 <= 106)) || ((i4 >= 126) && (i4 <= 127)) || ((i4 >= 144) && (i4 <= 145)) || ((i4 >= 63) && (i4 <= 64))) &&
			(((i5 >= 1) && (i5 <= 5)) || ((i5 >= 10) && (i5 <= 11)) || ((i5 >= 43) && (i5 <= 44)) || ((i5 >= 56) && (i5 <= 57)) || ((i5 >= 83) && (i5 <= 84)) || ((i5 >= 105) && (i5 <= 106)) || ((i5 >= 126) && (i5 <= 127)) || ((i5 >= 144) && (i5 <= 145)) || ((i5 >= 63) && (i5 <= 64))))
		{
			counter = 0;
			flag = 0;
			fc = 0.0;

			while (counter < SamplingSIZE)
			{
				pattern = Computing_difBSW_arbitray(state2, LFSR, NFSR, i1, i2, i3, i4, i5, z);         // to compute to see if the current ISD can generate desired KSD pattern

				if ((index == pattern) && (r == 0) && (flag == 0))                         // for BSW sampling, the index is 0
				{
					r++;
					
					fc = Computing_dfrateBSW_arbitrary(state2, LFSR, NFSR, i1, i2, i3, i4, i5, index, z);
//					printf("%u,%u,%u,%u,%u\n", i1, i2, i3, i4, i5);
					head = p1; p1->position[0] = i1; p1->position[1] = i2; p1->position[2] = i3; p1->position[3] = i4; p1->position[4] = i5; p1->fc = fc;  p1->next = NULL;
					p2 = p1;
//					printf("%u,%u,%u;\n", i1, i2, i3);
					flag = 1;

//					printf("%u,%u,%u: %.7f;\n", i1, i2, i3, fc);
					break;
					/*fc += 1.0;
					continue;*/
				}
				else if ((index == pattern) && (r > 0) && (flag == 0))
				{
					fc = Computing_dfrateBSW_arbitrary(state2, LFSR, NFSR, i1, i2, i3, i4,i5,index,z);
					p3 = (struct isd2 *)calloc(1, sizeof(struct isd2));
					if (NULL == p3)
					{
						printf("Error in calloc in Ksd2isd.\n");
						return 0;
					}
//					printf("%u,%u,%u,%u,%u\n", i1, i2, i3, i4, i5);
					p3->position[0] = i1;
					p3->position[1] = i2;
					p3->position[2] = i3;
					p3->position[3] = i4;
					p3->position[4] = i5;
//					printf("%u,%u,%u,%u,%u\n", p3->position[0], p3->position[1], p3->position[2], p3->position[3], p3->position[4]);
					p3->fc = fc;


					p3->next = NULL;
					p2->next = p3;
					p2 = p3;
					flag = 1;
					r++;

//					printf("%u,%u,%u: %.7f;\n", i1, i2, i3,fc);
					break;
					/*fc += 1.0;
					continue;*/
				}
//				else continue;
				counter = counter + 1;
			}

		}
		else continue;
	}

	ph = head;
	return(head);
}


void Isd2Pr_diversity(struct isd2 *isd2_table[4][4], double Pr[4][4])
{
	u8 i, j;
	u32 r;
	double sum;
	struct isd2 *pp1, *pp2;

	pp1 = NULL;
	pp2 = NULL;
	for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	{
		sum = 0.0;
		pp1 = isd2_table[i][j];
		pp2 = pp1;
		r = 0;

		do{
			sum = sum + pp2->fc;
			r++;
			pp2 = pp2->next;
		} while (pp2 != NULL);

		Pr[i][j] = (double)(sum) / (r);
	}
}

void output_Pr(double Pr[4][4])
{
	u8 i, j, r;

	r = 0;
	for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++){
		printf("i=%u,j=%u:", i, j);
		printf("%f ", Pr[i][j]);
		r++;
		if ((r % 4) == 0)printf("\n");
	}

}



void Output_isd2(struct isd2 *ph)                                          // output the link of isd2_table[i]
{
	struct isd2 *p;

	p = ph;

	if (ph != NULL)
	{
		do
		{
			printf("%u,%u: %.7f;\n", p->position[0], p->position[1],p->fc);
			p = p->next;
		} while (p != NULL);
	}
}

void Output_isd3(struct isd3 *ph)                                          // output the link of isd2_table2[i]
{
	struct isd3 *p;

	p = ph;

	if (ph != NULL)
	{
		do
		{
			printf("%u, %.7f;\n", p->position[0], p->fc);
			p = p->next;
		} while (p != NULL);
	}
}

void Output_isd4(struct isd4 *ph)                                          // output the link of isd2_table2[i]
{
	struct isd4 *p;

	p = ph;

	if (ph != NULL)
	{
		do
		{
			printf("%u,%u:%.7f;\n", p->position[0], p->position[1], p->fc);
			p = p->next;
		} while (p != NULL);
	}
}

void output_datalist(struct datalist *list)
{
	u32 i;
	u8 j;

	for (i = 0; i < ListSIZE; i++){
		for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
		printf("\n");
		for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
		printf("\n");
		printf("The ISnca_sum is %f\n", (list + i)->ISnca_sum);
		printf("The position is %u\n", (list + i)->position);
	}
}

void output_datalist2(struct datalist2 *list)
{
	struct datalist2 *p;
	u8 j;

	p = list;

	if (list != NULL)
	{
		do
		{
			for (j = 0; j < 2; j++)printf("%x", p->prefix[j]);
			printf("\n");
			for (j = 0; j < 23; j++)printf("%x", p->ISnca[j]);
			printf("\n");
			printf("The ISnca_sum is %f\n", p->ISnca_sum);
			printf("The position is %u\n", p->position);
			p = p->next;
		} while (p != NULL);
	}

}


struct datalist2 *SetA_generation_Self_contained2(u8 *LFSR, u8 *NFSR, struct datalist2 *list)        // generate the first set in the self-contained method with BSW sampling
{
	u8 z[2], j, k,zero[2],lfsr[40],nfsr[40];
	u32 i;
	struct ncastate state2[80];
	double sum;

	memset(z, 0, 2);
	memset(lfsr, 0, 40);
	memset(nfsr, 0, 40);
	memset(zero, 0, 2);
	memset(state2, 0, sizeof(struct ncastate) * 80);
	

//	printf("****************************************\n");
	// to ensure the continuous generation of 2-bit keystream prefix
	randomIV(LFSR, 40);
	randomIV(NFSR, 40);

	i = 0;
//	while (i < ListSIZE){
	while (i < 1){
		Extractstate_BSW(state2, LFSR, NFSR);
		Diff_back(state2, LFSR, NFSR);
	/*	if (i == ListSIZE - 190){
			memcpy(lfsr, LFSR, 40);
			memcpy(nfsr, NFSR, 40);
		}*/
		Keystreamgen(2, z, LFSR, NFSR);
		
		k = 0;
		j = 0;
		while (j < 23){
			if ((state2 + k)->flag == 1){
				(list + i)->ISnca[j] = (state2 + k)->state;
			    k++;
			    j++;
			}
			else k++;
		}

		for (j = 0; j < 2; j++)(list + i)->prefix[j] = z[j];
		for (sum = 0.0, j = 0; j < 23; j++)sum += ((list + i)->ISnca[j]) * pow(2, j);
		(list + i)->ISnca_sum = sum;
		(list + i)->position = i;
		(list + i)->next = list + i + 1;
		i++;
	}

	/*memcpy(LFSR, lfsr, 40);
	memcpy(NFSR, nfsr, 40);*/
	return(list);
}

u8 h(u8 x0, u8 x1, u8 x2, u8 x3, u8 x4)
{
	
	u8 t;

	t = x1 ^ x4 ^ (x0 & x3) ^ (x2 & x3) ^ (x3 & x4) ^ (x0 & x1 & x2) ^ (x0 & x2 & x3) ^ (x0 & x2 & x4) ^ (x1 & x2 & x4) ^ (x2 & x3 & x4);
	return(t);
}

u8 verify_stateBSW(struct datalist2 *list)
{
	u8 i;
	u8 num;

	for (num = 0, i = 6; i <= 7; i++){
		if (list->ISnca[i] == (list->ISnca[i - 6] ^ list->ISnca[i - 6 + 2] ^ list->ISnca[i - 6 + 4] ^ list->ISnca[i - 6 + 9] ^ h(list->ISnca[i - 6 + 11], list->ISnca[i - 6 + 13], list->ISnca[i - 6 + 15],
			list->ISnca[i - 6 + 17], list->ISnca[i - 6 + 8])))num++;
	}

	//	printf("The number is %u\n", num);

	return(num);
}

u8 verify_keystr(struct datalist2 *list, u8 inter[23])
{
	u8 i, LFSR[80], NFSR[80], z[2];

	memset(z, 0, 2);
	memset(LFSR, 0, 80);
	memset(NFSR, 0, 80);
	for (i = 1; i <= 5; i++)NFSR[i] = inter[i - 1];
	for (i = 10; i <= 11; i++)NFSR[i] = inter[i - 5];
	for (i = 31; i <= 32; i++)NFSR[i] = inter[i - 24];
	for (i = 43; i <= 44; i++)NFSR[i] = inter[i - 34];
	for (i = 56; i <= 57; i++)NFSR[i] = inter[i - 45];
	for (i = 63; i <= 64; i++)NFSR[i] = inter[i - 50];
	for (i = 3; i <= 4; i++)LFSR[i] = inter[i + 12];
	for (i = 25; i <= 26; i++)LFSR[i] = inter[i - 8];
	for (i = 46; i <= 47; i++)LFSR[i] = inter[i - 27];
	for (i = 64; i <= 65; i++)LFSR[i] = inter[i - 43];
	Keystreamgen(2, z, LFSR, NFSR);

	return(memcmp(list->prefix, z, 2));
}

u8 verify_ncakeystr4(struct datalist2_inner *list, u8 inter[24])
{
	u8 z[4];

	memset(z,0,4);
	Keystreamgen_nca6(z, inter);

	return(memcmp(list->prefix, z, 4));
}

u8 verify_ncakeystr(struct datalist2 *list, u8 inter[8])
{
	u8 z[2];

	memset(z, 0, 2);
	Keystreamgen_nca2(z, inter);

	return(memcmp(list->prefix, z, 2));
}

struct datalist2 *SetB_generation_Self_contained(u8 *LFSR, u8 *NFSR, struct datalist *list, struct isd2 *isd2_table[4])     // generate the second set in the self-contained method and test it
{
	u8 z[2], j, k, z1[2], index = 0, inter[19], flag = 0;
	u32 i, i1, r = 0, success = 0;
	struct ncastate state2[80];
	struct datalist *p1;
	struct datalist2 *p2, *head, *p3, *p4, *p5;         // p2,p3,p4,p5 and head are used to store the second set in the self-contained method
	struct isd2 *pp1, *pp2;
	double sum;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(state2, 0, sizeof(struct ncastate) * 80);
	memset(inter, 0, 19);
	pp1 = NULL;
	pp2 = NULL;

	randomIV(NFSR, 40);
	randomIV(LFSR, 40);
	printf("The randomly generated LFSR and NFSR are:");
	for (j = 0; j < 40; j++)printf("%u", LFSR[j]);
	printf("\n");
	for (j = 0; j < 40; j++)printf("%u", NFSR[j]);
	printf("\n");

	p1 = NULL;
	p1 = (struct datalist *)calloc(1, sizeof(struct datalist));
	if (p1 == NULL){
		printf("error in calloc in SetB_generation_Self_contained\n");
	}
	p2 = NULL;
	p2 = (struct datalist2 *)calloc(1, sizeof(struct datalist2));
	if (p2 == NULL){
		printf("error in calloc in SetB_generation_Self_contained*\n");
	}
	p3 = NULL;
	head = NULL;
	p4 = NULL;
	p5 = NULL;

	i1 = 0;
	do
	{
		Extractstate(state2, LFSR, NFSR);
		Diff_back(state2, LFSR, NFSR);

		k = 0;
		j = 0;
		while (j < 19){
			if ((state2 + k)->flag == 1){
				p1->ISnca[j] = (state2 + k)->state;
				k++;
				j++;
			}
			else k++;
		}

		Keystreamgen(2, z, LFSR, NFSR);
		printf("The generated keystream is:");
		for (j = 0; j < 2; j++)printf("%x", z[j]);
		printf("\n");
		for (j = 0; j < 2; j++)p1->prefix[j] = z[j];
		for (sum = 0.0, j = 0; j < 19; j++)sum += (p1->ISnca[j]) * pow(2, j);
		p1->ISnca_sum = sum;
		p1->position = i1;
		i1++;
		for (j = 0; j < 19; j++)inter[j] = p1->ISnca[j];

		for (i = 0; i < ListSIZE; i++){
			//		for (i = 0; i < 1; i++){
			printf("The %u th element is:", i);
			for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
			printf("\n");
			printf("The outside element is:");
			for (j = 0; j < 2; j++)printf("%x", p1->prefix[j]);
			printf("\n");
			for (j = 0; j < 2; j++)z1[j] = p1->prefix[j] ^ (list + i)->prefix[j];
			for (index = 0, j = 0; j < 2; j++)index += (u8)((z1[j]) * pow(2, j));
			printf("The KSD is %u \n", index);

			pp1 = isd2_table[index];                             // isd2_table is used to store the precomputation information
			pp2 = pp1;
			do
			{
				inter[pp2->position[0]] ^= 0x1;
				inter[pp2->position[1]] ^= 0x1;
				inter[pp2->position[2]] ^= 0x1;
				inter[pp2->position[3]] ^= 0x1;
				inter[pp2->position[4]] ^= 0x1;

				flag = verify_keystr(list + i, inter);
				if (flag == 0){
					r++;
					if (r == 1){
						head = p2;
						for (j = 0; j < 2; j++)p2->prefix[j] = (list + i)->prefix[j];
						for (j = 0; j < 19; j++)p2->ISnca[j] = inter[j];
						for (sum = 0.0, j = 0; j < 19; j++)sum += (p2->ISnca[j]) * pow(2, j);
						p2->ISnca_sum = sum;
						p2->position = r - 1;
						p2->next = NULL;
						p4 = p2;

						for (j = 0; j < 2; j++)printf("%x", p2->prefix[j]);
						printf("\n");
						for (j = 0; j < 19; j++)printf("%x", p2->ISnca[j]);
						printf("\n");
						printf("The checkup sum is %f\n", p2->ISnca_sum);
						printf("The new position is %u and the original position is %u\n", p2->position, i);
						printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
					}
					else
					{
						p5 = (struct datalist2 *)calloc(1, sizeof(struct datalist2));
						if (p5 == NULL){
							printf("error in calloc in SetB_generation_Self_contained i1=%u\n", i1);
						}
						for (j = 0; j < 2; j++)p5->prefix[j] = (list + i)->prefix[j];
						for (j = 0; j < 19; j++)p5->ISnca[j] = inter[j];
						for (sum = 0.0, j = 0; j < 19; j++)sum += (p5->ISnca[j]) * pow(2, j);
						p5->ISnca_sum = sum;
						p5->position = r - 1;

						p5->next = NULL;
						p4->next = p5;
						p4 = p5;

						for (j = 0; j < 2; j++)printf("%x", p5->prefix[j]);
						printf("\n");
						for (j = 0; j < 19; j++)printf("%x", p5->ISnca[j]);
						printf("\n");
						printf("The checkup sum is %f\n", p5->ISnca_sum);
						printf("The new position is %u and the original position is %u\n", p5->position, i);
						printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
					}
				}
				pp2 = pp2->next;

			} while (pp2 != NULL);
			printf("**********************************************\n");
		}
		success++;
	} while (success < 1);

	return(head);
}


void flip_inter(u8 inter[23], struct isd2 *pp2)
{
	u8 i1, i2, i3, i4, i5;

	i1 = pp2->position[0];
	i2 = pp2->position[1];
	i3 = pp2->position[2];
	i4 = pp2->position[3];
	i5 = pp2->position[4];

	if (((i1 >= 1) && (i1 <= 5)))i1 = i1 - 1;
	else if (((i1 >= 10) && (i1 <= 11)))i1 = i1 - 5;
	else if (((i1 >= 31) && (i1 <= 32)))i1 = i1 - 24;
	else if (((i1 >= 43) && (i1 <= 44)))i1 = i1 - 34;
	else if (((i1 >= 56) && (i1 <= 57)))i1 = i1 - 45;
	else if (((i1 >= 63) && (i1 <= 64)))i1 = i1 - 50;
	else if (((i1 >= 83) && (i1 <= 84)))i1 = i1 - 68;
	else if (((i1 >= 105) && (i1 <= 106)))i1 = i1 - 88;
	else if (((i1 >= 126) && (i1 <= 127)))i1 = i1 - 107;
	else if (((i1 >= 144) && (i1 <= 145)))i1 = i1 - 123;
	else i1 = i1;

	if (((i2 >= 1) && (i2 <= 5)))i2 = i2 - 1;
	else if (((i2 >= 10) && (i2 <= 11)))i2 = i2 - 5;
	else if (((i2 >= 31) && (i2 <= 32)))i2 = i2 - 24;
	else if (((i2 >= 43) && (i2 <= 44)))i2 = i2 - 34;
	else if (((i2 >= 56) && (i2 <= 57)))i2 = i2 - 45;
	else if (((i2 >= 63) && (i2 <= 64)))i2 = i2 - 50;
	else if (((i2 >= 83) && (i2 <= 84)))i2 = i2 - 68;
	else if (((i2 >= 105) && (i2 <= 106)))i2 = i2 - 88;
	else if (((i2 >= 126) && (i2 <= 127)))i2 = i2 - 107;
	else if (((i2 >= 144) && (i2 <= 145)))i2 = i2 - 123;
	else i2 = i2;

	if (((i3 >= 1) && (i3 <= 5)))i3 = i3 - 1;
	else if (((i3 >= 10) && (i3 <= 11)))i3 = i3 - 5;
	else if (((i3 >= 31) && (i3 <= 32)))i3 = i3 - 24;
	else if (((i3 >= 43) && (i3 <= 44)))i3 = i3 - 34;
	else if (((i3 >= 56) && (i3 <= 57)))i3 = i3 - 45;
	else if (((i3 >= 63) && (i3 <= 64)))i3 = i3 - 50;
	else if (((i3 >= 83) && (i3 <= 84)))i3 = i3 - 68;
	else if (((i3 >= 105) && (i3 <= 106)))i3 = i3 - 88;
	else if (((i3 >= 126) && (i3 <= 127)))i3 = i3 - 107;
	else if (((i3 >= 144) && (i3 <= 145)))i3 = i3 - 123;
	else i3 = i3;

	if (((i4 >= 1) && (i4 <= 5)))i4 = i4 - 1;
	else if (((i4 >= 10) && (i4 <= 11)))i4 = i4 - 5;
	else if (((i4 >= 31) && (i4 <= 32)))i4 = i4 - 24;
	else if (((i4 >= 43) && (i4 <= 44)))i4 = i4 - 34;
	else if (((i4 >= 56) && (i4 <= 57)))i4 = i4 - 45;
	else if (((i4 >= 63) && (i4 <= 64)))i4 = i4 - 50;
	else if (((i4 >= 83) && (i4 <= 84)))i4 = i4 - 68;
	else if (((i4 >= 105) && (i4 <= 106)))i4 = i4 - 88;
	else if (((i4 >= 126) && (i4 <= 127)))i4 = i4 - 107;
	else if (((i4 >= 144) && (i4 <= 145)))i4 = i4 - 123;
	else i4 = i4;

	if (((i5 >= 1) && (i5 <= 5)))i5 = i5 - 1;
	else if (((i5 >= 10) && (i5 <= 11)))i5 = i5 - 5;
	else if (((i5 >= 31) && (i5 <= 32)))i5 = i5 - 24;
	else if (((i5 >= 43) && (i5 <= 44)))i5 = i5 - 34;
	else if (((i5 >= 56) && (i5 <= 57)))i5 = i5 - 45;
	else if (((i5 >= 63) && (i5 <= 64)))i5 = i5 - 50;
	else if (((i5 >= 83) && (i5 <= 84)))i5 = i5 - 68;
	else if (((i5 >= 105) && (i5 <= 106)))i5 = i5 - 88;
	else if (((i5 >= 126) && (i5 <= 127)))i5 = i5 - 107;
	else if (((i5 >= 144) && (i5 <= 145)))i5 = i5 - 123;
	else i5 = i5;

	inter[i1] ^= 0x1;
	inter[i2] ^= 0x1;
	inter[i3] ^= 0x1;
	inter[i4] ^= 0x1;
	inter[i5] ^= 0x1;

}

void flip_ncainter(u8 inter[8], struct isd2 *pp2)
{
	u8 i1=0, i2=0;

	i1 = pp2->position[0];
	i2 = pp2->position[1];
		
	if ( (i1 !=0 ) && (i1 != 1) )inter[i1] ^= 0x1;
	else inter[i1] = inter[i1];
	if ((i2 != 0) && (i2 != 1))inter[i2] ^= 0x1;
	else inter[i2] = inter[i2];
		
}

void flip_ncainter_inner(u8 inter[12], struct isd3 *pp2)
{
	u8 i1 = 0;

	i1 = pp2->position[0];
//	printf("The i1=%u\n",i1);
	
	if (i1 != 0)inter[i1] ^= 0x1;
	else inter[i1] = inter[i1];
}

void flip_ncainter_inner4(u8 inter[24], struct isd4 *pp2)
{
	u8 i1 = 0, i2 = 0;

	i1 = pp2->position[0];
	i2 = pp2->position[1];
//	printf("The i1=%u\n",i1);

	if (i1 != 0)inter[i1] ^= 0x1;
	else inter[i1] = inter[i1];
	if (i2 != 0)inter[i2] ^= 0x1;
	else inter[i2] = inter[i2];
}

struct datalist2 *SetB_generation_Self_contained2(u8 *LFSR, u8 *NFSR, struct datalist2 *list, struct isd2 *isd2_table[4][4], double *counter)     // generate the second set in the self-contained method and test it
{
	u8 z[2], j, k, z1[2], index = 0, inter[19], flag = 0, zero[2];
	u32 i, i1, r = 0, success = 0;
	struct ncastate state2[80];
	struct datalist *p1;
	struct datalist2 *p2, *head, *p3, *p4, *p5;        // p2,p3,p4,p5 and head are used to store the second set in the self-contained method
	struct isd2 *pp1, *pp2;
	double sum, counter1;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(zero, 0, 2);
	memset(state2, 0, sizeof(struct ncastate) * 80);
	memset(inter, 0, 19);
	pp1 = NULL;
	pp2 = NULL;

	p1 = NULL;
	p1 = (struct datalist *)calloc(1, sizeof(struct datalist));
	if (p1 == NULL){
		printf("error in calloc in SetB_generation_Self_contained2\n");
	}
	p2 = NULL;
	p2 = (struct datalist2 *)calloc(1, sizeof(struct datalist2));
	if (p2 == NULL){
		printf("error in calloc in SetB_generation_Self_contained2*\n");
	}
	p3 = NULL;
	head = NULL;
	p4 = NULL;
	p4 = p2;
	p5 = NULL;
	for (j = 0; j < 2; j++)zero[j] = 0x1;

	i1 = 0;
	counter1 = 0;
	r = 0;

	//	while (success < 512){
	while (success < 524288){                            //524288 = 2^19
		randomIV(NFSR, 40);
		randomIV(LFSR, 40);
		Extractstate_BSW2(state2, LFSR, NFSR);           // generate the 0x1f keystream prefix
		/*for (j = 0; j < 40; j++)printf("%x", NFSR[j]);
		printf("\n");
		for (j = 0; j < 40; j++)printf("%x", LFSR[j]);
		printf("\n");*/
		Diff_back(state2, LFSR, NFSR);
		Keystreamgen(2, z, LFSR, NFSR);
		/*for (j = 0; j < 5; j++)printf("%x",z[j]);
		printf("\n");*/

		/*printf("The generated keystream is (should be 0x1f):");
		for (j = 0; j < 5; j++)printf("%x", z[j]);
		printf("\n");*/

		k = 0;
		j = 0;
		while (j < 19){
			if ((state2 + k)->flag == 1){
				p1->ISnca[j] = (state2 + k)->state;
				k++;
				j++;
			}
			else k++;
		}
		/*for (j = 0; j < 41; j++)printf("%x", p1->ISnca[j]);
		printf("\n");*/

		for (j = 0; j < 2; j++)p1->prefix[j] = z[j];
		for (sum = 0.0, j = 0; j < 19; j++)sum += (p1->ISnca[j]) * pow(2, j);
		p1->ISnca_sum = sum;
		p1->position = i1;
		i1++;
		for (j = 0; j < 19; j++)inter[j] = p1->ISnca[j];

		//		for (i = 0; i < ListSIZE; i++){
		for (i = 0; i < 1; i++){
			/*printf("The %u th element is:", i);
			for (j = 0; j < 5; j++)printf("%x", (list + i)->prefix[j]);
			printf("\n");
			printf("The outside element is (should be 0x1f):");
			for (j = 0; j < 5; j++)printf("%x", p1->prefix[j]);
			printf("\n");*/
			for (j = 0; j < 2; j++)z1[j] = p1->prefix[j] ^ (list + i)->prefix[j];
			for (index = 0, j = 0; j < 2; j++)index += (u8)((z1[j]) * pow(2, j));
			/*printf("The KSD is %u \n", index);*/

			/*printf("The set A is:\n");
			for (j = 0; j < 35; j++)printf("%x", (list + i)->ISnca[j]);
			printf("\n");*/
			pp1 = isd2_table[index][3];                             // isd2_table is used to store the precomputation information
			pp2 = pp1;
			do
			{
				for (j = 0; j < 19; j++)inter[j] = p1->ISnca[j];
				/*printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
				printf("The set B before fliping is:\n");
				for (j = 0; j < 35; j++)printf("%x", inter[j]);
				printf("\n");
				printf("The filped positions are: ");
				printf("%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2]);*/
				flip_inter(inter, pp2);
				/*printf("The set B after fliping is:\n");
				for (j = 0; j < 35; j++)printf("%x", inter[j]);
				printf("\n");
				printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");*/
				/*inter[pp2->position[0]] ^= 0x1;
				inter[pp2->position[1]] ^= 0x1;
				inter[pp2->position[2]] ^= 0x1;*/

				flag = verify_keystr(list + i, inter);
				/*printf("The flag is %u\n", flag);*/
				if (flag == 0){
					r++;
					/*printf("%d\n", memcmp((list + i)->ISnca, inter, 10));
					printf("%d\n", memcmp(((list + i)->ISnca) + 15, inter+15, 21) );*/
					//					printf("-----------------------------------------------\n");
					if (r == 1){
						//						if (memcmp((list + i)->ISnca, inter, 35) == 0){
						if ((memcmp((list + i)->ISnca, inter, 6) == 0) && (memcmp(((list + i)->ISnca) + 8, inter + 8, 11) == 0)) {
							(*counter)++;
							printf("Found!!!\n");
							printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 19; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							head = p2;
							for (j = 0; j < 2; j++)p2->prefix[j] = (list + i)->prefix[j];
							for (j = 0; j < 19; j++)p2->ISnca[j] = inter[j];
							for (sum = 0.0, j = 0; j < 19; j++)sum += (p2->ISnca[j]) * pow(2, j);
							p2->ISnca_sum = sum;
							p2->position = r - 1;
							p2->next = NULL;
							p4 = p2;

							for (j = 0; j < 2; j++)printf("%x", p2->prefix[j]);
							printf("\n");
							for (j = 0; j < 19; j++)printf("%x", p2->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", p2->ISnca_sum);
							printf("The new position is %u and the original position is %u\n", p2->position, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
						}
						//						else break;						

					}
					else
					{
						//						if (memcmp((list + i)->ISnca, inter, 41) == 0){
						if ((memcmp((list + i)->ISnca, inter, 6) == 0) && (memcmp(((list + i)->ISnca) + 8, inter + 8, 11) == 0)){
							(*counter)++;
							printf("Found!!!\n");
							printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 19; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							p5 = (struct datalist2 *)calloc(1, sizeof(struct datalist2));
							if (p5 == NULL){
								printf("error in calloc in SetB_generation_Self_contained2 i1=%u\n", i1);
							}
							for (j = 0; j < 2; j++)p5->prefix[j] = (list + i)->prefix[j];
							for (j = 0; j < 19; j++)p5->ISnca[j] = inter[j];
							for (sum = 0.0, j = 0; j < 19; j++)sum += (p5->ISnca[j]) * pow(2, j);
							p5->ISnca_sum = sum;
							p5->position = r - 1;

							p5->next = NULL;
							p4->next = p5;
							p4 = p5;

							for (j = 0; j < 2; j++)printf("%x", p5->prefix[j]);
							printf("\n");
							for (j = 0; j < 19; j++)printf("%x", p5->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", p5->ISnca_sum);
							printf("The new position is %u and the original position is %u\n", p5->position, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
						}
						//						else break;
					}
				}
				pp2 = pp2->next;
			} while (pp2 != NULL);
			//			printf("**********************************************\n");
		}
		success++;
	}   // corresponding to while(successs<)

	printf("The number of candidates is %u\n", r);
	printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");

	return(head);
}

u8 inter_keystr(u8 inter[19])
{
	u8 i, LFSR[40], NFSR[40], z[2], zero[2];

	memset(z, 0, 2);
	memset(zero, 0, 2);
	memset(LFSR, 0, 40);
	memset(NFSR, 0, 40);
	for (i = 1; i <= 2; i++)NFSR[i] = inter[i - 1];
	for (i = 4; i <= 5; i++)NFSR[i] = inter[i - 2];
	for (i = 7; i <= 8; i++)NFSR[i] = inter[i - 3];
	for (i = 21; i <= 22; i++)NFSR[i] = inter[i - 15];
	for (i = 31; i <= 33; i++)NFSR[i] = inter[i - 23];
	for (i = 3; i <= 4; i++)LFSR[i] = inter[11 + i - 3];
	for (i = 13; i <= 14; i++)LFSR[i] = inter[13 + i - 13];
	for (i = 17; i <= 18; i++)LFSR[i] = inter[15 + i - 17];
	for (i = 23; i <= 24; i++)LFSR[i] = inter[17 + i - 23];
	Keystreamgen(2, z, LFSR, NFSR);

	return(memcmp(zero, z, 2));
}

u8 inter_keystr_arbitrary(u8 inter[23], u8 keystr[2])
{
	u8 i, LFSR[80], NFSR[80], z[2];

	memset(z, 0, 2);
	memset(LFSR, 0, 80);
	memset(NFSR, 0, 80);

	for (i = 1; i <= 5; i++)NFSR[i] = inter[i - 1];
	for (i = 10; i <= 11; i++)NFSR[i] = inter[i - 5];
	for (i = 31; i <= 32; i++)NFSR[i] = inter[i - 24];
	for (i = 43; i <= 44; i++)NFSR[i] = inter[i - 34];
	for (i = 56; i <= 57; i++)NFSR[i] = inter[i - 45];
	for (i = 63; i <= 64; i++)NFSR[i] = inter[i - 50];
	for (i = 3; i <= 4; i++)LFSR[i] = inter[i + 12];
	for (i = 25; i <= 26; i++)LFSR[i] = inter[i - 8];
	for (i = 46; i <= 47; i++)LFSR[i] = inter[i - 27];
	for (i = 64; i <= 65; i++)LFSR[i] = inter[i - 43];
	Keystreamgen(2, z, LFSR, NFSR);

	return(memcmp(keystr, z, 2));
}

double ISD_BSW(u8 *LFSR, u8 *NFSR, struct isd2 *isd2_table)
{
	struct ncastate state2[80];
	u8 z[2], inter[19], k, j, temp[19], flag = 0;
	u32 i, r;
	double counter, fc;
	struct isd2 *pp1, *pp2;

	counter = 0.0;
	fc = 0;
	memset(state2, 0, sizeof(struct ncastate) * 80);
	memset(z, 0, 2);
	memset(inter, 0, 19);
	memset(temp, 0, 19);
	pp1 = NULL;
	pp2 = NULL;

	while (counter < SamplingSIZE)
	{
		randomIV(LFSR, 40);
		randomIV(NFSR, 40);
		Extractstate_BSW2(state2, LFSR, NFSR);
		Diff_back(state2, LFSR, NFSR);

		k = 0;
		j = 0;
		while (j < 19){
			if ((state2 + k)->flag == 1){
				inter[j] = (state2 + k)->state;
				k++;
				j++;
			}
			else k++;
		}
		//		memcpy(temp, inter, 41);
		Keystreamgen(2, z, LFSR, NFSR);

		pp1 = isd2_table;
		pp2 = pp1;
		r = 0;

		//		memcpy(inter, temp, 41);
		inter[pp2->position[0]] ^= 0x1;
		inter[pp2->position[1]] ^= 0x1;
		inter[pp2->position[2]] ^= 0x1;

		flag = inter_keystr(inter);
		if (flag == 0){
			fc += 1.0;
		}

		counter += 1.0;
	}

	printf("The ratio is %f\n", fc / SamplingSIZE);
	return((double)fc / SamplingSIZE);
}


double ISD_BSW_arbitrary(u8 *LFSR, u8 *NFSR, struct isd2 *isd2_table, u8 index, u8 keystr[2])
{
	struct ncastate state2[160];
	u8 z[2], z1[2], inter[23], k, j, temp[23], flag = 0, indexz[2], nfsr1[80], lfsr1[80];
	u32 i, r;
	double counter, fc;
	struct isd2 *pp1, *pp2;

	counter = 0.0;
	fc = 0;
	memset(state2, 0, sizeof(struct ncastate) * 160);
	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(indexz, 0, 2);
	memset(inter, 0, 23);
	memset(nfsr1, 0, 80);
	memset(lfsr1, 0, 80);
	memset(temp, 0, 23);
	pp1 = NULL;
	pp2 = NULL;
	for (j = 0; j < 2; j++)indexz[j] = ((index >> j) & 0x1) ^ keystr[j];

	while (counter < SamplingSIZE)
	{
		randomIV(LFSR, 80);
		randomIV(NFSR, 80);
		for (j = 0; j < 80; j++)nfsr1[j] = NFSR[j];
		for (j = 0; j < 80; j++)lfsr1[j] = LFSR[j];
		Extractstate_BSW_arbitrary(state2, LFSR, NFSR, keystr);
		Diff_back(state2, LFSR, NFSR);

		/*k = 0;
		j = 0;
		while (j < 35){
		if ((state2 + k)->flag == 1){
		inter[j] = (state2 + k)->state;
		k++;
		j++;
		}
		else k++;
		}*/
		//		memcpy(temp, inter, 35);
		Keystreamgen(2, z, LFSR, NFSR);

		pp1 = isd2_table;
		pp2 = pp1;
		r = 0;

		state2[pp2->position[0]].state ^= 0x1;
		state2[pp2->position[1]].state ^= 0x1;
		state2[pp2->position[2]].state ^= 0x1;
		state2[pp2->position[3]].state ^= 0x1;
		state2[pp2->position[4]].state ^= 0x1;

		Diff_back(state2, lfsr1, nfsr1);
		Keystreamgen(2, z1, lfsr1, nfsr1);

		if (memcmp(z1, indexz, 2) == 0){
			//				for (pattern = 0.0, j = 0; j < 5; j++)pattern += (z[j] ^ z1[j]) * pow(2, j);
			//				if ((u32)pattern == index)fc += 1.0;
			fc += 1.0;
		}

		//		memcpy(inter, temp, 41);
		/*inter[pp2->position[0]] ^= 0x1;
		inter[pp2->position[1]] ^= 0x1;
		inter[pp2->position[2]] ^= 0x1;*/

		/*flag = inter_keystr_arbitrary(inter,indexz);
		if (flag == 0){
		fc += 1.0;
		}*/

		counter += 1.0;
	}

	printf("The ratio is %f\n", fc / SamplingSIZE);
	return((double)fc / SamplingSIZE);
}


double ISD_BSW_arbitrary2(u8 *LFSR, u8 *NFSR, struct isd2 *isd2_table, u8 index, u8 keystr[2])
{
	struct ncastate state2[160];
	u8 z[2], inter[23], k, j, temp[23], flag = 0, indexz[2], i1, i2, i3, i4, i5;
	u32 i, r;
	double counter, fc;
	struct isd2 *pp1, *pp2;

	counter = 0.0;
	fc = 0;
	memset(state2, 0, sizeof(struct ncastate) * 160);
	memset(z, 0, 2);
	memset(indexz, 0, 2);
	memset(inter, 0, 23);
	memset(temp, 0, 23);
	pp1 = NULL;
	pp2 = NULL;
	for (j = 0; j < 2; j++)indexz[j] = ((index >> j) & 0x1) ^ keystr[j];

	while (counter < SamplingSIZE)
	{
		randomIV(LFSR, 80);
		randomIV(NFSR, 80);
		Extractstate_BSW_arbitrary(state2, LFSR, NFSR, keystr);
		Diff_back(state2, LFSR, NFSR);

		k = 0;
		j = 0;
		while (j < 23){
			if ((state2 + k)->flag == 1){
				inter[j] = (state2 + k)->state;
				k++;
				j++;
			}
			else k++;
		}
		Keystreamgen(2, z, LFSR, NFSR);

		pp1 = isd2_table;
		pp2 = pp1;
		r = 0;

		i1 = pp2->position[0];
		i2 = pp2->position[1];
		i3 = pp2->position[2];
		i4 = pp2->position[3];
		i5 = pp2->position[4];

		if (((i1 >= 1) && (i1 <= 5)))i1 = i1 - 1;
		else if (((i1 >= 10) && (i1 <= 11)))i1 = i1 - 5;
		else if (((i1 >= 31) && (i1 <= 32)))i1 = i1 - 24;
		else if (((i1 >= 43) && (i1 <= 44)))i1 = i1 - 34;
		else if (((i1 >= 56) && (i1 <= 57)))i1 = i1 - 45;
		else if (((i1 >= 63) && (i1 <= 64)))i1 = i1 - 50;
		else if (((i1 >= 83) && (i1 <= 84)))i1 = i1 - 68;
		else if (((i1 >= 105) && (i1 <= 106)))i1 = i1 - 88;
		else if (((i1 >= 126) && (i1 <= 127)))i1 = i1 - 107;
		else if (((i1 >= 144) && (i1 <= 145)))i1 = i1 - 123;
		else i1 = i1;

		if (((i2 >= 1) && (i2 <= 5)))i2 = i2 - 1;
		else if (((i2 >= 10) && (i2 <= 11)))i2 = i2 - 5;
		else if (((i2 >= 31) && (i2 <= 32)))i2 = i2 - 24;
		else if (((i2 >= 43) && (i2 <= 44)))i2 = i2 - 34;
		else if (((i2 >= 56) && (i2 <= 57)))i2 = i2 - 45;
		else if (((i2 >= 63) && (i2 <= 64)))i2 = i2 - 50;
		else if (((i2 >= 83) && (i2 <= 84)))i2 = i2 - 68;
		else if (((i2 >= 105) && (i2 <= 106)))i2 = i2 - 88;
		else if (((i2 >= 126) && (i2 <= 127)))i2 = i2 - 107;
		else if (((i2 >= 144) && (i2 <= 145)))i2 = i2 - 123;
		else i2 = i2;

		if (((i3 >= 1) && (i3 <= 5)))i3 = i3 - 1;
		else if (((i3 >= 10) && (i3 <= 11)))i3 = i3 - 5;
		else if (((i3 >= 31) && (i3 <= 32)))i3 = i3 - 24;
		else if (((i3 >= 43) && (i3 <= 44)))i3 = i3 - 34;
		else if (((i3 >= 56) && (i3 <= 57)))i3 = i3 - 45;
		else if (((i3 >= 63) && (i3 <= 64)))i3 = i3 - 50;
		else if (((i3 >= 83) && (i3 <= 84)))i3 = i3 - 68;
		else if (((i3 >= 105) && (i3 <= 106)))i3 = i3 - 88;
		else if (((i3 >= 126) && (i3 <= 127)))i3 = i3 - 107;
		else if (((i3 >= 144) && (i3 <= 145)))i3 = i3 - 123;
		else i3 = i3;

		if (((i4 >= 1) && (i4 <= 5)))i4 = i4 - 1;
		else if (((i4 >= 10) && (i4 <= 11)))i4 = i4 - 5;
		else if (((i4 >= 31) && (i4 <= 32)))i4 = i4 - 24;
		else if (((i4 >= 43) && (i4 <= 44)))i4 = i4 - 34;
		else if (((i4 >= 56) && (i4 <= 57)))i4 = i4 - 45;
		else if (((i4 >= 63) && (i4 <= 64)))i4 = i4 - 50;
		else if (((i4 >= 83) && (i4 <= 84)))i4 = i4 - 68;
		else if (((i4 >= 105) && (i4 <= 106)))i4 = i4 - 88;
		else if (((i4 >= 126) && (i4 <= 127)))i4 = i4 - 107;
		else if (((i4 >= 144) && (i4 <= 145)))i4 = i4 - 123;
		else i4 = i4;

		if (((i5 >= 1) && (i5 <= 5)))i5 = i5 - 1;
		else if (((i5 >= 10) && (i5 <= 11)))i5 = i5 - 5;
		else if (((i5 >= 31) && (i5 <= 32)))i5 = i5 - 24;
		else if (((i5 >= 43) && (i5 <= 44)))i5 = i5 - 34;
		else if (((i5 >= 56) && (i5 <= 57)))i5 = i5 - 45;
		else if (((i5 >= 63) && (i5 <= 64)))i5 = i5 - 50;
		else if (((i5 >= 83) && (i5 <= 84)))i5 = i5 - 68;
		else if (((i5 >= 105) && (i5 <= 106)))i5 = i5 - 88;
		else if (((i5 >= 126) && (i5 <= 127)))i5 = i5 - 107;
		else if (((i5 >= 144) && (i5 <= 145)))i5 = i5 - 123;
		else i5 = i5;

		inter[i1] ^= 0x1;
		inter[i2] ^= 0x1;
		inter[i3] ^= 0x1;
		inter[i4] ^= 0x1;
		inter[i5] ^= 0x1;

		flag = inter_keystr_arbitrary(inter, indexz);
		if (flag == 0){
			fc += 1.0;
		}

		counter += 1.0;
	}

	printf("The ratio is %f\n", fc / SamplingSIZE);
	return((double)fc / SamplingSIZE);
}


void TotalISD_BSW(u8 *LFSR, u8 *NFSR, struct isd2 *isd2_table[4])
{
	u32 i, r;
	double *fc;
	struct isd2 *pp1, *pp2;

	r = 0;
	fc = NULL;
	pp1 = NULL;
	pp2 = NULL;

	r = length_isd2(isd2_table[3]);
	fc = (double *)calloc(r + 1, sizeof(double));
	if (fc == NULL){
		printf("error in calloc in TotalISD_BSW\n");
	}

	pp1 = isd2_table[3];
	pp2 = pp1;
	i = 0;
	do{
		fc[i] = ISD_BSW(LFSR, NFSR, pp2);
		printf("The ISD is %u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
		printf("The %uth ISD ratio is %f\n", i, fc[i]);
		printf("----------------------------------------------\n");
		i++;
		pp2 = pp2->next;
	} while (pp2 != NULL);

}

void TotalISD_BSW_arbitrary(u8 *LFSR, u8 *NFSR, struct isd2 *isd2_table[4][4], u8 index, u8 keystr[2])
{
	u32 i, r;
	double *fc, sum;
	struct isd2 *pp1, *pp2;

	r = 0;
	fc = NULL;
	pp1 = NULL;
	pp2 = NULL;

	for (sum = 0.0, i = 0; i < 2; i++)sum += keystr[i] * pow(2, i);
	r = length_isd2(isd2_table[index][(u8)sum]);
	fc = (double *)calloc(r + 1, sizeof(double));
	if (fc == NULL){
		printf("error in calloc in TotalISD_BSW_arbitrary\n");
	}

	pp1 = isd2_table[index][(u8)sum];
	pp2 = pp1;
	i = 0;
	do{
//		fc[i] = ISD_BSW_arbitrary(LFSR, NFSR, pp2,index,keystr);
		fc[i] = ISD_BSW_arbitrary2(LFSR, NFSR, pp2, index, keystr);
		printf("The ISD is %u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
		printf("The %uth ISD ratio is %f\n", i, fc[i]);
		printf("----------------------------------------------\n");
		i++;
		pp2 = pp2->next;
	} while (pp2 != NULL);

	free(fc);
}

void verify_ISD_BSW(u8 *LFSR, u8 *NFSR, struct isd2 *isd2_table)
{
	struct ncastate state2[80];
	u8 z[2], inter[19], k, j, temp[19], flag = 0;
	u32 i, r;
	double counter, fc;
	struct isd2 *pp1, *pp2;

	counter = 0.0;
	fc = 0;
	memset(state2, 0, sizeof(struct ncastate) * 80);
	memset(z, 0, 2);
	memset(inter, 0, 19);
	memset(temp, 0, 19);
	pp1 = NULL;
	pp2 = NULL;

	while (counter < SamplingSIZE)
	{
		randomIV(LFSR, 40);
		randomIV(NFSR, 40);
		Extractstate_BSW2(state2, LFSR, NFSR);
		Diff_back(state2, LFSR, NFSR);

		k = 0;
		j = 0;
		while (j < 19){
			if ((state2 + k)->flag == 1){
				inter[j] = (state2 + k)->state;
				k++;
				j++;
			}
			else k++;
		}
		//		memcpy(temp, inter, 41);
		Keystreamgen(2, z, LFSR, NFSR);

		pp1 = isd2_table;
		pp2 = pp1;
		r = 0;

		//		memcpy(inter, temp, 41);
		inter[8] ^= 0x1;
		inter[3] ^= 0x1;
		inter[7] ^= 0x1;

		flag = inter_keystr(inter);
		if (flag == 0){
			fc += 1.0;
		}

		counter += 1.0;
	}

	printf("The ratio is %f\n", fc / SamplingSIZE);
	return((double)fc / SamplingSIZE);
}

u8 verify_inter_keystr(u8 inter[19], u8 keystr[2])
{
	u8 i, LFSR[40], NFSR[40], z[2];

	memset(z, 0, 2);
	memset(LFSR, 0, 40);
	memset(NFSR, 0, 40);
	for (i = 1; i <= 2; i++)NFSR[i] = inter[i - 1];
	for (i = 4; i <= 5; i++)NFSR[i] = inter[i - 2];
	for (i = 7; i <= 8; i++)NFSR[i] = inter[i - 3];
	for (i = 21; i <= 22; i++)NFSR[i] = inter[i - 15];
	for (i = 31; i <= 33; i++)NFSR[i] = inter[i - 23];
	for (i = 3; i <= 4; i++)LFSR[i] = inter[11 + i - 3];
	for (i = 13; i <= 14; i++)LFSR[i] = inter[13 + i - 13];
	for (i = 17; i <= 18; i++)LFSR[i] = inter[15 + i - 17];
	for (i = 23; i <= 24; i++)LFSR[i] = inter[17 + i - 23];
	Keystreamgen(2, z, LFSR, NFSR);

	if (memcmp(keystr, z, 2) == 0)return(1);
	else return(0);
}

u8 verify_inter(u8 *LFSR, u8 *NFSR)
{
	struct ncastate state2[80];
	u8 k, j, inter[19];
	u8 z[2], t;

	memset(z, 0, 2);
	memset(inter, 0, 19);

	randomIV(LFSR, 40);
	randomIV(NFSR, 40);
	Extractstate_BSW2(state2, LFSR, NFSR);
	//	Extractstate_BSW(state2, LFSR, NFSR);
	//	Extractstate(state2, LFSR, NFSR);
	Diff_back(state2, LFSR, NFSR);

	k = 0;
	j = 0;
	while (j < 19){
		if ((state2 + k)->flag == 1){
			inter[j] = (state2 + k)->state;
			k++;
			j++;
		}
		else k++;
	}
	Keystreamgen(2, z, LFSR, NFSR);
	for (j = 0; j < 2; j++)printf("%x", z[j]);
	printf("\n");

	t = verify_inter_keystr(inter, z);
	return(t);
}

u32 length_datalist2(struct datalist2 *ph)                   // find the length of the link struct datalist2
{
	struct datalist2 *p;
	int n;

	n = 0;
	p = ph;

	if (ph != NULL)
	{
		do
		{
			n++;
			p = p->next;
		} while (p != NULL);
	}

	//	return(n - 1);                                   // return the index of the last element  
	return(n);
}


u8 *SetB_generation_Self_contained3(u8 *LFSR, u8 *NFSR, struct datalist2 *list, struct isd2 *isd2_table[4][4], u8 *try)     // generate the second set in the self-contained method and test it
{
	u8 z[2], j, k, z1[2], index = 0, inter[23], flag = 0, zero[2];
	u32 i, i1, r = 0, success = 0, counter2=0;
	struct ncastate state2[160];
	struct datalist *p1;
//	struct datalist2 *p2, *head, *p3, *p4, *p5;             // p2,p3,p4,p5 and head are used to store the second set in the self-contained method
	struct isd2 *pp1, *pp2;
	double sum, counter1;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(zero, 0, 2);
	memset(state2, 0, sizeof(struct ncastate) * 160);
	memset(inter, 0, 23);
	pp1 = NULL;
	pp2 = NULL;

	p1 = NULL;
	p1 = (struct datalist *)calloc(1, sizeof(struct datalist));
	if (p1 == NULL){
		printf("error in calloc in SetB_generation_Self_contained3\n");
	}

	zero[0] = 0x0;
	zero[1] = 0x1;


	i1 = 0;
	counter1 = 0;
	r = 0;

//	while (success < 512){
//	while (success < 128){                                        //32768 = 2^17
	while (success < 1024){
//	while (success < 2097152){                                        //1048576 = 2^20
		randomIV(NFSR, 80);
		randomIV(LFSR, 80);
		Extractstate_BSW2(state2, LFSR, NFSR);           // generate the 0x1f keystream prefix
//		Extractstate_BSW(state2, LFSR, NFSR);           // generate the 0x00 keystream prefix
//		Extractstate_BSW_arbitrary(state2, LFSR, NFSR, zero);
		/*for (j = 0; j < 40; j++)printf("%x", NFSR[j]);
		printf("\n");
		for (j = 0; j < 40; j++)printf("%x", LFSR[j]);
		printf("\n");*/
		Diff_back(state2, LFSR, NFSR);
		Keystreamgen(2, z, LFSR, NFSR);
		/*for (j = 0; j < 3; j++)printf("%x",z[j]);
		printf("\n");*/

		/*printf("The generated keystream is (should be 0x1f):");
		for (j = 0; j < 5; j++)printf("%x", z[j]);
		printf("\n");*/

		k = 0;
		j = 0;
		while (j < 23){
			if ((state2 + k)->flag == 1){
				p1->ISnca[j] = (state2 + k)->state;
				k++;
				j++;
			}
			else k++;
		}
		/*for (j = 0; j < 41; j++)printf("%x", p1->ISnca[j]);
		printf("\n");*/

		for (j = 0; j < 2; j++)p1->prefix[j] = z[j];
		for (sum = 0.0, j = 0; j < 23; j++)sum += (p1->ISnca[j]) * pow(2, j);
		p1->ISnca_sum = sum;
		p1->position = i1;
		i1++;
		for (j = 0; j < 23; j++)inter[j] = p1->ISnca[j];

//		for (i = 0; i < ListSIZE; i++){
		for (i = 0; i < 1; i++){
			/*printf("The %u th element is:", i);
			for (j = 0; j < 5; j++)printf("%x", (list + i)->prefix[j]);
			printf("\n");
			printf("The outside element is (should be 0x1f):");
			for (j = 0; j < 5; j++)printf("%x", p1->prefix[j]);
			printf("\n");*/
			for (j = 0; j < 2; j++)z1[j] = p1->prefix[j] ^ (list + i)->prefix[j];
			for (index = 0, j = 0; j < 2; j++)index += (u8)((z1[j]) * pow(2, j));
			/*printf("The KSD is %u \n", index);*/

//			pp1 = isd2_table[index][7];                             // isd2_table is used to store the precomputation information
			pp1 = isd2_table[index][3];
//			pp1 = isd2_table[index][6];                             
			pp2 = pp1;
//			while ((pp2->fc) >= 0.3)pp2 = pp2->next;            // find the correct place to begin
			do
			{
				for (j = 0; j < 23; j++)inter[j] = p1->ISnca[j];
				/*printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
				printf("The set B before fliping is:\n");
				for (j = 0; j < 35; j++)printf("%x", inter[j]);
				printf("\n");
				printf("The filped positions are: ");
				printf("%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2]);*/
				flip_inter(inter, pp2);
				/*printf("The set B after fliping is:\n");
				for (j = 0; j < 35; j++)printf("%x", inter[j]);
				printf("\n");
				printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");*/
				/*inter[pp2->position[0]] ^= 0x1;
				inter[pp2->position[1]] ^= 0x1;
				inter[pp2->position[2]] ^= 0x1;*/
				counter2++;
				flag = verify_keystr(list + i, inter);
				/*printf("The flag is %u\n", flag);*/
				if (flag == 0){
					r++;
					for (sum = 0.0, j = 0; j < 23; j++)sum += inter[j] * pow(2, j);
					try[(u32)sum] ++;
					/*printf("%d\n", memcmp((list + i)->ISnca, inter, 10));
					printf("%d\n", memcmp(((list + i)->ISnca) + 15, inter+15, 21) );*/
//					printf("-----------------------------------------------\n");
					if (r == 1){
						/*for (sum = 0.0, j = 0; j < 28; j++)sum += inter[j] * pow(2, j);
						try[r - 1] = (float)sum; */
						/*head = p2;
						for (j = 0; j < 3; j++)p2->prefix[j] = (list + i)->prefix[j];
						for (j = 0; j < 28; j++)p2->ISnca[j] = inter[j];
						for (sum = 0.0, j = 0; j < 28; j++)sum += (p2->ISnca[j]) * pow(2, j);
						p2->ISnca_sum = sum;
						p2->position = r - 1;
						p2->next = NULL;
						p4 = p2;*/

//						if (memcmp((list + i)->ISnca, inter, 35) == 0){
						if ((memcmp((list + i)->ISnca, inter, 7) == 0) && (memcmp(((list + i)->ISnca) + 9, inter + 9, 14) == 0)) {

							printf("Found!!!\n");
							/*printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
//						else break;						

					}
					else
					{
						/*for (sum = 0.0, j = 0; j < 28; j++)sum += inter[j] * pow(2, j);
						try[r - 1] = (float)sum;*/

						/*for (j = 0; j < 3; j++)p5->prefix[j] = (list + i)->prefix[j];
						for (j = 0; j < 28; j++)p5->ISnca[j] = inter[j];
						for (sum = 0.0, j = 0; j < 28; j++)sum += (p5->ISnca[j]) * pow(2, j);
						p5->ISnca_sum = sum;
						p5->position = r - 1;

						p5->next = NULL;
						p4->next = p5;
						p4 = p5;*/

//						if (memcmp((list + i)->ISnca, inter, 41) == 0){
						if ((memcmp((list + i)->ISnca, inter, 7) == 0) && (memcmp(((list + i)->ISnca) + 9, inter + 9, 14) == 0)){

							printf("Found!!!\n");
							/*printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
//						else break;
					}
				}
				pp2 = pp2->next;
			} while (pp2 != NULL);


//			printf("**********************************************\n");
		}
		success++;
	}   // corresponding to while(successs<)

	/*printf("The counter1 is %u\n", counter2);
	printf("The counter2 is %u\n", r);*/
//	printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");

	return(try);
}


u8 *SetB_generation_Self_contained4(u8 *LFSR, u8 *NFSR, struct datalist2 *list, struct isd2 *isd2_table[4][4], u8 *try)     // generate the second set in the self-contained method and test it
{
	u8 z[2], j, k, z1[2], index = 0, inter[23], flag = 0, zero[2];
	u32 i, i1, r = 0, success = 0;
	struct ncastate state2[160];
	struct datalist *p1;
//	struct datalist2 *p2, *head, *p3, *p4, *p5;       // p2,p3,p4,p5 and head are used to store the second set in the self-contained method
	struct isd2 *pp1, *pp2;
	double sum, counter1;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(zero, 0, 2);
	memset(state2, 0, sizeof(struct ncastate) * 160);
	memset(inter, 0, 23);
	pp1 = NULL;
	pp2 = NULL;

	p1 = NULL;
	p1 = (struct datalist *)calloc(1, sizeof(struct datalist));
	if (p1 == NULL){
		printf("error in calloc in SetB_generation_Self_contained3\n");
	}

	for (j = 0; j < 2; j++)zero[j] = 0x1;

	i1 = 0;
	counter1 = 0;
	r = 0;

	//	while (success < 512){
//	while (success < 128){                                        //32768 = 2^19
	while (success < 1024){
//	while (success < 2048){
//	while (success < 2097152){                              //1048576 = 2^20
		randomIV(NFSR, 80);
		randomIV(LFSR, 80);
		Extractstate_BSW2(state2, LFSR, NFSR);           // generate the 0x1f keystream prefix
		/*for (j = 0; j < 40; j++)printf("%x", NFSR[j]);
		printf("\n");
		for (j = 0; j < 40; j++)printf("%x", LFSR[j]);
		printf("\n");*/
		Diff_back(state2, LFSR, NFSR);
		Keystreamgen(2, z, LFSR, NFSR);
		/*for (j = 0; j < 5; j++)printf("%x",z[j]);
		printf("\n");*/

		/*printf("The generated keystream is (should be 0x1f):");
		for (j = 0; j < 5; j++)printf("%x", z[j]);
		printf("\n");*/

		k = 0;
		j = 0;
		while (j < 23){
			if ((state2 + k)->flag == 1){
				p1->ISnca[j] = (state2 + k)->state;
				k++;
				j++;
			}
			else k++;
		}
		/*for (j = 0; j < 41; j++)printf("%x", p1->ISnca[j]);
		printf("\n");*/

		for (j = 0; j < 2; j++)p1->prefix[j] = z[j];
		for (sum = 0.0, j = 0; j < 23; j++)sum += (p1->ISnca[j]) * pow(2, j);
		p1->ISnca_sum = sum;
		p1->position = i1;
		i1++;
		for (j = 0; j < 23; j++)inter[j] = p1->ISnca[j];

//		for (i = 0; i < ListSIZE; i++){
		for (i = 0; i < 1; i++){
			/*printf("The %u th element is:", i);
			for (j = 0; j < 5; j++)printf("%x", (list + i)->prefix[j]);
			printf("\n");
			printf("The outside element is (should be 0x1f):");
			for (j = 0; j < 5; j++)printf("%x", p1->prefix[j]);
			printf("\n");*/
			for (j = 0; j < 2; j++)z1[j] = p1->prefix[j] ^ (list + i)->prefix[j];
			for (index = 0, j = 0; j < 2; j++)index += (u8)((z1[j]) * pow(2, j));
			/*printf("The KSD is %u \n", index);*/

			/*printf("The set A is:\n");
			for (j = 0; j < 35; j++)printf("%x", (list + i)->ISnca[j]);
			printf("\n");*/
			pp1 = isd2_table[index][3];                             // isd2_table is used to store the precomputation information
			pp2 = pp1;
			do
			{
				for (j = 0; j < 23; j++)inter[j] = p1->ISnca[j];
				/*printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
				printf("The set B before fliping is:\n");
				for (j = 0; j < 35; j++)printf("%x", inter[j]);
				printf("\n");
				printf("The filped positions are: ");
				printf("%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2]);*/
				flip_inter(inter, pp2);
				/*printf("The set B after fliping is:\n");
				for (j = 0; j < 35; j++)printf("%x", inter[j]);
				printf("\n");
				printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");*/
				/*inter[pp2->position[0]] ^= 0x1;
				inter[pp2->position[1]] ^= 0x1;
				inter[pp2->position[2]] ^= 0x1;*/

				flag = verify_keystr(list + i, inter);
				/*printf("The flag is %u\n", flag);*/
				if (flag == 0){
					r++;
					for (sum = 0.0, j = 0; j < 23; j++)sum += inter[j] * pow(2, j);
					try[(u32)sum] ++;
					/*printf("%d\n", memcmp((list + i)->ISnca, inter, 10));
					printf("%d\n", memcmp(((list + i)->ISnca) + 15, inter+15, 21) );*/
//					printf("-----------------------------------------------\n");
					if (r == 1){
						/*for (sum = 0.0, j = 0; j < 28; j++)sum += inter[j] * pow(2, j);
						try[r - 1] = (float)sum; */
						/*head = p2;
						for (j = 0; j < 3; j++)p2->prefix[j] = (list + i)->prefix[j];
						for (j = 0; j < 28; j++)p2->ISnca[j] = inter[j];
						for (sum = 0.0, j = 0; j < 28; j++)sum += (p2->ISnca[j]) * pow(2, j);
						p2->ISnca_sum = sum;
						p2->position = r - 1;
						p2->next = NULL;
						p4 = p2;*/

//						if (memcmp((list + i)->ISnca, inter, 35) == 0){
						if ((memcmp((list + i)->ISnca, inter, 7) == 0) && (memcmp(((list + i)->ISnca) + 9, inter + 9, 14) == 0)) {

							printf("Found!!!\n");
						/*	printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
//						else break;						

					}
					else
					{
						/*for (sum = 0.0, j = 0; j < 28; j++)sum += inter[j] * pow(2, j);
						try[r - 1] = (float)sum;*/

						/*for (j = 0; j < 3; j++)p5->prefix[j] = (list + i)->prefix[j];
						for (j = 0; j < 28; j++)p5->ISnca[j] = inter[j];
						for (sum = 0.0, j = 0; j < 28; j++)sum += (p5->ISnca[j]) * pow(2, j);
						p5->ISnca_sum = sum;
						p5->position = r - 1;

						p5->next = NULL;
						p4->next = p5;
						p4 = p5;*/

//						if (memcmp((list + i)->ISnca, inter, 41) == 0){
						if ((memcmp((list + i)->ISnca, inter, 7) == 0) && (memcmp(((list + i)->ISnca) + 9, inter + 9, 14) == 0)){

							printf("Found!!!\n");
							/*printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
//						else break;
					}
				}
				pp2 = pp2->next;
			} while (pp2 != NULL);
//			printf("**********************************************\n");
		}
		success++;
	}   // corresponding to while(successs<)

//	printf("The number of candidates is %u\n", r);
//	printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");

	return(try);
}

void Self_contained2_counting(u8 *try1)
{
	u32 i, num;

	for (num = 0, i = 0; i < 8388608; i++)
	if (try1[i] != 0)
	{
		num++;
	}
	printf("There are %u elements after one self_contained invoking\n", num);
}

//void Self_contained2_ncacounting(u8 *try1)
//{
//	u32 i, num;
//
//	for (num = 0, i = 0; i < 4096; i++)
//	if (try1[i] != 0)
//	{
//		num++;
//	}
////	printf("num=%u\n", num);
//}

void Self_contained2_ncacounting(u16 *try1)
{
	u32 i, num;

	for (num = 0, i = 0; i < 256; i++)
	if (try1[i] != 0)
	{
		num++;
	}
	//	printf("num=%u\n", num);
}

void Self_contained2_ncacounting4(u16 *try1)
{
	u32 i, num;

	for (num = 0, i = 0; i < 256; i++)
	if (try1[i] != 0)
	{
		num++;
	}
	//	printf("num=%u\n", num);
}

u32 Self_contained2_ncacounting2(u16 *try1)
{
	u32 i, num;

	for (num = 0, i = 0; i < 256; i++)
	if (try1[i] != 0)
	{
		num++;
	}
	return(num);
	//	printf("num=%u\n", num);
}

u32 Self_contained2_ncacounting4a(u16 *try1)
{
	u32 i, num;

	for (num = 0, i = 0; i < 256; i++)
	if (try1[i] != 0)
	{
		num++;
	}
	return(num);
	//	printf("num=%u\n", num);
}

void Self_contained2_counting2(u8 *try1)
{
	u32 i, num;

	for (num = 0, i = 0; i < 2097152; i++)
	if (try1[i] != 0)
	{
		num++;
	}
	printf("There are %u elements after one self_contained invoking\n", num);
}

//void Self_contained2_counting2(u8 *try1)
//{
//	u32 i, num;
//
//	for (num = 0, i = 0; i < 2097152; i++)
//	if (try1[i] != 0)
//	{
//		num++;
//	}
//	printf("There are %u elements after one self_contained invoking\n", num);
//}

void Merge(u8 *try, u8 *result)
{
	u32 i;

	for (i = 0; i < 2097152; i++)
	{
		if (try[i] != 0){

			result[i] = 1;

		}
		else continue;
	}
}

//void Merge_nca(u8 *try, u8 *result)
//{
//	u32 i;
//
//	for (i = 0; i < 4096; i++)
//	{
//		if (try[i] != 0){
//
//			result[i] = 1;
//
//		}
//		else continue;
//	}
//}

void Merge_nca(u16 *try, u16 *result)
{
	u32 i;

	for (i = 0; i < 256; i++)
	{
		if (try[i] != 0){

			result[i] = 1;

		}
		else continue;
	}
}

void Counting_merge(u8 *result)
{
	u32 i, num;

	for (num = 0, i = 0; i < 2097152; i++)
	if (result[i] != 0)
	{
		num++;
	}
	printf("There are %u elements after merging\n", num);
}

//u32 Counting_merge_nca(u8 *result)
//{
//	u32 i, num;
//
//	for (num = 0, i = 0; i < 4096; i++)
//	if (result[i] != 0)
//	{
//		num++;
//	}
////	printf("There are %u elements after merging nca\n", num);
//
//	return(num);
//}

u32 Counting_merge_nca(u16 *result)
{
	u32 i, num;

	for (num = 0, i = 0; i < 256; i++)
	if (result[i] != 0)
	{
		num++;
	}
	//	printf("There are %u elements after merging nca\n", num);

	return(num);
}

void Self_contained2_aftercounting(u8 *try1)
{
	u32 i, num;

	for (num = 0, i = 0; i < 8388608; i++)
	if (try1[i] != 0)
	{
		num++;
	}
	printf("There are %u elements after one intersection\n", num);
}


void Self_contained2_intersection(u8 *try1, u8 *try2)          // generate the intersection set of the candidates 
{
	u32 i, j, num;

	for (i = 0; i < 8388608; i++)
	if ((try1[i] != 0) && (try2[i] != 0))               // note that here the logical operation is &&, not ||
	{
//		try1[i] = try1[i] + try2[i];
		try1[i] = 1;
	}
	else
	{
		try1[i] = 0;
	}

	free(try2);

	for (num = 0, i = 0; i < 8388608; i++)
	if (try1[i] != 0)num++;

	printf("There are %u comment elements\n", num);
}

void Self_contained2_intersection2(u8 *try1, u8 *try2)      // generate the intersection set of the candidates 
{
	u32 i, j, num;

	for (i = 0; i < 8388608; i++)
	if ((try1[i] != 0) && (try2[i] != 0))               // note that here the logical operation is &&, not ||
	{
//		try1[i] = (u16)(try1[i] + try2[i]);
		try1[i] = 1;
	}
	else
	{
		try1[i] = 0;
	}

	free(try2);

	for (num = 0, i = 0; i < 8388608; i++)
	if (try1[i] != 0)num++;

	printf("There are %u comment elements\n", num);
}

void Self_contained2_intersection_new(u8 *try1, u8 *try2)      // generate the intersection set of the candidates 
{
	u32 i, j, num;

	for (i = 0; i < 8388608; i++)
	if ((try1[i] != 0) && (try2[i] != 0))               // note that here the logical operation is &&, not ||
	{
//		try1[i] = (u16)(try1[i] + try2[i]);
		try1[i] = 1;
	}
	else
	{
		try1[i] = 0;
	}

//	free(try2);

	for (num = 0, i = 0; i < 8388608; i++)
	if (try1[i] != 0)num++;

	printf("There are %u common elements\n", num);
}

//void Self_contained2_intersection_nca(u8 *try1, u8 *try2)      // generate the intersection set of the candidates 
//{
//	u32 i, j, num;
//
//	for (i = 0; i < 4096; i++)
//	if ((try1[i] != 0) && (try2[i] != 0))               // note that here the logical operation is &&, not ||
//	{
//		try1[i] = 1;
//	}
//	else
//	{
//		try1[i] = 0;
//	}
//
//
//	for (num = 0, i = 0; i < 4096; i++)
//	if (try1[i] != 0)num++;
//
////	printf("There are %u common elements in new\n", num);
//	printf("common %u in new\n", num);
//}


void Self_contained2_intersection_nca(u16 *try1, u16 *try2)      // generate the intersection set of the candidates 
{
	u32 i, j, num;

	for (i = 0; i < 256; i++)
	if ((try1[i] != 0) && (try2[i] != 0))               // note that here the logical operation is &&, not ||
	{
		try1[i] = 1;
	}
	else
	{
		try1[i] = 0;
	}


	for (num = 0, i = 0; i < 256; i++)
	if (try1[i] != 0)num++;

	//	printf("There are %u common elements in new\n", num);
//	printf("common %u in new\n", num);
}

void Self_contained2_intersection_nca4(u16 *try1, u16 *try2)      // generate the intersection set of the candidates 
{
	u32 i, j, num;

	for (i = 0; i < 256; i++)
	if ((try1[i] != 0) && (try2[i] != 0))               // note that here the logical operation is &&, not ||
	{
		try1[i] = 1;
	}
	else
	{
		try1[i] = 0;
	}


	for (num = 0, i = 0; i < 256; i++)
	if (try1[i] != 0)num++;

	//	printf("There are %u common elements in new\n", num);
	printf("common %u in new\n", num);
}


void Self_contained2_intersection_new2(u8 *try1, u8 *try2)      // generate the intersection set of the candidates 
{
	u32 i, j, num;

	for (i = 0; i < 2097152; i++)
	if ((try1[i] != 0) && (try2[i] != 0))               // note that here the logical operation is &&, not ||
	{
		try1[i] = 1;
	}
	else
	{
		try1[i] = 0;
	}

	
	for (num = 0, i = 0; i < 2097152; i++)
	if (try1[i] != 0)num++;

	printf("There are %u common elements\n", num);
}

//void Self_contained2_merging(u16 *try1, u16 *try2)          // generate the merging set of the candidates 
//{
//	u32 i, j, num;
//
//	for (i = 0; i < 8388608; i++)
//	if ((try1[i] != 0) || (try2[i] != 0))
//	{
//		try1[i] = (u16)(try1[i] + try2[i]);
//		try2[i] = 0;
//	}
//	/*else
//	{
//	try1[i] = 0;
//	}*/
//
//	free(try2);
//
//	for (num = 0, i = 0; i < 524288; i++)
//	if (try1[i] != 0)num++;
//
//	printf("There are %u merging elements\n", num);
//}

u8 check_candidate(u8 *try, double in)
{
	if (try[(u32)in] != 0)return(1);
	else return(0);
}

//u8 check_candidate1(u8 *try, double in)
//{
//	if (try[(u32)in] != 0)return(1);
//	else return(0);
//}

u8 check_candidate1(u16 *try, double in)
{
	if (try[(u32)in] != 0)return(1);
	else return(0);
}

//void counting_score(u16 *try1, double in)
//{
//	u32 i, num;
//
//	printf("The score of the correct element is %u\n", try1[(u32)in]);
//
//	num = 0;
//	for (i = 0; i < 524288; i++)
//	{
//		if (try1[i] >= try1[(u32)in])num++;
//	}
//
//	printf("The number of elements that have a higher score is %u\n", num);
//}


struct initialdata2 *Targetkeystream_generation(u8 *LFSR, u8 *NFSR, struct initialdata2 *list)        // generate the target keystream and its corresponding internal state
{
	u8 z[100], j, k, lfsr[40], nfsr[40];
	u32 i;
	struct ncastate state2[80];
	double sum;

	memset(z, 0, 100);
	memset(lfsr, 0, 40);
	memset(nfsr, 0, 40);
	memset(state2, 0, sizeof(struct ncastate) * 80);


	randomIV(LFSR, 40);
	randomIV(NFSR, 40);

	for (i = 0; i < 40; i++)list->lfsr[i] = LFSR[i];
	for (i = 0; i < 40; i++)list->nfsr[i] = NFSR[i];

	Keystreamgen(100, z, LFSR, NFSR);

	for (i = 0; i < 100; i++)list->keystream[i] = z[i];

	return(list);
}


struct initialdata2 *Targetkeystream_generation_inner(u8 *LFSR, u8 *NFSR, struct initialdata2 *list,u8 pattern[4])        // generate the target keystream and its corresponding internal state
{
	u8 z[100], j, k, lfsr[80], nfsr[80];
	u32 i;
	struct ncastate state2[160];
	double sum;
	u8 x0 = 0, x1 = 0, x2 = 0, x3 = 0, n63 = 0, n64 = 0, n65 = 0, n66 = 0;

	memset(z, 0, 100);
	memset(lfsr, 0, 80);
	memset(nfsr, 0, 80);
	memset(state2, 0, sizeof(struct ncastate) * 160);


	randomIV(LFSR, 80);
	randomIV(NFSR, 80);

	LFSR[3] = pattern[0] & 0x1;
	LFSR[4] = pattern[1] & 0x1;
	LFSR[5] = pattern[2] & 0x1;
	LFSR[6] = pattern[3] & 0x1;
	
	LFSR[25] = (pattern[0] >> 1) & 0x1;
	LFSR[26] = (pattern[1] >> 1) & 0x1;
	LFSR[27] = (pattern[2] >> 1) & 0x1;
	LFSR[28] = (pattern[3] >> 1) & 0x1;

	LFSR[46] = (pattern[0] >> 2) & 0x1;
	LFSR[47] = (pattern[1] >> 2) & 0x1;
	LFSR[48] = (pattern[2] >> 2) & 0x1;
	LFSR[49] = (pattern[3] >> 2) & 0x1;

	LFSR[64] = (pattern[0] >> 3) & 0x1;
	LFSR[65] = (pattern[1] >> 3) & 0x1;
	LFSR[66] = (pattern[2] >> 3) & 0x1;
	LFSR[67] = (pattern[3] >> 3) & 0x1;
	
	for (i = 0; i < 80; i++)list->lfsr[i] = LFSR[i];
	for (i = 0; i < 80; i++)list->nfsr[i] = NFSR[i];
	
	x0 = list->nfsr[1] ^ list->nfsr[2] ^ list->nfsr[4] ^ list->nfsr[10] ^ list->nfsr[31] ^ list->nfsr[43] ^ list->nfsr[56];
	x1 = list->nfsr[2] ^ list->nfsr[3] ^ list->nfsr[5] ^ list->nfsr[11] ^ list->nfsr[32] ^ list->nfsr[44] ^ list->nfsr[57];
	x2 = list->nfsr[3] ^ list->nfsr[4] ^ list->nfsr[6] ^ list->nfsr[12] ^ list->nfsr[33] ^ list->nfsr[45] ^ list->nfsr[58];
	x3 = list->nfsr[4] ^ list->nfsr[5] ^ list->nfsr[7] ^ list->nfsr[13] ^ list->nfsr[34] ^ list->nfsr[46] ^ list->nfsr[59];
	
	n63 = list->nfsr[63];
	n64 = list->nfsr[64];
	n65 = list->nfsr[65];
	n66 = list->nfsr[66];
	
	Keystreamgen(100, z, LFSR, NFSR);

	for (i = 0; i < 100; i++)list->keystream[i] = z[i];
	printf("The correct state is:");
	printf("x0=%x,x1=%x,x2=%x,x3=%x,n63=%x,n64=%x,n65=%x,n66=%x\n", x0,x1,x2,x3,n63,n64,n65,n66);
	return(list);
}

struct initialdata2 *Targetkeystream_generation2(u8 *LFSR, u8 *NFSR, struct initialdata2 *list, u8 prefix[2])        // generate the target keystream and its corresponding internal state
{
	u8 z[100], j, k, lfsr[80], nfsr[80];
	u32 i;
	struct ncastate state2[160];
	double sum;

	memset(z, 0, 100);
	memset(lfsr, 0, 80);
	memset(nfsr, 0, 80);
	memset(state2, 0, sizeof(struct ncastate) * 160);


	Extract_state(LFSR,NFSR,prefix);

	for (i = 0; i < 80; i++)list->lfsr[i] = LFSR[i];
	for (i = 0; i < 80; i++)list->nfsr[i] = NFSR[i];

	Keystreamgen(100, z, LFSR, NFSR);

	for (i = 0; i < 100; i++)list->keystream[i] = z[i];

	return(list);
}

void output_initialdata2(struct initialdata2 *list)
{
	u32 i;

	printf("The lfsr internal state is: ");
	for (i = 0; i < 40; i++)printf("%x", list->lfsr[i]);
	printf("\n");
	printf("The nfsr internal state is: ");
	for (i = 0; i < 40; i++)printf("%x", list->nfsr[i]);
	printf("\n");
	printf("The generated 100 keystream is: ");
	for (i = 0; i < 100; i++)printf("%x", list->keystream[i]);
	printf("\n");
}

void Prepare_data(struct initialdata2 *initiallist, struct datalist2 *list, u8 n)           // transform into the data structure used in our programm
{
	u32 i;
	struct ncastate state2[160];
	u8 z[2], z1[100], lfsr[80], nfsr[80],k,j;
	double sum;

	memset(z, 0, 2);
	memset(z1, 0, 100);
	memset(lfsr, 0, 80);
	memset(nfsr, 0, 80);
	memset(state2, 0, sizeof(struct ncastate) * 160);

	for (i = 0; i < 2; i++)list->prefix[i] = initiallist->keystream[n + i];
	list->position = n;

	for (i = 0; i < 2; i++)z[i] = initiallist->keystream[n + i];
	for (i = 0; i < 80; i++)lfsr[i] = initiallist->lfsr[i];
	for (i = 0; i < 80; i++)nfsr[i] = initiallist->nfsr[i];

	if (n != 0)Keystreamgen(n, z1, lfsr, nfsr);                            // run the beginning state up to the correct position
	Extractstate_BSW_arbitrary(state2, lfsr, nfsr, z);
	
	k = 0;
	j = 0;
	while (j < 23){
		if ((state2 + k)->flag == 1){
			list->ISnca[j] = (state2 + k)->state;
			k++;
			j++;
		}
		else k++;
	}

	for (sum = 0.0, j = 0; j < 23; j++)sum += (list->ISnca[j]) * pow(2, j);
	list->ISnca_sum = sum;
	list->next = NULL;
}

void Prepare_ncadata(struct initialdata2 *initiallist, struct datalist2 *list, u8 n)           // transform into the data structure used in our programm
{
	u32 i,j;
	u8 ns[8];
	u8 z[2], z1[100], lfsr[40], nfsr[40];
	double sum;

	memset(z, 0, 2);
	memset(z1, 0, 100);
	memset(ns, 0, 8);
	memset(lfsr, 0, 40);
	memset(nfsr, 0, 40);
	

	for (i = 0; i < 2; i++)list->prefix[i] = initiallist->keystream[n + i];
	list->position = n;

	for (i = 0; i < 2; i++)z[i] = initiallist->keystream[n + i];
	for (i = 0; i < 40; i++)lfsr[i] = initiallist->lfsr[i];
	for (i = 0; i < 40; i++)nfsr[i] = initiallist->nfsr[i];

	if (n != 0){
		Keystreamgen(n, z1, lfsr, nfsr);                            // run the beginning state up to the correct position
		Extractnca_arbitrary(ns, lfsr, nfsr, z);
	}
	else{
		Extractnca_arbitrary(ns, lfsr, nfsr, z);
	}

	for (j = 0; j < 8; j++)list->ISnca[j] = ns[j];
	
	for (sum = 0.0, j = 0; j < 8; j++)sum += (list->ISnca[j]) * pow(2, j);
	list->ISnca_sum = sum;
	list->next = NULL;
}

void Prepare_inner_ncadata(struct initialdata2 *initiallist, struct datalist2_inner *list, u8 n,u8 pattern[4])           // transform into the data structure used in our programm
{
	u32 i, j;
	u8 ns[24];
	u8 z[4], z1[100], lfsr[80], nfsr[80];
	double sum;

	memset(z, 0, 4);
	memset(z1, 0, 100);
	memset(ns, 0, 24);
	memset(lfsr, 0, 80);
	memset(nfsr, 0, 80);


	for (i = 0; i < 4; i++)list->prefix[i] = initiallist->keystream[n + i];
	list->position = n;

	for (i = 0; i < 4; i++)z[i] = initiallist->keystream[n + i];
	for (i = 0; i < 80; i++)lfsr[i] = initiallist->lfsr[i];
	for (i = 0; i < 80; i++)nfsr[i] = initiallist->nfsr[i];

	if (n != 0){
		Keystreamgen(n, z1, lfsr, nfsr);                            // run the beginning state up to the correct position
		Extractnca_arbitrary_inner4(ns, lfsr, nfsr, z,pattern);
	}
	else{
		Extractnca_arbitrary_inner4(ns, lfsr, nfsr, z,pattern);
	}

	list->ISnca[0] = ns[0];
	list->ISnca[1] = ns[1];
	list->ISnca[2] = ns[2];
	list->ISnca[3] = ns[3];

	list->ISnca[4] = ns[20];
	list->ISnca[5] = ns[21];
	list->ISnca[6] = ns[22];
	list->ISnca[7] = ns[23];

	for (sum = 0.0, j = 0; j < 8; j++)sum += (list->ISnca[j]) * pow(2, j);
	list->ISnca_sum = sum;
	list->next = NULL;
}



u8 *ISrecovery_Self_contained(u8 *LFSR, u8 *NFSR, struct datalist2 *list, struct isd2 *isd2_table[4][4], u8 *try)     // generate the internal state corresponding to the 2-bit keystream segment
{
	u8 z[2], j, k, z1[2], index = 0, inter[23], flag = 0, zero[2];
	u32 i, i1, r = 0, success = 0;
	struct ncastate state2[160];
	struct datalist *p1;
	struct isd2 *pp1, *pp2;
	double sum, counter1;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(zero, 0, 2);
	memset(state2, 0, sizeof(struct ncastate) * 160);
	memset(inter, 0, 23);
	pp1 = NULL;
	pp2 = NULL;

	p1 = NULL;
	p1 = (struct datalist *)calloc(1, sizeof(struct datalist));
	if (p1 == NULL){
		printf("error in calloc in ISrecovery_Self_contained\n");
	}

	for (j = 0; j < 2; j++)zero[j] = 0x1 ^ list->prefix[j];

	i1 = 0;
	counter1 = 0;
	r = 0;

	while (success < 1024){
//	while (success < 1100){
//	while (success < 2048){
		randomIV(NFSR, 80);
		randomIV(LFSR, 80);
		Extractstate_BSW_arbitrary(state2, LFSR, NFSR,zero);       // generate the desirable keystream prefix
		/*for (j = 0; j < 40; j++)printf("%x", NFSR[j]);
		printf("\n");
		for (j = 0; j < 40; j++)printf("%x", LFSR[j]);
		printf("\n");*/
		Diff_back(state2, LFSR, NFSR);
		Keystreamgen(2, z, LFSR, NFSR);

		k = 0;
		j = 0;
		while (j < 23){
			if ((state2 + k)->flag == 1){
				p1->ISnca[j] = (state2 + k)->state;
				k++;
				j++;
			}
			else k++;
		}
	
		for (j = 0; j < 2; j++)p1->prefix[j] = z[j];
		for (sum = 0.0, j = 0; j < 23; j++)sum += (p1->ISnca[j]) * pow(2, j);
		p1->ISnca_sum = sum;
		p1->position = 0;
		
		for (j = 0; j < 23; j++)inter[j] = p1->ISnca[j];

		for (i = 0; i < 1; i++){
			
			for (j = 0; j < 2; j++)z1[j] = (list + i)->prefix[j];
			for (index = 0, j = 0; j < 2; j++)index += (u8)((z1[j]) * pow(2, j));
			/*printf("The actual prefix is %u \n", index);*/

			/*printf("The set A is:\n");
			for (j = 0; j < 35; j++)printf("%x", (list + i)->ISnca[j]);
			printf("\n");*/
//			pp1 = isd2_table[index][3];                             // isd2_table is used to store the precomputation information
			pp1 = isd2_table[3][index];
			pp2 = pp1;
			do
			{
				for (j = 0; j < 23; j++)inter[j] = p1->ISnca[j];
				/*printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
				printf("The set B before fliping is:\n");
				for (j = 0; j < 35; j++)printf("%x", inter[j]);
				printf("\n");
				printf("The filped positions are: ");
				printf("%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2]);*/
				flip_inter(inter, pp2);
				/*printf("The set B after fliping is:\n");
				for (j = 0; j < 35; j++)printf("%x", inter[j]);
				printf("\n");
				printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");*/
				/*inter[pp2->position[0]] ^= 0x1;
				inter[pp2->position[1]] ^= 0x1;
				inter[pp2->position[2]] ^= 0x1;*/

				flag = verify_keystr(list + i, inter);
				/*printf("The flag is %u\n", flag);*/
				if (flag == 0){
					r++;
					for (sum = 0.0, j = 0; j < 23; j++)sum += inter[j] * pow(2, j);
					try[(u32)sum] ++;
					
					if (r == 1){
						printf("success=%u\n", success);
						if ((memcmp((list + i)->ISnca, inter, 7) == 0) && (memcmp(((list + i)->ISnca) + 9, inter + 9, 14) == 0)) {

							printf("Found!!!\n");
							/*	printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
//						else break;						

					}
					else
					{
						printf("success=%u\n", success);
						if ((memcmp((list + i)->ISnca, inter, 7) == 0) && (memcmp(((list + i)->ISnca) + 9, inter + 9, 14) == 0)){

							printf("Found!!!\n");
							/*printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
//						else break;
					}
				}
				pp2 = pp2->next;
			} while (pp2 != NULL);
//			printf("**********************************************\n");
		}
		success++;
	}   // corresponding to while(successs<)

	return(try);
}

u16 *ISrecovery_nca_Selfcontained(u8 *LFSR, u8 *NFSR, struct datalist2 *list, struct isd2 *isd2_table[4][4], u16 *try)     // generate the internal state corresponding to the 2-bit keystream segment
{
	u8 z[2], j, k, z1[2], index = 0, inter[8], flag = 0, zero[2];
	u32 i, i1, r = 0, success = 0, c = 0;
	u8 ns[8];
	struct datalist2 *p1;
	struct isd2 *pp1, *pp2;
	double sum, counter1;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(zero, 0, 2);
	memset(ns, 0, 8);
	memset(inter, 0, 8);
	pp1 = NULL;
	pp2 = NULL;

	p1 = NULL;
	p1 = (struct datalist2 *)calloc(1, sizeof(struct datalist2));
	if (p1 == NULL){
		printf("error in calloc in ISrecovery_nca_Selfcontained\n");
	}

//	for (j = 0; j < 2; j++)zero[j] = 0x1 ^ list->prefix[j];
	for (j = 0; j < 2; j++)zero[j] = 0x0 ^ list->prefix[j];

	i1 = 0;
	counter1 = 0;
	r = 0;
//	c = 5;
	c = 8;

	while (success < c*(NCAC)){

		randomIV(NFSR, 40);
		randomIV(LFSR, 40);
		Extractnca_arbitrary(ns, LFSR, NFSR, zero);       // generate the desirable keystream prefix
		Keystreamgen_nca2(z, ns);

		for (j = 0; j < 8; j++)p1->ISnca[j] = ns[j];
		
		for (j = 0; j < 2; j++)p1->prefix[j] = z[j];
		for (sum = 0.0, j = 0; j < 8; j++)sum += (p1->ISnca[j]) * pow(2, j);
		p1->ISnca_sum = sum;
		p1->position = 0;

		for (j = 0; j < 8; j++)inter[j] = p1->ISnca[j];

		for (i = 0; i < 1; i++){

			for (j = 0; j < 2; j++)z1[j] = (list + i)->prefix[j];
			for (index = 0, j = 0; j < 2; j++)index += (u8)((z1[j]) * pow(2, j));
			
//			pp1 = isd2_table[index][3];                             // isd2_table is used to store the precomputation information
//			pp1 = isd2_table[3][index];
			pp1 = isd2_table[0][index];
			pp2 = pp1;
			do
			{
				for (j = 0; j < 8; j++)inter[j] = p1->ISnca[j];
				flip_ncainter(inter, pp2);
				flag = verify_ncakeystr(list + i, inter);

				if (flag == 0){
					r++;
					for (sum = 0.0, j = 0; j < 8; j++)sum += inter[j] * pow(2, j);
					try[(u32)sum] ++;

					if (r == 1){
//						if ( memcmp((list + i)->ISnca + 2, inter + 2 , 10) == 0 ) {
						if (memcmp((list + i)->ISnca, inter, 8) == 0) {

//							printf("Fd!\n");
							/*	printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
//						else break;						

					}
					else
					{
//						if ( memcmp((list + i)->ISnca + 2, inter + 2, 10) == 0 ){
						if (memcmp((list + i)->ISnca, inter, 8) == 0) {

//							printf("Fd!\n");
							/*printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
//						else break;
					}
				}
				pp2 = pp2->next;
			} while (pp2 != NULL);
			
		}
		success++;
	}   // corresponding to while(successs<)

	return(try);
}


u16 *ISrecovery_innernca_Selfcontained(u8 *LFSR, u8 *NFSR, struct datalist2 *list, struct isd3 *isd2_table[4][4], u16 *try,u8 pattern[2])     // inner part function
{
	u8 z[2], j, k, z1[2], index = 0, inter[12], flag = 0, zero[2];
	u32 i, i1, r = 0, success = 0, c = 0;
	u8 ns[12];
	struct datalist2 *p1;
	struct isd3 *pp1, *pp2;
	double sum, counter1;

	memset(z, 0, 2);
	memset(z1, 0, 2);
	memset(zero, 0, 2);
	memset(ns, 0, 12);
	memset(inter, 0, 12);
	pp1 = NULL;
	pp2 = NULL;

	p1 = NULL;
	p1 = (struct datalist2 *)calloc(1, sizeof(struct datalist2));
	if (p1 == NULL){
		printf("error in calloc in ISrecovery_nca_Selfcontained\n");
	}

//	for (j = 0; j < 2; j++)zero[j] = 0x1 ^ list->prefix[j];
	flag = 0;
	for (i = 0; i < 4;i++)
	for (j = 0; j < 4;j++){
		if (isd2_table[i][j] != NULL){
			for (k = 0; k < 2; k++)zero[k] = ((i >> k) & 0x1) ^ list->prefix[k];
			i1 = i;
//			printf("i1=%u\n",i1);
			flag = 1;
			break;
		}
		if (flag == 1)break;
	}
//	printf("to 1\n");

	counter1 = 0;
	r = 0;
//	c = 5;
	c = 3;

	while (success < c*(INCAC)){

		randomIV(NFSR, 80);
		randomIV(LFSR, 80);
		Extractnca_arbitrary_inner(ns, LFSR, NFSR, zero, pattern);       // generate the desirable keystream prefix
		Keystreamgen_nca2(z, ns);
//		printf("to 2\n");
		for (j = 0; j < 12; j++)p1->ISnca[j] = ns[j];

		for (j = 0; j < 2; j++)p1->prefix[j] = z[j];
		for (sum = 0.0, j = 0; j < 12; j++)sum += (p1->ISnca[j]) * pow(2, j);
		p1->ISnca_sum = sum;
		p1->position = 0;

		for (j = 0; j < 12; j++)inter[j] = p1->ISnca[j];

		for (i = 0; i < 1; i++){

			for (j = 0; j < 2; j++)z1[j] = (list + i)->prefix[j];
			for (index = 0, j = 0; j < 2; j++)index += (u8)((z1[j]) * pow(2, j));
//			printf("to 2a\n");
//			pp1 = isd2_table[index][3];                             // isd2_table is used to store the precomputation information
//			pp1 = isd2_table[3][index];
			pp1 = isd2_table[i1][index];
//			printf("to 2d\n");
			pp2 = pp1;
			do
			{
//				printf("to 2e\n");
				for (j = 0; j < 12; j++)inter[j] = p1->ISnca[j];
				flip_ncainter_inner(inter, pp2);
//				printf("to 2f\n");
				flag = verify_ncakeystr(list + i, inter);
//				printf("to 3\n");
				if (flag == 0){
					r++;
					for (sum = 0.0, j = 0; j < 12; j++)sum += inter[j] * pow(2, j);
					try[(u32)sum] ++;

					if (r == 1){
//						if ( memcmp((list + i)->ISnca + 2, inter + 2 , 10) == 0 ) {
						if (memcmp((list + i)->ISnca, inter, 12) == 0) {

//							printf("Fd!\n");
							/*	printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
//						else break;						

					}
					else
					{
//						if ( memcmp((list + i)->ISnca + 2, inter + 2, 10) == 0 ){
						if (memcmp((list + i)->ISnca, inter, 12) == 0) {

//							printf("Fd!\n");
							/*printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
//						else break;
					}
				}
				pp2 = pp2->next;
			} while (pp2 != NULL);

		}
		success++;
	}   // corresponding to while(successs<)

//	printf("r=%u\n",r);
	free(p1);
	return(try);
}

u16 *ISrecovery_innernca_Selfcontained4(u8 *LFSR, u8 *NFSR, struct datalist2_inner *list, struct isd4 *isd2_table[16][16], u16 *try, u8 pattern[4])     // inner part function
{
	u8 z[4], j, k, z1[4], index = 0, inter[24], flag = 0, zero[4];
	u32 i, i1, r = 0, success = 0, c = 0;
	u8 ns[24];
	struct datalist2_inner *p1;
	struct isd4 *pp1, *pp2;
	double sum, counter1;

	memset(z, 0, 4);
	memset(z1, 0, 4);
	memset(zero, 0, 4);
	memset(ns, 0, 24);
	memset(inter, 0, 24);
	pp1 = NULL;
	pp2 = NULL;

	p1 = NULL;
	p1 = (struct datalist2_inner *)calloc(1, sizeof(struct datalist2_inner));
	if (p1 == NULL){
		printf("error in calloc in ISrecovery_nca_Selfcontained4\n");
	}

	//	for (j = 0; j < 2; j++)zero[j] = 0x1 ^ list->prefix[j];
	flag = 0;
	for (i = 0; i < 16; i++)
	for (j = 0; j < 16; j++){
		if (isd2_table[i][j] != NULL){
			for (k = 0; k < 4; k++)zero[k] = ((i >> k) & 0x1) ^ list->prefix[k];
			i1 = i;
//			printf("i1=%u\n",i1);
			flag = 1;
			break;
		}
		if (flag == 1)break;
	}
//	printf("to 1\n");

	counter1 = 0;
	r = 0;
//	c = 5;
	c = 14;

	while (success < c*(INCAC4)){

		randomIV(NFSR, 80);
		randomIV(LFSR, 80);
		Extractnca_arbitrary_inner4(ns, LFSR, NFSR, zero, pattern);       // generate the desirable keystream prefix
		Keystreamgen_nca6(z, ns);
//		printf("to 2\n");
		for (j = 0; j < 24; j++)p1->ISnca[j] = ns[j];
		/*list->ISnca[0] = ns[0];
		list->ISnca[1] = ns[1];
		list->ISnca[2] = ns[2];
		list->ISnca[3] = ns[3];

		list->ISnca[4] = ns[20];
		list->ISnca[5] = ns[21];
		list->ISnca[6] = ns[22];
		list->ISnca[7] = ns[23];*/

		for (j = 0; j < 4; j++)p1->prefix[j] = z[j];
		for (sum = 0.0, j = 0; j < 4; j++)sum += (p1->ISnca[j]) * pow(2, j);
		for (j = 20; j < 24; j++)sum += (p1->ISnca[j]) * pow(2, j-16);
		p1->ISnca_sum = sum;
		p1->position = 0;

		for (j = 0; j < 24; j++)inter[j] = p1->ISnca[j];

		for (i = 0; i < 1; i++){

			for (j = 0; j < 4; j++)z1[j] = (list + i)->prefix[j];
			for (index = 0, j = 0; j < 4; j++)index += (u8)((z1[j]) * pow(2, j));
//			printf("to 2a\n");
//			pp1 = isd2_table[index][3];                             // isd2_table is used to store the precomputation information
//			pp1 = isd2_table[3][index];
			pp1 = isd2_table[i1][index];
//			printf("to 2d\n");
			pp2 = pp1;
			do
			{
//				printf("to 2e\n");
				for (j = 0; j < 8; j++)inter[j] = p1->ISnca[j];
				flip_ncainter_inner4(inter, pp2);
//				printf("to 2f\n");
				flag = verify_ncakeystr4(list + i, inter);
//				printf("to 3\n");
				if (flag == 0){
					r++;
					for (sum = 0.0, j = 0; j < 4; j++)sum += inter[j] * pow(2, j);
					for (j = 20; j < 24; j++)sum += inter[j] * pow(2, j-16);
					try[(u32)sum] ++;

					if (r == 1){
//						if ( memcmp((list + i)->ISnca + 2, inter + 2 , 10) == 0 ) {
						if (memcmp((list + i)->ISnca, inter, 24) == 0) {

//							printf("Fd!\n");
							/*	printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
						//						else break;						

					}
					else
					{
//						if ( memcmp((list + i)->ISnca + 2, inter + 2, 10) == 0 ){
						if (memcmp((list + i)->ISnca, inter, 24) == 0) {

							//							printf("Fd!\n");
							/*printf("The flipped positions are: ");
							printf("%u,%u,%u,%u,%u\n", pp2->position[0], pp2->position[1], pp2->position[2], pp2->position[3], pp2->position[4]);
							for (j = 0; j < 2; j++)printf("%x", (list + i)->prefix[j]);
							printf("\n");
							for (j = 0; j < 23; j++)printf("%x", (list + i)->ISnca[j]);
							printf("\n");
							printf("The sum is %f\n", (list + i)->ISnca_sum);

							printf("The new position is %u and the original position is %u\n", r - 1, i);
							printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");*/
						}
						//						else break;
					}
				}
				pp2 = pp2->next;
			} while (pp2 != NULL);

		}
		success++;
	}   // corresponding to while(successs<)

//	printf("r=%u\n",r);
	free(p1);
	return(try);
}
#include "random.h"


void LFSR_clock(u8 *LFSR)
{
	u8 fd;
	u8 r;

	fd = 0x00;
	fd = LFSR[33] ^ LFSR[18] ^ LFSR[9] ^ LFSR[0];
	for(r=0; r<39; r++)LFSR[r] = LFSR[r+1];
	LFSR[39] = fd;
}


void NFSR_clock(u8 *NFSR, u8 LFSR_0)
{
	u8 fnd;
	u8 e;

	fnd = 0x00;
	fnd = LFSR_0 ^ NFSR[33] ^ NFSR[29] ^ NFSR[23] ^ NFSR[17] ^ NFSR[11] ^ NFSR[9] ^ ( NFSR[33] & NFSR[29] ) ^ ( NFSR[23] & NFSR[17] )
		^ ( NFSR[33] & NFSR[9] ) ^ ( NFSR[33] & NFSR[29] & NFSR[23]) ^ ( NFSR[29] & NFSR[23] & NFSR[17]) ^ ( NFSR[33] & NFSR[29] & NFSR[23] & NFSR[17]) ^ ( NFSR[29] & NFSR[23] & NFSR[17] & NFSR[11] & NFSR[9]);

	for(e=0; e<39; e++)NFSR[e] = NFSR[e+1];
	NFSR[39] = fnd;
}

u8 filter(u8 *LFSR,u8 *NFSR)
{
	u8 bot, out;

	bot = 0x00;
	out = 0x00;
	bot = LFSR[21] ^ (LFSR[1] & NFSR[22]) ^ (LFSR[21] & NFSR[22]) ^ (LFSR[1] & LFSR[21] & NFSR[22]);
	out = NFSR[1] ^ NFSR[7] ^ NFSR[15] ^ bot;
	return(out);

}

void LFSR_clockInmode(u8 *LFSR,u8 w)
{
	u8 fd;
	u8 r;

	fd = 0x00;
	fd = LFSR[33] ^ LFSR[18] ^ LFSR[9] ^ LFSR[0] ^ w;
	for(r=0; r<39; r++)LFSR[r] = LFSR[r+1];
	LFSR[39] = fd;
}


void NFSR_clockInmode(u8 *NFSR,u8 LFSR_0,u8 s)
{
	u8 fnd;
	u8 e;

	fnd = 0x00;
	fnd = s ^ LFSR_0 ^ NFSR[33] ^ NFSR[29] ^ NFSR[23] ^ NFSR[17] ^ NFSR[11] ^ NFSR[9] ^ (NFSR[33] & NFSR[29]) ^ (NFSR[23] & NFSR[17])
		^ (NFSR[33] & NFSR[9]) ^ (NFSR[33] & NFSR[29] & NFSR[23]) ^ (NFSR[29] & NFSR[23] & NFSR[17]) ^ (NFSR[33] & NFSR[29] & NFSR[23] & NFSR[17]) ^ (NFSR[29] & NFSR[23] & NFSR[17] & NFSR[11] & NFSR[9]);

	for(e=0; e<39; e++)NFSR[e] = NFSR[e+1];
	NFSR[39] = fnd;
}

void Initialization(u8 k[40], u8 iv[32], u8 *LFSR, u8 *NFSR)
{
	u32 t;
	u8 d;

	d = 0x00;

	memset(LFSR,0x00,40);
	memset(NFSR,0x00,40);
	memcpy(NFSR,k,40);
	memcpy(LFSR,iv,32);
	memset(LFSR+32,0x01,8);

//	for(t=0; t<40; t++)NFSR[t] = k[t];
//	for(t=0; t<32; t++)LFSR[t] = iv[t];
//	for(t=32; t<40; t++)LFSR[t] = 0x01;

	for(t=0; t<80; t++)
	{
		d = filter(LFSR,NFSR);
		NFSR_clockInmode(NFSR,LFSR[0],d);
		LFSR_clockInmode(LFSR,d);
		
	}
}


void Keystreamgen(u32 n, u8 *ks, u8 *LFSR, u8 *NFSR)
{
	u32 a;
	
	for(a=0; a<n; a++)
	{
		ks[a] = filter(LFSR,NFSR);
		NFSR_clock(NFSR,LFSR[0]);
		LFSR_clock(LFSR);
		
	}
}

u8 h_function(u8 x0, u8 x1, u8 x2)
{
	u8 bot;

	bot = 0x00;
	bot = x1 ^ (x0 & x2) ^ (x1 & x2) ^ (x0 & x1 & x2);
	return(bot);
}

u8 Nonlinear_masking(u8 *NFSR)
{
	u8 out;

	out = 0x00;
	out = NFSR[1] ^ NFSR[7] ^ NFSR[15];

	return(out);
}


/*  MT  */

//edited by Michael Trapp mt5918 CID: 01627245
//added content is in between */  MT */ and /*  (MT)  */

/*  (MT)  */

/* Copyright (C) 2000, Dr. Antonio Munjiza
 *
 * This code is provided as part of the book entitled "The Combined
 * Finite Discrete Element Method". It is distributed WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. Inclusion of a part or whole of this code into any other
 * commercial or research or other purpose code is not granted without author's
 * written explicit permission.
 * When results using whole or any part of this code
 * are published, Y code must be mentioned and acknowledgement to Dr Munjiza must be made.
 * Should you modify this source code, the Copyright (C) on the modified code
 * as a whole belongs to Dr. A. Munjiza regardless of the extent or nature
 * of modifications.
 * Copyright (C) to whole of any code containing any part of this code
 * also belongs to Dr. A.Munjiza.
 * Any code comprising any part of this source code
 * must be called Y program.
 * If you do not agree with this, you are not allowed to do
 * any modifications to any part of this source code or included
 * any part of it in any other program.
 */
 /* File  Yod.c */

#include "Yproto.h"

 /*  MT  */

void b64enc(FILE *fout, const UCHR *src, INS len)
{
	//encodes char array to base 64 and writes encoded chars to file

	//input
	//fout: output file
	//*src: input 1D char array containing ints or doubles
	//len:  length of input char array

	/*
	 * Base64 encoding/decoding (RFC1341)
	 * Copyright (c) 2005-2011, Jouni Malinen <j@w1.fi>
	 *
	 * This software may be distributed under the terms of the BSD license.
	 * See README for more details.
	 *
	*/
	const UCHR enc[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	const UCHR *end;
	UCHR *in;

	in = src;
	end = in + len;

	while (end - in > 2)
	{
		putc(enc[in[0] >> 2], fout);
		putc(enc[((in[0] & 0x03) << 4) | (in[1] >> 4)], fout);
		putc(enc[((in[1] & 0x0f) << 2) | (in[2] >> 6)], fout);
		putc(enc[in[2] & 0x3f], fout);
		in += 3;
	}

	if (end - in) 
	{
		putc(enc[in[0] >> 2], fout);
		if (end - in == 1)
		{
			putc(enc[(in[0] & 0x03) << 4], fout);
			putc('=', fout);
		}
		else
		{
			putc(enc[((in[0] & 0x03) << 4) | (in[1] >> 4)], fout);
			putc(enc[(in[1] & 0x0f) << 2], fout);
		}
		putc('=', fout);
	}
	return;
}

union INSToFLT {
	//converts int to float
	FLT f;
	INS i; //both i and f refer to the same location in memory
};

union INTToDBL {
	//converts int64 to double
	DBL d;
	INT i; //both i and f refer to the same location in memory
};

void d2d64(FILE *fout, DBL db10)
{
	//converts base 10 double to base 64 double
	union INTToDBL itd;
	itd.d = db10;
	putc((itd.i >> 56) & 0xFF, fout);
	putc((itd.i >> 48) & 0xFF, fout);
	putc((itd.i >> 40) & 0xFF, fout);
	putc((itd.i >> 32) & 0xFF, fout);
	putc((itd.i >> 24) & 0xFF, fout);
	putc((itd.i >> 16) & 0xFF, fout);
	putc((itd.i >> 8) & 0xFF, fout);
	putc(itd.i& 0xFF, fout);
}

void i2i64(FILE *fout, INS ib10)
{
	//converts base 10 int64 to base 64 int64 
	putc((ib10 >> 56) & 0xFF, fout);
	putc((ib10 >> 48) & 0xFF, fout);
	putc((ib10 >> 40) & 0xFF, fout);
	putc((ib10 >> 32) & 0xFF, fout);
	putc((ib10 >> 24) & 0xFF, fout);
	putc((ib10 >> 16) & 0xFF, fout);
	putc((ib10 >> 8) & 0xFF, fout);
	putc(ib10 & 0xFF, fout);
}

void f2f32(FILE *fout, FLT fb10)
{
	//converts base 10 float32 to base 64 float32 
	union INSToFLT itf;
	itf.f = fb10;
	putc((itf.i >> 24) & 0xFF, fout);
	putc((itf.i >> 16) & 0xFF, fout);
	putc((itf.i >> 8) & 0xFF, fout);
	putc(itf.i & 0xFF, fout);
}

void i2i32(FILE *fout, INS ib10)
{
	//base 10 int32 (short int, is) to base 64 int32
	putc((ib10 >> 24) & 0xFF, fout);
	putc((ib10 >> 16) & 0xFF, fout);
	putc((ib10 >> 8) & 0xFF, fout);
	putc(ib10 & 0xFF, fout);
}

void i2ui8(FILE *fout, UINT ib10)
{
	//base 10 int to base 2 string
	fprintf(fout, "%c", ib10 & 0xFF);
}

void writesheardisplacementb64(FILE *fout, YDJ ydj)
{
	const UCHR *bytes;
	INS header;

	header = sizeof(DBL) * ydj->njoint;
	//bytes = (UCHR *)&header;
	bytes = malloc(sizeof(INS));
	memcpy(&bytes[0], header, sizeof(INS));
	b64enc(fout, bytes, sizeof(INS));
	FREE(bytes);

	bytes = malloc(header);
	memcpy(&bytes[0], ydj->d1jsdc, header);
	b64enc(fout, bytes, header);
	FREE(bytes);
}

void writefailuremodeb64(FILE *fout, YDJ ydj)
{
	const UCHR *bytes;
	INS header;

	header = sizeof(DBL) * ydj->njoint;
	//bytes = (UCHR *)&header;
	bytes = malloc(sizeof(INS));
	memcpy(&bytes[0], header, sizeof(INS));
	b64enc(fout, bytes, sizeof(INS));
	FREE(bytes);

	//bytes = (UCHR *)&ydj->d1jfmd;
	bytes = malloc(header);
	memcpy(&bytes[0], ydj->d1jfmd, header);
	b64enc(fout, bytes, header);
	FREE(bytes);
}

void writefracturetype(FILE *fout, YDP ydp)
{
	const UCHR *bytes;
	INS header;

	header = sizeof(DBL) * ydp->nprop;
	//bytes = (UCHR *)&header;
	bytes = malloc(sizeof(INS));
	memcpy(&bytes[0], header, sizeof(INS));
	b64enc(fout, bytes, sizeof(INS));
	FREE(bytes);

	bytes = malloc(header);
	memcpy(&bytes[0], ydp->i1pefr, header);
	b64enc(fout, bytes, header);
	FREE(bytes);
}

void writecontactforceb64(FILE *fout, YDN ydn)
{
	const UCHR *bytes;
	INS header;
	DBL *cf;
	INS i, j, k;

	header = sizeof(DBL) * ydn->nnopo * ydn->nnodim;
	//bytes = (UCHR *)&header;
	bytes = malloc(sizeof(INS));
	memcpy(&bytes[0], header, sizeof(INS));
	b64enc(fout, bytes, sizeof(INS));
	FREE(bytes);

	cf = malloc(header);
	k = 0;
	for (i = 0; i != ydn->nnopo; i++)
	{
		for (j = 0; j != ydn->nnodim; j++)
		{
			cf[k] = ydn->d2nfcon[j][i];
			k++;
		}
		cf[k] = ydn->d2nfcon[j][i];
		k++;
	}

	bytes = malloc(header);
	memcpy(&bytes[0], cf, header);
	b64enc(fout, bytes, header);
	FREE(header);
	FREE(cf);
}

void writecoordsb64(FILE *fout, YDN ydn)
{

	UCHR *bytes;
	INS header;
	FLT *coords;
	INS i, j, k;

	header = sizeof(FLT) * (ydn->nnodim + 1) * ydn->nnopo;
	//bytes = (UCHR *)&header;
	bytes = malloc(sizeof(INS) * 100);
	memcpy(&bytes[0], &header, sizeof(INS));

	b64enc(fout, bytes, sizeof(INS));
	FREE(bytes);

	coords = malloc(header);
	k = 0;
	for (i = 0; i != ydn->nnopo; i++)
	{
		for (j = 0; j != ydn->nnodim; j++)
		{
			coords[k] = ydn->d2ncc[j][i];
			k++;
		}
		coords[k] = ydn->d2ncc[j][i];
		k++;
	}

	//bytes = (UCHR *)&coords;
	bytes = malloc(header);
	memcpy(&bytes[0], coords, header);
	b64enc(fout, bytes, header);
	FREE(bytes);
	FREE(coords);
}

void writeconnectivityb64(FILE *fout, YDE yde)
{
	const UCHR *bytes;
	INS header;
	INT *con;
	INS i, j, k;
	
	header = sizeof(INT) * yde->nelem * yde->nelno;
	//bytes = (UCHR *)&header;
	bytes = malloc(sizeof(INS));
	memcpy(&bytes[0], &header, sizeof(INS));
	b64enc(fout, bytes, sizeof(INS));
	FREE(bytes);

	con = malloc(header);
	k = 0;
	for (i = 0; i != yde->nelem; i++)
	{
		for (j = 0; j != yde->nelno; j++)
		{
			con[k] = yde->i2elto[j][i];
			k++;
		}
	}

	bytes = malloc(header);
	memcpy(&bytes[0], con, header);
	b64enc(fout, bytes, header);
	FREE(bytes);
	FREE(con);
}

void writeoffsetsb64(FILE *fout, YDE yde)
{
	UCHR *bytes;
	INS header;
	UINT *off;
	INS ielem;

	header = sizeof(UINT) * yde->nelem;
	//bytes = (UCHR *)&header;
	bytes = malloc(sizeof(INS));
	memcpy(&bytes[0], &header, sizeof(INS));
	b64enc(fout, bytes, sizeof(INS));
	FREE(bytes);

	off = malloc(sizeof(UINT) * yde->nelem);
	for (ielem = 0; ielem != yde->nelem; ielem++)
	{
		off[ielem] = (ielem + 1) * yde->nelno;
	}
	bytes = malloc(header);
	memcpy(&bytes[0], off, header);
	b64enc(fout, bytes, header);
	FREE(bytes);
	FREE(off);
}

void writetypesb64(FILE *fout, YDE yde)
{
	UCHR *bytes;
	INS header;
	UINT *types;
	INS ielem;

	header = sizeof(UINT) * yde->nelem;
	
	//bytes = (UCHR *)&header;
	bytes = malloc(sizeof(INS));
	memcpy(&bytes[0], &header, sizeof(INS));
	b64enc(fout, bytes, sizeof(INS));
	FREE(bytes);

	types = malloc(yde->nelem * sizeof(UINT));
	for (ielem = 0; ielem != yde->nelem; ielem++)
		types[ielem] = 5;

	//bytes = (UCHR *)&types;

	bytes = malloc(header);
	memcpy(&bytes[0], types, header);
	b64enc(fout, bytes, header);
	FREE(types);
	FREE(bytes);
}

static CHR *DBL_F[20] =	{ "%f","%1f","%2f","%3f","%4f","%5f","%6f","%7f","%8f","%9f",			   
 "%10f","%11f","%12f","%13f","%14f",															   
 "%15f","%16f" ,"%17v","%18f","%19f"															   
};																								   
																				   
static CHR *INT_LU[20] = { "%lu","%1lu","%2lu","%3lu","%4lu","%5lu","%6lu","%7lu","%8lu","%9lu",
 "%10lu","%11lu","%12lu","%13lu","%14lu",														   
 "%15lu","%16lu" ,"%17lu","%18lu","%19lu"														   
};																								   

#define FDBLw(file, x, ndigits){fprintf((file),DBL_F[ndigits],(x));}							   
#define UINTw(file, x, ndigits){fprintf((file),INT_LU[ndigits],(x));}							   
#define TAB {CHRwsp(fout); CHRwsp(fout);}														   
#define TAB2 {TAB; TAB;}																		   
#define TAB3 {TAB2; TAB; }																		   
#define TAB4 {TAB3; TAB; }																		   
#define TAB5 {TAB4; TAB; }																		   
/*  (MT)  */

static INT i1num[100];  /* numbers for space saving format     */ //What kind of space saving format? ~MT
static DBL d1num[100];  /* numbers for space saving format     */
static CHR c1code[500]; /* coded i1para in space saving format */

static void Yod2TRIELS (/* small strain softening triangle output */ nelem, fout, dcsizc, dcsizs, dcsizv, dcsizf, dpeks, dpela, dpemu, dpero, icoutp, iprop, d1nccx, d1nccy, d1ncix, d1nciy, d1nvcx, d1nvcy, i1elpr, i2elto, d2elcf) INT nelem; FILE *fout;DBL dcsizc;DBL dcsizs;DBL dcsizv;DBL dcsizf;DBL dpeks;DBL dpela;DBL dpemu;DBL dpero;INT icoutp;INT iprop;DBL *d1nccx;DBL *d1nccy;DBL *d1ncix;DBL *d1nciy;DBL *d1nvcx;DBL *d1nvcy;INT *i1elpr;INT **i2elto; DBL **d2elcf; {
	DBL voli, volc;
	DBL B[2][2];     /* left Cauchy-Green strain tensor */
	DBL D[2][2];     /* rate of deformation (stretching) tensor */
	DBL E[2][2];     /* strain tensor (small strains) */
	DBL F[2][2];     /* deformation gradient in global base */
	DBL F0[2][2];    /* initial local base */
	DBL FX[2][2];    /* current local base */
	DBL F0inv[2][2]; /* global base in initial local base */
	DBL FXinv[2][2]; /* global base in current local base */
	DBL L[2][2];     /* velocity gradient in global base */
	DBL LX[2][2];    /* vel. gradient in current local base = delta x/delta X */
	DBL T[2][2];     /* Cauchy stress */
	INT ielem;
	INT i, j, k;

	for (ielem = 0; ielem < nelem; ielem++)
	{
		if (i1elpr[ielem] == iprop)
		{ /* evaluate stress state */
			for (i = 1; i < 3; i++)
			{
				F0[0][i - 1] = d1ncix[(i2elto[i][ielem])] - d1ncix[(i2elto[0][ielem])];
				F0[1][i - 1] = d1nciy[(i2elto[i][ielem])] - d1nciy[(i2elto[0][ielem])];
				FX[0][i - 1] = d1nccx[(i2elto[i][ielem])] - d1nccx[(i2elto[0][ielem])];
				FX[1][i - 1] = d1nccy[(i2elto[i][ielem])] - d1nccy[(i2elto[0][ielem])];
				LX[0][i - 1] = d1nvcx[(i2elto[i][ielem])] - d1nvcx[(i2elto[0][ielem])];
				LX[1][i - 1] = d1nvcy[(i2elto[i][ielem])] - d1nvcy[(i2elto[0][ielem])];
			}
			YMATINV2(F0, F0inv, voli);
			YMATINV2(FX, FXinv, volc);
			for (i = 0; i < 2; i++)
			{
				for (j = 0; j < 2; j++)
				{
					F[i][j] = R0;
					L[i][j] = R0;
					for (k = 0; k < 2; k++)
					{
						F[i][j] = F[i][j] + FX[i][k] * F0inv[k][j];
						L[i][j] = L[i][j] + LX[i][k] * FXinv[k][j];
					}
				}
			}
			for (i = 0; i < 2; i++)
			{
				for (j = 0; j < 2; j++)
				{
					B[i][j] = R0;
					for (k = 0; k < 2; k++)
					{
						B[i][j] = B[i][j] + F[i][k] * F[j][k]; /* left Cauchy-Green strain */
					}
					D[i][j] = RP5 * (L[i][j] + L[j][i]); /* rate of deformation      */
					if (i == j)
					{
						E[i][j] = RP5 * (B[i][j] - R1); /* small strain             */
					}
					else
					{
						E[i][j] = RP5 * B[i][j];
					}
				}
			}
			for (i = 0; i < 2; i++) /* Cauchy stress */
			{
				for (j = 0; j < 2; j++)
				{
					T[i][j] = (R2 * dpemu * E[i][j]) * (voli / volc) + dpeks * D[i][j];
					if (i == j)
						T[i][j] = T[i][j] + dpela * (volc / voli - voli / volc);
				}
			}

			/* prepare output */
			i1num[0] = icoutp;										//number of character for each number ~MT
			i1num[1] = 17;											//size of this array
			d1num[3] = d1nccx[i2elto[0][ielem]] / dcsizc; 			//nodal coordinates current x of 0th node in ielem th element ~MT
			d1num[4] = d1nccx[i2elto[1][ielem]] / dcsizc;			//nodal coordinates current x of 1st node in ielem th element ~MT
			d1num[5] = d1nccx[i2elto[2][ielem]] / dcsizc;			//nodal coordinates current x of 2nd node in ielem th element ~MT
			d1num[6] = d1nccy[i2elto[0][ielem]] / dcsizc;			//nodal coordinates current y of 0th node in ielem th element ~MT
			d1num[7] = d1nccy[i2elto[1][ielem]] / dcsizc;			//nodal coordinates current y of 1st node in ielem th element ~MT
			d1num[8] = d1nccy[i2elto[2][ielem]] / dcsizc;			//nodal coordinates current y of 2nd node in ielem th element ~MT
			d1num[9] = (d1nvcx[i2elto[0][ielem]] +					//element velocities current x
				d1nvcx[i2elto[1][ielem]] +
				d1nvcx[i2elto[2][ielem]]) /
				(R3 * dcsizv);
			d1num[10] = (d1nvcy[i2elto[0][ielem]] +					//element velocities current y ~MT
				d1nvcy[i2elto[1][ielem]] +
				d1nvcy[i2elto[2][ielem]]) /
				(R3 * dcsizv);
			d1num[11] = T[0][0] / dcsizs;							//Cauchy ~MT
			d1num[12] = T[1][1] / dcsizs;
			d1num[13] = T[0][1] / dcsizs;
			d1num[14] = R0; /* elastic damage */
			d1num[15] = d2elcf[ielem][0] / dcsizf;					//contact force x ~MT
			d1num[16] = d2elcf[ielem][1] / dcsizf;					//contact force y ~MT
			for (i = 3; i < 17; i++)
			{
				d1num[i] = MAXIM((-R1), MINIM(d1num[i], R1)); //larger than -1?
			}

			//What do you mean: translate into INT? By which encoder?

			/* translate into INT */
			codeDBLtoINT(d1num, i1num);
			i1num[2] = YTE2TRIELS;
			codeINTtoCHR(c1code, i1num);
			CHRw(fout, c1code);
			CHRwcr(fout);
		}
	}
}

static void Yod2TRISOF(/* small strain softening triangle output */ nelem,fout,dcsizc, dcsizs, dcsizv,dpeks, dpela, dpemu, dpero,icoutp, iprop,d1nccx, d1nccy, d1ncix, d1nciy, d1nvcx,d1nvcy, d1sdel, i1elpr, i2elto) INT nelem;FILE *fout;DBL dcsizc;DBL dcsizs;DBL dcsizv;DBL dpeks;DBL dpela;DBL dpemu;DBL dpero;INT icoutp;INT iprop;DBL *d1nccx;DBL *d1nccy;DBL *d1ncix;DBL *d1nciy;DBL *d1nvcx;DBL *d1nvcy;DBL *d1sdel;INT *i1elpr;INT **i2elto; {
	DBL voli, volc;
	DBL B[2][2];     /* left Cauchy-Green strain tensor */
	DBL D[2][2];     /* rate of deformation (stretching) tensor */
	DBL E[2][2];     /* strain tensor (small strains) */
	DBL F[2][2];     /* deformation gradient in global base */
	DBL F0[2][2];    /* initial local base */
	DBL FX[2][2];    /* current local base */
	DBL F0inv[2][2]; /* global base in initial local base */
	DBL FXinv[2][2]; /* global base in current local base */
	DBL L[2][2];     /* velocity gradient in global base */
	DBL LX[2][2];    /* vel. gradient in current local base = delta x/delta X */
	DBL T[2][2];     /* Cauchy stress */
	INT ielem;
	INT i, j, k;

	for (ielem = 0; ielem < nelem; ielem++)
	{
		if (i1elpr[ielem] == iprop)
		{ /* evaluate stress state */
			for (i = 1; i < 3; i++)
			{
				F0[0][i - 1] = d1ncix[(i2elto[i][ielem])] - d1ncix[(i2elto[0][ielem])];
				F0[1][i - 1] = d1nciy[(i2elto[i][ielem])] - d1nciy[(i2elto[0][ielem])];
				FX[0][i - 1] = d1nccx[(i2elto[i][ielem])] - d1nccx[(i2elto[0][ielem])];
				FX[1][i - 1] = d1nccy[(i2elto[i][ielem])] - d1nccy[(i2elto[0][ielem])];
				LX[0][i - 1] = d1nvcx[(i2elto[i][ielem])] - d1nvcx[(i2elto[0][ielem])];
				LX[1][i - 1] = d1nvcy[(i2elto[i][ielem])] - d1nvcy[(i2elto[0][ielem])];
			}
			YMATINV2(F0, F0inv, voli);
			YMATINV2(FX, FXinv, volc);
			for (i = 0; i < 2; i++)
			{
				for (j = 0; j < 2; j++)
				{
					F[i][j] = R0;
					L[i][j] = R0;
					for (k = 0; k < 2; k++)
					{
						F[i][j] = F[i][j] + FX[i][k] * F0inv[k][j];
						L[i][j] = L[i][j] + LX[i][k] * FXinv[k][j];
					}
				}
			}
			for (i = 0; i < 2; i++)
			{
				for (j = 0; j < 2; j++)
				{
					B[i][j] = R0;
					for (k = 0; k < 2; k++)
					{
						B[i][j] = B[i][j] + F[i][k] * F[j][k]; /* left Cauchy-Green strain */
					}
					D[i][j] = RP5 * (L[i][j] + L[j][i]); /* rate of deformation      */
					if (i == j)
					{
						E[i][j] = RP5 * (B[i][j] - R1); /* small strain             */
					}
					else
					{
						E[i][j] = RP5 * B[i][j];
					}
				}
			}
			for (i = 0; i < 2; i++) /* Cauchy stress */
			{
				for (j = 0; j < 2; j++)
				{
					T[i][j] = (R1 - d1sdel[ielem]) *
						(R2 * dpemu * E[i][j]) * (voli / volc) +
						dpeks * D[i][j];
					if (i == j)
						T[i][j] = T[i][j] + dpela * (volc / voli - voli / volc);
				}
			}
			/* prepare output */
			i1num[0] = icoutp;
			i1num[1] = 15;
			d1num[3] = d1nccx[i2elto[0][ielem]] / dcsizc;
			d1num[4] = d1nccx[i2elto[1][ielem]] / dcsizc;
			d1num[5] = d1nccx[i2elto[2][ielem]] / dcsizc;
			d1num[6] = d1nccy[i2elto[0][ielem]] / dcsizc;
			d1num[7] = d1nccy[i2elto[1][ielem]] / dcsizc;
			d1num[8] = d1nccy[i2elto[2][ielem]] / dcsizc;
			d1num[9] = (d1nvcx[i2elto[0][ielem]] +
				d1nvcx[i2elto[1][ielem]] +
				d1nvcx[i2elto[2][ielem]]) /
				(R3 * dcsizv);
			d1num[10] = (d1nvcy[i2elto[0][ielem]] +
				d1nvcy[i2elto[1][ielem]] +
				d1nvcy[i2elto[2][ielem]]) /
				(R3 * dcsizv);
			d1num[11] = T[0][0] / dcsizs;
			d1num[12] = T[1][1] / dcsizs;
			d1num[13] = T[0][1] / dcsizs;
			d1num[14] = d1sdel[ielem];
			for (i = 3; i < 15; i++)
			{
				d1num[i] = MAXIM((-R1), MINIM(d1num[i], R1));
			}
			/* translate into INT */
			codeDBLtoINT(d1num, i1num);
			i1num[2] = YTE2TRIELS;
			codeINTtoCHR(c1code, i1num);
			CHRw(fout, c1code);
			CHRwcr(fout);
		}
	}
}

static void Yod2TRIRIG(/* small strain elastic triangle output */ nelem,fout,dcsizc, dcsizv,icoutp, iprop,d1nccx, d1nccy, d1nvcx, d1nvcy, i1elpr,i2elto) INT nelem;FILE *fout;DBL dcsizc;DBL dcsizv;INT icoutp;INT iprop;DBL *d1nccx;DBL *d1nccy;DBL *d1nvcx;DBL *d1nvcy;INT *i1elpr;INT **i2elto; {
	INT ielem;INT i;
	for (ielem = 0; ielem < nelem; ielem++)
	{
		if (i1elpr[ielem] == iprop)
		{ /* prepare output */
			i1num[0] = icoutp;
			i1num[1] = 15;
			d1num[3] = d1nccx[i2elto[0][ielem]] / dcsizc;
			d1num[4] = d1nccx[i2elto[1][ielem]] / dcsizc;
			d1num[5] = d1nccx[i2elto[2][ielem]] / dcsizc;
			d1num[6] = d1nccy[i2elto[0][ielem]] / dcsizc;
			d1num[7] = d1nccy[i2elto[1][ielem]] / dcsizc;
			d1num[8] = d1nccy[i2elto[2][ielem]] / dcsizc;
			d1num[9] = (d1nvcx[i2elto[0][ielem]] +
				d1nvcx[i2elto[1][ielem]] +
				d1nvcx[i2elto[2][ielem]]) /
				(R3 * dcsizv);
			d1num[10] = (d1nvcy[i2elto[0][ielem]] +
				d1nvcy[i2elto[1][ielem]] +
				d1nvcy[i2elto[2][ielem]]) /
				(R3 * dcsizv);
			d1num[11] = R0;
			d1num[12] = R0;
			d1num[13] = R0;
			d1num[14] = R0;
			for (i = 3; i < 15; i++)
			{
				d1num[i] = MAXIM((-R1), MINIM(d1num[i], R1)); //at least -1? Why?
			}
			/* translate into INT */
			codeDBLtoINT(d1num, i1num);
			i1num[2] = YTE2TRIELS;
			codeINTtoCHR(c1code, i1num);
			CHRw(fout, c1code);
			CHRwcr(fout);
		}
	}
}

static void Yod2JOINTSBROKEN(/* 2D joint output */ nelem,fout,dcsizc, dcsizv, dcsizf, dcsizd, dcsiza,icoutp, iprop,d1nccx, d1nccy, d1nvcx, d1nvcy, d1sdel,i1elpr, i2elto, d2tcs, dcsizs, i2eljp,i1elty, d1elsf, d2elcf,njoint, i1jtid, d1jknc, d1jksc, d1jnst,d1jsst, d1japc, d1japh, d1jsdc, d1jdlc,d1jphi, d1jfmd) INT nelem;FILE *fout;DBL dcsizc;DBL dcsizv;DBL dcsizf;DBL dcsizd;DBL dcsiza;INT icoutp;INT iprop;DBL *d1nccx;DBL *d1nccy;DBL *d1nvcx;DBL *d1nvcy;DBL *d1sdel;INT *i1elpr;INT **i2elto;DBL **d2tcs;DBL dcsizs;INT **i2eljp;INT *i1elty;DBL *d1elsf;DBL **d2elcf;INT njoint;INT *i1jtid;DBL *d1jknc;DBL *d1jksc;DBL *d1jnst;DBL *d1jsst;DBL *d1japc;DBL *d1japh;DBL *d1jsdc;DBL *d1jdlc;DBL *d1jphi;DBL *d1jfmd; {
	INT i, ielem, ijoint;
	for (ielem = 0; ielem < nelem; ielem++)
	{
		if (i1elty[ielem] == 1) /* output boundary */
		{
			i1num[0] = icoutp;
			i1num[1] = 19;
			d1num[3] = d1nccx[i2elto[0][ielem]] / dcsizc;
			d1num[4] = d1nccx[i2elto[1][ielem]] / dcsizc;
			d1num[5] = d1nccx[i2elto[2][ielem]] / dcsizc;
			d1num[6] = d1nccx[i2elto[3][ielem]] / dcsizc;
			d1num[7] = d1nccy[i2elto[0][ielem]] / dcsizc;
			d1num[8] = d1nccy[i2elto[1][ielem]] / dcsizc;
			d1num[9] = d1nccy[i2elto[2][ielem]] / dcsizc;
			d1num[10] = d1nccy[i2elto[3][ielem]] / dcsizc;
			d1num[11] = d2elcf[i2eljp[0][ielem]][0] / dcsizf;
			d1num[12] = d2elcf[i2eljp[0][ielem]][1] / dcsizf;
			d1num[13] = R0;
			d1num[14] = R0;
			d1num[15] = R0;
			d1num[16] = R0;
			d1num[17] = R0;
			d1num[18] = R0;
			for (i = 3; i < i1num[1]; i++)
			{
				d1num[i] = MAXIM((-R1), MINIM(d1num[i], R1));
			}
			/* translate into INT */
			codeDBLtoINT(d1num, i1num);
			i1num[2] = YTE2BOUNDS;
			codeINTtoCHR(c1code, i1num);
			CHRw(fout, c1code);
			CHRwcr(fout);
		}
		else if (i1elty[ielem] > 1) /* output fracture joints */
		{
			ijoint = ielem - (nelem - njoint);
			i1num[0] = icoutp;
			i1num[1] = 19;
			d1num[3] = d1nccx[i2elto[0][ielem]] / dcsizc;
			d1num[4] = d1nccx[i2elto[1][ielem]] / dcsizc;
			d1num[5] = d1nccx[i2elto[2][ielem]] / dcsizc;
			d1num[6] = d1nccx[i2elto[3][ielem]] / dcsizc;
			d1num[7] = d1nccy[i2elto[0][ielem]] / dcsizc;
			d1num[8] = d1nccy[i2elto[1][ielem]] / dcsizc;
			d1num[9] = d1nccy[i2elto[2][ielem]] / dcsizc;
			d1num[10] = d1nccy[i2elto[3][ielem]] / dcsizc;
			d1num[11] = d1jfmd[ijoint] / R2;
			d1num[12] = d1jnst[ijoint] / dcsizs;
			d1num[13] = d1jsst[ijoint] / dcsizs;
			d1num[14] = (DABS(d1elsf[i2eljp[0][ielem]]) + DABS(d1elsf[i2eljp[1][ielem]])) / R2 / dcsizf;
			d1num[15] = log10(d1japc[ijoint]) / log10(dcsiza);
			d1num[16] = DABS(d1jsdc[ijoint]) / dcsizd;
			d1num[17] = d1jdlc[ijoint] / dcsizd;
			d1num[18] = log10(d1japh[ijoint]) / log10(dcsiza);
			for (i = 3; i < i1num[1]; i++)
			{
				d1num[i] = MAXIM((-R1), MINIM(d1num[i], R1));
			}
			/* translate into INT */
			codeDBLtoINT(d1num, i1num);
			i1num[2] = YTE2JOINTS;
			codeINTtoCHR(c1code, i1num);
			CHRw(fout, c1code);
			CHRwcr(fout);
		}
	}
}

static void Yod2AESOURCE(/* acoustic emission source */ nelem,fout,dcsizc,icoutp,i2elto, d1elme, i1elyi,d1nccx, d1nccy) INT nelem;FILE *fout;DBL dcsizc;INT icoutp;INT **i2elto;DBL *d1elme;INT *i1elyi;DBL *d1nccx;DBL *d1nccy; {
	INT i, ielem;
	DBL dcsizm = 15.0;
	for (ielem = 0; ielem < nelem; ielem++)
	{
		if (i1elyi[ielem] >= 2)
		{
			i1num[0] = icoutp;
			i1num[1] = 10;
			d1num[3] = (d1nccx[i2elto[0][ielem]] + d1nccx[i2elto[1][ielem]] +
				d1nccx[i2elto[2][ielem]] + d1nccx[i2elto[3][ielem]]) /
				R4 / dcsizc;
			d1num[4] = (d1nccy[i2elto[0][ielem]] + d1nccy[i2elto[1][ielem]] +
				d1nccy[i2elto[2][ielem]] + d1nccy[i2elto[3][ielem]]) /
				R4 / dcsizc;
			d1num[5] = d1elme[ielem] / dcsizm;
			d1num[6] = R0;
			d1num[7] = R0;
			d1num[8] = R0;
			d1num[9] = R0;
			for (i = 3; i < 10; i++)
			{
				d1num[i] = MAXIM((-R1), MINIM(d1num[i], R1));
			}
			/* translate into INT */
			codeDBLtoINT(d1num, i1num);
			i1num[2] = YTE2AESOUR;
			codeINTtoCHR(c1code, i1num);
			CHRw(fout, c1code);
			CHRwcr(fout);
		}
	}
}

static void Yod2JOINTSINTACT(/* 2D joint output */ nelem,fout,dcsizc, dcsizv,dpefs, dpeft, dpegf, dpeks, dpepe,icoutp, iprop,d1nccx, d1nccy, d1nvcx, d1nvcy, d1sdel,i1elpr, i2elto) INT nelem;FILE *fout;DBL dcsizc;DBL dcsizv;DBL dpefs;DBL dpeft;DBL dpegf;DBL dpeks;DBL dpepe;INT icoutp;INT iprop;DBL *d1nccx;DBL *d1nccy;DBL *d1nvcx;DBL *d1nvcy;DBL *d1sdel;INT *i1elpr;INT **i2elto; {
	DBL small, o1, o2, s1, s2, op, sp, ot, st;
	DBL e1x, e1y, h;
	INT ielem, i, i0, i1, i2, i3;

	small = EPSILON;
	for (ielem = 0; ielem < nelem; ielem++)
	{
		if (i1elpr[ielem] == iprop)
		{
			i0 = i2elto[0][ielem];
			i1 = i2elto[1][ielem];
			i2 = i2elto[2][ielem];
			i3 = i2elto[3][ielem];
			e1x = RP5 * (d1nccx[i1] + d1nccx[i2] - d1nccx[i0] - d1nccx[i3]);
			e1y = RP5 * (d1nccy[i1] + d1nccy[i2] - d1nccy[i0] - d1nccy[i3]);
			h = SQRT(e1x * e1x + e1y * e1y);
			e1x = e1x / (h + small);
			e1y = e1y / (h + small);
			s1 = (d1nccy[i0] - d1nccy[i3]) * e1y + (d1nccx[i0] - d1nccx[i3]) * e1x;
			s2 = (d1nccy[i1] - d1nccy[i2]) * e1y + (d1nccx[i1] - d1nccx[i2]) * e1x;
			o1 = (d1nccy[i0] - d1nccy[i3]) * e1x - (d1nccx[i0] - d1nccx[i3]) * e1y;
			o2 = (d1nccy[i1] - d1nccy[i2]) * e1x - (d1nccx[i1] - d1nccx[i2]) * e1y;
			op = R2 * h * dpeft / dpepe;
			sp = R2 * h * dpefs / dpepe;
			ot = MAXIM((R2 * op), (R3 * dpegf / dpeft));
			st = MAXIM((R2 * sp), (R3 * dpegf / dpefs));

			/* prepare output */
			if ((((o1 + o2) / (R2 * ot)) > RP1) || (((s1 + s2) / (R2 * st)) > RP1))
			{
				i1num[0] = icoutp;
				i1num[1] = 14;
				d1num[3] = d1nccx[i2elto[0][ielem]] / dcsizc;
				d1num[4] = d1nccx[i2elto[1][ielem]] / dcsizc;
				d1num[5] = d1nccx[i2elto[2][ielem]] / dcsizc;
				d1num[6] = d1nccx[i2elto[3][ielem]] / dcsizc;
				d1num[7] = d1nccy[i2elto[0][ielem]] / dcsizc;
				d1num[8] = d1nccy[i2elto[1][ielem]] / dcsizc;
				d1num[9] = d1nccy[i2elto[2][ielem]] / dcsizc;
				d1num[10] = d1nccy[i2elto[3][ielem]] / dcsizc;
				d1num[11] = R1;
				d1num[12] = o1 / ot;
				d1num[13] = o2 / st;
				for (i = 3; i < 14; i++)
				{
					d1num[i] = MAXIM((-R1), MINIM(d1num[i], R1));
				}

				/* translate into INT */
				codeDBLtoINT(d1num, i1num);
				i1num[2] = YTE2JOINTS;
				codeINTtoCHR(c1code, i1num);
				CHRw(fout, c1code);
				CHRwcr(fout);
			}
		}
	}
}

static void Yod2CONTACTFORCE (/* output contact forces */ fcnf,nnopo, d1nfconx, d1nfcony, i1nopr) FILE *fcnf; INT nnopo;DBL *d1nfconx;DBL *d1nfcony;INT *i1nopr; 
{
	INT inode;
	DBL dsumfx, dsumfy, tmp, tmpy;

	dsumfx = R0;
	dsumfy = R0;
	tmp = 0.0;
	tmpy = 0.0;
	for (inode = 0; inode < nnopo; inode++)
	{
		dsumfx = dsumfx + DABS(d1nfconx[inode]);
		dsumfy = dsumfy + DABS(d1nfcony[inode]);

		if (i1nopr[inode] == 5)
		{
			tmp = tmp + d1nfconx[inode];
			tmpy = tmpy + d1nfcony[inode];
		}
	}

	DBLw(fcnf, dsumfx, 10);
	CHRwsp(fcnf);
	DBLw(fcnf, dsumfy, 10);
	CHRwsp(fcnf);
	DBLw(fcnf, tmp, 10);
	CHRwsp(fcnf);
	DBLw(fcnf, tmpy, 10);
	CHRwcr(fcnf);
}

/*  MT  */
DBL d1min(DBL *a, INT len)
{
	//finds min element of 1D double array
	DBL f = a[0]; //min element
	INT i;
	for (i = 1; i < len; ++i)
		if (f > a[i])
			f = a[i];
	return f;
}

INT i1min(INT *a, INT len)
{
	//finds min element of 1D int array
	INT f = a[0]; //min element
	INT i;
	for (i = 1; i < len; ++i)
		if (f > a[i])
			f = a[i];
	return f;
}

DBL d1max(DBL *a, INT len)
{
	//finds max element of 1D double array
	DBL f = a[0]; //max element
	INT i;
	for (i = 1; i < len; ++i)
		if (f < a[i])
			f = a[i];
	return f;
}

INT i1max(INT *a, INT len)
{
	//finds max element of 1D int array
	INT f = a[0];
	INT i;
	for (i = 1; i < len; ++i)
		if (f < a[i])
			f = a[i];
	return f;
}

DBL d2min(DBL **a, INT len1, INT len2)
{
	//finds min element of 2D double array
	DBL m = a[0][0]; //min element
	INT i;//loop indices
	INT j;

	for (j= 0; j < len1; j++)
		for (i = 0; i < len2; i++)
			if (m > a[j][i])
				m = a[j][i];
	return m;
}

INT i2min(INT **a, INT len1, INT len2)
{
	//finds min element of 2D int array
	INT m = a[0][0]; //max element
	INT i;
	INT j;
	for (j = 0; j < len1; j++)
		for (i = 0; i < len2; i++)
			if (m > a[j][i])
				m = a[j][i];
	return m;
}

DBL d2max(DBL **a, INT len1, INT len2)
{
	//finds max element of 2D double array
	DBL m = a[0][0];
	INT i;
	INT j;
	for (j = 0; j < len1; j++)
		for (i = 0; i < len2; i++)
			if (m < a[j][i])
				m = a[j][i];
	return m;
}

INT i2max(const INT **a, INT len1, INT len2)
{
	//finds max element of 2D int array
	INT m = a[0][0];
	INT i;
	INT j;
	for (j = 0; j < len1; j++)
		for (i = 0; i < len2; i++)
			if (m < a[j][i])
				m = a[j][i];
	return m;
}

/*  (MT)  */

void Yod(namep, yd) CHR *namep; YD yd; {
	YDC ydc = &(yd->ydc);
	YDE yde = &(yd->yde);
	YDI ydi = &(yd->ydi);
	YDN ydn = &(yd->ydn);
	YDO ydo = &(yd->ydo);
	YDP ydp = &(yd->ydp);
	YDJ ydj = &(yd->ydj);
	INT iprop, ihys, i;
	CHR namef[300];
	CHR cindex[50];
	static INT ncall = 0;
	FILE *fout = FILENULL;
	FILE *fcnf = FILENULL;
	DBL tmp, tmpy;

	/*  MT  */

	INT ielem, j; //loop index

	//file output flags
	INT ym; //ym file output (original output)
	INT vtu_leg; //vtu legacy file output
	INT vtu_xml; //vtu xml file ouput

	//data output flags
	INT app; //appended data output (only vtk xml)
	INT b10; //ascii data output
	INT b2; //binary data output (only vtk legacy files)
	INT b64; //b64 output (appended if app==1, inline if app==0)
	INT raw; //appended raw output (only vtk xml and if app==1)
	INT off; //offset of appended data (only vtk xml)
	ym = 0;
	vtu_leg = 0;
	vtu_xml = 1;

	b10 = 0;
	app = 0;
	b64 = 1; 


	/*  (MT)  */

	ncall = ncall + 1;

	/* output hystory */

	for (ihys = 0; ihys < (ydo->nohys); ihys++)
	{
		if ((ydo->d1ohyt[ihys]) == (ydc->dctime))
		{
			if ((ydo->f2ohyf[ihys]) == FILENULL)
			{
				CHRcpy(namef, namep);
				SINTw(cindex, ihys, 0);
				CHRcat(namef, "h");
				CHRcat(namef, cindex);
				ydo->f2ohyf[ihys] = fopen(namef, "a");
			}
			if ((ydo->f2ohyf[ihys]) != FILENULL)
			{
				tmp = (ydo->d1ohyt[ihys]) * (ydo->d1ohyc[ihys]);
				DBLw((ydo->f2ohyf[ihys]), tmp, 10);
				CHRwsp(ydo->f2ohyf[ihys]);
				tmp = (ydo->d1ohys[ihys]) * (ydo->d1ohyf[ihys]);
				DBLw((ydo->f2ohyf[ihys]), tmp, 10);
				CHRwcr(ydo->f2ohyf[ihys]);

				tmp = R0;
				tmpy = R0;
				for (i = 0; i < ydn->nnopo; i++)
				{
					if (ydn->i1nopr[i] == 1)
					{
						tmp = tmp + ydn->d2nfcon[0][i];
						tmpy = tmpy + ydn->d2nfcon[1][i];
					}
				}
				CHRwsp(ydo->f2ohyf[ihys]);
				DBLw((ydo->f2ohyf[ihys]), tmp, 10);
				CHRwsp(ydo->f2ohyf[ihys]);
				DBLw((ydo->f2ohyf[ihys]), tmpy, 10);
				CHRwcr(ydo->f2ohyf[ihys]);
			}
		}
	}

	if ((ncall > 100) || ((ydc->ncstep) >= (ydc->mcstep - 2)))
	{
		ncall = 0;
		for (ihys = 0; ihys < (ydo->nohys); ihys++)
		{
			if ((ydo->f2ohyf[ihys]) != FILENULL)
			{
				fclose(ydo->f2ohyf[ihys]);
			}
			ydo->f2ohyf[ihys] = FILENULL;
		}
	}

	/* output anymation */

	fout = FILENULL;
	fcnf = FILENULL;

	if ((ydc->ncstep % ydc->icoutf) == 0)
	{
		CHRcpynoext(namef, namep);
		SINTw(cindex, ydc->icouti, 0);
		CHRcat(namef, cindex);

		/*  MT  */
		if (ym)		
		/*  (MT)  */
		{			
			CHRcat(namef, ".ym");
			fout = fopen(namef, "w");
			if (fout != FILENULL)
			{
				CHRw(fout, "CODED");
				CHRwcr(fout);
				for (iprop = 0; iprop < ydp->nprop; iprop++)
				{
					if ((ydp->i1ptyp[iprop]) == (YTE2TRIELS)) //check if property type == plain stress triangle (=1) ~MT
					{
						Yod2TRIELS // small strain softening triangle output ~MT
						(
							yde->nelem,
							fout,
							ydc->dcsizc, ydc->dcsizs, ydc->dcsizv, ydc->dcsizf,
							ydp->d1peks[iprop], ydp->d1pela[iprop], ydp->d1pemu[iprop],
							ydp->d1pero[iprop],
							ydc->icoutp, iprop,
							ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2nci[0], ydn->d2nci[1], ydn->d2nvc[0],
							ydn->d2nvc[1],
							yde->i1elpr, yde->i2elto, yde->d2elcf
						);
					}
					else if ((ydp->i1ptyp[iprop]) == (YTE2TRISOF))
					{
						Yod2TRISOF(/* small strain elastic triangle output */
							yde->nelem,
							fout,
							ydc->dcsizc, ydc->dcsizs, ydc->dcsizv,
							ydp->d1peks[iprop], ydp->d1pela[iprop], ydp->d1pemu[iprop],
							ydp->d1pero[iprop],
							ydc->icoutp, iprop,
							ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2nci[0], ydn->d2nci[1], ydn->d2nvc[0],
							ydn->d2nvc[1], yde->d2elst[ydp->i1psde[iprop]], yde->i1elpr,
							yde->i2elto);
					}
					else if ((ydp->i1ptyp[iprop]) == (YTE2TRIRIG))
					{
						Yod2TRIRIG(/* rigid triangle output */
							yde->nelem,
							fout,
							ydc->dcsizc, ydc->dcsizv,
							ydc->icoutp, iprop,
							ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2nvc[0], ydn->d2nvc[1], yde->i1elpr,
							yde->i2elto);
					}
					else if ((ydp->i1ptyp[iprop]) == (YTE2JOINTS))
					{

						//Yod2JOINTSINTACT(  /* 2D joint output */
						/*yde->nelem ,
						fout       ,
						ydc->dcsizc       ,ydc->dcsizv       ,
						ydp->d1pefs[iprop],ydp->d1peft[iprop],ydp->d1pegf[iprop],
						ydp->d1peks[iprop],ydp->d1pepe[iprop],
						ydc->icoutp       ,iprop             ,
						ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nvc[0],ydn->d2nvc[1],
						yde->d2elst[ydp->i1psde[iprop]],
						yde->i1elpr,yde->i2elto
						);*/

						Yod2JOINTSBROKEN(/* broken joints */
							yde->nelem,
							fout,
							ydc->dcsizc, ydc->dcsizv, ydc->dcsizf, ydc->dcsizd, ydc->dcsiza,
							ydc->icoutp, iprop,
							ydn->d2ncc[0], ydn->d2ncc[1], ydn->d2nvc[0], ydn->d2nvc[1],
							yde->d2elst[ydp->i1psde[iprop]],
							yde->i1elpr, yde->i2elto, yde->d2tcs, ydc->dcsizs, yde->i2eljp,
							yde->i1elty, yde->d1elsf, yde->d2elcf,
							ydj->njoint, ydj->i1jtid, ydj->d1jknc, ydj->d1jksc, ydj->d1jnst,
							ydj->d1jsst, ydj->d1japc, ydj->d1japh, ydj->d1jsdc, ydj->d1jdlc,
							ydj->d1jphi, ydj->d1jfmd);

						Yod2AESOURCE(/* acoustic emission source */
							yde->nelem,
							fout,
							ydc->dcsizc,
							ydc->icoutp,
							yde->i2elto, yde->d1elme, yde->i1elyi,
							ydn->d2ncc[0], ydn->d2ncc[1]);
					}
				}
			}
		}
		/*  MT  */
		else if (vtu_leg)
		{
			CHRcat(namef, ".vtu");
			if (b10)
			{
				fout = fopen(namef, "w");
			}
			else
			{
				fout = fopen(namef, "wb"); //open in binary mode
			}
			if (fout != FILENULL)
			{
				CHRw(fout, "# vtk DataFile Version 4.2\n"); // older version to avoid paraview warnings
				CHRw(fout, namef); 
				CHRwcr(fout);
				if (b10)
				{
					CHRw(fout, "ASCII");
				}
				else
				{
					CHRw(fout, "BINARY");
				}
				CHRwcr(fout);
				CHRw(fout, "DATASET UNSTRUCTURED_GRID"); CHRwcr(fout);
				CHRw(fout, "POINTS"); CHRwsp(fout);
				fprintf(fout, "%d", ydn->nnopo); CHRwsp(fout);
				CHRw(fout, "float"); CHRwcr(fout);
				if (b10)
				{
					for (i = 0; i != ydn->nnopo; i++)
					{
						for (j = 0; j != ydn->nnodim; j++)
						{
							fprintf(fout, "%f", ydn->d2ncc[j][i]); CHRwsp(fout);
						}
						fprintf(fout, "0 ");
					}
				}
				else
				{
					for (i = 0; i != ydn->nnopo; i++)
					{
						for (j = 0; j != ydn->nnodim; j++)
						{
							d2d64(fout, ydn->d2ncc[j][i]);
						}
						d2d64(fout, 0);
					}
				}
				CHRwcr(fout);
				CHRw(fout, "CELLS"); CHRwsp(fout);
				UINTw(fout, yde->nelem, 0); CHRwsp(fout);
				UINTw(fout, yde->nelem * (yde->nelno + 1), 0); CHRwcr(fout);

				if (b10)
				{
					for (ielem = 0; ielem != yde->nelem; ielem++)
					{
						UINTw(fout, yde->nelno, 0); CHRwsp(fout);
						for (i = 0; i != yde->nelno; i++)
						{
							UINTw(fout, yde->i2elto[i][ielem], 0); CHRwsp(fout);
						}
						fprintf(fout, "\n");
					}
				}
				else
				{
					for (ielem = 0; ielem != yde->nelem; ielem++)
					{
						i2i32(fout, yde->nelno);
						for (i = 0; i != yde->nelno; i++)
						{
							i2i64(fout, yde->i2elto[i][ielem]);
						}
						fprintf(fout, "\n");
					}
				}

				CHRw(fout, "CELL_TYPES"); CHRwsp(fout);
				UINTw(fout, yde->nelem, 0); CHRwcr(fout);
				if (b10)
				{
					for (ielem = 0; ielem < yde->nelem; ielem++)
					{
						UINTw(fout, 5, 1); CHRwcr(fout);
					}
				}
				else
				{
					for (ielem = 0; ielem < yde->nelem; ielem++)
					{
						i2ui8(fout, 5);
					}
					CHRwcr(fout);
				}
			}
		}
		else if (vtu_xml)
		{
			//INT i = 17362929;
			//short *ptr = &i;
			//for(INT n = 0; n < 4; n++)
			//{
			//	printf("ptr: %d", *ptr);
			//	ptr++;
			//}

			CHRcat(namef, ".vtu");
			fout = fopen(namef, "w");

			if (fout != FILENULL)
			{
				//data offset
				off = 0;

				//header
				CHRw(fout, "<?xml version=\"1.0\"?>"); CHRwcr(fout);
				//CHRw(fout, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order = \"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"); 
				CHRw(fout, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order = \"LittleEndian\">");
				CHRwcr(fout);
				TAB; CHRw(fout, "<UnstructuredGrid>"); CHRwcr(fout);
				
				//piece
				TAB2; CHRw(fout, "<Piece NumberOfPoints=\"");
				UINTw(fout, ydn->nnopo, 0);
				CHRw(fout, "\" NumberOfCells=\"");
				UINTw(fout, yde->nelem, 0);
				CHRw(fout, "\">");
				CHRwcr(fout);

//				TAB3; CHRw(fout, "<PointData Scalars=\"Fracture type\" Vectors=\"Contact force\">"); CHRwcr(fout);
//				TAB4; CHRw(fout, "<DataArray type=\"Float64\" Name=\"Shear displacement\" format=\"");
//				if (app)
//				{
//					CHRw(fout, "appended");
//				}
//				else if (b10)
//				{
//					CHRw(fout, "ascii");
//				}
//				else
//				{
//					CHRw(fout, "binary");
//				}
//				//CHRw(fout, "\" RangeMin=\"");
//				//fprintf(fout, "%f", d1min(ydj->d1jsdc, ydj->njoint));
//				//CHRw(fout, "\" RangeMax=\"");
//				//fprintf(fout, "%f", d1max(ydj->d1jsdc, ydj->njoint));
//
//				if (app)
//				{
//					CHRw(fout, "offset=\"");
//					fprintf(fout, off);
//					CHRw(fout, "\"/>");
//					off += sizeof(INS);
//					off += sizeof(DBL) * ydj->njoint;
//				}
//				else
//				{
//					CHRw(fout, "\">"); CHRwcr(fout);
//					TAB5;
//					if (b10)
//					{
//						for (i = 0; i != ydj->njoint; i++)
//						{
//							fprintf(fout, "%f ", ydj->d1jsdc[i]);
//						}
//					}
//					else
//					{
//						writesheardisplacementb64(fout, ydj);
//					}
//					CHRwcr(fout);
//					TAB4; CHRw(fout, "</DataArray>");
//				}
//
//				CHRwcr(fout);
//				TAB4; CHRw(fout, "<DataArray type=\"Float64\" Name=\"Failure mode\" format=\"");
//				if (app)
//				{
//					CHRw(fout, "appended");
//				}
//				else if (b10)
//				{
//					CHRw(fout, "ascii");
//				}
//				else
//				{
//					CHRw(fout, "binary");
//				}
//				//CHRw(fout, "\" RangeMin=\"");
//				//fprintf(fout, "%f", d1min(ydj->d1jfmd, ydj->njoint));
//				//CHRw(fout, "\" RangeMax=\"");
//				//fprintf(fout, "%f", d1max(ydj->d1jfmd, ydj->njoint));
//
//				if (app)
//				{
//					CHRw(fout, "\" offset=\"");
//					fprintf(fout, off);
//					CHRw(fout, "\"/>");
//					off += sizeof(INS);
//					off += sizeof(DBL) * ydj->njoint;
//				}
//				else
//				{
//					CHRw(fout, "\">"); CHRwcr(fout);
//					TAB5;
//					if (b10)
//					{
//						for (i = 0; i != ydj->njoint; i++)
//						{
//							fprintf(fout, "%f ", ydj->d1jfmd[i]);
//						}
//					}
//					else if (b64)
//					{
//						writefailuremodeb64(fout, ydj);
//					
//					}
//					CHRwcr(fout);
//					TAB4; CHRw(fout, "</DataArray>");
//				}
//				CHRwcr(fout);
//
//				TAB4; CHRw(fout, "<DataArray type=\"Int64\" Name=\"Fracture type\" format=\"");
//				if (app)
//				{
//					CHRw(fout, "appended");
//				}
//				else if (b10)
//				{
//					CHRw(fout, "ascii");
//				}
//				else
//				{
//					CHRw(fout, "binary");
//				}
//				//CHRw(fout, "\" RangeMin=\"");
//				//fprintf("%d", i1min(ydp->i1pefr, ydp->nprop));
//				//CHRw(fout, "\" RangeMax=\"");
//				//fprintf("%d", i1max(ydp->i1pefr, ydp->nprop));
//				if (app)
//				{
//					CHRw(fout, "\" offset=\"");
//					fprintf(fout, off);
//					CHRw(fout, "\"/>");
//					CHRwcr(fout);
//					off += sizeof(INS);
//					off += sizeof(INT) * ydj->njoint;
//				}
//				else
//				{
//					CHRw(fout, "\">"); CHRwcr(fout);
//					TAB5;
//					if (b10)
//					{
//						for (i = 0; i != ydp->nprop; i++)
//						{
//							fprintf(fout, "%d ", ydp->i1pefr[i]);
//						}
//					}
//					else
//					{
//						writefracturetype(fout, ydp);
//					}
//					CHRwcr(fout);
//					TAB4; CHRw(fout, "</DataArray>");
//				}
//
//				CHRwcr(fout);
//				TAB4; CHRw(fout, "<DataArray type=\"Float64\" Name=\"Contact force\" NumberOfComponents=\"");
//				UINTw(fout, 3, 1);
//				CHRw(fout, "\" format=\"");
//				if (app)
//				{
//					CHRw(fout, "appended");
//				}
//				else if (b10)
//				{
//					CHRw(fout, "ascii");
//				}
//				else
//				{
//					CHRw(fout, "binary");
//				}
//				if (app)
//				{
//					CHRw(fout, "\" offset=\"");
//					fprintf(fout, off);
//					CHRw(fout, "\"/>");
//					off += sizeof(INS);
//					off += sizeof(INT) * yde->nelem;
//				}
//				else
//				{
//					CHRw(fout, "\">"); CHRwcr(fout);
//					TAB5;
//					if (b10)
//					{
//						i = 0;
//						while(i < ydj->njoint)
//						{
//							fprintf(fout, "%f ", ydj->d1jsdc[i]);
//							//ydn->d2nfcon
//							i++;
//							fprintf(fout, "%f ", ydj->d1jsdc[i]);
//							i++;
//							fprintf(fout, "%f ", 0.0);
//						}
//					}
//					else
//					{
//						writecontactforceb64(fout, ydn);
//					}
//					CHRwcr(fout);
//					TAB4; CHRw(fout, "</DataArray>");
//				}
//				CHRwcr(fout);
//				TAB3; CHRw(fout, "</PointData>"); CHRwcr(fout);

				TAB3; CHRw(fout, "<Points>"); CHRwcr(fout);
				TAB4; CHRw(fout, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"");

				if (app)
				{
					CHRw(fout, "appended");
				}
				else if (b10)
				{
					CHRw(fout, "ascii");
				}
				else
				{
					CHRw(fout, "binary");
				}
				CHRw(fout, "\" RangeMin=\"");
				fprintf(fout, "%f", d2min(ydn->d2ncc, ydn->nnodim, ydn->nnopo));
				CHRw(fout, "\" RangeMax=\"");
				fprintf(fout, "%f", d2max(ydn->d2ncc, ydn->nnodim, ydn->nnopo));

				if (app)
				{
					CHRw(fout, "\" offset=\"");
					UINTw(fout, off, 0);
					CHRw(fout, "\"/>"); CHRwcr(fout);
					off += sizeof(INS);
					off += sizeof(FLT) * (ydn->nnodim + 1) * ydn->nnopo;
				}
				else
				{
					CHRw(fout, "\">"); CHRwcr(fout);
					TAB5;
					if (b10)
					{
						for (i = 0; i != ydn->nnopo; i++)
						{
							for (j = 0; j != ydn->nnodim; j++)
							{
								fprintf(fout, "%f ", ydn->d2ncc[j][i]);
							}
							fprintf(fout, "0 ");
						}
					}
					else
					{
						writecoordsb64(fout, ydn);
					}
					CHRwcr(fout);
					TAB4; CHRw(fout, "</DataArray>"); CHRwcr(fout);
				}

				TAB3; CHRw(fout, "</Points>"); CHRwcr(fout);
				TAB3; CHRw(fout, "<Cells>"); CHRwcr(fout);

				TAB4; CHRw(fout, "<DataArray type=\"Int64\" Name=\"connectivity\" RangeMin=\"\" RangeMax=\"\" format=\"");
				if (app)
				{
					CHRw(fout, "appended\" offset=\"");
					UINTw(fout, off, 0);
					CHRw(fout, "\"/>"); CHRwcr(fout);
					off += sizeof(INS);
					off += sizeof(INT) * yde->nelem * yde->nelno;
				}
				else
				{
					if (b10)
					{
						CHRw(fout, "ascii\">"); CHRwcr(fout);
						TAB5;
						for (ielem = 0; ielem != yde->nelem; ielem++)
							for (i = 0; i != yde->nelno; i++)
								fprintf(fout, "%d ", yde->i2elto[i][ielem]);
					}
					else
					{
						CHRw(fout, "binary\">"); CHRwcr(fout);
						TAB5;
						writeconnectivityb64(fout, yde);
					}
					CHRwcr(fout);
					TAB4; CHRw(fout, "</DataArray>"); CHRwcr(fout);
				}

				TAB4; CHRw(fout, "<DataArray type=\"UInt32\" Name=\"offsets\" RangeMin=\"\" RangeMax=\"\" format=\"");

				if (app)
				{
					CHRw(fout, "appended\" offset=\"");
					INTw(fout, off, 0);
					CHRw(fout, "\"/>"); CHRwcr(fout);
					off += sizeof(INS);
					off += sizeof(UINT) * yde->nelem;
				}
				else
				{
					if (b10)
					{
						CHRw(fout, "ascii\">"); CHRwcr(fout);
						TAB5;
						for (ielem = 0; ielem != yde->nelem; ielem++)
							fprintf(fout, "%d ", (ielem + 1) * yde->nelno);
					}
					else
					{
						CHRw(fout, "binary\">"); CHRwcr(fout);
						TAB5;
						writeoffsetsb64(fout, yde);
					}
					CHRwcr(fout);
					TAB4; CHRw(fout, "</DataArray>"); CHRwcr(fout);
				}

				TAB4; CHRw(fout, "<DataArray type=\"UInt32\" Name=\"types\" RangeMin=\"\" RangeMax=\"\" format=\"");

				if (app)
				{
					CHRw(fout, "appended\" offset=\"");
					UINTw(fout, off, 0);
					CHRw(fout, "\"/>"); CHRwcr(fout);
					off += sizeof(INS);
					off += sizeof(UINT) * yde->nelem;
				}
				else
				{
					if (b10)
					{
						CHRw(fout, "ascii\">"); CHRwcr(fout);
						TAB5;
						for (ielem = 0; ielem != yde->nelem; ielem++)
							fprintf(fout, "%d ", 5);
					}
					else
					{
						CHRw(fout, "binary\">"); CHRwcr(fout);
						TAB5;
						writetypesb64(fout, yde);
					}
					CHRwcr(fout);
					TAB4; CHRw(fout, "</DataArray>"); CHRwcr(fout);
				}
				TAB3; CHRw(fout, "</Cells>"); CHRwcr(fout);
				TAB2; CHRw(fout, "</Piece>"); CHRwcr(fout);
				TAB; CHRw(fout, "</UnstructuredGrid>"); CHRwcr(fout);

				if (app)
				{
					TAB; CHRw(fout, "<AppendedData");
					if (b64)
					{
						CHRw(fout, " encoding = \"base64\"");
						CHRw(fout, ">"); CHRwcr(fout);
						TAB; CHRwsp(fout); CHRw(fout, "_");
						//writesheardisplacementb64(fout, ydj);
						//writefailuremodeb64(fout, ydj);
						//writefracturetype(fout, ydp);
						//writecontactforceb64(fout, ydj);

						writecoordsb64(fout, ydn);
						writeconnectivityb64(fout, yde);
						writeoffsetsb64(fout, yde);
						writetypesb64(fout, yde);
					}
					else
					{
						CHRw(fout, " encoding = \"raw\"");
						CHRw(fout, ">"); CHRwcr(fout);
						TAB; CHRwsp(fout); CHRw(fout, "_");

						//shear displacement
						i2i32(fout, sizeof(DBL) * ydj->njoint);
						for (i = 0; i != ydj->njoint; i++)
						{
							d2d64(fout, ydj->d1jsdc[i]);
						}

						//coords
						i2i32(fout, (INS) sizeof(DBL) * (ydn->nnodim + 1) * ydn->nnopo);
						for (i = 0; i != ydn->nnopo; i++)
						{
							for (j = 0; j != ydn->nnodim; j++)
							{
								d2d64(fout, ydn->d2ncc[j][i]);
							}
							d2d64(fout, 0.0);
						}

						//connectivity
						i2i32(fout, (INS) sizeof(INT) * yde->nelem * yde->nelno);
						for (ielem = 0; ielem != yde->nelem; ielem++)
							for (i = 0; i != yde->nelno; i++)
								i2i64(fout, yde->i2elto[i][ielem]);

						//connectivity offsets
						i2i32(fout, (INS) sizeof(INT) * yde->nelem);
						for (ielem = 0; ielem != yde->nelem; ielem++)
						{
							i2i64(fout, (ielem + 1) * yde->nelno);
						}

						//element type
						i2i32(fout, (INS)yde->nelem);
						for (ielem = 0; ielem != yde->nelem; ielem++)
							i2ui8(fout, 5);
					}
					CHRwcr(fout);
					TAB; CHRw(fout, "</AppendedData>"); CHRwcr(fout);
				}

				CHRw(fout, "</VTKFile>"); CHRwcr(fout);
			}
		}

		/*  (MT)  */

		fclose(fout);

		/* output total contact force */
		if (ydc->icouti == 0)
		{
			fcnf = fopen("contactforce.txt", "w");
		}
		else
		{
			fcnf = fopen("contactforce.txt", "a");
		}

		Yod2CONTACTFORCE
		(/* output contact forces */
			fcnf,
			ydn->nnopo, ydn->d2nfcon[0], ydn->d2nfcon[1], ydn->i1nopr);

		fclose(fcnf);

		ydc->icouti = ydc->icouti + 1;
	}
}

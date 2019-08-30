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
static INT  i1num[100];    /* numbers for space saving format     */
static DBL  d1num[100];    /* numbers for space saving format     */
static CHR c1code[500];    /* coded i1para in space saving format */

static void Yod2TRIELS(  /* small strain softening triangle output */
            nelem,
            fout,
            dcsizc,dcsizs,dcsizv,dcsizf,
            dpeks ,dpela ,dpemu ,dpero ,
            icoutp,iprop , 
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,
            i1elpr,i2elto,d2elcf
            )
  INT    nelem;
  FILE   *fout;
  DBL   dcsizc; DBL  dcsizs; DBL   dcsizv; DBL   dcsizf;
  DBL    dpeks; DBL   dpela; DBL    dpemu; DBL    dpero;
  INT   icoutp; INT   iprop;
  DBL  *d1nccx; DBL *d1nccy; DBL  *d1ncix; DBL  *d1nciy; DBL *d1nvcx;
  DBL  *d1nvcy;
  INT *i1elpr; INT **i2elto; DBL **d2elcf;
{ DBL voli,volc;
  DBL  B[2][2]; /* left Cauchy-Green strain tensor */
  DBL  D[2][2]; /* rate of deformation (stretching) tensor */
  DBL  E[2][2]; /* strain tensor (small strains) */
  DBL  F[2][2]; /* deformation gradient in global base */
  DBL F0[2][2]; /* initial local base */
  DBL FX[2][2]; /* current local base */
  DBL F0inv[2][2]; /* global base in initial local base */
  DBL FXinv[2][2]; /* global base in current local base */
  DBL  L[2][2]; /* velocity gradient in global base */
  DBL LX[2][2]; /* vel. gradient in current local base = delta x/delta X */
  DBL  T[2][2]; /* Cauchy stress */
  INT ielem;
  INT i,j,k;


  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { /* evaluate stress state */
      for(i=1;i<3;i++)
      { F0[0][i-1]=d1ncix[(i2elto[i][ielem])]-d1ncix[(i2elto[0][ielem])]; 
        F0[1][i-1]=d1nciy[(i2elto[i][ielem])]-d1nciy[(i2elto[0][ielem])];
        FX[0][i-1]=d1nccx[(i2elto[i][ielem])]-d1nccx[(i2elto[0][ielem])]; 
        FX[1][i-1]=d1nccy[(i2elto[i][ielem])]-d1nccy[(i2elto[0][ielem])];
        LX[0][i-1]=d1nvcx[(i2elto[i][ielem])]-d1nvcx[(i2elto[0][ielem])];
        LX[1][i-1]=d1nvcy[(i2elto[i][ielem])]-d1nvcy[(i2elto[0][ielem])];
      }
      YMATINV2(F0,F0inv,voli);
      YMATINV2(FX,FXinv,volc);
      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { F[i][j]=R0;
          L[i][j]=R0;
          for(k=0;k<2;k++)
          { F[i][j]=F[i][j]+FX[i][k]*F0inv[k][j];
            L[i][j]=L[i][j]+LX[i][k]*FXinv[k][j];
      } } }
      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { B[i][j]=R0;
          for(k=0;k<2;k++)
          { B[i][j]=B[i][j]+F[i][k]*F[j][k]; /* left Cauchy-Green strain */
          }
          D[i][j]=RP5*(L[i][j]+L[j][i]);     /* rate of deformation      */
          if(i==j)
          { E[i][j]=RP5*(B[i][j]-R1);        /* small strain             */
          }
          else
          { E[i][j]=RP5*B[i][j];
      } } }
      for(i=0;i<2;i++)     /* Cauchy stress */
      { for(j=0;j<2;j++)
        { T[i][j]=(R2*dpemu*E[i][j])*(voli/volc)+dpeks*D[i][j];
          if(i==j)T[i][j]=T[i][j]+dpela*(volc/voli-voli/volc);
      } } 
      /* prepare output */
      i1num[0]=icoutp;
      i1num[1]=17;
      d1num[3]=d1nccx[i2elto[0][ielem]]/dcsizc;
      d1num[4]=d1nccx[i2elto[1][ielem]]/dcsizc;
      d1num[5]=d1nccx[i2elto[2][ielem]]/dcsizc;
      d1num[6]=d1nccy[i2elto[0][ielem]]/dcsizc;
      d1num[7]=d1nccy[i2elto[1][ielem]]/dcsizc;
      d1num[8]=d1nccy[i2elto[2][ielem]]/dcsizc;
      d1num[9]=(d1nvcx[i2elto[0][ielem]]+
                d1nvcx[i2elto[1][ielem]]+
                d1nvcx[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[10]=(d1nvcy[i2elto[0][ielem]]+
                 d1nvcy[i2elto[1][ielem]]+
                 d1nvcy[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[11]=T[0][0]/dcsizs;
      d1num[12]=T[1][1]/dcsizs;
      d1num[13]=T[0][1]/dcsizs;
      d1num[14]=R0;             /* elastic damage */
	  d1num[15]=d2elcf[ielem][0]/dcsizf;
	  d1num[16]=d2elcf[ielem][1]/dcsizf;
      for(i=3;i<17;i++)
      { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));        
      }
      /* translate into INT */
      codeDBLtoINT(d1num,i1num); 
      i1num[2]=YTE2TRIELS;
      codeINTtoCHR(c1code,i1num);
      CHRw(fout,c1code); CHRwcr(fout);
} } }

static void Yod2TRISOF(  /* small strain softening triangle output */
            nelem,
            fout,
            dcsizc,dcsizs,dcsizv,
            dpeks ,dpela ,dpemu ,dpero ,
            icoutp,iprop , 
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,d1sdel,i1elpr,i2elto
            )
  INT    nelem;
  FILE   *fout;
  DBL   dcsizc; DBL  dcsizs;  DBL   dcsizv;
  DBL    dpeks; DBL   dpela; DBL   dpemu; DBL    dpero;
  INT   icoutp; INT   iprop;
  DBL  *d1nccx; DBL *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nvcx;
  DBL  *d1nvcy; DBL *d1sdel; INT *i1elpr; INT **i2elto;
{ DBL voli,volc;
  DBL  B[2][2]; /* left Cauchy-Green strain tensor */
  DBL  D[2][2]; /* rate of deformation (stretching) tensor */
  DBL  E[2][2]; /* strain tensor (small strains) */
  DBL  F[2][2]; /* deformation gradient in global base */
  DBL F0[2][2]; /* initial local base */
  DBL FX[2][2]; /* current local base */
  DBL F0inv[2][2]; /* global base in initial local base */
  DBL FXinv[2][2]; /* global base in current local base */
  DBL  L[2][2]; /* velocity gradient in global base */
  DBL LX[2][2]; /* vel. gradient in current local base = delta x/delta X */
  DBL  T[2][2]; /* Cauchy stress */
  INT ielem;
  INT i,j,k;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { /* evaluate stress state */
      for(i=1;i<3;i++)
      { F0[0][i-1]=d1ncix[(i2elto[i][ielem])]-d1ncix[(i2elto[0][ielem])]; 
        F0[1][i-1]=d1nciy[(i2elto[i][ielem])]-d1nciy[(i2elto[0][ielem])];
        FX[0][i-1]=d1nccx[(i2elto[i][ielem])]-d1nccx[(i2elto[0][ielem])]; 
        FX[1][i-1]=d1nccy[(i2elto[i][ielem])]-d1nccy[(i2elto[0][ielem])];
        LX[0][i-1]=d1nvcx[(i2elto[i][ielem])]-d1nvcx[(i2elto[0][ielem])];
        LX[1][i-1]=d1nvcy[(i2elto[i][ielem])]-d1nvcy[(i2elto[0][ielem])];
      }  
      YMATINV2(F0,F0inv,voli);
      YMATINV2(FX,FXinv,volc); 
      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { F[i][j]=R0;
          L[i][j]=R0;
          for(k=0;k<2;k++)
          { F[i][j]=F[i][j]+FX[i][k]*F0inv[k][j];
            L[i][j]=L[i][j]+LX[i][k]*FXinv[k][j];
      } } }
      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { B[i][j]=R0;
          for(k=0;k<2;k++)
          { B[i][j]=B[i][j]+F[i][k]*F[j][k]; /* left Cauchy-Green strain */
          }
          D[i][j]=RP5*(L[i][j]+L[j][i]);     /* rate of deformation      */
          if(i==j)
          { E[i][j]=RP5*(B[i][j]-R1);        /* small strain             */
          }
          else
          { E[i][j]=RP5*B[i][j];
      } } }
      for(i=0;i<2;i++)     /* Cauchy stress */
      { for(j=0;j<2;j++)
        { T[i][j]=(R1-d1sdel[ielem])*
                  (R2*dpemu*E[i][j])*(voli/volc)+dpeks*D[i][j];
          if(i==j)T[i][j]=T[i][j]+dpela*(volc/voli-voli/volc);
      } } 
      /* prepare output */
      i1num[0]=icoutp;
      i1num[1]=15;
      d1num[3]=d1nccx[i2elto[0][ielem]]/dcsizc;
      d1num[4]=d1nccx[i2elto[1][ielem]]/dcsizc;
      d1num[5]=d1nccx[i2elto[2][ielem]]/dcsizc;
      d1num[6]=d1nccy[i2elto[0][ielem]]/dcsizc;
      d1num[7]=d1nccy[i2elto[1][ielem]]/dcsizc;
      d1num[8]=d1nccy[i2elto[2][ielem]]/dcsizc;
      d1num[9]=(d1nvcx[i2elto[0][ielem]]+
                d1nvcx[i2elto[1][ielem]]+
                d1nvcx[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[10]=(d1nvcy[i2elto[0][ielem]]+
                 d1nvcy[i2elto[1][ielem]]+
                 d1nvcy[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[11]=T[0][0]/dcsizs;
      d1num[12]=T[1][1]/dcsizs;
      d1num[13]=T[0][1]/dcsizs;
      d1num[14]=d1sdel[ielem];
      for(i=3;i<15;i++)
      { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));        
      }
      /* translate into INT */
      codeDBLtoINT(d1num,i1num); 
      i1num[2]=YTE2TRIELS;
      codeINTtoCHR(c1code,i1num);
      CHRw(fout,c1code); CHRwcr(fout);
} } }

static void Yod2TRIRIG(  /* small strain elastic triangle output */
            nelem,
            fout,
            dcsizc,dcsizv,
            icoutp,iprop , 
            d1nccx,d1nccy,d1nvcx,d1nvcy,i1elpr,
            i2elto
            ) 
  INT    nelem;
  FILE   *fout;
  DBL   dcsizc; DBL   dcsizv;
  INT   icoutp; INT   iprop;
  DBL  *d1nccx; DBL *d1nccy; DBL *d1nvcx; DBL  *d1nvcy; INT *i1elpr;
  INT **i2elto;
{ INT ielem;
  INT i;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { /* prepare output */
      i1num[0]=icoutp;
      i1num[1]=15;
      d1num[3]=d1nccx[i2elto[0][ielem]]/dcsizc;
      d1num[4]=d1nccx[i2elto[1][ielem]]/dcsizc;
      d1num[5]=d1nccx[i2elto[2][ielem]]/dcsizc;
      d1num[6]=d1nccy[i2elto[0][ielem]]/dcsizc;
      d1num[7]=d1nccy[i2elto[1][ielem]]/dcsizc;
      d1num[8]=d1nccy[i2elto[2][ielem]]/dcsizc;
      d1num[9]=(d1nvcx[i2elto[0][ielem]]+
                d1nvcx[i2elto[1][ielem]]+
                d1nvcx[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[10]=(d1nvcy[i2elto[0][ielem]]+
                 d1nvcy[i2elto[1][ielem]]+
                 d1nvcy[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[11]=R0;
      d1num[12]=R0;
      d1num[13]=R0;
      d1num[14]=R0;
      for(i=3;i<15;i++)
      { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
      }
      /* translate into INT */
      codeDBLtoINT(d1num,i1num); 
      i1num[2]=YTE2TRIELS;
      codeINTtoCHR(c1code,i1num);
      CHRw(fout,c1code); CHRwcr(fout);
} } }

static void Yod2JOINTSBROKEN(  /* 2D joint output */
            nelem,
            fout,
            dcsizc,dcsizv,dcsizf,dcsizd,dcsiza,
            icoutp,iprop ,
            d1nccx,d1nccy,d1nvcx,d1nvcy,d1sdel,
            i1elpr,i2elto,d2tcs,dcsizs,i2eljp,
			i1elty,d1elsf,d2elcf,
			njoint,i1jtid,d1jknc,d1jksc,d1jnst,
			d1jsst,d1japc,d1japh,d1jsdc,d1jdlc,
			d1jphi,d1jfmd
            )
  INT   nelem;
  FILE  *fout;
  DBL  dcsizc; DBL   dcsizv; DBL  dcsizf; DBL  dcsizd; DBL   dcsiza;
  INT  icoutp; INT    iprop;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1nvcx; DBL *d1nvcy; DBL  *d1sdel;
  INT *i1elpr; INT **i2elto; DBL **d2tcs; DBL  dcsizs; INT **i2eljp;
  INT *i1elty; DBL  *d1elsf; DBL **d2elcf;
  INT  njoint; INT  *i1jtid; DBL *d1jknc; DBL *d1jksc; DBL  *d1jnst;
  DBL *d1jsst; DBL  *d1japc; DBL *d1japh; DBL *d1jsdc; DBL  *d1jdlc;
  DBL *d1jphi; DBL  *d1jfmd;
{ INT i,ielem,ijoint;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elty[ielem]==1)		  /* output boundary */
    { i1num[ 0]=icoutp;
	  i1num[ 1]=19;
      d1num[ 3]=d1nccx[i2elto[0][ielem]]/dcsizc;
      d1num[ 4]=d1nccx[i2elto[1][ielem]]/dcsizc;
      d1num[ 5]=d1nccx[i2elto[2][ielem]]/dcsizc;
      d1num[ 6]=d1nccx[i2elto[3][ielem]]/dcsizc;
      d1num[ 7]=d1nccy[i2elto[0][ielem]]/dcsizc;
      d1num[ 8]=d1nccy[i2elto[1][ielem]]/dcsizc;
      d1num[ 9]=d1nccy[i2elto[2][ielem]]/dcsizc;
      d1num[10]=d1nccy[i2elto[3][ielem]]/dcsizc;
      d1num[11]=d2elcf[i2eljp[0][ielem]][0]/dcsizf;
      d1num[12]=d2elcf[i2eljp[0][ielem]][1]/dcsizf;
      d1num[13]=R0;
	  d1num[14]=R0;
	  d1num[15]=R0;
	  d1num[16]=R0;
	  d1num[17]=R0;
	  d1num[18]=R0;
      for(i=3;i<i1num[1];i++)
      { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
      }
      /* translate into INT */
      codeDBLtoINT(d1num,i1num);
      i1num[2]=YTE2BOUNDS;
      codeINTtoCHR(c1code,i1num);
      CHRw(fout,c1code); CHRwcr(fout);
	}
	else if(i1elty[ielem]>1)	/* output fracture joints */
	{ ijoint=ielem-(nelem-njoint);
	  i1num[ 0]=icoutp;
	  i1num[ 1]=19;
	  d1num[ 3]=d1nccx[i2elto[0][ielem]]/dcsizc;
      d1num[ 4]=d1nccx[i2elto[1][ielem]]/dcsizc;
      d1num[ 5]=d1nccx[i2elto[2][ielem]]/dcsizc;
      d1num[ 6]=d1nccx[i2elto[3][ielem]]/dcsizc;
      d1num[ 7]=d1nccy[i2elto[0][ielem]]/dcsizc;
      d1num[ 8]=d1nccy[i2elto[1][ielem]]/dcsizc;
      d1num[ 9]=d1nccy[i2elto[2][ielem]]/dcsizc;
      d1num[10]=d1nccy[i2elto[3][ielem]]/dcsizc;
	  d1num[11]=d1jfmd[ijoint]/R2;
	  d1num[12]=d1jnst[ijoint]/dcsizs;
	  d1num[13]=d1jsst[ijoint]/dcsizs;
	  d1num[14]=(DABS(d1elsf[i2eljp[0][ielem]])+DABS(d1elsf[i2eljp[1][ielem]]))/R2/dcsizf;
	  d1num[15]=log10(d1japc[ijoint])/log10(dcsiza);
	  d1num[16]=DABS(d1jsdc[ijoint])/dcsizd;
	  d1num[17]=d1jdlc[ijoint]/dcsizd;
	  d1num[18]=log10(d1japh[ijoint])/log10(dcsiza);
      for(i=3;i<i1num[1];i++)
      { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
      }
      /* translate into INT */
      codeDBLtoINT(d1num,i1num);
      i1num[2]=YTE2JOINTS;
      codeINTtoCHR(c1code,i1num);
      CHRw(fout,c1code); CHRwcr(fout);
  } }
}

static void Yod2AESOURCE(	/* acoustic emission source */
			nelem,
			fout,
			dcsizc,
			icoutp,
			i2elto,d1elme,i1elyi,
			d1nccx,d1nccy
			)
  INT	 nelem;
  FILE   *fout;
  DBL   dcsizc;
  INT   icoutp;
  INT **i2elto; DBL *d1elme; INT *i1elyi;
  DBL  *d1nccx; DBL *d1nccy;
{ INT i,ielem;
  DBL dcsizm=15.0;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elyi[ielem]>=2)
	{ i1num[ 0]=icoutp;
	  i1num[ 1]=10;
	  d1num[ 3]=(d1nccx[i2elto[0][ielem]]+d1nccx[i2elto[1][ielem]]+
		d1nccx[i2elto[2][ielem]]+d1nccx[i2elto[3][ielem]])/R4/dcsizc;
      d1num[ 4]=(d1nccy[i2elto[0][ielem]]+d1nccy[i2elto[1][ielem]]+
		d1nccy[i2elto[2][ielem]]+d1nccy[i2elto[3][ielem]])/R4/dcsizc;
	  d1num[ 5]=d1elme[ielem]/dcsizm;
	  d1num[ 6]=R0;
	  d1num[ 7]=R0;
	  d1num[ 8]=R0;
	  d1num[ 9]=R0;
	  for(i=3;i<10;i++)
      { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
      }
	  /* translate into INT */
      codeDBLtoINT(d1num,i1num);
      i1num[2]=YTE2AESOUR;
      codeINTtoCHR(c1code,i1num);
      CHRw(fout,c1code); CHRwcr(fout);
  } }
}

static void Yod2JOINTSINTACT(  /* 2D joint output */
            nelem,
            fout,
            dcsizc,dcsizv,
            dpefs, dpeft, dpegf, dpeks, dpepe,
            icoutp,iprop , 
            d1nccx,d1nccy,d1nvcx,d1nvcy,d1sdel,
            i1elpr,i2elto
            ) 
  INT    nelem;
  FILE   *fout;
  DBL   dcsizc; DBL   dcsizv;
  DBL   dpefs; DBL   dpeft; DBL   dpegf; DBL   dpeks; DBL   dpepe;
  INT   icoutp; INT   iprop;
  DBL  *d1nccx; DBL *d1nccy; DBL *d1nvcx; DBL  *d1nvcy; DBL *d1sdel;
  INT  *i1elpr; INT **i2elto;
{ DBL small,o1,o2,s1,s2,op,sp,ot,st;
  DBL e1x,e1y,h;
  INT ielem,i,i0,i1,i2,i3;
  
  small=EPSILON; 
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { i0=i2elto[0][ielem];
      i1=i2elto[1][ielem];
      i2=i2elto[2][ielem];
      i3=i2elto[3][ielem];
      e1x=RP5*(d1nccx[i1]+d1nccx[i2]-d1nccx[i0]-d1nccx[i3]);
      e1y=RP5*(d1nccy[i1]+d1nccy[i2]-d1nccy[i0]-d1nccy[i3]);
      h=SQRT(e1x*e1x+e1y*e1y);
      e1x=e1x/(h+small);
      e1y=e1y/(h+small);
      s1=(d1nccy[i0]-d1nccy[i3])*e1y+(d1nccx[i0]-d1nccx[i3])*e1x;
      s2=(d1nccy[i1]-d1nccy[i2])*e1y+(d1nccx[i1]-d1nccx[i2])*e1x;
      o1=(d1nccy[i0]-d1nccy[i3])*e1x-(d1nccx[i0]-d1nccx[i3])*e1y;
      o2=(d1nccy[i1]-d1nccy[i2])*e1x-(d1nccx[i1]-d1nccx[i2])*e1y;
      op=R2*h*dpeft/dpepe;
      sp=R2*h*dpefs/dpepe;
      ot=MAXIM((R2*op),(R3*dpegf/dpeft));
      st=MAXIM((R2*sp),(R3*dpegf/dpefs));
      /* prepare output */
      if((((o1+o2)/(R2*ot))>RP1)||(((s1+s2)/(R2*st))>RP1))
      { i1num[0]=icoutp;
        i1num[1]=14;
        d1num[3]=d1nccx[i2elto[0][ielem]]/dcsizc;
        d1num[4]=d1nccx[i2elto[1][ielem]]/dcsizc;
        d1num[5]=d1nccx[i2elto[2][ielem]]/dcsizc;
        d1num[6]=d1nccx[i2elto[3][ielem]]/dcsizc;
        d1num[7]=d1nccy[i2elto[0][ielem]]/dcsizc;
        d1num[8]=d1nccy[i2elto[1][ielem]]/dcsizc;
        d1num[9]=d1nccy[i2elto[2][ielem]]/dcsizc;
        d1num[10]=d1nccy[i2elto[3][ielem]]/dcsizc;
        d1num[11]=R1;
        d1num[12]=o1/ot;
        d1num[13]=o2/st;
        for(i=3;i<14;i++)
        { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));        
        }
        /* translate into INT */
        codeDBLtoINT(d1num,i1num); 
        i1num[2]=YTE2JOINTS;
        codeINTtoCHR(c1code,i1num);
        CHRw(fout,c1code); CHRwcr(fout);  
} } } }

static void Yod2CONTACTFORCE( /* output contact forces */
		fcnf,
		nnopo,d1nfconx,d1nfcony,i1nopr
		)
  FILE *fcnf;
  INT  nnopo; DBL *d1nfconx; DBL *d1nfcony; INT *i1nopr;
{ INT inode;
  DBL dsumfx,dsumfy,tmp,tmpy;
  
  dsumfx=R0; dsumfy=R0;
  for(inode=0;inode<nnopo;inode++)
  { dsumfx=dsumfx+DABS(d1nfconx[inode]);
	dsumfy=dsumfy+DABS(d1nfcony[inode]);

        if(i1nopr[inode]==5)
        {
            tmp=tmp+d1nfconx[inode];
            tmpy=tmpy+d1nfcony[inode];
        }
  }

  DBLw(fcnf,dsumfx,10); CHRwsp(fcnf);
  DBLw(fcnf,dsumfy,10); CHRwsp(fcnf);
  DBLw(fcnf,tmp,10); CHRwsp(fcnf);
  DBLw(fcnf,tmpy,10); CHRwcr(fcnf);
}

void Yod(  
     namep,yd
     )
  CHR *namep; YD yd;
{ YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  YDJ ydj=&(yd->ydj);
  INT iprop,ihys,i;
  CHR namef[300];
  CHR cindex[50];
  static INT ncall=0;
  FILE *fout=FILENULL;
  FILE *fcnf=FILENULL;
  DBL tmp,tmpy;

  ncall=ncall+1;
  /* output hystory */
  for(ihys=0;ihys<(ydo->nohys);ihys++)
  { if((ydo->d1ohyt[ihys])==(ydc->dctime))
    { if((ydo->f2ohyf[ihys])==FILENULL)
      { CHRcpy(namef,namep);
        SINTw(cindex,ihys,0);
        CHRcat(namef,"h");
        CHRcat(namef,cindex); 
        ydo->f2ohyf[ihys]=fopen(namef,"a");
      }
      if((ydo->f2ohyf[ihys])!=FILENULL)
      { tmp=(ydo->d1ohyt[ihys])*(ydo->d1ohyc[ihys]);
        DBLw((ydo->f2ohyf[ihys]),tmp,10); 
        CHRwsp(ydo->f2ohyf[ihys]);
        tmp=(ydo->d1ohys[ihys])*(ydo->d1ohyf[ihys]);
        DBLw((ydo->f2ohyf[ihys]),tmp,10);
        CHRwcr(ydo->f2ohyf[ihys]);


        tmp=R0; tmpy=R0;
        for(i=0;i<ydn->nnopo;i++)
        {
        if(ydn->i1nopr[i]==1)
        {
            tmp=tmp+ydn->d2nfcon[0][i];
            tmpy=tmpy+ydn->d2nfcon[1][i];
        }
        }
        CHRwsp(ydo->f2ohyf[ihys]);
        DBLw((ydo->f2ohyf[ihys]),tmp,10);
        CHRwsp(ydo->f2ohyf[ihys]);
        DBLw((ydo->f2ohyf[ihys]),tmpy,10);
        CHRwcr(ydo->f2ohyf[ihys]);       
  } } 


}
  if((ncall>100)||((ydc->ncstep)>=(ydc->mcstep-2)))
  { ncall=0;
    for(ihys=0;ihys<(ydo->nohys);ihys++)
    { if((ydo->f2ohyf[ihys])!=FILENULL)
	  { fclose(ydo->f2ohyf[ihys]);
	  }
      ydo->f2ohyf[ihys]=FILENULL;
  } }
  /* output anymation */
  fout=FILENULL;
  fcnf=FILENULL;
  if((ydc->ncstep%ydc->icoutf)==0)
  { CHRcpynoext(namef,namep);
    SINTw(cindex,ydc->icouti,0);
    CHRcat(namef,cindex);
    CHRcat(namef,".ym");
    fout=fopen(namef,"w");
    if(fout!=FILENULL)
    { CHRw(fout,"CODED"); CHRwcr(fout);
      for(iprop=0;iprop<ydp->nprop;iprop++)
      { if((ydp->i1ptyp[iprop])==(YTE2TRIELS))
        { Yod2TRIELS(  
          yde->nelem ,
          fout       ,
          ydc->dcsizc,ydc->dcsizs,ydc->dcsizv,ydc->dcsizf,
          ydp->d1peks[iprop],ydp->d1pela[iprop],ydp->d1pemu[iprop],
          ydp->d1pero[iprop],
          ydc->icoutp       ,iprop             ,                     
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],
          yde->i1elpr,yde->i2elto,yde->d2elcf
          );
        }
        else if((ydp->i1ptyp[iprop])==(YTE2TRISOF))
        { Yod2TRISOF(  /* small strain elastic triangle output */
          yde->nelem ,
          fout       ,
          ydc->dcsizc       ,ydc->dcsizs       ,ydc->dcsizv       ,
          ydp->d1peks[iprop],ydp->d1pela[iprop],ydp->d1pemu[iprop],
          ydp->d1pero[iprop],
          ydc->icoutp       ,iprop             ,                     
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],yde->d2elst[ydp->i1psde[iprop]]    ,yde->i1elpr,
          yde->i2elto
          );
        }
        else if((ydp->i1ptyp[iprop])==(YTE2TRIRIG))
        { Yod2TRIRIG(  /* rigid triangle output */
          yde->nelem ,
          fout       ,
          ydc->dcsizc       ,ydc->dcsizv       ,
          ydc->icoutp       ,iprop             ,                    
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nvc[0],ydn->d2nvc[1],yde->i1elpr,
          yde->i2elto
          );
        }
        else if((ydp->i1ptyp[iprop])==(YTE2JOINTS))
        { //Yod2JOINTSINTACT(  /* 2D joint output */
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
          Yod2JOINTSBROKEN(	/* broken joints */
          yde->nelem ,
          fout       ,
          ydc->dcsizc,ydc->dcsizv,ydc->dcsizf,ydc->dcsizd,ydc->dcsiza,
          ydc->icoutp       ,iprop             ,                    
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nvc[0],ydn->d2nvc[1],
          yde->d2elst[ydp->i1psde[iprop]],
          yde->i1elpr,yde->i2elto,yde->d2tcs,ydc->dcsizs,yde->i2eljp,
		  yde->i1elty,yde->d1elsf,yde->d2elcf,
		  ydj->njoint,ydj->i1jtid,ydj->d1jknc,ydj->d1jksc,ydj->d1jnst,
		  ydj->d1jsst,ydj->d1japc,ydj->d1japh,ydj->d1jsdc,ydj->d1jdlc,
		  ydj->d1jphi,ydj->d1jfmd
          );
		  Yod2AESOURCE(	/* acoustic emission source */
		  yde->nelem ,
		  fout       ,
		  ydc->dcsizc,
		  ydc->icoutp,
		  yde->i2elto,yde->d1elme,yde->i1elyi,
		  ydn->d2ncc[0],ydn->d2ncc[1]
		  );
        }
      }
	  fclose(fout);

	  /* output total contact force */
	  if(ydc->icouti==0)
	  { fcnf=fopen("contactforce.txt","w");
	  }
	  else
	  { fcnf=fopen("contactforce.txt","a");
	  }
	  Yod2CONTACTFORCE( /* output contact forces */
	  fcnf,
	  ydn->nnopo,ydn->d2nfcon[0],ydn->d2nfcon[1],ydn->i1nopr
	  );
	  fclose(fcnf);

      ydc->icouti=ydc->icouti+1;
} } }

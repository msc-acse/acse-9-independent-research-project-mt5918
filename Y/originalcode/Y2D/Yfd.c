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
  /* File   Yfd.c */
#include "Yproto.h"

/**************GENERALISED NODAL FORCES***********/ 
static void Yfd2TRIELS(  /* small strain elastic triangle  */
			dctime,dcgrst,dcrmpt,
             nelem,njoint,
             iprop,
             dpeks, dpela, dpemu,
			 dpero, dpsem,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nfcx,
            d1nfcy,d1nmct,d1nvcx,d1nvcy,d1pnaf,
            d1pnap,d1pnat,i1elpr,i2eljp,i1elty,
			i1nopr,i2elto, d2tcs,d1jfpr
            )
  DBL  dctime; DBL   dcgrst; DBL  dcrmpt;
  INT   nelem; INT   njoint;
  INT   iprop;
  DBL   dpeks; DBL    dpela; DBL   dpemu;
  DBL   dpero; DBL    dpsem;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL  *d1nfcx;
  DBL *d1nfcy; DBL  *d1nmct; DBL *d1nvcx; DBL  *d1nvcy; DBL  *d1pnaf;
  DBL *d1pnap; DBL  *d1pnat; INT *i1elpr; INT **i2eljp; INT  *i1elty;
  INT *i1nopr; INT **i2elto; DBL **d2tcs; DBL  *d1jfpr;
{ DBL nx,ny,voli,volc;
  DBL  B[2][2];			/* left Cauchy-Green strain tensor */
  DBL  D[2][2];			/* rate of deformation (stretching) tensor */
  DBL  E[2][2];			/* strain tensor (small strains) */
  DBL  F[2][2];			/* deformation gradient in global base */
  DBL F0[2][2];			/* initial local base */
  DBL FX[2][2];			/* current local base */
  DBL F0inv[2][2];		/* global base in initial local base */
  DBL FXinv[2][2];		/* global base in current local base */
  DBL  L[2][2];			/* velocity gradient in global base */
  DBL LX[2][2];			/* vel. gradient in current local base = delta x/delta X */
  DBL  T[2][2];			/* Cauchy stress */
  DBL   dredf;			/* ramping reduction factor */
  INT ielem,ijoint,i,j,k,in,jn,kn,jnopr,knopr;	/* loop variables */

  /* staged loading */
  if(dctime<dcgrst)
  { dredf=R0;
  }
  else if(dctime<dcgrst+dcrmpt)
  { dredf=(dctime-dcgrst)/dcrmpt;
  }
  else
  { dredf=R1;
  }

  /* finite strain computation */
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { for(i=1;i<3;i++)
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
          E[i][j]=RP5*B[i][j];               /* small strain             */
          if(i==j)E[i][j]=E[i][j]-RP5;
      } }
      for(i=0;i<2;i++)     /* Cauchy stress T=2*mu*E+lamda*trace(E) */
      { for(j=0;j<2;j++)
        { T[i][j]=R2*dpemu*E[i][j]*(voli/volc)+dpeks*D[i][j];
        }
        T[i][i]=T[i][i]+dpela*(volc/voli-voli/volc);
      }
      d2tcs[0][ielem]=T[0][0];
      d2tcs[1][ielem]=T[0][1];
      d2tcs[2][ielem]=T[1][1];

      for(i=0;i<3;i++)      /* Nodal Forces */
      { j=i+1; if(j>2)j=0;
        k=j+1; if(k>2)k=0;
        in=i2elto[i][ielem];
        jn=i2elto[j][ielem];
        kn=i2elto[k][ielem];
        nx=d1nccy[kn]-d1nccy[jn];
        ny=d1nccx[jn]-d1nccx[kn];
        d1nmct[in]=d1nmct[in]+dpero*voli/R6;
        d1nfcx[in]=d1nfcx[in]+(T[0][0]*nx+T[0][1]*ny)/R2;
        d1nfcy[in]=d1nfcy[in]+(T[1][0]*nx+T[1][1]*ny)/R2;
        /* Nodal Forces due to edge force */
        jnopr=i1nopr[jn];
        knopr=i1nopr[kn];
        if( ((DABS(d1pnap[jnopr]))>EPSILON)&&
            ((DABS(d1pnap[knopr]))>EPSILON) )
        { d1nfcx[jn]=d1nfcx[jn]-
          dredf*d1pnap[jnopr]*d1pnaf[jnopr]*nx/R3-
          dredf*d1pnap[knopr]*d1pnaf[knopr]*nx/R6;
          d1nfcy[jn]=d1nfcy[jn]-
          dredf*d1pnap[jnopr]*d1pnaf[jnopr]*ny/R3-
          dredf*d1pnap[knopr]*d1pnaf[knopr]*ny/R6;
          d1nfcx[kn]=d1nfcx[kn]-
          dredf*d1pnap[jnopr]*d1pnaf[jnopr]*nx/R6-
          dredf*d1pnap[knopr]*d1pnaf[knopr]*nx/R3;
          d1nfcy[kn]=d1nfcy[kn]-
          dredf*d1pnap[jnopr]*d1pnaf[jnopr]*ny/R6-
          dredf*d1pnap[knopr]*d1pnaf[knopr]*ny/R3;
        }
        if( ((DABS(d1pnat[jnopr]))>EPSILON)&&
            ((DABS(d1pnat[knopr]))>EPSILON) )
        { d1nfcx[jn]=d1nfcx[jn]-
          dredf*d1pnat[jnopr]*d1pnaf[jnopr]*ny/R3-
          dredf*d1pnat[knopr]*d1pnaf[knopr]*ny/R6;
          d1nfcy[jn]=d1nfcy[jn]+
          dredf*d1pnat[jnopr]*d1pnaf[jnopr]*nx/R3+
          dredf*d1pnat[knopr]*d1pnaf[knopr]*nx/R6;
          d1nfcx[kn]=d1nfcx[kn]-
          dredf*d1pnat[jnopr]*d1pnaf[jnopr]*ny/R6-
          dredf*d1pnat[knopr]*d1pnaf[knopr]*ny/R3;
          d1nfcy[kn]=d1nfcy[kn]+
          dredf*d1pnat[jnopr]*d1pnaf[jnopr]*nx/R6+
          dredf*d1pnat[knopr]*d1pnaf[knopr]*nx/R3;
        }
		/* Nodal Forces due to internal fluid pressure */
		if(i2eljp[i][ielem]>0)
		{ if(i1elty[i2eljp[i][ielem]]>1)
		  { ijoint=i2eljp[i][ielem]-(nelem-njoint);
			d1nfcx[in]=d1nfcx[in]-dredf*d1jfpr[ijoint]*nx/R2;
			d1nfcy[in]=d1nfcy[in]-dredf*d1jfpr[ijoint]*ny/R2;
			d1nfcx[jn]=d1nfcx[jn]-dredf*d1jfpr[ijoint]*nx/R2;
			d1nfcy[jn]=d1nfcy[jn]-dredf*d1jfpr[ijoint]*ny/R2;
		} }
	  }
} } }

static void Yfd2TRISOF(  /* small strain softening triangle  */
            nelem,
            iprop,
             dpefs,dpeks,dpela ,dpemu , dpero,
             dpsem,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nfcx,
            d1nfcy,d1nmct,d1nvcx,d1nvcy,d1sdel,
            i1elpr,i2elto
            )
  INT    nelem;
  INT    iprop;
  DBL   dpefs; DBL    dpeks; DBL    dpela; DBL   dpemu; DBL   dpero;
  DBL   dpsem;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nfcx;
  DBL *d1nfcy; DBL  *d1nmct; DBL *d1nvcx; DBL  *d1nvcy; DBL *d1sdel; 
  INT *i1elpr; INT **i2elto;
{ DBL nx,ny,voli,volc,emax,dmax,dmas;
  DBL  B[2][2];		/* left Cauchy-Green strain tensor */
  DBL  D[2][2];		/* rate of deformation (stretching) tensor */
  DBL  E[2][2];		/* strain tensor (small strains) */
  DBL  F[2][2];		/* deformation gradient in global base */
  DBL F0[2][2];		/* initial local base */
  DBL FX[2][2];		/* current local base */
  DBL F0inv[2][2];	/* global base in initial local base */
  DBL FXinv[2][2];	/* global base in current local base */
  DBL  L[2][2];		/* velocity gradient in global base */
  DBL LX[2][2];		/* vel. gradient in current local base = delta x/delta X */
  DBL  T[2][2];		/* Cauchy stress */
  INT ielem;
  INT i,j,k;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { for(i=1;i<3;i++)
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
          E[i][j]=RP5*B[i][j];               /* small strain             */
          if(i==j)E[i][j]=E[i][j]-RP5;  
      } }
      emax=MAXIM(R1,SQRT(((B[0][0]+B[1][1])+
                    SQRT((B[0][0]-B[1][1])*(B[0][0]-B[1][1])+
                    R4*B[1][0]*B[0][1]))*RP5));
      if(d1sdel[ielem]>(-EPSILON))
      { d1sdel[ielem]=MAXIM(d1sdel[ielem],((emax-R1)/dpsem));
        dmax=d1sdel[ielem];
        dmax=MINIM(R1,dmax*dmax*dmax);
        dmas=R1;
      }
      else
      { dmax=0.0; dmas=R1/(dpefs*dpefs);
      }
      for(i=0;i<2;i++)     /* Cauchy stress T=2*mu*E+lamda*trace(E) */
      { for(j=0;j<2;j++)
        { T[i][j]=(R1-dmax)*
                  (R2*dpemu*E[i][j]*(voli/volc)+dpeks*D[i][j]);
          if(i==j)
          { if(volc<voli)
            { T[i][j]=T[i][j]+dpela*(volc/voli-voli/volc);
            }
            else
            { T[i][j]=T[i][j]+(R1-dmax)*dpela*(volc/voli-voli/volc);
      } } } }
      for(i=0;i<3;i++)      /* Nodal Forces due to stress */
      { j=i+1; if(j>2)j=0;
        k=j+1; if(k>2)k=0;
        nx=d1nccy[(i2elto[k][ielem])]-d1nccy[(i2elto[j][ielem])];
        ny=d1nccx[(i2elto[j][ielem])]-d1nccx[(i2elto[k][ielem])];
        d1nmct[(i2elto[i][ielem])]=
              d1nmct[(i2elto[i][ielem])]+dmas*dpero*voli/R6;
        d1nfcx[(i2elto[i][ielem])]=d1nfcx[(i2elto[i][ielem])]+
			  (T[0][0]*nx+T[0][1]*ny)/R2;
        d1nfcy[(i2elto[i][ielem])]=d1nfcy[(i2elto[i][ielem])]+
			  (T[1][0]*nx+T[1][1]*ny)/R2; 
      }
} } }

static void Yfd2JRCM(	/* mobilised joint roughness */
			JRC, K, PHIR, ratio, JRCm
  )
  DBL JRC; DBL K; DBL PHIR; DBL ratio; DBL *JRCm;
{ DBL delta1, delta2, delta;
  if(ratio<=0.3)
  { delta1=-PHIR/(JRC*K); delta2=0.0;
	delta=delta2-(0.3-ratio)/(0.3-0.0)*(delta2-delta1);
  }
  else if(ratio<=0.6)
  { delta1=0.0; delta2=0.75;
	delta=delta2-(0.6-ratio)/(0.6-0.3)*(delta2-delta1);
  }
  else if(ratio<=1.0)
  { delta1=0.75; delta2=1.0;
	delta=delta2-(1.0-ratio)/(1.0-0.6)*(delta2-delta1);
  }
  else if(ratio<=2.0)
  { delta1=1.0; delta2=0.85;
	delta=delta2-(2.0-ratio)/(2.0-1.0)*(delta2-delta1);
  }
  else if(ratio<=4.0)
  { delta1=0.85; delta2=0.7;
	delta=delta2-(4.0-ratio)/(4.0-2.0)*(delta2-delta1);
  }
  else if(ratio<=10.0)
  { delta1=0.7; delta2=0.5;
	delta=delta2-(10.0-ratio)/(10.0-4.0)*(delta2-delta1);
  }
  else if(ratio<=100.0)
  { delta1=0.5; delta2=0.0;
	delta=delta2-(100.0-ratio)/(100.0-10.0)*(delta2-delta1);
  }
  else
  { delta=0.0;
  }
  *JRCm=delta*JRC;
}

static void Yfd2JOINTS(  /* joint element  */
			dctime,dcgrst,dcrmpt,
             nelem,
			 iprop,
             dpeft, dpegt, dpegs,
			 dpeks, dpepe, dpicf,
			 dpcoh, dpefr, dpepf,
            d1nccx,d1nccy,d1nfcx,d1nfcy,
			d1nvcx,d1nvcy,d1nmct,
			d1sdel,i1elpr,i2elto,i2eljp,
			d2elfs,i1elty, d2tcs,d1eley,
			d1eles,d1elme,i1elyi,d0iedi,
			njoint,d1jkni,d1jknc,d1jksc,d1jnst,
			d1jsst,d1japi,d1japc,d1japh,d1japr,
			d1jsdc,d1jdlc,d1jsdp,d1jefl,d1jjrc,
			d1jjcs,d1jphi,d1jfmd,d1jfpr,d1jfet
            )
  DBL	dctime; DBL  dcgrst; DBL   dcrmpt;
  INT    nelem;
  INT    iprop;
  DBL    dpeft; DBL   dpegt; DBL    dpegs;
  DBL    dpeks; DBL   dpepe; DBL    dpicf; 
  DBL    dpcoh; DBL   dpefr; DBL    dpepf;
  DBL  *d1nccx; DBL *d1nccy; DBL  *d1nfcx; DBL  *d1nfcy;
  DBL  *d1nvcx; DBL *d1nvcy; DBL  *d1nmct;
  DBL  *d1sdel; INT *i1elpr; INT **i2elto; INT **i2eljp;
  DBL **d2elfs; INT *i1elty; DBL **d2tcs ; DBL  *d1eley; 
  DBL  *d1eles; DBL *d1elme; INT  *i1elyi; DBL  *d0iedi;
  INT   njoint; DBL	*d1jkni; DBL  *d1jknc; DBL  *d1jksc; DBL *d1jnst;
  DBL  *d1jsst; DBL *d1japi; DBL  *d1japc; DBL  *d1japh; DBL *d1japr;
  DBL  *d1jsdc; DBL *d1jdlc; DBL  *d1jsdp; DBL  *d1jefl; DBL  *d1jjrc;
  DBL  *d1jjcs; DBL *d1jphi; DBL  *d1jfmd; DBL  *d1jfpr; DBL  *d1jfet;
{ DBL dpefa=0.63;
  DBL dpefb=1.8;
  DBL dpefc=6.0;
  DBL dpefm=0.0;
  DBL dpefs;
  DBL small,sabs,o,s,o1,o2,s1,s2,op,sp,ot,st,z,sigma,tau;
  DBL e1x,e1y,h,area;
  INT ielem,jelem,integ,i0,i1,i2,i3,nfail,nsoft;
  INT i,iside,ijoint,inode;
  DBL nx,ny,Sxx,Sxy,Syy,nstr[2],sstr[2];
  DBL a0,ar,vm,kni,knn,knt,ktt,JRC,JCS,PHII,PHIR,JRCm,dm,M,K;
  DBL u1,u2,du,dv,up,uratio,as1,as2,an1,an2,ah,dsigma,blklen;
  DBL sigma1,sigma2,tau1,tau2,taup;
  DBL Ek;
  DBL dredf;	// ramping reduction factor
  DBL dprff;	// percentage of fluid filling

  /* staged loading */
  if(dctime<dcgrst)
  { dredf=R0;
  }
  else if(dctime<dcgrst+dcrmpt)
  { dredf=(dctime-dcgrst)/dcrmpt;
  }
  else
  { dredf=R1;
  }

  small=EPSILON; nsoft=0;
  for(ielem=0;ielem<nelem;ielem++)
  { ijoint=ielem-(nelem-njoint);
	if(i1elty[ielem]==0)	// unbroken joints
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
      ot=MAXIM((R2*op),(R3*dpegt/dpeft));
      nfail=0;
	  d1jfmd[ijoint]=R0;
      for(integ=0;integ<3;integ++)
      { dpefs=d2elfs[integ][ielem];
		sp=R2*h*dpefs/dpepe;
		st=MAXIM((R2*sp),(R3*dpegs/dpefs));
		if(integ==0)
        { o=o1; s=s1;
        }
        else if(integ==2)
        { o=o2; s=s2;
        } 
        else
        { o=RP5*(o1+o2); s=RP5*(s1+s2);
        }
		sabs=ABS(s);
        if((o>op)&&(sabs>sp))
		{ z=SQRT(((o-op)/ot)*((o-op)/ot)+((sabs-sp)/st)*((sabs-sp)/st));
		}
		else if(o>op)
		{ z=ABS((o-op)/(ot-op));
		}
		else if(sabs>sp)
		{ z=(sabs-sp)/(st-sp);
		}
        else
        { z=R0;
        }
		if(z>R0)	  // seismic energy by crack propagation
		{ // kinetic energy
		  Ek=R0;
		  for(iside=0;iside<2;iside++)
		  { jelem=i2eljp[iside][ielem];
			for(i=0;i<3;i++)
			{ inode=i2elto[i][jelem];
			  Ek+=d1nmct[inode]*(d1nvcx[inode]*d1nvcx[inode]+d1nvcy[inode]*d1nvcy[inode]);
		  } }
		  if(i1elyi[ielem]==0)
		  { // start recording, store initial kinetic energy
			d1eley[ielem]=Ek;
			d1eles[ielem]=R0;
			i1elyi[ielem]=1;
		  }
		  else if(i1elyi[ielem]==1)
		  { // seismic energy recording
			if(Ek>=d1eley[ielem])
			{ d1eles[ielem]=MAXIM(d1eles[ielem],Ek-d1eley[ielem]);
			  d1elme[ielem]=R2/R3*(log10(d1eles[ielem])-11.8);
			}
			else
			{ d1eley[ielem]=Ek;
			  d1eles[ielem]=R0;
			  d1elme[ielem]=-BEPSILON;
			  i1elyi[ielem]=0;
		  } }
		}
        if(z>=R1)
        { nfail=nfail+1;
		  if(nfail>1&&i1elty[ielem]==0)
          { // failure mode
			if(sabs<=sp)		  // mode I
			{ d1jfmd[ijoint]=R1;
			}
			else if(o<=op)		  // mode II
			{ d1jfmd[ijoint]=R2;
			}
			else				  // mixed mode I-II
			{ d1jfmd[ijoint]=R1+ATAN(((sabs-sp)/(st-sp))/((o-op)/(ot-op)));
			}
			// stop seismic recording -- report acoustic emission
			i1elyi[ielem]=2;
			// new crack
			i1elty[ielem]=3;
			(*d0iedi)=BEPSILON;
			// fluid begins to fill
			d1jfet[ijoint]=dctime;
          }
          z=R1;
        }
        z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*
		  exp(z*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))
		  *(dpefa*(R1-z)+dpefb*pow((R1-z),dpefc));
		/* normal stress */
		if(o<R0)
        { sigma=R2*o*dpeft/op;
        }
        else if(o>op)
        { sigma=dpeft*z; nsoft=nsoft+1;
        }
        else
        { sigma=(R2*o/op-(o/op)*(o/op))*dpeft;
        }
		if(sigma>dpeft)
		{ dpefs=-dpeft*dpicf+dpcoh;
		}
		else
		{ dpefs=-sigma*dpicf+dpcoh;
		}
		d2elfs[integ][ielem]=dpefs;
		/* shear stress */
        if((sigma>R0)&&(sabs>sp))
		{ tau=z*dpefs;
        }
        else if(sigma>R0)
        { tau=(R2*(sabs/sp)-(sabs/sp)*(sabs/sp))*dpefs;
        }
        else if(sabs>sp)
        { tau=z*(dpefs+dpefr*sigma)-dpefr*sigma;
        }
        else
        { tau=(R2*(sabs/sp)-(sabs/sp)*(sabs/sp))*dpefs;
        }
		if(s<R0) tau=-tau;
        if((integ==0)||(integ==2))
		{ area=h/R6;
		}
        else
		{ area=h/R3;
        }
        if(integ==0)  /* nodal forces */
        { area=h/R6;  /* area=h/6.0;  */
          d1nfcx[i0]=d1nfcx[i0]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i3]=d1nfcx[i3]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i0]=d1nfcy[i0]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i3]=d1nfcy[i3]+area*(tau*e1y+sigma*e1x);
        }
        else if(integ==1)
        { area=h/R3;  /* area=h/3.0;  */
          d1nfcx[i0]=d1nfcx[i0]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i3]=d1nfcx[i3]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i0]=d1nfcy[i0]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i3]=d1nfcy[i3]+area*(tau*e1y+sigma*e1x);
          d1nfcx[i1]=d1nfcx[i1]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i2]=d1nfcx[i2]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i1]=d1nfcy[i1]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i2]=d1nfcy[i2]+area*(tau*e1y+sigma*e1x);
        }
        else
        { area=h/R6; /* area=h/6.0;  */
          d1nfcx[i1]=d1nfcx[i1]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i2]=d1nfcx[i2]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i1]=d1nfcy[i1]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i2]=d1nfcy[i2]+area*(tau*e1y+sigma*e1x);
	} } }
	else if(i1elty[ielem]>1)	// pre-existing & new fracture joints
	{ /* joint constitutive model */
	  // meso opening and shear displacement
	  i0=i2elto[0][ielem];
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
	  o=RP5*(o1+o2); u2=RP5*(s1+s2);

	  // fluid pressure
	  if(dcrmpt>R0)
	  { dprff=(dctime-d1jfet[ijoint])/dcrmpt;
	  }
	  else
	  { dprff=R1;
	  }
	  if(dprff>R1) dprff=R1;
	  d1jfpr[ijoint]=dprff*dpepf;

	  // extract normal and shear stress
	  for(iside=0;iside<2;iside++)
	  { nx=d1nccy[i2elto[iside*2][ielem]]-d1nccy[i2elto[iside*2+1][ielem]];
		ny=d1nccx[i2elto[iside*2+1][ielem]]-d1nccx[i2elto[iside*2][ielem]];
		h=SQRT(nx*nx+ny*ny); nx=nx/(h+small); ny=ny/(h+small);
		Sxx=d2tcs[0][i2eljp[iside][ielem]];
		Sxy=d2tcs[1][i2eljp[iside][ielem]];
		Syy=d2tcs[2][i2eljp[iside][ielem]];
		nstr[iside]=Sxx*nx*nx+2*Sxy*nx*ny+Syy*ny*ny;
		sstr[iside]=(Syy-Sxx)*nx*ny+Sxy*(nx*nx-ny*ny);
	  }
	  sigma2=-RP5*(nstr[0]+nstr[1])-d1jfpr[ijoint];	// (+) compressive

	  if(o<R0)			  /* closed fracture */
	  { if(sigma2<=R0) sigma2=EPSILON;
		a0 =d1japi[ijoint];
		ar =d1japr[ijoint];
		kni=d1jkni[ijoint];
		u1 =d1jsdc[ijoint];
		up =d1jsdp[ijoint];
		JRC=d1jjrc[ijoint];
		JCS=d1jjcs[ijoint];
		as1=d1jdlc[ijoint];
		an1=d1japc[ijoint];
		sigma1=d1jnst[ijoint];
		tau1=d1jsst[ijoint];
		blklen=d1jefl[ijoint];
		vm=a0-ar;									  // maximum closure
		PHIR=atan(dpefr)*RAD2DEG;					  // residual friction angle (deg)
		K=log10(JCS*1e6/sigma2);					  // ratio between JCS and stress level
		du=u2-u1;									  // increment of shear displacement
		dsigma=sigma2-sigma1;						  // increment of normal stress

		// peak shear displacement
		up=0.0077*pow(blklen,0.45)*pow(sigma2/JCS/1e6,0.34)*cos(JRC*K*DEG2RAD);
		// up=blklen/500.0*pow(JRC/blklen,0.33);
		if(up<EPSILON) up=EPSILON;
		uratio=DABS(u2/up);
		
		// mobilised JRC
		M=JRC/12/K+0.70;
		if(uratio<R1)
		{ // pre-peak stage
		  Yfd2JRCM(JRC,K,PHIR,uratio,&JRCm);
		  dm=atan((R4*uratio-R1)*tan(R1/R3*JRC*K*DEG2RAD))*RAD2DEG;
		}
		else
		{ // post-peak stage
		  JRCm=JRC*pow(R1/uratio,0.381);
		  dm=R1/M*JRCm*K;
		}
		Yfd2JRCM(JRC,K,PHIR,uratio,&JRCm);
		dm=R1/M*JRCm*K;

		// roughness induced friction angle (deg)
		PHII=JRCm*K;

		// normal stiffness
		knn=pow(sigma2+kni*vm,2)/kni/pow(vm,2);

		// shear stress and stiffness
		if(JRC*K+PHIR>85.0)
		{ taup=sigma2*tan(85.0*DEG2RAD);
		}
		else
		{ taup=sigma2*tan((JRC*K+PHIR)*DEG2RAD);
		}
		tau2=sigma2*tan((JRCm*K+PHIR)*DEG2RAD);
		ktt=(tau2-tau1)/DABS(du);

		// shear-normal stiffness
		knt=knn*tan(dm*DEG2RAD);

		// dilatancy
		if(du*u1>=R0)
		{ as2=as1+DABS(du)*tan(dm*DEG2RAD);
		}
		else if(u1*u2>=R0)
		{ as2=as1-DABS(du)*tan(dm*DEG2RAD);
		}
		else
		{ as2=DABS(u2)*tan(dm*DEG2RAD);
		}

		// normal displacement
		if(du*u1>=R0)
		{ dv=(dsigma-knt*DABS(du))/knn;
		}
		else if(u1*u2>=R0)
		{ dv=(dsigma+knt*DABS(du))/knn;
		}
		else
		{ dv=(dsigma+knt*DABS(u1)-knt*DABS(u2))/knn;
		}
		an2=an1-dv;
		if(an2<ar) an2=ar;

		// hydraulic aperture
		if(uratio<0.75)
		{ ah=pow(an2*1e6,R2)/pow(JRC,2.5)/1e6;
		}
		else if(uratio>1.0)
		{ ah=pow(an2*1e6,0.5)*JRCm/1e6;
		}
		else
		{ ah=((1-(1-uratio)/0.25)*pow(an2*1e6,0.5)*JRCm+(1-uratio)/0.25*pow(an2*1e6,R2)/pow(JRC,2.5))/1e6;
		}
		if(ah>an2) ah=an2;

		// update current joint properties
		d1jphi[ijoint]=PHII;
		d1jknc[ijoint]=(sigma2+kni*vm)/vm;
		d1jksc[ijoint]=ktt;
		d1jdlc[ijoint]=as2;
		d1japc[ijoint]=an2;
		d1jsdc[ijoint]=u2;
		d1jnst[ijoint]=sigma2;
		d1jsst[ijoint]=tau2;
		d1japh[ijoint]=ah;

		// seismic energy by shear reactivation
		if(u2>=up)
		{ Ek=R0;
		  for(iside=0;iside<2;iside++)
		  { jelem=i2eljp[iside][ielem];
			for(i=0;i<3;i++)
			{ inode=i2elto[i][jelem];
			  Ek+=d1nmct[inode]*(d1nvcx[inode]*d1nvcx[inode]+d1nvcy[inode]*d1nvcy[inode]);
		  } }
		}
		if(i1elyi[ielem]==0)
		{ // start recording, store initial kinetic energy
		  d1eley[ielem]=Ek;
		  d1eles[ielem]=R0;
		  i1elyi[ielem]=3;
		}
		else if(i1elyi[ielem]==3)
		{ // seismic energy recording
		  if(Ek>=d1eley[ielem])
		  { d1eles[ielem]=MAXIM(d1eles[ielem],Ek-d1eley[ielem]);
			d1elme[ielem]=R2/R3*(log10(d1eles[ielem])-11.8);
		} }
	  }
	  else				/* open fracture */
	  { JRC=d1jjrc[ijoint];
		ah=pow(d1japi[ijoint]*1e6,R2)/pow(JRC,2.5)/1e6;
		d1jknc[ijoint]=d1jkni[ijoint];
		d1jksc[ijoint]=R0;
		d1japc[ijoint]=o+d1japi[ijoint];
		d1japh[ijoint]=o+ah;
		d1jsdc[ijoint]=u2;
		d1jnst[ijoint]=R0;
		d1jsst[ijoint]=R0;
	  }
} } }

/*********************PUBLIC*************************************/
void Yfd(  yde, ydj, ydn, ydp, ydi, ydc    /***  nodal forces  ***/
        )
  YDE yde; YDJ ydj;  YDN ydn; YDP ydp; YDI ydi; YDC ydc;
{ INT iprop,inopo,ielem,integ;
  
  /* initialise shear strength */
  if(yde->d2elfs==DBL2NULL)
  { yde->d2elfs=TalDBL2(yde->nelem,3);
	for(ielem=0;ielem<yde->nelem;ielem++)
    { for(integ=0;integ<3;integ++)
	  { yde->d2elfs[integ][ielem]=ydp->d1pcoh[yde->i1elpr[ielem]];
	} }
  }
  /* initialise acoustic emission parameters */
  if(yde->i1elyi==INT1NULL)
  { yde->i1elyi=TalINT1(yde->melem);
	yde->d1eley=TalDBL1(yde->melem);
	yde->d1eles=TalDBL1(yde->melem);
	yde->d1elme=TalDBL1(yde->melem);
	for(ielem=0;ielem<yde->nelem;ielem++)
	{ yde->i1elyi[ielem]=0;
	  yde->d1eley[ielem]=R0;
	  yde->d1eles[ielem]=R0;
	  yde->d1elme[ielem]=-BEPSILON;
	}
  }
  /* zero nodal forces and masses */
  for(inopo=0;inopo<ydn->nnopo;inopo++) 
  { ydn->d1nmct[inopo]=R0;
    if(ydn->nnodim>0)
    { ydn->d2nfc[0][inopo]=R0;
      ydn->d2nfcon[0][inopo]=R0;
    }
    if(ydn->nnodim>1)
    { ydn->d2nfc[1][inopo]=R0;
      ydn->d2nfcon[1][inopo]=R0;
    }
    if(ydn->nnodim>2)
    { ydn->d2nfc[2][inopo]=R0;
      ydn->d2nfcon[2][inopo]=R0;
    }
  }
  for(iprop=0;iprop<ydp->nprop;iprop++)
  { if((ydp->i1ptyp[iprop])==(YTE2TRIELS))
    { Yfd2TRIELS(   /* small strain elastic triangle  */
	  ydc->dctime,ydc->dcgrst,ydc->dcrmpt,
      yde->nelem, ydj->njoint,
      iprop,
      ydp->d1peks[iprop],ydp->d1pela[iprop],ydp->d1pemu[iprop],
	  ydp->d1pero[iprop],ydp->d1psem[iprop],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nfc[0],
      ydn->d2nfc[1],ydn->d1nmct  ,ydn->d2nvc[0],ydn->d2nvc[1],ydp->d1pnaf  ,
      ydp->d1pnap  ,ydp->d1pnat  ,yde->i1elpr  ,yde->i2eljp  ,yde->i1elty  ,
      ydn->i1nopr  ,yde->i2elto  ,yde->d2tcs   ,ydj->d1jfpr
      );
    }
    else if((ydp->i1ptyp[iprop])==(YTE2TRISOF))
    { Yfd2TRISOF(   /* small strain softening triangle  */
      yde->nelem,
      iprop,
      ydp->d1peks[iprop],ydp->d1pela[iprop],
      ydp->d1pemu[iprop],ydp->d1pero[iprop],ydp->d1psem[iprop],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nfc[0],
      ydn->d2nfc[1],ydn->d1nmct,ydn->d2nvc[0],ydn->d2nvc[1],
      yde->d2elst[ydp->i1psde[iprop]],
      yde->i1elpr,yde->i2elto
      );
    }
    else if((ydp->i1ptyp[iprop])==(YTE2JOINTS))
    { Yfd2JOINTS(  /* joint element  */
	  ydc->dctime,ydc->dcgrst,ydc->dcrmpt,
      yde->nelem,
	  iprop,
      ydp->d1peft[iprop],ydp->d1pegt[iprop],ydp->d1pegs[iprop],
	  ydp->d1peks[iprop],ydp->d1pepe[iprop],ydp->d1picf[iprop],
	  ydp->d1pcoh[iprop],ydp->d1pefr[iprop],ydp->d1pepf[iprop],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nfc[0],ydn->d2nfc[1],
	  ydn->d2nvc[0],ydn->d2nvc[1],ydn->d1nmct,
	  yde->d2elst[ydp->i1psde[iprop]],yde->i1elpr,yde->i2elto,yde->i2eljp,
	  yde->d2elfs,yde->i1elty,yde->d2tcs ,yde->d1eley,
	  yde->d1eles,yde->d1elme,yde->i1elyi,&ydi->diedi,
	  ydj->njoint,ydj->d1jkni,ydj->d1jknc,ydj->d1jksc,ydj->d1jnst,
	  ydj->d1jsst,ydj->d1japi,ydj->d1japc,ydj->d1japh,ydj->d1japr,
	  ydj->d1jsdc,ydj->d1jdlc,ydj->d1jsdp,ydj->d1jefl,ydj->d1jjrc,
	  ydj->d1jjcs,ydj->d1jphi,ydj->d1jfmd,ydj->d1jfpr,ydj->d1jfet
      );
    }
    /* increment load factor */
    ydp->d1pnaf[iprop]=ydp->d1pnaf[iprop]+ydp->d1pnai[iprop];
  }
}
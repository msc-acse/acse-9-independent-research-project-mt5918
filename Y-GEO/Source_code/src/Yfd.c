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
            nelem,
            iprop,
            npnfact,mprop,nprop,
            d3pnfac,
            i1ptyp,
            dpeks,dpela, dpemu, dpero ,
            dpsem,
            dpeem,dpenu,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nfcx,
            d1nfcy,d1nmct,d1nvcx,d1nvcy,d1pnaf,
            d1pnap,d1pnat,
            i1elpr,i1nopr,i2elto,
            nohys, dohyp, dctime,
            d1ohys, d1ohyt, d1ohyx, d1ohyy,
            i1ohyt, npnset, d1elfr,
            i1usan, d1peex, d1peey, d1pemx, d1pemy, d1peg,
            iuseis, dcstxx, dcstxy, dcstyy,
            dcsyxx, dcsyxy, dcsyyy, dcsrfy,
            i1pnfx, i1pnfy,
            i1pexc, i1nowe, iusehf,
            nsbar, d2elstr
            )
  INT    nelem;
  INT    iprop;
  INT   npnfact; INT    mprop; INT    nprop;
  INT npnset;
  DBL ***d3pnfac;
  DBL    dpeks; DBL   dpela; DBL    dpemu; DBL   dpero;
  DBL   dpsem; DBL   dpeem; DBL    dpenu;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nfcx;
  DBL *d1nfcy; DBL  *d1nmct; DBL *d1nvcx; DBL  *d1nvcy; DBL *d1pnaf;
  DBL *d1pnap; DBL  *d1pnat;
  INT *i1elpr; INT  *i1nopr; INT **i2elto;
  INT   nohys; DBL    dohyp; DBL   dctime;
  DBL *d1ohys; DBL  *d1ohyt; DBL  *d1ohyx; DBL *d1ohyy;
  INT *i1ohyt;
  DBL *d1elfr;
  
  INT i1usan;
  DBL d1peex; DBL d1peey; DBL d1pemx; DBL d1pemy; DBL d1peg;
  
  INT iuseis; DBL  dcstxx; DBL  dcstxy;  DBL  dcstyy;
  DBL  dcsyxx; DBL  dcsyxy;  DBL  dcsyyy; DBL dcsrfy; 
  INT *i1pnfx; INT *i1pnfy;
  INT *i1pexc; INT *i1nowe; INT iusehf;
  INT nsbar; DBL **d2elstr;
  
{ DBL nx,ny,voli,volc;
  DBL v0, v1, v2, rpx, rpy, r0x, r0y, r1x, r1y, r2x, r2y, stprev;
  DBL  V[3];
  DBL  B[2][2];     /* left Cauchy-Green strain tensor                        */
  DBL  D[2][2];     /* rate of deformation (stretching) tensor                */
  DBL  E[2][2];     /* strain tensor (small strains)                          */
  DBL  F[2][2];     /* deformation gradient in global base                    */
  DBL F0[2][2];     /* initial local base                                     */
  DBL FX[2][2];     /* current local base                                     */
  DBL F0inv[2][2];  /* global base in initial local base                      */
  DBL FXinv[2][2];  /* global base in current local base                      */
  DBL  L[2][2];     /* velocity gradient in global base                       */
  DBL LX[2][2];     /* vel. gradient in current local base = delta x/delta X  */
  DBL  T[2][2];     /* Cauchy stress                                          */
  INT ielem;
  INT i,j,k,in,jn,kn,jnopr,knopr;
  INT ihys;
  DBL **d2fact;
  DBL **d2time;
  
  DBL Tinsitu[2][2]; /* element in-situ stress */
  DBL cp; /* element p-wave velocity */
  DBL cs; /* element s-wave velocity */
  DBL mu;
  DBL lambda;
  
  /* Elastic constants for plane strain transverse isotropy */
  DBL d1peex_pstrain; 
  DBL d1peey_pstrain;
  DBL d1pemx_pstrain;
  DBL d1pemy_pstrain;
  
  DBL Delta;
  
  if(d1pnaf == DBL1NULL)
  { d1pnaf = TalDBL1(mprop);
  }
  for(i=0; i<nprop; i++)
  { d1pnaf[i]=R1;
  }

  if(d3pnfac != DBL3NULL)
  { if(d3pnfac[0][0][0] != -R1)
    { d2time = d3pnfac[0];
      d2fact = d3pnfac[1];
      for(j=0; j<npnset; j++)
      { for(i=1; i<npnfact; i++)
        { if((dctime>=d2time[j][i-1])&&(dctime<=d2time[j][i]))
          { d1pnaf[j]=d2fact[j][i-1]-((d2fact[j][i-1]-d2fact[j][i])*
                      ((dctime-d2time[j][i-1])/(d2time[j][i]-d2time[j][i-1])));
  } } } } }
  
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
      if(i1usan==1) /* Apply transversely isotropic elastic constitutive law (plane stress)*/
      { T[0][0] = (d1peex / (1 - d1pemx * d1pemy)) * (E[0][0] + d1pemy * E[1][1]) + dpeks * D[0][0];
        T[1][1] = (d1peey / (1 - d1pemx * d1pemy)) * (E[1][1] + d1pemx * E[0][0]) + dpeks * D[1][1];
        T[0][1] =  (2 * d1peg) * E[0][1] + dpeks * D[0][1];
        T[1][0] =  (2 * d1peg) * E[1][0] + dpeks * D[1][0];
      }
      else if(i1usan==2) /* Apply transversely isotropic elastic constitutive law (plane strain)*/
      { d1peex_pstrain = d1peex / (1 - d1pemx * d1pemy);
        d1peey_pstrain = d1peey / (1 - d1pemy * d1pemx);
	d1pemx_pstrain = (d1pemx + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	d1pemy_pstrain = (d1pemy + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	
	T[0][0] = (d1peex_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[0][0] + d1pemy_pstrain * E[1][1]) + dpeks * D[0][0];
        T[1][1] = (d1peey_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[1][1] + d1pemx_pstrain * E[0][0]) + dpeks * D[1][1];
        T[0][1] =  (2 * d1peg) * E[0][1] + dpeks * D[0][1];
        T[1][0] =  (2 * d1peg) * E[1][0] + dpeks * D[1][0];
      }
      else
      { if (i1ptyp == YTE2TRIELS)
        {
          for(i=0;i<2;i++)    
          {  
              for(j=0;j<2;j++)
              { T[i][j]=R2*dpemu*E[i][j]*(voli/volc)+dpeks*D[i][j];
              }
              T[i][i]=T[i][i]+dpela*(volc/voli-voli/volc);
          }
        }
        /* Plane Stress formulation based on E,nu */
        else if (i1ptyp == YTE2PLANESTRESS)
        {        
          T[0][0] = (dpeem / (1 - dpenu*dpenu)) * (E[0][0] + dpenu * E[1][1]) + dpeks * D[0][0];
          T[1][1] = (dpeem / (1 - dpenu*dpenu)) * (E[1][1] + dpenu * E[0][0]) + dpeks * D[1][1];
          T[0][1] = (dpeem / (1 + dpenu)) * E[0][1] + dpeks * D[0][1];
          T[1][0] = (dpeem / (1 + dpenu)) * E[1][0] + dpeks * D[1][0];
        }
        /* Plane Strain formulation based on E,nu */
        else if (i1ptyp == YTE2PLANESTRAIN)
        {
          T[0][0] = dpeem/((1+dpenu)*(1-2*dpenu)) * ( (1-dpenu)*E[0][0] +  dpenu*E[1][1]) + dpeks * D[0][0];
          T[1][1] = dpeem/((1+dpenu)*(1-2*dpenu)) * ( (1-dpenu)*E[1][1] +  dpenu*E[0][0]) + dpeks * D[1][1];
          T[0][1] = (dpeem / (1+dpenu)) * E[0][1] + dpeks * D[0][1];
          T[1][0] = (dpeem / (1+dpenu)) * E[1][0] + dpeks * D[1][0];
        }
      }
      
      /* output history states */
      r0x = d1ncix[(i2elto[0][ielem])];
      r0y = d1nciy[(i2elto[0][ielem])];
      r1x = d1ncix[(i2elto[1][ielem])];
      r1y = d1nciy[(i2elto[1][ielem])];
      r2x = d1ncix[(i2elto[2][ielem])];
      r2y = d1nciy[(i2elto[2][ielem])];
      for(ihys=0; ihys<nohys; ihys++)
      { 
	if(i1ohyt[ihys]==(YFLEE)) /* Strain energy */
	{
	  d1ohyt[ihys] = dctime;    /* output history time  */
          if(i1usan==1) /* Transversely isotropic elasticity */
	  { d1ohys[ihys] += 0.5 * ((E[0][0] * ((d1peex / (1 - d1pemx * d1pemy)) * (E[0][0] + d1pemy * E[1][1]))) +      /* output history state */
                                   (E[0][1] * ((1 / (2 * d1peg)) * E[0][1])) +
                                   (E[1][0] * ((1 / (2 * d1peg)) * E[0][1])) +
                                   (E[1][1] * ((d1peey / (1 - d1pemx * d1pemy)) * (E[1][1] + d1pemx * E[0][0])))) * (volc/2);    
	  }
	  else if(i1usan==2) /* Transversely isotropic elasticity (plane strain) */
	  { d1peex_pstrain = d1peex / (1 - d1pemx * d1pemy);
            d1peey_pstrain = d1peey / (1 - d1pemy * d1pemx);
	    d1pemx_pstrain = (d1pemx + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	    d1pemy_pstrain = (d1pemy + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	    d1ohys[ihys] += 0.5 * ((E[0][0] * ((d1peex_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[0][0] + d1pemy * E[1][1]))) +      /* output history state */
                                   (E[0][1] * ((1 / (2 * d1peg)) * E[0][1])) +
                                   (E[1][0] * ((1 / (2 * d1peg)) * E[0][1])) +
                                   (E[1][1] * ((d1peey_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[1][1] + d1pemx_pstrain * E[0][0])))) * (volc/2);    
	  }
	  else /* Isotropic elasticity */
	  { d1ohys[ihys] += 0.5 * ((E[0][0] * (2.0*dpemu*E[0][0]*(voli/volc) + dpela*(volc/voli-voli/volc))) +      /* output history state */
                                 (E[0][1] * (2.0*dpemu*E[0][1]*(voli/volc))) +
                                 (E[1][0] * (2.0*dpemu*E[1][0]*(voli/volc))) +
                                 (E[1][1] * (2.0*dpemu*E[1][1]*(voli/volc) + dpela*(volc/voli-voli/volc)))) * (volc/2);    
	  }
	}
	rpx=d1ohyx[ihys];  /* x coordinate of point P */
        rpy=d1ohyy[ihys];  /* y coordinate of point P */
        V2DCro(v0,(r1x-r0x),(r1y-r0y),(rpx-r0x),(rpy-r0y));
        V2DCro(v1,(r2x-r1x),(r2y-r1y),(rpx-r1x),(rpy-r1y));
        V2DCro(v2,(r0x-r2x),(r0y-r2y),(rpx-r2x),(rpy-r2y));

        if((v0>R0)&&(v1>R0)&&(v2>R0))    /* if point is inside the triangle */
        { if(i1ohyt[ihys]==(YFLDSXX))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            if((ABS(T[0][0]-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;    /* output history time  */
              d1ohys[ihys] = T[0][0];    /* output history state */
          } }
          else if(i1ohyt[ihys]==(YFLDSXY))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            if((ABS(T[0][1]-stprev))>=dohyp)  /* if((ABS(R1-T[0][1]/stprev))>dohyp) */
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = T[0][1];
          } }
          else if(i1ohyt[ihys]==(YFLDSYY))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            if((ABS(T[1][1]-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = T[1][1];
          } }
          else if(i1ohyt[ihys]==(YFLDSZZ))
          { d1ohyt[ihys] = dctime;
            d1ohys[ihys] = R0;
          }
          else if(i1ohyt[ihys]==(YFLDSZX))
          { d1ohyt[ihys] = dctime;
            d1ohys[ihys] = R0;
          }
          else if(i1ohyt[ihys]==(YFLDSZY))
          { d1ohyt[ihys] = dctime;
            d1ohys[ihys] = R0;
          }
          else if(i1ohyt[ihys]==(YFLDVEL))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            for(i=0; i<3; i++)
            { V[i] = SQRT((d1nvcx[(i2elto[i][ielem])]*d1nvcx[(i2elto[i][ielem])])
                         +(d1nvcy[(i2elto[i][ielem])]*d1nvcy[(i2elto[i][ielem])]));
            }
            if((ABS(((V[0]+V[1]+V[2])/R3)-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = (V[0]+V[1]+V[2])/R3;  /* average velocity of element */
          } }
          else if(i1ohyt[ihys]==(YFLDVEX))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            for(i=0; i<3; i++)
            { V[i] = d1nvcx[(i2elto[i][ielem])];
            }
            if((ABS(((V[0]+V[1]+V[2])/R3)-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = (V[0]+V[1]+V[2])/R3;  /* average velocity x of element */
          } }
          else if(i1ohyt[ihys]==(YFLDVEY))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            for(i=0; i<3; i++)
            { V[i] = d1nvcy[(i2elto[i][ielem])];
            }
            if((ABS(((V[0]+V[1]+V[2])/R3)-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = (V[0]+V[1]+V[2])/R3;  /* average velocity y of element */
            } 
          } 
      } }
       
      /* Store element stress tensor only if rebars are used */ 
      if(nsbar>0)
      { d2elstr[0][ielem]=T[0][0];
        d2elstr[1][ielem]=T[0][1];
        d2elstr[2][ielem]=T[1][0]; 
        d2elstr[3][ielem]=T[1][1];
      }

      /* Nodal Forces */
      for(i=0;i<3;i++)
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
	
	if (iuseis==1) /* Apply element in-situ stress */
        { Tinsitu[0][0]=-dcstxx-dcsyxx*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
          Tinsitu[0][1]=-dcstxy-dcsyxy*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
          Tinsitu[1][0]=-dcstxy-dcsyxy*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
          Tinsitu[1][1]=-dcstyy-dcsyyy*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
	  /* Nodal Forces due to in-situ stress */
	  d1nfcx[in]=d1nfcx[in]+(Tinsitu[0][0]*nx+Tinsitu[0][1]*ny)/R2;
	  d1nfcy[in]=d1nfcy[in]+(Tinsitu[1][0]*nx+Tinsitu[1][1]*ny)/R2;
	}
	
       /* Nodal Forces due to edge force*/
        jnopr=i1nopr[jn];
        knopr=i1nopr[kn];
        if( ((DABS(d1pnap[jnopr]))>EPSILON)&&
            ((DABS(d1pnap[knopr]))>EPSILON) )
        { d1nfcx[jn]=d1nfcx[jn]-
          d1pnap[jnopr]*d1pnaf[jnopr]*nx/R3-
          d1pnap[knopr]*d1pnaf[knopr]*nx/R6;
          d1nfcy[jn]=d1nfcy[jn]-
          d1pnap[jnopr]*d1pnaf[jnopr]*ny/R3-
          d1pnap[knopr]*d1pnaf[knopr]*ny/R6;
          d1nfcx[kn]=d1nfcx[kn]-
          d1pnap[jnopr]*d1pnaf[jnopr]*nx/R6-
          d1pnap[knopr]*d1pnaf[knopr]*nx/R3;
          d1nfcy[kn]=d1nfcy[kn]-
          d1pnap[jnopr]*d1pnaf[jnopr]*ny/R6-
          d1pnap[knopr]*d1pnaf[knopr]*ny/R3;
        }
        if( ((DABS(d1pnat[jnopr]))>EPSILON)&&
            ((DABS(d1pnat[knopr]))>EPSILON) )
        { d1nfcx[jn]=d1nfcx[jn]-
          d1pnat[jnopr]*d1pnaf[jnopr]*ny/R3-
          d1pnat[knopr]*d1pnaf[knopr]*ny/R6;
          d1nfcy[jn]=d1nfcy[jn]+
          d1pnat[jnopr]*d1pnaf[jnopr]*nx/R3+
          d1pnat[knopr]*d1pnaf[knopr]*nx/R6;
          d1nfcx[kn]=d1nfcx[kn]-
          d1pnat[jnopr]*d1pnaf[jnopr]*ny/R6-
          d1pnat[knopr]*d1pnaf[knopr]*ny/R3;
          d1nfcy[kn]=d1nfcy[kn]+
          d1pnat[jnopr]*d1pnaf[jnopr]*nx/R6+
          d1pnat[knopr]*d1pnaf[knopr]*nx/R3;
        }
        /* Nodal forces due to absorbing boundary condition */
        if (i1ptyp == YTE2TRIELS)
        { cp=SQRT((2*dpemu+dpela)/dpero);
          cs=SQRT(dpemu/dpero);
        }
        else if (i1ptyp == YTE2PLANESTRESS)
        { lambda=(dpeem*dpenu)/((1+dpenu)*(1-2*dpenu));
          mu=(dpeem)/(2*(1+dpenu));
          cp=SQRT((2*mu+lambda)/dpero);
          cs=SQRT(mu/dpero);
        }  
        else if (i1ptyp == YTE2PLANESTRAIN)
        { lambda=(dpeem*dpenu)/((1+dpenu)*(1-2*dpenu));
          mu=(dpeem)/(2*(1+dpenu));
          cp=SQRT((2*mu+lambda)/dpero);
          cs=SQRT(mu/dpero);
        }
        if ((i1pnfx[jnopr]==4)&&(i1pnfx[knopr]==4))
        { cp=SQRT((2*dpemu+dpela)/dpero);
          cs=SQRT(dpemu/dpero);
          d1nfcx[jn]=d1nfcx[jn]-dpero*cp*d1nvcx[jn]*ABS(nx)/R3-
                     dpero*cp*d1nvcx[kn]*ABS(nx)/R6-
                     dpero*cs*d1nvcx[jn]*ABS(ny)/R3-
                     dpero*cs*d1nvcx[kn]*ABS(ny)/R6;
          d1nfcx[kn]=d1nfcx[kn]-dpero*cp*d1nvcx[jn]*ABS(nx)/R6-
                     dpero*cp*d1nvcx[kn]*ABS(nx)/R3-
                     dpero*cs*d1nvcx[jn]*ABS(ny)/R6-
                     dpero*cs*d1nvcx[kn]*ABS(ny)/R3;
        }
        if ((i1pnfy[jnopr]==4)&&(i1pnfy[knopr]==4))
        { cp=SQRT((2*dpemu+dpela)/dpero);
          cs=SQRT(dpemu/dpero);
          d1nfcy[jn]=d1nfcy[jn]-dpero*cp*d1nvcy[jn]*ABS(ny)/R3-
                     dpero*cp*d1nvcy[kn]*ABS(ny)/R6-
                     dpero*cs*d1nvcy[jn]*ABS(nx)/R3-
                     dpero*cs*d1nvcy[kn]*ABS(nx)/R6;
          d1nfcy[kn]=d1nfcy[kn]-dpero*cp*d1nvcy[jn]*ABS(ny)/R6-
          dpero*cp*d1nvcy[kn]*ABS(ny)/R3-
          dpero*cs*d1nvcy[jn]*ABS(nx)/R6-
          dpero*cs*d1nvcy[kn]*ABS(nx)/R3;
        }
      }
      //! If element is "excavated" set nodal boundary condition to v_x = 0 and v_y = 0 (hardcoded bc ID)
      //! and translate by a constant vector (hardcoded, 200)
      if(i1pexc[i1elpr[ielem]]==1)
      { i1nopr[i2elto[0][ielem]] = 1; /* 1 is a hard-coded value */
        i1nopr[i2elto[1][ielem]] = 1; /* 1 is a hard-coded value */
	i1nopr[i2elto[2][ielem]] = 1; /* 1 is a hard-coded value */
	Delta = 200;
	d1nccx[(i2elto[0][ielem])]=d1ncix[(i2elto[0][ielem])]+Delta;
	d1nccx[(i2elto[1][ielem])]=d1ncix[(i2elto[1][ielem])]+Delta;
	d1nccx[(i2elto[2][ielem])]=d1ncix[(i2elto[2][ielem])]+Delta;
	d1nccy[(i2elto[0][ielem])]=d1nciy[(i2elto[0][ielem])]+Delta;
	d1nccy[(i2elto[1][ielem])]=d1nciy[(i2elto[1][ielem])]+Delta;
	d1nccy[(i2elto[2][ielem])]=d1nciy[(i2elto[2][ielem])]+Delta;
	if(iusehf==1)
	{ //! For hydrofrac: set nodes to be non-wettable
	  i1nowe[i2elto[0][ielem]]=3;
	  i1nowe[i2elto[1][ielem]]=3;
	  i1nowe[i2elto[2][ielem]]=3;
	}
      }
    }
  }
  FREE(d1pnaf);
}

static void Yfd2TRISOF(  /* small strain softening triangle  */
            nelem,
            iprop,
            dpeks,dpela ,dpemu ,dpero ,
            dpsem,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nfcx,
            d1nfcy,d1nmct,d1nvcx,d1nvcy,d1sdel,
            i1elpr,i2elto,d1elfs
            )
  INT    nelem;
  INT    iprop;
  DBL    dpeks; DBL    dpela; DBL   dpemu; DBL   dpero;
  DBL   dpsem;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nfcx;
  DBL *d1nfcy; DBL  *d1nmct; DBL *d1nvcx; DBL  *d1nvcy; DBL *d1sdel;
  INT *i1elpr; INT **i2elto; DBL *d1elfs;
{ DBL nx,ny,voli,volc,emax,dmax,dmas;
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
  DBL dpefs;

  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { dpefs=d1elfs[ielem];
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

static void Yfd2JOINTS(  /* joint element */
            nelem, iprop,
            dpeft, dpegf, dpegs,
            dpeco, dpefr, dpepe,
            d1nccx,d1nccy,d1nfcx,d1nfcy,d1nvcx,
            d1nvcy,
            i1elpr,i2elto,d1elfs,
            dctime, dcstec, nebrk,netbrk,i1ebrk,d2ecbrk,d2ecbrk_NEW,
            d1etbrk,d1elbrk,d1efe,
            nesft,netsft,i1esft,d2ecsft,d1etsft,i1esftf,d1ebrkf,
            d1eike,d1edke,d1nmct,d1etmke,
            i2elnext,i2eledge,
            i1pexc,
            dusaf,dpealp,
            dpecor,dpefrrd,dpeftr,dpegfr,dpegsr,
            i1nowe,i2noid, 
            iusefn,i1edfnf,ddfnft,ddfnco,ddfngf,ddfngs,iusehf,
            i1edft,d1etike,iusesm,dctwle,
            iusehy,d2eldmg
            )
  INT   nelem; INT   iprop;
  DBL   dpeft; DBL   dpegf; DBL   dpegs;
  DBL   dpeco; DBL   dpefr; DBL   dpepe;
  DBL *d1nccx; DBL *d1nccy; DBL *d1nfcx; DBL *d1nfcy;
  DBL *d1nvcx; DBL *d1nvcy;
  INT *i1elpr; INT **i2elto; DBL *d1elfs;
  DBL dctime; DBL dcstec;
  INT  *nebrk; INT  *netbrk; INT *i1ebrk;
  DBL  **d2ecbrk; DBL  **d2ecbrk_NEW; DBL *d1etbrk; DBL *d1elbrk; DBL *d1efe;
  INT  *nesft; INT  *netsft; INT *i1esft;
  DBL  **d2ecsft; DBL *d1etsft; INT *i1esftf; DBL *d1ebrkf;
  DBL *d1eike; DBL *d1edke; DBL *d1nmct; DBL *d1etmke;
  INT **i2elnext; INT **i2eledge;
  INT *i1pexc;
  DBL dusaf; DBL dpealp; DBL dpecor; DBL dpefrrd; DBL dpeftr; DBL dpegfr; DBL dpegsr;
  INT *i1nowe; INT **i2noid; 
  INT iusefn; INT *i1edfnf; DBL ddfnft; DBL ddfnco; DBL ddfngf; DBL ddfngs; INT iusehf;
  INT *i1edft; DBL *d1etike; INT iusesm; DBL dctwle;
  INT iusehy; DBL **d2eldmg;
  
{ DBL dpefa=0.63;
  DBL dpefb=1.8;
  DBL dpefc=6.0;
  DBL dpefm=0.0;
  DBL small,sabs,o,s,o1,o2,s1,s2,op,sp,ot,st,dmg,z,sigma,tau;
  DBL e1x,e1y,h,area;
  INT ielem,integ,i0,i1,i2,i3,nfail,el1,el2,el1edge,el2edge;
  INT nsoft;
  DBL dpefs;
  DBL joint_ke; /* joint kinetic energy (i.e. kinetic energy of the four joint nodes) */
  DBL delta_ke; /* difference between joint_ke and the joint kinetic energy calculated when it first yields */
  
  DBL beta;  /* joint orientation angle (0-180) */
  DBL gamma; /* relative angle between layering and joint element (0-90) */
  
  //! Assigning strength input values to "peak" values
  DBL dpeftp = dpeft;
  DBL dpegfp = dpegf;
  DBL dpegsp = dpegs;
  DBL C; /* exponent for power-law function */
  
  //DBL Fxi0;
  //DBL Fyi0;
  
  /*static FILE *out1=FILENULL;
  if(out1 == FILENULL)
  { out1=fopen("Yfd2JOINTS.txt", "a");
  }*/
  
  small=EPSILON; 
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    {// if((ielem%75)==0){i1elpr[ielem]=iprop-YIPROPMAX;}
      dpefs=d1elfs[ielem];
      i0=i2elto[0][ielem];
      i1=i2elto[1][ielem];
      i2=i2elto[2][ielem];
      i3=i2elto[3][ielem];
      
      //! Break random joints
      /*if((ncstep==1)&&((ielem%10)==0)&&(i1elpr[ielem]>=0)&&(i0!=i3))
      { i1elpr[ielem]=iprop-YIPROPMAX;
        i2noid[0][i0]=i1;
        i2noid[1][i1]=i0;
        i2noid[0][i2]=i3;
        i2noid[1][i3]=i2;
      }*/
      
      e1x=RP5*(d1nccx[i1]+d1nccx[i2]-d1nccx[i0]-d1nccx[i3]);
      e1y=RP5*(d1nccy[i1]+d1nccy[i2]-d1nccy[i0]-d1nccy[i3]);
      h=SQRT(e1x*e1x+e1y*e1y);
      e1x=e1x/(h+small);
      e1y=e1y/(h+small);
      s1=(d1nccy[i0]-d1nccy[i3])*e1y+(d1nccx[i0]-d1nccx[i3])*e1x;
      s2=(d1nccy[i1]-d1nccy[i2])*e1y+(d1nccx[i1]-d1nccx[i2])*e1x;
      o1=(d1nccy[i0]-d1nccy[i3])*e1x-(d1nccx[i0]-d1nccx[i3])*e1y;
      o2=(d1nccy[i1]-d1nccy[i2])*e1x-(d1nccx[i1]-d1nccx[i2])*e1y;
      //! Anisotropic fracture model
      if(dusaf>0.0)
      {
        //! Calculation of joint orientation (i.e. angle beta)
        //beta = atan( (((d1nccy[i1]+d1nccy[i2])/2)-((d1nccy[i0]+d1nccy[i3])/2))/(((d1nccx[i1]+d1nccx[i2])/2)-((d1nccx[i0]+d1nccx[i3])/2))) * 180/MYPI;
        beta = atan((d1nccy[i1]-d1nccy[i0])/(d1nccx[i1]-d1nccx[i0])) * 180/MYPI; // simplified formula
        if (beta < 0.0)
        { beta = beta + 180.0; }
      
        //! Calculation of relative angle between layering and joint element
        gamma = ABS(dpealp - beta);
        if(gamma > 90.0)
        { 
          gamma = 180.0 - gamma; 
        }
        
        //! Power-law variation with exponent C = IUSAF
        C = dusaf;
        dpeft = dpeftr + (dpeftp-dpeftr) * pow((gamma/90.0),C);
        dpegf = dpegfr + (dpegfp-dpegfr) * pow((gamma/90.0),C);
        dpegs = dpegsr + (dpegsp-dpegsr) * pow((gamma/90.0),C);
        // dpefs is updated below according to Mohr-Coulomb */        

      }
      
      op=R2*h*dpeft/dpepe;
      sp=R2*h*dpefs/dpepe;
      ot=MAXIM(EPSILON,(R3*dpegf/dpeft));
      st=MAXIM(EPSILON,(R3*dpegs/dpefs));
      //ot=MAXIM((R2*op),(R3*dpegf/dpeft));
      //st=MAXIM((R2*sp),(R3*dpegs/dpefs));
      
      //! Use "cohesive" DFN properties to calculate op, sp, ot, and st
      if((iusefn == 2) && (i1edfnf[ielem]==1))
      { //! If joint element belongs to DFN
        op=R2*h*ddfnft/dpepe;
        sp=R2*h*dpefs/dpepe;
        ot=MAXIM(EPSILON,(R3*ddfngf/ddfnft));
        st=MAXIM(EPSILON,(R3*ddfngs/dpefs));
        // dpefs is updated below according to Mohr-Coulomb 
      }
      
      nfail=0;
      nsoft=0;
      
      el1=i2elnext[0][ielem];     /* 1st element next to the joint */
      el2=i2elnext[1][ielem];     /* 2nd element next to the joint */
      el1edge=i2elnext[2][ielem]; /* edge number (0,1,2) of the 1st element (el1) */
      el2edge=i2elnext[3][ielem]; /* edge number (0,1,2) of the 2nd element (el2) */
      if((el1==-1)||(el2==-1)) continue; /* no need to do further computations if this is an external edge */
      
      //! Applying mixed DFN (i.e., DFN type 3) using the flag assigned to the joint element
      if(iusefn == 3)
      { if(i1edft[ielem] == 1) //! Broken-type DFN crack
      	{ i1elpr[ielem]=iprop-YIPROPMAX;
          d1ebrkf[ielem]=5.0;
          //! For hydrofrac
          i2noid[0][i0]=i1;
          i2noid[1][i1]=i0;
          i2noid[0][i2]=i3;
          i2noid[1][i3]=i2;
        }
        if(i1edft[ielem] == 2) //! Cohesive-type DFN crack
        { op=R2*h*ddfnft/dpepe;
          sp=R2*h*dpefs/dpepe;
          ot=MAXIM(EPSILON,(R3*ddfngf/ddfnft));
          st=MAXIM(EPSILON,(R3*ddfngs/dpefs));
          // dpefs is updated below according to Mohr-Coulomb
        }
      }
                  
      //! Performing excavation: set joint element state to broken if between at least one "excavated" element
      if((i1pexc[i1elpr[el1]]==1)||(i1pexc[i1elpr[el2]]==1))
      { i1elpr[ielem]=iprop-YIPROPMAX; 
        d1ebrkf[ielem]=4.0;
        if(iusehf==1)
	{ //! For hydrofrac
          i2noid[0][i0]=i1;
          i2noid[1][i1]=i0;
          i2noid[0][i2]=i3;
          i2noid[1][i3]=i2;
        }
      }
      else 
      {
      //Fxi0=0.0;
      //Fyi0=0.0;
      for(integ=0;integ<3;integ++)
      { if(integ==0)
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
        { dmg=SQRT(((o-op)/ot)*((o-op)/ot)+((sabs-sp)/st)*((sabs-sp)/st));
          nsoft=nsoft+1;
          if((nsoft>2)&&(i1esftf[ielem]==0))  
          { i1esft[*nesft]=ielem;
            d1etsft[*nesft]=dctime;
            d1etmke[ielem]=dctime;
            d2ecsft[0][*nesft]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
            d2ecsft[1][*nesft]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
            /* Kinetic energy of the joint as soon as it yields */ 
            d1eike[ielem]= 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                                d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                                d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                                d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
            d1etike[ielem] = dctime; //! initial time of KE monitoring window
            i1esftf[ielem]=1;
            (*nesft)++;
            (*netsft)++;
          }
        }
        else if(o>op)
        { dmg=(o-op)/ot;
          nsoft=nsoft+1;
          if((nsoft>2)&&(i1esftf[ielem]==0))  
          { 
            i1esft[*nesft]=ielem;
            d1etsft[*nesft]=dctime;
            d1etmke[ielem]=dctime;
            d2ecsft[0][*nesft]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
            d2ecsft[1][*nesft]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
            /* Kinetic energy of the joint as soon as it yields */ 
            d1eike[ielem]= 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                                d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                                d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                                d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
            d1etike[ielem] = dctime; //! initial time of KE monitoring window
            i1esftf[ielem]=2;
            (*nesft)++;
            (*netsft)++;
          }
        }
        else if(sabs>sp)
        { dmg=(sabs-sp)/st;
          nsoft=nsoft+1;
          if((nsoft>2)&&(i1esftf[ielem]==0))  
          { 
            i1esft[*nesft]=ielem;
            d1etsft[*nesft]=dctime;
            d1etmke[ielem]=dctime;
            d2ecsft[0][*nesft]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
            d2ecsft[1][*nesft]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
            /* Kinetic energy of the joint as soon as it yields */ 
            d1eike[ielem]= 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                                d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                                d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                                d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
            d1etike[ielem] = dctime; //! initial time of KE monitoring window
            i1esftf[ielem]=3;
            (*nesft)++;
            (*netsft)++;
          }
        }
        else
        { dmg=R0;
        }
        if(dmg>=R1) /* joint element broken */
        { nfail=nfail+1;
          if((nfail>1)&&(i1elpr[ielem]>=0))
          {
            i1elpr[ielem]=iprop-YIPROPMAX;
            i1ebrk[*nebrk]=ielem;
            d1etbrk[*nebrk]=dctime;
            //i1ebrkf[ielem]=i1esftf[ielem];
            if((o>=(op+ot))&&(sabs>=(sp+st))) // Mode 1 + mode 2 failure
            { d1ebrkf[ielem]=3.0; 
            }
            else if(o>=(op+ot)) // Mode 1 failure
            { d1ebrkf[ielem]=1.0; 
            }
            else if (sabs>=(sp+st)) // Mode 2 failure
            { d1ebrkf[ielem]=2.0; 
            }
            else 
            //{ d1ebrkf[ielem]=1.0+(sabs-sp)/(st-sp); } // Mode 1 + mode 2 failure (vectorial sum of s and o overcomes residual value)
            { d1ebrkf[ielem]=1.0+(sabs-sp)/(st); 
            } // Mode 1 + mode 2 failure (vectorial sum of s and o overcomes residual value)
            d2ecbrk[0][*nebrk]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
            d2ecbrk[1][*nebrk]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
            d2ecbrk_NEW[0][ielem]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
            d2ecbrk_NEW[1][ielem]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
            d1elbrk[*nebrk]=SQRT((d1nccx[i0]-d1nccx[i1])*(d1nccx[i0]-d1nccx[i1])+(d1nccy[i0]-d1nccy[i1])*(d1nccy[i0]-d1nccy[i1]));
            (*nebrk)++;
            (*netbrk)++;
            //! For hydrofrac
             i2noid[0][i0]=i1;
             i2noid[1][i1]=i0;
             i2noid[0][i2]=i3;
             i2noid[1][i3]=i2;
          }
          dmg=R1;
        }
        
        /* Specify intact/broken element edges */
        if(i0!=i3 && i1!=i2)
        { if(dmg==R1)
          { i2eledge[el1edge][el1]=-1;  /* -1: broken */
            i2eledge[el2edge][el2]=-1;  
          }
          else
          { i2eledge[el1edge][el1]=1;   /* 1: intact  */ 
            i2eledge[el2edge][el2]=1;
        } }
        
        /* Calculation of stress multiplier (z) from damage coefficient (dmg) */
        
        if(iusehy==0) /* Loading curve = unloading curve (classic formulation of Y-code) */
        { z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(dmg*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-dmg)+dpefb*pow((R1-dmg),dpefc));
        }
        else /* Hysteretic model with linear unloading */
        { if((integ==0) && (dmg>=d2eldmg[0][ielem])) /* Loading for integration point 0 */
          { d2eldmg[0][ielem]=dmg;
            z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[0][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[0][ielem])+dpefb*pow((R1-d2eldmg[0][ielem]),dpefc));
          }
          else if((integ==1) && (dmg>=d2eldmg[1][ielem])) /* Loading for integration point 1 */
          { d2eldmg[1][ielem]=dmg;
            z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[1][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[1][ielem])+dpefb*pow((R1-d2eldmg[1][ielem]),dpefc));
          }
          else if((integ==2) && (dmg>=d2eldmg[2][ielem]))  /* Loading for integration point 2 */
          { d2eldmg[2][ielem]=dmg;
            z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[2][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[2][ielem])+dpefb*pow((R1-d2eldmg[2][ielem]),dpefc));
          }
          if((integ==0) && (dmg<d2eldmg[0][ielem])) /* Unloading for integration point 0 */
          {  z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[0][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[0][ielem])+dpefb*pow((R1-d2eldmg[0][ielem]),dpefc));
             z=((op+dmg*ot)/(op+d2eldmg[0][ielem]*ot))*z;
          }
          else if((integ==1) && (dmg<d2eldmg[1][ielem])) /* Unloading for integration point 1 */
          {  z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[1][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[1][ielem])+dpefb*pow((R1-d2eldmg[1][ielem]),dpefc));
             z=((op+dmg*ot)/(op+d2eldmg[1][ielem]*ot))*z;
          }
          else if((integ==2) && (dmg<d2eldmg[2][ielem])) /* Unloading for integration point 2 */
          {  z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[2][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[2][ielem])+dpefb*pow((R1-d2eldmg[2][ielem]),dpefc));
             z=((op+dmg*ot)/(op+d2eldmg[2][ielem]*ot))*z;
          }
        }  
  
        if(o<R0)           /* normal stress */
        { sigma=R2*o*dpeft/op; /* sigma=R0; */
        }
        else if(o>op)
        { sigma=dpeft*z; nsoft=nsoft+1;
        }
        else
        { sigma=(R2*o/op-(o/op)*(o/op))*z*dpeft;
        }
        /* take into account Mohr-Coulomb   */
        if(dpeco>R0)
        { if(sigma>R0)
          { 
            if(dusaf>0) //! Anisotropic fracture model
            { 
              C = dusaf;
              dpefs = dpecor + (dpeco-dpecor) * pow((gamma/90.0),C); 
            }
            else
            { 
              dpefs = dpeco; 
            }
            //! Use "cohesive" DFN properties if joint element belongs to DFN
            //if((iusefn == 2)&&(i1edfnf[ielem]==1))
            if(((iusefn == 2)&&(i1edfnf[ielem]==1))||((iusefn == 3) && (i1edft[ielem] == 2)))
            { dpefs = ddfnco; }
          }
          else
          { 
            if(dusaf>0) //! Anisotropic fracture model
            { 
              C = dusaf; 
              dpefs = dpecor + (dpeco-dpecor) * pow((gamma/90.0),C) - sigma * (dpefrrd + (dpefr - dpefrrd) * pow((gamma/90.0),C));               
            }
            else
            { 
              dpefs = dpeco-sigma*dpefr; 
            }
            //! Use "cohesive" DFN properties if joint element belongs to DFN
            //if((iusefn == 2)&&(i1edfnf[ielem]==1))
            if(((iusefn == 2)&&(i1edfnf[ielem]==1))||((iusefn == 3) && (i1edft[ielem] == 2)))
            { dpefs = ddfnco-sigma*dpefr;}
        } }
        if((sigma>R0)&&(sabs>sp))           /* shear stress */
        { tau=z*dpefs;
        }
        else if(sigma>R0)
        { tau=(R2*(sabs/sp)-(sabs/sp)*(sabs/sp))*z*dpefs;
        }
        else if(sabs>sp)
        { tau=z*dpefs-dpefm*sigma;
        }
        else
        { tau=(R2*(sabs/sp)-(sabs/sp)*(sabs/sp))*(z*dpefs-dpefm*sigma);
        }
        d1elfs[ielem]=dpefs;    /* update fs to Mohr-Coulomb    */
        if(s<R0)tau=-tau;
        if(integ==0)  /* nodal forces */
        { area=h/R6; /* area=h/6.0; */
          d1nfcx[i0]=d1nfcx[i0]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i3]=d1nfcx[i3]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i0]=d1nfcy[i0]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i3]=d1nfcy[i3]+area*(tau*e1y+sigma*e1x);
          /* Element fracture energy */
          d1efe[ielem]=d1efe[ielem]+dcstec*((-area*(tau*e1x-sigma*e1y))*d1nvcx[i0]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i0]+
                        (area*(tau*e1x-sigma*e1y))*d1nvcx[i3]+(area*(tau*e1y+sigma*e1x))*d1nvcy[i3]);
          
          /*if(ielem==3)
          { Fxi0=Fxi0-area*(tau*e1x-sigma*e1y);
            Fyi0=Fyi0-area*(tau*e1y+sigma*e1x); }*/
        }
        else if(integ==1)
        { area=h/R3;  /* area=h/3.0; */
          d1nfcx[i0]=d1nfcx[i0]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i3]=d1nfcx[i3]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i0]=d1nfcy[i0]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i3]=d1nfcy[i3]+area*(tau*e1y+sigma*e1x);
          d1nfcx[i1]=d1nfcx[i1]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i2]=d1nfcx[i2]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i1]=d1nfcy[i1]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i2]=d1nfcy[i2]+area*(tau*e1y+sigma*e1x);
          /* Element fracture energy */
          d1efe[ielem]=d1efe[ielem]+dcstec*((-area*(tau*e1x-sigma*e1y))*d1nvcx[i0]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i0]+
                       (-area*(tau*e1x-sigma*e1y))*d1nvcx[i1]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i1]+
                       (+area*(tau*e1x-sigma*e1y))*d1nvcx[i2]+(+area*(tau*e1y+sigma*e1x))*d1nvcy[i2]+
                       (+area*(tau*e1x-sigma*e1y))*d1nvcx[i3]+(+area*(tau*e1y+sigma*e1x))*d1nvcy[i3]);

          /*if(ielem==3)
          { Fxi0=Fxi0-area*(tau*e1x-sigma*e1y);  
            Fyi0=Fyi0-area*(tau*e1y+sigma*e1x); }*/
          //if(ielem==3)
          //{ fprintf(out1,"%.6f \t %.6f \t %.6f \n",dctime,sigma,tau); }
        }
        else
        { area=h/R6; /* area=h/6.0; */
          d1nfcx[i1]=d1nfcx[i1]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i2]=d1nfcx[i2]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i1]=d1nfcy[i1]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i2]=d1nfcy[i2]+area*(tau*e1y+sigma*e1x);
          /* Element fracture energy */
          d1efe[ielem]=d1efe[ielem]+dcstec*((-area*(tau*e1x-sigma*e1y))*d1nvcx[i1]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i1]+
                       (+area*(tau*e1x-sigma*e1y))*d1nvcx[i2]+(+area*(tau*e1y+sigma*e1x))*d1nvcy[i2]);
        }
      }
      }
      /* If joint is yielded compute differential kinetic energy */
      if(i1esftf[ielem]>0)
      { joint_ke = 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                        d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                        d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                        d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
        delta_ke = joint_ke - d1eike[ielem];
        /* If the differential kinetic energy is greater than that from previous timestep update the joint differential kin energy */
        if(delta_ke > d1edke[ielem])
        { d1edke[ielem] = delta_ke; 
          d1etmke[ielem] = dctime;
        }
        if(iusesm == 1) //! Use maximum time window duration
        { if (dctime >= d1etike[ielem] + dctwle)
          { d1eike[ielem] = joint_ke;
            d1etike[ielem] = dctime; 
        } }
      }
      if(iusesm == 2) //! Record energy at the time of failure
      { if(d1ebrkf[ielem]>0)
        { joint_ke = 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                          d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                          d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                          d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
          d1edke[ielem] = joint_ke;
          d1etmke[ielem] = dctime;
      } }
    }
  }
}
/*********************PUBLIC*************************************/
void Yfd(   ydc,  yde,  ydn,  ydo,  ydpe, ydpn, ydpj, ydis, ydfn, ydhf, ydsm, ydsb    /***  nodal forces  ***/
        )
  YDC ydc; YDE yde; YDN ydn; YDO ydo; YDPE ydpe; YDPN ydpn; YDPJ ydpj; YDIS ydis; YDFN ydfn; YDHF ydhf; YDSM ydsm; YDSB ydsb;
{ INT iprop,inopo,jprop,i,j;
  INT ielem;
  static INT pmcstep=0; /* previous maximum number of time steps */
  INT ihys;
  
  INT s,k,r;
  DBL xi_0,yi_0,xi_1,yi_1;
  
  /* zero model strain energy */
  for(ihys=0;ihys<ydo->nohys;ihys++)
  { if(ydo->i1ohyt[ihys]==(YFLEE)) 
    { ydo->d1ohys[ihys] = 0.0;
    }
  }
  
    /* init. shear strength from joint database */
  if(ydc->ncstep==0 || ydc->ncstep==pmcstep)
  { pmcstep=ydc->mcstep;
    if(ydpj->npjset>0)
    { yde->d1elfs=TalDBL1(yde->melem);
      for(ielem=0;ielem<yde->melem;ielem++)
      { if(yde->i1elpr[ielem]>=ydpe->nprop)  /* joints   */
        { jprop=yde->i1elpr[ielem]-ydpe->nprop;
          if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
          { if(ydpj->d1pjfs[jprop]>R0)
            { yde->d1elfs[ielem]=ydpj->d1pjfs[jprop];
            }
            else
            { yde->d1elfs[ielem]=ydpj->d1pjco[jprop];  
        } } }
        else                                /* trians   */
        { iprop=yde->i1elpr[ielem];
	  // This part has to be commented out (still to figure out why :-)
	  /*if((ydpe->i1ptyp[iprop])==(YTE2TRISOF))
          { yde->d1elfs[ielem]=ydpj->d1pjfs[ydpe->i1pejp[iprop]];
          } */
  } } } }
  /* zero nodal forces and masses */
  for(inopo=0;inopo<ydn->nnopo;inopo++)
  { ydn->d1nmct[inopo]=R0;
    if(ydn->nnodim>0)ydn->d2nfc[0][inopo]=R0;
    if(ydn->nnodim>1)ydn->d2nfc[1][inopo]=R0;
    if(ydn->nnodim>2)ydn->d2nfc[2][inopo]=R0;
  }
  /* zero number of joint elements broken in the timestep */
  yde->netbrk=0;
  /* zero number of joint elements softened in the timestep */
  yde->netsft=0;
  /* Initializing i1ebrk */
  if(yde->i1ebrk == INT1NULL)
  { yde->i1ebrk=TalINT1(yde->melem);
    for (i=0;i<yde->melem;i++)
	{ yde->i1ebrk[i]=-1;
    }
  }
  /* Initializing i1esft */
  if(yde->i1esft == INT1NULL)
  { yde->i1esft=TalINT1(yde->melem);
    for (i=0;i<yde->melem;i++)
	{ yde->i1esft[i]=-1;
    }
  }
  /* Initializing d2ecbrk */
  if(yde->d2ecbrk==DBL2NULL)
  { yde->d2ecbrk=TalDBL2(ydn->mnodim,yde->melem);
    for(i=0;i<ydn->mnodim;i++)
    { for(j=0;j<yde->melem;j++)
        yde->d2ecbrk[i][j]=R0;
    }
  }
  /* Initializing d2ecbrk_NEW */
  if(yde->d2ecbrk_NEW==DBL2NULL)
  { yde->d2ecbrk_NEW=TalDBL2(ydn->mnodim,yde->melem);
    for(i=0;i<ydn->mnodim;i++)
    { for(j=0;j<yde->melem;j++)
        yde->d2ecbrk_NEW[i][j]=R0;
    }
  }
  /* Initializing d2ecsft */
  if(yde->d2ecsft==DBL2NULL)
  { yde->d2ecsft=TalDBL2(ydn->mnodim,yde->melem);
    for(i=0;i<ydn->mnodim;i++)
    { for(j=0;j<yde->melem;j++)
        yde->d2ecsft[i][j]=R0;
    }
  }
  /* initializing d1etbrk */
  if(yde->d1etbrk==DBL1NULL)
  { yde->d1etbrk=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
      yde->d1etbrk[i]=R0;
  }
  /* initializing d1etsft */
  if(yde->d1etsft==DBL1NULL)
  { yde->d1etsft=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
      yde->d1etsft[i]=R0;
  }
  /* initializing d1etmke */
  if(yde->d1etmke==DBL1NULL)
  { yde->d1etmke=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
      yde->d1etmke[i]=R0;
  }
  /* initializing d1elbrk */
  if(yde->d1elbrk==DBL1NULL)
  { yde->d1elbrk=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
      yde->d1elbrk[i]=R0;
  }
  /* Initializing d1efe */
  if(yde->d1efe==DBL1NULL)
  { yde->d1efe=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
	  yde->d1efe[i]=R0;
  }
  /* Initializing d1esftf */
  if(yde->i1esftf==INT1NULL)
  { yde->i1esftf=TalINT1(yde->melem);
    for(i=0;i<yde->melem;i++)
	  yde->i1esftf[i]=0;
  }
  /* Initializing d1ebrkf */
  if(yde->d1ebrkf==DBL1NULL)
  { yde->d1ebrkf=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
	  yde->d1ebrkf[i]=R0;
  }
  /* Initializing d1eike */
  if(yde->d1eike==DBL1NULL)
  { yde->d1eike=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
       yde->d1eike[i]=R0;
  }
  /* Initializing d1edke */
  if(yde->d1edke==DBL1NULL)
  { yde->d1edke=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
       yde->d1edke[i]=R0;
  }
  /* Initializing d1etike */
  if(yde->d1etike==DBL1NULL)
  { yde->d1etike=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
       yde->d1etike[i]=R0;
  }

  /* Initializing i1psup */
  if(ydpe->i1psup==INT1NULL)
  { ydpe->i1psup=TalINT1(ydpe->mprop);
    for (i=0;i<ydpe->mprop;i++)
    { yde->i1ebrk[i]=0; }
  }
  
//   /* Initializing i1pjhy with 0 if not specified in the input file */
//   if(ydpj->i1pjhy==INT1NULL)
//   { ydpj->i1pjhy=TalINT1(ydpj->mpjset);
//     for (i=0;i<ydpj->mpjset;i++)
//     { ydpj->i1pjhy[i]=0; }
//   }
  

  /* Initializing d2eldmg only if hysteretic joint model is used */
  if(ydpj->iusehy>0)
  { if(yde->d2eldmg==DBL2NULL)
    { yde->d2eldmg=TalDBL2(3,yde->melem);
      for(i=0;i<3;i++)
      { for(j=0;j<yde->melem;j++)
        yde->d2eldmg[i][j]=R0;
      }
    }
  }
  
  /* Initializing d2elstr only if rebars are used */
  if(ydsb->nsbar>0)
  { if(yde->d2elstr==DBL2NULL)
    { yde->d2elstr=TalDBL2(4,yde->melem);
      for(i=0;i<4;i++)
      { for(j=0;j<yde->melem;j++)
        yde->d2elstr[i][j]=R0;
      }
    }
  }
  
  /* Initializing d2nc0 (current coordinates at timestep 0, used to output displacements */
  if(ydn->d2nc0==DBL2NULL)
  { ydn->d2nc0=TalDBL2(4,ydn->mnopo);
    for(i=0;i<2;i++)
    { for(j=0;j<ydn->mnopo;j++)
      ydn->d2nc0[i][j]=ydn->d2ncc[i][j];
    }
  }
    
  /* At time step zero, if a mixed DFN type is used, find and assign crack type to all joint elements belonging to the DFN */
  if(ydc->ncstep == 0) 
  { if(ydfn->iusefn == 3)
    { for(jprop=0;jprop<ydpj->npjset;jprop++) // Loop over joint elements
      { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
        { for(ielem=0;ielem<yde->nelem;ielem++)
          { if(yde->i1elpr[ielem]==(jprop+(ydpe->nprop)))
            // Check if two edge nodes (of the joint element) belong to a DFN crack
       	    { xi_0 = ydn->d2nci[0][yde->i2elto[0][ielem]];  
      	      yi_0 = ydn->d2nci[1][yde->i2elto[0][ielem]];                
	      xi_1 = ydn->d2nci[0][yde->i2elto[1][ielem]];
	      yi_1 = ydn->d2nci[1][yde->i2elto[1][ielem]];
	      for(s=0; s<ydfn->mdfnfr; s++)
	      { for(k=0; k<ydfn->mdfnno; k++)
	        { if(ydfn->i2dfnn[k][s] >= 0)  
	          { if(xi_0 == ydn->d2nci[0][ydfn->i2dfnn[k][s]])
	            { if(yi_0 == ydn->d2nci[1][ydfn->i2dfnn[k][s]])
	              { for(r=0; r<ydfn->mdfnno; r++)
	                { if(ydfn->i2dfnn[r][s] >= 0)
		          { if(xi_1 == ydn->d2nci[0][ydfn->i2dfnn[r][s]])
	                    { if(yi_1 == ydn->d2nci[1][ydfn->i2dfnn[r][s]]) 
	                      { yde->i1edft[ielem] = ydfn->i1dfft[s]; //! Assign crack type to joint element (1 = broken, 2 = cohesive)
  } } } } } } } } } } } } } } }
  
  
   
  /* Apply support */
  for(iprop=0;iprop<ydpe->nprop;iprop++)
  { if(ydpe->i1psup[iprop]==1)
    { for(ielem=0;ielem<yde->nelem;ielem++)
      { if(yde->i1elpr[ielem]==iprop) //! Set initial coordinates equal to current coordinates (i.e., reset elastic deformation)
        { ydn->d2nci[0][(yde->i2elto[0][ielem])]=ydn->d2ncc[0][(yde->i2elto[0][ielem])];
          ydn->d2nci[0][(yde->i2elto[1][ielem])]=ydn->d2ncc[0][(yde->i2elto[1][ielem])];
          ydn->d2nci[0][(yde->i2elto[2][ielem])]=ydn->d2ncc[0][(yde->i2elto[2][ielem])];
          ydn->d2nci[1][(yde->i2elto[0][ielem])]=ydn->d2ncc[1][(yde->i2elto[0][ielem])];
          ydn->d2nci[1][(yde->i2elto[1][ielem])]=ydn->d2ncc[1][(yde->i2elto[1][ielem])];
          ydn->d2nci[1][(yde->i2elto[2][ielem])]=ydn->d2ncc[1][(yde->i2elto[2][ielem])];
      } } 
      ydpe->i1psup[iprop]=0;
   } }
   
  for(iprop=0;iprop<ydpe->nprop;iprop++)
  { if( (ydpe->i1ptyp[iprop])==(YTE2TRIELS) ||
        (ydpe->i1ptyp[iprop])==(YTE2PLANESTRESS) ||
        (ydpe->i1ptyp[iprop])==(YTE2PLANESTRAIN) )
    { Yfd2TRIELS(   /* small strain elastic triangle  */
      yde->nelem,
      iprop,
      ydpn->npnfact,ydpn->mpnset,ydpn->npnset,
      ydpn->d3pnfac,
      ydpe->i1ptyp[iprop],
      ydpe->d1peks[iprop],ydpe->d1pela[iprop],
      ydpe->d1pemu[iprop],ydpe->d1pero[iprop],ydpe->d1psem[iprop],
      ydpe->d1peem[iprop],ydpe->d1penu[iprop],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nfc[0],
      ydn->d2nfc[1],ydn->d1nmct  ,ydn->d2nvc[0],ydn->d2nvc[1],
      ydpn->d1pnaf , ydpn->d1pnap  ,ydpn->d1pnat  ,
      yde->i1elpr,ydn->i1nopr,yde->i2elto,
      ydo->nohys, ydo->dohyp , ydc->dctime,
      ydo->d1ohys, ydo->d1ohyt, ydo->d1ohyx, ydo->d1ohyy,
      ydo->i1ohyt, ydpn->npnset,yde->d1elfr,
      ydpe->i1usan[iprop], ydpe->d1peex[iprop], ydpe->d1peey[iprop],
      ydpe->d1pemx[iprop], ydpe->d1pemy[iprop], ydpe->d1peg[iprop],
      ydis->iuseis, ydis->dcstxx, ydis->dcstxy, ydis->dcstyy,
      ydis->dcsyxx, ydis->dcsyxy, ydis->dcsyyy, ydis->dcsrfy,
      ydpn->i1pnfx, ydpn->i1pnfy,
      ydpe->i1pexc, ydn->i1nowe, ydhf->iusehf,
      ydsb->nsbar,yde->d2elstr
      );
    }
    else if((ydpe->i1ptyp[iprop])==(YTE2TRISOF))
    { Yfd2TRISOF(   /* small strain softening triangle  */
      yde->nelem,
      iprop,      
      ydpe->d1peks[iprop],ydpe->d1pela[iprop],
      ydpe->d1pemu[iprop],ydpe->d1pero[iprop],ydpe->d1psem[iprop],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nfc[0],
      ydn->d2nfc[1],ydn->d1nmct,ydn->d2nvc[0],ydn->d2nvc[1],
      yde->d2elst[ydpe->i1psde[iprop]], 
      yde->i1elpr,yde->i2elto, yde->d1elfs
      );
    }
  }
  for(jprop=0;jprop<ydpj->npjset;jprop++)
  { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
    { Yfd2JOINTS(  /* joint element  */
      yde->nelem,(jprop+ydpe->nprop),
      ydpj->d1pjft[jprop], ydpj->d1pjgf[jprop], ydpj->d1pjgs[jprop],
      ydpj->d1pjco[jprop], ydpj->d1pjfr[jprop], ydpj->d1pjpe[jprop],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nfc[0],ydn->d2nfc[1],
      ydn->d2nvc[0],ydn->d2nvc[1],
      yde->i1elpr,yde->i2elto,yde->d1elfs,ydc->dctime,ydc->dcstec,&(yde->nebrk),&(yde->netbrk),yde->i1ebrk,yde->d2ecbrk,yde->d2ecbrk_NEW,
      yde->d1etbrk,yde->d1elbrk,yde->d1efe,&(yde->nesft),&(yde->netsft),yde->i1esft,yde->d2ecsft,
      yde->d1etsft,yde->i1esftf,yde->d1ebrkf,yde->d1eike,yde->d1edke,ydn->d1nmct,yde->d1etmke,
      yde->i2elnext,yde->i2eledge,
      ydpe->i1pexc,
      ydpj->d1usaf[jprop],ydpj->d1pjal[jprop],
      ydpj->d1pjcr[jprop],ydpj->d1pjfd[jprop],ydpj->d1pjtr[jprop],ydpj->d1pjgr[jprop],ydpj->d1pjsr[jprop],ydn->i1nowe,ydn->i2noid,
      ydfn->iusefn,yde->i1edfnf,ydfn->ddfnft,ydfn->ddfnco,ydfn->ddfngf,ydfn->ddfngs,ydhf->iusehf,
      yde->i1edft,yde->d1etike,ydsm->iusesm,ydsm->dctwle,
      ydpj->iusehy,yde->d2eldmg);
    }
  }
}


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
/* File  Yrb.c */
#include "Yproto.h"

/*****************  RBAR STEEL ELEMENTS *********************/
static void Yrb2CREATE(   /* created 1RBAR element  */
            nelem,  nsbar,  mrepo, n0repo, 
            d1nccx, d1nccy, d1ncix, d1nciy, d1nvcx, d1nvcy,
            d1sicAx,d1sicAy,d1sicBx,d1sicBy,d1rcigx,d1rcigy,
            d1rccgx,d1rccgy,d1rvcgx, d1rvcgy,
            d1riLcK,d1riLcE,d1riLcZ,
            i1elpr, i2elto, i1ptyp, i1rrpn, i1srpf,i1rmyel,
            d1spea,  
            i2relto,i1rprop, i1refbar, i1sbpr,
            d1rsctrx, d1rsctry, i2rbedn
            ) 
  INT    nelem; INT   nsbar; INT   mrepo; INT *n0repo;
  DBL  *d1nccx;  DBL *d1nccy;  DBL *d1ncix;   DBL *d1nciy;   DBL *d1nvcx;  DBL *d1nvcy; 
  DBL *d1sicAx;  DBL  *d1sicAy; DBL *d1sicBx; DBL *d1sicBy;  DBL *d1rcigx; DBL *d1rcigy; 
  DBL  *d1rccgx; DBL *d1rccgy;  DBL  *d1rvcgx; DBL *d1rvcgy; 
  DBL *d1riLcK;  DBL *d1riLcE;  DBL  *d1riLcZ;
  INT  *i1elpr;  INT **i2elto; INT *i1ptyp;   INT *i1rrpn;   INT  *i1srpf;
  INT *i1rmyel;  
  DBL *d1spea;
  INT **i2relto; INT *i1rprop; INT *i1refbar; INT *i1sbpr;   
  DBL *d1rsctrx; DBL *d1rsctry; INT **i2rbedn;
{ INT ielem,isbar,iprop,irpo,bezs,indk;
  DBL r0ix,r1ix,r2ix,r0iy,r1iy,r2iy,x0b,x1b,y0b,y1b;
  DBL v0ix,v1ix,v2ix,v0iy,v1iy,v2iy;
  DBL r0irx,r0iry,r1irx,r1iry,r2irx,r2iry,r0cx,r1cx,r2cx,r0cy,r1cy,r2cy;
  DBL rsxb, rsyb, rsxb0, rsyb0, nsxb0, nsyb0, VTemp;
  DBL p0ix, p0iy, p1ix, p1iy,p0cx, p0cy, p1cx, p1cy,p0vx,p0vy,p1vx,p1vy; 
  DBL dksi0,deta0,dzeta0,dksi1,deta1,dzeta1;
  DBL d0s,d1s,d2s,alfa0,beta0,alfa1,beta1,minXb,maxXb,minYb,maxYb;
  INT i,j,k=0,ncrepo=0,nrjoint=0;
  INT *i1rmyelT;

  DBL small=1.0e-8;    //1.0e-15;
  //DBL small=0.1; 
  i1rmyelT = INT1NULL;
  
  DBL elsize;
  
  /*static FILE *out1=FILENULL;
  if(out1 == FILENULL)
  {    out1=fopen("Yrb2CREATE.txt", "a");
  }*/
  
  /*Intersection Rbar triang Fem */
  irpo=(*n0repo);
  
  //printf("YRB2CREATE\n");

  for(isbar=0;isbar<nsbar;isbar++)
  { INT bfst = 0;
    for(ielem=0;ielem<nelem;ielem++)
    { //printf("YRB2CREATE-0\n");
      //printf("ielem: %d \t nelem: %d\n",ielem, nelem);
      iprop=i1elpr[ielem];
      //printf("YRB2CREATE-0-a\n");
      //printf("d1ncix[i2elto[0][ielem]]: %f\t d1nciy[i2elto[0][ielem]]: %f\t d1ncix[i2elto[1][ielem]]: %f\t d1nciy[i2elto[1][ielem]]: %f\n",d1ncix[i2elto[0][ielem]],d1nciy[i2elto[0][ielem]],d1ncix[i2elto[1][ielem]],d1nciy[i2elto[1][ielem]]);
      //printf("iprop: %d\n",iprop);
      //printf("i1ptyp[iprop]: %d\n",i1ptyp[iprop]);
      if(iprop<0) //! If ielem is a DFN joint (iprop = ipropj-YIPROPMAX)
      { continue; }
      if((i1ptyp[iprop])==(YTE2TRIELS) || (i1ptyp[iprop])==(YTE2TRIRIG) || 
         (i1ptyp[iprop])==(YTE2PLANESTRESS) || (i1ptyp[iprop])==(YTE2PLANESTRAIN))
      { //printf("YRB2CREATE-0-b\n");
	r0ix=d1ncix[i2elto[0][ielem]];
        r1ix=d1ncix[i2elto[1][ielem]];
        r2ix=d1ncix[i2elto[2][ielem]];
        r0iy=d1nciy[i2elto[0][ielem]];
        r1iy=d1nciy[i2elto[1][ielem]];           
        r2iy=d1nciy[i2elto[2][ielem]];

	//printf("Yrb2CREATE\n");

        x0b=d1sicAx[isbar];
        x1b=d1sicBx[isbar];
        y0b=d1sicAy[isbar];
        y1b=d1sicBy[isbar];
	
	//printf("Hello\n");
	//printf("isbar: %d \t ielem: %d\n",isbar,ielem);
	//printf("r0ix: %f \t r0iy: %f \t r1ix: %f \t r1iy: %f \t r2ix: %f \t r2iy: %f\n",r0ix,r0iy,r1ix,r1iy,r2ix,r2iy);
	//printf("x0b: %f \t y0b: %f \t x1b: %f \t y1b: %f \n",x0b,y0b,x1b,y1b);
	//printf("n0repo: |%d|\n",*n0repo);
	//printf("irpo: %d\n",irpo);
	
        rsxb=x1b-x0b;
        rsyb=y1b-y0b;
        
        //printf("rsxb: %f \t rsyb: %f \n",rsxb,rsyb);

        r0irx=r0ix-x0b;
        r0iry=r0iy-y0b;
        r1irx=r1ix-x0b;
        r1iry=r1iy-y0b;
        r2irx=r2ix-x0b;
        r2iry=r2iy-y0b;

        //printf("r0irx: %f \t r0iry: %f \n",r0irx,r0iry);
        //printf("r1irx: %f \t r1iry: %f \n",r1irx,r1iry);
        //printf("r2irx: %f \t r2iry: %f \n",r2irx,r2iry);
 
        V2DNor(rsxb0,rsyb0,rsxb,rsyb);
        //rsxb0=rsxb/(rsxb0*rsxb0+)
        //printf("rsxb0: %f \t rsyb0: %f \n",rsxb0,rsyb0);
        nsxb0=-rsyb0;
        nsyb0=rsxb0;

        V2DDot(d0s,nsxb0,nsyb0,r0irx,r0iry);
        V2DDot(d1s,nsxb0,nsyb0,r1irx,r1iry);
        V2DDot(d2s,nsxb0,nsyb0,r2irx,r2iry);

        //printf("d0s: %f \t d1s: %f \t d2s: %f \n",d0s,d1s,d2s);
         
	//printf("YRB2CREATE-1\n");
        /* if cross reference points*/
        if(((d0s>R0)&&(d1s>R0)&&(d2s>R0))||((d0s<R0)&&(d1s<R0)&&(d2s<R0)))
        {    bezs=0;
        }
        if(((d0s>R0)&&(d2s>R0)&&(d1s<R0))||((d0s<R0)&&(d2s<R0)&&(d1s>R0)))
        { alfa0=((ABS(d1s))/(ABS(d0s)+ABS(d1s)/*+small*/));
          beta0=1-alfa0;
          dksi0=alfa0;
          deta0=beta0;
          dzeta0=R0;
          alfa1=(ABS(d1s)/(ABS(d2s)+ABS(d1s)/*+small*/));
          beta1=1-alfa1;
          dksi1=R0;
          deta1=beta1;
          dzeta1=alfa1;
          if(bfst==0) bfst=1;
          bezs=1;
          i2rbedn[0][irpo]=i2elto[0][ielem];
          i2rbedn[1][irpo]=i2elto[1][ielem];
          i2rbedn[0][irpo+1]=i2elto[1][ielem];
          i2rbedn[1][irpo+1]=i2elto[2][ielem];
        }
        else if(((d0s>R0)&&(d1s>R0)&&(d2s<R0))||((d0s<R0)&&(d1s<R0)&&(d2s>R0)))
        { alfa0=(ABS(d2s)/(ABS(d0s)+ABS(d2s)/*+small*/));
          beta0=1-alfa0;
          dksi0=alfa0;
          deta0=R0;
          dzeta0=beta0;
          alfa1=(ABS(d2s)/(ABS(d1s)+ABS(d2s)/*+small*/));
          beta1=1-alfa1;
          dksi1=R0;
          deta1=alfa1;
          dzeta1=beta1;
          if(bfst==0) bfst=1;
          bezs=1;
          i2rbedn[0][irpo]=i2elto[2][ielem];
          i2rbedn[1][irpo]=i2elto[0][ielem];
          i2rbedn[0][irpo+1]=i2elto[1][ielem];
          i2rbedn[1][irpo+1]=i2elto[2][ielem];
        }
        else if(((d1s>R0)&&(d2s>R0)&&(d0s<R0))||((d1s<R0)&&(d2s<R0)&&(d0s>R0)))
        { alfa0=((ABS(d0s))/(ABS(d2s)+ABS(d0s)/*+small*/));
          beta0=1-alfa0;
          dksi0=beta0;
          deta0=R0;
          dzeta0=alfa0;
          alfa1=(ABS(d0s)/(ABS(d0s)+ABS(d1s)/*+small*/));
          beta1=1-alfa1;
          dksi1=beta1;
          deta1=alfa1;
          dzeta1=R0;
          if(bfst==0) bfst=1;
          bezs=1;
          i2rbedn[0][irpo]=i2elto[2][ielem];
          i2rbedn[1][irpo]=i2elto[0][ielem];
          i2rbedn[0][irpo+1]=i2elto[0][ielem];
          i2rbedn[1][irpo+1]=i2elto[1][ielem];
        } 
        //printf("YRB2CREATE-2\n");
        //printf("dksi0: %f \t deta0: %f \t dzeta0: %f \n",dksi0,deta0,dzeta0);
	//printf("dksi1: %f \t deta1: %f \t dzeta1: %f \n",dksi1,deta1,dzeta1);
        
        if (bezs==1)  
        { d1riLcK[irpo]=dksi0;
          d1riLcE[irpo]=deta0;
          d1riLcZ[irpo]=dzeta0;
          d1riLcK[irpo+1]=dksi1;
          d1riLcE[irpo+1]=deta1;
          d1riLcZ[irpo+1]=dzeta1;
          if(bfst==1)
          { i1srpf[isbar]=irpo;
            bfst=2;
          }
          i1rmyel[irpo]=ielem;
          i1rmyel[irpo+1]=ielem;
	  /*if(ielem==153)
	  { printf("r0ix: %f \t r0iy: %f \t r1ix: %f \t r1iy: %f \t r2ix: %f \t r2iy: %f\n",r0ix,r0iy,r1ix,r1iy,r2ix,r2iy); 
            printf("irpo: %d \t irpo+1: %d\n",irpo,irpo+1); }*/
          //fprintf(out1,"ielem: %d \n",ielem);
	  //fprintf(out1,"r0ix: %f \t r0iy: %f \t r1ix: %f \t r1iy: %f \t r2ix: %f \t r2iy: %f\n",r0ix,r0iy,r1ix,r1iy,r2ix,r2iy); 
          //fprintf(out1,"irpo: %d \t irpo+1: %d\n\n",irpo,irpo+1);
          
          i1rrpn[irpo]=irpo; 
          i1rrpn[irpo+1]=irpo+1;


         /* initial global coodinates Rbar elements */
          p0ix = r0ix*dksi0+r1ix*deta0+r2ix*dzeta0;
          p0iy = r0iy*dksi0+r1iy*deta0+r2iy*dzeta0;
          p1ix = r0ix*dksi1+r1ix*deta1+r2ix*dzeta1;
          p1iy = r0iy*dksi1+r1iy*deta1+r2iy*dzeta1;

          //printf("p0ix: %f \t p0iy: %f \t p1ix: %f  \t p1iy: %f\n\n",p0ix,p0iy,p1ix,p1iy);
              
          minXb=MINIM(x0b,x1b);
          maxXb=MAXIM(x0b,x1b);

          minYb=MINIM(y0b,y1b);
          maxYb=MAXIM(y0b,y1b);
              
          indk=1;

          if (x0b==x1b)
          { (p0ix=x0b)&&(p1ix=x1b);
            indk=1;
          }
          else if (y0b==y1b)
          { (p0iy=y0b)&&(p1iy=y1b);
            indk=1;
          }
          else if (indk=1);
          { if( (((p0ix>=minXb)&&(p0ix<=maxXb))&&((p0iy>=minYb)&&(p0iy<=maxYb)))&&
                (((p1ix>=minXb)&&(p1ix<=maxXb))&&((p1iy>=minYb)&&(p1iy<=maxYb))))
            { indk=0; }
          }
          
          d1rcigx[*n0repo]=p0ix;
          d1rcigy[*n0repo]=p0iy;
          d1rcigx[*n0repo+1]=p1ix;
          d1rcigy[*n0repo+1]=p1iy;

          /*relation betwen ref. point and bar*/ 
          i1refbar[*n0repo]=isbar;
          i1refbar[*n0repo+1]=isbar;

          /* properties nrepo*/
          i1rprop[*n0repo]=i1sbpr[isbar];
          i1rprop[*n0repo+1]=i1sbpr[isbar]; 
           
          

          /* current global coodinates Rbar elements */
       
          r0cx=d1nccx[i2elto[0][ielem]];
          r1cx=d1nccx[i2elto[1][ielem]];
          r2cx=d1nccx[i2elto[2][ielem]];
          r0cy=d1nccy[i2elto[0][ielem]];
          r1cy=d1nccy[i2elto[1][ielem]];           
          r2cy=d1nccy[i2elto[2][ielem]];
     
          p0cx = r0cx*dksi0+r1cx*deta0+r2cx*dzeta0;
          p0cy = r0cy*dksi0+r1cy*deta0+r2cy*dzeta0;
          p1cx = r0cx*dksi1+r1cx*deta1+r2cx*dzeta1;
          p1cy = r0cy*dksi1+r1cy*deta1+r2cy*dzeta1;

          d1rccgx[*n0repo]=p0cx;
          d1rccgy[*n0repo]=p0cy;
          d1rccgx[*n0repo+1]=p1cx;
          d1rccgy[*n0repo+1]=p1cy;
          
          /* initial velocity */

          v0ix=d1nvcx[i2elto[0][ielem]];
          v1ix=d1nvcx[i2elto[1][ielem]];
          v2ix=d1nvcx[i2elto[2][ielem]];
          v0iy=d1nvcy[i2elto[0][ielem]];
          v1iy=d1nvcy[i2elto[1][ielem]];           
          v2iy=d1nvcy[i2elto[2][ielem]];

          p0vx = v0ix*dksi0+v1ix*deta0+v2ix*dzeta0;
          p0vy = v0iy*dksi0+v1iy*deta0+v2iy*dzeta0;
          p1vx = v0ix*dksi1+v1ix*deta1+v2ix*dzeta1;
          p1vy = v0iy*dksi1+v1iy*deta1+v2iy*dzeta1;

          d1rvcgx[*n0repo]=p0vx;
          d1rvcgy[*n0repo]=p0vy;
          d1rvcgx[*n0repo+1]=p1vx;
          d1rvcgy[*n0repo+1]=p1vy;

          d1rsctrx[*n0repo]=d1rccgx[*n0repo+1];
          d1rsctry[*n0repo]=d1rccgy[*n0repo+1];
          d1rsctrx[*n0repo+1]=d1rccgx[*n0repo];
          d1rsctry[*n0repo+1]=d1rccgy[*n0repo];
        }
        //printf("YRB2CREATE-3\n");
        if ((bezs==1)&&(indk==0))  
        { //printf("YRB2CREATE-3-a\n");
          irpo=irpo+2;
          if (irpo>=mrepo)
          { //printf("YRB2CREATE-3-b\n");
	    CHRw(stderr,"Yjd: MREPO too small"); CHRwcr(stderr);
            exit(1);
          }
          //printf("YRB2CREATE-3-c\n");
          (*n0repo)=irpo;
	  //printf("YRB2CREATE-3-c-1\n");
        }
        else
        { //printf("YRB2CREATE-3-d\n");
	  irpo=irpo;
          (*n0repo)=irpo;
	  //printf("YRB2CREATE-3-e\n");
        } 
      } 
    }                 
  }
  //printf("n0repo: |%d|\n",n0repo);
  //printf("irpo: %d\n",irpo);
  //printf("YRB2CREATE-4\n");
  
  //! Determine the total number of reference points
  for(i=0; i1rmyel[i]>=0; i++)  
  { //fprintf(out1,"i: %d \t i1rmyel[i]: %d\n",i,i1rmyel[i]);
    ncrepo=ncrepo+1;
  }
  //fprintf(out1,"ncrepo: %d \n",ncrepo);
  //fprintf(out1,"ncrepo: %d \n",ncrepo);
  //! Create a temporary array to store the element IDs associated to each reference point
  i1rmyelT =TalINT1(ncrepo);
  for (i=0; i<ncrepo;i++)
  { i1rmyelT[i]=i1rmyel[i];
    //fprintf(out1,"i: %d \t i1rmyelT[i]: %d\n",i,i1rmyelT[i]);
  }
  //printf("YRB2CREATE-5\n");
  //! Determine the topology array for the 1D joint elements 
  for (i=0; i<ncrepo; i++) //! Loop over all reference points
  { /* Compute the average edge length of the associated triangle */ 
    //fprintf(out1,"First loop\n");
    r0ix=d1ncix[i2elto[0][i1rmyel[i]]];
    r1ix=d1ncix[i2elto[1][i1rmyel[i]]];
    r2ix=d1ncix[i2elto[2][i1rmyel[i]]];
    r0iy=d1nciy[i2elto[0][i1rmyel[i]]];
    r1iy=d1nciy[i2elto[1][i1rmyel[i]]];           
    r2iy=d1nciy[i2elto[2][i1rmyel[i]]];
    elsize=(SQRT((r0ix-r1ix)*(r0ix-r1ix)+(r0iy-r1iy)*(r0iy-r1iy))+SQRT((r0ix-r2ix)*(r0ix-r2ix)+(r0iy-r2iy)*(r0iy-r2iy))+SQRT((r2ix-r1ix)*(r2ix-r1ix)+(r2iy-r1iy)*(r2iy-r1iy)))/3;
        
    for (j=0; j<ncrepo; j++) //! Loop over all reference points
    { //fprintf(out1,"Second loop\n");
      if((i!=j) && //! If two reference points are distinct 
         //(ABS(d1rcigx[i]-d1rcigx[j])<small) && //! if the x distance is less than "small"
         //(ABS(d1rcigy[i]-d1rcigy[j])<small) && //! if the y distance is less than "small"
         (ABS(d1rcigx[i]-d1rcigx[j])<(elsize/5)) && //! if the x distance is less than "elsize/5" (Modification introduced to insert 1D joints also between "distinct" triangles)
         (ABS(d1rcigy[i]-d1rcigy[j])<(elsize/5)) && //! if the y distance is less than "elsize/5" (Modification introduced to insert 1D joints also between "distinct" triangles)
         (i1rmyelT[i]!=i1rmyelT[j]) && //! if the reference points belong to distinct triangles 
         (i1refbar[i]==i1refbar[j])) //! if the reference points belong to the same bar
      { //fprintf(out1,"Found, i: %d\n",k);
	i2relto[0][k]=i;
        i2relto[1][k]=j;
        i1rmyelT[i]=-1;
        i1rmyelT[j]=-1;
        //fprintf(out1,"i: %d \t j: %d \t k: %d \n",i,j,k);
        //fprintf(out1,"i2relto[0][k]: %d \t i2relto[1][k]: %d \n",i2relto[0][k],i2relto[1][k]);
        //fprintf(out1,"i1rmyelT[i]: %d \t i1rmyelT[j]: %d \n",i1rmyelT[i],i1rmyelT[j]);
        k=k+1;
         
      }
    }
  }
  //printf("YRB2CREATE-6\n");
  /* Setting the initial coordinates of rebar elements = current coordinates  */
  /* Added to Y-RC to avoid initial stress (i.e., at time step zero) in rebars when medium is pre-stressed (e.g., with in-situ stress) */  
  for (i=0; i<ncrepo; i++) //! Loop over reference points
  { d1rcigx[i]=d1rccgx[i];
    d1rcigy[i]=d1rccgy[i];
  }
  
  free (i1rmyelT);
  
  //printf("Bye YRB2CREATE\n");
  //printf("Hello Yrb2CREATE\n");
}

static void Yrb2CARENT( /*carent coordinates*/
            nelem,nsbar,n0repo,
            d1nccx,d1nccy,d1nvcx,d1nvcy,
            d1rccgx,d1rccgy,d1rvcgx,d1rvcgy,
            d1riLcK,d1riLcE,d1riLcZ,
            i2elto,i1rrpn,i1srpf,i1rmyel,
            d1rsctrx,d1rsctry
            ) 
  INT nelem; INT nsbar; INT *n0repo;
  DBL *d1nccx; DBL *d1nccy; DBL *d1nvcx; DBL *d1nvcy;
  DBL *d1rccgx; DBL *d1rccgy;  DBL *d1rvcgx; DBL *d1rvcgy; 
  DBL *d1riLcK;  DBL *d1riLcE;  DBL  *d1riLcZ;
  INT **i2elto; INT *i1rrpn;   INT  *i1srpf;  INT *i1rmyel; 
  DBL *d1rsctrx; DBL *d1rsctry;
{ INT imyelem,ip;
  DBL r0cx,r1cx,r2cx,r0cy,r1cy,r2cy;
  DBL v0cx,v1cx,v2cx,v0cy,v1cy,v2cy;
  DBL p0cx, p0cy, p1cx, p1cy; 
  DBL p0vx, p0vy, p1vx, p1vy; 
  DBL dksi0,deta0,dzeta0,dksi1,deta1,dzeta1;
 
  /*Intersection Rbar triang Fem */
  
  //printf("Yrb2CARENT\n");
  
  for(ip=0;ip<(*n0repo);ip=ip+2)
  { imyelem=i1rmyel[ip];
    { r0cx=d1nccx[i2elto[0][imyelem]];
      r1cx=d1nccx[i2elto[1][imyelem]];
      r2cx=d1nccx[i2elto[2][imyelem]];
      r0cy=d1nccy[i2elto[0][imyelem]];
      r1cy=d1nccy[i2elto[1][imyelem]];           
      r2cy=d1nccy[i2elto[2][imyelem]];

      v0cx=d1nvcx[i2elto[0][imyelem]];
      v1cx=d1nvcx[i2elto[1][imyelem]];
      v2cx=d1nvcx[i2elto[2][imyelem]];
      v0cy=d1nvcy[i2elto[0][imyelem]];
      v1cy=d1nvcy[i2elto[1][imyelem]];           
      v2cy=d1nvcy[i2elto[2][imyelem]];

      /* [mrdim][mrepo] reference point coordinate local  */
        
      dksi0=d1riLcK[ip];
      deta0=d1riLcE[ip];
      dzeta0=d1riLcZ[ip];
      dksi1=d1riLcK[ip+1];
      deta1=d1riLcE[ip+1];
      dzeta1=d1riLcZ[ip+1];

      /* current global coodinates Rbar elements */

      p0cx = r0cx*dksi0+r1cx*deta0+r2cx*dzeta0;
      p0cy = r0cy*dksi0+r1cy*deta0+r2cy*dzeta0;
      p1cx = r0cx*dksi1+r1cx*deta1+r2cx*dzeta1;
      p1cy = r0cy*dksi1+r1cy*deta1+r2cy*dzeta1;

      p0vx = v0cx*dksi0+v1cx*deta0+v2cx*dzeta0;
      p0vy = v0cy*dksi0+v1cy*deta0+v2cy*dzeta0;
      p1vx = v0cx*dksi1+v1cx*deta1+v2cx*dzeta1;
      p1vy = v0cy*dksi1+v1cy*deta1+v2cy*dzeta1;


      d1rccgx[ip]=p0cx;
      d1rccgy[ip]=p0cy;
      d1rccgx[ip+1]=p1cx;
      d1rccgy[ip+1]=p1cy;

      d1rvcgx[ip]=p0vx;
      d1rvcgy[ip]=p0vy;
      d1rvcgx[ip+1]=p1vx;
      d1rvcgy[ip+1]=p1vy;

      /* carent vector ref.points */
      d1rsctrx[ip]=d1rccgx[ip+1];
      d1rsctry[ip]=d1rccgy[ip+1];
      d1rsctrx[ip+1]=d1rccgx[ip];
      d1rsctry[ip+1]=d1rccgy[ip];
    }                   
  }
}  

static void Yrb2FORCE(  /* nodal force RBAR elements   */
            nelem, 
            iprop,
            nsbar, n0repo, 
            young,d1rcigx,d1rcigy,d1rccgx,d1rccgy,
            d1nfcx,d1nfcy, d1riLcK,d1riLcE,d1riLcZ,
            d1ncix,d1nciy,
            i1srpf, i1rprop, i1rrpn, i1rmyel, i1refbar,
            i2elto, d1spea,
            nohys, d1ohyx, d1ohyy, d1ohys, d1ohyt,dohyp,
            i1ohyt,
            dctime,ncstep,
            d1rbsig, d1rbfrc, d1rbstr,
            i1sbac
            ) 
   
  INT  nelem;
  INT iprop;
  INT  nsbar;   INT *n0repo;
  DBL  young;   DBL *d1rcigx; DBL *d1rcigy; DBL *d1rccgx; DBL *d1rccgy; 
  DBL *d1nfcx;  DBL *d1nfcy;  DBL *d1riLcK; DBL *d1riLcE; DBL *d1riLcZ;
  DBL *d1ncix;  DBL *d1nciy; 
  INT *i1srpf;  INT *i1rprop;   INT *i1rrpn;  INT *i1rmyel; INT *i1refbar;
  INT **i2elto; DBL *d1spea;
  INT nohys; DBL *d1ohyx; DBL *d1ohyy; DBL *d1ohys; DBL *d1ohyt; DBL dohyp; 
  INT *i1ohyt;
  DBL dctime; INT ncstep;  
  DBL *d1rbsig; DBL *d1rbfrc; DBL *d1rbstr;
  INT *i1sbac;

{ INT imyelem,ip,in,jn,kn,ihys;
  DBL drxi,dryi,dLi,drxc,dryc,dLc,dstrain,dsigma,dforce;
  DBL dksi0,deta0,dzeta0,dksi1,deta1,dzeta1,dpea;
  DBL r0x,r1x,r2x,r0y,r1y,r2y,rpx,rpy,v0,v1,v2,stprev;
      
  //printf("Yrb2FORCE\n");
  
 /*  force RBAR elements  */  
  for(ip=0;ip<(*n0repo);ip=ip+2)
  { if(i1rprop[ip]==iprop)
    { imyelem=i1rmyel[ip];
      drxi=d1rcigx[ip+1]-d1rcigx[ip];
      dryi=d1rcigy[ip+1]-d1rcigy[ip];
      dLi=SQRT(drxi*drxi+dryi*dryi);
      drxc=d1rccgx[ip+1]-d1rccgx[ip];
      dryc=d1rccgy[ip+1]-d1rccgy[ip];
      dLc=SQRT(drxc*drxc+dryc*dryc);

      dpea=d1spea[i1refbar[ip]];
      //printf("ncstep: %d \t ip: %d \t imyelem: %d  \t i1sbac[i1refbar[ip]]: %d \n",ncstep,ip,imyelem,i1sbac[i1refbar[ip]]);
      //printf("ncstep: %d \t ip: %d \t drxi: %.6f  \t dryi: %.6f \n",ncstep,ip,drxi,dryi);
      //printf("ncstep: %d \t ip: %d \t drxc: %.6f  \t dryc: %.6f \n",ncstep,ip,drxc,dryc);
      
      drxi=drxi/dLi;
      dryi=dryi/dLi;
      drxc=drxc/dLc;
      dryc=dryc/dLc;
      dstrain =(dLc-dLi)/dLi;
       
      /* Compute stress and force only if the rebar has been activated */
      if(i1sbac[i1refbar[ip]]==2)
      { dsigma=young*dstrain;
        dforce=dsigma*dpea;
      }
      else
      { dsigma=0.0;
        dforce=0.0;
      }      

      /* Storing dsigma, dforce and dstrain for output */
      d1rbsig[ip]=dsigma;
      d1rbfrc[ip]=dforce;
      d1rbstr[ip]=dstrain;
      
      //printf("ip: %d \t imyelem: %d  \t dsigma: %f \n",ip,imyelem,dsigma);
      
      in=i2elto[0][imyelem];
      jn=i2elto[1][imyelem];
      kn=i2elto[2][imyelem];
      dksi0=d1riLcK[ip];
      deta0=d1riLcE[ip];
      dzeta0=d1riLcZ[ip];
      dksi1=d1riLcK[ip+1];
      deta1=d1riLcE[ip+1];
      dzeta1=d1riLcZ[ip+1];
      d1nfcx[in]=d1nfcx[in]+dforce*drxc*dksi0;
      d1nfcx[jn]=d1nfcx[jn]+dforce*drxc*deta0;
      d1nfcx[kn]=d1nfcx[kn]+dforce*drxc*dzeta0;
      d1nfcy[in]=d1nfcy[in]+dforce*dryc*dksi0;
      d1nfcy[jn]=d1nfcy[jn]+dforce*dryc*deta0;
      d1nfcy[kn]=d1nfcy[kn]+dforce*dryc*dzeta0;
      d1nfcx[in]=d1nfcx[in]-dforce*drxc*dksi1;
      d1nfcx[jn]=d1nfcx[jn]-dforce*drxc*deta1;
      d1nfcx[kn]=d1nfcx[kn]-dforce*drxc*dzeta1;
      d1nfcy[in]=d1nfcy[in]-dforce*dryc*dksi1;
      d1nfcy[jn]=d1nfcy[jn]-dforce*dryc*deta1;
      d1nfcy[kn]=d1nfcy[kn]-dforce*dryc*dzeta1;

      /* output history states */
//       r0x = d1ncix[in];
//       r0y = d1nciy[in];
//       r1x = d1ncix[jn];
//       r1y = d1nciy[jn];
//       r2x = d1ncix[kn];
//       r2y = d1nciy[kn];
//       for(ihys=0; ihys<nohys; ihys++)
//       { rpx=d1ohyx[ihys];	/* x coordinate of point P */
//         rpy=d1ohyy[ihys];	/* y coordinate of point P */
//         V2DCro(v0,(r1x-r0x),(r1y-r0y),(rpx-r0x),(rpy-r0y));
//         V2DCro(v1,(r2x-r1x),(r2y-r1y),(rpx-r1x),(rpy-r1y));
//         V2DCro(v2,(r0x-r2x),(r0y-r2y),(rpx-r2x),(rpy-r2y));
// 
//         if((v0>R0)&&(v1>R0)&&(v2>R0))		/* if point is inside the triangle */
//         { if(i1ohyt[ihys]==(YFLDSBAR))
//           { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
//             if(/*(ABS(dforce-stprev))>=dohyp ||*/ !(ncstep%icoutf) )
//             { d1ohyt[ihys] = dctime;		/* output history time  */
//               d1ohys[ihys] = dsigma;		/* output history state */  
//             }
//           }
//         }
//       }
    }
  }
} 

static void Yrb1RJOINT(  /* 1D joint element  */
            nbrjointrb, nelem,
            iprop,
            young, dcstec,
            d1rccgx,d1rccgy,
            d1nccx,d1nccy,d1nfcx,d1nfcy,d1nvcx,
            d1nvcy,d1peftc,d1pefsc,d1pegfs,d1pepec,d1spea,
            d1sdiam,d1smmdiam,d1crlcr,i1sbpr,
            d1sfc,d1mpsfc,
            d1sfy, d1epssh, d1sfu, d1epsu, d1sfbr, d1epsbr, 
            i1elpr,i1rprop,i2elto,i2relto,mrepo,i1refbar,
            i1rmyel,
            i1myjoint,
            d1rsctrx,d1rsctry,d1rvcgx,d1rvcgy,
            d1pela,d1pemu,nprop,
            dctime, ncstep,
            i1ptyp,
            d1rjsig, d1rjtau,
            i1petyp, d1peem, d1peex, d1peey, i1usan,
            i1sbty, d1elfs,
            i1pexc, i1sbac,
            d1stkn, d1stkt,
            d1rjslpnor, d1rjdelnor,
            d1styns, d1strns,
            i1rjstnor,
            d2elstr,
            d1stcoh,d1stfri,
            d2rjfrv,i2rbedn
            ) 
  INT nbrjointrb; INT nelem;
  INT iprop;
  DBL  young; DBL  dcstec;
  DBL *d1rccgx; DBL *d1rccgy; 
  DBL *d1nccx; DBL *d1nccy; DBL *d1nfcx; DBL *d1nfcy; DBL *d1nvcx;
  DBL *d1nvcy; DBL *d1peftc; DBL *d1pefsc; DBL *d1pegfs; DBL *d1pepec;DBL *d1spea;
  DBL *d1sdiam; DBL *d1smmdiam; DBL *d1crlcr; INT *i1sbpr;
  DBL *d1sfc;   DBL *d1mpsfc;
  DBL *d1sfy; DBL *d1epssh; DBL *d1sfu; DBL *d1epsu; DBL *d1sfbr; DBL *d1epsbr; 
  INT *i1elpr; INT *i1rprop; INT **i2elto; INT **i2relto; INT mrepo; INT *i1refbar;
  INT *i1rmyel;
  INT *i1myjoint;
  DBL *d1rsctrx; DBL *d1rsctry; DBL *d1rvcgx; DBL *d1rvcgy; 
  DBL *d1pela; DBL *d1pemu; INT nprop;
  DBL dctime; INT ncstep;
  INT *i1ptyp;
  DBL *d1rjsig; DBL *d1rjtau; 
  INT *i1petyp; DBL *d1peem; DBL *d1peex; DBL *d1peey; INT *i1usan;
  INT *i1sbty; DBL *d1elfs;
  INT *i1pexc; INT *i1sbac;
  DBL *d1stkn; DBL *d1stkt;
  DBL *d1rjslpnor; DBL *d1rjdelnor;
  DBL *d1styns; DBL *d1strns;
  INT *i1rjstnor;
  DBL **d2elstr;
  DBL *d1stcoh; DBL *d1stfri;
  DBL **d2rjfrv; INT **i2rbedn;
{ INT i0,i1,i2,i3,i0r,i1r;
  INT i,imyelemj;
  DBL drx,dry,rx,ry,slp,slpnor,delta,syld,sang,slpnor0;
  DBL e1x,e1y,alfa,h,o,op,s,sp,st;  
  DBL dsigma,dtau=0,eps=0,veps=0,dnforce,dtforce,eps0; 
  DBL fy,epsy,epssh,epsu,sigmash,sigmasu,sigmabr,epsbr;
  DBL youngcon,dpemucon,dpelacon,dpeftcon,dpefscon;
  DBL k1e,k2e,k3e,k4e,k5e,k6e,eb,koef;
  DBL di0ri1, di0i1, di1ri2, di2i3;
  DBL k1s,k2s; 
  DBL pois;
  DBL cralfa;
  DBL reduc;
  DBL stif,sma,lc=0,di,mod=0;
  DBL small=0.000001,snorpl=0.0;
  DBL beta=1.0, ka=2.0;
  
  //DBL rebar_normal_penalty=210e9;
  //DBL rebar_tangential_penalty=210e9;
  
  DBL length; // average length of two rebars

  /*parameters in moment when crack open*/
  static DBL* d1slpin=DBL1NULL; //slip when crack open in tension
  static DBL* d1deltain=DBL1NULL; //delta when crack open in tension
  static DBL* d1slpnorpr=DBL1NULL;//previous normalized slip
  static DBL* d1deltapr=DBL1NULL;//previous delta
  static INT* indic=INT1NULL; //indicator if crack is open
  static DBL* d1sigmacr=DBL1NULL;//stress in steeel when crack open

  /*parameters for slip - strain relation */
  static DBL* d1snormax=DBL1NULL; //maximal slip

  /* parameters for stress - strain relation */
  static DBL* d1epspr=DBL1NULL; //previous strain 
  static DBL* d1vepspr=DBL1NULL; //previous velocity in bar direction
  static DBL* d1sigmapr=DBL1NULL; //previous stress
  static DBL* d1epsv=DBL1NULL;  //strain where velocity change sign
  static DBL* d1sigmav=DBL1NULL; //stress where velocity change sign
  static DBL* d1epss=DBL1NULL;  //strain where stress change sign
  static DBL* d1sigmas=DBL1NULL;  //stress where stress change sign
  static DBL* d1epsmax=DBL1NULL;  //maximum strain where velocity change sign
  static DBL* d1epssvmin=DBL1NULL; //strain where veloc. change sign when minimum stress
  static DBL* d1sigmavmin=DBL1NULL;  //minimum stress where veloc. change sign
  static DBL* d1sigmasmax=DBL1NULL;  //stress where stress change sign when max strain
  static DBL* d1epssmax=DBL1NULL;  //maximum strain where stress change sign
  static DBL* d1sigmak=DBL1NULL;  //loc. parametar for Fe5
  static DBL* d1epsk=DBL1NULL;  //loc. parametar for Fe5
  static FILE *contra=FILENULL; //brisi
  
  /* Variables to store opening, sliding, slip, and delta at the time of rebar activation */
  static DBL* d1oac=DBL1NULL; // opening when rebar is activated
  static DBL* d1sac=DBL1NULL; // sliding when rebar is activated
  static DBL* d1slpac=DBL1NULL; // half slip when rebar is activated
  static DBL* d1delac=DBL1NULL; // delta when rebar is activated
  static INT* icrac=INT1NULL; // flag if crack has an activated bar
  
  DBL normal_strain; // 1D joint strain in the direction of the rebar 
  DBL normal_displacement; // 1D joint displacement in the direction of the rebar 
  DBL yield_normal_strain; // 1D joint strain in the direction of the rebar at yielding
  DBL yield_normal_displacement; // 1D joint displacement in the direction of the rebar at yielding
  DBL yield_normal_stress; // 1D joint stress in the direction of the rebar at yielding
  DBL rupture_normal_strain ; // 1D joint strain in the direction of the rebar at rupture
  DBL rupture_normal_displacement ; // 1D joint displacement in the direction of the rebar at rupture
  DBL tangential_strain; // 1D joint strain in the direction perpendicular to the rebar
  DBL tangential_displacement; // 1D joint displacement in the direction perpendicular to the rebar
  DBL plastic_strain; // 1D joint plastic strain component in the direction of the rebar
  DBL plastic_displacement; // 1D joint plastic displacement component in the direction of the rebar
  DBL elastic_stress; //! 1D joint elastic stress component in the direction of the rebar
  DBL plastic_stress; //! 1D joint plastic stress component in the direction of the rebar
  
  static DBL* d1plstr=DBL1NULL; // plastic strain of 1D joint (used only if i1sbty[i1refbar[i0r]]==6)
  
  DBL rebar_sigman_el0; // Stress in the first triangle in the direction perpendicular to the average rebar direction
  DBL rebar_sigman_el1; // Stress in the second triangle in the direction perpendicular to the average rebar direction
  DBL rebar_sigman_average; // Average normal stress to the rebar
  DBL interface_shear_strength; // Shear strength of rebar-rock interface
  //DBL slip_normal_strain; // 1D joint strain in the direction of the rebar at slip 
  DBL interface_shear_strength_force; // Force resisting pullout
  DBL pullout_force; // Pullout force
  
  //printf("Yrb1RJOINT\n");
  
  /*static FILE *out2=FILENULL;
  if(out2 == FILENULL)
  { out2=fopen("Yrb1RJOINT.txt", "a");
  }*/
  
  /* function for relation slip-strain */
  #define FE1s(k1s,k2s,snor) (snor/k1s-(k2s/k1s))
  /* function for relation stress-strain */
  #define FE1(eps,sigmav,epsv,young) (sigmav+young*(eps-epsv))
  #define FE2(eps, fy, epsy,small) (fy+small*(eps-epsy))
  #define FE3(k1e,k2e,k3e,eps) (k1e*(eps*eps)+k2e*eps+k3e)
  #define FE4(fy,koef,eb,eps,epss,sigmas) (-fy*(koef-(koef*(koef-1.0))/((-eb/fy)*(eps-epss)+koef-1.0))+sigmas)
  #define FE5(fy,koef,eb,eps,epsk,sigmak) (sigmak+fy*(koef-(koef*(koef-1.0))/((-eb/fy)*(epsk-eps)+koef-1.0)))
  #define FE6(k4e,k5e,k6e,eps) (k4e*(eps*eps)+k5e*eps+k6e) 

  if(d1slpin==DBL1NULL)
  { d1slpin=TalDBL1(nbrjointrb);
    d1deltain=TalDBL1(nbrjointrb);
    d1slpnorpr=TalDBL1(nbrjointrb);
    d1deltapr=TalDBL1(nbrjointrb);
    indic=TalINT1(nbrjointrb);
    d1sigmacr=TalDBL1(nbrjointrb);
    for (i=0; i<nbrjointrb; i++)
    { d1slpin[i]=R0;
      d1deltain[i]=R0;
      d1slpnorpr[i]=R0;
      d1deltapr[i]=R0;
      indic[i]=0;
      d1sigmacr[i]=R0;
    }
  }
  if(d1vepspr==DBL1NULL)
  { d1snormax=TalDBL1(nbrjointrb);
    for(i=0; i<nbrjointrb; i++)
    { d1snormax[i]=R0;
    }
  }
  if(d1vepspr==DBL1NULL)
  { d1epspr=TalDBL1(nbrjointrb);
    d1vepspr=TalDBL1(nbrjointrb);
    d1sigmapr=TalDBL1(nbrjointrb);
    d1epsv = TalDBL1(nbrjointrb);
    d1sigmav=TalDBL1(nbrjointrb);
    d1epss=TalDBL1(nbrjointrb);
    d1sigmas=TalDBL1(nbrjointrb);
    d1epsmax=TalDBL1(nbrjointrb);
    d1epssvmin=TalDBL1(nbrjointrb);
    d1sigmavmin=TalDBL1(nbrjointrb);
    d1sigmasmax=TalDBL1(nbrjointrb);
    d1epssmax=TalDBL1(nbrjointrb);
    d1sigmak=TalDBL1(nbrjointrb);
    d1epsk=TalDBL1(nbrjointrb);
    for(i=0; i<nbrjointrb; i++)
    { d1vepspr[i]=R0;
      d1sigmapr[i]=R0;
      d1epsv[i]=R0;
      d1sigmav[i]=R0;
      d1epss[i]=R0;
      d1sigmas[i]=R0;
      d1epsmax[i]=R0;
      d1epssvmin[i]=R0;
      d1sigmavmin[i]=R0;
      d1sigmasmax[i]=R0;
      d1epssmax[i]=R0;
      d1sigmak[i]=R0;
      d1epsk[i]=R0;
    }
  }
  
  if(d1oac==DBL1NULL)
  { d1oac=TalDBL1(nbrjointrb);
    d1sac=TalDBL1(nbrjointrb);
    d1slpac=TalDBL1(nbrjointrb);
    d1delac=TalDBL1(nbrjointrb);
    icrac=TalINT1(nbrjointrb);
    for (i=0; i<nbrjointrb; i++)
    { d1oac[i]=R0;
      d1sac[i]=R0;
      d1slpac[i]=R0;
      d1delac[i]=R0;
      icrac[i]=0;
    }
  }
  
  if(d1plstr==DBL1NULL)
  { d1plstr=TalDBL1(nbrjointrb);
    for (i=0; i<nbrjointrb; i++)
    { d1plstr[i]=R0; } 
  }  
    
  //fprintf(out2,"Yrb1RJOINT \n");
  //! Main loop
  for(i=0;i2relto[0][i]>=0; i++) 
  { //fprintf(out2,"ncstep: %d \n",ncstep);
    //fprintf(out2,"i: %d \t i2relto[0][i]: %d \t i2relto[1][i]: %d \t iprop: %d \t i1myjoint[i]: %d\n",i,i2relto[0][i],i2relto[1][i],iprop,i1myjoint[i]);
    //fprintf(out2,"i1rprop[i2relto[0][i]]: %d \t i1rprop[i2relto[1][i]: %d \n\n",i1rprop[i2relto[0][i]],i1rprop[i2relto[1][i]]);
    if(i1rprop[i2relto[0][i]]==iprop && i1rprop[i2relto[1][i]]==iprop && i1myjoint[i]>=R0)
    { //fprintf(out2,"ncstep: %d \n",ncstep);
      //fprintf(out2,"i: %d \t i1myjoint[i]: %d \n",i,i1myjoint[i]);
      //printf("Yrb1RJOINT\n");
      //printf("ncstep: %d \n",ncstep);
      //printf("i: %d \t i1myjoint[i]: %d \n",i,i1myjoint[i]);
      imyelemj=i1myjoint[i]; //! ID of joint with the given 1D joint 
      //i0=i2elto[0][imyelemj]; //! ID of first node of joint 
      //i1=i2elto[1][imyelemj]; //! ID of second node of joint 
      //i2=i2elto[2][imyelemj]; //! ID of third node of joint 
      //i3=i2elto[3][imyelemj]; //! ID of fourth node of joint 
      i0r=i2relto[0][i]; //! ID of first reference point of 1D joint 
      i1r=i2relto[1][i]; //! ID of second reference point of 1D joint 
      //printf("i0: %d \t i1: %d \t i2: %d \t i3: %d\n",i0,i1,i2,i3);
      //printf("i0r: %d \t i1r: %d \n",i0r,i1r);
      //printf("i2elto[0][i1rmyel[i0r]]: %d\t i2elto[1][i1rmyel[i0r]]: %d\t i2elto[2][i1rmyel[i0r]]: %d\n",i2elto[0][i1rmyel[i0r]],i2elto[1][i1rmyel[i0r]],i2elto[2][i1rmyel[i0r]]);
      //printf("i2elto[0][i1rmyel[i1r]]: %d\t i2elto[1][i1rmyel[i1r]]: %d\t i2elto[2][i1rmyel[i1r]]: %d\n",i2elto[0][i1rmyel[i1r]],i2elto[1][i1rmyel[i1r]],i2elto[2][i1rmyel[i1r]]);
      drx=(d1rccgx[i1r]-d1rccgx[i0r]); // x-component of 1D joint in system of reference centered in i0r
      dry=(d1rccgy[i1r]-d1rccgy[i0r]); // y-component of 1D joint in system of reference centered in i0r
      
      //printf("i2rbedn[0][i0r]: %d\t i2rbedn[1][i0r]: %d\n",i2rbedn[0][i0r],i2rbedn[1][i0r]);
      //printf("i2rbedn[0][i1r]: %d\t i2rbedn[1][i1r]: %d\n\n",i2rbedn[0][i1r],i2rbedn[1][i1r]);
      /* Load node IDS using edges instead of joint elements */
      i0=i2rbedn[0][i0r]; //! ID of first node of first edge 
      i1=i2rbedn[1][i0r]; //! ID of second node of first edge
      i2=i2rbedn[0][i1r]; //! ID of first node of second edge
      i3=i2rbedn[1][i1r]; //! ID of second node of second edge
      
      
      //fprintf(out2,"d1rccgx[i0r]: %.6f \t d1rccgy[i0r]: %.6f \t d1rccgx[i1r]: %.6f \t d1rccgy[i1r]: %.6f \n",d1rccgx[i0r],d1rccgy[i0r],d1rccgx[i1r],d1rccgy[i1r]);
      //fprintf(out2,"drx: %.6f \t dry: %.6f \n",drx,dry);
      rx=((d1rccgx[i0r]-d1rsctrx[i0r])+(d1rsctrx[i1r]-d1rccgx[i1r]))/2.0; //! Average direction of two rebars
      ry=((d1rccgy[i0r]-d1rsctry[i0r])+(d1rsctry[i1r]-d1rccgy[i1r]))/2.0; //! Average direction of two rebars
      length=SQRT(rx*rx+ry*ry); //! Average length of two rebars
      rx=rx/(SQRT(rx*rx+ry*ry)+EPSILON); //! Normalized direction of rebar 
      ry=ry/(SQRT(rx*rx+ry*ry)+EPSILON); //! Normalized direction of rebar
      //fprintf(out2,"rx: %.6f \t ry: %.6f \n",rx,ry);
      di0ri1=SQRT((d1rccgy[i0r]-d1nccy[i1])*(d1rccgy[i0r]-d1nccy[i1])+
                  (d1rccgx[i0r]-d1nccx[i1])*(d1rccgx[i0r]-d1nccx[i1])); //! Distance between i0r and i1
      di0i1=SQRT((d1nccy[i0]-d1nccy[i1])*(d1nccy[i0]-d1nccy[i1])+
                 (d1nccx[i0]-d1nccx[i1])*(d1nccx[i0]-d1nccx[i1])); //! Distance between i0 and i1 
      di1ri2=SQRT((d1rccgy[i1r]-d1nccy[i2])*(d1rccgy[i1r]-d1nccy[i2])+
                  (d1rccgx[i1r]-d1nccx[i2])*(d1rccgx[i1r]-d1nccx[i2])); //! Distance between i1r and i2
      di2i3=SQRT((d1nccy[i2]-d1nccy[i3])*(d1nccy[i2]-d1nccy[i3])+
                 (d1nccx[i2]-d1nccx[i3])*(d1nccx[i2]-d1nccx[i3])); //! Distance between i2 and i3
      //fprintf(out2,"di0ri1: %.6f \t di0i1: %.6f \t di1ri2: %.6f \t di2i3: %.6f\n",di0ri1,di0i1,di1ri2,di2i3);
      slp=(drx*rx+dry*ry)/2.0; //! Half slip = 0.5 * (component of 1D joint along average bar direction)    
      delta=(-drx*ry+dry*rx)/2.0; //! Component of 1D joint in the direction perpendicular to the average bar direction
      /* Normalize slip by rebar diameter and multiply by Kfc = (UCS_concrete (in MPa))/(20)^(2/3) */
      slpnor=(slp/d1sdiam[i1refbar[i0r]])*pow((d1mpsfc[i1sbpr[i1refbar[i0r]]]/20.0),(2.0/3.0));
      //fprintf(out2,"slp: %.6f \t delta: %.6f \t slpnor: %.6f \t d1sdiam: %.6f \t d1mpsfc: %.6f\n",slp,delta,slpnor,d1sdiam[i1refbar[i0r]],d1mpsfc[i1sbpr[i1refbar[i0r]]]);      
      e1x=RP5*(d1nccx[i1]+d1nccx[i2]-d1nccx[i0]-d1nccx[i3]); //! x-component of joint axis 
      e1y=RP5*(d1nccy[i1]+d1nccy[i2]-d1nccy[i0]-d1nccy[i3]); //! y-component of joint axis
      h=SQRT(e1x*e1x+e1y*e1y); //! Length of joint axis
      e1x=e1x/(h+EPSILON); //! x-component of normalized joint axis
      e1y=e1y/(h+EPSILON); //! y-component of normalized joint axis
      //fprintf(out2,"e1x: %.6f \t e1y: %.6f \n",e1x,e1y); 
      if((e1y*rx-e1x*ry)>1)
      { alfa = 0.0; }
      else if ((e1y*rx-e1x*ry)<-1)
      { alfa = MYPI; }
      else
      { alfa=acos(e1y*rx-e1x*ry); } // Angle (in radians) betwen average direction of rebar and crack normal
      //fprintf(out2,"alfa: %.6f \n",alfa);
      o=(d1rccgy[i0r]-d1rccgy[i1r])*e1x-(d1rccgx[i0r]-d1rccgx[i1r])*e1y; //! Opening between two reference points of 1D joint in the joint reference system 
      s=(d1rccgy[i0r]-d1rccgy[i1r])*e1y+(d1rccgx[i0r]-d1rccgx[i1r])*e1x; //! Sliding  between two reference points of 1D joint in the joint reference system  
      //fprintf(out2,"o: %.6f \t s: %.6f \n",o,s); 
      
      /* Getting Young's modulus of element with i0r */
      if (i1petyp[i1elpr[i1rmyel[i0r]]]==YTE2TRIELS) /* Elastic element with Lame's parameters */
      { dpemucon=d1pemu[i1elpr[i1rmyel[i0r]]];
        dpelacon=d1pela[i1elpr[i1rmyel[i0r]]];
        youngcon=dpemucon*(3*dpelacon+2*dpemucon)/(dpelacon+dpemucon);
        pois=(dpelacon/(2.0*(dpelacon+dpemucon))); //! Poisson's ratio of "concrete"
        youngcon=youngcon/(1.0-pois*pois); //! Young's modulus of "concrete"
      }
      else if (i1petyp[i1elpr[i1rmyel[i0r]]]==YTE2PLANESTRESS) /* Plane stress elastic element */
      { youngcon=d1peem[i1elpr[i1rmyel[i0r]]]; }
      else /* Plane stress elastic element */
      { youngcon=d1peem[i1elpr[i1rmyel[i0r]]]; }
      
      if(i1usan[i1elpr[i1rmyel[i0r]]]==1)
      { youngcon=0.5*(d1peex[i1elpr[i1rmyel[i0r]]]+d1peey[i1elpr[i1rmyel[i0r]]]); }/* Transversely isotropic elastic element */
      
      //dpemucon=d1pemu[i1elpr[i1rmyel[i0r]]]; //! Second Lame's parameter of triangle with i0r 
      //dpelacon=d1pela[i1elpr[i1rmyel[i0r]]]; //! First Lame's parameter of triangle with i0r
      //fprintf(out2,"d1pemu: %.6f \t dpelacon: %.6f \n",dpelacon,dpelacon);
      
      if(i1elpr[imyelemj]>=0) //! Check if joint is not broken
      { dpeftcon=d1peftc[i1elpr[imyelemj]-nprop]; //! Tensile strength of "concrete"
        //dpefscon=d1pefsc[i1elpr[imyelemj]-nprop]; //! Shear strength of "concrete" 
      }
      else
      { dpeftcon=d1peftc[i1elpr[imyelemj]+YIPROPMAX-nprop]; //! Tensile strength of "concrete"
        //dpefscon=d1pefsc[i1elpr[imyelemj]+YIPROPMAX-nprop]; //! Shear strength of "concrete" 
      }
      
      dpefscon=d1elfs[imyelemj]; //! Shear strength of "concrete" (it varies element by element due to Mohr-Coulomb)
      //fprintf(out2,"dpeftcon: %.6f \t dpefscon: %.6f \n",dpeftcon,dpefscon);
      //printf("dpeftcon: %.6f \t dpefscon: %.6f \n",dpeftcon,dpefscon);
      //youngcon=dpemucon*(3*dpelacon+2*dpemucon)/(dpelacon+dpemucon);
      //pois=(dpelacon/(2.0*(dpelacon+dpemucon))); //! Poisson's ratio of "concrete"
      //youngcon=youngcon/(1.0-pois*pois); //! Young's modulus of "concrete"
      //fprintf(out2,"youngcon: %.6f \t pois: %.6f \n",youngcon,pois);
      //fprintf(out2,"youngcon: %.6f \n",youngcon);
      
      //fprintf(out2,"i1elpr[imyelemj]: %d \n",i1elpr[imyelemj]);
      if(i1elpr[imyelemj]>=0) //! Check if joint is not broken
      { //fprintf(out2,"Not broken\n");
        op=R2*h*d1peftc[i1elpr[imyelemj]-nprop]/d1pepec[i1elpr[imyelemj]-nprop]; //! Peak opening of non-broken joint
        //sp=R2*h*d1pefsc[i1elpr[imyelemj]-nprop]/d1pepec[i1elpr[imyelemj]-nprop]; //! Peak sliding of non-broken joint
        //printf("sp: %.6f \t d1pefsc[i1elpr[imyelemj]-nprop]: %.6f \n",sp,d1pefsc[i1elpr[imyelemj]-nprop]);
        sp=R2*h*dpefscon/d1pepec[i1elpr[imyelemj]-nprop]; //! Peak sliding of non-broken joint
        //printf("sp: %.6f \t dpefscon: %.6f \n",sp,dpefscon);
        if(i1ptyp[i1elpr[imyelemj]-nprop]==YTE2JOINTS) //! Notice that "-nprop" was missing in the original version! wtf!
        { //st=MAXIM(EPSILON,(R3*d1pegfs[i1elpr[imyelemj]-nprop]/d1pefsc[i1elpr[imyelemj]-nprop])); //! Residual sliding of joint
          //printf("st: %.6f \t d1pefsc[i1elpr[imyelemj]-nprop]: %.6f \n",st,d1pefsc[i1elpr[imyelemj]-nprop]);
          st=MAXIM(EPSILON,(R3*d1pegfs[i1elpr[imyelemj]-nprop]/dpefscon)); //! Residual sliding of joint
          //printf("st: %.6f \t dpefscon: %.6f \n",st,dpefscon);
          //fprintf(out2,"st: %.6f \n",R3*d1pegfs[i1elpr[imyelemj]-nprop]/d1pefsc[i1elpr[imyelemj]-nprop]);
          //fprintf(out2,"st: %.6f \n",R3*d1pegfs[i1elpr[imyelemj]-nprop]/d1pefsc[i1elpr[imyelemj]-nprop]);  
	}
        else if (i1ptyp[i1elpr[imyelemj]-nprop]==YTE2JOINTSCON) //! Notice that "-nprop" was missing in the original version! wtf!
        { st=MAXIM((R2*sp),(R2*d1pegfs[i1elpr[imyelemj]-nprop]/d1pefsc[i1elpr[imyelemj]-nprop])); }
      }
      else 
      { //fprintf(out2,"Broken\n");
        op=R2*h*d1peftc[i1elpr[imyelemj]+YIPROPMAX-nprop]/d1pepec[i1elpr[imyelemj]+YIPROPMAX-nprop]; //! Peak opening of broken joint
        sp=R2*h*d1pefsc[i1elpr[imyelemj]+YIPROPMAX-nprop]/d1pepec[i1elpr[imyelemj]+YIPROPMAX-nprop]; //! Peak sliding of broken joint
        sp=R2*h*dpefscon/d1pepec[i1elpr[imyelemj]+YIPROPMAX-nprop]; //! Peak sliding of broken joint
        if(i1ptyp[i1elpr[imyelemj]+YIPROPMAX-nprop]==YTE2JOINTS) //! Notice that "-nprop" was missing in the original version! wtf!
        { st=MAXIM((EPSILON),(R3*d1pegfs[i1elpr[imyelemj]+YIPROPMAX-nprop]/d1pefsc[i1elpr[imyelemj]+YIPROPMAX-nprop]));
          st=MAXIM((EPSILON),(R3*d1pegfs[i1elpr[imyelemj]+YIPROPMAX-nprop]/d1pefsc[i1elpr[imyelemj]+YIPROPMAX-nprop]));
        }
        else if (i1ptyp[i1elpr[imyelemj]+YIPROPMAX-nprop]==YTE2JOINTSCON) //! Notice that "-nprop" was missing in the original version! wtf!
        { st=MAXIM((R2*sp),(R2*d1pegfs[i1elpr[imyelemj]+YIPROPMAX-nprop]/dpefscon));
        } 
      }
      //printf("i0r: %d \t i1refbar[i0r]: %d \t i1sbty[i1refbar[i0r]]: %d \n",i0r,i1refbar[i0r],i1sbty[i1refbar[i0r]]);
      
//       fprintf(out2,"d1pepec: %.6f \n",d1pepec[i1elpr[imyelemj]-nprop]);
//       fprintf(out2,"d1pegfs: %.6f \t d1pefsc: %.6f\n",d1pegfs[i1elpr[imyelemj]-nprop],d1pefsc[i1elpr[imyelemj]-nprop]);
      //fprintf(out2,"op: %.6f \t sp: %.6f \t st: %.6f \n",op,sp,st);
      
      //! Check if joint opened 
      if((o>op||ABS(s)>sp) && indic[i]==R0) //! Joint opens for the first time (i.e., joint @ 1D joint reference point is softening)
      { d1slpin[i]=d1slpnorpr[i]; //! Get normalized slip from previous time step (index i refers to  1D joint ID)
        d1deltain[i]=d1deltapr[i]; //! Get delta from previous time step (index i refers to  1D joint ID)
        d1sigmacr[i]=d1sigmapr[i]; //! Get stress in 1D joint element from previous time step  (index i refers to  1D joint ID)                                                                          
        indic[i]=1; //! Set the flag for that 1D joint to opened crack
      }
      if(indic[i]==R0) //! Joint is still elastic 
      { d1slpnorpr[i]=slpnor; //! Store normalized slip
        d1deltapr[i]=delta; //! Store delta 
      }
      //fprintf(out2,"young: %.6f \n\n\n",young);
      
      /* 1D joint formulation wih normal and tangential stiffnesses, infinitely elastic */
      if(i1sbty[i1refbar[i0r]]==5)
      { 
        //printf("1D joint type 5\n");
        if((i1sbac[i1refbar[i0r]]>0) && icrac[i]==0) //! If 1D joint is activate for the first time  
        { d1oac[i]=o; //! Store opening at the time of rebar activation
          d1sac[i]=s; //! Store sliding at the time of rebar activation
          d1slpac[i]=slp; //! Store half slip at the time of rebar activation 
          d1delac[i]=delta; //! Store delta at the time of rebar activation
          icrac[i]=1; //! Set the flag for this 1D joint element to be activated
          //printf("Activated joint element\n");
          //printf("i: %d \t d1oac[i]: %.10f \t d1sac[i]: %.10f \t d1slpac[i]: %.10f \t d1delac[i]: %.10f \n",i,d1oac[i],d1sac[i],d1slpac[i],d1delac[i]);
        }
        
        //printf("i: %d \t o: %.10f \t s: %.10f \t slp: %.10f \t delta[i]: %.10f \n",i,o,s,slp,delta);
        
        o=o-d1oac[i]; //! Current opening - opening at the time of rebar activation
        s=s-d1sac[i]; //! Current sliding - sliding at the time of rebar activation
        slp=slp-d1slpac[i]; //! Current half slip - half slip at the time of rebar activation
        delta=delta-d1delac[i]; //! Current delta - delta at the time of rebar activation
        
        //printf("i: %d \t o: %.10f \t s: %.10f \t slp: %.10f \t delta[i]: %.10f \n",i,o,s,slp,delta);
        
        /* Compute normal stress: normal_stiffness * slip */  
        dsigma=d1stkn[i1sbpr[i1refbar[i0r]]]*slp; //! Units of d1stkn are those of a stress / length
        //normal_strain = slp/(length+EPSILON);
        //dsigma=d1stkn[i1sbpr[i1refbar[i0r]]]*normal_strain; //! Units of d1stkn are those of a stress 
        /* Compute tangential stress: tangential_stiffness * delta */
        dtau=d1stkt[i1sbpr[i1refbar[i0r]]]*delta; //! Units of d1stkt are those of a stress / length
        //tangential_strain = delta/(length+EPSILON);
        //dtau=d1stkt[i1sbpr[i1refbar[i0r]]]*tangential_strain; //! Units of d1stkt are those of a stress
       
        //printf("i: %d \t d1stkn: %f \t d1stkt: %f \t disgma: %.10f \t dtau: %.10f \n",i,d1stkn[i1sbpr[i1refbar[i0r]]],d1stkt[i1sbpr[i1refbar[i0r]]],dsigma,dtau);
      }
       /* 1D joint formulation wih normal and tangential stiffnesses (in [stress/length], elastic - perfectly plastic  */
      else if(i1sbty[i1refbar[i0r]]==6)
      { 
        //printf("1D joint type 6\n");
        if((i1sbac[i1refbar[i0r]]>0) && icrac[i]==0) //! If 1D joint is activated for the first time  
        { d1oac[i]=o; //! Store opening at the time of rebar activation
          d1sac[i]=s; //! Store sliding at the time of rebar activation
          d1slpac[i]=slp; //! Store half slip at the time of rebar activation 
          d1delac[i]=delta; //! Store delta at the time of rebar activation
          icrac[i]=1; //! Set the flag for this 1D joint element to be activated
          d1plstr[i]=0.0; //! Set the plastic displacement of the 1D joint to zero (axial direction)
        }
        
        o=o-d1oac[i]; //! Current opening - opening at the time of rebar activation
        s=s-d1sac[i]; //! Current sliding - sliding at the time of rebar activation
        slp=slp-d1slpac[i]; //! Current half slip - half slip at the time of rebar activation
        delta=delta-d1delac[i]; //! Current delta - delta at the time of rebar activation
        
        /* Compute normal stress (elastic - perfectly plastic */  
        
        yield_normal_stress = d1styns[i1sbpr[i1refbar[i0r]]]; //! Normal stress of 1D joint at yielding point
        rupture_normal_displacement = d1strns[i1sbpr[i1refbar[i0r]]]; //! Normal displacement of 1D joint at rupture
        normal_displacement = slp; //! Displacement in the direction of the rebar (extension)
        yield_normal_displacement = yield_normal_stress/d1stkn[i1sbpr[i1refbar[i0r]]]; //! Normal displacement of 1D joint at yielding
        plastic_displacement = d1plstr[i];
        
        if((normal_displacement < (-rupture_normal_displacement)) || (normal_displacement > (rupture_normal_displacement))) //! broken
        { dsigma=0.0;
          i1rprop[i2relto[0][i]]=iprop-YIPROPMAX; //! Remove 1D joint
          i1rprop[i2relto[1][i]]=iprop-YIPROPMAX; //! Remove 1D joint
        }
        else if (normal_displacement >= plastic_displacement) //! Tension
        { elastic_stress = (normal_displacement - plastic_displacement) * d1stkn[i1sbpr[i1refbar[i0r]]];
          plastic_stress = yield_normal_stress;
          if(elastic_stress<plastic_stress)
          { dsigma = elastic_stress; 
            if(i1rjstnor[i]>1) //! 1D joint state = yielded in the past
            { i1rjstnor[i]=3;
            }
          }
          else
          { dsigma = plastic_stress; 
            i1rjstnor[i] = 2; //! 1D joint state = yielded
          }
        }
        else //! Compression
        { elastic_stress = (normal_displacement - plastic_displacement) * d1stkn[i1sbpr[i1refbar[i0r]]];
          plastic_stress = -yield_normal_stress;
          if(elastic_stress>plastic_stress)
          { dsigma = elastic_stress; 
            if(i1rjstnor[i]>1) //! 1D joint state = yielded in the past
            { i1rjstnor[i]=3;
            }
          }
          else
          { dsigma = plastic_stress; 
            i1rjstnor[i] = 2; //! 1D joint state = yielded
          }
        }
        
        /* Update plastic_displacement of 1D joint */
        d1plstr[i] = normal_displacement - dsigma / d1stkn[i1sbpr[i1refbar[i0r]]];
        
        //fprintf(out2,"dctime: %.6f \t normal_strain: %.6f \t dsigma: %.6f \n",dctime,normal_strain,dsigma);
        //fprintf(out2,"%.6f \t %.6f \t %.6f \n",dctime,normal_displacement,dsigma); 

        /* old formulation for monotonic loading only */
        /*if(ABS(normal_strain)<=yield_normal_strain) //! elastic 
        { dsigma=d1stkn[i1sbpr[i1refbar[i0r]]]*normal_strain; 
          i1rjstnor[i]=1; //! 1D joint state = elastic
        }  
        else if((ABS(normal_strain)>yield_normal_strain)&&(ABS(normal_strain)<=rupture_normal_strain)) //! perfectly-plastic 
        { dsigma=(slp/(ABS(slp)+EPSILON))*yield_normal_stress;
          i1rjstnor[i]=2; //! 1D joint state = yielded
        } 
        else //! broken
        { dsigma=0.0;
          i1rprop[i2relto[0][i]]=iprop-YIPROPMAX; //! Remove 1D joint
          i1rprop[i2relto[1][i]]=iprop-YIPROPMAX; //! Remove 1D joint
        }*/ 
        
        /* Compute tangential stress (elastic without interaction with the normal stress */
        tangential_displacement = delta;
        dtau=d1stkt[i1sbpr[i1refbar[i0r]]]*tangential_displacement; //! Units of d1stkt are those of a stress/length
       
      }
      /* 1D joint formulation wih normal and tangential stiffnesses */
      /* + maximum tensile stress limited by frictional resistance of rebar-rock interface */
      else if(i1sbty[i1refbar[i0r]]==7)
      { 
        //printf("1D joint type 7\n");
        if((i1sbac[i1refbar[i0r]]>0) && icrac[i]==0) //! If 1D joint is activated for the first time  
        { d1oac[i]=o; //! Store opening at the time of rebar activation
          d1sac[i]=s; //! Store sliding at the time of rebar activation
          d1slpac[i]=slp; //! Store half slip at the time of rebar activation 
          d1delac[i]=delta; //! Store delta at the time of rebar activation
          icrac[i]=1; //! Set the flag for this 1D joint element to be activated
        }
        
        o=o-d1oac[i]; //! Current opening - opening at the time of rebar activation
        s=s-d1sac[i]; //! Current sliding - sliding at the time of rebar activation
        slp=slp-d1slpac[i]; //! Current half slip - half slip at the time of rebar activation
        delta=delta-d1delac[i]; //! Current delta - delta at the time of rebar activation
        
        /* Compute normal (axial) stress */  
        
        yield_normal_stress = d1styns[i1sbpr[i1refbar[i0r]]];
        rupture_normal_strain = d1strns[i1sbpr[i1refbar[i0r]]];
        normal_strain = slp/(length+EPSILON); //! Strain in the direction of the rebar
        yield_normal_strain = yield_normal_stress/d1stkn[i1sbpr[i1refbar[i0r]]]; //! Normal strain of 1D joint at yielding
        
        /* Stress normal to the average rebar direction in first and second triangular element */
        rebar_sigman_el0=d2elstr[0][i1rmyel[i0r]]*ry*ry+d2elstr[3][i1rmyel[i0r]]*rx*rx-2*d2elstr[1][i1rmyel[i0r]]*rx*ry;
        rebar_sigman_el1=d2elstr[0][i1rmyel[i1r]]*ry*ry+d2elstr[3][i1rmyel[i1r]]*rx*rx-2*d2elstr[1][i1rmyel[i1r]]*rx*ry;
        rebar_sigman_average = (rebar_sigman_el0+rebar_sigman_el1)/2; //! Average normal stress
        
        if(rebar_sigman_average>=0) //! Tension, cohesion only
        { interface_shear_strength=d1stcoh[i1sbpr[i1refbar[i0r]]]; } 
        else //! Compression, Mohr-Coulomb shear strength criterion
        { interface_shear_strength=d1stcoh[i1sbpr[i1refbar[i0r]]]-d1stfri[i1sbpr[i1refbar[i0r]]]*rebar_sigman_average; } 
        
        interface_shear_strength_force=interface_shear_strength*length*MYPI*d1sdiam[i1refbar[i0r]]; //! Force resisting pullout
                
        //slip_normal_strain = interface_shear_strength/d1stkn[i1sbpr[i1refbar[i0r]]]; //! Normal strain of 1D joint at slip

        if(normal_strain >= 0) /* tension (ie rebar is being pulled out */
        { dsigma=d1stkn[i1sbpr[i1refbar[i0r]]]*normal_strain;
          i1rjstnor[i]=1; //! 1D joint state = elastic
          pullout_force=dsigma*d1spea[i1refbar[i0r]]; //! Pullout force
          if(pullout_force>interface_shear_strength_force) //! Rebar is slipping
          { dsigma=interface_shear_strength_force/d1spea[i1refbar[i0r]]; //! Equivalent axial stress
            i1rjstnor[i]=4; //! 1D joint state = slipped
          }
        }
        else /* compression */
        {
          dsigma=d1stkn[i1sbpr[i1refbar[i0r]]]*normal_strain;
          i1rjstnor[i]=1; //! 1D joint state = elastic
        }

        /* Compute tangential stress (elastic without interaction with the normal stress */
        tangential_strain = delta/(length+EPSILON);
        dtau=d1stkt[i1sbpr[i1refbar[i0r]]]*tangential_strain; //! Units of d1stkt are those of a stress
        
        //printf("dsigma: %.6f\n",dsigma);
	//printf("interface_shear_strength: %.6f\n",interface_shear_strength);
        //printf("interface_shear_strength_force: %.6f\n",interface_shear_strength_force);
	//printf("i1rjstnor: %d\n",i1rjstnor[i]);
      }
      /* 1D joint formulation from Y-RC (i1sbty[i1refbar[i0r]] == 0 or 1)*/
      else
      { /* Compute stress in 1D joint when the joint state is elastic:
        dsigma = (E_steel/E_concrete)*(sigma_joint +/- tau_joint*(sin alpha / cos alfa) 
        (see eq. (14) in Zivaljic et al 2012) */
        //fprintf(out2,"ncstep: %d \t i: %d \t dsigma: %.6f \t young: %.6f \t youngcon: %.6f \t o: %.6f \t dpeftcon: %.6f \t op: %.6f \n",ncstep,i, dsigma,young,youngcon,o,dpeftcon,op);
        //fprintf(out2,"ABS(s): %.6f \t sp: %.6f \t dpefscon: %.6f \t e1x: %.6f \t rx: %.6f \t e1y: %.6f \t ry: %.6f \t alfa: %.6f \t cos(alfa): %.6f\n",ABS(s),sp, dpefscon,e1x,rx,e1y,ry,alfa,cos(alfa));
        if(indic[i]==R0) 
        { if (s>R0 && o>R0) //! Tensile state with sliding greater than zero
          { dsigma=(young/youngcon)*((R2*o/op-(o/op)*(o/op))*dpeftcon-
                   (R2*ABS(s)/sp-(ABS(s)/sp)*(ABS(s)/sp))*dpefscon*(e1x*rx+e1y*ry)/(cos(alfa)+EPSILON));
          /*if(ncstep>58600 && ncstep<58650)
	  { fprintf(out2,"ncstep: %d \t i: %d \t dsigma: %.6f \t 1 \n",ncstep,i, dsigma); }*/
          }
          else if (s<=R0 && o>R0) //! Tensile state with sliding less than zero
          { dsigma=(young/youngcon)*((R2*o/op-(o/op)*(o/op))*dpeftcon+
                   (R2*ABS(s)/sp-(ABS(s)/sp)*(ABS(s)/sp))*dpefscon*(e1x*rx+e1y*ry)/(cos(alfa)+EPSILON));
          /*if(ncstep>58000 && ncstep<58650)
	  {fprintf(out2,"ncstep: %d \t i: %d \t dsigma: %.6f \t 2 \n",ncstep,i, dsigma); }*/
          }
          else if(s>R0 && o<=R0) //! Compression state with sliding greater than zero (TODO: brackets are missing?)
          { dsigma=(young/youngcon)*(R2*o*dpeftcon/op)-
                   (R2*ABS(s)/sp-(ABS(s)/sp)*(ABS(s)/sp))*dpefscon*(e1x*rx+e1y*ry)/(cos(alfa)+EPSILON);
            /*if(ncstep>58000 && ncstep<58650)
	    { fprintf(out2,"ncstep: %d \t i: %d \t dsigma: %.6f \t 3 \n",ncstep,i, dsigma); }*/
          }
          else //! Compression state with sliding less than than zero (TODO: brackets are missing?)
          { dsigma=(young/youngcon)*(R2*o*dpeftcon/op)+
                   (R2*ABS(s)/sp-(ABS(s)/sp)*(ABS(s)/sp))*dpefscon*(e1x*rx+e1y*ry)/(cos(alfa)+EPSILON);
            /*if(ncstep>58000 && ncstep<58650)
	    { fprintf(out2,"ncstep: %d \t i: %d \t dsigma: %.6f \t 4 \n",ncstep,i, dsigma); }*/
          }
        }
        /* Compute stress in 1D joint when the joint state is either softened or broken */
        else 
        { /* The bar is elastic-plastic-brittle */
          if(i1sbty[i1refbar[i0r]]==1)
          { slpnor=slpnor-d1slpin[i]; //! Current normalized slip - normalized slip when joint @ 1D joint reference point is softening
            delta=delta-d1deltain[i]; //! Current Delta - delta when joint @ 1D joint reference point is softening
            slp=slpnor/pow((d1mpsfc[i1sbpr[i1refbar[i0r]]]/20.0),(2.0/3.0))*d1sdiam[i1refbar[i0r]]; //! Un-normalized slpnor
            fy=d1sfy[i1sbpr[i1refbar[i0r]]]; // Yield stress of steel bar

            /* Determine shear stress */
            sma=MYPI*pow(d1sdiam[i1refbar[i0r]],4.0)/64.0; // Area moment of inertia of bar: I_2 = (pi/4.0) * r^4
            //! Foundation stiffness for concrete: k = 150 * UCS_concrete^(0.85)* D_bar (eq. 31)
            stif=150.0*pow((d1mpsfc[i1sbpr[i1refbar[i0r]]]),0.85)/d1smmdiam[i1refbar[i0r]]*
                 (d1sfc[i1sbpr[i1refbar[i0r]]]/d1mpsfc[i1sbpr[i1refbar[i0r]]])*
                 (d1smmdiam[i1refbar[i0r]]/d1sdiam[i1refbar[i0r]]);
            //! Compute the length of the curvature-influenzing zone
            lc=(3.0*MYPI)/4.0*pow((4.0*young*sma/(stif*d1sdiam[i1refbar[i0r]])),0.25); //! Value derived from modelling the bar as beam resting on an elastic foundation
            //! Empirical non-dimensional parameter (Eq. 29)
            di=(1+150*slp/d1sdiam[i1refbar[i0r]])*(delta/d1sdiam[i1refbar[i0r]]);
            if(di>0.02)
            { lc=lc*(1.0+3.0*pow((di-0.02),0.8)); //! Eq. (28)
            }
            dtau=((young*sma*4.0)/(d1sdiam[i1refbar[i0r]]*d1sdiam[i1refbar[i0r]]*MYPI))*
                 ((34.9091*delta)/(lc*lc*lc)); //! Shear stress in the  1D joint element (eq. 34)
            mod=SQRT(MAXIM(EPSILON,1.0-3.0*pow((dtau/fy),2))); //! Reducing factor for the yield stress, fy (eq. 35, Von Mises criterion + isotropic hardenind)

            /* Determine critical strain and stress values of steel */
            epsy=(fy/young)*mod;  // Yield strain of bar  
            epssh=d1epssh[i1sbpr[i1refbar[i0r]]]; // Hardening strain of bar (from input) 
            epsu=d1epsu[i1sbpr[i1refbar[i0r]]];  // Ultimate strain of bar (from input)
            epsbr=d1epsbr[i1sbpr[i1refbar[i0r]]];  // Breakage strain of bar (from input)
            fy=d1sfy[i1sbpr[i1refbar[i0r]]]*mod; // Yield stress of bar (reduced by shear stress factor) 
            sigmash=(fy+(epssh-epsy)*small); // Hardening stress of bar (small hardening slope)
            sigmasu=d1sfu[i1sbpr[i1refbar[i0r]]]*mod;  // Utimate stress of bar (from input) reduced by shear stress factor
            sigmabr=d1sfbr[i1sbpr[i1refbar[i0r]]]*mod;  // Breakage stress of bar (from input) reduced by shear stress factor 

            /* Compute reduction factor to account for adjacent cracks, eq. 38 */ 
            /*cralfa=1.0-exp(-pow((0.065*d1crlcr[i1refbar[i0r]]/d1sdiam[i1refbar[i0r]]+0.5),3)); 
            if(cralfa>=0.087*d1crlcr[i1refbar[i0r]]/d1sdiam[i1refbar[i0r]])
            { cralfa=0.087*d1crlcr[i1refbar[i0r]]/d1sdiam[i1refbar[i0r]];
            }*/
            cralfa=1.0;
 
            /* Calculation of critical adimensional slip values based on the relation by Soltani and Maekawa (2008)*/
            sang=(epsy*(ka+3500.0*epsy))*cralfa; //! If epsilon_y < epsilon < epsilon_sh (eq. 19) (ka = 2.0)                                                                             
            syld=sang-sang*d1sigmacr[i]/(epsy*young); //! Yielding slip ?                                                             
            eps0=epssh+(0.06-syld/2.0)/(0.013*(sigmasu-fy)* 
              (d1mpsfc[i1sbpr[i1refbar[i0r]]]/d1sfc[i1sbpr[i1refbar[i0r]]])); //! (eq. 22)
            slpnor0=(0.047*(sigmasu-fy)*(d1mpsfc[i1sbpr[i1refbar[i0r]]]/d1sfc[i1sbpr[i1refbar[i0r]]])*(eps0-epssh))*cralfa+syld; //! If epsilon_sh < epsilon < epsilon_0 (eq. 21) 
       
            /* Calculation of strain from slip */
       
            /* Compute maximum strain values (differential normalized slip is greater than any previous value (i.e., loading conditions)) 
               Piecewise bi-linear relationship */
            if(slpnor>d1snormax[i])  
            { d1snormax[i]=slpnor; //! Store differential normalized slip value 
              if(d1snormax[i]<=syld) //! Differential normalized slip value is <= yielding slip  
              { k1s=sang/epsy; 
                k2s=-k1s*d1sigmacr[i]/young; //! d1sigmacr is the stress in the bar at the time step previous to softening (i.e., when crack opens)
                d1epsmax[i]=FE1s(k1s,k2s,d1snormax[i]); //! Compute and store strain (highest value reached so far)
              }
              else if (d1snormax[i]<=syld+small) //! Differential normalized slip value is <= yielding slip + small
              { k1s=(syld+small)/(epssh-epsy);
                k2s=-k1s*epsy+syld;
                d1epsmax[i]=FE1s(k1s,k2s,d1snormax[i]); //! Compute and store strain (highest value reached so far)
              }
              else if (d1snormax[i]<=slpnor0) //! Differential normalized slip value is <= slip_0 value (See Fig. 9)
              { k1s=(0.047*(sigmasu-fy)*(d1mpsfc[i1sbpr[i1refbar[i0r]]]/
                    d1sfc[i1sbpr[i1refbar[i0r]]]))*cralfa;
                k2s=(-0.047*(sigmasu-fy)*(d1mpsfc[i1sbpr[i1refbar[i0r]]]/
                    d1sfc[i1sbpr[i1refbar[i0r]]])*epssh)*cralfa+syld+small;
                d1epsmax[i]=FE1s(k1s,k2s,d1snormax[i]); //! Compute and store strain (highest value reached so far)
              }
              else //!  Differential normalized slip value is > slip_0 value (See Fig. 9)
              { k1s=(0.007*(sigmasu-fy)*(d1mpsfc[i1sbpr[i1refbar[i0r]]]/
                    d1sfc[i1sbpr[i1refbar[i0r]]]))*cralfa;
                k2s=(-0.007*(sigmasu-fy)*(d1mpsfc[i1sbpr[i1refbar[i0r]]]/
                    d1sfc[i1sbpr[i1refbar[i0r]]])*eps0)*cralfa+slpnor0;
                d1epsmax[i]=FE1s(k1s,k2s,d1snormax[i]); //! Compute and store strain (highest value reached so far)
              }
            } 
        
            /* Compute strain value */ 
            if(d1snormax[i]<=syld ) //! Bar has not yielded yet (elastic behaviour)                                                       
            {  k1s=sang/epsy;                                                                              
               k2s=-k1s*d1sigmacr[i]/young;
               eps=FE1s(k1s,k2s,slpnor); //! Compute strain value
            }
            else // Bar has yielded (see eq. 24 and 25, beta = 1.0) 
            { k1s=(epssh*sang+d1epsmax[i]*sang+epsy*(d1snormax[i]-syld)*(1.0+beta))/
                  (epsy*(d1epsmax[i]+epssh));
              k2s=(-d1sigmacr[i]*sang*(epssh+d1epsmax[i])-young*(epssh*sang*
                  (d1epsmax[i]-epsy)+d1epsmax[i]*d1epsmax[i]*sang+d1epsmax[i]*epsy*
                  (beta*(d1snormax[i]-syld)-sang)+epssh*epsy*
                  (syld-d1snormax[i])))/(epsy*young*(d1epsmax[i]+epssh)); 
              eps=FE1s(k1s,k2s,slpnor); //! Compute strain value
            }

            veps=(eps-d1epspr[i])/dcstec; //! Strain rate: (current strain - previous strain) / time step size
            if (veps*d1vepspr[i]<R0) //! Current strain rate x previous strain rate < 0 (i.e., change from loading to unloading or viceversa)
            { d1epsv[i]=d1epspr[i]; //! Storing strain from previous time step (i.e., strain at the time of strain rate inversion)
              d1sigmav[i]=d1sigmapr[i]; //! Storing stress from previous time step (i.e., stress at the time of strain rate inversion)
            }
        
            // Steel material model (See Fig. 15)
            eb=(-young/R6)*log(10*(d1epsmax[i]-epsy)); //! E_b in eq. 40
            koef=young/(young-eb); //! a in eq. 40
            if(d1sigmav[i]<=d1sigmavmin[i]) 
            { d1sigmavmin[i]=d1sigmav[i]; //! Storing minimum value of stress at the time of strain inversion reached so far
              d1epssvmin[i]=d1epsv[i]; //! Storing associated minimum value of strain at the time of strain inversion reached so far
            }
      
            if((d1sigmav[i]<R0 && ABS(d1sigmav[i]-FE4(fy,koef,eb,d1epspr[i],d1epssmax[i],d1sigmasmax[i]))<EPSILON ) ||
               (d1sigmav[i]<R0 && ABS(FE1(eps, d1sigmavmin[i],d1epssvmin[i],young)-d1sigmav[i])>EPSILON))
            { d1sigmak[i]=d1sigmav[i]+fy; //! Local parameter for FE5
              d1epsk[i]=d1epsv[i]+epsy; //! Local parameter for FE5
            }
        
            //! Check if bar is broken 
            if(eps>epsbr || mod<=SQRT(EPSILON)|| delta>=3*d1sdiam[i1refbar[i0r]] )   
            { dsigma=R0;
              i1rprop[i2relto[0][i]]=iprop-YIPROPMAX;
              i1rprop[i2relto[1][i]]=iprop-YIPROPMAX;
            }
 
            /* Calculation of stress from strain */

            k1e=(sigmash-sigmasu)/((epssh-epsu)*(epssh-epsu));
            k2e=2*epsu*(sigmasu-sigmash)/((epssh-epsu)*(epssh-epsu));
            k3e=(epssh*epssh*sigmasu-2*epssh*epsu*sigmasu+epsu*epsu*sigmash)/((epssh-epsu)*(epssh-epsu));
            k4e=-sigmasu/(epsbr*epsbr-2*epsbr*epsu+epsu*epsu);
            k5e=2*epsu*sigmasu/(epsbr*epsbr-2*epsbr*epsu+epsu*epsu);
            k6e=epsbr*sigmasu*(epsbr-2*epsu)/(epsbr*epsbr-2*epsbr*epsu+epsu*epsu);
        
            /* Unloading */
            if(veps<R0 && i1rprop[i2relto[1][i]]==iprop)
            { if((d1sigmapr[i]<R0) && (FE1(eps, d1sigmav[i],d1epsv[i],young)<FE4(fy,koef,eb,eps,d1epssmax[i],d1sigmasmax[i])))
              { dsigma = FE4(fy,koef,eb,eps,d1epssmax[i],d1sigmasmax[i]);
              }
              else 
              { dsigma = FE1(eps, d1sigmav[i],d1epsv[i],young);
              }
            }
            /* Loading */
            if (veps>=R0 && i1rprop[i2relto[1][i]]==iprop)
            { if(eps>d1epsk[i] &&
                 FE1(eps,d1sigmav[i],d1epsv[i],young)>FE5(fy,koef,eb,eps, d1epsk[i],d1sigmak[i])
                 && d1epsmax[i]>epsy
                 && d1sigmavmin[i]<0
                 && ((FE5(fy,koef,eb,eps, d1epsk[i],d1sigmak[i])<FE2(eps, fy, epsy, small) && (eps<=epssh ))
                     || (FE5(fy,koef,eb,eps, d1epsk[i],d1sigmak[i])<FE3(k1e,k2e,k3e,eps)  &&(eps<epsu))
                     || (FE5(fy,koef,eb,eps, d1epsk[i],d1sigmak[i])< FE6(k4e,k5e,k6e,eps) &&(eps>epsu)))) 
              { dsigma = FE5(fy,koef,eb,eps, d1epsk[i],d1sigmak[i]);
              }
              else if(eps>epsy && eps<=epssh && FE1(eps, d1sigmav[i],d1epsv[i],young)>FE2(eps, fy, epsy, small)) 
              { dsigma = FE2(eps, fy, epsy, small);
              }
              else if(eps>epssh && eps<=epsu && FE1(eps, d1sigmav[i],d1epsv[i],young)>FE3(k1e,k2e,k3e,eps)) 
              { dsigma = FE3(k1e,k2e,k3e,eps);
              }
              else if (eps>epsu && eps<=epsbr && FE1(eps, d1sigmav[i],d1epsv[i],young)>FE6(k4e,k5e,k6e,eps))  
              { dsigma = FE6(k4e,k5e,k6e,eps);
              }
              else 
              { dsigma = FE1(eps, d1sigmav[i],d1epsv[i],young);
              }
            }
     
            //dsigma = 0.0;
          }
          else /* The bar is infinitely elastic */
          { slpnor=slpnor-d1slpin[i]; //! Current normalized slip - normalized slip when joint @ 1D joint reference point is softening
            delta=delta-d1deltain[i]; //! Current Delta - delta when joint @ 1D joint reference point is softening
            slp=slpnor/pow((d1mpsfc[i1sbpr[i1refbar[i0r]]]/20.0),(2.0/3.0))*d1sdiam[i1refbar[i0r]]; //! Un-normalized slpnor
            fy=d1sfy[i1sbpr[i1refbar[i0r]]]; // Yield stress of steel bar

            /* Determine shear stress */
            sma=MYPI*pow(d1sdiam[i1refbar[i0r]],4.0)/64.0; // Area moment of inertia of bar: I_2 = (pi/4.0) * r^4
            //! Foundation stiffness for concrete: k = 150 * UCS_concrete^(0.85)* D_bar (eq. 31)
            stif=150.0*pow((d1mpsfc[i1sbpr[i1refbar[i0r]]]),0.85)/d1smmdiam[i1refbar[i0r]]*
                 (d1sfc[i1sbpr[i1refbar[i0r]]]/d1mpsfc[i1sbpr[i1refbar[i0r]]])*
                 (d1smmdiam[i1refbar[i0r]]/d1sdiam[i1refbar[i0r]]);
            //! Compute the length of the curvature-influenzing zone
            lc=(3.0*MYPI)/4.0*pow((4.0*young*sma/(stif*d1sdiam[i1refbar[i0r]])),0.25); //! Value derived from modelling the bar as beam resting on an elastic foundation
            //! Empirical non-dimensional parameter (Eq. 29)
            di=(1+150*slp/d1sdiam[i1refbar[i0r]])*(delta/d1sdiam[i1refbar[i0r]]);
            if(di>0.02)
            { lc=lc*(1.0+3.0*pow((di-0.02),0.8)); //! Eq. (28)
            }
            dtau=((young*sma*4.0)/(d1sdiam[i1refbar[i0r]]*d1sdiam[i1refbar[i0r]]*MYPI))*
                 ((34.9091*delta)/(lc*lc*lc)); //! Shear stress in the  1D joint element (eq. 34)
            mod=SQRT(MAXIM(EPSILON,1.0-3.0*pow((dtau/fy),2))); //! Reducing factor for the yield stress, fy (eq. 35, Von Mises criterion + isotropic hardenind)

            /* Determine critical strain and stress values of steel */
            epsy=(fy/young)*mod;  // Yield strain of bar  
            epssh=d1epssh[i1sbpr[i1refbar[i0r]]]; // Hardening strain of bar (from input) 
            epsu=d1epsu[i1sbpr[i1refbar[i0r]]];  // Ultimate strain of bar (from input)
            epsbr=d1epsbr[i1sbpr[i1refbar[i0r]]];  // Breakage strain of bar (from input)
            fy=d1sfy[i1sbpr[i1refbar[i0r]]]*mod; // Yield stress of bar (reduced by shear stress factor) 
            sigmash=(fy+(epssh-epsy)*small); // Hardening stress of bar (small hardening slope)
            sigmasu=d1sfu[i1sbpr[i1refbar[i0r]]]*mod;  // Utimate stress of bar (from input) reduced by shear stress factor
            sigmabr=d1sfbr[i1sbpr[i1refbar[i0r]]]*mod;  // Breakage stress of bar (from input) reduced by shear stress factor 

            /* Compute reduction factor to account for adjacent cracks, eq. 38 */ 
            /*cralfa=1.0-exp(-pow((0.065*d1crlcr[i1refbar[i0r]]/d1sdiam[i1refbar[i0r]]+0.5),3)); 
            if(cralfa>=0.087*d1crlcr[i1refbar[i0r]]/d1sdiam[i1refbar[i0r]])
            { cralfa=0.087*d1crlcr[i1refbar[i0r]]/d1sdiam[i1refbar[i0r]];
            }*/
            cralfa=1.0;
 
            /* Calculation of critical adimensional slip values based on the relation by Soltani and Maekawa (2008)*/
            sang=(epsy*(ka+3500.0*epsy))*cralfa; //! If epsilon_y < epsilon < epsilon_sh (eq. 19) (ka = 2.0)                                                                             
            syld=sang-sang*d1sigmacr[i]/(epsy*young); //! Yielding slip ?                                                             
            eps0=epssh+(0.06-syld/2.0)/(0.013*(sigmasu-fy)* 
                 (d1mpsfc[i1sbpr[i1refbar[i0r]]]/d1sfc[i1sbpr[i1refbar[i0r]]])); //! (eq. 22)
            slpnor0=(0.047*(sigmasu-fy)*(d1mpsfc[i1sbpr[i1refbar[i0r]]]/d1sfc[i1sbpr[i1refbar[i0r]]])*(eps0-epssh))*cralfa+syld; //! If epsilon_sh < epsilon < epsilon_0 (eq. 21) 
       
            /* Calculation of strain from slip */
       
            /* Compute maximum strain values (differential normalized slip is greater than any previous value (i.e., loading conditions)) 
               Piecewise bi-linear relationship */
            if(slpnor>d1snormax[i])  
            { d1snormax[i]=slpnor; //! Store differential normalized slip value 
              if(d1snormax[i]<=syld) //! Differential normalized slip value is <= yielding slip  
              { k1s=sang/epsy; 
                k2s=-k1s*d1sigmacr[i]/young; //! d1sigmacr is the stress in the bar at the time step previous to softening (i.e., when crack opens)
                d1epsmax[i]=FE1s(k1s,k2s,d1snormax[i]); //! Compute and store strain (highest value reached so far)
              }
              else if (d1snormax[i]<=syld+small) //! Differential normalized slip value is <= yielding slip + small
              { k1s=(syld+small)/(epssh-epsy);
                k2s=-k1s*epsy+syld;
                d1epsmax[i]=FE1s(k1s,k2s,d1snormax[i]); //! Compute and store strain (highest value reached so far)
              }
              else if (d1snormax[i]<=slpnor0) //! Differential normalized slip value is <= slip_0 value (See Fig. 9)
              { k1s=(0.047*(sigmasu-fy)*(d1mpsfc[i1sbpr[i1refbar[i0r]]]/
                     d1sfc[i1sbpr[i1refbar[i0r]]]))*cralfa;
                k2s=(-0.047*(sigmasu-fy)*(d1mpsfc[i1sbpr[i1refbar[i0r]]]/
                     d1sfc[i1sbpr[i1refbar[i0r]]])*epssh)*cralfa+syld+small;
                d1epsmax[i]=FE1s(k1s,k2s,d1snormax[i]); //! Compute and store strain (highest value reached so far)
              }
              else //!  Differential normalized slip value is > slip_0 value (See Fig. 9)
              { k1s=(0.007*(sigmasu-fy)*(d1mpsfc[i1sbpr[i1refbar[i0r]]]/
                     d1sfc[i1sbpr[i1refbar[i0r]]]))*cralfa;
                k2s=(-0.007*(sigmasu-fy)*(d1mpsfc[i1sbpr[i1refbar[i0r]]]/
                     d1sfc[i1sbpr[i1refbar[i0r]]])*eps0)*cralfa+slpnor0;
                d1epsmax[i]=FE1s(k1s,k2s,d1snormax[i]); //! Compute and store strain (highest value reached so far)
              }
            } 
        
            k1s=sang/epsy;                                                                              
            k2s=-k1s*d1sigmacr[i]/young;
            eps=FE1s(k1s,k2s,slpnor); //! Compute strain value

            veps=(eps-d1epspr[i])/dcstec; //! Strain rate: (current strain - previous strain) / time step size
            if (veps*d1vepspr[i]<R0) //! Current strain rate x previous strain rate < 0 (i.e., change from loading to unloading or viceversa)
            { d1epsv[i]=d1epspr[i]; //! Storing strain from previous time step (i.e., strain at the time of strain rate inversion)
              d1sigmav[i]=d1sigmapr[i]; //! Storing stress from previous time step (i.e., stress at the time of strain rate inversion)
            }
        
            // Steel material model (See Fig. 15)
            eb=(-young/R6)*log(10*(d1epsmax[i]-epsy)); //! E_b in eq. 40
            koef=young/(young-eb); //! a in eq. 40
            if(d1sigmav[i]<=d1sigmavmin[i]) 
            { d1sigmavmin[i]=d1sigmav[i]; //! Storing minimum value of stress at the time of strain inversion reached so far
              d1epssvmin[i]=d1epsv[i]; //! Storing associated minimum value of strain at the time of strain inversion reached so far
            }
      
            if((d1sigmav[i]<R0 && ABS(d1sigmav[i]-FE4(fy,koef,eb,d1epspr[i],d1epssmax[i],d1sigmasmax[i]))<EPSILON ) ||
               (d1sigmav[i]<R0 && ABS(FE1(eps, d1sigmavmin[i],d1epssvmin[i],young)-d1sigmav[i])>EPSILON))
            { d1sigmak[i]=d1sigmav[i]+fy; //! Local parameter for FE5
              d1epsk[i]=d1epsv[i]+epsy; //! Local parameter for FE5
            }
        
            /* Calculation of stress from strain */
        
            dsigma = FE1(eps, d1sigmav[i],d1epsv[i],young);

          }  
        }
      }
      /* Update total forces */
            
      /* Compute reduction factor for shear stress */
      
      /* if (ABS(delta)<(2*d1sdiam[i1refbar[i0r]]) && eps<d1epsu[i1sbpr[i1refbar[i0r]]])
      { reduc=1.0; }
      else if(ABS(delta)>(2*d1sdiam[i1refbar[i0r]]) && eps>d1epsu[i1sbpr[i1refbar[i0r]]])
      { reduc=(1.0-(ABS(delta)-2*d1sdiam[i1refbar[i0r]])/(3*d1sdiam[i1refbar[i0r]]-2*d1sdiam[i1refbar[i0r]]))*
              (1.0-(eps-d1epsu[i1sbpr[i1refbar[i0r]]])/(epsbr-d1epsu[i1sbpr[i1refbar[i0r]]])); }
      else if(ABS(delta)>(2*d1sdiam[i1refbar[i0r]]))
      { reduc=(1.0-(ABS(delta)-2*d1sdiam[i1refbar[i0r]])/(3*d1sdiam[i1refbar[i0r]]-2*d1sdiam[i1refbar[i0r]])); }
      else 
      { reduc=(1.0-(eps-d1epsu[i1sbpr[i1refbar[i0r]]])/(d1epsbr[i1sbpr[i1refbar[i0r]]]-d1epsu[i1sbpr[i1refbar[i0r]]])); }
      if (ABS(delta)>(3*d1sdiam[i1refbar[i0r]]) && eps>d1epsbr[i1sbpr[i1refbar[i0r]]])
      { reduc=0.0; }
     
      dtau=dtau*reduc; */
      
      //fprintf(out2,"%.6f \t %.6f \t %.6f \n",dctime,dsigma,dtau); 
      
      /* Set normal and tangential stress to zero if the 1D joint is 
         of type 5, 6 or 7 and has not been actiavated yet */
      if(((i1sbty[i1refbar[i0r]]==5)||(i1sbty[i1refbar[i0r]]==6)||(i1sbty[i1refbar[i0r]]==7))&&(i1sbac[i1refbar[i0r]]==0))
      { dsigma=0.0;
        dtau=0.0;
        //printf("Non-active\n");
      }
     
      //printf("dsigma: %.6f\n",dsigma);
      /* Store normal and tangential stress for output */
      d1rjsig[i]=dsigma;
      d1rjtau[i]=dtau;
      /* Store normal and tangential strains for output (only if i1sbty[i1refbar[i0r]] == 5, 6 or 7) */
      if((i1sbty[i1refbar[i0r]]==5)||(i1sbty[i1refbar[i0r]]==6)||(i1sbty[i1refbar[i0r]]==7))
      { d1rjslpnor[i]=normal_displacement;
        d1rjdelnor[i]=tangential_displacement;
      }
      else
      { d1rjslpnor[i]=0.0;
        d1rjdelnor[i]=0.0;
      }
      
      /* Check if at least one of the two adjacent triangles has been excavated */
      if((i1pexc[i1elpr[i1rmyel[i0r]]]==1)||(i1pexc[i1elpr[i1rmyel[i1r]]]==1))
      { i1rprop[i2relto[0][i]]=iprop-YIPROPMAX; //! Set first reference point of 1D joint to broken
        i1rprop[i2relto[1][i]]=iprop-YIPROPMAX; //! Set second reference point of 1D joint to broken
      }
      else /* Update forces only if the triangles have not been excavated */
      { //fprintf(out2,"dsigma: %.6f \n",dsigma);
        dnforce=dsigma*d1spea[i1refbar[i0r]]; //! Normal force
        dtforce=dtau*d1spea[i1refbar[i0r]]; //! Shear force
        
        d2rjfrv[0][i]=dnforce*rx+dtforce*rx; //! x-component of force
        d2rjfrv[1][i]=dnforce*ry+dtforce*ry; //! y-component of force

        d1nfcx[i0]=d1nfcx[i0] + dnforce*rx*di0ri1/di0i1; 
        d1nfcx[i1]=d1nfcx[i1] + dnforce*rx*(di0i1-di0ri1)/di0i1;
        d1nfcx[i2]=d1nfcx[i2] - dnforce*rx*(di2i3-di1ri2)/di2i3;
        d1nfcx[i3]=d1nfcx[i3] - dnforce*rx*(di1ri2/di2i3);
        d1nfcy[i0]=d1nfcy[i0] + dnforce*ry*di0ri1/di0i1;
        d1nfcy[i1]=d1nfcy[i1] + dnforce*ry*(di0i1-di0ri1)/di0i1;
        d1nfcy[i2]=d1nfcy[i2] - dnforce*ry*(di2i3-di1ri2)/di2i3;
        d1nfcy[i3]=d1nfcy[i3] - dnforce*ry*(di1ri2/di2i3);
 
        d1nfcx[i0]=d1nfcx[i0] - dtforce*ry*di0ri1/di0i1; 
        d1nfcx[i1]=d1nfcx[i1] - dtforce*ry*(di0i1-di0ri1)/di0i1;
        d1nfcx[i2]=d1nfcx[i2] + dtforce*ry*(di2i3-di1ri2)/di2i3;
        d1nfcx[i3]=d1nfcx[i3] + dtforce*ry*(di1ri2/di2i3);
        d1nfcy[i0]=d1nfcy[i0] + dtforce*rx*di0ri1/di0i1;
        d1nfcy[i1]=d1nfcy[i1] + dtforce*rx*(di0i1-di0ri1)/di0i1;
        d1nfcy[i2]=d1nfcy[i2] - dtforce*rx*(di2i3-di1ri2)/di2i3;
        d1nfcy[i3]=d1nfcy[i3] - dtforce*rx*(di1ri2/di2i3);
	
	

        /* Update control points for stress-strain */
        if(indic[i]==1) //! If crack is open
        { if(dsigma*d1sigmapr[i]<R0) //! If goinf from tension to compression or viceversa  
          { d1epss[i]=eps;  //! Store strain when stress changes sign
            d1sigmas[i]=dsigma; //! Store stress when stress changes sign
          }
          if(d1epss[i]>=d1epssmax[i]) 
          { d1sigmasmax[i]=d1sigmas[i]; //! Store maximum stress when stress changes sign
            d1epssmax[i]=d1epss[i]; //! Store maximum strain when stress changes sign
          }
          d1epspr[i]=eps; //! Store strain as previous strain
          d1vepspr[i]=veps;  //! Store strain rate as previous strain rate
        }
        d1sigmapr[i]=dsigma; //! Store stress as previous stress
      }
    }
  } 
  #undef FE1
  #undef FE2
  #undef FE3
  #undef FE4
  #undef FE5
  #undef FE6
  #undef FE1s
} 


/*********************PUBLIC********************************************************/
void Yrb(ydc,yde,ydn,ydpe,ydpj,ydps,ydr,ydsb,ydo   /***  RBAR elements  ***/
        )
  YDC ydc; YDE yde; YDN ydn; YDPE ydpe; YDPJ ydpj; YDPS ydps; YDR ydr; YDSB ydsb; YDO ydo; 
{ INT i=0,j,k,iprop,indic=0,temp;
  INT i0,i1,i2,ir,ij1,ij2;
  DBL rx0,rx1,rx2,rxr,ry0,ry1,ry2,ryr;
  DBL small=1.0e-9;
  INT isbar;

  /* initializing d1rbsig */
  if(ydr->d1rbsig==DBL1NULL)
  { ydr->d1rbsig=TalDBL1(ydr->mrepo);
    for(i=0;i<ydr->mrepo;i++)
    { ydr->d1rbsig[i]=R0;
    }
  }
  /* initializing d1rbfrc */
  if(ydr->d1rbfrc==DBL1NULL)
  { ydr->d1rbfrc=TalDBL1(ydr->mrepo);
    for(i=0;i<ydr->mrepo;i++)
    { ydr->d1rbfrc[i]=R0;
    }
  }
  /* initializing d1rbstr */
  if(ydr->d1rbstr==DBL1NULL)
  { ydr->d1rbstr=TalDBL1(ydr->mrepo);
    for(i=0;i<ydr->mrepo;i++)
    { ydr->d1rbstr[i]=R0;
    }
  }
  /* initializing d1rjsig */
  if(ydr->d1rjsig==DBL1NULL)
  { ydr->d1rjsig=TalDBL1(ydr->mrepo);
    for(i=0;i<ydr->mrepo;i++)
    { ydr->d1rjsig[i]=R0;
    }
  }
  /* initializing d1rjtau */
  if(ydr->d1rjtau==DBL1NULL)
  { ydr->d1rjtau=TalDBL1(ydr->mrepo);
    for(i=0;i<ydr->mrepo;i++)
    { ydr->d1rjtau[i]=R0;
    }
  }
  /* initializing d1rjslpnor */
  if(ydr->d1rjslpnor==DBL1NULL)
  { ydr->d1rjslpnor=TalDBL1(ydr->mrepo);
    for(i=0;i<ydr->mrepo;i++)
    { ydr->d1rjslpnor[i]=R0;
    }
  }
  /* initializing d1rjdelnor */
  if(ydr->d1rjdelnor==DBL1NULL)
  { ydr->d1rjdelnor=TalDBL1(ydr->mrepo);
    for(i=0;i<ydr->mrepo;i++)
    { ydr->d1rjdelnor[i]=R0;
    }
  }
  /* initializing i1rjstnor */
  if(ydr->i1rjstnor==INT1NULL)
  { ydr->i1rjstnor=TalINT1(ydr->mrepo);
    for(i=0;i<ydr->mrepo;i++)
    { ydr->i1rjstnor[i]=1;
    }
  }
  /* initializing d2rjfrv only if rebars are used */
  if(ydr->d2rjfrv==DBL2NULL)
  { ydr->d2rjfrv=TalDBL2(2,ydr->mrepo);
    for(i=0;i<2;i++)
    { for(j=0;j<ydr->mrepo;j++)
      ydr->d2rjfrv[i][j]=R0;
    }
  }
  /* initializing i2rbedn only if rebars are used */
  if(ydr->i2rbedn==INT2NULL)
  { ydr->i2rbedn=TalINT2(2,ydr->mrepo);
    for(i=0;i<2;i++)
    { for(j=0;j<ydr->mrepo;j++)
      ydr->i2rbedn[i][j]=-1;
    }
  }
  
  
  //printf("ydr->nrepo: %d\n",ydr->nrepo);
  //printf("ydr->nrepo: |%d|\n",ydr->nrepo);
  
  if(ydsb->nsbar>0)
  { if(ydsb->isfirst==0)
    { ydsb->isfirst=1;
      Yrb2CREATE(          /* Create rebar elements  */
       yde->nelem,ydsb->nsbar,   ydr->mrepo,   &(ydr->nrepo),
       ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],ydn->d2nvc[1], 
       ydsb->d2sic[0],ydsb->d2sic[1],ydsb->d2sic[2],ydsb->d2sic[3],ydr->d2rcig[0],ydr->d2rcig[1], 
       ydr->d2rccg[0],ydr->d2rccg[1],ydr->d2rvcg[0],ydr->d2rvcg[1], 
       ydr->d2riLc[0],ydr->d2riLc[1],ydr->d2riLc[2],
       yde->i1elpr,yde->i2elto,ydpe->i1ptyp,ydr->i1rrpn,ydsb->i1srpf, 
       ydr->i1rmyel,ydsb->d1spea,ydr->i2relto,ydr->i1rprop,ydr->i1refbar,
       ydsb->i1sbpr,    
       ydr->d2rsctr[0], ydr->d2rsctr[1],ydr->i2rbedn
      );
      i=i;
    } 
    else
    { /* If rebars are activated, set the initial coordinates of reference points equal to current coordinates */
      /* This is needed to have stress-free rebars at the time of their installation */
      for(isbar=0;isbar<ydsb->nsbar;isbar++) //! Loop over rebars
      { if(ydsb->i1sbac[isbar]==1) //! If rebar is activated
        { for (i=0; i<ydr->nrepo; i++) //! Loop over reference points
          { if(ydr->i1refbar[i]==isbar) //! Check if reference point belongs to that rebar
            { ydr->d2rcig[0][i]=ydr->d2rccg[0][i];
              ydr->d2rcig[1][i]=ydr->d2rccg[1][i];
            } 
          }
          ydsb->i1sbac[isbar]=2; //! Flag rebar with "2"
        }
      }  
      
      Yrb2CARENT(          /* Compute current coordinates of rebar elements  */
      yde->nelem,ydsb->nsbar,&(ydr->nrepo),
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nvc[0],ydn->d2nvc[1],
      ydr->d2rccg[0],ydr->d2rccg[1],ydr->d2rvcg[0],ydr->d2rvcg[1], 
      ydr->d2riLc[0],ydr->d2riLc[1],ydr->d2riLc[2],
      yde->i2elto,ydr->i1rrpn,ydsb->i1srpf, 
      ydr->i1rmyel, 
      ydr->d2rsctr[0], ydr->d2rsctr[1]
      );
    }
    
    
    //printf("ydr->nrepo: %d\n",ydr->nrepo);
    //printf("ydr->d2rcig[0][720]: %0.6f \t ydr->d2rcig[1][720]: %0.6f\n",ydr->d2rcig[0][720],ydr->d2rcig[1][720]);
    //printf("ydr->d2rcig[0][721]: %0.6f \t ydr->d2rcig[1][721]: %0.6f\n",ydr->d2rcig[0][721],ydr->d2rcig[1][721]);
    //printf("ydr->i1refbar[720]: %d \t ydr->i1refbar[721]: %d\n",ydr->i1refbar[720],ydr->i1refbar[721]);
    
    //printf("ydr->d2rcig[0][48]: %0.6f \t ydr->d2rcig[1][48]: %0.6f\n",ydr->d2rcig[0][48],ydr->d2rcig[1][48]);
    //printf("ydr->d2rcig[0][49]: %0.6f \t ydr->d2rcig[1][49]: %0.6f\n",ydr->d2rcig[0][49],ydr->d2rcig[1][49]);
    //printf("ydr->i1refbar[48]: %d \t ydr->i1refbar[49]: %d\n",ydr->i1refbar[48],ydr->i1refbar[49]);
    
    //printf("ydr->d2rcig[0][744]: %0.6f \t ydr->d2rcig[1][744]: %0.6f\n",ydr->d2rcig[0][744],ydr->d2rcig[1][744]);
    //printf("ydr->d2rcig[0][745]: %0.6f \t ydr->d2rcig[1][745]: %0.6f\n",ydr->d2rcig[0][745],ydr->d2rcig[1][745]);
    //printf("ydr->i1refbar[744]: %d \t ydr->i1refbar[745]: %d\n",ydr->i1refbar[744],ydr->i1refbar[745]);
    
    for (iprop=0;iprop<ydps->nprop;iprop++)
    { Yrb2FORCE(          /* Compute nodal forces due to rebar elements */
      yde->nelem,
      iprop,
      ydsb->nsbar, &(ydr->nrepo),  
      ydps->d1young[iprop],   
      ydr->d2rcig[0], ydr->d2rcig[1], ydr->d2rccg[0], ydr->d2rccg[1], 
      ydn->d2nfc[0],  ydn->d2nfc[1],  ydr->d2riLc[0], ydr->d2riLc[1], ydr->d2riLc[2],
      ydn->d2nci[0],  ydn->d2nci[1],
      ydsb->i1srpf,   ydr->i1rprop,      ydr->i1rrpn,    ydr->i1rmyel, ydr->i1refbar,   
      yde->i2elto,    ydsb->d1spea,
      ydo->nohys,     ydo->d1ohyx,    ydo->d1ohyy, ydo->d1ohys,ydo->d1ohyt,ydo->dohyp,
      ydo->i1ohyt,
      ydc->dctime, ydc->ncstep,
      ydr->d1rbsig, ydr->d1rbfrc, ydr->d1rbstr,
      ydsb->i1sbac
      ); 
    }
   
    
    /* Allocate memory and initialize array of 1D joint elements */
    if(ydr->i1myjoint==INT1NULL)
    { ydr->i1myjoint=TalINT1(ydr->mrepo);
      
      for(j=0; j<ydr->mrepo;j++)
      { ydr->i1myjoint[j]=-1; }
            
      for (j=0;ydr->i2relto[0][j]>=0; j++)
      { //printf("j: %d \t ydr->i2relto[0][j]: %d\n",j,ydr->i2relto[0][j]);
        i0=yde->i2elto[0][ydr->i1rmyel[ydr->i2relto[0][j]]];
        i1=yde->i2elto[1][ydr->i1rmyel[ydr->i2relto[0][j]]];
        i2=yde->i2elto[2][ydr->i1rmyel[ydr->i2relto[0][j]]];
        ir=ydr->i2relto[0][j];
        
        /* Y-RC used initial coordinates to determine 1D joint elements */
        /*rx0=ydn->d2nci[0][i0];
        rx1=ydn->d2nci[0][i1];
        rx2=ydn->d2nci[0][i2];
        rxr=ydr->d2rcig[0][ir];
        ry0=ydn->d2nci[1][i0];
        ry1=ydn->d2nci[1][i1];
        ry2=ydn->d2nci[1][i2];
        ryr=ydr->d2rcig[1][ir];*/
        
        /* Using current coordinates to determine the 1D joint elements */
        /* Doing so because the initial coordinates of the reference points are set equal to the current coordinates, see Yrb2CREATE */
        rx0=ydn->d2ncc[0][i0];
        rx1=ydn->d2ncc[0][i1];
        rx2=ydn->d2ncc[0][i2];
        rxr=ydr->d2rccg[0][ir];
        ry0=ydn->d2ncc[1][i0];
        ry1=ydn->d2ncc[1][i1];
        ry2=ydn->d2ncc[1][i2];
        ryr=ydr->d2rccg[1][ir];

        //printf("i0: %d \t i1: %d \t i2: %d \t ir: %d\n",i0,i1,i2,ir);
        //printf("rx0: %f \t rx1: %f \t rx2: %f \t rxr: %f\n",rx0,rx1,rx2,rxr);
        //printf("ry0: %f \t ry1: %f \t ry2: %f \t ryr: %f\n",ry0,ry1,ry2,ryr);

        if(ABS((rx0-rxr)*(ry1-ryr)-(rx1-rxr)*(ry0-ryr))<small)
        { ij1=i0;
          ij2=i1;
        }
        else if (ABS((rx1-rxr)*(ry2-ryr)-(rx2-rxr)*(ry1-ryr))<small)
        { ij1=i1;
          ij2=i2;
        }
        else if (ABS((rx0-rxr)*(ry2-ryr)-(rx2-rxr)*(ry0-ryr))<small)
        { ij1=i0;
          ij2=i2;
        }
        
        for(k=0;(indic!=1 && k<yde->nelem); k++)
        { if((yde->i2elto[0][k]==ij1 || yde->i2elto[1][k]==ij1 || yde->i2elto[2][k]==ij1  || yde->i2elto[3][k]==ij1) &&
             (yde->i2elto[0][k]==ij2 || yde->i2elto[1][k]==ij2 || yde->i2elto[2][k]==ij2  || yde->i2elto[3][k]==ij2) &&
             (yde->i2elto[3][k]!=(-1)))
          { ydr->i1myjoint[j]=k;
            indic=1;
          }
        }
        
        indic=0;
        if(ij1==yde->i2elto[0][k-1] || ij1==yde->i2elto[1][k-1]) {}
        else
        { temp=ydr->i2relto[0][j];
          ydr->i2relto[0][j]=ydr->i2relto[1][j];
          ydr->i2relto[1][j]=temp;
        }
        //printf("j: %d \t i1myjoint[j]: %d \t ydr->i2relto[0][j]: %d \t ydr->i2relto[1][j]: %d\n",j,ydr->i1myjoint[j],ydr->i2relto[0][j],ydr->i2relto[1][j]);
      }
      
      /* Determine the total number of 1D joint elements */
      for(j=0;ydr->i2relto[0][j]>=R0; j++)
      { ydr->nbrjointrb=ydr->nbrjointrb+1;
        //printf("ydr->nbrjointrb: %d\n",ydr->nbrjointrb);
      }
      //printf("ydr->nbrjointrb: %d\n\n",ydr->nbrjointrb);
    }
    
    /* Compute forces due to 1D joint elements */
    for (iprop=0;iprop<ydps->nprop;iprop++)
    { Yrb1RJOINT(  
      ydr->nbrjointrb, yde->nelem,
      iprop,
      ydps->d1young[iprop],ydc->dcstec,   
      ydr->d2rccg[0], ydr->d2rccg[1],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nfc[0],ydn->d2nfc[1],ydn->d2nvc[0],
      ydn->d2nvc[1],ydpj->d1pjft,ydpj->d1pjco,ydpj->d1pjgs, ydpj->d1pjpe,ydsb->d1spea, 
      ydsb->d1sdiam,ydsb->d1smmdiam,ydsb->d1crlcr, ydsb->i1sbpr,
      ydps->d1sfc, ydps->d1mpsfc, 
      ydps->d1sfy, ydps->d1epssh, ydps->d1sfu, ydps->d1epsu, ydps->d1sfbr, ydps->d1epsbr,  
      yde->i1elpr,ydr->i1rprop,yde->i2elto,ydr->i2relto,ydr->mrepo, ydr->i1refbar, 
      ydr->i1rmyel, 
      ydr->i1myjoint,
      ydr->d2rsctr[0],ydr->d2rsctr[1],ydr->d2rvcg[0],ydr->d2rvcg[1],
      ydpe->d1pela, ydpe->d1pemu, 
      ydpe->nprop,
      ydc->dctime, ydc->ncstep,
      ydpj->i1ptyp,
      ydr->d1rjsig, ydr->d1rjtau,
      ydpe->i1ptyp,ydpe->d1peem, ydpe->d1peex, ydpe->d1peey, ydpe->i1usan,
      ydsb->i1sbty, yde->d1elfs,
      ydpe->i1pexc, ydsb->i1sbac,
      ydps->d1stkn,ydps->d1stkt,
      ydr->d1rjslpnor, ydr->d1rjdelnor,
      ydps->d1styns,ydps->d1strns,
      ydr->i1rjstnor,
      yde->d2elstr,
      ydps->d1stcoh,ydps->d1stfri,
      ydr->d2rjfrv,ydr->i2rbedn
      );
    }
  }
}
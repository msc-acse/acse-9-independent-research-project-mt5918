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
/* File   Yfrd.c */
#include "Yproto.h"
/**************FRACTURE ELEMENTS***********/
static void Yfrd2TRIANGLE(  /* fracture triangle  */
             melem,mnopo ,nelest,nnopst,
            iprop ,
            n0elem,n0nopo,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,d1elfs,d1sdel,i1elpr,i1nobf,
            i1nopr,i2elto
            )
  INT    melem; INT   mnopo; INT  nelest; INT  nnopst;
  INT    iprop;
  INT  *n0elem; INT  *n0nopo;
  DBL  *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL *d1nciy; DBL *d1nvcx;
  DBL  *d1nvcy; DBL  *d1elfs; DBL *d1sdel; INT *i1elpr; INT *i1nobf;
  INT  *i1nopr; INT **i2elto;
{ DBL vx[3], vy[3], x[3], y[3], xi[3], yi[3];
  DBL xc,yc,vxc,vyc,ratio,volc,voli;
  INT nelem, nnopo, inopo;
  INT i,ielem;
  INT i1t[12]={0,2,4,8, 1,6,5,9,  11,10,3,7 };
  DBL dpefs;

  nelem=(*n0elem);
  nnopo=(*n0nopo);
  for(ielem=0;ielem<nelest;ielem++)
  { if(i1elpr[ielem]==iprop)
    { dpefs=d1elfs[ielem];
      if(d1sdel[ielem]>(RP15))
      { d1sdel[ielem]=-R1;
        for(i=0;i<3;i++)
        { inopo=i2elto[i][ielem];
          x[i]=d1nccx[inopo];
          y[i]=d1nccy[inopo];
          xi[i]=d1ncix[inopo];
          yi[i]=d1nciy[inopo];
          vx[i]=d1nvcx[inopo];
          vy[i]=d1nvcy[inopo];
          i1nobf[inopo]=1;
        }
        xc=(x[0]+x[1]+x[2])/R3;
        yc=(y[0]+y[1]+y[2])/R3;
        vxc=(vx[0]+vx[1]+vx[2])/R3;
        vyc=(vy[0]+vy[1]+vy[2])/R3;
        for(i=0;i<3;i++)
        { x[i]=x[i]-xc;
          y[i]=y[i]-yc;
          vx[i]=vx[i]-vxc;
          vy[i]=vy[i]-vyc;
        }
        volc=(x[1]-x[0])*(y[2]-y[0])-(y[1]-y[0])*(x[2]-x[0]);
        voli=(xi[1]-xi[0])*(yi[2]-yi[0])-(yi[1]-yi[0])*(xi[2]-xi[0]);
        ratio=dpefs*SQRT(MINIM(R1,(voli/volc)));
        for(i=0;i<3;i++)
        { if(nnopo>=mnopo)
          { CHRw(stderr,"Ymd: MNOPO too small"); CHRwcr(stderr);
            exit(1);
          }
          inopo=i2elto[i][ielem];
          d1nccx[nnopo]=xc+x[i]*ratio;
          d1nccy[nnopo]=yc+y[i]*ratio;
          d1ncix[nnopo]=d1nccx[nnopo];
          d1nciy[nnopo]=d1nccy[nnopo];
          d1nvcx[nnopo]=vxc;
          d1nvcy[nnopo]=vyc;
          i1nopr[nnopo]=i1nopr[inopo];
          i1nobf[nnopo]=1;
          i2elto[i][ielem]=nnopo;
          nnopo=nnopo+1;
  } } } }
  (*n0nopo)=nnopo;
  (*n0elem)=nelem;
}
/*********************PUBLIC********************************************************/
void Yfrd(   ydc, yde, ydi, ydn, ydpe, ydpn, ydpj, ydpm    /***  mesh elements  ***/
         )
  YDC ydc; YDE yde; YDI ydi;  YDN ydn; YDPE ydpe; YDPN ydpn; YDPJ ydpj; YDPM ydpm;
{ INT nelest, nnopst;
  INT iprop;

  nelest=yde->nelem;
  nnopst=ydn->nnopo;
  for(iprop=0;iprop<ydpe->nprop;iprop++)
  { if((ydpe->i1ptyp[iprop]==YTE2TRISOF)&&
       (ydpe->i1pefr[iprop]>0))
    { Yfrd2TRIANGLE(  /* fracture triangle  */
      yde->melem   ,ydn->mnopo   ,nelest       ,nnopst     ,
      iprop ,
      &(yde->nelem),&(ydn->nnopo),
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
      ydn->d2nvc[1],
      yde->d1elfs,
      yde->d2elst[ydpe->i1psde[iprop]],
      yde->i1elpr,ydn->i1nobf,
      ydn->i1nopr,yde->i2elto
      );
  } }
  if((ydn->nnopo)!=nnopst)ydi->diedi=R5*ydi->diezon;
}

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
/* File   Ysd.c */
#include "Yproto.h"
/**********************SOLVERS*****************************************/
static void Ysd2MEC(  nnopo,  /* mechanical solver for 2D nodes with x,y d.o.f.  */
                      dcstec,dcgray, dctime,npnfact,ibc,
                      d3pnfac,
                      d1nccx,d1nccy,d1nfcx,d1nfcy,d1nmct,
                      d1nvcx,d1nvcy,d1pnax,d1pnay,d1pnxx,
                      d1pnxy,d1pnyx,d1pnyy,i1nopr,
											i1pnfx,i1pnfy
                    )
  INT   nnopo;
  DBL  dcstec;  DBL  dcgray;  DBL  dctime;  INT   npnfact;  INT ibc;
  DBL ***d3pnfac;
  DBL *d1nccx;  DBL *d1nccy;  DBL *d1nfcx;  DBL *d1nfcy;  DBL *d1nmct;
  DBL *d1nvcx;  DBL *d1nvcy;  DBL *d1pnax;  DBL *d1pnay;  DBL *d1pnxx;
  DBL *d1pnxy;  DBL *d1pnyx;  DBL *d1pnyy;  INT *i1nopr;
	INT *i1pnfx;	INT *i1pnfy;
{ INT inopo,i;
  DBL T[2][2];
  DBL Tinv[2][2];
  DBL volc,vX,vY,aX,aY,fX,fY,fact;
  DBL *d1time;
  DBL *d1fact;

  T[0][0]=d1pnxx[ibc];   /* local curvlinear base in global coordinates */
  T[1][0]=d1pnxy[ibc];
  T[0][1]=d1pnyx[ibc];
  T[1][1]=d1pnyy[ibc];

  aX=d1pnax[ibc];
  aY=d1pnay[ibc];
  fact=R1;
  if(d3pnfac != DBL3NULL)
  { if(d3pnfac[0][0][0] != -R1)
    { d1time = d3pnfac[0][ibc];
      d1fact = d3pnfac[1][ibc];
      for(i=1; i<npnfact; i++)
      { if((dctime>=d1time[i-1])&&(dctime<=d1time[i]))
        { fact=d1fact[i-1]-((d1fact[i-1]-d1fact[i])*
               ((dctime-d1time[i-1])/(d1time[i]-d1time[i-1])));
  } } } }
  aX=aX*fact;
  aY=aY*fact;

  YMATINV2(T,Tinv,volc);
  for(inopo=0;inopo<nnopo;inopo++)
  { if((i1nopr[inopo]==ibc)&&(d1nmct[inopo]>EPSILON))
    {  /* velocity in local curvlinear coordinates */
       vX=Tinv[0][0]*d1nvcx[inopo]+Tinv[0][1]*d1nvcy[inopo];
       vY=Tinv[1][0]*d1nvcx[inopo]+Tinv[1][1]*d1nvcy[inopo];
       fX=Tinv[0][0]*d1nfcx[inopo]+
          Tinv[0][1]*(d1nfcy[inopo]+dcgray*d1nmct[inopo]);
       fY=Tinv[1][0]*d1nfcx[inopo]+
          Tinv[1][1]*(d1nfcy[inopo]+dcgray*d1nmct[inopo]);
       /* solve equations in local coordinates */
       if(i1pnfx[ibc]==1)      /* supplied force        */
       { vX=vX+((aX+fX)/d1nmct[inopo])*dcstec;
       }
       else if(i1pnfx[ibc]==2) /* supplied acceleration */
       { vX=vX+(aX+fX/d1nmct[inopo])*dcstec;
       }
       else if(i1pnfx[ibc]==3) /* supplied velocity     */
       { vX=aX;
       }
       else if(i1pnfx[ibc]==4) /* supplied force with absorbing boundary */
       { vX=vX+((aX+fX)/d1nmct[inopo])*dcstec;
       }
       if(i1pnfy[ibc]==1)
       { vY=vY+((aY+fY)/d1nmct[inopo])*dcstec;
       }
       else if(i1pnfy[ibc]==2)
       { vY=vY+(aY+fY/d1nmct[inopo])*dcstec;
       }
       else if(i1pnfy[ibc]==3)
       { vY=aY;
       }
       else if(i1pnfy[ibc]==4) /* supplied force with absorbing boundary */
       { vY=vY+((aY+fY)/d1nmct[inopo])*dcstec;
       }
       /* velocity back to global coordinates */
       d1nvcx[inopo]=T[0][0]*vX+T[0][1]*vY;
       d1nvcy[inopo]=T[1][0]*vX+T[1][1]*vY;
       /* new coordinates */
       d1nccx[inopo]=d1nccx[inopo]+d1nvcx[inopo]*dcstec;
       d1nccy[inopo]=d1nccy[inopo]+d1nvcy[inopo]*dcstec;
    }
  }
}
static void Ysd2TRIRIG( /* mechanical solver rigid triangles  */
             nelem,
            dcstec, dpero, iprop,
            d1nccx,d1nccy,d1nfcx,d1nfcy,d1nmct,
            d1nvcx,d1nvcy,i1elpr,i2elto
            )
  INT   nelem;
  DBL  dcstec; DBL   dpero; INT   iprop;
  DBL *d1nccx; DBL *d1nccy; DBL *d1nfcx; DBL  *d1nfcy; DBL *d1nmct;
  DBL *d1nvcx; DBL *d1nvcy; INT *i1elpr; INT **i2elto;
{ INT i,ielem;
  DBL mas,iner,xc,yc,vxc,vyc,fxc,fyc,omega,c,s,momf,momv;
  DBL fx[3];
  DBL fy[3];
  DBL x[3];
  DBL y[3];
  DBL xnew;
  DBL ynew;
  DBL vx[3];
  DBL vy[3];
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { for(i=0;i<3;i++)
      { x[i]=d1nccx[i2elto[i][ielem]];
        y[i]=d1nccy[i2elto[i][ielem]];
        vx[i]=d1nvcx[i2elto[i][ielem]];
        vy[i]=d1nvcy[i2elto[i][ielem]];
        fx[i]=d1nfcx[i2elto[i][ielem]];
        fy[i]=d1nfcy[i2elto[i][ielem]];
      }
      xc=(x[0]+x[1]+x[2])/R3;
      yc=(y[0]+y[1]+y[2])/R3;
      vxc=(vx[0]+vx[1]+vx[2])/R3;
      vyc=(vy[0]+vy[1]+vy[2])/R3;
      fxc=(fx[0]+fx[1]+fx[2])/R3;
      fyc=(fy[0]+fy[1]+fy[2])/R3;
      for(i=0;i<3;i++)
      { x[i]=x[i]-xc;
        y[i]=y[i]-yc;
        vx[i]=vx[i]-vxc;
        vy[i]=vy[i]-vyc;
      }
      momf=R0;  momv=R0;
      for(i=0;i<3;i++)
      { momf=momf+x[i]*fy[i]-y[i]*fx[i];
        momv=momv+x[i]*vy[i]-y[i]*vx[i];
      }
      mas=dpero*((x[1]-x[0])*(y[2]-y[0])-
                 (y[1]-y[0])*(x[2]-x[0]))/R6;
      iner=mas*(x[0]*x[0]+y[0]*y[0]
               +x[1]*x[1]+y[1]*y[1]
               +x[2]*x[2]+y[2]*y[2]);
      momv=momv*mas;
      vxc=vxc+fxc*dcstec/mas;
      vyc=vyc+fyc*dcstec/mas;
      xc=xc+dcstec*vxc;
      yc=yc+dcstec*vyc;
      omega=momv/iner+momf*dcstec/iner;
      c=COS(dcstec*omega);
      s=SIN(dcstec*omega);
      for(i=0;i<3;i++)
      { xnew=x[i]*c-y[i]*s;
        ynew=y[i]*c+x[i]*s;
        d1nccx[i2elto[i][ielem]]=xc+xnew;
        d1nccy[i2elto[i][ielem]]=yc+ynew;
        d1nvcx[i2elto[i][ielem]]=vxc-ynew*omega;
        d1nvcy[i2elto[i][ielem]]=vyc+xnew*omega;
        d1nmct[i2elto[i][ielem]]=mas;
      }
    }
  }
}

/********* PUBLIC ********************/
void Ysd(ydc,yde, ydn, ydo, ydpe, ydpn   /*** explicit solver of equations   ***/
         )
  YDC ydc; YDE yde; YDN ydn; YDO ydo; YDPE ydpe; YDPN ydpn;
{ DBL Ek,stprev;
  INT iprop,i,ihys, ibc;

  /* Solve Y elem. and nodes */
  for(ibc=0;ibc < ydpn->npnset;ibc++)		/* For Property Nodes ONLY */
  { /* YTN2MEC = -1 2D mechanical x,y d.o.f. node */
    Ysd2MEC(/* mechanical solver for 2D nodes with x,y d.o.f.  */
    ydn->nnopo,
    ydc->dcstec,ydc->dcgray, ydc->dctime, ydpn->npnfact, ibc,
    ydpn->d3pnfac,
    ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nfc[0],ydn->d2nfc[1],ydn->d1nmct,
    ydn->d2nvc[0],ydn->d2nvc[1],ydpn->d1pnax,ydpn->d1pnay,ydpn->d1pnxx,
    ydpn->d1pnxy,ydpn->d1pnyx,ydpn->d1pnyy,ydn->i1nopr,
    ydpn->i1pnfx, ydpn->i1pnfy
    );
  }
  for(iprop=0;iprop<ydpe->nprop;iprop++)
  { if((ydpe->i1ptyp[iprop])==(YTE2TRIRIG))		/* YTE2TRIRIG = -2 rigid triangle*/
    { Ysd2TRIRIG( /* mechanical solver rigid triangles  */
      yde->nelem,
      ydc->dcstec,ydpe->d1pero[iprop], iprop,
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nfc[0],ydn->d2nfc[1],ydn->d1nmct,
      ydn->d2nvc[0],ydn->d2nvc[1],yde->i1elpr,yde->i2elto
      );
  } }
  for(ihys=0;ihys<ydo->nohys;ihys++)  /* get history variables */
  { if(ydo->i1ohyt[ihys]==(YFLEK))
    { Ek=R0;
      for(i=0;i<ydn->nnopo;i++)
      {  Ek=Ek+RP5*(ydn->d2nvc[0][i]*ydn->d2nvc[0][i]+
                ydn->d2nvc[1][i]*ydn->d2nvc[1][i])*ydn->d1nmct[i];
      }
      stprev=MAXIM((EPSILON),(ABS(ydo->d1ohys[ihys])));
      if((ABS(R1-Ek/stprev))>=ydo->dohyp)
      { ydo->d1ohys[ihys]=Ek;
        ydo->d1ohyt[ihys]=ydc->dctime;
  } } }
}

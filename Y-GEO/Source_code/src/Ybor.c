/* Copyright (C) 2009, Queen Mary, University of London, Imperial College London
 * Implemented by Tomas Lukas
 * This code is provided as part of the book entitled "The Combined
 * Finite Discrete Element Method". It is distributed WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. Inclusion of a part or whole of this code into any other
 * commercial or research or other purpose code is not granted without author's
 * written explicit permission.
 */

/* File  Ybor.c */
#include "Yproto.h"

DBL Ybor2factor(dctime,npaf,Ts,Te,d1pres)
  DBL  dctime; INT  npaf; DBL  Ts; DBL  Te;
  DBL *d1pres;
{ DBL factor;
  DBL Tstep;
  INT Nstart;
  INT Nend;
  DBL Tstart;
  DBL Tend;

  /* find pressure factor */
  if(dctime == Ts)
  { factor = d1pres[0];
  }
  else if(dctime == Te)
  { factor = d1pres[npaf-1];
  }
  else if((dctime > Ts) && (dctime < Te))
  { Tstep = (Te - Ts)/(npaf-1);
    Tend = Ts;
    Nend = 0;
    while(Tend < dctime)
    { Tend = Tend + Tstep;
      Nend = Nend + 1;
    }
    Tstart = Tend - Tstep;
    Nstart = Nend - 1;
    if(dctime == Tend)
    { factor = d1pres[Nend];
    }
    else
    { factor = d1pres[Nstart]-
             ((d1pres[Nstart]-d1pres[Nend])*((dctime-Tstart)/Tstep));
  } }
  else
  { factor = R0;
  }
  return factor;
}

static void YborIntersection(r,xc,yc,
                             x0,y0,x1,y1,
                             dxi,dyi)
  DBL    r; DBL   xc; DBL yc;
  DBL   x0; DBL   y0; DBL x1; DBL y1;
  DBL *dxi; DBL *dyi;
{ DBL x0t,y0t,x1t,y1t;  /* temp. (local) coordinates    */
  DBL dx,dy,dr,D,Discr;
  DBL dx1,dx2,dy1,dy2;
  INT sgndy=1;

  /* equations are derived for circle of radius r and center (0, 0) */
  /* and line determined by two points (x0,y0) and (x1,y1)          */

  /* => set coordinate system to the center of circle               */
  x0t=x0-xc;
  y0t=y0-yc;
  x1t=x1-xc;
  y1t=y1-yc;

  dx=x1t-x0t;
  dy=y1t-y0t;
  dr=SQRT((dx*dx)+(dy*dy));
  if(dy<R0) sgndy=-1;
  D=(x0t*y1t)-(x1t*y0t);
  Discr=SQRT(((r*r)*(dr*dr))-(D*D));

  dx1=((D*dy)+(sgndy*dx*Discr))/(dr*dr);
  dx2=((D*dy)-(sgndy*dx*Discr))/(dr*dr);

  dy1=((-D*dx)+(DABS(dy)*Discr))/(dr*dr);
  dy2=((-D*dx)-(DABS(dy)*Discr))/(dr*dr);

  /* transform results back to global coord. system */
  dx1=dx1+xc;
  dy1=dy1+yc;
  dx2=dx2+xc;
  dy2=dy2+yc;

  dx=DABS(dx);
  dy=DABS(dy);
  if(dx>dy)
  { if(x1>x0)
    { if((x0<dx1)&&(dx1<x1))
      { *dxi=dx1;
        *dyi=dy1;
      }
      else
      { *dxi=dx2;
        *dyi=dy2;
    } }
    else
    { if((x1<dx1)&&(dx1<x0))
      { *dxi=dx1;
        *dyi=dy1;
      }
      else
      { *dxi=dx2;
        *dyi=dy2;
  } } }
  else
  { if(y1>y0)
    { if((y0<dy1)&&(dy1<y1))
      { *dxi=dx1;
        *dyi=dy1;
      }
      else
      { *dxi=dx2;
        *dyi=dy2;
    } }
    else
    { if((y1<dy1)&&(dy1<y0))
      { *dxi=dx1;
        *dyi=dy1;
      }
      else
      { *dxi=dx2;
        *dyi=dy2;
  } } }
}

static void Ybor2JOINTS(  /* joints on edges only  */
            nelem,dctime,
            jprop,
            nborh,nbpaf,d1bprs,dbbuf,
            d1bcax,d1bcay,d1bcbx,d1bcby,
            d1brad,d1bpaf,d1bpts,d1bpte,d1bvdt,
            d1nccx,d1nccy,d1ncix,d1nciy,
            d1nfcx,d1nfcy,d1nmct,
            i1elpr,i1nopr,i2elto
            )
  INT    nelem; DBL   dctime;
  INT    jprop;
  INT    nborh; INT     nbpaf; DBL  *d1bprs; DBL    dbbuf;
  DBL  *d1bcax; DBL  *d1bcay;  DBL  *d1bcbx; DBL  *d1bcby;
  DBL  *d1brad; DBL   *d1bpaf; DBL   *d1bpts; DBL   *d1bpte; DBL  *d1bvdt;
  DBL  *d1nccx; DBL  *d1nccy; DBL  *d1ncix; DBL  *d1nciy;
  DBL  *d1nfcx; DBL  *d1nfcy; DBL  *d1nmct;
  INT  *i1elpr; INT  *i1nopr; INT **i2elto;
{ INT ielem,iborh;
  INT i0,i1;
  DBL nx,ny;
  DBL L;
  DBL Ts,Te;        /* start and end time       */
  DBL u0,v0,u1,v1;
  DBL delta0,delta1,fac0,fac1;
  DBL dpres0,dpres1;
  DBL *d1rad;
  DBL *Lb;          /* lenght of each borhole   */
  DBL *rbx;
  DBL *rby;
  DBL *e1x;
  DBL *e1y;
  DBL *e2x;
  DBL *e2y;
  DBL *d1xmax;
  DBL *d1xmin;
  DBL *d1ymax;
  DBL *d1ymin;

  d1xmax=TalDBL1(nborh);
  d1xmin=TalDBL1(nborh);
  d1ymax=TalDBL1(nborh);
  d1ymin=TalDBL1(nborh);
  d1rad=TalDBL1(nborh);
  Lb=TalDBL1(nborh);
  rbx=TalDBL1(nborh);
  rby=TalDBL1(nborh);
  e1x=TalDBL1(nborh);
  e1y=TalDBL1(nborh);
  e2x=TalDBL1(nborh);
  e2y=TalDBL1(nborh);
  for(iborh=0;iborh<nborh;iborh++)
  { d1xmax[iborh]=R0;
    d1xmin[iborh]=R0;
    d1ymax[iborh]=R0;
    d1ymin[iborh]=R0;
    rbx[iborh]=R0;
    rby[iborh]=R0;
    Lb[iborh]=R0;
    e1x[iborh]=R0;
    e1y[iborh]=R0;
    e2x[iborh]=R0;
    e2y[iborh]=R0;
  }
  for(iborh=0;iborh<nborh;iborh++)
  { d1xmax[iborh]=MAXIM(d1bcax[iborh],d1bcbx[iborh])+d1brad[iborh]+dbbuf+EPSILON;
    d1xmin[iborh]=MINIM(d1bcax[iborh],d1bcbx[iborh])-d1brad[iborh]-dbbuf-EPSILON;
    d1ymax[iborh]=MAXIM(d1bcay[iborh],d1bcby[iborh])+d1brad[iborh]+dbbuf+EPSILON;
    d1ymin[iborh]=MINIM(d1bcay[iborh],d1bcby[iborh])-d1brad[iborh]-dbbuf-EPSILON;
    d1rad[iborh]=d1brad[iborh]+EPSILON;
    /* prepare local coord. system u,v (origin at point B; A lies on v axis)  */
    rbx[iborh]=d1bcax[iborh]-d1bcbx[iborh];
    rby[iborh]=d1bcay[iborh]-d1bcby[iborh];
    UnVec(&e2x[iborh],&e2y[iborh],rbx[iborh],rby[iborh]);
    Rot90(&e1x[iborh],&e1y[iborh],e2x[iborh],e2y[iborh]);
    V2DTranToLoc(&u0,&v0,rbx[iborh],rby[iborh],
                 e1x[iborh],e1y[iborh],e2x[iborh],e2y[iborh]);
    Lb[iborh]=v0;                   /* length of iborh  */
  }

  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==jprop)
    { if((i2elto[0][ielem]==i2elto[3][ielem])&& /* edge=joint with jprop*/
         (i2elto[1][ielem]==i2elto[2][ielem]))
      { i0=i2elto[0][ielem];    /* 1st node */
        i1=i2elto[1][ielem];    /* 2nd node */
        for(iborh=0; iborh<nborh; iborh++)
        { if(((d1xmin[iborh]<d1nccx[i0])&&(d1nccx[i0]<d1xmax[iborh])) &&
             ((d1ymin[iborh]<d1nccy[i0])&&(d1nccy[i0]<d1ymax[iborh])))
          { /* 1st node - transform to local coord. system  */
            V2DTranToLoc(&u0,&v0,(d1nccx[i0]-d1bcbx[iborh]),
                         (d1nccy[i0]-d1bcby[iborh]),e1x[iborh],e1y[iborh],
                         e2x[iborh],e2y[iborh]);
            delta0=DABS(u0/d1rad[iborh]);
            /* 2nd node - transform to local coord. system  */
            V2DTranToLoc(&u1,&v1,(d1nccx[i1]-d1bcbx[iborh]),
                         (d1nccy[i1]-d1bcby[iborh]),e1x[iborh],e1y[iborh],
                         e2x[iborh],e2y[iborh]);
            delta1=DABS(u1/d1rad[iborh]);
            if((delta0<=R1) && (delta1<=R1))    /* both nodes inside    */
            { L=Lb[iborh]-v0;
              Ts=d1bpts[iborh]+(L/d1bvdt[iborh]);   /* start time - i0  */
              Te=d1bpte[iborh]+(L/d1bvdt[iborh]);   /* end time - i0    */
              fac0=Ybor2factor(dctime,nbpaf,Ts,Te,d1bpaf);
              L=Lb[iborh]-v1;
              Ts=d1bpts[iborh]+(L/d1bvdt[iborh]);   /* start time - i1  */
              Te=d1bpte[iborh]+(L/d1bvdt[iborh]);   /* end time - i1    */
              fac1=Ybor2factor(dctime,nbpaf,Ts,Te,d1bpaf);
              dpres0=d1bprs[iborh]*fac0;
              dpres1=d1bprs[iborh]*fac1;
              nx=d1nccy[i1]-d1nccy[i0];	/* normal vector perpendicular  */
              ny=d1nccx[i0]-d1nccx[i1]; /* to the edge == length        */
              d1nfcx[i0]=d1nfcx[i0]-dpres0*nx/R3-dpres1*nx/R6;
              d1nfcy[i0]=d1nfcy[i0]-dpres0*ny/R3-dpres1*ny/R6;
              d1nfcx[i1]=d1nfcx[i1]-dpres0*nx/R6-dpres1*nx/R3;
              d1nfcy[i1]=d1nfcy[i1]-dpres0*ny/R6-dpres1*ny/R3;
  } } } } } } 

  FREE(d1xmax);
  FREE(d1xmin);
  FREE(d1ymax);
  FREE(d1ymin);
  FREE(d1rad);
  FREE(Lb);
  FREE(rbx);
  FREE(rby);
  FREE(e1x);
  FREE(e1y);
  FREE(e2x);
  FREE(e2y);
}
static void Yinterfluid2JOINTS(  /* broken joints or edge only  */
            nelem,dctime,
            jprop,
            nsour,nsdim,nspaf,nssaf,
            d1scsx,d1scsy,d1spaf,d1ssaf,
            d1spts,d1spte,d1svpr,d1sprs,
            d1ssir,dsbuf,
            d1nccx,d1nccy,d1ncix,d1nciy,
            d1nfcx,d1nfcy,d1nmct,
            i1elpr,i1nopr,i2elto
            )
  INT    nelem; DBL   dctime;
  INT    jprop;
  INT    nsour; INT    nsdim; INT    nspaf; INT    nssaf;
  DBL  *d1scsx; DBL  *d1scsy; DBL  *d1spaf; DBL  *d1ssaf;
  DBL  *d1spts; DBL  *d1spte; DBL  *d1svpr; DBL  *d1sprs;
  DBL  *d1ssir; DBL    dsbuf;
  DBL  *d1nccx; DBL  *d1nccy; DBL  *d1ncix; DBL  *d1nciy;
  DBL  *d1nfcx; DBL  *d1nfcy; DBL  *d1nmct;
  INT  *i1elpr; INT  *i1nopr; INT **i2elto;
{ INT ielem,isour,irep,nrep;
  INT i0,i1;
  DBL nx,ny;
  DBL drad0,drad1,dradx,drady;
  DBL delta0,delta1,fac0,fac1;
  DBL dpres0,dpres1;
  DBL dxi,dyi;  /* intersection point (circle-line)         */
  DBL dL;       /* length of pressure-ends at intersection  */
  DBL dL2;      /* length of the edge of joint              */
  DBL dA,dB;    /* dA=force at node inside; dB=node outside */
  DBL *d1Rad;   /* radii for each source                    */
  DBL *d1pa;    /* pressure amplitudes for each source      */
  DBL *d1xmax;
  DBL *d1xmin;
  DBL *d1ymax;
  DBL *d1ymin;

  d1Rad=TalDBL1(nsour);
  d1pa=TalDBL1(nsour);
  d1xmax=TalDBL1(nsour);
  d1xmin=TalDBL1(nsour);
  d1ymax=TalDBL1(nsour);
  d1ymin=TalDBL1(nsour);
  for(isour=0;isour<nsour;isour++)
  { d1Rad[isour]=R0;
    d1pa[isour]=R0;
    d1xmax[isour]=R0;
    d1xmin[isour]=R0;
    d1ymax[isour]=R0;
    d1ymin[isour]=R0;
  }
  for(isour=0;isour<nsour;isour++)
  { if((dctime>=d1spts[isour]) && (dctime<=d1spte[isour]))
    { d1Rad[isour]=d1svpr[isour]*(dctime-d1spts[isour]);
      fac0=Ybor2factor(dctime,nspaf,d1spts[isour],d1spte[isour],d1spaf);
      d1pa[isour]=d1sprs[isour]*fac0;
      d1xmax[isour]=d1scsx[isour]+d1Rad[isour]+dsbuf+EPSILON;
      d1xmin[isour]=d1scsx[isour]-d1Rad[isour]-dsbuf-EPSILON;
      d1ymax[isour]=d1scsy[isour]+d1Rad[isour]+dsbuf+EPSILON;
      d1ymin[isour]=d1scsy[isour]-d1Rad[isour]-dsbuf-EPSILON;
  } }

  for(ielem=0;ielem<nelem;ielem++)
  { if(((i1elpr[ielem]+YIPROPMAX)==jprop) ||    /* broken joint */
       ((i2elto[0][ielem]==i2elto[3][ielem])&&      /* or edge  */
        (i2elto[1][ielem]==i2elto[2][ielem])))
    { nrep=0;
      if(i1elpr[ielem]==jprop) nrep=1;      /* edge with jprop  */
      if((i1elpr[ielem]+YIPROPMAX)==jprop) nrep=2;
      for(irep=0;irep<nrep;irep++)
      { i0=i2elto[0+2*irep][ielem];
        i1=i2elto[1+2*irep][ielem];
        for(isour=0;isour<nsour;isour++)
        { if(((d1xmin[isour]<d1nccx[i0])&&(d1nccx[i0]<d1xmax[isour])) &&
             ((d1ymin[isour]<d1nccy[i0])&&(d1nccy[i0]<d1ymax[isour])))
          { /* 1st node */
            dradx=d1nccx[i0]-d1scsx[isour];
            drady=d1nccy[i0]-d1scsy[isour];
            drad0=SQRT((dradx*dradx)+(drady*drady));
            delta0=(drad0-d1ssir[isour])/(d1Rad[isour]-d1ssir[isour]+EPSILON);
            if(delta0<R0) delta0=R0;
            /* 2nd node */
            dradx=d1nccx[i1]-d1scsx[isour];
            drady=d1nccy[i1]-d1scsy[isour];
            drad1=SQRT((dradx*dradx)+(drady*drady));
            delta1=(drad1-d1ssir[isour])/(d1Rad[isour]-d1ssir[isour]+EPSILON);
            if(delta1<R0) delta1=R0;
            if((delta0<=R1) && (delta1<=R1))    /* both nodes inside    */
            { fac0=Ybor2factor(delta0,nssaf,R0,R1,d1ssaf);
              fac1=Ybor2factor(delta1,nssaf,R0,R1,d1ssaf);
              dpres0=d1pa[isour]*fac0;
              dpres1=d1pa[isour]*fac1;
              nx=d1nccy[i1]-d1nccy[i0];	/* normal vector perpendicular  */
              ny=d1nccx[i0]-d1nccx[i1]; /* to the edge == length        */
              d1nfcx[i0]=d1nfcx[i0]-dpres0*nx/R3-dpres1*nx/R6;
              d1nfcy[i0]=d1nfcy[i0]-dpres0*ny/R3-dpres1*ny/R6;
              d1nfcx[i1]=d1nfcx[i1]-dpres0*nx/R6-dpres1*nx/R3;
              d1nfcy[i1]=d1nfcy[i1]-dpres0*ny/R6-dpres1*ny/R3;
            }
            else if((delta0<R1) && (delta1>R1))  /* i0 inside           */
            { fac0=Ybor2factor(delta0,nssaf,R0,R1,d1ssaf);
              dpres0=d1pa[isour]*fac0;
              YborIntersection(d1Rad[isour],d1scsx[isour],d1scsy[isour],
                               d1nccx[i0],d1nccy[i0],d1nccx[i1],d1nccy[i1],
                               &dxi,&dyi);
              dL=SQRT(((dxi-d1nccx[i0])*(dxi-d1nccx[i0]))+
                      ((dyi-d1nccy[i0])*(dyi-d1nccy[i0])));
              nx=d1nccy[i1]-d1nccy[i0];	/* normal vector perpendicular  */
              ny=d1nccx[i0]-d1nccx[i1]; /* to the edge == length        */
              dL2=SQRT((nx*nx)+(ny*ny));
              dA=(dpres0*dL/R2)*(R1-(dL/(R3*dL2)));
              dB=dpres0*dL*dL/(R6*dL2);
              nx=nx/dL2;                /* unit normal vector           */
              ny=ny/dL2;
              d1nfcx[i0]=d1nfcx[i0]-dA*nx;
              d1nfcy[i0]=d1nfcy[i0]-dA*ny;
              d1nfcx[i1]=d1nfcx[i1]-dB*nx;
              d1nfcy[i1]=d1nfcy[i1]-dB*ny;
            }
            else if((delta0>R1) && (delta1<R1))  /* i1 inside           */
            { fac1=Ybor2factor(delta1,nssaf,R0,R1,d1ssaf);
              dpres1=d1pa[isour]*fac1;
              YborIntersection(d1Rad[isour],d1scsx[isour],d1scsy[isour],
                               d1nccx[i0],d1nccy[i0],d1nccx[i1],d1nccy[i1],
                               &dxi,&dyi);
              dL=SQRT(((dxi-d1nccx[i1])*(dxi-d1nccx[i1]))+
                      ((dyi-d1nccy[i1])*(dyi-d1nccy[i1])));
              nx=d1nccy[i1]-d1nccy[i0];	/* normal vector perpendicular  */
              ny=d1nccx[i0]-d1nccx[i1]; /* to the edge == length        */
              dL2=SQRT((nx*nx)+(ny*ny));
              dA=(dpres0*dL/R2)*(R1-(dL/(R3*dL2)));
              dB=dpres0*dL*dL/(R6*dL2);
              nx=nx/dL2;                /* unit normal vector           */
              ny=ny/dL2;
              d1nfcx[i0]=d1nfcx[i0]-dB*nx;
              d1nfcy[i0]=d1nfcy[i0]-dB*ny;
              d1nfcx[i1]=d1nfcx[i1]-dA*nx;
              d1nfcy[i1]=d1nfcy[i1]-dA*ny;
  } } } } } }

  FREE(d1Rad);
  FREE(d1pa);
  FREE(d1xmax);
  FREE(d1xmin);
  FREE(d1ymax);
  FREE(d1ymin);
}

/*********************PUBLIC********************************************************/
void Ybor( ydc, yde, ydn, ydb, yds,            /*** pressure on borholes; source ***/
           ydpe, ydpj, ydpn)
 YDC ydc; YDE yde; YDN ydn; YDB ydb; YDS yds;
 YDPE ydpe; YDPJ ydpj; YDPN ydpn;
{ INT jprop,iborh,isour;
  DBL tmin,tmax,rmax;
  DBL vmin;

  if(ydb->nborh>0)
  { tmin=ydb->d1bpts[0];
    tmax=ydb->d1bpte[0];
    rmax=ydb->d1brad[0];
    vmin=ydb->d1bvdt[0];
    for(iborh=0; iborh<ydb->nborh; iborh++)
    { tmin=MINIM(tmin,ydb->d1bpts[iborh]);
      tmax=MAXIM(tmax,ydb->d1bpte[iborh]);
      rmax=MAXIM(rmax,ydb->d1brad[iborh]);
      vmin=MINIM(vmin,ydb->d1bvdt[iborh]);
    }
    tmax=tmax+((ydb->dblmax+rmax)/vmin);	/* max end time */
    if((ydc->dctime>=tmin) && (ydc->dctime<=tmax))
    { for(jprop=0;jprop<ydpj->npjset;jprop++)
      { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
        { Ybor2JOINTS(   /* apply pressure on edge of borehole  */
          yde->nelem,ydc->dctime,
          (jprop+ydpe->nprop),
          ydb->nborh,ydb->nbpaf,ydb->d1bprs,ydb->dbbuf,
          ydb->d2bca[0],ydb->d2bca[1],ydb->d2bcb[0],ydb->d2bcb[1],
          ydb->d1brad,ydb->d1bpaf,ydb->d1bpts,ydb->d1bpte,ydb->d1bvdt,
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],
          ydn->d2nfc[0],ydn->d2nfc[1],ydn->d1nmct,
          yde->i1elpr,ydn->i1nopr,yde->i2elto
          );
  } } } }

  if(yds->nsour>0)
  { tmin=R0;
    tmax=R0;
    for(isour=0;isour<yds->nsour;isour++)
    { tmin=MINIM(tmin,yds->d1spts[isour]);
      tmax=MAXIM(tmax,yds->d1spte[isour]);
    }
    if((ydc->dctime>=tmin) && (ydc->dctime<=tmax))
    { for(jprop=0;jprop<ydpj->npjset;jprop++)
      { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
        { Yinterfluid2JOINTS(   /* broken joints or edge only  */
          yde->nelem,ydc->dctime,
          (jprop+ydpe->nprop),
          yds->nsour,yds->nsdim,yds->nspaf,yds->nssaf,
          yds->d2scs[0],yds->d2scs[1],yds->d1spaf,yds->d1ssaf,
          yds->d1spts,yds->d1spte,yds->d1svpr,yds->d1sprs,
          yds->d1ssir,yds->dsbuf,
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],
          ydn->d2nfc[0],ydn->d2nfc[1],ydn->d1nmct,
          yde->i1elpr,ydn->i1nopr,yde->i2elto
          );
  } } } }
}



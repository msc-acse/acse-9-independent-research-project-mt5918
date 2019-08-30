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
/* File  Ycd.c */
#include "Yproto.h"
#define G2GETFORS 0
#define G2GETPRES 1

static void Ycdvel2TRI( /* find maximum squered velocity */
             nelem,
             iprop,
            d0mxv2,
            d1nvcx,d1nvcy,i1elpr,i2elto
            ) 
  INT   nelem;
  INT   iprop;
  DBL *d0mxv2;
  DBL *d1nvcx; DBL *d1nvcy; INT *i1elpr; INT **i2elto;
{ INT ielem;
  DBL x0,x1,x2,y0,y1,y2;

  if(nelem<1)return;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { /* find max velocity */
      x0=d1nvcx[i2elto[0][ielem]];
      x1=d1nvcx[i2elto[1][ielem]];
      x2=d1nvcx[i2elto[2][ielem]];
      y0=d1nvcy[i2elto[0][ielem]];
      y1=d1nvcy[i2elto[1][ielem]];
      y2=d1nvcy[i2elto[2][ielem]];
      (*d0mxv2)=MAXIM((*d0mxv2),(x0*x0+y0*y0));
      (*d0mxv2)=MAXIM((*d0mxv2),(x1*x1+y1*y1));
      (*d0mxv2)=MAXIM((*d0mxv2),(x2*x2+y2*y2));
} } }

static void Ycdcor2TRI(  /* find element coordinates */
             nelem,
             iprop,
            d1erad,d1eccx,d1eccy,d1nccx,d1nccy,
            i1nobf,i1elpr,i2elto
            ) 
  INT   nelem;
  INT   iprop;
  DBL *d1erad; DBL *d1eccx; DBL *d1eccy; DBL  *d1nccx; DBL *d1nccy; 
  INT *i1elpr; INT *i1nobf; INT **i2elto;
{ INT ielem,i;
  DBL x0,x1,x2,y0,y1,y2,xc,yc;

  if(nelem<1)return;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { for(i=0;i<3;i++)
      { if(i1nobf[i2elto[i][ielem]]>0)
        { x0=d1nccx[i2elto[0][ielem]];
          x1=d1nccx[i2elto[1][ielem]];
          x2=d1nccx[i2elto[2][ielem]];
          y0=d1nccy[i2elto[0][ielem]];
          y1=d1nccy[i2elto[1][ielem]];
          y2=d1nccy[i2elto[2][ielem]];
          xc=(x0+x1+x2)/R3;
          yc=(y0+y1+y2)/R3;
          d1eccx[ielem]=xc;
          d1eccy[ielem]=yc;
          d1erad[ielem]=(xc-x0)*(xc-x0)+(yc-y0)*(yc-y0);
          d1erad[ielem]=MAXIM(d1erad[ielem],
          ((xc-x1)*(xc-x1)+(yc-y1)*(yc-y1)));
          d1erad[ielem]=MAXIM(d1erad[ielem],
          ((xc-x2)*(xc-x2)+(yc-y2)*(yc-y2)));
          d1erad[ielem]=SQRT(d1erad[ielem]);
          break;
} } } } }
 
static void Ycdprocess(
     melem  ,micoup,nelem , nprop,
     dcstec ,diezon,d0iedi,
     i0iecff,
     d1nccx ,d1nccy,d1nccz,d1nvcx,d1nvcy,
     d1nvcz ,i1elcf,i1elpr,i1iecn,i1iect,
     i1nobf,i1ptyp ,i2elto,i2eljp,i1elty
     )
  INT    melem; INT micoup; INT   nelem; INT   nprop;
  DBL   dcstec; DBL diezon; DBL *d0iedi;
  INT *i0iecff;
  DBL  *d1nccx; DBL *d1nccy; DBL *d1nccz;  DBL  *d1nvcx; DBL *d1nvcy;
  DBL  *d1nvcz; INT *i1elcf; INT *i1elpr;  INT  *i1iecn; INT *i1iect;
  INT  *i1nobf; INT *i1ptyp; INT **i2elto; INT **i2eljp; INT *i1elty;
{ INT ncelx,ncely,ncelz;  /* total number of x, y, z cells   */
  INT nelemd;             /* twice total number of elements  */
  DBL diam;               /* maximum diameter of element     */
  DBL dmxv2;              /* maximum displacement increment  */
  INT icoup,jcoup;        /* couple                          */
  INT icouco;             /* couple's contactor              */
  INT icouta;             /* couple's target                 */
  INT ielem;                /* element                       */
  INT ielemx;               /* element assigned to cell x    */
  INT ielemy;               /* element assigned to cell y    */
  INT ielemz;               /* element assigned to cell z    */
  INT ihx,ihy,ihz;          /* x, y, z  head of a   list     */
  INT iminx,iminy,iminz;    /* space boundaries              */
  INT imaxx,imaxy,imaxz;    /* space boundaries              */
  INT iprop;                /* element property id           */  
  INT ix,iy,iz;             /* x, y, z cell                  */
  INT jelemx;               /* element assigned to cell x    */
  DBL *d1eccx;     /* element coordinate current x           */
  DBL *d1eccy;     /* element coordinate current y           */
  DBL *d1eccz;     /* element coordinate current z           */
  DBL *d1erad;     /* element radius squered                 */
  INT *i1cnx;      /* contactor next  x                      */
  INT *i1cny;      /* contactor next  y                      */
  INT *i1cnz;      /* contactor next  z                      */
  INT *i1cfz;      /* contactor first z                      */
  INT *i1eccx;     /* element coordinate current x           */
  INT *i1eccy;     /* element coordinate current y           */
  INT *i1eccz;     /* element coordinate current z           */
  INT  i1heax[5];  /* heads of 5 connected lists for x cells */
  INT  i1heay[5];  /* heads of 5 connected lists for y cells */
  INT  i1heaz[2];  /* heads of 2 connected lists for z cells */
  INT *i2cfx[2];   /* contactor first x                      */
  INT *i2cfy[2];   /* contactor first y                      */
  INT i,j;		   /* loop variables                         */
  INT isintern;	   /* contact via internal faces	         */

  /* Process contact detection if specified by input data */
  if(micoup>0)
  { nelemd=2*nelem;
    /* initialise data if needed i.e. if(i0iecff ==default) */
    if((*i0iecff)==-2)
    { for(ielem=0;ielem<melem;ielem++)
      { i1elcf[ielem]=-1;
      }
      for(icoup=0;icoup<micoup;icoup++)
      { i1iecn[icoup]=icoup+1;
      }
      (*i0iecff)=0;
      i1iecn[micoup-1]=-1;
    }
    /* find maximum  velocity */
    dmxv2=R0;
    for(iprop=0;iprop<nprop;iprop++)   
    { if(((i1ptyp[iprop])==(YTE2TRIELS))||
         ((i1ptyp[iprop])==(YTE2TRIRIG)))
      { Ycdvel2TRI( 
         nelem,
         iprop,
        &dmxv2,
        d1nvcx,d1nvcy,i1elpr,i2elto 
        );
    } } 
    (*d0iedi)=(*d0iedi)+R2*dcstec*SQRT(dmxv2);
    /* if not time for contact detection, return */
    if((*d0iedi)>=diezon)
    { /* unmark new couples */
      for(icoup=0;icoup<micoup;icoup++)
      { if(i1iect[icoup]<0) i1iect[icoup]=-1-i1iect[icoup];
      }
      /* find element coordinates */
      (*d0iedi)=R0;
      d1erad=TalDBL1(nelem);    /* element radius to the power 2     */
      d1eccz=TalDBL1(nelem);    /* element coordinate current z      */
      d1eccy=TalDBL1(nelem);    /* element coordinate current y      */
      d1eccx=TalDBL1(nelem);    /* element coordinate current x      */
      i1eccz=TalINT1(nelem);    /* element coordinate current z      */
      i1eccy=TalINT1(nelem);    /* element coordinate current y      */
      i1eccx=TalINT1(nelem);    /* element coordinate current x      */
      for(ielem=0;ielem<nelem;ielem++)
      { d1erad[ielem]=R0;
        d1eccx[ielem]=R0;
        d1eccy[ielem]=R0;
        d1eccz[ielem]=R0;
      }
      for(iprop=0;iprop<nprop;iprop++)
      { if(((i1ptyp[iprop])==(YTE2TRIELS))||
           ((i1ptyp[iprop])==(YTE2TRIRIG)))
        { Ycdcor2TRI( 
          nelem ,
          iprop ,
          d1erad,d1eccx,d1eccy,d1nccx,
          d1nccy,i1nobf,i1elpr,i2elto
          );
      } }
      /* find maximum diameter */
      diam=d1erad[0];
      for(ielem=1;ielem<nelem;ielem++)
      { diam=MAXIM(diam,d1erad[ielem]);        
      }
      if(diam<EPSILON)diam=BEPSILON;
      diam=R2*(diam+diezon);
      /* intigerise element coordinates */
      for(ielem=0;ielem<nelem;ielem++)
      { i1eccx[ielem]=(INT)(d1eccx[ielem]/diam);
        i1eccy[ielem]=(INT)(d1eccy[ielem]/diam);
        i1eccz[ielem]=(INT)(d1eccz[ielem]/diam);
      }
      /* find space boundaries */
      iminx=i1eccx[0];
      imaxx=i1eccx[0];
      iminy=i1eccy[0];
      imaxy=i1eccy[0];
      iminz=i1eccz[0];
      imaxz=i1eccz[0];
      for(ielem=1;ielem<nelem;ielem++)
      { iminx=MINIM(iminx,i1eccx[ielem]);
        imaxx=MAXIM(imaxx,i1eccx[ielem]);
        iminy=MINIM(iminy,i1eccy[ielem]);
        imaxy=MAXIM(imaxy,i1eccy[ielem]);
        iminz=MINIM(iminz,i1eccz[ielem]);
        imaxz=MAXIM(imaxz,i1eccz[ielem]);
      }
      iminx=iminx-1;
      iminy=iminy-1;
      iminz=iminz-1;
      imaxx=imaxx+2;
      imaxy=imaxy+2;
      imaxz=imaxz+1;
      /* normalise coordinates */
      for(ielem=0;ielem<nelem;ielem++)
      { i1eccx[ielem]=i1eccx[ielem]-iminx;
        i1eccy[ielem]=i1eccy[ielem]-iminy;
        i1eccz[ielem]=i1eccz[ielem]-iminz;
      }
      ncelx=imaxx-iminx;
      ncely=imaxy-iminy;
      ncelz=imaxz-iminz;
      /* allocate memory  */
      i1cfz=TalINT1(ncelz);    /* contactor first z                   */
      i2cfy[0]=TalINT1(ncely); /* contactor first y  (iz-1)           */
      i2cfy[1]=TalINT1(ncely); /* contactor first y  (iz  )           */
      i2cfx[0]=TalINT1(ncelx); /* contactor first x  (iz-1,iy-1)      */
      i2cfx[1]=TalINT1(ncelx); /* contactor first x  (iz-1,iy  )      */
      i1cnz=TalINT1(nelem);    /* contactor next z                    */
      i1cny=TalINT1(nelem);    /* contactor next y                    */
      i1cnx=TalINT1(nelem);    /* contactor next x                    */
      /* assume no contactors at any cell */
      for(iz=0;iz<ncelz;iz++)
      { i1cfz[iz]=-1;  
      }
      for(ihy=0;ihy<2;ihy++)
      { for(iy=0;iy<ncely;iy++)
        { i2cfy[ihy][iy]=-1;
      } }  
      for(ihx=0;ihx<2;ihx++)
      { for(ix=0;ix<ncelx;ix++)
        { i2cfx[ihx][ix]=-1;
      } }
      /* assign all contactors to z-cells */
      for(ielem=0;ielem<nelem;ielem++) 
      { if(d1erad[ielem]>R0)
        { i1cnz[ielem]=i1cfz[i1eccz[ielem]];
          i1cfz[i1eccz[ielem]]=ielem;
      } }
      /* scan all loaded z cells */
      for(ielem=0;ielem<nelem;ielem++)
      { iz=i1eccz[ielem];
        if(i1cfz[iz]<nelem)
        { i1heaz[0]=i1cfz[iz];
          i1heaz[1]=i1cfz[iz-1];
          if(i1heaz[1]>nelem) i1heaz[1]=i1heaz[1]-nelemd;
          i1cfz[iz]=i1cfz[iz]+nelemd;
          /* load elements from cells iz & iz-1  onto y cells */
          for(ihz=0;ihz<2;ihz++)
          { ielemz=i1heaz[ihz];
            while(ielemz>=0)
            { i1cny[ielemz]=i2cfy[ihz][(i1eccy[ielemz])];
              i2cfy[ihz][(i1eccy[ielemz])]=ielemz;
              ielemz=i1cnz[ielemz];
          } }
          /* scan all loaded y cells */
          ielemz=i1heaz[0];
          while(ielemz>=0)
          { iy=i1eccy[ielemz];
            if(i2cfy[0][iy]<nelem)
            { i1heay[0]=i2cfy[0][iy];
              i1heay[1]=i2cfy[0][iy-1];
              i1heay[2]=i2cfy[1][iy+1];
              i1heay[3]=i2cfy[1][iy];
              i1heay[4]=i2cfy[1][iy-1];
              if(i1heay[1]>nelem)i1heay[1]=i1heay[1]-nelemd;
              i2cfy[0][iy]=i2cfy[0][iy]+nelemd;
              /* load elements from y cells onto x cells */
              ihx=0;
              for(ihy=0;ihy<5;ihy++)
              { if(ihy>0)ihx=1;
                ielemy=i1heay[ihy];
                while(ielemy>=0)
                { i1cnx[ielemy]=i2cfx[ihx][(i1eccx[ielemy])];
                  i2cfx[ihx][(i1eccx[ielemy])]=ielemy;
                  ielemy=i1cny[ielemy];
              } }
              /* scan all loaded x cells */
              ielemy=i1heay[0];
              while(ielemy>=0)
              { ix=i1eccx[ielemy];
                if(i2cfx[0][ix]<nelem)
                { i1heax[0]=i2cfx[0][ix];
                  i1heax[1]=i2cfx[0][ix-1];
                  i1heax[2]=i2cfx[1][ix+1];
                  i1heax[3]=i2cfx[1][ix];
                  i1heax[4]=i2cfx[1][ix-1];
                  if(i1heax[1]>nelem)i1heax[1]=i1heax[1]-nelemd;
                  i2cfx[0][ix]=i2cfx[0][ix]+nelemd;
                  /* detect contacts for cell (ix,iy,iz) */
                  ielemx=i1heax[0];
                  while(ielemx>=0)
                  { for(ihx=0;ihx<5;ihx++)
                    { jelemx=i1heax[ihx];
                      while(jelemx>=0)
                      { if((ihx!=0)||(ielemx>jelemx))
                        { icouco=MAXIM(ielemx,jelemx);
                          icouta=MINIM(ielemx,jelemx);
                          /* exclude old couple */
                          icoup=i1elcf[icouco];
                          jcoup=icoup;
					      isintern=0;
						  if(i1elty[i2eljp[0][icouta]]==0 && i1elty[i2eljp[1][icouta]]==0
						  && i1elty[i2eljp[2][icouta]]==0 && i1elty[i2eljp[0][icouco]]==0
						  && i1elty[i2eljp[1][icouco]]==0 && i1elty[i2eljp[2][icouco]]==0)
						  { isintern=1;			// both are internal elements
						  }
						  else
						  { /* discontinuous mesh */
							for(i=0;i<3;i++)
						    { if(i1elty[i2eljp[i][icouco]]==0)
							  { for(j=0;j<3;j++)
								{ if(i1elty[i2eljp[j][icouta]]==0)
								  { if(i2eljp[i][icouco]==i2eljp[j][icouta])
									{ isintern=1; break;	  // connect not by fracture interface
							} } } } }
							/* continuous mesh */
							for(i=0;i<3;i++)
							{ for(j=0;j<3;j++)
							  { if(i2elto[i][icouco]==i2elto[j][icouta])
								{ isintern=1; break;
							} } }
						  }
						  /* contact via fracture or boundary */
						  if(isintern==0)
						  { while((icoup>=0)&&(icouta>=0))  // eliminate those uncontacted elements which previously contact
                            { if(i1iect[icoup]==icouta)
                              { if(((d1eccx[icouco]-d1eccx[icouta])*
                                    (d1eccx[icouco]-d1eccx[icouta])+
                                    (d1eccy[icouco]-d1eccy[icouta])*
                                    (d1eccy[icouco]-d1eccy[icouta])+
                                    (d1eccz[icouco]-d1eccz[icouta])*
                                    (d1eccz[icouco]-d1eccz[icouta]))>
                                    ((d1erad[icouco]+d1erad[icouta]+diezon)*
                                    (d1erad[icouco]+d1erad[icouta]+diezon)))
                                { i1iecn[jcoup]=i1iecn[icoup];
                                  if(jcoup==icoup) i1elcf[icouco]=i1iecn[icoup];
                                  i1iecn[icoup]=(*i0iecff);
                                  (*i0iecff)=icoup;
                                }
                                icouta=-1;
                              }
                              jcoup=icoup;
                              icoup=i1iecn[icoup];
                            }
                            if(icouta>=0)  /* add new couple if close */
                            { if(((d1eccx[icouco]-d1eccx[icouta])*
                                  (d1eccx[icouco]-d1eccx[icouta])+
                                  (d1eccy[icouco]-d1eccy[icouta])*
                                  (d1eccy[icouco]-d1eccy[icouta])+
                                  (d1eccz[icouco]-d1eccz[icouta])*
                                  (d1eccz[icouco]-d1eccz[icouta]))<
                                  ((d1erad[icouco]+d1erad[icouta]+diezon)*
                                  (d1erad[icouco]+d1erad[icouta]+diezon)))
                              { icoup=(*i0iecff);
                                if(icoup>=0)
                                { (*i0iecff)=i1iecn[icoup];
                                  i1iecn[icoup]=i1elcf[icouco];
                                  i1elcf[icouco]=icoup;
                                  i1iect[icoup]=-1-icouta;/* mark new */
                                }
                                else
                                { CHRw(stderr,"Ycd: to small meicoc  ");
                                  CHRwcr(stderr); exit(1);
						} } } } }
                        jelemx=i1cnx[jelemx];
                    } }
                    ielemx=i1cnx[ielemx];
                  } }
                  ielemy=i1cny[ielemy];
                }
                /* unload elements from x cells */
                ihx=0;
                for(ihy=0;ihy<5;ihy++)
                { if(ihy>0)ihx=1;
                  ielemy=i1heay[ihy];
                  while(ielemy>=0)
                  { i2cfx[ihx][(i1eccx[ielemy])]=-1;
                    ielemy=i1cny[ielemy];
              } } }
              ielemz=i1cnz[ielemz];
            }
            /* unload elements from y cells */
            for(ihz=0;ihz<2;ihz++)
            { ielemz=i1heaz[ihz];
              while(ielemz>=0)
              { i2cfy[ihz][(i1eccy[ielemz])]=-1;
                ielemz=i1cnz[ielemz];
        } } } }
      /* free memory  */
      FREE(i1cnx);
      FREE(i1cny);
      FREE(i1cnz);
      FREE(i2cfx[1]);
      FREE(i2cfx[0]);
      FREE(i2cfy[1]);  
      FREE(i2cfy[0]);
      FREE(i1cfz);
      FREE(i1eccx);   
      FREE(i1eccy);
      FREE(i1eccz);
      FREE(d1eccx);   
      FREE(d1eccy);
      FREE(d1eccz);
      FREE(d1erad);
} } }

/*********************PUBLIC********************************************************/
void Ycd(  ydc,  yde, ydj, ydi,  ydn, ydp   /***  contact detection  ***/
        )
  YDC ydc; YDE yde; YDJ ydj; YDI ydi; YDN ydn; YDP ydp;
{ INT i=0;
  if(ydn->nnopo>0)
  if(ydn->nnodim>1)
  { Ycdprocess( /* Y elements contact detection */
    yde->melem  ,ydi->micoup,yde->nelem ,ydp->nprop ,
    ydc->dcstec ,ydi->diezon,&ydi->diedi,
    &ydi->iiecff,
    ydn->d2ncc[0] ,ydn->d2ncc[1],DBL1NULL,ydn->d2nvc[0],ydn->d2nvc[1],
    DBL1NULL ,yde->i1elcf,yde->i1elpr,ydi->i1iecn,ydi->i1iect,
    ydn->i1nobf ,ydp->i1ptyp ,yde->i2elto, yde->i2eljp, yde->i1elty);
  }
}

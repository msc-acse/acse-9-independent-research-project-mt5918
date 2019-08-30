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
/* File   Ymd.c */
#include "Yproto.h"  
/**************GRAINS G2 R=0.5**********/
/* File   Ymd.c */
#include "Yproto.h"  
/**************MESH GRAINS***********/ 

/**************JOINT ELEMENTS***********/ 
static void Yjd2TRIANGLE(  /* create joints */
             melem, mnopo,nelest,nnopst,
             iprop,ipropj,
            n0elem,n0nopo,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,d1sdel,
            i1elpr,i1jnef,i1jnen,i1nobf,i1nopr,
            i2elto,i2eljp,i1elty
            )
  INT    melem; INT   mnopo; INT  nelest; INT nnopst;
  INT    iprop; INT  ipropj;
  INT  *n0elem; INT *n0nopo;
  DBL  *d1nccx; DBL *d1nccy; DBL *d1ncix; DBL *d1nciy; DBL *d1nvcx;
  DBL  *d1nvcy; DBL *d1sdel; 
  INT  *i1elpr; INT *i1jnef; INT *i1jnen; INT *i1nobf; INT *i1nopr;
  INT **i2elto; INT **i2eljp; INT *i1elty;
{ INT nelem,nnopo,nelecur;
  INT i,j,k,in,jn,kn,ijnew,ielem,jelem,kelem;
  DBL nx,ny;

  nelem=(*n0elem);
  nnopo=(*n0nopo);
  for(ielem=0;ielem<nelest;ielem++)
  { if(i1elpr[ielem]==iprop)
    { for(i=0;i<3;i++)
      { j=i+1; if(j>2)j=0;
        in=i2elto[i][ielem];
        jn=i2elto[j][ielem];
        kn=MAXIM(in,jn);
        if(nnopo>=mnopo)      /* create new node */
        { CHRw(stderr,"Yjd: MNOPO too small"); CHRwcr(stderr);
          exit(1);
        }
        d1nccx[nnopo]=d1nccx[in];
        d1nccy[nnopo]=d1nccy[in];
        d1ncix[nnopo]=d1ncix[in];
        d1nciy[nnopo]=d1nciy[in];
        d1nvcx[nnopo]=d1nvcx[in];
        d1nvcy[nnopo]=d1nvcy[in];
        i1nopr[nnopo]=i1nopr[in];
        i1nobf[nnopo]=1;
        nnopo=nnopo+1;
        /* check if joint element already existent */
        ijnew=i1jnef[kn];
        while(ijnew>=0)
        { if((i2elto[2][ijnew]==in)&&(i2elto[3][ijnew]==jn))
          { i1elpr[ijnew]=MAXIM(i1elpr[ijnew],ipropj);
            if(i<2)
            { i2elto[2][ijnew]=nnopo-1;
              i2elto[3][ijnew]=nnopo;
            }
            else
            { i2elto[2][ijnew]=nnopo-1;
              i2elto[3][ijnew]=nnopo-3;
            }
            i2eljp[1][ijnew]=ielem;
			ijnew=-100;
          }
          else
          { ijnew=i1jnen[ijnew-nelest];
        } }
        /* create new mid-edge joint node */
        if(ijnew>(-10))
        { if(nelem>=melem)
          { CHRw(stderr,"Yjd: MELEM too small"); CHRwcr(stderr);
            exit(1);
          }
          if(i<2)
          { i2elto[1][nelem]=nnopo;
            i2elto[0][nelem]=nnopo-1;
          }
          else
          { i2elto[1][nelem]=nnopo-3;
            i2elto[0][nelem]=nnopo-1;
          }
          i2elto[2][nelem]=jn;
          i2elto[3][nelem]=in;
		  i2eljp[0][nelem]=ielem;
          i1elpr[nelem]=ipropj;
          i1jnen[nelem-nelest]=i1jnef[kn];
		  i1elty[nelem]=0;
          i1jnef[kn]=nelem;
          nelem=nelem+1;
      } }
      for(i=0;i<3;i++) /* detach element */
      { i2elto[i][ielem]=nnopo-3+i;
      }
  } }
  /* reorganise joint elements */	  // added by Qinghua Lei
  nelecur=nelem;	// current number of elements
  for(ielem=nelest;ielem<nelem;ielem++)
  { if(d1sdel!=DBL1NULL) d1sdel[ielem]=R0;
	if(i2elto[2][ielem]<nnopst && i2elto[2][ielem]>-1)
	{ // checking joint whether overlaps with another one
	  for(jelem=ielem+1;jelem<nelecur;jelem++)
	  { if(i2elto[2][jelem]<nnopst)
		{ if(DABS(d1ncix[i2elto[0][ielem]]-d1ncix[i2elto[1][jelem]])<EPSILON
		  && DABS(d1ncix[i2elto[1][ielem]]-d1ncix[i2elto[0][jelem]])<EPSILON
		  && DABS(d1nciy[i2elto[0][ielem]]-d1nciy[i2elto[1][jelem]])<EPSILON
		  && DABS(d1nciy[i2elto[1][ielem]]-d1nciy[i2elto[0][jelem]])<EPSILON)
		  { // sort the order -- upper side first
			nx=d1nciy[i2elto[1][ielem]]-d1nciy[i2elto[0][ielem]];
			ny=d1ncix[i2elto[0][ielem]]-d1ncix[i2elto[1][ielem]];
			if((ny>EPSILON)||(DABS(ny)<=EPSILON&&nx>R0))
			{ i2elto[2][ielem]=i2elto[0][ielem];
			  i2elto[3][ielem]=i2elto[1][ielem];
			  i2elto[0][ielem]=i2elto[0][jelem];
			  i2elto[1][ielem]=i2elto[1][jelem];
			  i2eljp[1][ielem]=i2eljp[0][ielem];
			  i2eljp[0][ielem]=i2eljp[0][jelem];
			}
			else
			{ i2elto[2][ielem]=i2elto[0][jelem];
			  i2elto[3][ielem]=i2elto[1][jelem];
			  i2eljp[1][ielem]=i2eljp[0][jelem];
			}
			i1elty[ielem]=2;
			break;
	  } } }
	  if(i1elty[ielem]==2)		// joint elements on pre-exisitng fractures
	  { // i1elpr[ielem]=i1elpr[ielem]-YIPROPMAX;
		for(kelem=jelem;kelem<nelecur;kelem++)
		{ for(i=0;i<2;i++)
		  { i2elto[2*i][kelem]=i2elto[2*i][kelem+1];
			i2elto[2*i+1][kelem]=i2elto[2*i+1][kelem+1];
			i2eljp[i][kelem]=i2eljp[i][kelem+1];
		  }
		  i1elpr[kelem]=i1elpr[kelem+1];
		  i1elty[kelem]=i1elty[kelem+1];
		}
		for(i=0;i<2;i++)
		{ i2elto[2*i][nelecur]=-1;
		  i2elto[2*i+1][nelecur]=-1;
		  i2eljp[i][nelecur]=-1;
		}
		i1elty[nelecur]=-1;
		nelecur=nelecur-1;
	  }
	  else if(i1elty[ielem]==0)	  // joint elements on model boundaries
	  { i2elto[3][ielem]=i2elto[0][ielem];
		i2elto[2][ielem]=i2elto[1][ielem];
		i1elty[ielem]=1;
		i2eljp[1][ielem]=i2eljp[0][ielem];
  } } }
  nelem=nelecur;
  /* link joints to triangles */
  for(ielem=nelest;ielem<nelem;ielem++)		  // boundaries
  { if(i2elto[3][ielem]==i2elto[0][ielem]
	&& i2elto[2][ielem]==i2elto[1][ielem])
	{ jelem=i2eljp[0][ielem];
	  for(j=0;j<3;j++)
	  { k=j+1; if(k>2)k=0;
		if(i2elto[0][ielem]==i2elto[j][jelem]
		&& i2elto[1][ielem]==i2elto[k][jelem])
		{ i2eljp[j][jelem]=ielem; break;
	} } }
	else
	{ for(i=0;i<2;i++)
	  { jelem=i2eljp[i][ielem];
		for(j=0;j<3;j++)
		{ k=j+1; if(k>2)k=0;
		  if(i2elto[i*2][ielem]==i2elto[j][jelem]
		  && i2elto[i*2+1][ielem]==i2elto[k][jelem])
		  { if(i1elty[ielem]>0)				// fracture joints
			{ i2eljp[j][jelem]=ielem;
			}
			else							// unbroken joints
			{ i2eljp[j][jelem]=ielem;
		    }
			break;
	} } } }
  }
  /* update number of nodes and elements */
  (*n0nopo)=nnopo;
  (*n0elem)=nelem;
}

/**************MESH ELEMENTS***********/ 
static void Ymd2TRIANGLE(  /* mesh triangle  */
             melem,mnopo ,nelest,nnopst,
             ipemb,iprop ,
            n0elem,n0nopo,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,d1sdel,i1elpr,i1mnnf,i1mnnn,
            i1nobf,i1nopr,i2elto
            ) 
  INT    melem; INT   mnopo; INT  nelest; INT  nnopst;
  INT    ipemb; INT   iprop;
  INT *n0elem; INT  *n0nopo;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nvcx;
  DBL *d1nvcy; DBL  *d1sdel; INT  *i1elpr; INT *i1mnnf; INT *i1mnnn;
  INT *i1nobf; INT *i1nopr; INT **i2elto;
{ INT nelem, nnopo;
  DBL x,y;
  INT i,j,in,jn,innew,ielem,jelem;
  INT i1t[9]={3,0,2, 4,1,0, 5,2,1};
  INT i1new[6];

  nelem=(*n0elem);
  nnopo=(*n0nopo);
  for(ielem=0;ielem<nelest;ielem++)
  { if(i1elpr[ielem]==iprop)
    { for(i=0;i<3;i++)
      { i1new[i+3]=i2elto[i][ielem];
        /* set mid-edge node */
        j=i+1; if(j>2)j=0;
        in=MAXIM(i2elto[i][ielem],i2elto[j][ielem]);
        jn=MINIM(i2elto[i][ielem],i2elto[j][ielem]);
        x=(d1nccx[in]+d1nccx[jn])/R2;
        y=(d1nccy[in]+d1nccy[jn])/R2;
        /* check if mid-edge node already existent */
        innew=i1mnnf[in];
        while(innew>=0)
        { if(ABS(x-d1nccx[innew])+ABS(y-d1nccy[innew])<EPSILON)
          { i1new[i]=innew;
            i1nobf[innew]=0;  /* not a boundary node */
            break;
          }
          innew=i1mnnn[innew-nnopst];
        }
        /* if non-existent, create new mid-edge node */
        if(innew<0)
        { if(nnopo>mnopo)
          { CHRw(stderr,"Ymd: MNOPO too small"); CHRwcr(stderr);
            exit(1);
          }
          i1new[i]=nnopo;
          d1nccx[nnopo]=x;
          d1nccy[nnopo]=y;
          d1ncix[nnopo]=(d1ncix[in]+d1ncix[jn])/R2;
          d1nciy[nnopo]=(d1nciy[in]+d1nciy[jn])/R2;
          d1nvcx[nnopo]=(d1nvcx[in]+d1nvcx[jn])/R2;
          d1nvcy[nnopo]=(d1nvcy[in]+d1nvcy[jn])/R2;
          i1nopr[nnopo]=MINIM(i1nopr[in],i1nopr[jn]);
          i1nobf[nnopo]=1;      /* assume that it is a boundary */
          i1mnnn[nnopo-nnopst]=i1mnnf[in];
          i1mnnf[in]=nnopo;
          nnopo=nnopo+1;
      } }
      /* create new elements */
      i=0;
      for(jelem=0;jelem<3;jelem++)
      { for(in=0;in<3;in++)
        { i2elto[in][nelem]=i1new[i1t[i]];
          i=i+1;
        }
        i1elpr[nelem]=i1elpr[ielem];
        if(d1sdel!=DBL1NULL)d1sdel[nelem]=d1sdel[ielem];
        nelem=nelem+1;     
        if(nelem>melem)
        { CHRw(stderr,"Ymd: MELEM too small"); CHRwcr(stderr);
          exit(1);
      } }
      for(in=0;in<3;in++)
      { i2elto[in][ielem]=i1new[in];
      }
      if(ipemb<1)
      { for(i=0;i<6;i++)
        { i1nobf[i1new[i]]=1;
      } } 
  } } 
  (*n0nopo)=nnopo;
  (*n0elem)=nelem;
}

static void Ymdicsreted2TRIANGLE(  /* mesh triangle discrete */
             melem,mnopo ,nelest,nnopst,
            iprop ,
            n0elem,n0nopo,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,i1elpr,i1nobf,i1nopr,i2elto
            ) 
  INT    melem; INT   mnopo; INT  nelest; INT  nnopst;
  INT    iprop;
  INT *n0elem; INT  *n0nopo;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL *d1nciy; DBL  *d1nvcx;
  DBL *d1nvcy; INT  *i1elpr; INT *i1nobf; INT *i1nopr; INT **i2elto;
{ INT nelem, nnopo;
  INT i,ielem,inopo;
  nelem=(*n0elem);
  nnopo=(*n0nopo);
  for(ielem=0;ielem<nelest;ielem++)
  { if(i1elpr[ielem]==iprop)
    { for(i=0;i<3;i++)
      { inopo=i2elto[i][ielem];
        d1nccx[nnopo]=d1nccx[inopo];
        d1nccy[nnopo]=d1nccy[inopo];
        d1ncix[nnopo]=d1ncix[inopo];
        d1nciy[nnopo]=d1nciy[inopo];
        d1nvcx[nnopo]=d1nvcx[inopo];
        d1nvcy[nnopo]=d1nvcy[inopo];
        i1nopr[nnopo]=i1nopr[inopo];
        i1nobf[nnopo]=1;
        if(nnopo>mnopo)
        { CHRw(stderr,"Ymd: MNOPO too small"); CHRwcr(stderr);
          exit(1);
        }
        i2elto[i][ielem]=nnopo;
        nnopo=nnopo+1;
  } } }
  (*n0nopo)=nnopo;
  (*n0elem)=nelem;
}

static void Ymdskewed2TRIANGLE(  /* mesh triangle  skewed */
             melem,mnopo ,nelest,nnopst,
             ipemb,iprop ,
            n0elem,n0nopo,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,d1sdel,i1elpr,i1mnnf,i1mnnn,
            i1nobf,i1nopr,i2elto
            ) 
  INT    melem; INT   mnopo; INT  nelest; INT  nnopst;
  INT    ipemb; INT   iprop;
  INT *n0elem; INT  *n0nopo;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nvcx;
  DBL *d1nvcy; DBL  *d1sdel; INT  *i1elpr; INT *i1mnnf; INT *i1mnnn;
  INT *i1nobf; INT *i1nopr; INT **i2elto;
{ INT nelem, nnopo;
  DBL x,y;
  INT i,j,k,in,jn,innew,ielem,jelem,kelem;
  INT i1t[18]={3,0,6,  0,4,6, 4,1,6,  1,5,6,  5,2,6,  2,3,6}; 
  INT i1new[7];

  nelem=(*n0elem);
  nnopo=(*n0nopo);
  for(ielem=0;ielem<nelest;ielem++)
  { if(i1elpr[ielem]==iprop)
    { for(i=0;i<3;i++)
      { i1new[i+3]=i2elto[i][ielem];
        /* set mid-edge node */
        j=i+1; if(j>2)j=0;
        in=MAXIM(i2elto[i][ielem],i2elto[j][ielem]);
        jn=MINIM(i2elto[i][ielem],i2elto[j][ielem]);
        x=(d1nccx[in]+d1nccx[jn])/R2;
        y=(d1nccy[in]+d1nccy[jn])/R2;
        /* check if mid-edge node already existent */
        innew=i1mnnf[in];
        while(innew>=0)
        { if(ABS(x-d1nccx[innew])+ABS(y-d1nccy[innew])<EPSILON)
          { i1new[i]=innew;
            i1nobf[innew]=0;  /* not a boundary node */
            break;
          }
          innew=i1mnnn[innew-nnopst];
        }
        /* if non-existent, create new mid-edge node */
        if(innew<0)
        { if(nnopo>mnopo)
          { CHRw(stderr,"Ymd: MNOPO too small"); CHRwcr(stderr);
            exit(1);
          }
          i1new[i]=nnopo;
          d1nccx[nnopo]=x;
          d1nccy[nnopo]=y;
          d1ncix[nnopo]=(d1ncix[in]+d1ncix[jn])/R2;
          d1nciy[nnopo]=(d1nciy[in]+d1nciy[jn])/R2;
          d1nvcx[nnopo]=(d1nvcx[in]+d1nvcx[jn])/R2;
          d1nvcy[nnopo]=(d1nvcy[in]+d1nvcy[jn])/R2;
          i1nopr[nnopo]=MINIM(i1nopr[in],i1nopr[jn]);
          i1nobf[nnopo]=1;      /* assume that it is a boundary */
          i1mnnn[nnopo-nnopst]=i1mnnf[in];
          i1mnnf[in]=nnopo;
          nnopo=nnopo+1;
      } }
      /* create new central node */
      if(nnopo>mnopo)
      { CHRw(stderr,"Ymd: MNOPO too small"); CHRwcr(stderr);
        exit(1);
      }
      i=i1new[3];
      j=i1new[4];
      k=i1new[5];
      i1new[6]=nnopo;
      d1nccx[nnopo]=(d1nccx[i]+d1nccx[j]+d1nccx[k])/R3;
      d1nccy[nnopo]=(d1nccy[i]+d1nccy[j]+d1nccy[k])/R3;
      d1ncix[nnopo]=(d1ncix[i]+d1ncix[j]+d1ncix[k])/R3;
      d1nciy[nnopo]=(d1nciy[i]+d1nciy[j]+d1nciy[k])/R3;
      d1nvcx[nnopo]=(d1nvcx[i]+d1nvcx[j]+d1nvcx[k])/R3;
      d1nvcy[nnopo]=(d1nvcy[i]+d1nvcy[j]+d1nvcy[k])/R3; 
      i1nopr[nnopo]=MINIM(i1nopr[k],MINIM(i1nopr[i],i1nopr[j]));
      i1nobf[nnopo]=0;      /* assume that it is not a boundary */
      nnopo=nnopo+1;
      /* create new elements */
      i=0;
      for(jelem=0;jelem<6;jelem++)
      { if(jelem==0)
        { kelem=ielem;
        }
        else
        { kelem=nelem;
          if(nelem>melem)
          { CHRw(stderr,"Ymd: MELEM too small"); CHRwcr(stderr);
            exit(1);
          }
          nelem=nelem+1;
        }
        for(in=0;in<3;in++)
        { i2elto[in][kelem]=i1new[i1t[i]];
          i=i+1;
        }
        i1elpr[kelem]=i1elpr[ielem];
        if(d1sdel!=DBL1NULL)d1sdel[kelem]=d1sdel[ielem];    
      }
      if(ipemb<1)
      { for(i=0;i<6;i++)
        { i1nobf[i1new[i]]=1;
      } } 
  } } 
  (*n0nopo)=nnopo;
  (*n0elem)=nelem;
}

/*********************PUBLIC********************************************************/
void Ymd(   ydc, yde, ydi, ydn, ydp    /***  mesh elements  ***/
        )
  YDC ydc; YDE yde; YDI ydi;  YDN ydn; YDP ydp;
{ INT nelest;		// actual number of elements
  INT nnopst;		// actual number of nodal points
  INT iprop;		// loop indicator for properties
  INT inopo;		// loop indicator for nodal points
  INT imesh;		// loop indicator for mesh
  INT imestyp;		// 
  INT *i1mnnf;		// mesh node new first for each old node
  INT *i1mnnn;		// mesh node new next for each new node
  INT *i1jnef;		// joint node element first for each old node
  INT *i1jnen;		// joint node element next for each new joint
  DBL *d1sdel;		// 

  for(imesh=0;imesh<10;imesh++)
  { nelest=yde->nelem;
    nnopst=ydn->nnopo;
    i1mnnf=INT1NULL;
    i1mnnn=INT1NULL;
    i1jnef=INT1NULL;
    i1jnen=INT1NULL;
    for(iprop=0;iprop<ydp->nprop;iprop++)
    { imestyp=ydp->i1pemn[iprop]%10;
      ydp->i1pemn[iprop]=ydp->i1pemn[iprop]/10;
      if((ydp->i1ptyp[iprop]==YTE2TRIELS)||
         (ydp->i1ptyp[iprop]==YTE2TRIRIG))
      { d1sdel=yde->d2elst[ydp->i1psde[iprop]];
        if(imestyp==1)
        { if(i1mnnf==INT1NULL)
          { ydi->diedi=ydi->diezon+ydi->diezon;
            i1mnnf=TalINT1(nnopst);
            i1mnnn=TalINT1(3*nelest);
            for(inopo=0;inopo<nnopst;inopo++)
            { i1mnnf[inopo]=-1;
            }
          }
          Ymd2TRIANGLE(  /* mesh triangle  */
          yde->melem   ,ydn->mnopo   ,nelest       ,nnopst     ,
          ydp->i1pemb[iprop],iprop ,
          &(yde->nelem),&(ydn->nnopo),
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],     d1sdel,yde->i1elpr,i1mnnf     ,i1mnnn     ,
          ydn->i1nobf,ydn->i1nopr,yde->i2elto
          );
        }
        else if(imestyp==2)
        { if(i1mnnf==INT1NULL)
          { ydi->diedi=ydi->diezon+ydi->diezon;
            i1mnnf=TalINT1(nnopst);
            i1mnnn=TalINT1(3*nelest);
            for(inopo=0;inopo<nnopst;inopo++)
            { i1mnnf[inopo]=-1;
            }
          }
          Ymdskewed2TRIANGLE(  /* skew triangle  */
          yde->melem   ,ydn->mnopo   ,nelest       ,nnopst     ,
          ydp->i1pemb[iprop],iprop ,
          &(yde->nelem),&(ydn->nnopo),
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],     d1sdel,yde->i1elpr,i1mnnf     ,i1mnnn     ,
          ydn->i1nobf,ydn->i1nopr,yde->i2elto
          );
        }
        else if(imestyp==3)
        { if(i1jnef==INT1NULL)
          { ydi->diedi=ydi->diezon+ydi->diezon;
            i1jnef=TalINT1(nnopst);
            i1jnen=TalINT1(3*nelest);
            for(inopo=0;inopo<nnopst;inopo++)
            { i1jnef[inopo]=-1;
            }
          }
          Yjd2TRIANGLE(  /* create joints  */
          yde->melem   ,ydn->mnopo   ,nelest       ,nnopst     ,
          iprop ,ydp->i1pejp[iprop],
          &(yde->nelem),&(ydn->nnopo),
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],d1sdel,
          yde->i1elpr,i1jnef,i1jnen,ydn->i1nobf,ydn->i1nopr,
          yde->i2elto,yde->i2eljp,yde->i1elty
          );
        }
        else if(imestyp==4)
        { ydi->diedi=ydi->diezon+ydi->diezon;
          Ymdicsreted2TRIANGLE(  /* separate triangle  */
          yde->melem   ,ydn->mnopo   ,nelest       ,nnopst     ,
          iprop ,
          &(yde->nelem),&(ydn->nnopo),
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],yde->i1elpr,ydn->i1nobf,ydn->i1nopr,yde->i2elto
          );
	    }
	  }
    }
    /* free memory */
    FREE(i1mnnn);
    FREE(i1mnnf);
    FREE(i1jnef);
    FREE(i1jnen);
  }
}
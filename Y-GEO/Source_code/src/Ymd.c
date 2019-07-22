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
static void Yjd2TRIANGLE(  /* mesh triangle  */
             melem,mnopo ,nelest,nnopst,
             iprop, i2pmij, mpmrow, nprop,
            n0elem,n0nopo,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,d1sdel,
            i1elpr, i1elprtmp, i1jnef,i1jnen,i1next, i1nobf,i1nopr,i1nowe,
            i2elto,
	    d1elfr, iusefn, mdfnfr, mdfnno, i2dfnn, d1dffr, d1dfpe, d1elpe, d1dfpt, d1elpt,
            i2elnext, i1edfnf
            )
  INT    melem; INT   mnopo; INT  nelest; INT  nnopst;
  INT    iprop; INT **i2pmij; INT mpmrow; INT nprop;
  INT *n0elem; INT  *n0nopo;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nvcx;
  DBL *d1nvcy; DBL  *d1sdel;
  INT *i1elpr; INT  *i1elprtmp; INT *i1jnef; INT *i1jnen; INT *i1next; INT *i1nobf; INT *i1nopr; INT *i1nowe;
  INT **i2elto;
  DBL *d1elfr; INT iusefn; INT mdfnfr; INT mdfnno; INT **i2dfnn; DBL *d1dffr;
  DBL *d1elpe; DBL *d1dfpe; DBL *d1elpt; DBL *d1dfpt; 
  INT **i2elnext; INT *i1edfnf; 
  
{ INT nelem, nnopo;
  INT i,j,in,jn,kn,ijnew,ielem,irow;
  INT ipropj = -1;
 
  INT fracture_flag; /* Flag that identifies a joint (i.e. element edge) element that belongs to the DFN */
  DBL xi_0, yi_0, xi_1, yi_1;   
  INT k,r,s;
  
  nelem=(*n0elem);
  nnopo=(*n0nopo);
  for(ielem=0;ielem<nelest;ielem++)
  { if(i1elprtmp[ielem]==iprop)
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
        i1nowe[nnopo]=i1nowe[in];
	i1nobf[nnopo]=1;
        nnopo=nnopo+1;
	
	fracture_flag = 0;
	// Check if the two nodes belong to a crack
	if(iusefn > 0)
	{ xi_0=d1ncix[in];
	  yi_0=d1nciy[in];
	  xi_1=d1ncix[jn];
	  yi_1=d1nciy[jn];
	  for(s=0; s<mdfnfr; s++)
	  { for(k=0; k<mdfnno; k++)
	    { if(i2dfnn[k][s] >= 0)  
	      { if(xi_0 == d1ncix[i2dfnn[k][s]])
	        { if(yi_0 == d1nciy[i2dfnn[k][s]])
	          { for(r=0; r<mdfnno; r++)
	            { if(i2dfnn[r][s] >= 0)
		      { if(xi_1 == d1ncix[i2dfnn[r][s]])
	                { if(yi_1 == d1nciy[i2dfnn[r][s]]) 
	                  { fracture_flag = 1;
			    d1elfr[ielem] = d1dffr[s]; //! Assign DFN friction to triangle
			    d1elpe[ielem] = d1dfpe[s]; //! Assign DFN normal penalty to triangle
			    d1elpt[ielem] = d1dfpt[s]; //! Assign DFN tangential penalty to triangle
		          }
         } } } } } } } } }
	
        /* check if joint element already existent */
        ijnew=i1jnef[kn];
        while(ijnew>=0)
        { if((i2elto[2][ijnew]==in)&&(i2elto[3][ijnew]==jn))
          { i2elnext[1][ijnew]=ielem; /* element# next to the joint    */
            i2elnext[3][ijnew]=i;     /* edge# (0,1,2) of that element */
            if(i1elpr[ielem] == i1elpr[i1next[ijnew]])
            { for(irow=0; irow<mpmrow;irow++)
              { if( (i2pmij[0][irow] == i1elpr[ielem]) &&
                    (i2pmij[1][irow] == i1elpr[ielem]) )
                { i1elpr[ijnew]=i2pmij[2][irow]+nprop;
		
		  if(fracture_flag==1) 
		  { if(iusefn==1) //! Break elements only when iusefn==1
		    { i1elpr[ijnew]=ipropj-YIPROPMAX;
		    }
                    //! For hydrofrac and cohesive DFN
		    i1edfnf[ijnew] = 1;
                  }  
	    } } } 
            else if(i1elpr[ielem] != i1elpr[i1next[ijnew]])
            { for(irow=0; irow<mpmrow;irow++)
              { if( (i2pmij[0][irow] == i1elpr[ielem]) ||
                    (i2pmij[1][irow] == i1elpr[ielem]) )
                { if( (i2pmij[0][irow] == i1elpr[i1next[ijnew]]) ||
                      (i2pmij[1][irow] == i1elpr[i1next[ijnew]]) )
                  { i1elpr[ijnew]=i2pmij[2][irow]+nprop;
		  
		     if(fracture_flag==1)
		     { if(iusefn==1) //! Break elements only when iusefn==1
		       { i1elpr[ijnew]=ipropj-YIPROPMAX; }
                       //! For hydrofrac and cohesive DFN
                       i1edfnf[ijnew] = 1;
                     }
            } } } }
            if(i<2)
            { i2elto[2][ijnew]=nnopo-1;
              i2elto[3][ijnew]=nnopo;
            }
            else
            { i2elto[2][ijnew]=nnopo-1;
              i2elto[3][ijnew]=nnopo-3;
            }
            ijnew=-100;
          } 
          else
          { ijnew=i1jnen[ijnew-nelest];
        } }
        /* create new mid-edge joint node  */
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

          for(irow=0; irow<mpmrow;irow++)
          {
            if(i2pmij[0][irow] == i2pmij[1][irow])
              if(i2pmij[0][irow] == i1elpr[ielem])
                ipropj=i2pmij[2][irow]+nprop;
          }
          i1elpr[nelem]=ipropj;
          i1next[nelem]=ielem;
          i2elnext[0][nelem]=ielem; /* element# next to the joint    */
          i2elnext[2][nelem]=i;     /* edge# (0,1,2) of that element */
          i1jnen[nelem-nelest]=i1jnef[kn];
          i1jnef[kn]=nelem;
	  
	  if(fracture_flag==1) 
          { if(iusefn==1) //! Break elements only when iusefn==1
	    { i1elpr[nelem]=ipropj-YIPROPMAX; }
            //! For hydrofrac and cohesive DFN
	    i1edfnf[nelem] = 1;
	  }
  
	  nelem=nelem+1;
	} }
      for(i=0;i<3;i++) /* detach element */
      { i2elto[i][ielem]=nnopo-3+i;
      }
  } }
  for(ielem=nelest;ielem<nelem;ielem++)
  { if(d1sdel!=DBL1NULL)d1sdel[ielem]=R0;
    for(i=2;i<4;i++)
    if(i2elto[2][ielem]<nnopst)
    { i2elto[3][ielem]=i2elto[0][ielem];
      i2elto[2][ielem]=i2elto[1][ielem];
  } }
  (*n0nopo)=nnopo;
  (*n0elem)=nelem;
}


/**************MESH ELEMENTS***********/
static void Ymd2TRIANGLE(  /* mesh triangle  */
             melem,mnopo ,nelest,nnopst,
             ipemb,iprop ,
            n0elem,n0nopo,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,d1sdel,i1elpr, i1elprtmp, i1mnnf,i1mnnn,
            i1nobf,i1nopr,i1nowe,i2elto
            )
  INT    melem; INT   mnopo; INT  nelest; INT  nnopst;
  INT    ipemb; INT   iprop;
  INT *n0elem; INT  *n0nopo;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nvcx;
  DBL *d1nvcy; DBL  *d1sdel; INT  *i1elpr; INT *i1elprtmp; INT *i1mnnf; INT *i1mnnn;
  INT *i1nobf; INT *i1nopr; INT *i1nowe; INT **i2elto;
{ INT nelem, nnopo;
  DBL x,y;
  INT i,j,in,jn,innew,ielem,jelem;
  INT i1t[9]={3,0,2, 4,1,0, 5,2,1};
  INT i1new[6];
 
  nelem=(*n0elem);
  nnopo=(*n0nopo);
  for(ielem=0;ielem<nelest;ielem++)
  { if(i1elprtmp[ielem]==iprop)
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
          i1nowe[nnopo]=MINIM(i1nowe[in],i1nowe[jn]);
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
            d1nvcy,i1elpr,i1nobf,i1nopr,i1nowe,i2elto
            )
  INT    melem; INT   mnopo; INT  nelest; INT  nnopst;
  INT    iprop;
  INT *n0elem; INT  *n0nopo;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL *d1nciy; DBL  *d1nvcx;
  DBL *d1nvcy; INT  *i1elpr; INT *i1nobf; INT *i1nopr; INT *i1nowe; INT **i2elto;
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
        i1nowe[nnopo]=i1nowe[inopo];
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
            d1nvcy,d1sdel,i1elpr, i1elprtmp, i1mnnf,i1mnnn,
            i1nobf,i1nopr,i1nowe,i2elto
            )
  INT    melem; INT   mnopo; INT  nelest; INT  nnopst;
  INT    ipemb; INT   iprop;
  INT *n0elem; INT  *n0nopo;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nvcx;
  DBL *d1nvcy; DBL  *d1sdel; INT  *i1elpr; INT *i1elprtmp; INT *i1mnnf; INT *i1mnnn;
  INT *i1nobf; INT *i1nopr; INT *i1nowe; INT **i2elto;
{ INT nelem, nnopo;
  DBL x,y;
  INT i,j,k,in,jn,innew,ielem,jelem,kelem;
  INT i1t[18]={3,0,6,  0,4,6, 4,1,6,  1,5,6,  5,2,6,  2,3,6};
  INT i1new[7];

  nelem=(*n0elem);
  nnopo=(*n0nopo);
  for(ielem=0;ielem<nelest;ielem++)
  { if(i1elprtmp[ielem]==iprop)
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
          i1nowe[nnopo]=MINIM(i1nowe[in],i1nowe[jn]);
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
      i1nowe[nnopo]=MINIM(i1nowe[k],MINIM(i1nowe[i],i1nowe[j]));
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
void Ymd(   ydc, yde, ydi, ydn, ydpe, ydpn, ydpm, ydfn    /***  mesh elements  ***/
        )
  YDC ydc; YDE yde; YDI ydi;  YDN ydn; YDPE ydpe; YDPN ydpn; YDPM ydpm; YDFN ydfn;
{ INT nelest, nnopst;
  INT iprop;
  INT icom, i,j, ielem, irow;
  INT inopo;
  INT imesh,imestyp;
  INT nProps;         /* max of number Pr sets to be meshed together in one line of i2pmset  */
  INT *i1props;       /* List of Property sets to be meshed together                         */
  INT *i1mnnf;        /* mesh node new first for each old node                               */
  INT *i1mnnn;        /* mesh node new next for each new node                                */
  INT *i1jnef;        /* joint node element first for each old node                          */
  INT *i1jnen;        /* joint node element next for each new joint                          */
  DBL *d1sdel;
  INT *i1next;        /* stores index of the element next to joint                           */
  INT *i1isMeshing;   /* indicates if meshing is allowed for the selected Pr Set,
                         only if prtype = YTE2TRIELS || YTE2TRIRIG || YTE2PLANESTRESS ||  
                         YTE2PLANESTRAIN */
  INT i0, i1, i2, i3;
  
  /* initializing the array of element indexes next to joints and flags of broken element edges */
  if(yde->i2elnext==INT2NULL)
  { yde->i2elnext=TalINT2(4,yde->melem);
    for (i=0;i<4;i++)
    { for(j=0; j<yde->melem; j++)
      { yde->i2elnext[i][j]=-1;
    } }
    yde->i2eledge=TalINT2(yde->nelno,yde->melem);
    for (i=0;i<yde->nelno;i++)
    { for(j=0; j<yde->melem; j++)
      { yde->i2eledge[i][j]=-1;   /* -1: broken (+1: intact) */
    } }
  }
  
  /* initializing the array of node indexes next to each node in clockwise and counter-clockwise order*/
  if(ydn->i2noid==INT2NULL)
  { ydn->i2noid=TalINT2(2,ydn->mnopo);
    for (i=0;i<2;i++)
    { for(j=0; j<ydn->mnopo; j++)
      { ydn->i2noid[i][j]=-1;
    } }
  }
  
  //! Allocating memory and initializing i1edfnf with 0 only if a DFN is defined
  if(ydfn->iusefn > 0)
  { if(yde->i1edfnf==INT1NULL)
    { yde->i1edfnf=TalINT1(yde->melem);
      for(ielem=0;ielem<yde->melem;ielem++)
      { yde->i1edfnf[ielem]=0;
      } 
  } }
    
  if(ydpe->i1pejp==INT1NULL)
  { ydpe->i1pejp=TalINT1(ydpe->mprop);
  }
  i1isMeshing=INT1NULL;
  i1isMeshing=TalINT1(ydpm->mpmcom);
  for(icom=0;icom<ydpm->mpmcom;icom++)
  { i1isMeshing[icom]=-1;  
  }
  for(icom=0;icom<ydpm->mpmcom;icom++)
  { nProps=ydpm->i2pmset[0][icom];
    i1props=INT1NULL;
    if(i1props==INT1NULL)
    { i1props=TalINT1(nProps);
    }
    for(i=0; i<nProps; i++)
    { i1props[i]=ydpm->i2pmset[2+i][icom];
      if( (ydpe->i1ptyp[i1props[i]]==YTE2TRIELS) ||
          (ydpe->i1ptyp[i1props[i]])==(YTE2PLANESTRESS) ||
          (ydpe->i1ptyp[i1props[i]])==(YTE2PLANESTRAIN) ||
          (ydpe->i1ptyp[i1props[i]]==YTE2TRIRIG) )
      { i1isMeshing[icom]=1;
      }
      else
      { i1isMeshing[icom]=-1;
      }
    }
    FREE(i1props);
  }
  //! Allocating memory and initializing d1elfr with -1.0 only if a DFN is defined
  if(ydfn->iusefn > 0)
  { if(yde->d1elfr==DBL1NULL)
    { yde->d1elfr=TalDBL1(yde->melem);
      for(ielem=0;ielem<yde->melem;ielem++)
      { yde->d1elfr[ielem]=-1.0; } 
  } }
  //! Allocating memory and initializing d1elpe with -1.0 only if a DFN is defined
  if(ydfn->iusefn > 0)
  { if(yde->d1elpe==DBL1NULL)
    { yde->d1elpe=TalDBL1(yde->melem);
      for(ielem=0;ielem<yde->melem;ielem++)
      { yde->d1elpe[ielem]=-1.0;
      } 
  } }
  //! Allocating memory and initializing d1elpt with -1.0 only if a DFN is defined
  if(ydfn->iusefn > 0)
  { if(yde->d1elpt==DBL1NULL)
    { yde->d1elpt=TalDBL1(yde->melem);
      for(ielem=0;ielem<yde->melem;ielem++)
      { yde->d1elpt[ielem]=-1.0;
      } 
  } }
  //! Allocating memory and initializing i1edft with -1 only if a DFN type = 3 is defined
  if(ydfn->iusefn == 3)
  { if(yde->i1edft==INT1NULL)
    { yde->i1edft=TalINT1(yde->melem);
      for(ielem=0;ielem<yde->melem;ielem++)
      { yde->i1edft[ielem]=-1;
      } 
  } }
    
  if(ydc->ncstep==0)
  { for(ielem=0;ielem<yde->nelem;ielem++)
    { for(i=0;i<3;i++)
      { yde->i2eltost[i][ielem]=yde->i2elto[i][ielem];
      }
    }
  }

  for(imesh=0;imesh<10;imesh++)
  { nelest=yde->nelem;
    nnopst=ydn->nnopo;
    i1mnnf=INT1NULL;
    i1mnnn=INT1NULL;
    i1jnef=INT1NULL;
    i1jnen=INT1NULL;
    i1next=INT1NULL;
    for(icom=0;icom<ydpm->mpmcom;icom++)
    { iprop = ydpm->i2pmset[2][icom];
      imestyp=ydpm->i2pmset[1][icom]%10;
      ydpm->i2pmset[1][icom]=ydpm->i2pmset[1][icom]/10;
      if(i1isMeshing[icom]>0)
      { for(irow=0; irow<ydpm->mpmrow; irow++)
        { if(ydpm->i2pmij[0][irow] == ydpm->i2pmij[1][irow])
          { ydpe->i1pejp[ydpm->i2pmij[0][irow]] = ydpm->i2pmij[2][irow];
        } }
        d1sdel=yde->d2elst[ydpe->i1psde[iprop]];
        if(imestyp>0)    /* construct i1elprtmp only when refining! */
        { for(i=0; i<ydpm->i2pmset[0][icom]; i++)
          { for(ielem=0; ielem<yde->nelem; ielem++)
            { if(yde->i1elpr[ielem] == ydpm->i2pmset[2+i][icom])
              {  yde->i1elprtmp[ielem]=icom+prArbit;
        } } } }
        if(imestyp==1)
        { if(i1mnnf==INT1NULL)
          { ydi->diedi=ydi->diezon+ydi->diezon;
            i1mnnf=TalINT1(nnopst);
            i1mnnn=TalINT1(3*nelest);
            for(inopo=0;inopo<nnopst;inopo++)
            { i1mnnf[inopo]=-1;
            }
          }
          iprop = icom+prArbit;
          Ymd2TRIANGLE(  /* mesh triangle  */
          yde->melem   ,ydn->mnopo   ,nelest       ,nnopst     ,
          ydpe->i1pemb[ydpm->i2pmset[2][icom]],iprop ,
          &(yde->nelem),&(ydn->nnopo),
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],     d1sdel,yde->i1elpr, yde->i1elprtmp, i1mnnf     ,i1mnnn     ,
          ydn->i1nobf,ydn->i1nopr,ydn->i1nowe,yde->i2elto
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
          iprop = icom+prArbit;
          Ymdskewed2TRIANGLE(  /* skew triangle  */
          yde->melem   ,ydn->mnopo   ,nelest       ,nnopst     ,
          ydpe->i1pemb[ydpm->i2pmset[2][icom]],iprop ,
          &(yde->nelem),&(ydn->nnopo),
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],     d1sdel,yde->i1elpr, yde->i1elprtmp, i1mnnf     ,i1mnnn     ,
          ydn->i1nobf,ydn->i1nopr,ydn->i1nowe,yde->i2elto
          );
        }
        else if(imestyp==3)
        { if(i1jnef==INT1NULL)
          { ydi->diedi=ydi->diezon+ydi->diezon;
            i1jnef=TalINT1(nnopst);
            i1jnen=TalINT1(3*nelest);
            i1next=TalINT1(yde->melem);
            for(inopo=0;inopo<nnopst;inopo++)
            { i1jnef[inopo]=-1;
            }
            for(ielem=0;ielem<nelest;ielem++)
            { i1next[ielem]=-1;
            }
          }
          iprop = icom+prArbit;
          Yjd2TRIANGLE(  /* create joints  */
          yde->melem   ,ydn->mnopo   ,nelest       ,nnopst     ,
          iprop, ydpm->i2pmij, ydpm->mpmrow, ydpe->nprop,
          &(yde->nelem),&(ydn->nnopo),
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],d1sdel,
          yde->i1elpr, yde->i1elprtmp, i1jnef,i1jnen, i1next, ydn->i1nobf,ydn->i1nopr,ydn->i1nowe,
          yde->i2elto,
          yde->d1elfr, ydfn->iusefn, ydfn->mdfnfr, ydfn->mdfnno, ydfn->i2dfnn, ydfn->d1dffr, ydfn->d1dfpe, yde->d1elpe, ydfn->d1dfpt, yde->d1elpt,
          yde->i2elnext,yde->i1edfnf);
        }
        else if(imestyp==4)
        { ydi->diedi=ydi->diezon+ydi->diezon;
          Ymdicsreted2TRIANGLE(  /* separate triangle  */
          yde->melem   ,ydn->mnopo   ,nelest       ,nnopst     ,
          iprop ,
          &(yde->nelem),&(ydn->nnopo),
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],yde->i1elpr,ydn->i1nobf,ydn->i1nopr,ydn->i1nowe,yde->i2elto
          );
        }
      }
    }
    /* free memory */
    FREE(i1mnnn);
    FREE(i1mnnf);
    FREE(i1jnef);
    FREE(i1jnen);
    FREE(i1next);
  }
  FREE(i1isMeshing);
}

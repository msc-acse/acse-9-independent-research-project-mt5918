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
/* File   Yid.c */
#include "Yproto.h"
/**********************GENERALISED INTERACTION FORCES*****************************************/
static void Yid2TRITRI(              /* Triangle to Triangle */
             nelem,nicoup,dctime, /* nicoup added by RG the 29/03/06 */ 
            dcstec,icfmty,diezon,
			iprop,  jprop,
            d1iesl,d1nccx,d1nccy,d1nfcx,d1nfcy,
            d1nvcx,d1nvcy,d1pefr,d1pepc,d1pbif,
			i1elcf,i1elpr,i2eljp,i2elto,i1elty,
			d1elsf,d2elcf,i1icff,i1iecn,i1iect,
			d2sldis,njoint,d1jkni,d1jknc,d1jksc,
			d1jphi,d1jnst,d1jfpr
            )
  INT    nelem;
  INT  *nicoup; DBL   dctime; /* nicoup added by RG the 29/03/06 */
  DBL   dcstec; INT   icfmty; DBL   diezon;
  INT    iprop; INT    jprop;
  DBL  *d1iesl; DBL  *d1nccx; DBL  *d1nccy; DBL  *d1nfcx; DBL *d1nfcy;
  DBL  *d1nvcx; DBL  *d1nvcy; DBL  *d1pefr; DBL  *d1pepc; DBL *d1pbif;
  INT  *i1elcf; INT  *i1elpr; INT **i2eljp; INT **i2elto; INT *i1elty;
  DBL  *d1elsf; DBL **d2elcf; INT  *i1icff; INT  *i1iecn; INT *i1iect;
  DBL **d2sldis; INT  njoint; DBL  *d1jkni; DBL  *d1jknc; DBL  *d1jksc;
  DBL  *d1jphi; DBL  *d1jnst; DBL  *d1jfpr;
{ INT kprop,icontact,ielem,jelem,icoup,jcoup,it,jt,in,jn,ie,je,ip,jp,np;
  DBL a0,a1,a2,b0,b1,b2,c0,c1,c2,n0,n1,n2;
  DBL pen,penT,tmp,dmin2,smin,smax,PHIR,PHII;
  DBL resfric,totfric,saver,aver0,aver1,aver2;
  DBL vrelx, vrely;
  DBL vx[2][3];
  DBL vy[2][3];
  DBL vrel,vreln;
  DBL fn,fna,fnb;   /* normal forces      */
  DBL ft,fta,ftb;	/* tangential forces  */
  DBL small=EPSILON;
  DBL nsmall=-EPSILON;
  DBL big=BEPSILON;
  DBL zone2;
  DBL p[10];
  DBL s[10];
  DBL fx[3];
  DBL fy[3];
  DBL vol[2];
  DBL rx[2][3];
  DBL ry[2][3];
  DBL nx[2][3];
  DBL ny[2][3];
  DBL d[2][3][3];
  INT i2to[2][3];
  DBL *d1sldis;     /* sliding distance */
  INT icoupID;
  INT ijoint,jjoint;
  DBL defsr;		/* effective stress ratio */

  d1sldis = NULL;

  /* interaction calculation */
  zone2=(R4*diezon*diezon);
  pen=MINIM(d1pepc[iprop],d1pepc[jprop]);
  penT=RP1*pen;
  if(iprop!=jprop)	  // residual friction coefficient
  { resfric=MINIM(d1pbif[iprop],d1pbif[jprop]);
  }
  else if(d1pefr==DBL1NULL)
  { resfric=R0;
  }
  else
  { resfric=MINIM(d1pefr[iprop],d1pefr[jprop]);
  }
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { kprop=jprop;
    }
    else if(i1elpr[ielem]==jprop)
    { kprop=iprop;
    }
    else
    { kprop=-1;
    }
    if(kprop>=0)
    { icoup=i1elcf[ielem];
      jcoup=-1;
      while(icoup>=0)
      { if(i1iect[icoup]<0) i1iect[icoup]=-1-i1iect[icoup];
        jelem=i1iect[icoup];
        icontact=-2;
        if(i1elpr[jelem]==kprop)
        { icontact=-1;
          jt=jelem;
          for(it=0;it<2;it++)
          { for(in=0;in<3;in++)
            { i2to[it][in]=i2elto[in][jt];
            }
            jt=ielem;
          }
          for(it=0;it<2;it++)
          { for(in=0;in<3;in++)
            { rx[it][in]=d1nccx[i2to[it][in]];
              ry[it][in]=d1nccy[i2to[it][in]];
              vx[it][in]=d1nvcx[i2to[it][in]];  /* velocities for friction */
              vy[it][in]=d1nvcy[i2to[it][in]];
          } }
          for(it=0;it<2;it++)
          { vol[it]=(rx[it][1]-rx[it][0])*(ry[it][2]-ry[it][0])-
                    (ry[it][1]-ry[it][0])*(rx[it][2]-rx[it][0]);
            for(ie=0;ie<3;ie++)
            { je=ie+1; if(je>2)je=0;
              nx[it][ie]=ry[it][je]-ry[it][ie];
              ny[it][ie]=rx[it][ie]-rx[it][je];
          } }
          for(it=0;it<2;it++)
          { jt=it+1; if(jt>1)jt=0;
            for(in=0;in<3;in++)
            { for(ie=0;ie<3;ie++)
              { d[it][in][ie]=((rx[jt][ie]-rx[it][in])*nx[jt][ie]+
                (ry[jt][ie]-ry[it][in])*ny[jt][ie])/vol[jt];
          } } }
          dmin2=big;
		  icoupID=-1;
          for(it=0;it<2;it++)
          { jt=it+1; if(jt>1)jt=0;
            for(in=0;in<3;in++)
            { fx[in]=R0; fy[in]=R0;
            }
            n0=(nx[jt][0]*nx[jt][0]+ny[jt][0]*ny[jt][0])/
               (vol[jt]*vol[jt]);
            n1=(nx[jt][1]*nx[jt][1]+ny[jt][1]*ny[jt][1])/
               (vol[jt]*vol[jt]);
            n2=(nx[jt][2]*nx[jt][2]+ny[jt][2]*ny[jt][2])/
               (vol[jt]*vol[jt]);
            for(in=0;in<3;in++)
            { icoupID=icoupID+1;
              d1sldis=d2sldis[icoupID];

			  jn=in+1; if(jn>2)jn=0;
              a0=d[it][in][0];
              a1=d[it][in][1];
              a2=d[it][in][2];
              b0=d[it][jn][0];
              b1=d[it][jn][1];
              b2=d[it][jn][2];
              c0=d[jt][0][in];
              c1=d[jt][1][in];
              c2=d[jt][2][in];
              /* check if contact */
              if((((c0>nsmall)&&(c1>nsmall)&&(c2>nsmall))||
                  ((c0<small)&&(c1<small)&&(c2<small)))||
                 (((a0<small)&&(b0<small))||((a1<small)&&(b1<small))||
                  ((a2<small)&&(b2<small))))
              { if((a0<=a1)&&(a0<=a2))
                { dmin2=MINIM(dmin2,(a0*a0/n0));
                }
                else if((a1<=a0)&&(a1<=a2))
                { dmin2=MINIM(dmin2,(a1*a1/n1));
                }
                else
                { dmin2=MINIM(dmin2,(a2*a2/n2));
                }
              }
              else
              { icontact=it;
                /* domain of contact */
                smin=R0; smax=R1;
                if((a0<R0)&&(b0>small))smin=MAXIM(smin,(a0/(a0-b0)));
                if((a1<R0)&&(b1>small))smin=MAXIM(smin,(a1/(a1-b1)));
                if((a2<R0)&&(b2>small))smin=MAXIM(smin,(a2/(a2-b2)));
                if((a0>small)&&(b0<R0))smax=MINIM(smax,(a0/(a0-b0)));
                if((a1>small)&&(b1<R0))smax=MINIM(smax,(a1/(a1-b1)));
                if((a2>small)&&(b2<R0))smax=MINIM(smax,(a2/(a2-b2)));
                if(smax>smin)
                { s[0]=smin;
                  p[0]=MINIM((a0+smin*(b0-a0)),(a1+smin*(b1-a1)));
                  p[0]=MINIM(p[0],(a2+smin*(b2-a2)));
                  np=1;
                  /* intermediate points */
                  tmp=b0-a0+a1-b1;
                  if((DABS(tmp))>small)
                  { tmp=(a1-a0)/tmp;
                    if((tmp>smin)&&(tmp<smax)&&
                       ((a0+tmp*(b0-a0))<(a2+tmp*(b2-a2))))
                    { s[np]=tmp;
                      p[np]=a0+tmp*(b0-a0);
                      np=np+1;
                  } }
                  tmp=b0-a0+a2-b2;
                  if((DABS(tmp))>small)
                  { tmp=(a2-a0)/tmp;
                    if((tmp>smin)&&(tmp<smax)&&
                       ((a0+tmp*(b0-a0))<(a1+tmp*(b1-a1))))
                    { s[np]=tmp;
                      p[np]=a0+tmp*(b0-a0);
                      np=np+1;
                  } }
                  tmp=b1-a1+a2-b2;
                  if((DABS(tmp))>small)
                  { tmp=(a2-a1)/tmp;
                    if((tmp>smin)&&(tmp<smax)&&
                       ((a1+tmp*(b1-a1))<(a0+tmp*(b0-a0))))
                    { s[np]=tmp;
                      p[np]=a1+tmp*(b1-a1);
                      np=np+1;
                  } }
                  s[np]=smax;
                  p[np]=MINIM((a0+smax*(b0-a0)),(a1+smax*(b1-a1)));
                  p[np]=MINIM(p[np],(a2+smax*(b2-a2)));
                  np=np+1;
                  /* order intermediate points */
                  for(ip=0;ip<(np-1);ip++)
                  { for(jp=(ip+1);jp<np;jp++)
                    { if(s[ip]>s[jp])
                      { tmp=s[jp]; s[jp]=s[ip]; s[ip]=tmp;
                        tmp=p[jp]; p[jp]=p[ip]; p[ip]=tmp;
                  } } }

                  /* calculate relative velocity and sliding distance */
                  saver = R0;
                  for(ip=0;ip<np;ip++)
                  { saver=saver+s[ip];
                  }
                  saver=saver/((DBL)np);
                  aver0=((a0*(R1-saver))+(b0*saver));
                  aver1=((a1*(R1-saver))+(b1*saver));
                  aver2=((a2*(R1-saver))+(b2*saver));
                  vrelx=(aver0*vx[jt][2]+aver1*vx[jt][0]+aver2*vx[jt][1])-
                        ((R1-saver)*vx[it][in]+saver*vx[it][jn]);
                  vrely=(aver0*vy[jt][2]+aver1*vy[jt][0]+aver2*vy[jt][1])-
                        ((R1-saver)*vy[it][in]+saver*vy[it][jn]);
                  vrel=-vrelx*ny[it][in]+vrely*nx[it][in];
				  vreln=vrelx*nx[it][in]-vrely*ny[it][in];

                  d1sldis[icoup]=d1sldis[icoup]+vrel*dcstec;      /* sliding distance */

                  /* integrate normal force */
                  fn=p[0]*(s[1]-s[0])+p[np-1]*(s[np-1]-s[np-2]);
                  fnb=p[0]*(s[1]-s[0])*(s[1]+R2*s[0])+
                      p[np-1]*(s[np-1]-s[np-2])*(s[np-2]+R2*s[np-1]);
                  for(ip=1;ip<(np-1);ip++)
                  { fn=fn+p[ip]*(s[ip+1]-s[ip-1]);
                    fnb=fnb+p[ip]*(
                    (s[ip]-s[ip-1])*(s[ip-1]+R2*s[ip])+
                    (s[ip+1]-s[ip])*(s[ip+1]+R2*s[ip]));
                  }
                  fnb=fnb*pen*RP5;
                  fn=fn*pen*RP15;
                  fna=fn-fnb;

				  /* effective stress ratio - for friction calculation */
				  defsr=R1;
				  if(i2eljp[in][jelem]>=0)
				  { if(i1elty[i2eljp[in][jelem]]>1)
					{ // contact joint ids
					  jjoint=i2eljp[in][jelem]-(nelem-njoint);
					  if(d1jnst[jjoint]>R0)
					  { defsr=R1-MINIM(d1jfpr[jjoint]/d1jnst[jjoint],R1);
					  }
					  else
					  { defsr=R1;
					  }
				  } }

				  /* tangential forces */
				  if(icfmty==0)
				  { /* Coulomb friction model */
					ft=d1sldis[icoup]*penT;
					if(DABS(ft)>(resfric*DABS(fn*defsr)))
					{ ft=ft*DABS(resfric*fn*defsr/ft);
					  d1sldis[icoup]=ft/penT;
				  } }
				  else if(icfmty==1)
				  { /* Barton-Bandis friction model */
					if((i1elty[i2eljp[in][jelem]]>0)&&(i1elty[i2eljp[in][jelem]]<3))
					{ // pre-existing fractures -- friction + dilation
					  PHIR=atan(resfric)*RAD2DEG;
					  PHII=d1jphi[jjoint];
					  for(ip=0;ip<3;ip++)
					  { if(i1elty[i2eljp[ip][ielem]]>1)
						{ ijoint=i2eljp[ip][ielem]-(nelem-njoint);
						  PHII=MINIM(PHII,d1jphi[ijoint]);
					  } }
					  if((PHIR+PHII)>85.0)
					  { totfric=tan(85.0*DEG2RAD);
					  }
					  else
					  { totfric=tan((PHIR+PHII)*DEG2RAD);
					} }
					else if(i1elty[i2eljp[in][jelem]]>=3)
					{ // new fractures -- only residual friction
					  totfric=resfric;
					}
					else
					{ totfric=R0;
					}
					if(vrel!=R0)
					{ ft=vrel*DABS(totfric*fn*defsr/vrel);
					}
					else
					{ ft=R0;
					}
				  }
                  fta=ft*fna/fn;
                  ftb=ft*fnb/fn;

                  /* update total force */
                  fx[in]=fx[in]-fna*nx[it][in]-fta*ny[it][in];
                  fy[in]=fy[in]-fna*ny[it][in]+fta*nx[it][in];
                  fx[jn]=fx[jn]-fnb*nx[it][in]-ftb*ny[it][in];
                  fy[jn]=fy[jn]-fnb*ny[it][in]+ftb*nx[it][in];

				  /* store frictional force */
				  d1elsf[ielem]=d1elsf[ielem]+ft*SQRT(nx[it][in]*nx[it][in]+ny[it][in]*ny[it][in]);
				  d1elsf[jelem]=d1elsf[jelem]-ft*SQRT(nx[it][in]*nx[it][in]+ny[it][in]*ny[it][in]);

				  /* store boundary contact force */
				  /*if(i1elpr[ielem]!=i1elpr[jelem])
				  { d2elcf[ielem][0]=d2elcf[ielem][0]-fn*nx[it][in];
					d2elcf[ielem][1]=d2elcf[ielem][1]-fn*ny[it][in];
					d2elcf[jelem][0]=d2elcf[jelem][0]+fn*nx[it][in];
					d2elcf[jelem][1]=d2elcf[jelem][1]+fn*ny[it][in];
				  }*/
				  if(i1elpr[ielem]!=i1elpr[jelem])
				  { d2elcf[ielem][0]=d2elcf[ielem][0]-fn*nx[it][in]-ft*ny[it][in];
					d2elcf[ielem][1]=d2elcf[ielem][1]-fn*ny[it][in]+ft*nx[it][in];
					d2elcf[jelem][0]=d2elcf[jelem][0]+fn*nx[it][in]+ft*ny[it][in];
					d2elcf[jelem][1]=d2elcf[jelem][1]+fn*ny[it][in]-ft*nx[it][in];
				  }
            } } }
            if(icontact==it) /* update nodal forces  */
            { for(in=0;in<3;in++)
              { d1nfcx[i2to[it][in]]=d1nfcx[i2to[it][in]]+fx[in];
                d1nfcy[i2to[it][in]]=d1nfcy[i2to[it][in]]+fy[in];
                ie=in+1; if(ie>2)ie=0;
                for(jn=0;jn<3;jn++)
                { d1nfcx[i2to[jt][in]]=d1nfcx[i2to[jt][in]]-fx[jn]*d[it][jn][ie];
                  d1nfcy[i2to[jt][in]]=d1nfcy[i2to[jt][in]]-fy[jn]*d[it][jn][ie];
        } } } } }
        /* remove the couple if too far from each other  */
        if((icontact==(-1))&&(dmin2>zone2))
        { if(jcoup<0)
          { i1elcf[ielem]=i1iecn[icoup];
            i1iecn[icoup]=*i1icff;
            *i1icff=icoup;
            icoup=i1elcf[ielem];
			d1sldis[icoup]=R0;
          }
          else
          { i1iecn[jcoup]=i1iecn[icoup];
            i1iecn[icoup]=*i1icff;
            *i1icff=icoup;
            icoup=i1iecn[jcoup];
			d1sldis[icoup]=R0;
        } }
        else
        { jcoup=icoup;
          icoup=i1iecn[icoup];
  } } } }
}

/*********************PUBLIC***********************/
void Yid(  ydc,  yde, ydj, ydi,  ydn, ydp     /***  nodal forces  ***/
        )
  YDC ydc; YDE yde; YDJ ydj; YDI ydi; YDN ydn; YDP ydp;
{ INT iprop,jprop,i,j;
  
  ydi->nicoup=0; /* added by RG the 29/03/06 */
  if(ydi->micoup>0)
  { if(ydi->d2sldis==DBL2NULL)    /* initializing the array of sliding distances */
    { ydi->d2sldis=TalDBL2(ydi->mistate,ydi->micoup);
      if(ydi->d2sldis!=DBL2NULL)
      { for(i=0;i<ydi->micoup;i++)
        { for(j=0;j<ydi->mistate;j++)
          { ydi->d2sldis[j][i]=R0;
    } } } }
	if(yde->d1elsf==DBL1NULL)
	{ yde->d1elsf=TalDBL1(yde->nelem);
	}
	if(yde->d2elcf==DBL2NULL)
	{ yde->d2elcf=TalDBL2(yde->nelem,2);
	}
	for(i=0;i<yde->nelem;i++)
	{ yde->d1elsf[i]=R0;
	  yde->d2elcf[i][0]=R0;
	  yde->d2elcf[i][1]=R0;
	}
	for(iprop=0;iprop<ydp->nprop;iprop++)
    { if(((ydp->i1ptyp[iprop])==(YTE2TRIELS))||
         ((ydp->i1ptyp[iprop])==(YTE2TRISOF))||
         ((ydp->i1ptyp[iprop])==(YTE2TRIRIG)))
      { for(jprop=iprop;jprop<ydp->nprop;jprop++)
        { if(((ydp->i1ptyp[jprop])==(YTE2TRIELS))||
             ((ydp->i1ptyp[jprop])==(YTE2TRISOF))||
             ((ydp->i1ptyp[jprop])==(YTE2TRIRIG)))
          { Yid2TRITRI( /* Triangle to Triange */
            yde->nelem ,
            &(ydi->nicoup),ydc->dctime,/* added by RG the 29/03/06 */
            ydc->dcstec,ydc->icfmty,ydi->diezon,
			iprop,jprop,
            ydi->d1iesl,ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nfcon[0],ydn->d2nfcon[1],
			ydn->d2nvc[0],ydn->d2nvc[1],ydp->d1pefr,ydp->d1pepc,ydp->d1pbif,
			yde->i1elcf,yde->i1elpr,yde->i2eljp,yde->i2elto,yde->i1elty,
			yde->d1elsf,yde->d2elcf,&(ydi->iiecff),ydi->i1iecn,ydi->i1iect,
			ydi->d2sldis,ydj->njoint,ydj->d1jkni,ydj->d1jknc,ydj->d1jksc,
			ydj->d1jphi,ydj->d1jnst,ydj->d1jfpr
            );
  } } } } }
}
/* This Yad module was developed for the topological diagnosis of 
 * interconnected fracture patterns and the initial configuration
 * of fracture apertures
 * (1) the topological analysis is based on a binary-tree search
 *     algorithm
 * (2) the initial apertures can be assigned correlated to trace
 *     lengths or determined empirically by fracture roughenss 
 *     coefficient
 */	// written by Qinghua Lei

/* File   Yad.c */
#include "Yproto.h"

static void Yad2JointDataExtraction(	  /* joint data extraction */
			njoint,i2elto,i1elty,i1jtid,d1ncix,d1nciy,
			i1fj2j,d2fjcx,d2fjcy,d1fjln
  )
  INT  njoint; INT **i2elto; INT  *i1elty; INT *i1jtid; DBL *d1ncix; DBL *d1nciy;
  INT *i1fj2j; DBL **d2fjcx; DBL **d2fjcy; DBL *d1fjln;
{ INT ijoint,ifrjt,inode;
  
  /* extract fracture joint information */
  ifrjt=0;
  for(ijoint=0;ijoint<njoint;ijoint++)
  { if(i1elty[i1jtid[ijoint]]>1)
	{ // coordinates
	  for(inode=0;inode<4;inode++)
	  { d2fjcx[ifrjt][inode]=d1ncix[i2elto[inode][i1jtid[ijoint]]];
		d2fjcy[ifrjt][inode]=d1nciy[i2elto[inode][i1jtid[ijoint]]];
	  }
	  // length of joints
	  V2DLen(d1fjln[ifrjt],d2fjcx[ifrjt][1]-d2fjcx[ifrjt][0],
		d2fjcy[ifrjt][1]-d2fjcy[ifrjt][0]);
	  // index in element list
	  i1fj2j[ifrjt]=ijoint;
	  ifrjt++;
  } }
}

static void Yad2JointConnectivity(	  /* joint connectivity analysis */
			 nfrjt,d2fjcx,d2fjcy,
			d2fjbx,d2fjnm,i2fjcn
  )
  INT    nfrjt; DBL **d2fjcx; DBL **d2fjcy;
  DBL **d2fjbx; DBL **d2fjnm; INT **i2fjcn;
{ DBL xmin, xmax, ymin, ymax, zmin, zmax;
  DBL nx,ny,len;
  INT inode,jnode,ijoint,jjoint;

  /* bounding box */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { xmin=BEPSILON; xmax=-BEPSILON;
	ymin=BEPSILON; ymax=-BEPSILON;
	zmin=BEPSILON; zmax=-BEPSILON;
	for(inode=0;inode<4;inode++)
	{ xmin=MINIM(d2fjcx[ijoint][inode],xmin);
	  xmax=MAXIM(d2fjcx[ijoint][inode],xmax);
	  ymin=MINIM(d2fjcy[ijoint][inode],ymin);
	  ymax=MAXIM(d2fjcy[ijoint][inode],ymax);
	}
	d2fjbx[ijoint][0]=xmin; d2fjbx[ijoint][1]=xmax;
	d2fjbx[ijoint][2]=ymin; d2fjbx[ijoint][3]=ymax;
  }
  /* upward unit normal */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { nx=d2fjcy[ijoint][0]-d2fjcy[ijoint][1];
	ny=d2fjcx[ijoint][1]-d2fjcx[ijoint][0];
	V2DLen(len,nx,ny);
	d2fjnm[ijoint][0]=nx/len;
	d2fjnm[ijoint][1]=ny/len;
  }
  /* connectivity analysis */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { i2fjcn[ijoint][0]=-1;
	i2fjcn[ijoint][1]=-1;
  }
  for(ijoint=0;ijoint<nfrjt-1;ijoint++)
  { for(jjoint=ijoint+1;jjoint<nfrjt;jjoint++)
	{ // bounding box test
	  if(d2fjbx[ijoint][0]>d2fjbx[jjoint][1]+EPSILON
	  || d2fjbx[ijoint][1]<d2fjbx[jjoint][0]-EPSILON
	  || d2fjbx[ijoint][2]>d2fjbx[jjoint][3]+EPSILON
	  || d2fjbx[ijoint][3]<d2fjbx[jjoint][2]-EPSILON)
	  { continue;
	  }
	  // connectivity analysis
	  for(inode=0;inode<2;inode++)
	  { for(jnode=0;jnode<2;jnode++)
		 { // check overlapping
		   if(DABS(d2fjcx[ijoint][inode]-d2fjcx[jjoint][jnode])<EPSILON
		   && DABS(d2fjcy[ijoint][inode]-d2fjcy[jjoint][jnode])<EPSILON)
		   { // connectivity of the 1st joint
			 if(i2fjcn[ijoint][inode]==-1)
			 { i2fjcn[ijoint][inode]=jjoint;	// store neighbour
			 }
			 else if(i2fjcn[ijoint][inode]>=0)
			 { i2fjcn[ijoint][inode]=-2;		// fracture intersection
			 }
			 // connectivity of the 2nd joint
			 if(i2fjcn[jjoint][jnode]==-1)
			 { i2fjcn[jjoint][jnode]=ijoint;	// store neighbour
			 }
			 else if(i2fjcn[jjoint][jnode]>=0)
			 { i2fjcn[jjoint][jnode]=-2;		// fracture intersection
			 }
  } } } } }
}

static void Yad2CountChildren(	  /* count real children of a tree-node */
			i1jtcn,ijoint,i1fj2s,nreal
			)
  INT **i1jtcn; INT ijoint; INT *i1fj2s; INT *nreal;
{ INT inode,jjoint;
  (*nreal)=0;
  if(ijoint!=-1)
  { for(inode=0;inode<2;inode++)
	{ jjoint=i1jtcn[ijoint][inode];
	  if(jjoint>=0 && i1fj2s[jjoint]==-1)
	  { i1fj2s[jjoint]=-2;
		(*nreal)++;
  } } }
}

static void Yad2FindChildren(	/* find real children of a tree-node */
			i1jtcn,ijoint,sectindx,
			nreal,i1fj2s,i1child
			)
  INT **i1jtcn; INT ijoint; INT sectindx;
  INT nreal; INT *i1fj2s; INT *i1child;
{ INT inode, jjoint, ichild;
  ichild=0;
  for(inode=0;inode<2;inode++)
  { jjoint=i1jtcn[ijoint][inode];
	if(jjoint>=0 && i1fj2s[jjoint]<0)
	{ i1child[ichild]=jjoint;
	  i1fj2s[jjoint]=sectindx;
	  ichild++;
  } }
}

static void Yad2LinkJoints2Sections(	  /* identify isolated sections */
			i2fjcn,nfrjt,i1fj2s,nsect)
  INT **i2fjcn; INT nfrjt; INT *i1fj2s; INT *nsect;
{ INT ijoint,jjoint;
  INT *i1child,*i1parent,*i1subchild;
  INT iparent,nparent,ichild,nchild,ireal,*nreal;

  /* initialisation */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { i1fj2s[ijoint]=-1;
  }
  (*nsect)=0;

  /* binary tree search */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { if(i1fj2s[ijoint]==-1)
	{ // a new section
	  i1fj2s[ijoint]=*nsect;
	  nchild=2; nparent=1;
	  i1parent=TalINT1(nparent);
	  *i1parent=ijoint;
	  while(nchild>0)
	  { // count number of real children
		nchild=0;
		nreal=TalINT1(nparent);
		for(iparent=0;iparent<nparent;iparent++)
		{ jjoint=i1parent[iparent];
		  Yad2CountChildren(i2fjcn,jjoint,i1fj2s,&nreal[iparent]);
		  nchild+=nreal[iparent];
		}
		// find and mark real children
		if(nchild>0)
		{ ichild=0;
		  i1child=TalINT1(nchild);
		  for(iparent=0;iparent<nparent;iparent++)
		  { if(nreal[iparent]>0)
			{ jjoint=i1parent[iparent];
			  i1subchild=TalINT1(nreal[iparent]);
			  Yad2FindChildren(i2fjcn,jjoint,*nsect,nreal[iparent],i1fj2s,i1subchild);
			  for(ireal=0;ireal<nreal[iparent];ireal++)
			  { i1child[ichild]=i1subchild[ireal];
				ichild++;
			  }
			  FREE(i1subchild);
		  } }
		  FREE(i1parent);
		  nparent=nchild;
		  i1parent=TalINT1(nparent);
		  i1parent=i1child;
	  } }
	  (*nsect)++;
	  FREE(i1parent);
  } }
}

static void Yad2SectionDataExtraction(	  /* extract isolated section data */
			d1fjln,i2fjcn,nfrjt,i1fj2s,nsect,
			i1scnj,d1scln,i1scbj
			)
  DBL *d1fjln; INT **i2fjcn; INT nfrjt; INT *i1fj2s; INT nsect;
  INT *i1scnj; DBL  *d1scln;  INT *i1scbj;
{ INT inode,ijoint,isect;
  
  /* initialisation */
  for(isect=0;isect<nsect;isect++)
  { i1scnj[isect] = 0;
	d1scln[isect] = R0;
  }
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { i1scbj[ijoint] = 0;
  }

  /* length calculation */
  for(isect=0;isect<nsect;isect++)
  { for(ijoint=0;ijoint<nfrjt;ijoint++)
	{ if(i1fj2s[ijoint]==isect)
	  { d1scln[isect]+=d1fjln[ijoint];
		i1scnj[isect]++;
  } } }

  /* boundary joint identification */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { for(inode=0;inode<2;inode++)
	{ if(i2fjcn[ijoint][inode]<0)
	  { i1scbj[ijoint]=1;
	  }
  } }
}

static void Yad2LinkSections2Fractures(	  /* group sections to form fractures */
			d2fjcx,d2fjcy,i2fjcn,d2fjnm,d2fjbx,nfrjt,
			i1fj2s,i1scnj, nsect,i1scbj,i1sc2f,nfrac
			)
  DBL **d2fjcx; DBL **d2fjcy; INT **i2fjcn; DBL **d2fjnm; DBL **d2fjbx; INT  nfrjt;
  INT  *i1fj2s; INT  *i1scnj; INT    nsect; INT  *i1scbj; INT  *i1sc2f; INT *nfrac;
{ DBL **d2bdbx;			// bounding boxes of sections
  DBL xmin, xmax, ymin, ymax;
  INT isect,jsect,ksect,ijoint,jjoint,inode,jnode,lowID,uppID;
  DBL vecdot;
  INT connflag,isolfrac;

  /* initialisation */
  for(isect=0;isect<nsect;isect++)
  { i1sc2f[isect]=-1;
  }
  (*nfrac)=0;
  d2bdbx=TalDBL2(nsect,4);

  /* only one section */
  if(nsect==1)
  { i1sc2f[0]=0;
	*nfrac=1;
	return;
  }

  /* bounding box */
  for(isect=0;isect<nsect;isect++)
  { xmin=BEPSILON; xmax=-BEPSILON;
	ymin=BEPSILON; ymax=-BEPSILON;
	for(ijoint=0;ijoint<nfrjt;ijoint++)
	{ for(inode=0;inode<2;inode++)
	  { xmin=MINIM(d2fjcx[ijoint][inode],xmin);
		xmax=MAXIM(d2fjcx[ijoint][inode],xmax);
		ymin=MINIM(d2fjcy[ijoint][inode],ymin);
		ymax=MAXIM(d2fjcy[ijoint][inode],ymax);
	} }
	d2bdbx[isect][0]=xmin; d2bdbx[isect][1]=xmax;
	d2bdbx[isect][2]=ymin; d2bdbx[isect][3]=ymax;
  }

  /* section connectivity */
  for(isect=0;isect<nsect-1;isect++)
  { isolfrac=TRUE;
	for(jsect=isect+1;jsect<nsect;jsect++)
	{ // bounding box test
	  if(d2bdbx[isect][0]>d2bdbx[jsect][1]+EPSILON
	  || d2bdbx[isect][1]<d2bdbx[jsect][0]-EPSILON
	  || d2bdbx[isect][2]>d2bdbx[jsect][3]+EPSILON
	  || d2bdbx[isect][3]<d2bdbx[jsect][2]-EPSILON)
	  { continue;
	  }
	  // connected state indicator
	  connflag=FALSE;
	  // recognise connected joint element pairs
	  for(ijoint=0;ijoint<nfrjt;ijoint++)
	  { // not boundary element of the ith section
		if(i1fj2s[ijoint]!=isect||i1scbj[ijoint]==0)
		{ continue;
		}
		for(jjoint=0;jjoint<nfrjt;jjoint++)
		{ // not boundary element of the jth section
		  if(i1fj2s[jjoint]!=jsect||i1scbj[jjoint]==0)
		  { continue;
		  }
		  // bounding box test
		  if(d2fjbx[ijoint][0]>d2fjbx[jjoint][1]+EPSILON
		  || d2fjbx[ijoint][1]<d2fjbx[jjoint][0]-EPSILON
		  || d2fjbx[ijoint][2]>d2fjbx[jjoint][3]+EPSILON
		  || d2fjbx[ijoint][3]<d2fjbx[jjoint][2]-EPSILON)
		  { continue;
		  }
		  // check coplanarity
		  V2DDot(vecdot,d2fjnm[ijoint][0],d2fjnm[ijoint][1],
			d2fjnm[jjoint][0],d2fjnm[jjoint][1]);
		  if(DABS(vecdot)<R1-ANGLETOLER) continue;
		  // check connectivity
		  for(inode=0;inode<2;inode++)
		  { if(i2fjcn[ijoint][inode]>=0) continue;
			for(jnode=0;jnode<2;jnode++)
			{ if(i2fjcn[jjoint][jnode]>=0) continue;
			  if(DABS(d2fjcx[ijoint][inode]-d2fjcx[jjoint][jnode])<EPSILON
			  && DABS(d2fjcy[ijoint][inode]-d2fjcy[jjoint][jnode])<EPSILON)
			  { connflag=TRUE;
				isolfrac=FALSE;
				// update joint connectivity data
				i2fjcn[ijoint][inode]=jjoint;
				i2fjcn[jjoint][jnode]=ijoint;
				break;
			} }
			if(connflag==TRUE) break;
		  }
		  if(connflag==TRUE) break;
		}
		if(connflag==TRUE) break;
	  }
	  if(connflag==TRUE)
	  { // the two sections are connected
		if(i1sc2f[isect]==-1 && i1sc2f[jsect]==-1)
		{ i1sc2f[isect]=*nfrac;
		  i1sc2f[jsect]=*nfrac;
		  (*nfrac)++;
		}
		else if(i1sc2f[isect]==-1)
		{ i1sc2f[isect]=i1sc2f[jsect];
		}
		else if(i1sc2f[jsect]==-1)
		{ i1sc2f[jsect]=i1sc2f[isect];
		}
		else if(i1sc2f[isect]!=i1sc2f[jsect])
		{ lowID=MINIM(i1sc2f[isect],i1sc2f[jsect]);
		  uppID=MAXIM(i1sc2f[isect],i1sc2f[jsect]);
		  i1sc2f[isect]=lowID; i1sc2f[jsect]=lowID;
		  for(ksect=0;ksect<nsect;ksect++)
		  { if(i1sc2f[ksect]==uppID)
			{ i1sc2f[ksect]=lowID;
			}
			else if(i1sc2f[ksect]>uppID)
			{ i1sc2f[ksect]-=1;
		  }	}
		  (*nfrac)--;
	} } }
	// it is an isolated fracture
	if(isolfrac==TRUE && i1sc2f[isect]<0)
	{ i1sc2f[isect]=*nfrac;
	  (*nfrac)++;
	}
	// the last fracture is isolated
	if(connflag==FALSE && isect==nsect-2 && i1sc2f[isect+1]<0)
	{ i1sc2f[isect+1]=*nfrac;
	  (*nfrac)++;
  } }

  /* T type intersection */
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { for(inode=0;inode<2;inode++)
	{ if(i2fjcn[ijoint][inode]==-2)
	  { i2fjcn[ijoint][inode]=-1;
  } } }

  /* memory clean */
  FREE(d2bdbx);
}

static void Yad2FractureDataExtraction(	  /* extract fracture data */
			d1scln,nsect,i1sc2f,nfrac,d1frln
			)
  DBL *d1scln; INT nsect; INT *i1sc2f; INT nfrac; DBL *d1frln;
{ INT isect, ifrac;
  
  /* initialisation */
  for(ifrac=0;ifrac<nfrac;ifrac++)
  { d1frln[ifrac] = 0;
  }

  /* length calculation */
  for(isect=0;isect<nsect;isect++)
  { d1frln[i1sc2f[isect]]+=d1scln[isect];
  }
}

static void Yad2LinkJoints2Fractures(	/* link joint elements to fractures */
			i1fj2s,i2fjcn, nfrjt,i1sc2f,nsect,nfrac,
			i1fj2f,i1frhd,i1frnx
			)
  INT *i1fj2s; INT **i2fjcn; INT   nfrjt; INT *i1sc2f; INT nsect; INT nfrac;
  INT *i1fj2f; INT  *i1frhd; INT *i1frnx;
{ INT ijoint,ifrac,iside;	  // loop variables
  INT *i1frns;				  // number of sections in each fracture
  INT IsEnd;				  // indicator of the end of fracture
  INT preID, curID, nxtID;	  // joint ids

  /* initialisation */
  i1frns = TalINT1(nfrac);
  for(ifrac=0;ifrac<nfrac;ifrac++)
  { i1frns[ifrac]=0;
	i1frhd[ifrac]=-1;
  }
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { i1fj2f[ijoint]=-1;
	i1frnx[ijoint]=-1;
  }
  
  /* link joints to fractures	*/
  for(ijoint=0;ijoint<nfrjt;ijoint++)
  { i1fj2f[ijoint]=i1sc2f[i1fj2s[ijoint]];
  }

  /* count number of joints in each fracture */
  for(ifrac=0;ifrac<nfrac;ifrac++)
  { for(ijoint=0;ijoint<nfrjt;ijoint++)
	{ if(i1fj2f[ijoint]==ifrac)
	  { i1frns[ifrac]++;
  } } }

  /* sort joint order */
  for(ifrac=0;ifrac<nfrac;ifrac++)
  { for(ijoint=0;ijoint<nfrjt;ijoint++)
	{ // if not belonging to this fracture
	  if(i1fj2f[ijoint]!=ifrac) continue;
	  // find header joint
	  if(i2fjcn[ijoint][0]==-1||i2fjcn[ijoint][1]==-1)
	  { if(i1frns[ifrac]==1)
		{ i1frhd[ifrac]=ijoint;
		  i1frnx[ijoint]=-1;
		  break;
		}
		preID=-1;
		curID=ijoint;
		i1frhd[ifrac]=ijoint;
		IsEnd=FALSE;
		while(IsEnd==FALSE)
		{ for(iside=0;iside<2;iside++)
		  { nxtID=i2fjcn[curID][iside];
			if(nxtID!=preID)
			{ if(nxtID<0) IsEnd=TRUE;
			  i1frnx[curID]=nxtID;
			  preID=curID;
			  curID=nxtID;
			  break;
		} } }
		// finish
		break;
  } } }

  /* memory clean */
  FREE(i1frns);
}

static void Yad2GeoFractureJointConfiguration(	/* configure fracture joints geologically */
			i1jtid,i1fj2j,d1fjln,i1fj2s,d1scln,
			 nfrac,i1fj2f,d1frln,i1frhd,i1frnx,
			i1elpr,d1pefr,d1pegt,d1pela,
			d1pemu,d1pjrc,d1pjcs,d1pjsl,
			d1jkni,d1jknc,d1jksc,d1jnst,d1jsst,
			d1japi,d1japc,d1japh,d1japr,
			d1jsdc,d1jdlc,d1jsdp,d1jefl,d1jjrc,
			d1jjcs,d1jphi
			)
  INT *i1jtid; INT *i1fj2j; DBL *d1fjln; INT *i1fj2s; DBL *d1scln;
  INT   nfrac; INT *i1fj2f; DBL *d1frln; INT *i1frhd; INT *i1frnx;
  INT *i1elpr; DBL *d1pefr; DBL *d1pegt; DBL *d1pela;
  DBL *d1pemu; DBL *d1pjrc; DBL *d1pjcs; DBL *d1pjsl;
  DBL *d1jkni; DBL *d1jknc; DBL *d1jksc; DBL *d1jnst; DBL *d1jsst;
  DBL *d1japi; DBL *d1japc; DBL *d1japh; DBL *d1japr;
  DBL *d1jsdc; DBL *d1jdlc; DBL *d1jsdp; DBL *d1jefl; DBL *d1jjrc;
  DBL *d1jjcs; DBL *d1jphi;
{ INT ifrjt,ifrac,ijoint;	  // loop variable
  DBL curvx[2],curvy[2];	  // normalised curvilinear coordinate
  DBL apimax;				  // maximum initial aperture
  DBL aprmax;				  // maximum residual aperture
  DBL vm;					  // maximum normal closure
  INT propID;				  // joint element property id
  DBL E,v,mu,lamda,Gf;		  // mechanical parameters
  DBL JRC0,JCS0,L0,JCSn,JRCn; // joint parameters

  /* calculate apertues of joint elements */
  for(ifrac=0;ifrac<nfrac;ifrac++)
  { ifrjt=i1frhd[ifrac];
	curvx[0]=RP5; curvy[0]=R0;
	while(ifrjt!=-1)
	{ // index of joint element
	  ijoint=i1fj2j[ifrjt];

	  // effective fracture length (i.e. block size)
	  d1jefl[ijoint]=d1scln[i1fj2s[ifrjt]];

	  // mechanical properties
	  propID=i1elpr[i1jtid[ijoint]];
	  mu=d1pemu[propID];
	  lamda=d1pela[propID];
	  Gf=d1pegt[propID];
	  E=mu*(3*lamda+2*mu)/(lamda+mu);
	  v=lamda/2/(lamda+mu);
	  JRC0=d1pjrc[propID];
	  JCS0=d1pjcs[propID];
	  L0=d1pjsl[propID];

	  // scale-dependent JRC and JCS
	  JRCn=JRC0*pow(d1jefl[ijoint]/L0,-0.02*JRC0);
	  JCSn=JCS0*pow(d1jefl[ijoint]/L0,-0.03*JRC0);
	  d1jjrc[ijoint]=JRCn;
	  d1jjcs[ijoint]=JCSn;
	  d1jphi[ijoint]=R0;

	  // maximum initial aperture
	  apimax=SQRT(8*Gf*(1-v*v)/E/MYPI*d1frln[ifrac]);

	  // maximum residual aperture
	  vm=(-0.1032-0.0074*JRCn+1.1350*pow(JCSn/apimax*1e-3,-0.2510))*1e-3;
	  if(vm<R0) vm=apimax*R9/R10;
	  aprmax=apimax-vm;
	  aprmax=MAXIM(apimax/R10,aprmax);

	  // local curvilinear coordinate
	  curvx[1]=curvx[0]-d1fjln[ifrjt]/d1frln[ifrac];

	  // elliptical aperture distribution
	  curvy[1]=SQRT(R1-curvx[1]*curvx[1]);

	  // initial and current aperture
	  // d1japi[ijoint]=RP5*apimax*(curvy[0]+curvy[1]);
	  d1japi[ijoint]=MYPI/R4*apimax;
	  d1japc[ijoint]=d1japi[ijoint];

	  // hydraulic aperture
	  d1japh[ijoint]=pow(d1japi[ijoint]*1e6,R2)/pow(JRCn,2.5)/1e6;;

	  // residual aperture
	  // d1japr[ijoint]=RP5*aprmax*(curvy[0]+curvy[1]);
	  d1japr[ijoint]=MYPI/R4*aprmax;

	  // initial and current joint normal stiffness
	  d1jkni[ijoint]=(-7.15+1.75*JRCn+0.02*JCSn/d1japi[ijoint]*1e-3)*1e9;
	  d1jknc[ijoint]=d1jkni[ijoint];

	  // current joint shear stiffness
	  d1jksc[ijoint]=R0;

	  // normal and shear stress
	  d1jnst[ijoint]=R0; d1jsst[ijoint]=R0;

	  // current and peak shear displacement, and dilation
	  d1jsdc[ijoint]=R0; d1jdlc[ijoint]=R0;
	  d1jsdp[ijoint]=d1jefl[ijoint]/500.0*pow(JRCn/d1jefl[ijoint],0.33);

	  // next joint
	  curvx[0]=curvx[1]; curvy[0]=curvy[1];
	  ifrjt=i1frnx[ifrjt];
  } }
}

static void Yad2EmpFractureJointConfiguration(	/* configure fracture joints empirically */
			i1jtid,i1fj2j,d1fjln,i1fj2s,d1scln,
			 nfrac,i1fj2f,d1frln,i1frhd,i1frnx,
			i1elpr,d1pjrc,d1pjcs,d1pjsl,
			d1jkni,d1jknc,d1jksc,d1jnst,d1jsst,
			d1japi,d1japc,d1japh,d1japr,
			d1jsdc,d1jdlc,d1jsdp,d1jefl,d1jjrc,
			d1jjcs,d1jphi
			)
  INT *i1jtid; INT *i1fj2j; DBL *d1fjln; INT *i1fj2s; DBL *d1scln;
  INT   nfrac; INT *i1fj2f; DBL *d1frln; INT *i1frhd; INT *i1frnx;
  INT *i1elpr; DBL *d1pjrc; DBL *d1pjcs; DBL *d1pjsl;
  DBL *d1jkni; DBL *d1jknc; DBL *d1jksc; DBL *d1jnst; DBL *d1jsst;
  DBL *d1japi; DBL *d1japc; DBL *d1japh; DBL *d1japr;
  DBL *d1jsdc; DBL *d1jdlc; DBL *d1jsdp; DBL *d1jefl; DBL *d1jjrc;
  DBL *d1jjcs; DBL *d1jphi;
{ INT ifrjt,ifrac,ijoint;	  // loop variable
  DBL a0;					  // initial aperture
  DBL ar;					  // residual aperture
  DBL vm;					  // maximum normal closure
  INT propID;				  // joint element property id
  DBL JRC0,JCS0,L0,JCSn,JRCn; // joint parameters

  /* calculate apertues of joint elements */
  for(ifrac=0;ifrac<nfrac;ifrac++)
  { ifrjt=i1frhd[ifrac];
	while(ifrjt!=-1)
	{ // index of joint element
	  ijoint=i1fj2j[ifrjt];

	  // effective fracture length (i.e. block size)
	  d1jefl[ijoint]=d1scln[i1fj2s[ifrjt]];

	  // mechanical properties
	  propID=i1elpr[i1jtid[ijoint]];
	  JRC0=d1pjrc[propID];
	  JCS0=d1pjcs[propID];
	  L0=d1pjsl[propID];

	  // scale-dependent JRC and JCS
	  JRCn=JRC0*pow(d1jefl[ijoint]/L0,-0.02*JRC0);
	  JCSn=JCS0*pow(d1jefl[ijoint]/L0,-0.03*JRC0);
	  // JRCn=JRC0; JCSn=JCS0;
	  d1jjrc[ijoint]=JRCn;
	  d1jjcs[ijoint]=JCSn;
	  d1jphi[ijoint]=R0;

	  // initial aperture
	  a0=JRC0/50/1e3;

	  // closure
	  vm=(-0.1032-0.0074*JRCn+1.1350*pow(JCSn/a0*1e-3,-0.2510))*1e-3;
	  if(vm<R0) vm=a0*R9/R10;

	  // residual aperture
	  ar=a0-vm;

	  // initial and current aperture
	  d1japi[ijoint]=a0;
	  d1japc[ijoint]=a0;

	  // hydraulic aperture
	  d1japh[ijoint]=pow(a0*1e6,R2)/pow(JRCn,2.5)/1e6;

	  // residual aperture
	  d1japr[ijoint]=ar;

	  // initial and current joint normal stiffness
	  d1jkni[ijoint]=(-7.15+1.75*JRCn+0.02*JCSn/d1japi[ijoint]*1e-3)*1e9;
	  d1jknc[ijoint]=d1jkni[ijoint];

	  // current joint shear stiffness
	  d1jksc[ijoint]=R0;

	  // normal and shear stress
	  d1jnst[ijoint]=R0; d1jsst[ijoint]=R0;

	  // current and peak shear displacement, and dilation
	  d1jsdc[ijoint]=R0; d1jdlc[ijoint]=R0;
	  d1jsdp[ijoint]=d1jefl[ijoint]/500.0*pow(JRCn/d1jefl[ijoint],0.33);

	  // next joint
	  ifrjt=i1frnx[ifrjt];
  } }
}

static void Yad2UnfracturedJointConfiguration(	/* configure unfractured joints */
			njoint,i1elty,i1elpr,i2elto,
			d1ncix,d1nciy,
			d1pjrc,d1pjcs,d1pjsl,
			i1jtid,d1jkni,d1jknc,d1jksc,
			d1jnst,d1jsst,d1japi,
			d1japc,d1japh,d1japr,d1jsdc,
			d1jdlc,d1jsdp,d1jefl,d1jjrc,
			d1jjcs,d1jphi
			)
  INT  njoint; INT *i1elty; INT *i1elpr; INT **i2elto;
  DBL *d1ncix; DBL *d1nciy;
  DBL *d1pjrc; DBL *d1pjcs; DBL *d1pjsl;
  INT *i1jtid; DBL *d1jkni; DBL *d1jknc; DBL  *d1jksc;
  DBL *d1jnst; DBL *d1jsst; DBL  *d1japi;
  DBL *d1japc; DBL *d1japh; DBL *d1japr; DBL  *d1jsdc;
  DBL *d1jdlc; DBL *d1jsdp; DBL *d1jefl; DBL  *d1jjrc;
  DBL *d1jjcs; DBL *d1jphi;
{ INT ijoint,ielem,propID;
  DBL JRC0,JCS0;

  for(ijoint=0;ijoint<njoint;ijoint++)
  { ielem=i1jtid[ijoint];
	if(i1elty[ielem]==0)
	{ V2DLen(d1jefl[ijoint],d1ncix[i2elto[0][ielem]]-d1ncix[i2elto[1][ielem]],
		d1nciy[i2elto[0][ielem]]-d1nciy[i2elto[1][ielem]]);
	  propID=i1elpr[ielem];
	  JRC0=d1pjrc[propID];
	  JCS0=d1pjcs[propID];
	  d1jefl[ijoint]=d1pjsl[propID];
	  d1jjrc[ijoint]=JRC0;
	  d1jjcs[ijoint]=JCS0;
	  d1jphi[ijoint]=R0;
	  d1japi[ijoint]=d1jjrc[ijoint]/50.0*1e-3;
	  d1japc[ijoint]=d1japi[ijoint];
	  d1japh[ijoint]=R0;
	  d1japr[ijoint]=MAXIM(d1japi[ijoint]/R10,(-0.296-0.0056*d1jjrc[ijoint]+
		2.241*pow(d1jjcs[ijoint]/d1japi[ijoint]*1e-3,-0.245))*1e-3);
	  d1jkni[ijoint]=(-7.15+1.75*d1jjrc[ijoint]+
		0.02*d1jjcs[ijoint]/d1japi[ijoint]*1e-3)*1e9;
	  d1jknc[ijoint]=d1jkni[ijoint];
	  d1jksc[ijoint]=R0; d1jnst[ijoint]=R0;
	  d1jsst[ijoint]=R0; d1jsdc[ijoint]=R0;
	  d1jdlc[ijoint]=R0;
	  d1jsdp[ijoint]=d1jefl[ijoint]/500.0*pow(d1pjrc[propID]/d1jefl[ijoint],0.33);
	}
  }
}

static void Yad2JointConfiguration(	  /* configure joint element properties */
			iciaty,
			njoint,i2elto,i1elpr,i1elty,
			i1jtid,d1jkni,d1jknc,d1jksc,
			d1jnst,d1jsst,d1japi,
			d1japc,d1japh,d1japr,d1jsdc,
			d1jdlc,d1jsdp,d1jefl,d1jjrc,
			d1jjcs,d1jphi,
			d1ncix,d1nciy,
			d1pefr,d1pegt,d1pela,d1pemu,
			d1pjrc,d1pjcs,d1pjsl
			)
  INT   iciaty;
  INT   njoint; INT **i2elto; INT *i1elpr; INT *i1elty;
  INT  *i1jtid; DBL  *d1jkni; DBL *d1jknc; DBL *d1jksc;
  DBL  *d1jnst; DBL  *d1jsst; DBL *d1japi;
  DBL  *d1japc; DBL  *d1japh; DBL *d1japr; DBL *d1jsdc;
  DBL  *d1jdlc; DBL  *d1jsdp; DBL *d1jefl; DBL *d1jjrc;
  DBL  *d1jjcs; DBL  *d1jphi;
  DBL  *d1ncix; DBL  *d1nciy;
  DBL  *d1pefr; DBL  *d1pegt; DBL *d1pela;
  DBL  *d1pemu; DBL  *d1pjrc; DBL *d1pjsl;
{ INT   ijoint;				// joint loop variable
  INT    nfrjt;				// number of fracture joints
  INT  *i1fj2j;				// fracture joint to joint element indicator
  DBL **d2fjcx,**d2fjcy;	// fracture joint coordinates
  DBL  *d1fjln;				// fracture joint length
  INT **i2fjcn;				// fracture joint connectivity matrix
  DBL **d2fjbx;				// fracture joint bounding boxes
  DBL **d2fjnm;				// fracture joint upward unit normal
  INT    nsect;				// number of isolated sections
  INT  *i1fj2s;				// fracture joint to isolated section indicator
  INT  *i1scnj;				// number of joints in each isolated section
  INT  *i1scbj;				// isolated section boundary joints
  DBL  *d1scln;				// isolated section length
  INT  *i1sc2f;				// isolated section to fracture indicator
  INT    nfrac;				// number of fractures
  INT  *i1fj2f;				// fracture joint to fracture indicator
  DBL  *d1frln;				// length of each fracture
  INT  *i1frhd;				// fracture head joint
  INT  *i1frnx;				// fracture next joint

  /******************************************/
  /*	  Identify fracture geometries		*/
  /******************************************/
  /* count fracture joints */
  nfrjt=0;
  for(ijoint=0;ijoint<njoint;ijoint++)
  { if(i1elty[i1jtid[ijoint]]>1)
	{ nfrjt++;
  } }

  /* joint data extraction */
  i1fj2j=TalINT1(nfrjt);
  d1fjln=TalDBL1(nfrjt);
  d2fjcx=TalDBL2(nfrjt,4);
  d2fjcy=TalDBL2(nfrjt,4);
  Yad2JointDataExtraction(
	njoint,i2elto,i1elty,i1jtid,d1ncix,d1nciy,
	i1fj2j,d2fjcx,d2fjcy,d1fjln);
  
  /* joint connectivity analysis */
  d2fjbx=TalDBL2(nfrjt,4);
  d2fjnm=TalDBL2(nfrjt,2);
  i2fjcn=TalINT2(nfrjt,2);
  Yad2JointConnectivity(
	nfrjt ,d2fjcx,d2fjcy,
	d2fjbx,d2fjnm,i2fjcn);

  /* joint and isolated section linkage */
  i1fj2s=TalINT1(nfrjt);
  Yad2LinkJoints2Sections(
	i2fjcn,nfrjt,i1fj2s,&nsect);

  /* isolated section data extraction */
  i1scnj=TalINT1(nsect);
  d1scln=TalDBL1(nsect);
  i1scbj=TalINT1(nfrjt);
  Yad2SectionDataExtraction(
	d1fjln,i2fjcn,nfrjt ,i1fj2s,nsect,
	i1scnj,d1scln,i1scbj);

  /* section and fracture linkage */
  i1sc2f=TalINT1(nsect);
  Yad2LinkSections2Fractures(
	d2fjcx,d2fjcy,i2fjcn,d2fjnm,d2fjbx,nfrjt,
	i1fj2s,i1scnj, nsect,i1scbj,i1sc2f,&nfrac);

  /* fracture data extraction */
  d1frln=TalDBL1(nfrac);
  Yad2FractureDataExtraction(
	d1scln,nsect,i1sc2f,nfrac,d1frln);

  /* joint element and fracture linkage */
  i1fj2f=TalINT1(nfrjt); i1frhd=TalINT1(nfrac); i1frnx=TalINT1(nfrjt);
  Yad2LinkJoints2Fractures(
	i1fj2s,i2fjcn,nfrjt ,i1sc2f,nsect,nfrac,
	i1fj2f,i1frhd,i1frnx);
  
  /******************************************/
  /*	Configure fracture joint elements	*/
  /******************************************/
  if(iciaty==0)
  { /* roughness-correlated */
	Yad2EmpFractureJointConfiguration(
	  i1jtid,i1fj2j,d1fjln,i1fj2s,d1scln,
	  nfrac ,i1fj2f,d1frln,i1frhd,i1frnx,
	  i1elpr,d1pjrc,d1pjcs,d1pjsl,
	  d1jkni,d1jknc,d1jksc,d1jnst,d1jsst,
	  d1japi,d1japc,d1japh,d1japr,
	  d1jsdc,d1jdlc,d1jsdp,d1jefl,d1jjrc,
	  d1jjcs,d1jphi);
  }
  else if(iciaty==1)
  { /* length-correlated */ // for fractures, not for bedding planes bounded by two materials
	Yad2GeoFractureJointConfiguration(
	  i1jtid,i1fj2j,d1fjln,i1fj2s,d1scln,
	  nfrac ,i1fj2f,d1frln,i1frhd,i1frnx,
	  i1elpr,d1pefr,d1pegt,d1pela,
	  d1pemu,d1pjrc,d1pjcs,d1pjsl,
	  d1jkni,d1jknc,d1jksc,d1jnst,d1jsst,
	  d1japi,d1japc,d1japh,d1japr,
	  d1jsdc,d1jdlc,d1jsdp,d1jefl,d1jjrc,
	  d1jjcs,d1jphi);
  }

  /******************************************/
  /*	Configure unbroken joint elements	*/
  /******************************************/
  Yad2UnfracturedJointConfiguration(
	njoint,i1elty,i1elpr,i2elto,
	d1ncix,d1nciy,
	d1pjrc,d1pjcs,d1pjsl,
	i1jtid,d1jkni,d1jknc,d1jksc,
	d1jnst,d1jsst,d1japi,
	d1japc,d1japh,d1japr,d1jsdc,
	d1jdlc,d1jsdp,d1jefl,d1jjrc,
	d1jjcs,d1jphi);

  /* memory clean */
  FREE(d2fjcx);
  FREE(d2fjcy);
  FREE(d1fjln);
  FREE(i2fjcn);
  FREE(d2fjbx);
  FREE(d2fjnm);
  FREE(i1fj2s);
  FREE(i1scnj);
  FREE(i1scbj);
  FREE(d1scln);
  FREE(i1sc2f);
  FREE(d1frln);
  FREE(i1fj2f);
  FREE(i1frhd);
  FREE(i1frnx);
}

/*********************PUBLIC**********************/
void Yad( ydc, yde, ydj, ydn, ydp		/***  aperture configuration ***/
		)
  YDC ydc; YDE yde; YDJ ydj; YDN ydn; YDP ydp;
{ INT ielem,ijoint;		// loop variables
  INT njoint;			// number of joint elements
  INT ntrele;			// number of triangle elements
  
  if(ydc->ncstep==0)
  { /* initialisation of joint data structure */
	njoint=0;
	for(ielem=0;ielem<yde->nelem;ielem++)
	{ if(yde->i1elty[ielem]!=-1)
	  { njoint++;
	} }
	ntrele=(yde->nelem)-njoint;
	ydj->njoint=njoint;
	ydj->i1jtid=TalINT1(njoint);
	ydj->d1jkni=TalDBL1(njoint);
	ydj->d1jknc=TalDBL1(njoint);
	ydj->d1jksc=TalDBL1(njoint);
	ydj->d1jnst=TalDBL1(njoint);
	ydj->d1jsst=TalDBL1(njoint);
	ydj->d1japi=TalDBL1(njoint);
	ydj->d1japc=TalDBL1(njoint);
	ydj->d1japh=TalDBL1(njoint);
	ydj->d1japr=TalDBL1(njoint);
	ydj->d1jsdc=TalDBL1(njoint);
	ydj->d1jdlc=TalDBL1(njoint);
	ydj->d1jsdp=TalDBL1(njoint);
	ydj->d1jefl=TalDBL1(njoint);
	ydj->d1jjrc=TalDBL1(njoint);
	ydj->d1jjcs=TalDBL1(njoint);
	ydj->d1jphi=TalDBL1(njoint);
	ydj->d1jfmd=TalDBL1(njoint);
	ydj->d1jfpr=TalDBL1(njoint);
	ydj->d1jfet=TalDBL1(njoint);
	for(ijoint=0;ijoint<njoint;ijoint++)
	{ ydj->i1jtid[ijoint]=ntrele+ijoint;
	  ydj->d1jkni[ijoint]=BEPSILON;
	  ydj->d1jknc[ijoint]=BEPSILON;
	  ydj->d1jksc[ijoint]=R0;
	  ydj->d1jnst[ijoint]=R0;
	  ydj->d1jsst[ijoint]=R0;
	  ydj->d1japi[ijoint]=R0;
	  ydj->d1japc[ijoint]=R0;
	  ydj->d1japh[ijoint]=R0;
	  ydj->d1japr[ijoint]=R0;
	  ydj->d1jsdc[ijoint]=R0;
	  ydj->d1jdlc[ijoint]=R0;
	  ydj->d1jsdp[ijoint]=R0;
	  ydj->d1jefl[ijoint]=R0;
	  ydj->d1jjrc[ijoint]=R0;
	  ydj->d1jjcs[ijoint]=R0;
	  ydj->d1jphi[ijoint]=R0;
	  ydj->d1jfmd[ijoint]=R0;
	  ydj->d1jfpr[ijoint]=R0;
	  ydj->d1jfet[ijoint]=R0;
	}

	/* configuration of joint element properties */
	Yad2JointConfiguration(	/* configure joint properties */
	  ydc->iciaty,
	  njoint,yde->i2elto,yde->i1elpr,yde->i1elty,
	  ydj->i1jtid,ydj->d1jkni,ydj->d1jknc,ydj->d1jksc,
	  ydj->d1jnst,ydj->d1jsst,ydj->d1japi,ydj->d1japc,
	  ydj->d1japh,ydj->d1japr,ydj->d1jsdc,ydj->d1jdlc,
	  ydj->d1jsdp,ydj->d1jefl,ydj->d1jjrc,ydj->d1jjcs,
	  ydj->d1jphi,
	  ydn->d2nci[0],ydn->d2nci[1],
	  ydp->d1pefr,ydp->d1pegt,ydp->d1pela,ydp->d1pemu,
	  ydp->d1pjrc,ydp->d1pjcs,ydp->d1pjsl);
  }
}
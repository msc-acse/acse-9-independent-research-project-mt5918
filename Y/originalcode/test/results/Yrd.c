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
/* File  Yrd.c */
#include "Yproto.h"

static void Yrdc(ydc,finp,name)				/* read control data */
  YDC ydc; FILE *finp; CHR *name;
{ CHR *namep;
  
  namep=name+8;
  if(CHRcmp(namep,"MCSTEP",6)==0)
  { INTr(finp,&(ydc->mcstep));
  }
  else if(CHRcmp(namep,"NCSTEP",6)==0)
  { INTr(finp,&(ydc->ncstep));
  }
  else if(CHRcmp(namep,"DCGRAY",6)==0)
  { DBLr(finp,&(ydc->dcgray));
  }
  else if(CHRcmp(namep,"DCSIZC",6)==0)
  { DBLr(finp,&(ydc->dcsizc));
  }
  else if(CHRcmp(namep,"DCSIZF",6)==0)
  { DBLr(finp,&(ydc->dcsizf));
  }
  else if(CHRcmp(namep,"DCSIZS",6)==0)
  { DBLr(finp,&(ydc->dcsizs));
  }
  else if(CHRcmp(namep,"DCSIZV",6)==0)
  { DBLr(finp,&(ydc->dcsizv));
  }
  else if(CHRcmp(namep,"DCSIZD",6)==0) 
  { DBLr(finp,&(ydc->dcsizd));
  }
  else if(CHRcmp(namep,"DCSIZA",6)==0) 
  { DBLr(finp,&(ydc->dcsiza));
  }
  else if(CHRcmp(namep,"DCSTEC",6)==0)
  { DBLr(finp,&(ydc->dcstec));
  }
  else if(CHRcmp(namep,"DCTIME",6)==0)
  { DBLr(finp,&(ydc->dctime));
  }
  else if(CHRcmp(namep,"DCGRST",6)==0) 
  { DBLr(finp,&(ydc->dcgrst));
  }
  else if(CHRcmp(namep,"DCRMPT",6)==0) 
  { DBLr(finp,&(ydc->dcrmpt));
  }
  else if(CHRcmp(namep,"ICOUTF",6)==0)
  { INTr(finp,&(ydc->icoutf));
  }
  else if(CHRcmp(namep,"ICOUTI",6)==0)
  { INTr(finp,&(ydc->icouti));
  }
  else if(CHRcmp(namep,"ICSAVF",6)==0)
  { INTr(finp,&(ydc->icsavf));
  }
  else if(CHRcmp(namep,"ICOUTP",6)==0)
  { INTr(finp,&(ydc->icoutp));
  }
  else if(CHRcmp(namep,"ICFMTY",6)==0)
  { INTr(finp,&(ydc->icfmty));
  }
  else if(CHRcmp(namep,"ICIATY",6)==0)
  { INTr(finp,&(ydc->iciaty));
  }
  else
  { CHRw(stderr,"Yrdc: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(103);
} }
 
static void Yrdd(yd) /* default values */
   YD yd;
{ YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDJ ydj=&(yd->ydj);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);

  /* Set Control to default   */
  ydc->mcstep=0;
  ydc->ncstep=0;
  ydc->finp=(FILE*)NULL;
  ydc->fcheck=(FILE*)NULL;
  ydc->dcgray=R0;
  ydc->dcsizc=R1;
  ydc->dcsizf=R1;
  ydc->dcsizv=R1;
  ydc->dcstec=R1;
  ydc->dctime=R0;
  ydc->icoutf=0;
  ydc->icouti=0;
  ydc->icsavf=0;
  ydc->icoutp=2;
  ydc->icfmty=0;
  ydc->iciaty=0;

  /* Set Elements to default */
  yde->melem=0; yde->nelem=0;
  yde->melst=0; yde->nelst=0;
  yde->melno=0; yde->nelno=0;
  yde->i1elcf=INT1NULL; 
  yde->i1elpr=INT1NULL;
  yde->i1elty=INT1NULL;
  yde->d2elst=DBL2NULL;
  yde->i2elto=INT2NULL;
  yde->i2eljp=INT2NULL;
  yde->d2tcs =DBL2NULL;
  yde->d1elsf=DBL1NULL;
  yde->d2elcf=DBL1NULL;
  yde->d2elfs=DBL2NULL;
  yde->d1eley=DBL1NULL;
  yde->d1eles=DBL1NULL;
  yde->d1elme=DBL1NULL;
  yde->i1elyi=INT1NULL;

  /* Set Joints to default */	// added by Qinghua
  ydj->njoint=0;
  ydj->i1jtid=INT1NULL;
  ydj->d1jkni=DBL1NULL;
  ydj->d1jknc=DBL1NULL;
  ydj->d1jksc=DBL1NULL;
  ydj->d1jnst=DBL1NULL;
  ydj->d1jsst=DBL1NULL;
  ydj->d1japi=DBL1NULL;
  ydj->d1japc=DBL1NULL;
  ydj->d1japh=DBL1NULL;
  ydj->d1japr=DBL1NULL;
  ydj->d1jsdc=DBL1NULL;
  ydj->d1jdlc=DBL1NULL;
  ydj->d1jsdp=DBL1NULL;
  ydj->d1jefl=DBL1NULL;
  ydj->d1jjrc=DBL1NULL;
  ydj->d1jjcs=DBL1NULL;
  ydj->d1jphi=DBL1NULL;
  ydj->d1jfmd=DBL1NULL;
  ydj->d1jfpr=DBL1NULL;
  ydj->d1jfet=DBL1NULL;

  /* Set Interaction to default */
  ydi->micoup=0; ydi->nicoup=0;  
  ydi->iiecff=-2;
  ydi->diedi=BEPSILON;
  ydi->diezon=R0;
  ydi->d1iesl=DBL1NULL;
  ydi->i1iecn=INT1NULL;
  ydi->i1iect=INT1NULL;
  ydi->mistate=6;
  ydi->d2sldis=DBL2NULL;

  /* Set Nodes to default */
  ydn->mnodim=0; ydn->nnodim=0;
  ydn->mnopo=0;  ydn->nnopo=0;
  ydn->d1nmct=DBL1NULL;
  ydn->d2ncc=DBL2NULL;
  ydn->d2nci=DBL2NULL;
  ydn->d2nfc=DBL2NULL;
  ydn->d2nfcon=DBL2NULL;
  ydn->d2nvc=DBL2NULL;
  ydn->i1nobf=INT1NULL;
  ydn->i1nopr=INT1NULL;

  /* Set Output to default  */
  ydo->mohys=0;  ydo->nohys=0;
  ydo->dohyp=0.05;  /* 5% accuracy */
  ydo->d1ohyc=DBL1NULL;
  ydo->d1ohyf=DBL1NULL;
  ydo->d1ohys=DBL1NULL;
  ydo->d1ohyt=DBL1NULL;
  ydo->d1ohyx=DBL1NULL;
  ydo->d1ohyy=DBL1NULL;
  ydo->i1ohyt=INT1NULL;

  /* Set Properties to default  */
  ydp->mprop=0; ydp->nprop=0;
  ydp->d1pefr=DBL1NULL;
  ydp->d1peca=DBL1NULL;
  ydp->d1pecl=DBL1NULL;
  ydp->d1peft=DBL1NULL;
  ydp->d1pegt=DBL1NULL;
  ydp->d1pegs=DBL1NULL;
  ydp->d1peks=DBL1NULL;
  ydp->d1pela=DBL1NULL;
  ydp->d1pemu=DBL1NULL;
  ydp->d1pepe=DBL1NULL;
  ydp->d1pepc=DBL1NULL;
  ydp->d1pepf=DBL1NULL;
  ydp->d1pbif=DBL1NULL;
  ydp->d1pera=DBL1NULL;
  ydp->d1pero=DBL1NULL;
  ydp->d1pevi=DBL1NULL;
  ydp->d1picf=DBL1NULL;
  ydp->d1pcoh=DBL1NULL;
  ydp->d1pnaf=DBL1NULL;
  ydp->d1pnai=DBL1NULL;
  ydp->d1pnap=DBL1NULL;
  ydp->d1pnat=DBL1NULL;
  ydp->d1pnax=DBL1NULL;
  ydp->d1pnay=DBL1NULL;
  ydp->d1pnaz=DBL1NULL;
  ydp->d1pnxx=DBL1NULL;
  ydp->d1pnxy=DBL1NULL;
  ydp->d1pnxz=DBL1NULL;
  ydp->d1pnyx=DBL1NULL;
  ydp->d1pnyy=DBL1NULL;
  ydp->d1pnyz=DBL1NULL;
  ydp->d1pnzx=DBL1NULL;
  ydp->d1pnzy=DBL1NULL;
  ydp->d1pnzz=DBL1NULL;
  ydp->d1psem=DBL1NULL;
  ydp->d1pjrc=DBL1NULL;
  ydp->d1pjcs=DBL1NULL;
  ydp->d1pjsl=DBL1NULL;
  ydp->i1pecn=INT1NULL;
  ydp->i1pecp=INT1NULL;
  ydp->i1pect=INT1NULL;
  ydp->i1pefr=INT1NULL;
  ydp->i1pejp=INT1NULL;
  ydp->i1pemb=INT1NULL;
  ydp->i1pemn=INT1NULL;
  ydp->i1pnfx=INT1NULL;
  ydp->i1pnfy=INT1NULL;
  ydp->i1pnfz=INT1NULL;
  ydp->i1pnib=INT1NULL;
  ydp->i1psde=INT1NULL;
  ydp->i1ptyp=INT1NULL;
}

static void Yrde(yde,finp,name)			/* read elements data */
  YDE    yde; FILE *finp; CHR *name; 
{ CHR *namep;
  INT i,j;
  INT iinit=-1;

  namep=name+8;
  if(CHRcmp(namep,"MELEM",5)==0)
  { INTr(finp,&(yde->melem));
  }
  else if(CHRcmp(namep,"NELEM",5)==0)
  { INTr(finp,&(yde->nelem));
  }
  else if(CHRcmp(namep,"MELST",5)==0)
  { INTr(finp,&(yde->melst));
  }
  else if(CHRcmp(namep,"NELST",5)==0)
  { INTr(finp,&(yde->nelst));
  }
  else if(CHRcmp(namep,"MELNO",5)==0)
  { INTr(finp,&(yde->melno));
  }
  else if(CHRcmp(namep,"NELNO",5)==0)
  { INTr(finp,&(yde->nelno));
  }
  else if(CHRcmp(namep,"I1ELCF",6)==0)
  { TformINT1(finp,iinit,yde->melem,&(yde->i1elcf));
  }
  else if(CHRcmp(namep,"I1ELPR",6)==0)
  { TformINT1(finp,iinit,yde->melem,&(yde->i1elpr));
  }
  else if(CHRcmp(namep,"I1ELTY",6)==0)
  { TformINT1(finp,iinit,yde->melem,&(yde->i1elty));
  }
  else if(CHRcmp(namep,"D2ELST",6)==0)
  { TformDBL2(finp,R0,yde->melst,yde->melem,&(yde->d2elst));
  }
  else if(CHRcmp(namep,"I2ELTO",6)==0)
  { TformINT2(finp,iinit,yde->melno,yde->melem,&(yde->i2elto));
    yde->d2tcs=TalDBL2(4,yde->melem);
	for(i=0;i<yde->melem;i++)
	{ for(j=0; j<4; j++) 
	  { yde->d2tcs[j][i]=R0;
	  }
	}
  }
  else if(CHRcmp(namep,"D2ELFS",6)==0)
  { TformDBL2(finp,R0,yde->melno,yde->melem,&(yde->d2elfs));
  }
  else if(CHRcmp(namep,"I2ELJP",6)==0)
  { TformINT2(finp,iinit,yde->melno,yde->melem,&(yde->i2eljp));
  }
  else if(CHRcmp(namep,"D1ELEY",6)==0)
  { TformDBL1(finp,R0,yde->melem,&(yde->d1eley));
  }
  else if(CHRcmp(namep,"D1ELES",6)==0)
  { TformDBL1(finp,R0,yde->melem,&(yde->d1eles));
  }
  else if(CHRcmp(namep,"D1ELME",6)==0)
  { TformDBL1(finp,-BEPSILON,yde->melem,&(yde->d1elme));
  }
  else if(CHRcmp(namep,"I1ELYI",6)==0)
  { TformINT1(finp,iinit,yde->melem,&(yde->i1elyi));
  }
  else
  { CHRw(stderr,"Yrde: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
} }

static void Yrdj(ydj,finp,name)			/* read joint elements data */
  YDJ    ydj; FILE *finp; CHR *name; 
{ CHR *namep;
  INT iinit=-1;

  namep=name+8;
  if(CHRcmp(namep,"NJOINT",6)==0)
  { INTr(finp,&(ydj->njoint));
  }
  else if(CHRcmp(namep,"I1JTID",6)==0)
  { TformINT1(finp,iinit,ydj->njoint,&(ydj->i1jtid));
  }
  else if(CHRcmp(namep,"D1JKNI",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jkni));
  }
  else if(CHRcmp(namep,"D1JKNC",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jknc));
  }
  else if(CHRcmp(namep,"D1JKSC",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jksc));
  }
  else if(CHRcmp(namep,"D1JNST",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jnst));
  }
  else if(CHRcmp(namep,"D1JSST",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jsst));
  }
  else if(CHRcmp(namep,"D1JAPI",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1japi));
  }
  else if(CHRcmp(namep,"D1JAPC",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1japc));
  }
  else if(CHRcmp(namep,"D1JAPH",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1japh));
  }
  else if(CHRcmp(namep,"D1JAPR",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1japr));
  }
  else if(CHRcmp(namep,"D1JSDC",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jsdc));
  }
  else if(CHRcmp(namep,"D1JDLC",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jdlc));
  }
  else if(CHRcmp(namep,"D1JSDP",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jsdp));
  }
  else if(CHRcmp(namep,"D1JEFL",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jefl));
  }
  else if(CHRcmp(namep,"D1JJRC",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jjrc));
  }
  else if(CHRcmp(namep,"D1JJCS",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jjcs));
  }
  else if(CHRcmp(namep,"D1JPHI",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jphi));
  }
  else if(CHRcmp(namep,"D1JFMD",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jfmd));
  }
  else if(CHRcmp(namep,"D1JFPR",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jfpr));
  }
  else if(CHRcmp(namep,"D1JFET",6)==0)
  { TformDBL1(finp,R0,ydj->njoint,&(ydj->d1jfet));
  }
  else
  { CHRw(stderr,"Yrdj: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
} }

static void Yrdi(ydi,finp,name)				/* read interaction data */
  YDI   ydi; FILE *finp; CHR *name; 
{ CHR *namep;
  
  namep=name+8;
  if(CHRcmp(namep,"MICOUP",6)==0)			// maximum number of contacting couples of finite elements
  { INTr(finp,&(ydi->micoup));
  }
  else if(CHRcmp(namep,"NICOUP",6)==0)		// actual number of contacting couples of finite elements
  { INTr(finp,&(ydi->nicoup));
  }
  else if(CHRcmp(namep,"IIECFF",6)==0)		// internal variable for contact
  { INTr(finp,&(ydi->iiecff));
  }
  else if(CHRcmp(namep,"DIEDI",5)==0)		// internal variable to trigger contact detection
  { DBLr(finp,&(ydi->diedi));
  }
  else if(CHRcmp(namep,"DIEZON",6)==0)		// buffer size around each element for contact detection purposes
  { DBLr(finp,&(ydi->diezon));
  }
  else if(CHRcmp(namep,"D1IESL",6)==0)		// array for storing sliding distance between couples in contact
  /* size should be melem != micoup */
  { TformDBL1(finp,R0,ydi->micoup,&(ydi->d1iesl)); 
  }
  else if(CHRcmp(namep,"I1IECN",6)==0)		// array to store next contacting couple in the list
  { TformINT1(finp,-1,ydi->micoup,&(ydi->i1iecn)); 
  }
  else if(CHRcmp(namep,"I1IECT",6)==0)		// array to store target finite element for each contacting couple
  { TformINT1(finp,-1,ydi->micoup,&(ydi->i1iect)); 
  }
  else
  { CHRw(stderr,"Yrdi: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
} }

static void Yrdn(ydn,finp,name)				/* read nodes data */
  YDN    ydn; FILE *finp; CHR *name; 
{ CHR *namep;
  INT i;
  
  namep=name+8;
  if(CHRcmp(namep,"MNODIM",5)==0)			// maximum number of freedom degrees per node
  { INTr(finp,&(ydn->mnodim));
  }
  else if(CHRcmp(namep,"NNODIM",5)==0)		// actual number of freedom degrees per node
  { INTr(finp,&(ydn->nnodim));
  }
  else if(CHRcmp(namep,"MNOPO",5)==0)		// maximum number of nodes
  { INTr(finp,&(ydn->mnopo));
  }
  else if(CHRcmp(namep,"NNOPO",5)==0)		// actual number of nodes
  { INTr(finp,&(ydn->nnopo));
  }
  else if(CHRcmp(namep,"D1NMCT",6)==0)		// array of current mass for each node
  { TformDBL1(finp,R0,ydn->mnopo,&(ydn->d1nmct));
  }
  else if(CHRcmp(namep,"D2NCC",5)==0)		// array of current coordinates of nodes
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2ncc)); 
  }
  else if(CHRcmp(namep,"D2NCI",5)==0)		// array of initial coordinates of nodes
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nci)); 
  }
  else if(CHRcmp(namep,"D2NFC",5)==0)		// array of the current nodal forces
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nfc));
    ydn->d2nfcon=TalDBL2(2,ydn->mnopo);
	for(i=0;i<ydn->mnopo;i++)
	{ ydn->d2nfcon[0][i]=ydn->d2nfc[0][i];
	  ydn->d2nfcon[1][i]=ydn->d2nfc[1][i];
	}
  }
  else if(CHRcmp(namep,"D2NVC",5)==0)		// array of current velocities
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nvc)); 
  }
  else if(CHRcmp(namep,"I1NOBF",6)==0)		// flag array for indicating whether a node is on the boundary
  { TformINT1(finp,1,ydn->mnopo,&(ydn->i1nobf));
  }
  else if(CHRcmp(namep,"I1NOPR",6)==0)		// array of property set ID
  { TformINT1(finp,-1,ydn->mnopo,&(ydn->i1nopr));
  }
  else
  { CHRw(stderr,"Yrdn: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
} }

static void Yrdo(ydo,finp,name)					/* read output data */
  YDO    ydo; FILE *finp; CHR *name; 
{ INT i; CHR *namep;
  
  namep=name+8;
  if(CHRcmp(namep,"MOHYS",5)==0) 
  {  INTr(finp,&(ydo->mohys))
     i=(ydo->mohys)*sizeof(FILE*);
     if(i>0)ydo->f2ohyf=(FILE**)MALLOC(i);
     for(i=0;i<(ydo->mohys);i++)
     { ydo->f2ohyf[i]=FILENULL;
     }
  }
  else if(CHRcmp(namep,"NOHYS",5)==0) 
  { INTr(finp,&(ydo->nohys));
  }
  else if(CHRcmp(namep,"DOHYP",5)==0) 
  { DBLr(finp,&(ydo->dohyp));
  }
  else if(CHRcmp(namep,"D1OHYS",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohys));
  }
  else if(CHRcmp(namep,"D1OHYC",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyc));
  }
  else if(CHRcmp(namep,"D1OHYF",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyf));
  }
  else if(CHRcmp(namep,"D1OHYT",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyt));
  }
  else if(CHRcmp(namep,"D1OHYX",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyx));
  }
  else if(CHRcmp(namep,"D1OHYY",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyy));
  }
  else if(CHRcmp(namep,"D1OHYZ",6)==0) 
  { TformDBL1(finp,R0,ydo->mohys,&(ydo->d1ohyz));
  }
  else if(CHRcmp(namep,"I1OHYT",6)==0) 
  { TformINT1(finp,-1,ydo->mohys,&(ydo->i1ohyt));
  }
  else
  { CHRw(stderr,"Yrdo: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
} }

static void Yrdp(ydp,finp,name)					/* read properties data */
  YDP    ydp; FILE *finp; CHR *name; 
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MPROP",5)==0)				// maximum number of properties sets
  {  INTr(finp,&(ydp->mprop));
  }
  else if(CHRcmp(namep,"NPROP",5)==0)			// actual number of properties sets
  {  INTr(finp,&(ydp->nprop));
  }
  else if(CHRcmp(namep,"D1PECA",6)==0)			// 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1peca)); 
  }
  else if(CHRcmp(namep,"D1PECL",6)==0)			// 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pecl)); 
  }
  else if(CHRcmp(namep,"D1PEFT",6)==0)			// tensile strength
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1peft));
  }
  else if(CHRcmp(namep,"D1PEGT",6)==0)			// mode I energy release rate
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pegt));
  }
  else if(CHRcmp(namep,"D1PEGS",6)==0)			// mode II energy release rate
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pegs));
  }
  else if(CHRcmp(namep,"D1PEKS",6)==0)			// viscous damping
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1peks));
  }
  else if(CHRcmp(namep,"D1PELA",6)==0)			// 1st lame elastic constant
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pela));
  }
  else if(CHRcmp(namep,"D1PEMU",6)==0)			// 2nd lame elastic constant
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pemu));
  }
  else if(CHRcmp(namep,"D1PEPE",6)==0)			// element penalty term
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepe));
  }
  else if(CHRcmp(namep,"D1PEPC",6)==0)			// contact penalty term
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepc));
  }
  else if(CHRcmp(namep,"D1PEPF",6)==0)			// pore fluid pressure
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pepf));
  }
  else if(CHRcmp(namep,"D1PBIF",6)==0)			// bedding interface friction coeff.
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pbif));
  }
  else if(CHRcmp(namep,"D1PCOH",6)==0)			// cohesion
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pcoh));
  }
  else if(CHRcmp(namep,"D1PICF",6)==0)			// internal friction coefficient
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1picf));
  }
  else if(CHRcmp(namep,"D1PERA",6)==0)			// 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pera));
  }
  else if(CHRcmp(namep,"D1PERO",6)==0)			// density
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pero));
  }
  else if(CHRcmp(namep,"D1PEVI",6)==0)			// 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pevi));
  }
  else if(CHRcmp(namep,"D1PEFR",6)==0)			// joint friction coefficient
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pefr));
  }
  else if(CHRcmp(namep,"D1PNAF",6)==0)			// amplitude factor
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnaf));
  }
  else if(CHRcmp(namep,"D1PNAI",6)==0)			// increment of the amplitude factor each time step
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnai));
  }
  else if(CHRcmp(namep,"D1PNAP",6)==0)			// amplitude of load applied as element surface pressure
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnap)); 
  }
  else if(CHRcmp(namep,"D1PNAT",6)==0)			// amplitude of load applied as element surface traction
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnat)); 
  }
  else if(CHRcmp(namep,"D1PNAX",6)==0)			// amplitude of force or velocity in local x direction
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnax)); 
  }
  else if(CHRcmp(namep,"D1PNAY",6)==0)			// amplitude of force or velocity in local y direction
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnay)); 
  }
  else if(CHRcmp(namep,"D1PNAZ",6)==0)			//amplitude of force or velocity in local z direction
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnaz)); 
  }
  else if(CHRcmp(namep,"D1PNXX",6)==0)			// x component of x axis of local coordinate system
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnxx)); 
  }
  else if(CHRcmp(namep,"D1PNXY",6)==0)			// y component of x axis of local coordinate system
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnxy)); 
  }
  else if(CHRcmp(namep,"D1PNXZ",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnxz)); 
  }
  else if(CHRcmp(namep,"D1PNYX",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnyx)); 
  }
  else if(CHRcmp(namep,"D1PNYY",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnyy)); 
  }
  else if(CHRcmp(namep,"D1PNYZ",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnyz)); 
  }
  else if(CHRcmp(namep,"D1PNZX",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnzx)); 
  }
  else if(CHRcmp(namep,"D1PNZY",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnzy)); 
  }
  else if(CHRcmp(namep,"D1PNZZ",6)==0) 
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pnzz)); 
  }
  else if(CHRcmp(namep,"D1PSEM",6)==0)
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1psem)); 
  }
  else if(CHRcmp(namep,"D1PJRC",6)==0)
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pjrc)); 
  }
  else if(CHRcmp(namep,"D1PJCS",6)==0)
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pjcs)); 
  }
  else if(CHRcmp(namep,"D1PJSL",6)==0)
  { TformDBL1(finp,R0,ydp->mprop,&(ydp->d1pjsl)); 
  }
  else if(CHRcmp(namep,"I1PECN",6)==0)
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pecn));
  }
  else if(CHRcmp(namep,"I1PECP",6)==0) 
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pecp));
  }
  else if(CHRcmp(namep,"I1PECT",6)==0)
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pect));
  }
  else if(CHRcmp(namep,"I1PEFR",6)==0)
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pefr));
  }  
  else if(CHRcmp(namep,"I1PEJP",6)==0)
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pejp));
  }
  else if(CHRcmp(namep,"I1PEMB",6)==0)			// flag for boundary nodes marking
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pemb));
  }
  else if(CHRcmp(namep,"I1PEMN",6)==0)			// number of successive mesh refinement
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pemn));
  }
  else if(CHRcmp(namep,"I1PNFX",6)==0)			// boundary condition type: 1-force, 2-acceleration, 3-velocity
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pnfx));
  }
  else if(CHRcmp(namep,"I1PNFY",6)==0) 
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pnfy));
  }
  else if(CHRcmp(namep,"I1PNFZ",6)==0) 
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pemn));
  }
  else if(CHRcmp(namep,"I1PNIB",6)==0) 
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1pnib));
  }
  else if(CHRcmp(namep,"I1PTYP",6)==0)			// type of each property
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1ptyp));
  }
  else if(CHRcmp(namep,"I1PSDE",6)==0)			// ID of elastic damage state variable
  { TformINT1(finp,-1,ydp->mprop,&(ydp->i1psde));
  }
  else
  { CHRw(stderr,"Yrdp: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(104);
} }

/***************************PUBLIC***********************************/
INT Yrd(namep,yd)
   CHR *namep; YD yd;
{ INT icount;
  CHR name[300];
  YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDJ ydj=&(yd->ydj);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);

  if((ydc->finp)==(FILENULL))
  { Yrdd(yd);
    ydc->finp=fopen(namep,"r");
    ydc->fcheck=fopen("Ytmp","w");
  }
  if(((ydc->finp)==(FILENULL))||((ydc->fcheck)==(FILENULL)))
  { CHRw(stderr,"Yrd: Could not open input file - usage -i inputfile"); 
    CHRwcr(stderr);   
    return 0;
  }
 
  SETLINEBUF(ydc->fcheck);
  CHRr(ydc->finp,name);
  while(FILEND(ydc->finp)==0) 
  { if(CHRcmp(name, "$YSTOP",6)==0)
    { CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
      return 0;
    }
    else if(CHRcmp(name, "$YDOIT",6)==0)
    { CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
      return 1;
    }
    else if(CHRcmp(name,"/*",2)==0)    /* read and ignore comments */
    { icount=0;
      do
      { CHRr(ydc->finp,name); icount++;
        if(icount>100)
        { CHRw(stderr,"Yrd: too long comment near - ");
          CHRw(stderr,name);
          CHRwcr(stderr);      
          return 0;
      } }
	  while((FILEND(ydc->finp)==0)&&(CHRcmp(name,"*/",2)!=0));
    }
    else if(CHRcmp(name, "/YD/",4)==0)
    { CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
      if(CHRcmp(name,"/YD/YDC/",8)==0)			/* read control data       */
      { Yrdc(ydc,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDE/",8)==0)		/* read data elements      */
      { Yrde(yde,ydc->finp,name);
      }
	  else if(CHRcmp(name,"/YD/YDJ/",8)==0)		/* read data joints        */
	  { Yrdj(ydj,ydc->finp,name);
	  }
      else if(CHRcmp(name,"/YD/YDI/",8)==0)		/* read data interaction   */
      { Yrdi(ydi,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDN/",8)==0)		/* read data nodes          */
      { Yrdn(ydn,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDO/",8)==0)		/* read data output         */
      { Yrdo(ydo,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDP/",8)==0)		/* read data properties     */
      { Yrdp(ydp,ydc->finp,name);
      }
      else
      { CHRw(stderr,"Yrd: unknown name: ");
        CHRw(stderr,name); 
        CHRwcr(stderr);   
        return 0;
    } }
    else
    { CHRw(stderr,"Yrd: unknown name: ");
      CHRw(stderr,name); 
      CHRwcr(stderr);
      return 0;
    }
    CHRr(ydc->finp,name);
  }
  fclose(ydc->finp);
  fclose(ydc->fcheck);
  return 0;
}

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
static  CHR *cdig="0123456789";

static void Yrdc(ydc,finp,name)
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
  else if(CHRcmp(namep,"DCSTEC",6)==0)
  { DBLr(finp,&(ydc->dcstec));
  }
  else if(CHRcmp(namep,"DCTIME",6)==0)
  { DBLr(finp,&(ydc->dctime));
  }
  else if(CHRcmp(namep,"ICOUTF",6)==0)
  { INTr(finp,&(ydc->icoutf));
  }
  else if(CHRcmp(namep,"ICOUTRF",7)==0)
  { INTr(finp,&(ydc->icoutrf));
  }
  else if(CHRcmp(namep,"ICOUTNF",7)==0)
  { INTr(finp,&(ydc->icoutnf));
  }
  else if(CHRcmp(namep,"ICOUTI",6)==0)
  { INTr(finp,&(ydc->icouti));
  }
  else if(CHRcmp(namep,"ICOUTP",6)==0)
  { INTr(finp,&(ydc->icoutp));
  }
  else if(CHRcmp(namep,"ICRESF",6)==0)
  { INTr(finp,&(ydc->icresf));
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
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDB ydb=&(yd->ydb);
  YDS yds=&(yd->yds);
  YDO ydo=&(yd->ydo);
  YDPE ydpe=&(yd->ydpe);
  YDPN ydpn=&(yd->ydpn);
  YDPJ ydpj=&(yd->ydpj);
  YDPM ydpm=&(yd->ydpm);
  YDFN ydfn=&(yd->ydfn);
  YDIS ydis=&(yd->ydis);
  YDHF ydhf=&(yd->ydhf);
  YDSM ydsm=&(yd->ydsm);
  YDR ydr=&(yd->ydr);
  YDSB ydsb=&(yd->ydsb);
  YDPS ydps=&(yd->ydps);
  
  
  /* Set Control to default   */
  ydc->mcstep=0; ydc->ncstep=0;
  ydc->finp=(FILE*)NULL;
  ydc->fcheck=(FILE*)NULL;
  ydc->dcgray=R0;
  ydc->dcsizc=R1;
  ydc->dcsizf=R1;
  ydc->dcsizv=R1;
  ydc->dcstec=R1;
  ydc->dctime=R0;
  ydc->icoutf=0;
  ydc->icoutrf=0;
  ydc->icoutnf=0;
  ydc->icouti=0;
  ydc->icoutp=2;
  ydc->icresf=0;
  
  /* Set Elements to default    */
  yde->melem=0; yde->nelem=0;
  yde->melst=0; yde->nelst=0;
  yde->melno=0; yde->nelno=0;
  yde->nelemst=0;
  yde->nebrk=0; yde->netbrk=0;
  yde->nesft=0; yde->netsft=0;
  yde->i1elcf=INT1NULL;
  yde->i1elpr=INT1NULL;
  yde->i1elprtmp=INT1NULL;
  yde->d1elfs=DBL1NULL;
  yde->d2elst=DBL2NULL;
  yde->i2elto=INT2NULL;
  yde->i2eltost=INT2NULL;
  yde->i1ebrk=INT1NULL;
  yde->i1esft=INT1NULL;
  yde->d2ecbrk=DBL2NULL;
  yde->d2ecbrk_NEW=DBL2NULL;
  yde->d1etbrk=DBL1NULL;
  yde->d1elbrk=DBL1NULL;
  yde->d1efe=DBL1NULL;
  yde->d2ecsft=DBL2NULL;
  yde->d1etsft=DBL1NULL;
  yde->i1esftf=INT1NULL;
  yde->d1ebrkf=DBL1NULL;
  yde->d1etmke=DBL1NULL;
  yde->d1eike=DBL1NULL;           
  yde->d1edke=DBL1NULL;            
  yde->d1elfr=DBL1NULL; 
  yde->d1elpe=DBL1NULL; 
  yde->d1elpt=DBL1NULL;
  yde->i2elnext=INT2NULL;
  yde->i2eledge=INT2NULL;
  yde->i1edfnf=INT1NULL; 
  yde->i1edft=INT1NULL; 
  yde->d1etike=DBL1NULL;
  yde->d2eldmg=DBL2NULL;
  yde->d2elstr=DBL2NULL;
  
  /* Set Interaction to default */
  ydi->micoup=0; ydi->nicoup=0;
  ydi->iiecff=-2;
  ydi->diedi=BEPSILON;
  ydi->diezon=R0;
  ydi->d1iesl=DBL1NULL;
  ydi->i1iecn=INT1NULL;
  ydi->i1iect=INT1NULL;
  ydi->d2sldis=DBL2NULL;
  /* Set Nodes to default  */
  ydn->mnodim=0;  ydn->nnodim=0;
  ydn->mnopo=0;  ydn->nnopo=0;
  ydn->d1nmct=DBL1NULL;
  ydn->d2ncc=DBL2NULL;
  ydn->d2nci=DBL2NULL;
  ydn->d2nfc=DBL2NULL;
  ydn->d2nvc=DBL2NULL;
  ydn->i1nobf=INT1NULL;
  ydn->i1nopr=INT1NULL;
  ydn->i1nowe=INT1NULL; //! for hydrofrac
  ydn->i2noid=INT2NULL; //! for hydrofrac
  ydn->d1nfp=DBL1NULL;
  ydn->d2nc0=DBL2NULL;
  /* Set Boreholes to default */
  ydb->mborh=0; ydb->nborh=0;
  ydb->mbdim=0; ydb->nbdim=0;
  ydb->nbpaf=0;
  ydb->d2bca=DBL2NULL;
  ydb->d2bcb=DBL2NULL;
  ydb->d1brad=DBL1NULL;
  ydb->d1bpaf=DBL1NULL;
  ydb->d1bpts=DBL1NULL;
  ydb->d1bpte=DBL1NULL;
  ydb->d1bvdt=DBL1NULL;
  ydb->d1bprs=DBL1NULL;
  ydb->dblmax=R0;
  ydb->dbbuf=R0;
  /* Set Sources to default */
  yds->msour=0; yds->nsour=0;
  yds->msdim=0; yds->nsdim=0;
  yds->nspaf=0;
  yds->nssaf=0;
  yds->d2scs=DBL2NULL;
  yds->d1spaf=DBL1NULL;
  yds->d1ssaf=DBL1NULL;
  yds->d1spts=DBL1NULL;
  yds->d1spte=DBL1NULL;
  yds->d1svpr=DBL1NULL;
  yds->d1sprs=DBL1NULL;
  yds->d1ssir=DBL1NULL;
  yds->dsbuf=R0;
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
  /* element property database	*/
  ydpe->mprop=0; ydpe->nprop=0;
  ydpe->d1peca=DBL1NULL;
  ydpe->d1pecl=DBL1NULL;
  ydpe->d1peks=DBL1NULL;
  ydpe->d1pela=DBL1NULL;
  ydpe->d1pemu=DBL1NULL;
  ydpe->d1pepe=DBL1NULL;
  //ydpe->d1pept=DBL1NULL;
  //ydpe->d1pefr=DBL1NULL;
  ydpe->d1pera=DBL1NULL;
  ydpe->d1pero=DBL1NULL;
  ydpe->d1pevi=DBL1NULL;
  
  ydpe->mperow=0;
  ydpe->d2peint=DBL2NULL;
  
  ydpe->d1peem=DBL1NULL;
  ydpe->d1penu=DBL1NULL;

  ydpe->d1psem=DBL1NULL;

  ydpe->i1pecn=INT1NULL;
  ydpe->i1pecp=INT1NULL;
  ydpe->i1pect=INT1NULL;
  ydpe->i1pefr=INT1NULL;
  ydpe->i1pejp=INT1NULL;
  ydpe->i1pemb=INT1NULL;
  ydpe->i1pemn=INT1NULL;

  ydpe->i1pnib=INT1NULL;
  ydpe->i1psde=INT1NULL;
  ydpe->i1ptyp=INT1NULL;
  
  ydpe->i1pexc=INT1NULL;
  
  ydpe->i1usan=INT1NULL;
  ydpe->d1peex=DBL1NULL;
  ydpe->d1peey=DBL1NULL;
  ydpe->d1pemx=DBL1NULL;
  ydpe->d1pemy=DBL1NULL;
  ydpe->d1peg=DBL1NULL;
  
  ydpe->i1psup=INT1NULL;
    
  /* joint property database	*/
  ydpj->mpjset=0; ydpj->npjset=0;
  ydpj->d1pjfs=DBL1NULL;
  ydpj->d1pjft=DBL1NULL;
  ydpj->d1pjgf=DBL1NULL;
  ydpj->d1pjgs=DBL1NULL; 
  ydpj->d1pjco=DBL1NULL;
  ydpj->d1pjfr=DBL1NULL;
  ydpj->d1pjpe=DBL1NULL;
  ydpj->i1ptyp=INT1NULL;
  ydpj->i1psde=INT1NULL;
  ydpj->d1usaf=DBL1NULL;
  ydpj->d1pjcr=DBL1NULL;
  ydpj->d1pjfd=DBL1NULL;
  ydpj->d1pjtr=DBL1NULL;
  ydpj->d1pjgr=DBL1NULL;
  ydpj->d1pjsr=DBL1NULL;
  ydpj->d1pjal=DBL1NULL;
  ydpj->iusehy=0;
  
  /* mesh property database	*/
  ydpm->mpmcom=0;
  ydpm->mpmcol=0;
  ydpm->i2pmset=INT2NULL;
  ydpm->mpmrow=0;
  ydpm->i2pmij=INT2NULL;
  /* node property database	*/
  ydpn->mpnset=0; ydpn->npnset=0;
  ydpn->mpnfact=0, ydpn->npnfact=0;
  ydpn->d3pnfac=DBL3NULL;

  ydpn->i1pnfx=INT1NULL;
  ydpn->i1pnfy=INT1NULL;
  ydpn->i1pnfz=INT1NULL;

  ydpn->d1pnaf=DBL1NULL;
  ydpn->d1pnap=DBL1NULL;
  ydpn->d1pnat=DBL1NULL;
  ydpn->d1pnax=DBL1NULL;
  ydpn->d1pnay=DBL1NULL;
  ydpn->d1pnaz=DBL1NULL;

  ydpn->d1pnxx=DBL1NULL;
  ydpn->d1pnxy=DBL1NULL;
  ydpn->d1pnxz=DBL1NULL;
  ydpn->d1pnyx=DBL1NULL;
  ydpn->d1pnyy=DBL1NULL;
  ydpn->d1pnyz=DBL1NULL;
  ydpn->d1pnzx=DBL1NULL;
  ydpn->d1pnzy=DBL1NULL;
  ydpn->d1pnzz=DBL1NULL;
  
  /* Set discrete fracture network values to default */
  ydfn->iusefn=0;
  ydfn->mdfnfr=0;
  ydfn->mdfnno=0;
  ydfn->i2dfnn=INT2NULL;
  ydfn->d1dffr=DBL1NULL;
  ydfn->d1dfpe=DBL1NULL;
  ydfn->d1dfpt=DBL1NULL;
  ydfn->ddfnft=0;
  ydfn->ddfnco=0;
  ydfn->ddfngf=0;
  ydfn->ddfngs=0;
  ydfn->i1dfft=INT1NULL;
  
  /* set in-situ stress parameters to zero */
  ydis->iuseis=0;
  ydis->dcstxx=R0;
  ydis->dcstxy=R0;
  ydis->dcstyy=R0;
  ydis->dcsyxx=R0;
  ydis->dcsyxy=R0;
  ydis->dcsyyy=R0;
  ydis->dcsrfy=R0;
  
  /* set hydro-frac parameters to default */
  ydhf->iusehf=0;
  ydhf->ihftyp=1;               
  ydhf->dhfflp=R0;               
  ydhf->dhfflq=R0;
  ydhf->hfarow=0;
  ydhf->d2hfaf=DBL2NULL;
  ydhf->fluvol=R0;       
  ydhf->flupres=R0;
  ydhf->flurho0=R0;              
  ydhf->flupres0=R0;            
  ydhf->flubulk=R0;
  ydhf->fradim=2;
  ydhf->ihfmsin=0;
  ydhf->gravacc=R0;
  ydhf->d2wtlev=DBL2NULL;
  
  /* set seismic monitoring parameters to default */
  ydsm->iusesm=0;
  ydsm->dctwle=1;
  
  /* Set reference points to default Y-RC */    
  ydr->mnodim=0;  ydr->nnodim=0;  
  ydr->mrdim=0;  ydr->nrdim=0;
  ydr->mrldm=0;  ydr->nrldm=0;
  ydr->mrepo=0;  ydr->nrepo=0;
  ydr->nbrjointrb=0;
  ydr->d2rcig=DBL2NULL;
  ydr->d2rccg=DBL2NULL;
  ydr->d2rccl=DBL2NULL;
  ydr->d2rvcg=DBL2NULL;
  ydr->d2riLc=DBL2NULL;
  ydr->d2rsctr=DBL2NULL;
  ydr->i2relto=INT2NULL;
  ydr->i1rmyel=INT1NULL;
  ydr->i1rrpn=INT1NULL;
  ydr->i1rprop=INT1NULL;
  ydr->i1refbar=INT1NULL;
  ydr->i1myjoint=INT1NULL;
  ydr->d1rbsig=DBL1NULL;
  ydr->d1rbfrc=DBL1NULL;
  ydr->d1rbstr=DBL1NULL;
  ydr->d1rjsig=DBL1NULL;
  ydr->d1rjtau=DBL1NULL;
  ydr->d1rjslpnor=DBL1NULL;
  ydr->d1rjdelnor=DBL1NULL;
  ydr->i1rjstnor=INT1NULL;
  ydr->d2rjfrv=DBL2NULL;
  ydr->i2rbedn=INT2NULL;
  
  /* Set Steel Nodes/bar to default  Y-RC */   
  ydsb->msdim=0;  ydsb->nsdim=0;
  ydsb->msbar=0;  ydsb->nsbar=0;  ydsb->isfirst=0;
  ydsb->d2sic=DBL2NULL;
  ydsb->i1srpf=INT1NULL;
  ydsb->i1sbpr=INT1NULL;
  ydsb->d1spea=DBL1NULL;
  ydsb->d1sdiam=DBL1NULL;
  ydsb->d1crlcr=DBL1NULL;
  ydsb->d1smmdiam=DBL1NULL;
  ydsb->i1sbty=INT1NULL;
  ydsb->i1sbac=INT1NULL;
  
  /* Set Steel Property to default Y-RC */
  ydps->mprop=0;  ydps->nprop=0;
  ydps->d1young=DBL1NULL;
  ydps->d1sfc=DBL1NULL;
  ydps->d1mpsfc=DBL1NULL;
  ydps->d1sfy=DBL1NULL;
  ydps->d1epssh=DBL1NULL;
  ydps->d1sfu=DBL1NULL;
  ydps->d1epsu=DBL1NULL;
  ydps->d1sfbr=DBL1NULL;
  ydps->d1epsbr=DBL1NULL;
  ydps->d1stkn=DBL1NULL;
  ydps->d1stkt=DBL1NULL;
  ydps->d1styns=DBL1NULL;
  ydps->d1strns=DBL1NULL;
  ydps->d1stcoh=DBL1NULL;
  ydps->d1stfri=DBL1NULL;
}

static void Yrde(yde,finp,name)
  YDE    yde; FILE *finp; CHR *name;
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MELEM",5)==0)
  { INTr(finp,&(yde->melem));
  }
  else if(CHRcmp(namep,"NELEM",5)==0)
  { INTr(finp,&(yde->nelem));
    yde->nelemst=yde->nelem;
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
  { TformINT1(finp,-1,yde->melem,&(yde->i1elcf));
  }
  else if(CHRcmp(namep,"I1ELPR",6)==0)
  { TformINT1(finp,-1,yde->melem,&(yde->i1elpr));
    TformINT1(finp,-1,yde->melem,&(yde->i1elprtmp));
  }
  else if(CHRcmp(namep,"D2ELST",6)==0)
  { TformDBL2(finp,R0,yde->melst,yde->melem,&(yde->d2elst));
  }
  else if(CHRcmp(namep,"I2ELTO",6)==0)
  { TformINT2(finp,-1,yde->melno,yde->melem,&(yde->i2elto));
    TformINT2(finp,-1,yde->melno,yde->melem,&(yde->i2eltost));
  }
  else
  { CHRw(stderr,"Yrde: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
} }

static void Yrdi(ydi,finp,name)
  YDI   ydi; FILE *finp; CHR *name;
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MICOUP",6)==0)
  { INTr(finp,&(ydi->micoup));
  }
  else if(CHRcmp(namep,"NICOUP",6)==0)
  { INTr(finp,&(ydi->nicoup));
  }
  else if(CHRcmp(namep,"IIECFF",6)==0)
  { INTr(finp,&(ydi->iiecff));
  }
  else if(CHRcmp(namep,"DIEDI",5)==0)
  { DBLr(finp,&(ydi->diedi));
  }
  else if(CHRcmp(namep,"DIEZON",6)==0)
  { DBLr(finp,&(ydi->diezon));
  }
  else if(CHRcmp(namep,"D1IESL",6)==0)
  { TformDBL1(finp,R0,ydi->micoup,&(ydi->d1iesl));
  }
  else if(CHRcmp(namep,"I1IECN",6)==0)
  { TformINT1(finp,-1,ydi->micoup,&(ydi->i1iecn));
  }
  else if(CHRcmp(namep,"I1IECT",6)==0)
  { TformINT1(finp,-1,ydi->micoup,&(ydi->i1iect));
  }
  else if(CHRcmp(namep,"MISTATE",7)==0)
  {  INTr(finp,&(ydi->mistate));
  }
  else if(CHRcmp(namep,"D2SLDIS",7)==0)
  { TformDBL2(finp,R0,ydi->mistate, ydi->micoup,&(ydi->d2sldis));
  }
  else
  { CHRw(stderr,"Yrdi: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
} }

static void Yrdn(ydn,finp,name)
   YDN    ydn; FILE *finp; CHR *name;
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MNODIM",5)==0)
  { INTr(finp,&(ydn->mnodim));
  }
  else if(CHRcmp(namep,"NNODIM",5)==0)
  { INTr(finp,&(ydn->nnodim));
  }
  else if(CHRcmp(namep,"MNOPO",5)==0)
  { INTr(finp,&(ydn->mnopo));
  }
  else if(CHRcmp(namep,"NNOPO",5)==0)
  { INTr(finp,&(ydn->nnopo));
    ydn->nnopst=ydn->nnopo;
  }
  else if(CHRcmp(namep,"D1NMCT",6)==0)
  { TformDBL1(finp,R0,ydn->mnopo,&(ydn->d1nmct));
  }
  else if(CHRcmp(namep,"D2NCC",5)==0)
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2ncc));
  }
  else if(CHRcmp(namep,"D2NCI",5)==0)
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nci));
  }
  else if(CHRcmp(namep,"D2NFC",5)==0)
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nfc));
  }
  else if(CHRcmp(namep,"D2NVC",5)==0)
  { TformDBL2(finp,R0,ydn->mnodim,ydn->mnopo,&(ydn->d2nvc));
  }
  else if(CHRcmp(namep,"I1NOBF",6)==0)
  { TformINT1(finp,-1,ydn->mnopo,&(ydn->i1nobf));
  }
  else if(CHRcmp(namep,"I1NOPR",6)==0)
  { TformINT1(finp,-1,ydn->mnopo,&(ydn->i1nopr));
  }
  else if(CHRcmp(namep,"I1NOWE",6)==0)
  { TformINT1(finp,-1,ydn->mnopo,&(ydn->i1nowe));
  }
  else
  { CHRw(stderr,"Yrdn: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(105);
} }

static void Yrdb(ydb,finp,name)
   YDB    ydb; FILE *finp; CHR *name;
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MBORH",5)==0)
  {  INTr(finp,&(ydb->mborh));
  }
  else if(CHRcmp(namep,"NBORH",5)==0)
  {  INTr(finp,&(ydb->nborh));
  }
  else if(CHRcmp(namep,"MBDIM",5)==0)
  {  INTr(finp,&(ydb->mbdim));
  }
  else if(CHRcmp(namep,"NBDIM",5)==0)
  {  INTr(finp,&(ydb->nbdim));
  }
  else if(CHRcmp(namep,"NBPAF",5)==0)
  {  INTr(finp,&(ydb->nbpaf));
  }
  else if(CHRcmp(namep,"D2BCA",5)==0)
  { TformDBL2(finp,R0,ydb->mbdim,ydb->mborh,&(ydb->d2bca));
  }
  else if(CHRcmp(namep,"D2BCB",5)==0)
  { TformDBL2(finp,R0,ydb->mbdim,ydb->mborh,&(ydb->d2bcb));
  }
  else if(CHRcmp(namep,"D1BRAD",6)==0)
  { TformDBL1(finp,R0,ydb->mborh,&(ydb->d1brad));
  }
  else if(CHRcmp(namep,"D1BPAF",6)==0)
  { TformDBL1(finp,R0,ydb->nbpaf,&(ydb->d1bpaf));
  }
  else if(CHRcmp(namep,"D1BPTS",6)==0)
  { TformDBL1(finp,R0,ydb->mborh,&(ydb->d1bpts));
  }
  else if(CHRcmp(namep,"D1BPTE",6)==0)
  { TformDBL1(finp,R0,ydb->mborh,&(ydb->d1bpte));
  }
  else if(CHRcmp(namep,"D1BVDT",6)==0)
  { TformDBL1(finp,R0,ydb->mborh,&(ydb->d1bvdt));
  }
  else if(CHRcmp(namep,"D1BPRS",6)==0)
  { TformDBL1(finp,R0,ydb->mborh,&(ydb->d1bprs));
  }
  else if(CHRcmp(namep,"DBLMAX",6)==0)
  {  DBLr(finp,&(ydb->dblmax));
  }
  else if(CHRcmp(namep,"DBBUF",5)==0)
  {  DBLr(finp,&(ydb->dbbuf));
  }
  else
  { CHRw(stderr,"Yrdb: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(105);
} }

static void Yrds(yds,finp,name)
   YDS yds; FILE *finp; CHR *name;
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MSOUR",5)==0)
  {  INTr(finp,&(yds->msour));
  }
  else if(CHRcmp(namep,"NSOUR",5)==0)
  {  INTr(finp,&(yds->nsour));
  }
  else if(CHRcmp(namep,"MSDIM",5)==0)
  {  INTr(finp,&(yds->msdim));
  }
  else if(CHRcmp(namep,"NSDIM",5)==0)
  {  INTr(finp,&(yds->nsdim));
  }
  else if(CHRcmp(namep,"NSPAF",5)==0)
  {  INTr(finp,&(yds->nspaf));
  }
  else if(CHRcmp(namep,"NSSAF",5)==0)
  {  INTr(finp,&(yds->nssaf));
  }
  else if(CHRcmp(namep,"D2SCS",5)==0)
  { TformDBL2(finp,R0,yds->msdim,yds->msour,&(yds->d2scs));
  }
  else if(CHRcmp(namep,"D1SPAF",6)==0)
  { TformDBL1(finp,R0,yds->nspaf,&(yds->d1spaf));
  }
  else if(CHRcmp(namep,"D1SSAF",6)==0)
  { TformDBL1(finp,R0,yds->nssaf,&(yds->d1ssaf));
  }
  else if(CHRcmp(namep,"D1SPTS",6)==0)
  { TformDBL1(finp,R0,yds->msour,&(yds->d1spts));
  }
  else if(CHRcmp(namep,"D1SPTE",6)==0)
  { TformDBL1(finp,R0,yds->msour,&(yds->d1spte));
  }
  else if(CHRcmp(namep,"D1SVPR",6)==0)
  { TformDBL1(finp,R0,yds->msour,&(yds->d1svpr));
  }
  else if(CHRcmp(namep,"D1SPRS",6)==0)
  { TformDBL1(finp,R0,yds->msour,&(yds->d1sprs));
  }
  else if(CHRcmp(namep,"D1SSIR",6)==0)
  { TformDBL1(finp,R0,yds->msour,&(yds->d1ssir));
  }
  else if(CHRcmp(namep,"DSBUF",5)==0)
  { DBLr(finp,&(yds->dsbuf));
  }
  else
  { CHRw(stderr,"Yrds: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(105);
} }

static void Yrdo(ydo,finp,name)
   YDO    ydo; FILE *finp; CHR *name;
{ INT i; CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MOHYS",5)==0)
  { INTr(finp,&(ydo->mohys))
    i=(ydo->mohys)*sizeof(FILE*);
    if(i>0)ydo->f2ohyf=(FILE**)MALLOC(i);
    for(i=0;i<(ydo->mohys);i++)
    { ydo->f2ohyf[i]=FILENULL;
  } }
  else if(CHRcmp(namep,"NOHYS",5)==0)
  {  INTr(finp,&(ydo->nohys));
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



static void Yrdpe(ydpe,finp,name)
  YDPE    ydpe; FILE *finp; CHR *name;
{ CHR *namep;

  namep=name+9;
  if(CHRcmp(namep,"MPROP",5)==0)
  { INTr(finp,&(ydpe->mprop));
  }
  else if(CHRcmp(namep,"NPROP",5)==0)
  { INTr(finp,&(ydpe->nprop));
  }
  else if(CHRcmp(namep,"D1PECA",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1peca));
  }
  else if(CHRcmp(namep,"D1PECL",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pecl));
  }
  else if(CHRcmp(namep,"D1PEKS",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1peks));
  }
  else if(CHRcmp(namep,"D1PELA",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pela));
  }
  else if(CHRcmp(namep,"D1PEMU",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pemu));
  }
  else if(CHRcmp(namep,"D1PEEM",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1peem));
  }
  else if(CHRcmp(namep,"D1PENU",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1penu));
  }
  else if(CHRcmp(namep,"D1PEPE",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pepe));
  }
  /*else if(CHRcmp(namep,"D1PEPT",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pept));
  }
  else if(CHRcmp(namep,"D1PEFR",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pefr));
  }*/
  else if(CHRcmp(namep,"D1PERA",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pera));
  }
  else if(CHRcmp(namep,"D1PERO",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pero));
  }
  else if(CHRcmp(namep,"D1PEVI",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pevi));
  }
  else if(CHRcmp(namep,"D1PSEM",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1psem));
  }
  else if(CHRcmp(namep,"I1PECN",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1pecn));
  }
  else if(CHRcmp(namep,"I1PECP",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1pecp));
  }
  else if(CHRcmp(namep,"I1PECT",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1pect));
  }
  else if(CHRcmp(namep,"I1PEFR",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1pefr));
  }
  else if(CHRcmp(namep,"I1PEJP",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1pejp));
  }
  else if(CHRcmp(namep,"I1PEMB",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1pemb));
  }
  else if(CHRcmp(namep,"I1PEMN",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1pemn));
  }
  else if(CHRcmp(namep,"I1PNIB",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1pnib));
  }
  else if(CHRcmp(namep,"I1PSDE",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1psde));
  }
  else if(CHRcmp(namep,"I1PTYP",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1ptyp));
  }
  else if(CHRcmp(namep,"I1PEXC",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1pexc));
  }
  else if(CHRcmp(namep,"I1USAN",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1usan));
  }
  else if(CHRcmp(namep,"D1PEEX",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1peex));
  }
  else if(CHRcmp(namep,"D1PEEY",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1peey));
  }
  else if(CHRcmp(namep,"D1PEMX",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pemx));
  }
  else if(CHRcmp(namep,"D1PEMY",6)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1pemy));
  }
  else if(CHRcmp(namep,"D1PEG",5)==0)
  { TformDBL1(finp,R0,ydpe->mprop,&(ydpe->d1peg));
  }
  else if(CHRcmp(namep,"I1PSUP",6)==0)
  { TformINT1(finp,-1,ydpe->mprop,&(ydpe->i1psup));
  }
  else if(CHRcmp(namep,"MPEROW",6)==0)
  { INTr(finp,&(ydpe->mperow));
  }
  else if(CHRcmp(namep,"D2PEINT",7)==0)
  { TformDBL2(finp,0.0,5,ydpe->mperow,&(ydpe->d2peint));
  } 
  else
  { CHRw(stderr,"Yrdpe: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
} }

static void Yrdpj(ydpj,finp,name)
YDPJ ydpj; FILE *finp; CHR *name;
{ CHR *namep;
  namep=name+9;

  if(CHRcmp(namep,"MPJSET",6)==0)
  { INTr(finp,&(ydpj->mpjset));
  }
  else if(CHRcmp(namep,"NPJSET",6)==0)
  { INTr(finp,&(ydpj->npjset));
  }
  else if(CHRcmp(namep,"D1PJFS",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjfs));
  }
  else if(CHRcmp(namep,"D1PJFT",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjft));
  }
  else if(CHRcmp(namep,"D1PJGF",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjgf));
  }
  else if(CHRcmp(namep,"D1PJGS",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjgs));
  }
  else if(CHRcmp(namep,"D1PJCO",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjco));
  }
  else if(CHRcmp(namep,"D1PJFR",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjfr));
  }
  else if(CHRcmp(namep,"D1PJPE",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjpe));
  }
  else if(CHRcmp(namep,"I1PSDE",6)==0)
  { TformINT1(finp,-1,ydpj->mpjset,&(ydpj->i1psde));
  }
  else if(CHRcmp(namep,"I1PTYP",6)==0)
  { TformINT1(finp,-1,ydpj->mpjset,&(ydpj->i1ptyp));
  }
  else if(CHRcmp(namep,"D1USAF",6)==0)
  { TformDBL1(finp,-1,ydpj->mpjset,&(ydpj->d1usaf));
  } 
  else if(CHRcmp(namep,"D1PJAL",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjal));
  }
  else if(CHRcmp(namep,"D1PJCR",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjcr));
  }
  else if(CHRcmp(namep,"D1PJFD",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjfd));
  }
  else if(CHRcmp(namep,"D1PJTR",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjtr));
  }
  else if(CHRcmp(namep,"D1PJGR",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjgr));
  }
  else if(CHRcmp(namep,"D1PJSR",6)==0)
  { TformDBL1(finp,R0,ydpj->mpjset,&(ydpj->d1pjsr));
  }
  else if(CHRcmp(namep,"IUSEHY",6)==0)
  { INTr(finp,&(ydpj->iusehy));
  }
  else
  { CHRw(stderr,"Yrdpj: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
  }
}
 
static void Yrdpm(ydpm,finp,name)
YDPM ydpm; FILE *finp; CHR *name;
{ CHR *namep;
  namep=name+9;

  if(CHRcmp(namep,"MPMCOM",6)==0)
  { INTr(finp,&(ydpm->mpmcom));
  }
  else if(CHRcmp(namep,"MPMCOL",6)==0)
  { INTr(finp,&(ydpm->mpmcol));
  }
  else if(CHRcmp(namep,"I2PMSET",7)==0)
  { TformINT2(finp,-1,ydpm->mpmcol,ydpm->mpmcom,&(ydpm->i2pmset));
  }
  else if(CHRcmp(namep,"MPMROW",6)==0)
  { INTr(finp,&(ydpm->mpmrow));
  }
  else if(CHRcmp(namep,"I2PMIJ",6)==0)
  { TformINT2(finp,-1,3,ydpm->mpmrow,&(ydpm->i2pmij));
  }
  else
  { CHRw(stderr,"Yrdpm: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
  }
}

static void Yrdpn(ydpn,finp,name)
YDPN ydpn; FILE *finp; CHR *name;
{ CHR *namep;

  namep=name+9;
  if(CHRcmp(namep,"MPNSET",6)==0)
  { INTr(finp,&(ydpn->mpnset));
  }
  else if(CHRcmp(namep,"NPNSET",6)==0)
  { INTr(finp,&(ydpn->npnset));
  }
  else if(CHRcmp(namep,"MPNFACT",7)==0)
  { INTr(finp,&(ydpn->mpnfact));
  }
  else if(CHRcmp(namep,"NPNFACT",7)==0)
  { INTr(finp,&(ydpn->npnfact));
  }
  else if(CHRcmp(namep,"D3PNFAC",7)==0)
  { TformDBL3(finp,-R1,2,ydpn->mpnset,ydpn->mpnfact,&(ydpn->d3pnfac));
  }
  else if(CHRcmp(namep,"I1PNFX",6)==0)
  { TformINT1(finp,-1,ydpn->mpnset,&(ydpn->i1pnfx));
  }
  else if(CHRcmp(namep,"I1PNFY",6)==0)
  { TformINT1(finp,-1,ydpn->mpnset,&(ydpn->i1pnfy));
  }
  else if(CHRcmp(namep,"I1PNFZ",6)==0)
  { TformINT1(finp,-1,ydpn->mpnset,&(ydpn->i1pnfz));
  }
  else if(CHRcmp(namep,"D1PNAF",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnaf));
  }
  else if(CHRcmp(namep,"D1PNAT",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnat));
  }
  else if(CHRcmp(namep,"D1PNAP",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnap));
  }
  else if(CHRcmp(namep,"D1PNAX",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnax));
  }
  else if(CHRcmp(namep,"D1PNAY",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnay));
  }
  else if(CHRcmp(namep,"D1PNAZ",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnaz));
  }
  else if(CHRcmp(namep,"D1PNXX",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnxx));
  }
  else if(CHRcmp(namep,"D1PNXY",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnxy));
  }
  else if(CHRcmp(namep,"D1PNXZ",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnxz));
  }
  else if(CHRcmp(namep,"D1PNYX",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnyx));
  }
  else if(CHRcmp(namep,"D1PNYY",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnyy));
  }
  else if(CHRcmp(namep,"D1PNYZ",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnyz));
  }
  else if(CHRcmp(namep,"D1PNZX",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnzx));
  }
  else if(CHRcmp(namep,"D1PNZY",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnzy));
  }
  else if(CHRcmp(namep,"D1PNZZ",6)==0)
  { TformDBL1(finp,R0,ydpn->mpnset,&(ydpn->d1pnzz));
  }
  else
  { CHRw(stderr,"Yrdpn: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
  }
}

static void Yrdfn(ydfn,finp,name)
YDFN ydfn; FILE *finp; CHR *name;
{ CHR *namep;
  namep=name+9;

  if(CHRcmp(namep,"IUSEFN",6)==0)
  { INTr(finp,&(ydfn->iusefn));
  }
  else if(CHRcmp(namep,"MDFNFR",6)==0)
  { INTr(finp,&(ydfn->mdfnfr));
  }
  else if(CHRcmp(namep,"MDFNNO",6)==0)
  { INTr(finp,&(ydfn->mdfnno));
  }
  else if(CHRcmp(namep,"I2DFNN",6)==0)
  { TformINT2(finp,-1,ydfn->mdfnno,ydfn->mdfnfr,&(ydfn->i2dfnn));
  }
  else if(CHRcmp(namep,"D1DFFR",6)==0)
  { TformDBL1(finp,R0,ydfn->mdfnfr,&(ydfn->d1dffr));
  }
  else if(CHRcmp(namep,"D1DFPE",6)==0)
  { TformDBL1(finp,R0,ydfn->mdfnfr,&(ydfn->d1dfpe));
  }
  else if(CHRcmp(namep,"D1DFPT",6)==0)
  { TformDBL1(finp,R0,ydfn->mdfnfr,&(ydfn->d1dfpt));
  }
  else if(CHRcmp(namep,"DDFNFT",6)==0)
  { DBLr(finp,&(ydfn->ddfnft));
  }
  else if(CHRcmp(namep,"DDFNCO",6)==0)
  { DBLr(finp,&(ydfn->ddfnco));
  }
  else if(CHRcmp(namep,"DDFNGF",6)==0)
  { DBLr(finp,&(ydfn->ddfngf));
  }
  else if(CHRcmp(namep,"DDFNGS",6)==0)
  { DBLr(finp,&(ydfn->ddfngs));
  }
  else if(CHRcmp(namep,"I1DFFT",6)==0)
  { TformINT1(finp,-1,ydfn->mdfnfr,&(ydfn->i1dfft));
  }
  else
  { CHRw(stderr,"Yrdpm: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
  }
}
/* read parameters for the in-situ stress */
static void Yrdis(ydis,finp,name)
YDIS ydis; FILE *finp; CHR *name;
{ CHR *namep;
  namep=name+9;
  if(CHRcmp(namep,"IUSEIS",5)==0)
  { INTr(finp,&(ydis->iuseis));
  }
  else if(CHRcmp(namep,"DCSTXX",6)==0)
  { DBLr(finp,&(ydis->dcstxx));
  }
  else if(CHRcmp(namep,"DCSTXY",6)==0)
  { DBLr(finp,&(ydis->dcstxy));
  }
  else if(CHRcmp(namep,"DCSTYY",6)==0)
  { DBLr(finp,&(ydis->dcstyy));
  }
  else if(CHRcmp(namep,"DCSYXX",6)==0)
  { DBLr(finp,&(ydis->dcsyxx));
  }
  else if(CHRcmp(namep,"DCSYXY",6)==0)
  { DBLr(finp,&(ydis->dcsyxy));
  }
  else if(CHRcmp(namep,"DCSYYY",6)==0)
  { DBLr(finp,&(ydis->dcsyyy));
  }
  else if(CHRcmp(namep,"DCSRFY",6)==0)
  { DBLr(finp,&(ydis->dcsrfy));
  }
  else
  { CHRw(stderr,"Yrdis: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
} }

/* read parameters for hydro-frac */
static void Yrdhf(ydhf,finp,name)
YDHF ydhf; FILE *finp; CHR *name;
{ CHR *namep;
  namep=name+9;
  if(CHRcmp(namep,"IUSEHF",6)==0)
  { INTr(finp,&(ydhf->iusehf));
  }
  else if(CHRcmp(namep,"IHFTYP",6)==0)
  { INTr(finp,&(ydhf->ihftyp));
  }
  else if(CHRcmp(namep,"DHFFLP",6)==0)
  { DBLr(finp,&(ydhf->dhfflp));
  }
  else if(CHRcmp(namep,"DHFFLQ",6)==0)
  { DBLr(finp,&(ydhf->dhfflq));
  }
  else if(CHRcmp(namep,"HFAROW",6)==0)
  { INTr(finp,&(ydhf->hfarow));
  }
  else if(CHRcmp(namep,"D2HFAF",6)==0)
  { TformDBL2(finp,-R1,2,ydhf->hfarow,&(ydhf->d2hfaf));
  }
  else if(CHRcmp(namep,"FLUPRES",7)==0)
  { DBLr(finp,&(ydhf->flupres));
  }
  else if(CHRcmp(namep,"FLURHO0",7)==0)
  { DBLr(finp,&(ydhf->flurho0));
  }
  else if(CHRcmp(namep,"FLUPRES0",8)==0)
  { DBLr(finp,&(ydhf->flupres0));
  }
  else if(CHRcmp(namep,"FLUBULK",7)==0)
  { DBLr(finp,&(ydhf->flubulk));
  }
  else if(CHRcmp(namep,"FRADIM",6)==0)
  { INTr(finp,&(ydhf->fradim));
  }
  else if(CHRcmp(namep,"GRAVACC",7)==0)
  { DBLr(finp,&(ydhf->gravacc));
  }
  else if(CHRcmp(namep,"D2WTLEV",7)==0)
  { TformDBL2(finp,-R1,2,2,&(ydhf->d2wtlev));
  }
  else
  { CHRw(stderr,"Yrdhf: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
} }

/* read parameters for seismic monitoring */
static void Yrdsm(ydsm,finp,name)
YDSM ydsm; FILE *finp; CHR *name;
{ CHR *namep;
  namep=name+9;
  if(CHRcmp(namep,"IUSESM",6)==0)
  { INTr(finp,&(ydsm->iusesm));
  }
  else if(CHRcmp(namep,"DCTWLE",6)==0)
  { DBLr(finp,&(ydsm->dctwle));
  }
  else
  { CHRw(stderr,"Yrdsm: unknown name: ");
    CHRw(stderr,name);
    CHRwcr(stderr);
    exit(104);
} }


/* read parameters for reference points Y-RC */
static void Yrdr(ydr,finp,name)                   
   YDR    ydr; FILE *finp; CHR *name; 
{ CHR *namep;

  namep=name+8;
  if(CHRcmp(namep,"MNODIM",6)==0) 
  {  INTr(finp,&(ydr->mnodim));
  }
  else if(CHRcmp(namep,"NNODIM",6)==0) 
  {  INTr(finp,&(ydr->nnodim));
  }
  else if(CHRcmp(namep,"MRDIM",5)==0) 
  {  INTr(finp,&(ydr->mrdim));
  }
  else if(CHRcmp(namep,"NRDIM",5)==0) 
  {  INTr(finp,&(ydr->nrdim));
  }
  else if(CHRcmp(namep,"MRLDM",5)==0) 
  {  INTr(finp,&(ydr->mrldm));
  }
  else if(CHRcmp(namep,"NRLDM",5)==0) 
  {  INTr(finp,&(ydr->nrldm));
  }
  else if(CHRcmp(namep,"MREPO",5)==0) 
  {  INTr(finp,&(ydr->mrepo));
  }
  else if(CHRcmp(namep,"NREPO",5)==0) 
  {  INTr(finp,&(ydr->nrepo));
  }
  else if(CHRcmp(namep,"NBRJOINTRB",10)==0) 
  {  INTr(finp,&(ydr->nbrjointrb));
  }
  else if(CHRcmp(namep,"D2RCIG",6)==0) 
  { TformDBL2(finp,R0,ydr->mnodim,ydr->mrepo,&(ydr->d2rcig)); 
  }
  else if(CHRcmp(namep,"D2RCCG",6)==0) 
  { TformDBL2(finp,R0,ydr->mnodim,ydr->mrepo,&(ydr->d2rccg)); 
  }
  else if(CHRcmp(namep,"D2RCCL",6)==0) 
  { TformDBL2(finp,R0,ydr->mrldm,ydr->mrepo,&(ydr->d2rccl)); 
  }
  else if(CHRcmp(namep,"D2RVCG",6)==0) 
  { TformDBL2(finp,R0,ydr->mnodim,ydr->mrepo,&(ydr->d2rvcg)); 
  }
  else if(CHRcmp(namep,"D2RILC",6)==0) 
  { TformDBL2(finp,R0,ydr->mrldm,ydr->mrepo,&(ydr->d2riLc)); 
  }
  else if(CHRcmp(namep,"D2RSCTR",7)==0) 
  { TformDBL2(finp,R0,2,ydr->mrepo,&(ydr->d2rsctr)); 
  }
  else if(CHRcmp(namep,"I2RELTO",7)==0) 
  { TformINT2(finp,-1,2,ydr->mrepo,&(ydr->i2relto)); 
  }
  else if(CHRcmp(namep,"I1RMYEL",7)==0) 
  { TformINT1(finp,-1,ydr->mrepo,&(ydr->i1rmyel));
  }
  else if(CHRcmp(namep,"I1RRPN",6)==0) 
  { TformINT1(finp,-1,ydr->mrepo,&(ydr->i1rrpn));
  }
  else if(CHRcmp(namep,"I1RPROP",7)==0) 
  { TformINT1(finp,-1,ydr->mrepo,&(ydr->i1rprop));
  }
  else if(CHRcmp(namep,"I1REFBAR",8)==0) 
  { TformINT1(finp,-1,ydr->mrepo,&(ydr->i1refbar));
  }
  else if(CHRcmp(namep,"I1MYJOINT",9)==0) 
  { TformINT1(finp,-1,ydr->mrepo,&(ydr->i1myjoint));
  }
  else
  { CHRw(stderr,"Yrdr: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
} }

/* read parameters for steel elements Y-RC */
static void Yrdsb(ydsb,finp,name)                   
   YDSB    ydsb; FILE *finp; CHR *name; 
{ CHR *namep;

  namep=name+9;
  if(CHRcmp(namep,"MSDIM",5)==0) 
  {  INTr(finp,&(ydsb->msdim));
  }
  else if(CHRcmp(namep,"NSDIM",5)==0) 
  {  INTr(finp,&(ydsb->nsdim));
  }
  else if(CHRcmp(namep,"MSBAR",5)==0) 
  {  INTr(finp,&(ydsb->msbar));
  }
  else if(CHRcmp(namep,"NSBAR",5)==0) 
  {  INTr(finp,&(ydsb->nsbar));
  }
  else if(CHRcmp(namep,"ISFIRST",7)==0) 
  {  INTr(finp,&(ydsb->isfirst));
  }
  else if(CHRcmp(namep,"D2SIC",5)==0) 
  { TformDBL2(finp,R0,ydsb->msdim,ydsb->msbar,&(ydsb->d2sic)); 
  }
  else if(CHRcmp(namep,"I1SRPF",6)==0) 
  { TformINT1(finp,-1,ydsb->nsbar,&(ydsb->i1srpf));
  }
  else if(CHRcmp(namep,"I1SBPR",6)==0) 
  { TformINT1(finp,-1,ydsb->msbar,&(ydsb->i1sbpr));
  }
  else if(CHRcmp(namep,"D1SPEA",6)==0) 
  { TformDBL1(finp,-1,ydsb->msbar,&(ydsb->d1spea)); 
  }
  else if(CHRcmp(namep,"D1SMMDIAM",9)==0) 
  { TformDBL1(finp,-1,ydsb->msbar,&(ydsb->d1smmdiam)); 
  }
  else if(CHRcmp(namep,"D1SDIAM",7)==0) 
  { TformDBL1(finp,-1,ydsb->msbar,&(ydsb->d1sdiam)); 
  }
  else if(CHRcmp(namep,"D1CRLCR",7)==0) 
  { TformDBL1(finp,-1,ydsb->msbar,&(ydsb->d1crlcr)); 
  }
  else if(CHRcmp(namep,"I1SBTY",6)==0) 
  { TformINT1(finp,-1,ydsb->msbar,&(ydsb->i1sbty));
  }
  else if(CHRcmp(namep,"I1SBAC",6)==0) 
  { TformINT1(finp,-1,ydsb->msbar,&(ydsb->i1sbac));
  }
  else
  { CHRw(stderr,"Yrdsb: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
} }

/* read steel properties Y-RC */
static void Yrdps(ydps,finp,name)                    
   YDPS    ydps; FILE *finp; CHR *name; 
{ CHR *namep;

  namep=name+9;
  if(CHRcmp(namep,"MPROP",5)==0) 
  {  INTr(finp,&(ydps->mprop));
  }
  else if(CHRcmp(namep,"NPROP",5)==0) 
  {  INTr(finp,&(ydps->nprop));
  }
  else if(CHRcmp(namep,"D1YOUNG",7)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1young)); 
  }
  else if(CHRcmp(namep,"D1MPSFC",7)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1mpsfc)); 
  }
  else if(CHRcmp(namep,"D1SFC",5)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1sfc)); 
  }
  else if(CHRcmp(namep,"D1SFY",5)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1sfy)); 
  }
  else if(CHRcmp(namep,"D1EPSSH",7)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1epssh)); 
  }
  else if(CHRcmp(namep,"D1SFU",5)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1sfu)); 
  }
  else if(CHRcmp(namep,"D1EPSU",6)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1epsu)); 
  }
  else if(CHRcmp(namep,"D1SFBR",6)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1sfbr)); 
  }
  else if(CHRcmp(namep,"D1EPSBR",7)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1epsbr)); 
  }
  else if(CHRcmp(namep,"D1STKN",6)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1stkn)); 
  }
  else if(CHRcmp(namep,"D1STKT",6)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1stkt)); 
  }
  else if(CHRcmp(namep,"D1STYNS",7)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1styns)); 
  }
  else if(CHRcmp(namep,"D1STRNS",7)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1strns)); 
  }
  else if(CHRcmp(namep,"D1STCOH",7)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1stcoh)); 
  }
  else if(CHRcmp(namep,"D1STFRI",7)==0) 
  { TformDBL1(finp,-1,ydps->mprop,&(ydps->d1stfri)); 
  }
  else
  { CHRw(stderr,"Yrdps: unknown name: ");
    CHRw(stderr,name); 
    CHRwcr(stderr);
    exit(105);
} }



/***************************PUBLIC***********************************/
INT Yrd(namep,yd)
   CHR *namep; YD yd;
{ INT icount;
  CHR name[300];
  YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDB ydb=&(yd->ydb);
  YDS yds=&(yd->yds);
  YDO ydo=&(yd->ydo);
  YDPE ydpe=&(yd->ydpe);
  YDPN ydpn=&(yd->ydpn);
  YDPJ ydpj=&(yd->ydpj);
  YDPM ydpm=&(yd->ydpm);
  YDFN ydfn=&(yd->ydfn);
  YDIS ydis=&(yd->ydis);
  YDHF ydhf=&(yd->ydhf);
  YDSM ydsm=&(yd->ydsm);
  YDR ydr=&(yd->ydr);         /* Y-RC */
  YDSB ydsb=&(yd->ydsb);      /* Y-RC */
  YDPS ydps=&(yd->ydps);      /* Y-RC */
  
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
        }
      }while((FILEND(ydc->finp)==0)&&(CHRcmp(name,"*/",2)!=0));
    }
    else if(CHRcmp(name, "/YD/",4)==0)
    { CHRw(ydc->fcheck,name); CHRwcr(ydc->fcheck);
      if(CHRcmp(name,"/YD/YDC/",8)==0)        /* read control data                  */
      { Yrdc(ydc,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDE/",8)==0)   /* read data elements                 */
      { Yrde(yde,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDI/",8)==0)   /* read data interaction              */
      { Yrdi(ydi,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDN/",8)==0)   /* read data nodes                    */
      { Yrdn(ydn,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDB/",8)==0)   /* read data boreholes                */
      { Yrdb(ydb,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDS/",8)==0)   /* read data sources                  */
      { Yrds(yds,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDO/",8)==0)   /* read data output                   */
      { Yrdo(ydo,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDPE/",9)==0)  /* read data properties elements      */
      { Yrdpe(ydpe,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDPN/",9)==0)  /* read data properties nodes         */
      { Yrdpn(ydpn,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDPJ/",9)==0)  /* read data properties joints        */
      { Yrdpj(ydpj,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDPM/",9)==0)  /* read data properties meshing       */
      { Yrdpm(ydpm,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDFN/",9)==0)  /* read data properties DFN           */
      { Yrdfn(ydfn,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDIS/",9)==0)  /* read in-situ stress parameters     */
      { Yrdis(ydis,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDHF/",9)==0)  /* read hydro-frac parameters         */
      { Yrdhf(ydhf,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDSM/",9)==0)  /* read seismic monitoring parameters */
      { Yrdsm(ydsm,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDR/",8)==0)  /* read referent points  Y-RC           */
      { Yrdr(ydr,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDSB/",9)==0)  /* read steel element  Y-RC             */
      { Yrdsb(ydsb,ydc->finp,name);
      }
      else if(CHRcmp(name,"/YD/YDPS/",9)==0)  /* read steel property  Y-RC           */
      { Yrdps(ydps,ydc->finp,name);
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
    ydo->foebrk=FILENULL;
    ydo->foesft=FILENULL;
    ydo->fohyfr=FILENULL;
    ydo->fofrac=FILENULL;
    CHRr(ydc->finp,name);
  }

  fclose(ydc->finp);
  fclose(ydc->fcheck);
  return 0;
}

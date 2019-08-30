/**********************************************************************/
/** Copyright (C) 2008,                                              **/
/** Queen Mary University of London (QMUL) & Imperial College        **/
/** of Science, Technology and Medicine (ICSTM). All rights reserved.**/
/** Implemented for you by Dr Jiansheng Xiang                        **
 
* This code is part of the Virtual Geoscience Workbench (VGW) developed
* jointly by ICSTM and QMUL through two related parallel projects at 
* ICSTM and QMUL respectively funded by EPSRC. 
*
* This code is provided by copyright holders under the GNU Lesser 
* General Public License (LGPL). It is open source code; you can 
* redistribute it and/or modify it under the terms of the GNU Lesser 
* General Public License version 3.  
*  
* This code is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty 
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
* the GNU Lesser General Public License for more details,
* http://www.gnu.org/licenses/lgpl.txt. 
*  
* You should have received a copy of the GNU Lesser General Public 
* License along with this code; if not, write to 
* Dr Jiansheng Xiang Prof Antonio Munjiza or Dr John-Paul Latham 
* j.xiang@imperial.ac.uk a.munjiza@qmul.ac.uk 
* or j.p.latham@imperial.ac.uk 
* ******************************************************************* */ 
#include "Yd.h"

/****************** Control ***************************************/
static void Ywdc(ydc,fout)
  YDC ydc; FILE *fout;
{ CHRw(fout,"     /*   Control   */");
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/MCSTEP");
  CHRwsp(fout);
  INTw(fout,ydc->mcstep,8);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/NCSTEP");
  CHRwsp(fout);
  INTw(fout,ydc->ncstep+1,8);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCGRAY");
  CHRwsp(fout);
  DBLw(fout,ydc->dcgray,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSIZC");
  CHRwsp(fout);
  DBLw(fout,ydc->dcsizc,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSIZF");
  CHRwsp(fout);
  DBLw(fout,ydc->dcsizf,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSIZS");
  CHRwsp(fout);
  DBLw(fout,ydc->dcsizs,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSIZV");
  CHRwsp(fout);
  DBLw(fout,ydc->dcsizv,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSIZD");
  CHRwsp(fout);
  DBLw(fout,ydc->dcsizd,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSIZA");
  CHRwsp(fout);
  DBLw(fout,ydc->dcsiza,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCSTEC");
  CHRwsp(fout);
  DBLw(fout,ydc->dcstec,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCTIME");
  CHRwsp(fout);
  DBLw(fout,ydc->dctime,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCGRST");
  CHRwsp(fout);
  DBLw(fout,ydc->dcgrst,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/DCRMPT");
  CHRwsp(fout);
  DBLw(fout,ydc->dcrmpt,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/ICOUTF");
  CHRwsp(fout);
  INTw(fout,ydc->icoutf,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/ICOUTI");
  CHRwsp(fout);
  INTw(fout,ydc->icouti,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/ICSAVF");
  CHRwsp(fout);
  INTw(fout,ydc->icsavf,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/ICOUTP");
  CHRwsp(fout);
  INTw(fout,ydc->icoutp,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/ICFMTY");
  CHRwsp(fout);
  INTw(fout,ydc->icfmty,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDC/ICIATY");
  CHRwsp(fout);
  INTw(fout,ydc->iciaty,5);
  CHRwcr(fout);
}

/****************** Elements ***************************************/
static void Ywde(yde,fout)
  YDE yde; FILE *fout;
{ INT icount;
  INT jcount;
  
  CHRw(fout,"     /*   Elements   */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDE/MELEM");
  CHRwsp(fout);
  INTw(fout,yde->melem,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/NELEM");
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/MELST");
  CHRwsp(fout);
  INTw(fout,yde->melst,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/NELST");
  CHRwsp(fout);
  INTw(fout,yde->nelst,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/MELNO");
  CHRwsp(fout);
  INTw(fout,yde->melno,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/NELNO");
  CHRwsp(fout);
  INTw(fout,yde->nelno,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/D2ELST");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/I1ELCF");
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);
  for(icount=0;icount<yde->nelem;icount++)
  { INTwn(fout,-1,4);
	CHRwsp(fout);
  }
  CHRwcr(fout);
  
  CHRw(fout,"/YD/YDE/I1ELTY");
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);
  for(icount=0;icount<yde->nelem;icount++)
  { INTwn(fout,yde->i1elty[icount],4);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/I1ELPR");
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);
  for(icount=0;icount<yde->nelem;icount++)
  { INTw(fout,yde->i1elpr[icount],4);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/I2ELTO");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,yde->melno,5);
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);
  for(icount=0;icount<yde->nelem;icount++)
  { for(jcount=0;jcount<yde->melno;jcount++)
    { INTwn(fout,yde->i2elto[jcount][icount],6);
	  CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/D2ELFS");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,yde->melno,5);
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);
  for(icount=0;icount<yde->nelem;icount++)
  { for(jcount=0;jcount<yde->melno;jcount++)
    { DBLwa(fout,yde->d2elfs[jcount][icount]);
	  CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDE/I2ELJP");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,yde->nelno,5);
  CHRwsp(fout);
  INTw(fout,yde->nelem,5);
  CHRwcr(fout);
  for(icount=0;icount<yde->nelem;icount++)
  { for(jcount=0;jcount<yde->nelno;jcount++)
    { INTwn(fout,yde->i2eljp[jcount][icount],6);
	  CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout);
}

/****************** Joints *************************************/
static void Ywdj(ydj,fout)
  YDJ ydj; FILE *fout;
{ INT icount;
  
  CHRw(fout,"     /*   Joints   */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDJ/NJOINT");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/I1JTID");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { INTwn(fout,ydj->i1jtid[icount],4);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JKNI");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jkni[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JKNC");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jknc[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JKSC");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jksc[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JNST");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jnst[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JSST");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jsst[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JAPI");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1japi[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JAPC");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1japc[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JAPH");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1japh[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JAPR");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1japr[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JSDC");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jsdc[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JDLC");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jdlc[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JSDP");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jsdp[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JEFL");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jefl[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JJRC");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jjrc[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JJCS");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jjcs[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JPHI");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jphi[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JFMD");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jfmd[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JFPR");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jfpr[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDJ/D1JFET");
  CHRwsp(fout);
  INTw(fout,ydj->njoint,5);
  CHRwcr(fout);
  for(icount=0;icount<ydj->njoint;icount++)
  { DBLwa(fout,ydj->d1jfet[icount]);
	CHRwsp(fout);
  }
  CHRwcr(fout);
}

/****************** Interactions *************************************/
static void Ywdi(ydi,fout)
  YDI ydi; FILE *fout;
{ CHRw(fout,"     /*   Interactions     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDI/MICOUP");
  CHRwsp(fout);
  INTw(fout,ydi->micoup,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/NICOUP");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/IIECFF");
  CHRwsp(fout);
  INTwn(fout,-2,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/DIEDI");
  CHRwsp(fout);
  DBLw(fout,200.0,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/DIEZON");
  CHRwsp(fout);
  DBLw(fout,ydi->diezon,18);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/D1IESL");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/I1IECN");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDI/I1IECT");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);
}

/****************** Nodes ******************************************/
static void Ywdn(ydn,fout)
  YDN ydn; FILE *fout;
{ INT icount,jcount;
  
  CHRw(fout,"     /*   Nodes     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDN/MNODIM");
  CHRwsp(fout);
  INTw(fout,ydn->mnodim,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/NNODIM");
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/MNOPO");
  CHRwsp(fout);
  INTw(fout,ydn->mnopo,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/NNOPO");
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D2NCC");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLwa(fout,ydn->d2ncc[jcount][icount]);
      CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D2NCI");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLwa(fout,ydn->d2nci[jcount][icount]);
      CHRwsp(fout);    
    }
    CHRwcr(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D2NFC");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLwa(fout,ydn->d2nfc[jcount][icount]);
      CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D1NMCT");
  CHRwsp(fout);
  INTw(fout,0,5);
  CHRwcr(fout);
  for(icount=0;icount<0;icount++)
  { DBLwa(fout,ydn->d1nmct[icount]);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/D2NVC");
  CHRwsp(fout);
  INTw(fout,21,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnodim,5);
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { for(jcount=0;jcount<ydn->nnodim;jcount++)
    { DBLwa(fout,ydn->d2nvc[jcount][icount]);
      CHRwsp(fout);
    }
    CHRwcr(fout);
  }
  CHRwcr(fout); 

  CHRw(fout,"/YD/YDN/I1NOBF");
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { INTw(fout,ydn->i1nobf[icount],4);   
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDN/I1NOPR");
  CHRwsp(fout);
  INTw(fout,ydn->nnopo,5);
  CHRwcr(fout);
  for(icount=0;icount<ydn->nnopo;icount++)
  { INTw(fout,ydn->i1nopr[icount],4);   
    CHRwsp(fout);
  }
  CHRwcr(fout);
}

/****************** Output ******************************************/
static void Ywdo(ydo,fout)
  YDO ydo; FILE *fout;
{ INT icount;
  
  CHRw(fout,"     /*   Output     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDO/MOHYS");
  CHRwsp(fout);
  INTw(fout,ydo->mohys,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/NOHYS");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/DOHYP");
  CHRwsp(fout);
  DBLwa(fout,ydo->dohyp);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYS");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLwa(fout,ydo->d1ohys[icount]);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYC");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLwa(fout,ydo->d1ohyc[icount]);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYF");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLwa(fout,ydo->d1ohyf[icount]);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYT");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLwa(fout,ydo->d1ohyt[icount]);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYX");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLwa(fout,ydo->d1ohyx[icount]);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/D1OHYY");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { DBLwa(fout,ydo->d1ohyy[icount]);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDO/I1OHYT");
  CHRwsp(fout);
  INTw(fout,ydo->nohys,5);
  CHRwcr(fout);
  for(icount=0;icount<ydo->nohys;icount++)
  { INTw(fout,ydo->i1ohyt[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);
}

/****************** Properties **************************************/
static void Ywdp(ydp,fout)
  YDP ydp; FILE *fout;
{ INT icount;
  
  CHRw(fout,"     /*   Properties     */");
  CHRwcr(fout);
  CHRw(fout,"/YD/YDP/MPROP");
  CHRwsp(fout);
  INTw(fout,ydp->mprop,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/NPROP");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEFT");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1peft[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);
  
  CHRw(fout,"/YD/YDP/D1PEGT");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pegt[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEGS");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pegs[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);
  
  CHRw(fout,"/YD/YDP/D1PEKS");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1peks[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PELA");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pela[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);
  
  CHRw(fout,"/YD/YDP/D1PEMU");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pemu[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEPE");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pepe[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEPC");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pepc[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEPF");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pepf[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PBIF");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pbif[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PCOH");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pcoh[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PICF");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1picf[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PEFR");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pefr[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PERO");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pero[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PJRC");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pjrc[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PJCS");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pjcs[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);
  
  CHRw(fout,"/YD/YDP/D1PJSL");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pjsl[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PNAF");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pnaf[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PNAI");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pnai[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PNAP");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pnap[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PNAT");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pnat[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PNAX");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pnax[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PNAY");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pnay[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PNXX");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pnxx[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PNXY");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pnxy[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PNYX");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pnyx[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PNYY");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1pnyy[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/D1PSEM");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { DBLw(fout,ydp->d1psem[icount],18);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PEFR");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1pefr[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PEJP");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1pejp[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PEMB");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1pemb[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PEMN");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1pemn[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PNFX");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1pnfx[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PNFY");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1pnfy[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PSDE");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1psde[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);

  CHRw(fout,"/YD/YDP/I1PTYP");
  CHRwsp(fout);
  INTw(fout,ydp->nprop,5);
  CHRwcr(fout);
  for(icount=0;icount<ydp->nprop;icount++)
  { INTw(fout,ydp->i1ptyp[icount],5);
    CHRwsp(fout);
  }
  CHRwcr(fout);
}

/****************** Mesh ******************************************/
static void Ywgid(yde,ydn,fout)
  YDE yde;YDN ydn; FILE *fout;
{ INT icount,jcount;
  
  CHRw(fout,"MESH    dimension 2 ElemType Triangle  Nnode ");
  INTw(fout,yde->nelno, 4);
  CHRwsp(fout);
  CHRwcr(fout);
  CHRw(fout,"Coordinates ");
  CHRwcr(fout);

  for(icount=0;icount<ydn->nnopo;icount++)
  { INTw(fout,icount+1,5);
    CHRwsp(fout);
    DBLwa(fout,ydn->d2ncc[0][icount]);
    CHRwsp(fout);
    DBLwa(fout,ydn->d2ncc[1][icount]);
    CHRwsp(fout);
    CHRwcr(fout);
  }

  CHRw(fout,"End Coordinates ");
  CHRwcr(fout);
  CHRw(fout,"Elements ");
  CHRwcr(fout);

  for(icount=0;icount<yde->nelem;icount++)
  { INTw(fout,icount+1,5);
    CHRwsp(fout);
    for(jcount=0;jcount<yde->nelno;jcount++)
    { INTw(fout,yde->i2elto[jcount][icount]+1,10);
      CHRwsp(fout);
	}
    INTw(fout,1,2);
    CHRwcr(fout);
  }
  CHRw(fout,"End Elements ");
  CHRwcr(fout);
}

/***************************PUBLIC***********************************/
void Ywd(namep,yd)
  CHR *namep; YD yd;
{ YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDJ ydj=&(yd->ydj);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDP ydp=&(yd->ydp);
  CHR namef[300];
  CHR cindex[50];
  FILE *fp=FILENULL;
  FILE *fgid=FILENULL;

  if((ydc->icsavf)!=0&&(ydc->ncstep)!=0&&
	((((ydc->ncstep+1)%ydc->icsavf)==0)||(ydc->ncstep==(ydc->mcstep-1))))
  { CHRcpynoext(namef,namep);
    SINTw(cindex,ydc->icouti-1,0);
	CHRcat(namef,cindex);
	CHRcat(namef,"_restart.y");
	fp=fopen(namef,"w");
	CHRcpynoext(namef,namep);
    CHRcat(namef,".msh");
    fgid=fopen(namef,"w");
    if((fp!=FILENULL)&&(fgid!=FILENULL))
    { Ywdc(ydc,fp);			/* Write Control variables into the input file		*/
	  Ywde(yde,fp);			/* Write Elements variables into the input file		*/
	  Ywdj(ydj,fp);			/* Write Joint variables into the input file		*/
	  Ywdi(ydi,fp);			/* Write Interaction variables into the input file	*/
	  Ywdn(ydn,fp);			/* Write Nodes variables into the input file		*/
	  Ywdp(ydp,fp);			/* Write Properties variables into the input file	*/
	  Ywdo(ydo,fp);			/* Write Output variables into the input file		*/
	  Ywgid(yde,ydn,fgid);  /* Write GiD mesh file								*/
	  CHRw(fp,"$YDOIT"); CHRwcr(fp);
	  CHRw(fp,"$YSTOP"); CHRwcr(fp);
	}
	fclose(fp);
	fclose(fgid);
} }

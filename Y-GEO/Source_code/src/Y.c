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
/* File   Y.c */
#include "Yproto.h"
#include <time.h>
#include "paraview.h" //! PARAVIEW

main(argc, argv)
  INT argc; char **argv;
{ CHR c1name[300];         /* name of the problem i.e. input file */
  struct YD_struct yd;     /* Y database                          */
  YDC ydc=&(yd.ydc);       /* Y control database                  */
  YDE yde=&(yd.yde);       /* Y element database                  */
  YDI ydi=&(yd.ydi);       /* Y interaction database              */
  YDN ydn=&(yd.ydn);       /* Y node database                     */
  YDB ydb=&(yd.ydb);       /* Y borehole database                 */
  YDS yds=&(yd.yds);       /* Y source (inter. fluid) database    */
  YDO ydo=&(yd.ydo);       /* Y output database                   */
  YDPE ydpe=&(yd.ydpe);    /* Y property database  for elements   */
  YDPN ydpn=&(yd.ydpn);    /* Y property database  for nodes (BC) */
  YDPJ ydpj=&(yd.ydpj);    /* Y property database  for joints     */
  YDPM ydpm=&(yd.ydpm);    /* Y property database  for meshing    */
  YDFN ydfn=&(yd.ydfn);    /* Y property database for DFN         */
  YDIS ydis=&(yd.ydis);    /* Y in-situ stress parameters         */
  YDHF ydhf=&(yd.ydhf);    /* Y hydro-frac parameters             */
  YDSM ydsm=&(yd.ydsm);    /* Y seismic monitoring parameters     */
  YDR ydr=&(yd.ydr);       /* Y reference points database Y-RC    */
  YDSB ydsb=&(yd.ydsb);    /* Y steel-bar element database Y-RC   */
  YDPS ydps=&(yd.ydps);    /* Y property database  for steel Y-RC */
  
  time_t beginTime, endTime;
  time(&beginTime);
  
  /* get name of the problem */
  if(argv[1]!=NULL)
  { CHRcpy(c1name,argv[1]); 

/* process data */
    CHRw(stdout,"Copyright (C) 2011, Dr. Antonio Munjiza \n");
    CHRw(stdout,"Y-Geo is developed and maintained by the Geomechanics_Group@UofT. \n");
    CHRw(stdout,"Y-Geo is based on the original Y code that is provided as part of the book entitled \n");
    CHRw(stdout,"The Combined Finite Discrete Element Method. Under supervision of Dr. G. \n");
    CHRw(stdout,"Grasselli (with some help from Prof A. Munjiza) UofT students Mr. O. K. Mahabadi \n");
    CHRw(stdout,"and Mr. A. Lisjak and others have made and will be making all past, current and \n");
    CHRw(stdout,"future modifications to both the 2D and 3D original source codes - with \n");
    CHRw(stdout,"understanding that original copyright remains. \n");
    CHRw(stdout,"Y-Geo is distributed WITHOUT ANY WARRANTY; without even the implied \n");
    CHRw(stdout,"warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. \n");
    CHRw(stdout,"Inclusion of a part or whole of this code into any other commercial or research or \n");
    CHRw(stdout,"other purpose code is not granted without author's written explicit permission. \n");
    CHRw(stdout,"When results using whole or any part of this code are published, Y-Geo code must be \n");
    CHRw(stdout,"mentioned and acknowledgement to the Geomechanics_Group@UoIT and Prof. A. \n");
    CHRw(stdout,"Munjiza must be made. \n");
    CHRw(stdout,"Should you modify this source code, the Copyright (C) on the modified code \n");
    CHRw(stdout,"as a whole belongs to Dr. A. Munjiza regardless of the extent or nature \n");
    CHRw(stdout,"of modifications. \n");
    CHRw(stdout,"Copyright (C) to whole of any code containing any part of this code \n");
    CHRw(stdout,"also belongs to Dr. A.Munjiza. \n");
    CHRw(stdout,"Any code comprising any part of this source code \n");
    CHRw(stdout,"must be called Y program.\n");
    CHRw(stdout,"If you do not agree with this, you are not allowed to do \n");
    CHRw(stdout,"any modifications to any part of this source code or include \n");
    CHRw(stdout,"any part of it in any other program. \n\n");
    CHRw(stdout,"Mahabadi, O., Lisjak, A., Munjiza, A., and Grasselli, G. (2012).  \n");
    CHRw(stdout,"Y-Geo: New Combined Finite-Discrete Element Numerical Code \n");
    CHRw(stdout,"for Geomechanical Applications. \n");
    CHRw(stdout,"International Journal of Geomechanics. 12, SPECIAL ISSUE:  \n");
    CHRw(stdout,"Advances in Modeling Rock Engineering Problems, \n");
    CHRw(stdout,"doi: http://dx.doi.org/10.1061/(ASCE)GM.1943-5622.0000216 \n\n");
//Mahabadi, O., Lisjak, A., Munjiza, A., and Grasselli, G. (2012). 
//Y-Geo: New Combined Finite-Discrete Element Numerical Code for Geomechanical Applications.
//International Journal of Geomechanics. 12, SPECIAL ISSUE: Advances in Modeling Rock Engineering Problems, 
//(doi: http://dx.doi.org/10.1061/(ASCE)GM.1943-5622.0000216)

    ydc->finp=FILENULL; ydc->fcheck=FILENULL;
    while(Yrd(c1name,&yd)>0)                 /* Process while any input */
    { CHRw(stdout,"NEW INPUT: "); CHRw(stdout, c1name); CHRwcr(stdout);
      for(ydc->ncstep=ydc->ncstep;ydc->ncstep<ydc->mcstep;ydc->ncstep++)
      { Ymd(ydc,yde,ydi,ydn,ydpe,ydpn,ydpm,ydfn);                      /* mesh elements            */
        Yfd(ydc,yde,ydn,ydo,ydpe,ydpn, ydpj,ydis,ydfn,ydhf,ydsm,ydsb); /* nodal forces             */
        Yhd(ydc,yde,ydn,ydo,ydpe,ydpn, ydpj,ydis,ydhf,ydfn);           /* hydrofrac                */
        Ybor(ydc,yde,ydn,ydb,yds,ydpe,ydpj,ydpn);                      /* borholes, inter. fluid   */
        Ycd(ydc,yde,ydi,ydn,ydpe,ydpn);                                /* contact detection        */
        Yrb(ydc,yde,ydn,ydpe,ydpj,ydps,ydr,ydsb,ydo);                  /* RBAR steel elements  Y-RC*/
        Yid(ydc,yde,ydi,ydn,ydo,ydpe,ydpn, ydpj,ydpm,ydfn);            /* interaction              */
        Ysd(ydc,yde,ydn,ydo,ydpe,ydpn);                                /* solve equations          */
	Yod(c1name,&yd);                                               /* output results           */
        Yfrd(ydc,yde,ydi,ydn,ydpe,ydpn,ydpj,ydpm);                     /* fracture                 */
        ydc->dctime=ydc->dctime+ydc->dcstec;                           /* update time              */
      }
    }
    pv_free();
    time(&endTime);
    printf("\nStart time: ");
    printf(ctime(&beginTime));
    printf("End time:   ");
    printf(ctime(&endTime));
    printf("Total run time (seconds) = %.1lf \n\n", difftime(endTime, beginTime));
    CHRw(stderr,"   ***** Y HAS ORDERLY FINISHED *****");
    CHRwcr(stderr);
    CHRw(stderr,"Press a key to continue");
    CHRwcr(stderr);
    getchar();
  }
  else
  { CHRw(stdout,"Double click the data file to run the program.\n");
  } 
}

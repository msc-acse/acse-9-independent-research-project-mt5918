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

main(argc, argv)
  INT argc; char **argv;
{ CHR c1name[300];         /* name of the problem i.e. input file */
  struct YD_struct yd;     /* Y database                          */
  YDC ydc=&(yd.ydc);       /* Y control database                  */
  YDE yde=&(yd.yde);       /* Y element database                  */
  YDJ ydj=&(yd.ydj);	   /* Y joint database					  */
  YDI ydi=&(yd.ydi);       /* Y interaction database              */
  YDN ydn=&(yd.ydn);       /* Y node database                     */
  YDO ydo=&(yd.ydo);       /* Y output database                   */
  YDP ydp=&(yd.ydp);       /* Y property database                 */

  /* get name of the problem */
  if(argv[1]!=NULL)
  { CHRcpy(c1name,argv[1]);
	
    /* process data */
	CHRw(stdout,"Copyright (C) 2000, Dr. Antonio Munjiza. \n");									  
	CHRw(stdout,"This code is provided as part of the book entitled The Combined \n");
	CHRw(stdout,"Finite Discrete Element Method. It is distributed WITHOUT ANY WARRANTY; \n");
	CHRw(stdout,"without even the implied warranty of MERCHANTABILITY or FITNESS FOR A \n");
	CHRw(stdout,"PARTICULAR PURPOSE. Inclusion of a part or whole of this code into any other \n");
	CHRw(stdout,"commercial or research or other purpose code is not granted without author's \n");
	CHRw(stdout,"written explicit permission. \n");
	CHRw(stdout,"When results using whole or any part of this code \n");
	CHRw(stdout,"are published, Y code must be mentioned and acknowledgement to \n");
	CHRw(stdout,"Dr Munjiza must be made. \n");
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
	
    while(Yrd(c1name,&yd)>0)                 /* Process while any input */
    { CHRw(stdout,"NEW INPUT"); CHRwcr(stdout);
      for(ydc->ncstep=ydc->ncstep;ydc->ncstep<ydc->mcstep;ydc->ncstep++)
      {
		  Ymd(ydc,yde,ydi,ydn,ydp);							/* mesh elements           */
		  Yad(ydc,yde,ydj,ydn,ydp);							/* aperture configuration  */ // added by Qinghua
	    Yfd(yde,ydj,ydn,ydp,ydi,ydc);						/* nodal forces            */
        Ycd(ydc,yde,ydj,ydi,ydn,ydp);						/* contact detection       */
        Yid(ydc,yde,ydj,ydi,ydn,ydp);						/* interaction             */
        Ysd(ydc,yde,ydn,ydo,ydp);							/* solve equations         */
		Yod(c1name,&yd);									/* output results          */
        Yfrd(ydc,yde,ydi,ydn,ydp);							/* fracture                */
        ydc->dctime=ydc->dctime+ydc->dcstec;				/* update time             */
		Ywd(c1name,&yd);									/* write restart file	   */ // added by Qinghua
		if((ydc->ncstep%ydc->icoutf)==0)					/* prompt progress		   */
 		{ INTw(stderr,ydc->ncstep/ydc->icoutf,5); CHRwsp(stderr); CHRw(stderr,"/");
		  INTw(stderr,ydc->mcstep/ydc->icoutf,5); CHRwcr(stderr);
		}
    } }
    CHRw(stderr,"   ***** Y HAS ORDERLY FINISHED *****");
    CHRwcr(stderr);
  }
  else
  { CHRw(stdout,"Double click the data file to run the program.\n");
  }
}
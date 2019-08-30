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
/* file Yd.h  Y data base description */
#include "Ytypes.h"
#ifndef FRAMEINCL
#include "frame.h"
#define FRAMEINCL
#endif
#ifndef YTYPESINCL
#include "Ytypes.h"
#define YTYPESINCL
#endif

typedef struct YDC_struct *YDC; 
struct YDC_struct
{ INT  mcstep, ncstep;      /* maximum/current number of time steps   */
  FILE *finp, *fcheck;
  DBL  dcgray;				/* gravity y							  */
  DBL  dcsizc;				/* size coord.							  */
  DBL  dcsizf;				/* size force							  */
  DBL  dcsizs;				/* size stress							  */
  DBL  dcsizv;				/* size velocity						  */
  DBL  dcsizd;              /* size displacement					  */
  DBL  dcsiza;              /* size aperture						  */
  DBL  dcstec;              /* current time step size                 */
  DBL  dctime;              /* current time                           */
  DBL  dcgrst;				/* gravity settling time				  */
  DBL  dcrmpt;				/* ramping time for loading application   */
  INT  icsavf;				/* restart file output frequency		  */
  INT  icoutf;              /* write output frequency                 */
  INT  icouti;              /* current write output No                */
  INT  icoutp;              /* output precision - digits per number   */
  INT  icfmty;				/* fric. model 0:Coulomb 1:Barton-Bandis  */
  INT  iciaty;				/* init. aper. type 0:roughness,1:length  */
};

typedef struct YDE_struct *YDE;
struct YDE_struct
{ INT melem, nelem; /* maximum (actual) number of elements            */
  INT melst, nelst; /* maximum (actual) number of elemen. states var. */
  INT melno, nelno; /* maximum (actual) number of elemen. nodes       */
  INT  *i1elcf;     /*[melem]    contacting couple first		 */ 
  INT  *i1elpr;     /*[melem]    element property				 */
  INT  *i1elty;		/*[melem]	  element type					 */
					/*3:new crack; 2:pre-existing; 1:boundary;	 */
					/*0:unbroken; -1:non-joint elemnt			 */
  DBL **d2elst;     /*[melst][melem]    - element state          */
  INT **i2elto;     /*[melno][melem]    - element topology       */
  INT **i2eljp;     /*[melno][melem]    - element-joint relation */
  DBL **d2tcs ;		/*[3][melem]		 - element stress tensor */
  DBL  *d1elsf;		/*[melem]			 - element sliding force */
  DBL **d2elcf;		/*[melem][2]		 - element contact force */
  DBL **d2elfs;		/*[melem][3] - shear strength at integ. point*/
  DBL  *d1eley;		/*[melem]	  - kinetic energy when yielding */
  DBL  *d1eles;		/*[melem]	  - seismic energy				 */
  DBL  *d1elme;		/*[melem]	  - seismic event magnitude		 */
  INT  *i1elyi;		/*[melem]	  - yielding indicator			 */
};

typedef struct YDJ_struct *YDJ;						// added by Qinghua
struct YDJ_struct
{ INT  njoint;		  /* number of joint elements					  */
  INT *i1jtid;		  /* joint element index in element topology list */
  DBL *d1jkni;		  /* [mjelem] joint initial normal stiffness	  */
  DBL *d1jknc;		  /* [mjelem] joint current normal stiffness	  */
  DBL *d1jksc;		  /* [mjelem] joint current shear stiffness		  */
  DBL *d1jnst;		  /* [mjelem] joint normal stress				  */
  DBL *d1jsst;		  /* [mjelem] joint shear stress				  */
  DBL *d1japi;		  /* [mjelem] joint initial normal aperture		  */
  DBL *d1japc;		  /* [mjelem] joint current normal aperture		  */
  DBL *d1japh;		  /* [mjelem] joint hydraulic aperture			  */
  DBL *d1japr;		  /* [mjelem] joint residual normal aperture	  */
  DBL *d1jsdc;		  /* [mjelem] joint current shear displacement	  */
  DBL *d1jdlc;		  /* [mjelem] joint current shear dilation		  */
  DBL *d1jsdp;		  /* [mjelem] joint peak shear displacement		  */
  DBL *d1jefl;		  /* [mjelem] joint effective fracture length	  */
  DBL *d1jjrc;		  /* [mjelem] joint roughness coefficient		  */
  DBL *d1jjcs;		  /* [mjelem] joint compressive strength (MPa)	  */
  DBL *d1jphi;		  /* [mjelem] joint asperity angle (deg)		  */
  DBL *d1jfmd;		  /* [mjelem] joint (if new crack) failure mode   */
  DBL *d1jfpr;		  /* [mjelem] joint fluid pressure				  */
  DBL *d1jfet;		  /* [mjelem] joint fluid entering time			  */
};

typedef struct YDI_struct *YDI;
struct YDI_struct
{ INT micoup, nicoup;	/*maximum possible number of contacting couples */
  INT    iiecff;		/*interaction element contact. couple free first*/
  DBL    diedi;			/*travel since last detection					*/
  DBL    diezon;		/*buffer zone size								*/
  DBL   *d1iesl;		/*[mcoup] contact sliding                       */
  INT   *i1iecn;		/*[mcoup] couple next                           */
  INT   *i1iect;		/*[mcoup] couple target                         */
  INT   mistate;        /* number of states for d2sldis                 */
  DBL **d2sldis;		/*[mistate][mcoup] sliding distance of couples  */
};

typedef struct YDN_struct *YDN;
struct YDN_struct
{ INT mnodim, nnodim;   /* max(actual) nodal dimensions number          */
  INT mnopo, nnopo;     /* maximum (actual) number of nodal points      */
  DBL   *d1nmct;        /* [mnopo] nodal mass current translation       */
  DBL  **d2ncc;         /* [mnodim][mnopo] nodal coordinate current     */
  DBL  **d2nci;         /* [mnodim][mnopo] nodal coordinate initial     */
  DBL  **d2nfc;         /* [mnodim][mnopo] nodal force current          */
  DBL  **d2nfcon;       /* [mnodim][mnopo] nodal force due to contact   */
  DBL  **d2nvc;         /* [mnodim][mnopo] nodal velocity current       */
  INT  *i1nobf;         /* [mnopo] nodal boundary >0 is boundary        */
  INT  *i1nopr;         /* [mnopo] nodal property                       */
};

typedef struct YDO_struct *YDO;
struct YDO_struct
{ INT mohys, nohys; /* maximum (actual) number of hystory variables      */
  DBL     dohyp;  /* output hystory accuracy                             */
  DBL   *d1ohyf;  /*[mohys] output hystory factor to scale state         */
  DBL   *d1ohyc;  /*[mohys] output hystory factor to scale time          */
  DBL   *d1ohys;  /*[mohys] output hystory state                         */
  DBL   *d1ohyt;  /*[mohys] output hystory time                          */
  DBL   *d1ohyx;  /*[mohys] output history x coordinate of the point     */
  DBL   *d1ohyy;  /*[mohys] output history y coordinate of the point     */
  DBL   *d1ohyz;  /*[mohys] output history z coordinate of the point     */
  FILE **f2ohyf;  /*[mohys] output history files                        */
  INT   *i1ohyt;  /*[mohys] output hystory type, i.e. which variable     */
};

typedef struct YDP_struct *YDP;
struct YDP_struct
{ INT mprop, nprop; /* maximum (actual) number of properties             */
  DBL   *d1peca;  /*[mprop] child age - procreation                      */
  DBL   *d1pecl;  /*[mprop] child life - interval for procreation        */
  DBL   *d1pefs;  /*[mprop] ultimate shear stress at joint			     */
  DBL   *d1peft;  /*[mprop] ultimate tensile stress at joint             */
  DBL   *d1pegt;  /*[mprop] ultimate fracture energy mode I	             */
  DBL   *d1pegs;  /*[mprop] ultimate fracture energy mode II             */
  DBL   *d1peks;  /*[mprop] dpeks=2hbeta*sqrt(E*ro) in 2D or 3D,0<beta<1 */
  DBL   *d1pela;  /*[mprop] property lamda - Lame elastic constant       */
  DBL   *d1pemu;  /*[mprop] property mu    - Lame elastic constant       */
  DBL   *d1pepe;  /*[mprop] property penalty parameter (element)         */
  DBL   *d1pepc;  /*[mprop] property penalty parameter (contact)         */
  DBL   *d1pepf;  /*[mprop] pore fluid pressure				             */
  DBL   *d1pbif;  /*[mprop] bedding interface friction coefficient       */
  DBL   *d1pera;  /*[mprop] property radius of sphere                    */
  DBL   *d1pero;  /*[mprop] property ro    - density                     */
  DBL   *d1pevi;  /*[mprop] viscosity for  granular flow                 */
  DBL   *d1pefr;  /*[mprop] coefficient of friction                      */
  DBL   *d1picf;  /*[mprop] internal friction angle                      */
  DBL   *d1pcoh;  /*[mprop] the cohesion constant                        */
  DBL   *d1pnaf;  /*[mprop] amplitude factor all ampltd multp by it      */
  DBL   *d1pnai;  /*[mprop] amplitude factor increment each time step    */
  DBL   *d1pnap;  /*[mprop] amplitude of element surface pressure        */
  DBL   *d1pnat;  /*[mprop] amplitude of element surface traction        */
  DBL   *d1pnax;  /*[mprop] amplitude of force/velocity x                */
  DBL   *d1pnay;  /*[mprop] amplitude of force/velocity y                */
  DBL   *d1pnaz;  /*[mprop] amplitude of force/velocity z                */
  DBL   *d1pnxx;  /*[mprop] direction of local x                         */
  DBL   *d1pnxy;  /*[mprop] direction of local x                         */
  DBL   *d1pnxz;  /*[mprop] direction of local x                         */
  DBL   *d1pnyx;  /*[mprop] direction of local y                         */
  DBL   *d1pnyy;  /*[mprop] direction of local y                         */
  DBL   *d1pnyz;  /*[mprop] direction of local y                         */
  DBL   *d1pnzx;  /*[mprop] direction of local z                         */
  DBL   *d1pnzy;  /*[mprop] direction of local z                         */
  DBL   *d1pnzz;  /*[mprop] direction of local z                         */
  DBL   *d1psem;  /*[mprop] maximum tensile stretch                      */
  DBL   *d1pjrc;  /*[mprop] joint roughness coefficient					 */
  DBL   *d1pjcs;  /*[mprop] joint compressive strength					 */
  DBL   *d1pjsl;  /*[mprop] joint sample length							 */
  INT   *i1pecn;  /*[mprop] property No to be assigned to child's nodes  */
  INT   *i1pecp;  /*[mprop] permanent property to be assigned to child   */
  INT   *i1pect;  /*[mprop] temporary property to be assigned to child   */
  INT   *i1pefr;  /*[mprop] if >, fracture                               */
  INT   *i1pejp;  /*[mprop] joint property; if<0, no joints              */
  INT   *i1pemb;  /*[mprop] mark boundary nodes 1 yes 0 no               */
  INT   *i1pemn;  /*[mprop] number of mesh refinements                   */
  INT   *i1pnfx;  /*[mprop] fixity x direction 1 force; 2 acc. 3 vel.    */
  INT   *i1pnfy;  /*[mprop] fixity y direction 1 force; 2 acc. 3 vel.    */
  INT   *i1pnfz;  /*[mprop] fixity z direction 1 force; 2 acc. 3 vel.    */  
  INT   *i1pnib;  /*[mprop] 1 if borhole regardles boundary, else 0      */  
  INT   *i1ptyp;  /*[mprop] property type                                */
  INT   *i1psde;  /*[mprop] state damage elastic id                      */
};

typedef struct YD_struct *YD;
struct YD_struct
{ struct YDC_struct ydc;   /* control structure                  */
  struct YDE_struct yde;   /* element description structure      */
  struct YDJ_struct ydj;   /* joint element structure		     */
  struct YDI_struct ydi;   /* interaction structure              */
  struct YDN_struct ydn;   /* node description structure         */
  struct YDO_struct ydo;   /* output description structure       */
  struct YDP_struct ydp;   /* property description structure     */
};
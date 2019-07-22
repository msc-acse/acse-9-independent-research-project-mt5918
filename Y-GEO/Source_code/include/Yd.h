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
{ INT  mcstep, ncstep;       /* maximum/current number of time steps                  */
  FILE *finp, *fcheck;

  DBL  dcgray;               /* gravity y                                              */
  DBL  dcsizc;               /* size coord.                                            */
  DBL  dcsizf;               /* size force                                             */
  DBL  dcsizs;               /* size stress                                            */
  DBL  dcsizv;               /* size velocity                                          */
  DBL  dcstec;               /* current time step size                                 */
  DBL  dctime;               /* current time                                           */
  INT  icoutf;               /* write output frequency                                 */
  INT  icoutrf;              /* write output reduced frequency                         */
  INT  icoutnf;              /* threshold number of fractures beyond which use icoutrf */
  INT  icouti;               /* current write output No                                */
  INT  icoutp;               /* output precision - digits per number                   */
  INT  icresf;               /* output frequency of restart data                       */
};

typedef struct YDE_struct *YDE;
struct YDE_struct
{ INT melem, nelem;          /* maximum (actual) number of elements                   */
  INT melst, nelst;          /* maximum (actual) number of elemen. states var.        */
  INT melno, nelno;          /* maximum (actual) number of elemen. nodes              */
  INT nelemst;               /* actual number of elements @ the beginning             */
  INT nebrk;                 /* actual number of joint elements broken                */
  INT netbrk;                /* number of joint elements broken in one time step      */
  INT nesft;                 /* number of joint elements softened                     */
  INT netsft;                /* number of joint elements softened in one time step    */

  INT   *i1elcf;             /*[melem]    contacting couple first                     */
  INT   *i1elpr;             /*[melem]    element property                            */
  INT   *i1elprtmp;          /*[melem]    element property  used for meshing          */
  DBL   *d1elfs;             /*[melem]    shear strength at joint (Mohr-Coulomb)      */
  DBL   **d2elst;            /*[melst][melem]    - element state                      */
  INT   **i2elto;            /*[melno][melem]    - element topology                   */
  INT   **i2eltost;          /*[melno][melem]    - element topology @ the beginning   */
  INT   *i1ebrk;             /*[melem]         broken joint elements                  */
  DBL   **d2ecbrk;           /*[mnodim][melem] joint element breakage coordinates     */
  DBL   **d2ecbrk_NEW;           /*[mnodim][melem] joint element breakage coordinates     */
  DBL   *d1etbrk;            /*[melem]         joint element breakage time            */
  DBL   *d1elbrk;            /*[melem]         length of broken joint elements        */
  DBL   *d1efe;              /*[melem]         joint element fracture energy          */
  INT   *i1esft;             /*[melem]         softened joint elements                */
  DBL   **d2ecsft;           /*[mnodim][melem] joint element yielding coordinates     */
  DBL   *d1etsft;            /*[melem]         joint element yielding time            */
  INT   *i1esftf;            /*[melem]         softened element flag                  */
  DBL   *d1ebrkf;            /*[melem]         broken element flag                    */
  DBL   *d1eike;             /*[melem]    joint element initial kinetic energy        */
  DBL   *d1edke;             /*[melem]    joint element differential kinetic energy   */
  DBL   *d1etmke;            /*[melem]    joint element maximum kinetic energy time   */
  DBL   *d1elfr;             /*[melem] friction angle of elements next to DFN cracks  */
  DBL   *d1elpe;             /*[melem] normal penalty of elements next to DFN cracks  */
  DBL   *d1elpt;             /*[melem] tangent penalty of elements next to DFN cracks */
  INT   **i2elnext;          /*[4][melem] element & edge indeces next to joints       */
  INT   **i2eledge;          /*[nelno][melem] flags for intact/broken element edges   */
  INT   *i1edfnf;            /*[melem] flag for joint element belonging to a DFN      */
  INT   *i1edft;             /*[melem] flag for triangles next to DFN crack type 3    */
  INT   *d1etike;            /*[melem] initial time of KE monitoring window           */
  
  DBL   **d2eldmg;           /*[3][melem] maximum damage coefficient of joint elements*/
  
  DBL   **d2elstr;           /*[4][melem] stress tensor of elements                   */
};

typedef struct YDI_struct *YDI;
struct YDI_struct
{ INT micoup, nicoup;        /* maximum possible number of contacting couples         */
  INT    iiecff;             /* interaction element contact. couple free first        */

  DBL    diedi;              /* travel since last detection                           */
  DBL    diezon;             /* buffer zone size                                      */
  DBL   *d1iesl;             /*[mcoup] contact sliding                                */
  INT   *i1iecn;             /*[mcoup] couple next                                    */
  INT   *i1iect;             /*[mcoup] couple target                                  */
  INT   mistate;             /* number of states for d2sldis                          */
  DBL   **d2sldis;           /*[mistate][mcoup] sliding distance btw  couples         */
};

typedef struct YDN_struct *YDN;
struct YDN_struct
{ INT mnodim, nnodim;        /* max(actual) nodal dimensions number                   */
  INT mnopo, nnopo;          /* maximum (actual) number of nodal points               */
  INT nnopst;                /* actual number of nodal points at the beginning        */

  DBL   *d1nmct;             /* [mnopo] nodal mass current translation                */
  DBL  **d2ncc;              /* [mnodim][mnopo] nodal coordinate current              */
  DBL  **d2nci;              /* [mnodim][mnopo] nodal coordinate initial              */
  DBL  **d2nfc;              /* [mnodim][mnopo] nodal force current                   */
  DBL  **d2nvc;              /* [mnodim][mnopo] nodal velocity current                */
  INT  *i1nobf;              /* [mnopo] nodal boundary >0 is boundary                 */
  INT  *i1nopr;              /* [mnopo] nodal boundary condition                      */
  INT  *i1nowe;              /* [mnopo] hydrofrac boundary, >0 node is wet            */
  INT  **i2noid;             /* [2][mnopo] hydrofrac, IDs of cw and ccw nodes         */
  DBL   *d1nfp;              /* [mnopo] hydrofrac, nodal fluid pressure               */
  DBL  **d2nc0;              /* [mnodim][mnopo] nodal coordinate at time step 0       */
};

typedef struct YDB_struct *YDB;
struct YDB_struct
{ INT mborh, nborh;          /* maximum (actual) number of boreholes                  */
  INT mbdim, nbdim;          /* maximum (actual) number of dimensions of boreholes    */
  INT nbpaf;                 /* number of amplitude factors                           */

  DBL   **d2bca;             /* [mbdim][mborh] coordinates of point A in borehole     */
  DBL   **d2bcb;             /* [mbdim][mborh] coordinates of point B in borehole     */
  DBL   *d1brad;             /* [mborh] radii of boreholes                            */
  DBL   *d1bpaf;             /* [nbpf] amplitude factor of pressure amplitude         */
  DBL   *d1bpts;             /* [mborh] start time of pressure load on boreholes      */
  DBL   *d1bpte;             /* [mborh] end time of pressure load on boreholes        */
  DBL   *d1bvdt;             /* [mborh] velocity of detonation on boreholes           */
  DBL   *d1bprs;             /* [mborh] amplitudes of pressure for each borehole      */
  DBL   dblmax;              /* max length of all boreholes                           */
  DBL   dbbuf;               /* buffer (max. size of element)                         */
};

typedef struct YDS_struct *YDS;
struct YDS_struct
{ INT msour, nsour;          /* maximum (actual) number of sources                    */
  INT msdim, nsdim;          /* maximum (actual) number of dimensions of sources      */
  INT nspaf;                 /* number of pressure amplitude factors (p=p(t) )        */
  INT nssaf;                 /* number of pressure amplitude factors (s=s(r) )        */

  DBL   **d2scs;             /* [msdim][msour] coordinates of sources                 */
  DBL   *d1spaf;             /* [nspaf] pressure amplitude factors (p=p(t) )          */
  DBL   *d1ssaf;             /* [nssaf] pressure amplitude factors (s=s(r) )          */
  DBL   *d1spts;             /* [msour] start time of pressure load                   */
  DBL   *d1spte;             /* [msour] end time of pressure load                     */
  DBL   *d1svpr;             /* [msour] velocity of pressure propagation              */
  DBL   *d1sprs;             /* [msour] amplitudes of pressures (p=p(t) )             */
  DBL   *d1ssir;             /* [msour] source initial radius                         */
  DBL   dsbuf;               /* buffer (max. size of element)                         */
};

typedef struct YDO_struct *YDO;
struct YDO_struct
{ INT mohys, nohys;          /* maximum (actual) number of history variables          */

  DBL     dohyp;             /* output history accuracy                               */

  DBL   *d1ohyf;             /*[mohys] output history factor to scale state           */
  DBL   *d1ohyc;             /*[mohys] output history factor to scale time            */
  DBL   *d1ohys;             /*[mohys] output history state                           */
  DBL   *d1ohyt;             /*[mohys] output history time                            */
  DBL   *d1ohyx;             /*[mohys] output history x coordinate of the point       */
  DBL   *d1ohyy;             /*[mohys] output history y coordinate of the point       */
  DBL   *d1ohyz;             /*[mohys] output history z coordinate of the point       */
  FILE  **f2ohyf;            /*[mohys] output history files                           */
  INT   *i1ohyt;             /*[mohys] output history type, i.e. which variable       */
  FILE  **f2orsf;            /*[...] restart data files                               */
  FILE  *foebrk;             /* output time, position and energy of broken elements   */
  FILE  *foesft;             /* output time, position and energy of softened elements */
  FILE  *fohyfr;             /* time, flow rate and fluid pressure for hydrofrac      */
  FILE  *fofrac;             /* total fracture area in the model                      */
};

/* Y Database Properties - Elements */
typedef struct YDPE_struct *YDPE;
struct YDPE_struct
{ INT mprop, nprop;          /* maximum (actual) number of properties                 */

  DBL   *d1peca;             /*[mprop] child age - procreation                        */
  DBL   *d1pecl;             /*[mprop] child life - interval for procreation          */
  DBL   *d1peks;             /*[mprop] dpeks=2hbeta*sqrt(E*ro) in 2D or 3D,0<beta<1   */
  DBL   *d1pela;             /*[mprop] property lamda - Lame elastic constant         */
  DBL   *d1pemu;             /*[mprop] property mu    - Lame elastic constant         */
  DBL   *d1pepe;             /*[mprop] contact penalty parameter for Y-RC             */
  //DBL   *d1pept;             /*[mprop] tangential penalty param., (0.1*ydpj->d1pjpe)  */
  //DBL   *d1pefr;             /*[mprop] Coloumb friction                               */
  DBL   *d1pera;             /*[mprop] property radius of sphere                      */
  DBL   *d1pero;             /*[mprop] property ro    - density                       */
  DBL   *d1pevi;             /*[mprop] viscosity for  granular flow                   */
  
  INT    mperow;             /* Possible combinations of Pr sets for phi, p_n, p_t    */
  DBL   **d2peint;           /*[mprop][5] Coulomb friction, normal penalty, tangential penalty between property sets */
  
  DBL   *d1peem;             /*[mprop] Young's modulus                                */
  DBL   *d1penu;             /*[mprop] Poisson's ratio                                */

  DBL   *d1psem;             /*[mprop] maximum tensile stretch                        */

  INT   *i1pecn;             /*[mprop] property No to be assigned to child's nodes    */
  INT   *i1pecp;             /*[mprop] permanent property to be assigned to child     */
  INT   *i1pect;             /*[mprop] temporary property to be assigned to child     */
  INT   *i1pefr;             /*[mprop] if >, fracture                                 */
  INT   *i1pejp;             /*[mprop] joint property; if<0, no joints                */
  INT   *i1pemb;             /*[mprop] mark boundary nodes 1 yes 0 no                 */
  INT   *i1pemn;             /*[mprop] number of mesh refinements                     */

  INT   *i1pnib;             /*[mprop] 1 if borhole regardles boundary, else 0        */
  INT   *i1ptyp;             /*[mprop] property type                                  */ 
  INT   *i1psde;             /*[mprop] state damage elastic id                        */
  
  INT   *i1pexc;             /*[mprop] excavation flag                                */
  
  INT   *i1usan;             /*[mprop] flag for using anisotropic elasticity (=1:use) */
  DBL   *d1peex;             /*[mprop] Young's modulus E_x                            */
  DBL   *d1peey;             /*[mprop] Young's modulus E_y                            */
  DBL   *d1pemx;             /*[mprop] Poisson's ratio mu_xy                          */
  DBL   *d1pemy;             /*[mprop] Poisson's ratio mu_yx                          */
  DBL   *d1peg;              /*[mprop] Shear modulus G                                */
  
  INT   *i1psup;             /*[mprop] excavation support flag                        */ 
};

/* Y Database Properties for Joints */
typedef struct YDPJ_struct *YDPJ;
struct YDPJ_struct
{ INT mpjset, npjset;        /* maximum (actual) number of Joint Property sets       */

  DBL *d1pjfs;               /* [mpjset] ultimate shear strength at joint            */
  DBL *d1pjft;               /* [mpjset] ultimate tensile strength at joint          */
  DBL *d1pjgf;               /* [mpjset] ultimate fracture energy at joint in tension*/
  DBL *d1pjgs;               /* [mpjset] ultimate fracture energy at joint in shear  */
  DBL *d1pjco;               /* [mpjset] cohesion at joint                           */
  DBL *d1pjfr;               /* [mpjset] friction at joint                           */

  DBL *d1pjpe;               /* [mpjset] fracture penalty parameter                  */

  INT *i1psde;               /* [mpjset] state damage elastic id                     */
  INT *i1ptyp;               /* [mpjset] Property Joint Type                         */
  
  DBL *d1usaf;               /* [mpjset] flag for using anisotropic fracture (=1or2) */
  DBL *d1pjcr;               /* [mpjset] reduced cohesion at joint                   */
  DBL *d1pjfd;               /* [mpjset] reduced friction at joint                   */
  DBL *d1pjtr;               /* [mpjset] reduced ultimate tensile strength at joint  */
  DBL *d1pjgr;               /* [mpjset] reduced ultimate fracture energy at joint in tension */
  DBL *d1pjsr;               /* [mpjset] reduced ultimate fracture energy at joint in shear */
  DBL *d1pjal;               /* [mpjset] layering orientation (0-180)                */
 
  INT iusehy;                /* flag for using hysteretic fracture model             */
};

/* Y Database Properties - Nodes, Nodal Properties or Boundary Conditions for nodes  */
typedef struct YDPN_struct *YDPN;
struct YDPN_struct
{ INT mpnset, npnset;        /* maximum (actual) number of node property sets        */
  INT mpnfact, npnfact;      /* maximum (actual) number of factors                   */
  DBL ***d3pnfac;            /*[2][mbc][mfact] time and amplitude factor             */

  INT   *i1pnfx;             /*[mpnset] fixity x direction 1 force; 2 acc. 3 vel.    */
  INT   *i1pnfy;             /*[mpnset] fixity y direction 1 force; 2 acc. 3 vel.    */
  INT   *i1pnfz;             /*[mpnset] fixity z direction 1 force; 2 acc. 3 vel.    */

  DBL   *d1pnaf;             /*[mpnset] amplitude factor all ampltd multp by it      */
  DBL   *d1pnap;             /*[mpnset] amplitude of element surface pressure        */
  DBL   *d1pnat;             /*[mpnset] amplitude of element surface traction        */
  DBL   *d1pnax;             /*[mpnset] amplitude of force/velocity x                */
  DBL   *d1pnay;             /*[mpnset] amplitude of force/velocity y                */
  DBL   *d1pnaz;             /*[mpnset] amplitude of force/velocity z                */

  DBL   *d1pnxx;             /*[mpnset] direction of local x                         */
  DBL   *d1pnxy;             /*[mpnset] direction of local x                         */
  DBL   *d1pnxz;             /*[mpnset] direction of local x                         */
  DBL   *d1pnyx;             /*[mpnset] direction of local y                         */
  DBL   *d1pnyy;             /*[mpnset] direction of local y                         */
  DBL   *d1pnyz;             /*[mpnset] direction of local y                         */
  DBL   *d1pnzx;             /*[mpnset] direction of local z                         */
  DBL   *d1pnzy;             /*[mpnset] direction of local z                         */
  DBL   *d1pnzz;             /*[mpnset] direction of local z                         */
};

/* Y Database Properties - Meshing */
typedef struct YDPM_struct *YDPM;
struct YDPM_struct
{ INT mpmcom;                /* combination of pr. sets (rows in I2PMSET)            */
  INT mpmcol;                /* Number of columns  in I2PMSET                        */
  INT **i2pmset;             /*[mpmcol][mpmcom] mesh pr.sets i & j together          */
  INT mpmrow;                /* Number of rows in I2PMIJ                             */
  INT **i2pmij;              /* Meshing combinations defining joints btw 2 pr. sets  */
};

/* Y Database Properties - Meshing */
typedef struct YDFN_struct *YDFN;
struct YDFN_struct
{ INT iusefn;                /* flag for using discrete fracture network (=1: use)        */
  INT mdfnfr;                /* maximum number of fractures in the DFN                    */
  INT mdfnno;                /* maximum number of nodes per fracture                      */
  INT **i2dfnn;              /* [mdfnfr][mdfnno] node IDs of each fracture                */
  DBL *d1dffr;               /* [mdfnfr] friction coefficient of each fracture            */
  DBL *d1dfpe;               /* [mdfnfr] normal penalty of each fracture                  */
  DBL *d1dfpt;               /* [mdfnfr] tangential penalty of each fracture              */
  DBL ddfnft;                /* Tensile strength of DFN fractures (when iusefn==2)        */
  DBL ddfnco;                /* Cohesion of DFN fractures (when iusefn==2)                */ 
  DBL ddfngf;                /* Mode I fracture energy of DFN fractures (when iusefn==2)  */ 
  DBL ddfngs;                /* Mode II fracture energy of DFN fractures (when iusefn==2) */
  INT *i1dfft;               /* [mdfnfr] type of DFN fracture (when iusefn==3, 1=broken, 2=cohesive */
};

typedef struct YDIS_struct *YDIS;
struct YDIS_struct
{ 
  INT  iuseis;               /* flag for using in-situ stress (=1: use)               */
  DBL  dcstxx;               /* in-situ stress tensor xx component                    */
  DBL  dcstxy;               /* in-situ stress tensor xy component                    */ 
  DBL  dcstyy;               /* in-situ stress tensor yy component                    */ 
  DBL  dcsyxx;               /* in-situ stress tensor xx component y gradient         */	
  DBL  dcsyxy;               /* in-situ stress tensor xy component y gradient         */	
  DBL  dcsyyy;               /* in-situ stress tensor yy component y gradient         */
  DBL  dcsrfy;               /* y coordinate of the (flat) topographic surface        */
};

typedef struct YDHF_struct *YDHF;
struct YDHF_struct
{ 
  INT  iusehf;               /* flag for using in-situ stress (=1: use)               */
  INT  ihftyp;               /* hydro-frac input type, 1 = pressure, 2 = flow rate    */
  DBL  dhfflp;               /* input fluid pressure                                  */
  DBL  dhfflq;               /* input flow rate                                       */
  INT  hfarow;               /* number of amplitude factors                           */
  DBL **d2hfaf;              /* pressure or flow rate amplitude factor vs time        */
  DBL  fluvol;               /* Fluid volume                                          */
  DBL  flupres;              /* Fluid pressure                                        */
  DBL  flumass;              /* Fluid mass                                            */
  DBL  flurho0;              /* Fluid density at reference pressure p_0               */
  DBL  flupres0;             /* Reference fluid pressure p_0                          */
  DBL  flubulk;              /* Bulk modulus                                          */
  INT  fradim;               /* Dimensionality of fractures (2 = 2D, 3 = 3D)          */
  INT  ihfmsin;              /* Flag for fluid mass initialization                    */
  DBL  gravacc;              /* Gravitational acceleration for hydrostatic pressure   */
  DBL **d2wtlev;             /* Upstream and downstream x-coord and water level       */
};

typedef struct YDSM_struct *YDSM;
struct YDSM_struct
{ 
  INT  iusesm;               /* flag for using alternative monitoring (=1: use)       */
  DBL  dctwle;               /* maximum duration of monitoring window                 */
};

/* reference points */
typedef struct YDR_struct *YDR;
struct YDR_struct
{ INT mnodim, nnodim; /* max(actual) nodal dimensions number                   */
  INT mrdim, nrdim;   /* max(actual) reference point dimensions number         */
  INT mrldm, nrldm;   /* max(actual) local reference point dim. number         */
  INT mrepo, nrepo;   /* maximum (actual) number of ref. nodal points          */
  INT nbrjointrb;
  DBL  **d2rcig;  /* [mnodim][mrepo] coordinate initial global                 */
  DBL  **d2rccg;  /* [mnodim][mrepo] coordinate current global                 */
  DBL  **d2rccl;  /* [mrldm][mrepo]  coordinate current local                  */
  DBL  **d2rvcg;  /* [mnodim][mrepo] velocity current global                   */
  DBL  **d2riLc;  /* [mrldm][mrepo] referent coordinate initial local          */ 
  DBL   **d2rsctr;/* [2][mrepo]carent vector ref.points                        */ 
  INT  **i2relto; /* [2][mrepo]  1D joint element topology                     */
  INT   *i1rmyel; /* [mrepo]         my element                                */
  INT   *i1rrpn;  /* [mrepo]  ref. point next on sing list                     */
  INT   *i1rprop;  /*[mrepo]  properties of ref. points                        */
  INT   *i1refbar; /*[mrepo]  ref. points - bar                                */
  INT   *i1myjoint;   /* [mrepo]  steel joint                                  */
  DBL   *d1rbsig;   /* [mrepo]  stress in rebar                                */
  DBL   *d1rbfrc;   /* [mrepo]  force in rebar                                 */
  DBL   *d1rbstr;   /* [mrepo]  strain in rebar                                */
  DBL   *d1rjsig;   /* [mrepo]  normal (axial) stress in 1D joint              */
  DBL   *d1rjtau;   /* [mrepo]  tangential (shear) stress in 1D joint          */
  DBL   *d1rjslpnor;/* [mrepo]  normal (axial) strain in 1D joint              */
  DBL   *d1rjdelnor;/* [mrepo]  tangential (shear) strain in 1D joint          */
  INT   *i1rjstnor; /* [mrepo]  axial state of 1D joint  (1 = ela, 2 = yielded */
  DBL   **d2rjfrv;  /* [2][mrepo] force vector in 1D joint                     */
  INT   **i2rbedn;  /* [2][mrepo] edge nodes associated with reference point   */
};

/* steel reinforcement elements*/
typedef struct YDSB_struct *YDSB;
struct YDSB_struct
{ INT msdim, nsdim; /* max(actual) reference point dimensions number    */
  INT msbar, nsbar; /* maximum (actual) number of bar elements          */
  INT isfirst;      /* iffirst step for steel bar    yes/no             */   

  DBL    **d2sic;    /*[msdim][msbar] initial steel bar coordinate        */
  INT    *i1srpf;    /*[nsbar]        reference point first               */ 
  INT    *i1sbpr;    /*[msbar]        bar element property                */
  DBL    *d1spea;    /*[msbar]        bar area                            */
  DBL    *d1sdiam;   /*[msbar]        bar diametar                        */ 
  DBL    *d1smmdiam; /*[msbar]        bar diametar  in (mm)               */
  DBL    *d1crlcr;   /*[msbar]        interval of discrete cracks in (mm) */
  INT    *i1sbty;    /*[msbar]        type of bar, 0 = elastic, 1 = plastic (Y-RC formulation), 5 = elastic (Andrea's formulation) */
  INT    *i1sbac;    /*[msbar]        0 = inactivated, 1 = activated      */
};

/* steel reinforcement property*/
typedef struct YDPS_struct *YDPS;
struct YDPS_struct
{ INT mprop, nprop; /* max(actual) number of steel property                                        */
  DBL    *d1young;    /*[mprop] property young modulus of elasticity                               */
  DBL    *d1sfc;      /*[mprop] concrete compressive strength                                      */
  DBL    *d1mpsfc;    /*[mprop] concrete compressive strength in MPa                               */
  DBL    *d1sfy;      /*[mprop] yielding strength of steel                                         */
  DBL    *d1epssh;    /*[mprop] strain - hardening point of steel                                  */
  DBL    *d1sfu;      /*[mprop] ultimate strength  of steel                                        */
  DBL    *d1epsu;     /*[mprop] ultimate strain  of steel                                          */
  DBL    *d1sfbr;     /*[mprop] break strength  of steel                                           */
  DBL    *d1epsbr;    /*[mprop] break strain  of steel                                             */
  DBL    *d1stkn;     /*[mprop] normal stiffness of 1D joint element (used if i1sbty = 5 or 6)     */
  DBL    *d1stkt;     /*[mprop] tangential stiffness of 1D joint element (used if i1sbty = 5 or 6) */
  DBL    *d1styns;    /*[mprop] yield stress in the axial direction  (used if i1sbty = 6)          */
  DBL    *d1strns;    /*[mprop] rupture strain in the axial direction  (used if i1sbty = 6)        */
  DBL    *d1stcoh;    /*[mprop] cohesion of rebar-rock interface                                   */
  DBL    *d1stfri;    /*[mprop] friction coefficient of rebar-rock interface                       */
};



typedef struct YD_struct *YD;
struct YD_struct
{ struct YDC_struct ydc;     /* control structure                                    */
  struct YDE_struct yde;     /* element description structure                        */
  struct YDI_struct ydi;     /* interaction structure                                */
  struct YDN_struct ydn;     /* node description structure                           */
  struct YDB_struct ydb;     /* borehole description structure                       */
  struct YDS_struct yds;     /* inter. fluid (source) description structure          */
  struct YDO_struct ydo;     /* output description structure                         */
  struct YDPE_struct ydpe;   /* property description structure for elements          */
  struct YDPN_struct ydpn;   /* property description structure for nodes             */
  struct YDPJ_struct ydpj;   /* property description structure for joints            */
  struct YDPM_struct ydpm;   /* property description structure for meshing           */
  struct YDFN_struct ydfn;   /* discrete fracture network structure                  */
  struct YDIS_struct ydis;   /* in-situ stress parameters                            */
  struct YDHF_struct ydhf;   /* hydro-fracturing parameters                          */
  struct YDSM_struct ydsm;   /* seismic monitoring parameters                        */
  struct YDR_struct ydr;     /* reference points description structure Y-RC          */
  struct YDSB_struct ydsb;   /* reinforcement description structure Y-RC             */
  struct YDPS_struct ydps;   /* property description structure for steel Y-RC        */
};


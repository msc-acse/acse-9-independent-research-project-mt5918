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
/* File  Yod.c */
#include "Yproto.h"
#include "paraview.h" //! PARAVIEW
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>


static INT  i1num[100];    /* numbers for space saving format     */
static DBL  d1num[100];    /* numbers for space saving format     */
static CHR c1code[500];    /* coded i1para in space saving format */

static void Yod2TRIELS(  /* small strain elastic triangle output */
            nelem,
            fout,
            dcsizc,dcsizs,dcsizv,
            dpeks ,dpela ,dpemu ,dpero ,
            icoutp,iprop ,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,
            i1elpr,i2elto
            )
  INT    nelem;
  FILE   *fout;
  DBL   dcsizc; DBL  dcsizs;  DBL   dcsizv;
  DBL    dpeks; DBL   dpela; DBL   dpemu; DBL    dpero;
  INT   icoutp; INT   iprop;
  DBL  *d1nccx; DBL *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nvcx;
  DBL  *d1nvcy;
  INT *i1elpr; INT **i2elto;
{ DBL voli,volc;
  DBL  B[2][2]; /* left Cauchy-Green strain tensor */
  DBL  D[2][2]; /* rate of deformation (stretching) tensor */
  DBL  E[2][2]; /* strain tensor (small strains) */
  DBL  F[2][2]; /* deformation gradient in global base */
  DBL F0[2][2]; /* initial local base */
  DBL FX[2][2]; /* current local base */
  DBL F0inv[2][2]; /* global base in initial local base */
  DBL FXinv[2][2]; /* global base in current local base */
  DBL  L[2][2]; /* velocity gradient in global base */
  DBL LX[2][2]; /* vel. gradient in current local base = delta x/delta X */
  DBL  T[2][2]; /* Cauchy stress */
  INT ielem;
  INT i,j,k;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { /* evaluate stress state */
      for(i=1;i<3;i++)
      { F0[0][i-1]=d1ncix[(i2elto[i][ielem])]-d1ncix[(i2elto[0][ielem])];
        F0[1][i-1]=d1nciy[(i2elto[i][ielem])]-d1nciy[(i2elto[0][ielem])];
        FX[0][i-1]=d1nccx[(i2elto[i][ielem])]-d1nccx[(i2elto[0][ielem])];
        FX[1][i-1]=d1nccy[(i2elto[i][ielem])]-d1nccy[(i2elto[0][ielem])];
        LX[0][i-1]=d1nvcx[(i2elto[i][ielem])]-d1nvcx[(i2elto[0][ielem])];
        LX[1][i-1]=d1nvcy[(i2elto[i][ielem])]-d1nvcy[(i2elto[0][ielem])];
      }
      YMATINV2(F0,F0inv,voli);
      YMATINV2(FX,FXinv,volc);
      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { F[i][j]=R0;
          L[i][j]=R0;
          for(k=0;k<2;k++)
          { F[i][j]=F[i][j]+FX[i][k]*F0inv[k][j];
            L[i][j]=L[i][j]+LX[i][k]*FXinv[k][j];
      } } }
      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { B[i][j]=R0;
          for(k=0;k<2;k++)
          { B[i][j]=B[i][j]+F[i][k]*F[j][k]; /* left Cauchy-Green strain */
          }
          D[i][j]=RP5*(L[i][j]+L[j][i]);     /* rate of deformation      */
          if(i==j)
          { E[i][j]=RP5*(B[i][j]-R1);        /* small strain             */
          }
          else
          { E[i][j]=RP5*B[i][j];
      } } }
      for(i=0;i<2;i++)     /* Cauchy stress */
      { for(j=0;j<2;j++)
        { T[i][j]=(R2*dpemu*E[i][j])*(voli/volc)+dpeks*D[i][j];
          if(i==j)T[i][j]=T[i][j]+dpela*(volc/voli-voli/volc);
      } }
      /* prepare output */
      i1num[0]=icoutp;
      i1num[1]=15;
      d1num[3]=d1nccx[i2elto[0][ielem]]/dcsizc;
      d1num[4]=d1nccx[i2elto[1][ielem]]/dcsizc;
      d1num[5]=d1nccx[i2elto[2][ielem]]/dcsizc;
      d1num[6]=d1nccy[i2elto[0][ielem]]/dcsizc;
      d1num[7]=d1nccy[i2elto[1][ielem]]/dcsizc;
      d1num[8]=d1nccy[i2elto[2][ielem]]/dcsizc;
      d1num[9]=(d1nvcx[i2elto[0][ielem]]+
                d1nvcx[i2elto[1][ielem]]+
                d1nvcx[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[10]=(d1nvcy[i2elto[0][ielem]]+
                 d1nvcy[i2elto[1][ielem]]+
                 d1nvcy[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[11]=T[0][0]/dcsizs;
      d1num[12]=T[1][1]/dcsizs;
      d1num[13]=T[0][1]/dcsizs;
      d1num[14]=R0;             /* elastic damage */
      for(i=3;i<15;i++)
      { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
      }
      /* translate into INT */
      codeDBLtoINT(d1num,i1num);
      i1num[2]=YTE2TRIELS;
      codeINTtoCHR(c1code,i1num);
      CHRw(fout,c1code); CHRwcr(fout);
} } }
static void Yod2TRISOF(  /* small strain softening triangle output */
            nelem,
            fout,
            dcsizc,dcsizs,dcsizv,
            dpeks ,dpela ,dpemu ,dpero ,
            icoutp,iprop ,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nvcx,
            d1nvcy,d1sdel,i1elpr,i2elto
            )
  INT    nelem;
  FILE   *fout;
  DBL   dcsizc; DBL  dcsizs;  DBL   dcsizv;
  DBL    dpeks; DBL   dpela; DBL   dpemu; DBL    dpero;
  INT   icoutp; INT   iprop;
  DBL  *d1nccx; DBL *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nvcx;
  DBL  *d1nvcy; DBL *d1sdel; INT *i1elpr; INT **i2elto;
{ DBL voli,volc;
  DBL  B[2][2]; /* left Cauchy-Green strain tensor */
  DBL  D[2][2]; /* rate of deformation (stretching) tensor */
  DBL  E[2][2]; /* strain tensor (small strains) */
  DBL  F[2][2]; /* deformation gradient in global base */
  DBL F0[2][2]; /* initial local base */
  DBL FX[2][2]; /* current local base */
  DBL F0inv[2][2]; /* global base in initial local base */
  DBL FXinv[2][2]; /* global base in current local base */
  DBL  L[2][2]; /* velocity gradient in global base */
  DBL LX[2][2]; /* vel. gradient in current local base = delta x/delta X */
  DBL  T[2][2]; /* Cauchy stress */
  INT ielem;
  INT i,j,k;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { /* evaluate stress state */
      for(i=1;i<3;i++)
      { F0[0][i-1]=d1ncix[(i2elto[i][ielem])]-d1ncix[(i2elto[0][ielem])];
        F0[1][i-1]=d1nciy[(i2elto[i][ielem])]-d1nciy[(i2elto[0][ielem])];
        FX[0][i-1]=d1nccx[(i2elto[i][ielem])]-d1nccx[(i2elto[0][ielem])];
        FX[1][i-1]=d1nccy[(i2elto[i][ielem])]-d1nccy[(i2elto[0][ielem])];
        LX[0][i-1]=d1nvcx[(i2elto[i][ielem])]-d1nvcx[(i2elto[0][ielem])];
        LX[1][i-1]=d1nvcy[(i2elto[i][ielem])]-d1nvcy[(i2elto[0][ielem])];
      }
      YMATINV2(F0,F0inv,voli);
      YMATINV2(FX,FXinv,volc);
      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { F[i][j]=R0;
          L[i][j]=R0;
          for(k=0;k<2;k++)
          { F[i][j]=F[i][j]+FX[i][k]*F0inv[k][j];
            L[i][j]=L[i][j]+LX[i][k]*FXinv[k][j];
      } } }
      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { B[i][j]=R0;
          for(k=0;k<2;k++)
          { B[i][j]=B[i][j]+F[i][k]*F[j][k]; /* left Cauchy-Green strain */
          }
          D[i][j]=RP5*(L[i][j]+L[j][i]);     /* rate of deformation      */
          if(i==j)
          { E[i][j]=RP5*(B[i][j]-R1);        /* small strain             */
          }
          else
          { E[i][j]=RP5*B[i][j];
      } } }
      for(i=0;i<2;i++)     /* Cauchy stress */
      { for(j=0;j<2;j++)
        { T[i][j]=(R1-d1sdel[ielem])*
                  (R2*dpemu*E[i][j])*(voli/volc)+dpeks*D[i][j];
          if(i==j)T[i][j]=T[i][j]+dpela*(volc/voli-voli/volc);
      } }
      /* prepare output */
      i1num[0]=icoutp;
      i1num[1]=15;
      d1num[3]=d1nccx[i2elto[0][ielem]]/dcsizc;
      d1num[4]=d1nccx[i2elto[1][ielem]]/dcsizc;
      d1num[5]=d1nccx[i2elto[2][ielem]]/dcsizc;
      d1num[6]=d1nccy[i2elto[0][ielem]]/dcsizc;
      d1num[7]=d1nccy[i2elto[1][ielem]]/dcsizc;
      d1num[8]=d1nccy[i2elto[2][ielem]]/dcsizc;
      d1num[9]=(d1nvcx[i2elto[0][ielem]]+
                d1nvcx[i2elto[1][ielem]]+
                d1nvcx[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[10]=(d1nvcy[i2elto[0][ielem]]+
                 d1nvcy[i2elto[1][ielem]]+
                 d1nvcy[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[11]=T[0][0]/dcsizs;
      d1num[12]=T[1][1]/dcsizs;
      d1num[13]=T[0][1]/dcsizs;
      d1num[14]=d1sdel[ielem];
      for(i=3;i<15;i++)
      { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
      }
      /* translate into INT */
      codeDBLtoINT(d1num,i1num);
      i1num[2]=YTE2TRIELS;
      codeINTtoCHR(c1code,i1num);
      CHRw(fout,c1code); CHRwcr(fout);
} } }
static void Yod2TRIRIG(  /* rigid triangle output */
            nelem,
            fout,
            dcsizc,dcsizv,
            icoutp,iprop ,
            d1nccx,d1nccy,d1nvcx,d1nvcy,i1elpr,
            i2elto
            )
  INT    nelem;
  FILE   *fout;
  DBL   dcsizc; DBL   dcsizv;
  INT   icoutp; INT   iprop;
  DBL  *d1nccx; DBL *d1nccy; DBL *d1nvcx; DBL  *d1nvcy; INT *i1elpr;
  INT **i2elto;
{ INT ielem;
  INT i;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { /* prepare output */
      i1num[0]=icoutp;
      i1num[1]=15;
      d1num[3]=d1nccx[i2elto[0][ielem]]/dcsizc;
      d1num[4]=d1nccx[i2elto[1][ielem]]/dcsizc;
      d1num[5]=d1nccx[i2elto[2][ielem]]/dcsizc;
      d1num[6]=d1nccy[i2elto[0][ielem]]/dcsizc;
      d1num[7]=d1nccy[i2elto[1][ielem]]/dcsizc;
      d1num[8]=d1nccy[i2elto[2][ielem]]/dcsizc;
      d1num[9]=(d1nvcx[i2elto[0][ielem]]+
                d1nvcx[i2elto[1][ielem]]+
                d1nvcx[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[10]=(d1nvcy[i2elto[0][ielem]]+
                 d1nvcy[i2elto[1][ielem]]+
                 d1nvcy[i2elto[2][ielem]])/(R3*dcsizv);
      d1num[11]=R0;
      d1num[12]=R0;
      d1num[13]=R0;
      d1num[14]=R0;
      for(i=3;i<15;i++)
      { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
      }
      /* translate into INT */
      codeDBLtoINT(d1num,i1num);
      i1num[2]=YTE2TRIELS;
      codeINTtoCHR(c1code,i1num);
      CHRw(fout,c1code); CHRwcr(fout);
} } }
static void Yod2JOINTSBROKEN(  /* 2D joint output */
            nelem,
            fout,
            dcsizc,dcsizv,
            icoutp,iprop ,
            d1nccx,d1nccy,d1nvcx,d1nvcy,d1sdel,
            i1elpr,i2elto
            )
  INT    nelem;
  FILE   *fout;
  DBL   dcsizc; DBL   dcsizv;
  INT   icoutp; INT   iprop;
  DBL  *d1nccx; DBL *d1nccy; DBL *d1nvcx; DBL  *d1nvcy; DBL *d1sdel;
  INT  *i1elpr; INT **i2elto;
{ INT ielem;
  INT i;
  INT ipropc;
  for(ielem=0;ielem<nelem;ielem++)
  { ipropc=i1elpr[ielem];
    if(ipropc==iprop)
    { if((i2elto[0][ielem]==i2elto[3][ielem])&&
         (i2elto[1][ielem]==i2elto[2][ielem]))
      { ipropc=ipropc-YIPROPMAX;
    } }
    if((ipropc<0)&&((ipropc+YIPROPMAX)==iprop))
    { /* prepare output */
      i1num[0]=icoutp;
      i1num[1]=12;
      d1num[3]=d1nccx[i2elto[0][ielem]]/dcsizc;
      d1num[4]=d1nccx[i2elto[1][ielem]]/dcsizc;
      d1num[5]=d1nccx[i2elto[2][ielem]]/dcsizc;
      d1num[6]=d1nccx[i2elto[3][ielem]]/dcsizc;
      d1num[7]=d1nccy[i2elto[0][ielem]]/dcsizc;
      d1num[8]=d1nccy[i2elto[1][ielem]]/dcsizc;
      d1num[9]=d1nccy[i2elto[2][ielem]]/dcsizc;
      d1num[10]=d1nccy[i2elto[3][ielem]]/dcsizc;
      d1num[11]=-R1;
      if(d1sdel!=DBL1NULL)d1num[11]=d1sdel[ielem];
      for(i=3;i<15;i++)
      { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
      }
      /* translate into INT */
      codeDBLtoINT(d1num,i1num);
      i1num[2]=YTE2JOINTS;
      codeINTtoCHR(c1code,i1num);
      CHRw(fout,c1code); CHRwcr(fout);
} } }

static void Yod2JOINTSINTACT(  /* 2D joint output - d1peks taken out! */
            nelem,
            fout,
            dcsizc,dcsizv,
            dpeft, dpegf, dpepe,
            icoutp,iprop ,
            d1nccx,d1nccy,d1nvcx,d1nvcy,d1sdel,
            i1elpr,i2elto,d1elfs
            )
  INT    nelem;
  FILE   *fout;
  DBL   dcsizc; DBL   dcsizv;
  DBL   dpeft; DBL   dpegf; DBL   dpepe;
  INT   icoutp; INT   iprop;
  DBL  *d1nccx; DBL *d1nccy; DBL *d1nvcx; DBL  *d1nvcy; DBL *d1sdel;
  INT  *i1elpr; INT **i2elto; DBL *d1elfs;
{ DBL small,o1,o2,s1,s2,op,sp,ot,st;
  DBL e1x,e1y,h;
  INT ielem,i,i0,i1,i2,i3;
  DBL dpefs;

  small=EPSILON;
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { dpefs=d1elfs[ielem];
      i0=i2elto[0][ielem];
      i1=i2elto[1][ielem];
      i2=i2elto[2][ielem];
      i3=i2elto[3][ielem];
      e1x=RP5*(d1nccx[i1]+d1nccx[i2]-d1nccx[i0]-d1nccx[i3]);
      e1y=RP5*(d1nccy[i1]+d1nccy[i2]-d1nccy[i0]-d1nccy[i3]);
      h=SQRT(e1x*e1x+e1y*e1y);
      e1x=e1x/(h+small);
      e1y=e1y/(h+small);
      s1=(d1nccy[i0]-d1nccy[i3])*e1y+(d1nccx[i0]-d1nccx[i3])*e1x;
      s2=(d1nccy[i1]-d1nccy[i2])*e1y+(d1nccx[i1]-d1nccx[i2])*e1x;
      o1=(d1nccy[i0]-d1nccy[i3])*e1x-(d1nccx[i0]-d1nccx[i3])*e1y;
      o2=(d1nccy[i1]-d1nccy[i2])*e1x-(d1nccx[i1]-d1nccx[i2])*e1y;
      op=R2*h*dpeft/dpepe;
      sp=R2*h*dpefs/dpepe;
      ot=MAXIM((R2*op),(R3*dpegf/dpeft));
      st=MAXIM((R2*sp),(R3*dpegf/dpefs));
     /* prepare output */
      if((((o1+o2)/(R2*ot))>RP1)||(((s1+s2)/(R2*st))>RP1))
      { i1num[0]=icoutp;
        i1num[1]=14;
        d1num[3]=d1nccx[i2elto[0][ielem]]/dcsizc;
        d1num[4]=d1nccx[i2elto[1][ielem]]/dcsizc;
        d1num[5]=d1nccx[i2elto[2][ielem]]/dcsizc;
        d1num[6]=d1nccx[i2elto[3][ielem]]/dcsizc;
        d1num[7]=d1nccy[i2elto[0][ielem]]/dcsizc;
        d1num[8]=d1nccy[i2elto[1][ielem]]/dcsizc;
        d1num[9]=d1nccy[i2elto[2][ielem]]/dcsizc;
        d1num[10]=d1nccy[i2elto[3][ielem]]/dcsizc;
        d1num[11]=R1;
        d1num[12]=o1/ot;
        d1num[13]=o2/st;
        for(i=3;i<14;i++)
        { d1num[i]=MAXIM((-R1),MINIM(d1num[i],R1));
        }
        /* translate into INT */
        codeDBLtoINT(d1num,i1num);
        i1num[2]=YTE2JOINTS;
        codeINTtoCHR(c1code,i1num);
        CHRw(fout,c1code); CHRwcr(fout);
} } } }

//! PARAVIEW OUTPUT FUNCTIONS (from paraview.c by AJ)

/* Global stress / strain variables */
static pv_tensor_t* pv_stress;
static pv_tensor_t* pv_principal_stress;
static pv_tensor_t* pv_strain;
static pv_vector_t* pv_principal_stress_direction_1;
static pv_vector_t* pv_principal_stress_direction_2;

/* Global nodal variables */
static pv_vector_t* pv_nodal_force;
static pv_vector_t* pv_nodal_current_coord;
static pv_vector_t* pv_nodal_velocity;
static pv_vector_t* pv_nodal_displacement;

static float *pv_nodal_fluid_pressure;

/* Global element variables */
static int* pv_element_property_id;

/* Global joint variables */
static float *pv_broken_joint_mode_of_failure;
static float *pv_broken_joint_sliding;
static float *pv_broken_joint_opening;
static float *pv_broken_joint_area;
static int* pv_yielded_joint_mode_of_failure;

static pv_vector_t* pv_broken_joint_coordinate;
static pv_vector_t* pv_yielded_joint_coordinate;

static int pv_num_broken_joints;
static int pv_num_yielded_joints;

/* Global rebars variables */

static pv_vector_t* pv_rebar_coordinate;
static float *pv_rebar_area;
static float *pv_rebar_stress;
static float *pv_rebar_force;
static float *pv_rebar_strain;
static int *pv_rebar_activation_flag;

static int pv_num_rebars;

/* Global 1D joint variables */

static pv_vector_t* pv_1D_joint_coordinate;
static float *pv_1D_joint_normal_stress;
static float *pv_1D_joint_tangential_stress;
static float *pv_1D_joint_normal_strain;
static float *pv_1D_joint_tangential_strain;
static int   *pv_1D_joint_state;
//static pv_vector_t* pv_1D_joint_force;

static int pv_num_unbroken_1D_joints;
static int pv_num_1D_joints;

/* Connectivity tables */
static int* pv_element_to_nodes;
static int* pv_joint_to_nodes;

/* Other important variables to hold onto */
static int pv_num_elements;
static int pv_num_joints;
static const char* pv_file_prefix;
static double pv_realtime_per_timestep;

/* De-Munjiza Mangler lookup table  - Map munjiza mangling to sanity */
static int* pv_munjiza_to_sane_table;

/* Count of Munjiza joints we have seen so far */
static int pv_munjiza_joints_seen;

/* Count of Munjiza elements we have seen so far */
static int pv_munjiza_elements_seen;

/* Number of dead nodes in the model, paraview_node_id = node_id - num_dead_nodes */
static int pv_num_dead_nodes;

/* Acoustic Information */
static int pv_num_acoustic_events;
static pv_vector_t* pv_acoustic_coordinate;
static float *pv_acoustic_energy;
static float *pv_acoustic_magnitude;
static float *pv_acoustic_fracture_energy;
static float *pv_acoustic_event_mode;

/* Principal directions information */
static pv_vector_t* pv_element_centroid_coordinate;
static pv_vector_t* pv_principal_stress_direction_1_andrea;
static pv_vector_t* pv_principal_stress_direction_2_andrea;




void pv_register_munjiza_joint(int munjiza_joint)
{
    pv_munjiza_to_sane_table[munjiza_joint] = pv_munjiza_joints_seen;
    pv_munjiza_joints_seen++;
}

void pv_register_munjiza_element(int munjiza_element)
{
    pv_munjiza_to_sane_table[munjiza_element] = pv_munjiza_elements_seen;
    pv_munjiza_elements_seen++;
}

int translate_munjiza_element_to_sane(int munjiza_element)
{
    return pv_munjiza_to_sane_table[munjiza_element];
}

int translate_munjiza_joint_to_sane(int munjiza_joint)
{
    return pv_munjiza_to_sane_table[munjiza_joint];
}

void pv_init(int num_elements, int upper_bound_num_joints, int num_rebars, int upper_bound_num_1D_joints, int num_dead_nodes, int upper_bound_acoustic_events, const char* file_prefix, double realtime_per_timestep)
{

    /* Munjiza Translator */
    //pv_munjiza_joints_seen = 0;
    pv_munjiza_elements_seen = 0;
    pv_num_dead_nodes = num_dead_nodes;

    /* Remember variables */
    pv_num_elements = num_elements;
    //pv_num_joints = num_joints;
    pv_file_prefix = file_prefix;
    pv_realtime_per_timestep = realtime_per_timestep;

    /* Allocate element arrays */
    pv_stress                        = malloc( num_elements * sizeof(pv_tensor_t) );
    if(pv_stress == NULL) { exit(1); }

    pv_principal_stress              = malloc( num_elements * sizeof(pv_tensor_t) );
    if(pv_principal_stress == NULL) { exit(1); }

    pv_strain                        = malloc( num_elements * sizeof(pv_tensor_t) );
    if(pv_strain == NULL) { exit(1); }

    //pv_principal_stress_direction_1  = malloc( num_elements * sizeof(pv_vector_t) );
    //if(pv_principal_stress_direction_1 == NULL) { exit(1); }

    //pv_principal_stress_direction_2  = malloc( num_elements * sizeof(pv_vector_t) );
    //if(pv_principal_stress_direction_2 == NULL) { exit(1); }

     
    /* Allocate node arrays */
    pv_nodal_force                   = malloc( 3 * num_elements * sizeof(pv_vector_t) );
    if(pv_nodal_force == NULL) { exit(1); }

    pv_nodal_current_coord           = malloc( 3 * num_elements * sizeof(pv_vector_t) );
    if(pv_nodal_current_coord == NULL) { exit(1); }

    pv_nodal_velocity                = malloc( 3 * num_elements * sizeof(pv_vector_t) );
    if(pv_nodal_velocity == NULL) { exit(1); }

    pv_nodal_displacement            = malloc( 3 * num_elements * sizeof(pv_vector_t) );
    if(pv_nodal_displacement == NULL) { exit(1); }
    
    pv_nodal_fluid_pressure          = malloc( 3 * num_elements * sizeof(float) );
    if(pv_nodal_fluid_pressure == NULL) { exit(1); }
     
    /* Allocate element properties */
    pv_element_property_id           = malloc( num_elements * sizeof(int) );
    if(pv_element_property_id == NULL) { exit(1); }

     
    /* Allocate joint variables */
    pv_broken_joint_mode_of_failure  = malloc( upper_bound_num_joints * sizeof(float) );
    if(pv_broken_joint_mode_of_failure == NULL) { exit(1); }

    pv_broken_joint_sliding  = malloc( upper_bound_num_joints * sizeof(float) );
    if(pv_broken_joint_sliding == NULL) { exit(1); }
    
    pv_broken_joint_opening  = malloc( upper_bound_num_joints * sizeof(float) );
    if(pv_broken_joint_opening == NULL) { exit(1); }
    
    pv_broken_joint_area  = malloc( upper_bound_num_joints * sizeof(float) );
    if(pv_broken_joint_area == NULL) { exit(1); }
   
    pv_yielded_joint_mode_of_failure  = malloc( upper_bound_num_joints * sizeof(int) );
    if(pv_yielded_joint_mode_of_failure == NULL) { exit(1); }

    pv_broken_joint_coordinate  = malloc( 4 * upper_bound_num_joints * sizeof(pv_vector_t) );
    if(pv_broken_joint_coordinate == NULL) { exit(1); }

    pv_yielded_joint_coordinate  = malloc( 4 * upper_bound_num_joints * sizeof(pv_vector_t) );
    if(pv_yielded_joint_coordinate == NULL) { exit(1); }

    pv_num_broken_joints = 0;
    pv_num_yielded_joints = 0;
    
    /* Allocate rebars variables */
    pv_rebar_coordinate  = malloc(2 * num_rebars * sizeof(pv_vector_t) );
    if(pv_rebar_coordinate == NULL) { exit(1); }
    
    pv_rebar_area  = malloc( num_rebars * sizeof(float) );
    if(pv_rebar_area == NULL) { exit(1); }
    
    pv_rebar_stress  = malloc( num_rebars * sizeof(float) );
    if(pv_rebar_stress == NULL) { exit(1); }
    
    pv_rebar_force  = malloc( num_rebars * sizeof(float) );
    if(pv_rebar_force == NULL) { exit(1); }
    
    pv_rebar_strain  = malloc( num_rebars * sizeof(float) );
    if(pv_rebar_strain == NULL) { exit(1); }
    
    pv_rebar_activation_flag  = malloc( num_rebars * sizeof(int) );
    if(pv_rebar_activation_flag == NULL) { exit(1); }
    
    //pv_num_rebars = 0;
    
     /* Allocate 1D joints variables */
    pv_1D_joint_coordinate  = malloc(2 * upper_bound_num_1D_joints * sizeof(pv_vector_t) );
    if(pv_1D_joint_coordinate == NULL) { exit(1); }
    
    pv_1D_joint_normal_stress  = malloc(upper_bound_num_1D_joints * sizeof(float) );
    if(pv_1D_joint_normal_stress == NULL) { exit(1); }
    
    pv_1D_joint_tangential_stress  = malloc(upper_bound_num_1D_joints * sizeof(float) );
    if(pv_1D_joint_tangential_stress == NULL) { exit(1); }
    
    pv_1D_joint_normal_strain  = malloc(upper_bound_num_1D_joints * sizeof(float) );
    if(pv_1D_joint_normal_strain == NULL) { exit(1); }
    
    pv_1D_joint_tangential_strain  = malloc(upper_bound_num_1D_joints * sizeof(float) );
    if(pv_1D_joint_tangential_strain == NULL) { exit(1); }
   
    pv_1D_joint_state  = malloc(upper_bound_num_1D_joints * sizeof(int) );
    if(pv_1D_joint_state == NULL) { exit(1); }
    
    //pv_1D_joint_force = malloc(4 * upper_bound_num_1D_joints * sizeof(pv_vector_t) );
    //if(pv_1D_joint_force == NULL) { exit(1); }
   
    //pv_num_1D_joints = 0;
    pv_num_unbroken_1D_joints = 0;

    /* Allocate acoustic emission data */
    pv_acoustic_coordinate  = malloc( upper_bound_acoustic_events * sizeof(pv_vector_t) );
    if(pv_acoustic_coordinate == NULL) { exit(1); }

    pv_acoustic_energy  = malloc( upper_bound_acoustic_events * sizeof(float) );
    if(pv_acoustic_energy == NULL) { exit(1); }

    pv_acoustic_magnitude  = malloc( upper_bound_acoustic_events * sizeof(float) );
    if(pv_acoustic_magnitude == NULL) { exit(1); }

    pv_acoustic_fracture_energy  = malloc( upper_bound_acoustic_events * sizeof(float) );
    if(pv_acoustic_fracture_energy == NULL) { exit(1); }

    pv_acoustic_event_mode  = malloc( upper_bound_acoustic_events * sizeof(int) );
    if(pv_acoustic_event_mode == NULL) { exit(1); }

    pv_num_acoustic_events = 0;
    
    /* Allocate principal direction data */
    
    pv_element_centroid_coordinate  = malloc( num_elements * sizeof(pv_vector_t) );
    if(pv_element_centroid_coordinate == NULL) { exit(1); }
    
    pv_principal_stress_direction_1_andrea  = malloc( num_elements * sizeof(pv_vector_t) );
    if(pv_principal_stress_direction_1_andrea == NULL) { exit(1); }

    pv_principal_stress_direction_2_andrea  = malloc( num_elements * sizeof(pv_vector_t) );
    if(pv_principal_stress_direction_2_andrea == NULL) { exit(1); }
    
     
    /* Allocate connectivity table */
    pv_element_to_nodes              = malloc( num_elements * 3 * sizeof(int) );
    if(pv_element_to_nodes == NULL) { exit(1); }

    //pv_joint_to_nodes                = malloc( num_joints * 4 * sizeof(int) );
    //if(pv_joint_to_nodes == NULL) { exit(1); }

    /* Allocate Munjiza translator */
    pv_munjiza_to_sane_table         = malloc( (upper_bound_num_joints + num_elements)*sizeof(int) );
    if(pv_munjiza_to_sane_table == NULL) { exit(1); }
}

void pv_free()
{
    free(pv_stress);
    free(pv_principal_stress);
    free(pv_strain);
    //free(pv_principal_stress_direction_1);
    //free(pv_principal_stress_direction_2);
    free(pv_nodal_force);
    free(pv_nodal_current_coord);
    free(pv_nodal_velocity);
    free(pv_nodal_displacement);
    free(pv_nodal_fluid_pressure);
    free(pv_element_property_id);
    free(pv_broken_joint_mode_of_failure);
    free(pv_broken_joint_sliding);
    free(pv_broken_joint_opening);
    free(pv_broken_joint_area);
    free(pv_yielded_joint_mode_of_failure);
    free(pv_broken_joint_coordinate);
    free(pv_yielded_joint_coordinate);
    free(pv_element_to_nodes);
    free(pv_munjiza_to_sane_table);
    free(pv_acoustic_coordinate);
    free(pv_acoustic_energy);
    free(pv_acoustic_magnitude);
    free(pv_acoustic_fracture_energy);
    free(pv_acoustic_event_mode);
    free(pv_element_centroid_coordinate);
    free(pv_principal_stress_direction_1_andrea);
    free(pv_principal_stress_direction_2_andrea);
    free(pv_rebar_coordinate);
    free(pv_rebar_area);
    free(pv_rebar_stress);
    free(pv_rebar_force);
    free(pv_rebar_strain);
    free(pv_rebar_activation_flag);
    free(pv_1D_joint_coordinate);
    free(pv_1D_joint_normal_stress);
    free(pv_1D_joint_tangential_stress);
    free(pv_1D_joint_normal_strain);
    free(pv_1D_joint_tangential_strain);
    free(pv_1D_joint_state);
    //free(pv_1D_joint_force);
}

void pv_write_stress(int munjiza_element_id, pv_tensor_t tensor)
{
    int element_id;
    element_id =  translate_munjiza_element_to_sane(munjiza_element_id);

    if(element_id >= pv_num_elements) { printf("Element ID %d out of range", element_id); exit(1); }
    pv_stress[element_id] = tensor;
}

void pv_write_principal_stress(int munjiza_element_id, pv_tensor_t tensor)
{
    int element_id;
    element_id =  translate_munjiza_element_to_sane(munjiza_element_id);

    if(element_id >= pv_num_elements) { printf("Element ID %d out of range", element_id); exit(1); }
    pv_principal_stress[element_id] = tensor;
}

void pv_write_strain(int munjiza_element_id, pv_tensor_t tensor)
{
    int element_id;
    element_id =  translate_munjiza_element_to_sane(munjiza_element_id);

    if(element_id >= pv_num_elements) { printf("Element ID %d out of range", element_id); exit(1); }
    pv_strain[element_id] = tensor;
}

/*void pv_write_principal_stress_direction_1(int munjiza_element_id, pv_vector_t vector)
{
    int element_id;
    element_id =  translate_munjiza_element_to_sane(munjiza_element_id);

    if(element_id >= pv_num_elements) { printf("Element ID %d out of range", element_id); exit(1); }
    pv_principal_stress_direction_1[element_id] = vector;
}*/

/*void pv_write_principal_stress_direction_2(int munjiza_element_id, pv_vector_t vector)
{
    int element_id;
    element_id =  translate_munjiza_element_to_sane(munjiza_element_id);

    if(element_id >= pv_num_elements) { printf("Element ID %d out of range", element_id); exit(1); }
    pv_principal_stress_direction_2[element_id] = vector;
}*/

void pv_write_nodal_force(int munjiza_node_id, pv_vector_t vector)
{
    int node_id;
    node_id = munjiza_node_id - pv_num_dead_nodes; 


    if(node_id >= 3*pv_num_elements) { printf(" pv_write_nodal_force - Node ID %d out of range", node_id); exit(1); }
    pv_nodal_force[node_id] = vector;
}

void pv_write_nodal_current_coord(int munjiza_node_id, pv_vector_t vector)
{
    int node_id;
    node_id = munjiza_node_id - pv_num_dead_nodes; 

    if(node_id >= 3*pv_num_elements) { printf(" pv_write_nodal_current_coord - Node ID %d out of range", node_id); exit(1); }
    pv_nodal_current_coord[node_id] = vector;
}

void pv_write_nodal_velocity(int munjiza_node_id, pv_vector_t vector)
{
    int node_id;
    node_id = munjiza_node_id - pv_num_dead_nodes; 

    if(node_id >= 3*pv_num_elements) { printf(" pv_write_nodal_velocity - Node ID %d out of range", node_id); exit(1); }
    pv_nodal_velocity[node_id] = vector;
}

void pv_write_nodal_displacement(int munjiza_node_id, pv_vector_t vector)
{
    int node_id;
    node_id = munjiza_node_id - pv_num_dead_nodes; 

    if(node_id >= 3*pv_num_elements) { printf(" pv_write_nodal_displacement - Node ID %d out of range", node_id); exit(1); }
    pv_nodal_displacement[node_id] = vector;
}

void pv_write_nodal_fluid_pressure(int munjiza_node_id, float fluid_pressure)
{
    int node_id;
    node_id = munjiza_node_id - pv_num_dead_nodes; 

    if(node_id >= 3*pv_num_elements) { printf(" pv_write_nodal_fluid_pressure - Node ID %d out of range", node_id); exit(1); }
    pv_nodal_fluid_pressure[node_id] = fluid_pressure;
}

void pv_write_element_property_id(int munjiza_element_id, int property_id)
{
    int element_id;
    element_id =  translate_munjiza_element_to_sane(munjiza_element_id);

    if(element_id >= pv_num_elements) { printf("Element ID %d out of range", element_id); exit(1); }
    pv_element_property_id[element_id] = property_id;
}

/// void pv_write_joint_state(int munjiza_joint_id, int joint_state)
/// {
///     int joint_id;
///     joint_id = translate_munjiza_joint_to_sane(munjiza_joint_id);
/// 
///     pv_joint_state[joint_id] = joint_state;
/// }
/// 
/// void pv_write_joint_mode_of_failure(int munjiza_joint_id, int mode_of_failure)
/// {
///     int joint_id;
///     joint_id = translate_munjiza_joint_to_sane(munjiza_joint_id);
/// 
///     pv_joint_mode_of_failure[joint_id] = mode_of_failure;
/// }
/// 
/// void pv_write_joint_breakage_kinetic_energy(int munjiza_joint_id, float energy)
/// {
///     int joint_id;
///     joint_id = translate_munjiza_joint_to_sane(munjiza_joint_id);
/// 
///     pv_joint_breakage_kinetic_energy[joint_id] = energy;
/// }

void pv_write_element_connectivity(int munjiza_element_id, int node1, int node2, int node3)
{
    int element_id;
    element_id =  translate_munjiza_element_to_sane(munjiza_element_id);

    pv_element_to_nodes[3*element_id    ] = node1 - pv_num_dead_nodes;
    pv_element_to_nodes[3*element_id + 1] = node2 - pv_num_dead_nodes;
    pv_element_to_nodes[3*element_id + 2] = node3 - pv_num_dead_nodes;
}

//void pv_write_joint_connectivity(int munjiza_joint_id, int node1, int node2, int node3, int node4)
//{
//    int joint_id;
//    joint_id = translate_munjiza_joint_to_sane(munjiza_joint_id);
//
//    pv_joint_to_nodes[4*joint_id    ] = node1 - pv_num_dead_nodes;
//    pv_joint_to_nodes[4*joint_id + 1] = node2 - pv_num_dead_nodes;
//    pv_joint_to_nodes[4*joint_id + 2] = node3 - pv_num_dead_nodes;
//    pv_joint_to_nodes[4*joint_id + 3] = node4 - pv_num_dead_nodes;
//}

void pv_write_broken_joint_mode_of_failure(int broken_joint_id, float mode_of_failure)
{
     if(broken_joint_id >= pv_num_broken_joints) { printf(" pv_write_broken_joint_mode_of_failure - Out of range."); exit(1); }
     pv_broken_joint_mode_of_failure[broken_joint_id] = mode_of_failure;
}

void pv_write_broken_joint_sliding(int broken_joint_id, float sliding)
{
     if(broken_joint_id >= pv_num_broken_joints) { printf(" pv_write_broken_joint_sliding - Out of range."); exit(1); }
     pv_broken_joint_sliding[broken_joint_id] = sliding;
}

void pv_write_broken_joint_opening(int broken_joint_id, float opening)
{
     if(broken_joint_id >= pv_num_broken_joints) { printf(" pv_write_broken_joint_opening - Out of range."); exit(1); }
     pv_broken_joint_opening[broken_joint_id] = opening;
}

void pv_write_broken_joint_area(int broken_joint_id, float area)
{
     if(broken_joint_id >= pv_num_broken_joints) { printf(" pv_write_broken_joint_area - Out of range."); exit(1); }
     pv_broken_joint_area[broken_joint_id] = area;
}

void pv_write_yielded_joint_mode_of_failure(int yielded_joint_id, int mode_of_failure)
{ 
     if(yielded_joint_id >= pv_num_yielded_joints) { printf(" pv_write_yielded_joint_mode_of_failure - Out of range."); exit(1); }
     pv_yielded_joint_mode_of_failure[yielded_joint_id] = mode_of_failure;
}

void pv_write_broken_joint_coordinates(int broken_joint_id, pv_vector_t pt0,  pv_vector_t pt1,  pv_vector_t pt2,  pv_vector_t pt3)
{
     if(broken_joint_id >= pv_num_broken_joints) { printf(" pv_write_broken_joint_coordinates - Out of range."); exit(1); }
     pv_broken_joint_coordinate[4*broken_joint_id    ] = pt0;
     pv_broken_joint_coordinate[4*broken_joint_id + 1] = pt1;
     pv_broken_joint_coordinate[4*broken_joint_id + 2] = pt2;
     pv_broken_joint_coordinate[4*broken_joint_id + 3] = pt3;
}

void pv_write_yielded_joint_coordinates(int yielded_joint_id, pv_vector_t pt0,  pv_vector_t pt1,  pv_vector_t pt2,  pv_vector_t pt3)
{
     if(yielded_joint_id >= pv_num_yielded_joints) { printf(" pv_write_yielded_joint_coordinates - Out of range."); exit(1); }
     pv_yielded_joint_coordinate[4*yielded_joint_id    ] = pt0;
     pv_yielded_joint_coordinate[4*yielded_joint_id + 1] = pt1;
     pv_yielded_joint_coordinate[4*yielded_joint_id + 2] = pt2;
     pv_yielded_joint_coordinate[4*yielded_joint_id + 3] = pt3;
}

void pv_write_rebar_coordinates(int rebar_id, pv_vector_t pt0,  pv_vector_t pt1)
{
     if(rebar_id >= pv_num_rebars) { printf(" pv_write_rebar_coordinates - Out of range."); exit(1); }
     pv_rebar_coordinate[2*rebar_id    ] = pt0;
     pv_rebar_coordinate[2*rebar_id + 1] = pt1;
}

void pv_write_rebar_area(int rebar_id, float area)
{ 
     if(rebar_id >= pv_num_rebars) { printf(" pv_write_rebar_area - Out of range."); exit(1); }
     pv_rebar_area[rebar_id] = area;
}

void pv_write_rebar_stress(int rebar_id, float stress)
{ 
     if(rebar_id >= pv_num_rebars) { printf(" pv_write_rebar_stress - Out of range."); exit(1); }
     pv_rebar_stress[rebar_id] = stress;
}

void pv_write_rebar_force(int rebar_id, float force)
{ 
     if(rebar_id >= pv_num_rebars) { printf(" pv_write_rebar_force - Out of range."); exit(1); }
     pv_rebar_force[rebar_id] = force;
}

void pv_write_rebar_strain(int rebar_id, float strain)
{ 
     if(rebar_id >= pv_num_rebars) { printf(" pv_write_rebar_strain - Out of range."); exit(1); }
     pv_rebar_strain[rebar_id] = strain;
}

void pv_write_rebar_activation_flag(int rebar_id, int rebar_activation_flag)
{ 
     if(rebar_id >= pv_num_rebars) { printf(" pv_write_rebar_activation_flag - Out of range."); exit(1); }
     pv_rebar_activation_flag[rebar_id] = rebar_activation_flag;
}

void pv_write_1D_joint_coordinates(int oneD_joint_id, pv_vector_t pt0,  pv_vector_t pt1)
{
     if(oneD_joint_id >= pv_num_unbroken_1D_joints) { printf(" pv_write_1D_joint_coordinates - Out of range."); exit(1); }
     pv_1D_joint_coordinate[2*oneD_joint_id    ] = pt0;
     pv_1D_joint_coordinate[2*oneD_joint_id + 1] = pt1;
}

void pv_write_1D_joint_normal_stress(int oneD_joint_id, float stress)
{ 
     if(oneD_joint_id >= pv_num_unbroken_1D_joints) { printf(" pv_write_1D_joint_normal_stress - Out of range."); exit(1); }
     pv_1D_joint_normal_stress[oneD_joint_id] = stress;
}

void pv_write_1D_joint_tangential_stress(int oneD_joint_id, float stress)
{ 
     if(oneD_joint_id >= pv_num_unbroken_1D_joints) { printf(" pv_write_1D_joint_tangential_stress - Out of range."); exit(1); }
     pv_1D_joint_tangential_stress[oneD_joint_id] = stress;
}

void pv_write_1D_joint_normal_strain(int oneD_joint_id, float strain)
{ 
     if(oneD_joint_id >= pv_num_unbroken_1D_joints) { printf(" pv_write_1D_joint_normal_strain - Out of range."); exit(1); }
     pv_1D_joint_normal_strain[oneD_joint_id] = strain;
}

void pv_write_1D_joint_tangential_strain(int oneD_joint_id, float strain)
{ 
     if(oneD_joint_id >= pv_num_unbroken_1D_joints) { printf(" pv_write_1D_joint_tangential_strain - Out of range."); exit(1); }
     pv_1D_joint_tangential_strain[oneD_joint_id] = strain;
}

void pv_write_1D_joint_state(int oneD_joint_id, int state)
{ 
     if(oneD_joint_id >= pv_num_unbroken_1D_joints) { printf(" pv_write_1D_joint_state - Out of range."); exit(1); }
     pv_1D_joint_state[oneD_joint_id] = state;
}

/*void pv_write_1D_joint_force(int oneD_joint_id, pv_vector_t vector)
{

    if(oneD_joint_id >= pv_num_unbroken_1D_joints) { printf(" pv_write_1D_joint_force - Out of range"); exit(1); }
    pv_1D_joint_force[oneD_joint_id] = vector;
}*/

void pv_write_acoustic_coordinate(int acoustic_node_id, pv_vector_t coord)
{
    pv_acoustic_coordinate[acoustic_node_id] = coord;
}

void pv_write_acoustic_energy(int acoustic_node_id, float energy)
{
    pv_acoustic_energy[acoustic_node_id] = energy;
}

void pv_write_acoustic_magnitude(int acoustic_node_id, float magnitude)
{
    pv_acoustic_magnitude[acoustic_node_id] = magnitude;
}

void pv_write_acoustic_fracture_energy(int acoustic_node_id, float energy)
{
    pv_acoustic_fracture_energy[acoustic_node_id] = energy;
}

void pv_write_acoustic_event_mode(int acoustic_node_id, float mode)
{
    pv_acoustic_event_mode[acoustic_node_id] = mode;
}

void pv_set_num_broken_joints(int count)
{
    pv_num_broken_joints = count;
}

void pv_set_num_yielded_joints(int count)
{
    pv_num_yielded_joints = count;
}

void pv_set_num_unbroken_1D_joints(int count)
{
    pv_num_unbroken_1D_joints = count;
}

void pv_set_num_acoustic_events(int count)
{
    pv_num_acoustic_events = count;
}

void pv_set_num_rebars(int count)
{
    pv_num_rebars = count;
}

/* Principal direction */
void pv_write_element_centroid_coordinate(int munjiza_element_id, pv_vector_t coord)
{
    int element_id;
    element_id =  translate_munjiza_element_to_sane(munjiza_element_id);
    
    if(element_id >= pv_num_elements) { printf("Element ID %d out of range", element_id); exit(1); }
    pv_element_centroid_coordinate[element_id] = coord;
}
void pv_write_principal_stress_direction_1_andrea(int munjiza_element_id, pv_vector_t vector)
{
    int element_id; 
    element_id =  translate_munjiza_element_to_sane(munjiza_element_id);

    if(element_id >= pv_num_elements) { printf("Principal direction point ID %d out of range", element_id); exit(1); }
    pv_principal_stress_direction_1_andrea[element_id] = vector;
}
void pv_write_principal_stress_direction_2_andrea(int munjiza_element_id, pv_vector_t vector)
{
    int element_id;
    element_id =  translate_munjiza_element_to_sane(munjiza_element_id);

    if(element_id >= pv_num_elements) { printf("Principal direction point ID %d out of range", element_id); exit(1); }
    pv_principal_stress_direction_2_andrea[element_id] = vector;
}


/* Write output file for the timestep.  Called when the timestep is done. */
void pv_write_timestep(int timestep)
{
    FILE* vtp_file;
    FILE* broken_joint_vtp_file;
    FILE* yielded_joint_vtp_file;
    FILE* acoustic_vtp_file;
    FILE* pvd_file;
    FILE* principal_direction_vtp_file;
    FILE* rebar_vtp_file;
    FILE* oneD_joint_vtp_file;
    
    char tmp_filename[1024];

    int i, j, k;

    snprintf(tmp_filename, 1024, "%s_%d.vtp", pv_file_prefix, timestep);
//printf("DEBUG: %s", tmp_filename);

    /* Open the .vtp file */
    vtp_file = fopen(tmp_filename, "w");
    if(vtp_file == NULL) { perror("Could not open VTP output file"); exit(1); }

    fprintf(vtp_file, "<?xml version='1.0'?>\n<VTKFile type='PolyData' version='0.1' byte_order='LittleEndian'>");

   fprintf(vtp_file, "<PolyData>\n<Piece NumberOfPoints='%d' NumberOfVerts='0' NumberOfLines='0' NumberOfStrips='0' NumberOfPolys='%d'>\n", 3*pv_num_elements, pv_num_elements);

   /* Write out the points */
   fprintf(vtp_file, "<Points>\n<DataArray type='Float32' NumberOfComponents='3' format='ascii'>\n");

   for(i = 0; i < 3*pv_num_elements; ++i)
   {
       fprintf(vtp_file, "%.6f %.6f %.6f ", (double)pv_nodal_current_coord[i].x, (double)pv_nodal_current_coord[i].y, (double)pv_nodal_current_coord[i].z );
   }

   fprintf(vtp_file, "\n</DataArray>\n</Points>\n");

   /* Write out the data associated with the points */
   fprintf(vtp_file, "<PointData Vectors='velocity force displacement'>\n");

       /* Write out velocity */
   fprintf(vtp_file, "<DataArray type='Float32' Name='velocity' NumberOfComponents='3' format='ascii'>\n");

   for(i = 0; i < 3*pv_num_elements; ++i)
   {
       fprintf(vtp_file, "%.6f %.6f %.6f ", (double)pv_nodal_velocity[i].x, (double)pv_nodal_velocity[i].y, (double)pv_nodal_velocity[i].z );
   }
   fprintf(vtp_file, "\n</DataArray>\n");

       /* Write out force */
   fprintf(vtp_file, "<DataArray type='Float32' Name='force' NumberOfComponents='3' format='ascii'>\n");

   for(i = 0; i < 3*pv_num_elements; ++i)
   {
       fprintf(vtp_file, "%.6f %.6f %.6f ", (double)pv_nodal_force[i].x, (double)pv_nodal_force[i].y, (double)pv_nodal_force[i].z );
   }
   fprintf(vtp_file, "\n</DataArray>\n");

       /* Write out displacement */
   fprintf(vtp_file, "<DataArray type='Float32' Name='displacement' NumberOfComponents='3' format='ascii'>\n");

   for(i = 0; i < 3*pv_num_elements; ++i)
   {
       fprintf(vtp_file, "%.6f %.6f %.6f ", (double)pv_nodal_displacement[i].x, (double)pv_nodal_displacement[i].y, (double)pv_nodal_displacement[i].z );
   }
   fprintf(vtp_file, "\n</DataArray>\n");

      /* Write out nodal fluid pressure */
   fprintf(vtp_file, "<DataArray type='Float32' Name='fluid pressure' NumberOfComponents='1' format='ascii'>\n");

   for(i = 0; i < 3*pv_num_elements; ++i)
   {
       fprintf(vtp_file, "%.6f ", (double)pv_nodal_fluid_pressure[i]);
   }
   fprintf(vtp_file, "\n</DataArray>\n");

   
   /* Done writing point data */
   fprintf(vtp_file, "</PointData>\n");


   /* Write out the element data */
   fprintf(vtp_file, "<CellData Tensors='stress strain principal_stress' Vectors='principal_stress_direction_1 principal_stress_direction_2' Scalars='property_id'>\n");

       /* Write out stress */
   fprintf(vtp_file, "<DataArray type='Float32' Name='stress' NumberOfComponents='9' format='ascii'>\n");

   for(i = 0; i < pv_num_elements; ++i)
   {
       fprintf(vtp_file, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ", (double)pv_stress[i].a11,  (double)pv_stress[i].a12,  (double)pv_stress[i].a13,  (double)pv_stress[i].a21,  (double)pv_stress[i].a22,  (double)pv_stress[i].a23,  (double)pv_stress[i].a31,  (double)pv_stress[i].a32,  (double)pv_stress[i].a33);  
   }
   fprintf(vtp_file, "\n</DataArray>\n");

       /* Write out strain */
   fprintf(vtp_file, "<DataArray type='Float32' Name='strain' NumberOfComponents='9' format='ascii'>\n");

   for(i = 0; i < pv_num_elements; ++i)
   {
       fprintf(vtp_file, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ", (double)pv_strain[i].a11,  (double)pv_strain[i].a12,  (double)pv_strain[i].a13,  (double)pv_strain[i].a21,  (double)pv_strain[i].a22,  (double)pv_strain[i].a23,  (double)pv_strain[i].a31,  (double)pv_strain[i].a32,  (double)pv_strain[i].a33);  
   }
   fprintf(vtp_file, "\n</DataArray>\n");

       /* Write out principal_stress */
   fprintf(vtp_file, "<DataArray type='Float32' Name='principal_stress' NumberOfComponents='9' format='ascii'>\n");

   for(i = 0; i < pv_num_elements; ++i)
   {
       fprintf(vtp_file, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ", (double)pv_principal_stress[i].a11,  (double)pv_principal_stress[i].a12,  (double)pv_principal_stress[i].a13,  (double)pv_principal_stress[i].a21,  (double)pv_principal_stress[i].a22,  (double)pv_principal_stress[i].a23,  (double)pv_principal_stress[i].a31,  (double)pv_principal_stress[i].a32,  (double)pv_principal_stress[i].a33);  
   }
   fprintf(vtp_file, "\n</DataArray>\n");

   //    /* Write out principal_stress_direction_1 */
   //fprintf(vtp_file, "<DataArray type='Float32' Name='principal_stress_direction_1' NumberOfComponents='3' format='ascii'>\n");

   //for(i = 0; i < pv_num_elements; ++i)
   //{
   //    fprintf(vtp_file, "%.6f %.6f %.6f ", (double)pv_principal_stress_direction_1[i].x, (double)pv_principal_stress_direction_1[i].y, (double)pv_principal_stress_direction_1[i].z );
   //}
   //fprintf(vtp_file, "\n</DataArray>\n");

   //    /* Write out principal_stress_direction_2 */
   //fprintf(vtp_file, "<DataArray type='Float32' Name='principal_stress_direction_2' NumberOfComponents='3' format='ascii'>\n");

   //for(i = 0; i < pv_num_elements; ++i)
   //{
   //    fprintf(vtp_file, "%.6f %.6f %.6f ", (double)pv_principal_stress_direction_2[i].x, (double)pv_principal_stress_direction_2[i].y, (double)pv_principal_stress_direction_2[i].z );
   //}
   //fprintf(vtp_file, "\n</DataArray>\n");

       /* Write out property_id */
   fprintf(vtp_file, "<DataArray type='Int32' Name='property_id' NumberOfComponents='1' format='ascii'>\n");

   for(i = 0; i < pv_num_elements; ++i)
   {
       fprintf(vtp_file, "%d ", pv_element_property_id[i]);
   }
   fprintf(vtp_file, "\n</DataArray>\n");

   /* Done writing element data */
   fprintf(vtp_file, "\n</CellData>\n");

   /* Write out the polygons */
   fprintf(vtp_file, "<Polys>\n");

   /* Write connectivity table */
   fprintf(vtp_file, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");

   for(i = 0; i < pv_num_elements; ++i)
   {
       fprintf(vtp_file, "%d %d %d ", pv_element_to_nodes[3*i],  pv_element_to_nodes[3*i + 1], pv_element_to_nodes[3*i + 2]);  
   }

   fprintf(vtp_file, "\n</DataArray>\n");

   /* Write offsets table */
   fprintf(vtp_file, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");

   for(i = 3; i <= 3*pv_num_elements; i = i + 3)
   {
       fprintf(vtp_file, "%d ", i); 
   }

   fprintf(vtp_file, "\n</DataArray>\n");


   /* Done writing polygons */
   fprintf(vtp_file, "</Polys>\n");

   /* Done writing file */
   fprintf(vtp_file, "</Piece>\n</PolyData>\n</VTKFile>");
  
    if(fflush(vtp_file)) { perror("Could not fflush VTP output file"); exit(1); }
    if(fclose(vtp_file)) { perror("Could not close VTP output file"); exit(1); }

    /*
     *
     *
     *
     *
     */

    /* Open the broken joint .vtp file */
    snprintf(tmp_filename, 1024, "%s_broken_joints_%d.vtp", pv_file_prefix, timestep);

    /* Open the .vtp file */
    broken_joint_vtp_file = fopen(tmp_filename, "w");
    if(broken_joint_vtp_file == NULL) { perror("Could not open broken joints VTP output file"); exit(1); }

    fprintf(broken_joint_vtp_file, "<?xml version='1.0'?>\n<VTKFile type='PolyData' version='0.1' byte_order='LittleEndian'>");
    fprintf(broken_joint_vtp_file, "<PolyData>\n<Piece NumberOfPoints='%d' NumberOfVerts='0' NumberOfLines='0' NumberOfStrips='0' NumberOfPolys='%d'>\n", 8*pv_num_broken_joints, 2*pv_num_broken_joints);

    /* Write out the points */
    fprintf(broken_joint_vtp_file, "<Points>\n<DataArray type='Float32' NumberOfComponents='3' format='ascii'>\n");

    for(i = 0; i < pv_num_broken_joints; ++i)
    {
        /* Write two polygons (8 points) : 0011 and 2233 */
        /* Each joint has 4 nodes */
        for(j = 0; j < 4; ++j)
        {
            for(k = 0; k < 2; ++k)
            {
                fprintf(broken_joint_vtp_file, "%.6f %.6f %.6f ", (double)pv_broken_joint_coordinate[4*i + j].x, (double)pv_broken_joint_coordinate[4*i + j].y, (double)pv_broken_joint_coordinate[4*i + j].z); 
            }
        }
    }

    fprintf(broken_joint_vtp_file, "\n</DataArray>\n</Points>\n");

    /* Write out the data associated with the points */
    fprintf(broken_joint_vtp_file, "<PointData>\n");
    
    /* Done writing point data */
    fprintf(broken_joint_vtp_file, "</PointData>\n");


    /* Write out the element data */
    fprintf(broken_joint_vtp_file, "<CellData Scalars='broken_mode_of_failure'>\n");
 //   fprintf(broken_joint_vtp_file, "<CellData Scalars='broken_mode_of_failure broken_sliding broken_opening broken_area'>\n");
 
        /* Write out property_id */
    fprintf(broken_joint_vtp_file, "<DataArray type='Float32' Name='broken_mode_of_failure' NumberOfComponents='1' format='ascii'>\n");

    for(i = 0; i < pv_num_broken_joints; ++i)
    {
        for(k = 0; k < 2; ++k)
        {
            fprintf(broken_joint_vtp_file, "%.6f ", (double)pv_broken_joint_mode_of_failure[i]); 
        }
    }
    fprintf(broken_joint_vtp_file, "\n</DataArray>\n");
    
    fprintf(broken_joint_vtp_file, "<DataArray type='Float32' Name='broken_sliding' NumberOfComponents='1' format='ascii'>\n");
    
    for(i = 0; i < pv_num_broken_joints; ++i)
    {
        for(k = 0; k < 2; ++k)
        {
            fprintf(broken_joint_vtp_file, "%.9f ", (double)pv_broken_joint_sliding[i]); 
        }
    }
    fprintf(broken_joint_vtp_file, "\n</DataArray>\n");
    
    fprintf(broken_joint_vtp_file, "<DataArray type='Float32' Name='broken_opening' NumberOfComponents='1' format='ascii'>\n");
    
    for(i = 0; i < pv_num_broken_joints; ++i)
    {
        for(k = 0; k < 2; ++k)
        {
            fprintf(broken_joint_vtp_file, "%.9f ", (double)pv_broken_joint_opening[i]); 
        }
    }
    fprintf(broken_joint_vtp_file, "\n</DataArray>\n");    
    
    fprintf(broken_joint_vtp_file, "<DataArray type='Float32' Name='broken_area' NumberOfComponents='1' format='ascii'>\n");
    
    for(i = 0; i < pv_num_broken_joints; ++i)
    {
        for(k = 0; k < 2; ++k)
        {
            fprintf(broken_joint_vtp_file, "%.9f ", (double)pv_broken_joint_area[i]); 
        }
    }
    fprintf(broken_joint_vtp_file, "\n</DataArray>\n"); 
	
    /* Done writing element data */
    fprintf(broken_joint_vtp_file, "\n</CellData>\n");

    /* Write out the polygons */
    fprintf(broken_joint_vtp_file, "<Polys>\n");

    /* Write connectivity table */
    fprintf(broken_joint_vtp_file, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");

    for(i = 0; i < 8*pv_num_broken_joints; ++i)
    {
        fprintf(broken_joint_vtp_file, "%d ", i); 
    }

    fprintf(broken_joint_vtp_file, "\n</DataArray>\n");

    /* Write offsets table */
    fprintf(broken_joint_vtp_file, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");

    for(i = 4; i <= 4*2*pv_num_broken_joints; i = i + 4)
    {
        fprintf(broken_joint_vtp_file, "%d ", i); 
    }

    fprintf(broken_joint_vtp_file, "\n</DataArray>\n");

    /* Done writing polygons */
    fprintf(broken_joint_vtp_file, "</Polys>\n");

    /* Done writing file */
    fprintf(broken_joint_vtp_file, "</Piece>\n</PolyData>\n</VTKFile>");

    if(fclose(broken_joint_vtp_file)) { perror("Could not close broken joint VTP output file"); exit(1); }

    /*
     *
     *
     *
     *
     */

    /* Open the yielded joint .vtp file */
    snprintf(tmp_filename, 1024, "%s_yielded_joints_%d.vtp", pv_file_prefix, timestep);

    /* Open the .vtp file */
    yielded_joint_vtp_file = fopen(tmp_filename, "w");
    if(yielded_joint_vtp_file == NULL) { perror("Could not open yielded joints VTP output file"); exit(1); }

    fprintf(yielded_joint_vtp_file, "<?xml version='1.0'?>\n<VTKFile type='PolyData' version='0.1' byte_order='LittleEndian'>");
    fprintf(yielded_joint_vtp_file, "<PolyData>\n<Piece NumberOfPoints='%d' NumberOfVerts='0' NumberOfLines='0' NumberOfStrips='0' NumberOfPolys='%d'>\n", 8*pv_num_yielded_joints, 2*pv_num_yielded_joints);

    /* Write out the points */
    fprintf(yielded_joint_vtp_file, "<Points>\n<DataArray type='Float32' NumberOfComponents='3' format='ascii'>\n");

    for(i = 0; i < pv_num_yielded_joints; ++i)
    {
        /* Write two polygons (8 points) : 0011 and 2233 */
        /* Each joint has 4 nodes */
        for(j = 0; j < 4; ++j)
        {
            for(k = 0; k < 2; ++k)
            {
                fprintf(yielded_joint_vtp_file, "%.6f %.6f %.6f ", (double)pv_yielded_joint_coordinate[4*i + j].x, (double)pv_yielded_joint_coordinate[4*i + j].y, (double)pv_yielded_joint_coordinate[4*i + j].z); 
            }
        }
    }

    fprintf(yielded_joint_vtp_file, "\n</DataArray>\n</Points>\n");

    /* Write out the data associated with the points */
    fprintf(yielded_joint_vtp_file, "<PointData>\n");
    
    /* Done writing point data */
    fprintf(yielded_joint_vtp_file, "</PointData>\n");


    /* Write out the element data */
    fprintf(yielded_joint_vtp_file, "<CellData Scalars='yielded_mode_of_failure'>\n");

        /* Write out property_id */
    fprintf(yielded_joint_vtp_file, "<DataArray type='Int32' Name='yielded_mode_of_failure' NumberOfComponents='1' format='ascii'>\n");

    for(i = 0; i < pv_num_yielded_joints; ++i)
    {
        for(k = 0; k < 2; ++k)
        {
            fprintf(yielded_joint_vtp_file, "%d ", pv_yielded_joint_mode_of_failure[i]); 
        }
    }

    fprintf(yielded_joint_vtp_file, "\n</DataArray>\n");

    /* Done writing element data */
    fprintf(yielded_joint_vtp_file, "\n</CellData>\n");

    /* Write out the polygons */
    fprintf(yielded_joint_vtp_file, "<Polys>\n");

    /* Write connectivity table */
    fprintf(yielded_joint_vtp_file, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");

    for(i = 0; i < 8*pv_num_yielded_joints; ++i)
    {
        fprintf(yielded_joint_vtp_file, "%d ", i); 
    }

    fprintf(yielded_joint_vtp_file, "\n</DataArray>\n");

    /* Write offsets table */
    fprintf(yielded_joint_vtp_file, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");

    for(i = 4; i <= 4*2*pv_num_yielded_joints; i = i + 4)
    {
        fprintf(yielded_joint_vtp_file, "%d ", i); 
    }

    fprintf(yielded_joint_vtp_file, "\n</DataArray>\n");

    /* Done writing polygons */
    fprintf(yielded_joint_vtp_file, "</Polys>\n");

    /* Done writing file */
    fprintf(yielded_joint_vtp_file, "</Piece>\n</PolyData>\n</VTKFile>");

    if(fclose(yielded_joint_vtp_file)) { perror("Could not close yielded joint VTP output file"); exit(1); }


    /*
     *
     *
     *
     *
     */
    
    /* Open the rebar .vtp file */
    snprintf(tmp_filename, 1024, "%s_rebars_%d.vtp", pv_file_prefix, timestep);

    /* Open the .vtp file */
    rebar_vtp_file = fopen(tmp_filename, "w");
    if(rebar_vtp_file == NULL) { perror("Could not open rebar VTP output file"); exit(1); }

    fprintf(rebar_vtp_file, "<?xml version='1.0'?>\n<VTKFile type='PolyData' version='0.1' byte_order='LittleEndian'>");
    fprintf(rebar_vtp_file, "<PolyData>\n<Piece NumberOfPoints='%d' NumberOfVerts='0' NumberOfLines='0' NumberOfStrips='0' NumberOfPolys='%d'>\n", 4*pv_num_rebars, pv_num_rebars);

    /* Write out the points */
    fprintf(rebar_vtp_file, "<Points>\n<DataArray type='Float32' NumberOfComponents='3' format='ascii'>\n");

    for(i = 0; i < pv_num_rebars; ++i)
    {
        /* Write one polygons (4 points) : 0011 */
        /* Each rebar has 2 nodes */
        for(j = 0; j < 2; ++j)
        {   
            for(k = 0; k < 2; ++k)
            {
                fprintf(rebar_vtp_file, "%.6f %.6f %.6f ", (double)pv_rebar_coordinate[2*i + j].x, (double)pv_rebar_coordinate[2*i + j].y, (double)pv_rebar_coordinate[2*i + j].z); 
            }
        }
    }

    fprintf(rebar_vtp_file, "\n</DataArray>\n</Points>\n");

    /* Write out the data associated with the points */
    fprintf(rebar_vtp_file, "<PointData>\n");
    
    /* Done writing point data */
    fprintf(rebar_vtp_file, "</PointData>\n");

    /* Write out the element data */

    fprintf(rebar_vtp_file, "<CellData Scalars='rebar_diameter'>\n");

    /* Write out rebar_area */
    fprintf(rebar_vtp_file, "<DataArray type='Float32' Name='rebar_area' NumberOfComponents='1' format='ascii'>\n");
    for(i = 0; i < pv_num_rebars; ++i)
    {
        fprintf(rebar_vtp_file, "%.6f ", (double)pv_rebar_area[i]);
    }
    fprintf(rebar_vtp_file, "\n</DataArray>\n");
    
    /* Write out rebar_stress */
    fprintf(rebar_vtp_file, "<DataArray type='Float32' Name='rebar_stress' NumberOfComponents='1' format='ascii'>\n");
    for(i = 0; i < pv_num_rebars; ++i)
    {
        fprintf(rebar_vtp_file, "%.6f ", (double)pv_rebar_stress[i]);

    }
    fprintf(rebar_vtp_file, "\n</DataArray>\n");
    
    /* Write out rebar_force */
    fprintf(rebar_vtp_file, "<DataArray type='Float32' Name='rebar_force' NumberOfComponents='1' format='ascii'>\n");
    for(i = 0; i < pv_num_rebars; ++i)
    {
        fprintf(rebar_vtp_file, "%.6f ", (double)pv_rebar_force[i]);
    }
    fprintf(rebar_vtp_file, "\n</DataArray>\n");
    
    /* Write out rebar_strain */
    fprintf(rebar_vtp_file, "<DataArray type='Float32' Name='rebar_strain' NumberOfComponents='1' format='ascii'>\n");
    for(i = 0; i < pv_num_rebars; ++i)
    {
        fprintf(rebar_vtp_file, "%.6f ", (double)pv_rebar_strain[i]);
    }
    fprintf(rebar_vtp_file, "\n</DataArray>\n");
    
    /* Write out rebar_activation_flag */
    fprintf(rebar_vtp_file, "<DataArray type='Int32' Name='rebar_activation_flag' NumberOfComponents='1' format='ascii'>\n");
    for(i = 0; i < pv_num_rebars; ++i)
    {
        fprintf(rebar_vtp_file, "%d ", pv_rebar_activation_flag[i]);
    }
    fprintf(rebar_vtp_file, "\n</DataArray>\n");
    
    
    /* Done writing element data */
    fprintf(rebar_vtp_file, "\n</CellData>\n");

    

    /* Write out the polygons */
    fprintf(rebar_vtp_file, "<Polys>\n");

    /* Write connectivity table */
    fprintf(rebar_vtp_file, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");

    for(i = 0; i < 4*pv_num_rebars; ++i)
    {
        fprintf(rebar_vtp_file, "%d ", i); 
    }

    fprintf(rebar_vtp_file, "\n</DataArray>\n");

    /* Write offsets table */
    fprintf(rebar_vtp_file, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");

    for(i = 4; i <= 4*1*pv_num_rebars; i = i + 4)
    {
        fprintf(rebar_vtp_file, "%d ", i); 
    }

    fprintf(rebar_vtp_file, "\n</DataArray>\n");

    /* Done writing polygons */
    fprintf(rebar_vtp_file, "</Polys>\n");

    /* Done writing file */
    fprintf(rebar_vtp_file, "</Piece>\n</PolyData>\n</VTKFile>");

    if(fclose(rebar_vtp_file)) { perror("Could not close rebar VTP output file"); exit(1); }
    
     /*
     *
     *
     *
     *
     */
    
    /* Open the 1D_joint .vtp file */
    snprintf(tmp_filename, 1024, "%s_1D_joints_%d.vtp", pv_file_prefix, timestep);

    /* Open the .vtp file */
    oneD_joint_vtp_file = fopen(tmp_filename, "w");
    if(oneD_joint_vtp_file == NULL) { perror("Could not open 1D joint VTP output file"); exit(1); }

    fprintf(oneD_joint_vtp_file, "<?xml version='1.0'?>\n<VTKFile type='PolyData' version='0.1' byte_order='LittleEndian'>");
    fprintf(oneD_joint_vtp_file, "<PolyData>\n<Piece NumberOfPoints='%d' NumberOfVerts='0' NumberOfLines='0' NumberOfStrips='0' NumberOfPolys='%d'>\n", 4*pv_num_unbroken_1D_joints, pv_num_unbroken_1D_joints);

    /* Write out the points */
    fprintf(oneD_joint_vtp_file, "<Points>\n<DataArray type='Float32' NumberOfComponents='3' format='ascii'>\n");

    for(i = 0; i < pv_num_unbroken_1D_joints; ++i)
    {
        /* Write one polygons (4 points) : 0011 */
        /* Each 1D joint has 2 nodes */
        for(j = 0; j < 2; ++j)
        {
            for(k = 0; k < 2; ++k)
            {
                fprintf(oneD_joint_vtp_file, "%.6f %.6f %.6f ", (double)pv_1D_joint_coordinate[2*i + j].x, (double)pv_1D_joint_coordinate[2*i + j].y, (double)pv_1D_joint_coordinate[2*i + j].z); 
            }
        }
    }

    fprintf(oneD_joint_vtp_file, "\n</DataArray>\n</Points>\n");

    /* Write out the data associated with the points */
    fprintf(oneD_joint_vtp_file, "<PointData>\n");
    //fprintf(oneD_joint_vtp_file, "<PointData Vectors='force'>\n");
    
   // /* Write out force */
   //fprintf(vtp_file, "<DataArray type='Float32' Name='force' NumberOfComponents='3' format='ascii'>\n");

   //for(i = 0; i < pv_num_unbroken_1D_joints; ++i)
   //{
   //    /* Write out the same force four times, one for each point */
   //    for(j = 0; j < 4; ++j)
   //    { 
   //       fprintf(oneD_joint_vtp_file, "%.6f %.6f %.6f ", (double)pv_1D_joint_force[i].x, (double)pv_1D_joint_force[i].y, (double)pv_1D_joint_force[i].z );
   //    }
   //}
   //fprintf(vtp_file, "\n</DataArray>\n");
    
    /* Done writing point data */
    fprintf(oneD_joint_vtp_file, "</PointData>\n");

    /* Write out the element data */

    fprintf(oneD_joint_vtp_file, "<CellData Scalars='1D_joint_stress'>\n");

    /* Write out 1D joint normal stress */
    fprintf(oneD_joint_vtp_file, "<DataArray type='Float32' Name='1D_joint_axial_stress' NumberOfComponents='1' format='ascii'>\n");
    for(i = 0; i < pv_num_unbroken_1D_joints; ++i)
    {
        fprintf(oneD_joint_vtp_file, "%.6f ", (double)pv_1D_joint_normal_stress[i]);
    }
    fprintf(oneD_joint_vtp_file, "\n</DataArray>\n");
    
    /* Write out 1D joint tangential stress */
    fprintf(oneD_joint_vtp_file, "<DataArray type='Float32' Name='1D_joint_shear_stress' NumberOfComponents='1' format='ascii'>\n");
    for(i = 0; i < pv_num_unbroken_1D_joints; ++i)
    {
        fprintf(oneD_joint_vtp_file, "%.6f ", (double)pv_1D_joint_tangential_stress[i]);
    }
    fprintf(oneD_joint_vtp_file, "\n</DataArray>\n");
    
    /* Write out 1D joint normal strain */
    fprintf(oneD_joint_vtp_file, "<DataArray type='Float32' Name='1D_joint_axial_displacement' NumberOfComponents='1' format='ascii'>\n");
    for(i = 0; i < pv_num_unbroken_1D_joints; ++i)
    {
        fprintf(oneD_joint_vtp_file, "%.6f ", (double)pv_1D_joint_normal_strain[i]);
    }
    fprintf(oneD_joint_vtp_file, "\n</DataArray>\n");
    
    /* Write out 1D joint tangential strain */
    fprintf(oneD_joint_vtp_file, "<DataArray type='Float32' Name='1D_joint_shear_displacement' NumberOfComponents='1' format='ascii'>\n");
    for(i = 0; i < pv_num_unbroken_1D_joints; ++i)
    {
        fprintf(oneD_joint_vtp_file, "%.6f ", (double)pv_1D_joint_tangential_strain[i]);
    }
    fprintf(oneD_joint_vtp_file, "\n</DataArray>\n");
    
    /* Write out 1D joint state */
    fprintf(oneD_joint_vtp_file, "<DataArray type='Int32' Name='1D_joint_state' NumberOfComponents='1' format='ascii'>\n");
    for(i = 0; i < pv_num_unbroken_1D_joints; ++i)
    {
        fprintf(oneD_joint_vtp_file, "%d ", pv_1D_joint_state[i]);
    }
    fprintf(oneD_joint_vtp_file, "\n</DataArray>\n");
      
    /* Done writing element data */
    fprintf(oneD_joint_vtp_file, "\n</CellData>\n");

    

    /* Write out the polygons */
    fprintf(oneD_joint_vtp_file, "<Polys>\n");

    /* Write connectivity table */
    fprintf(oneD_joint_vtp_file, "<DataArray type='Int32' Name='connectivity' format='ascii'>\n");

    for(i = 0; i < 4*pv_num_unbroken_1D_joints; ++i)
    {
        fprintf(oneD_joint_vtp_file, "%d ", i); 
    }

    fprintf(oneD_joint_vtp_file, "\n</DataArray>\n");

    /* Write offsets table */
    fprintf(oneD_joint_vtp_file, "<DataArray type='Int32' Name='offsets' format='ascii'>\n");

    for(i = 4; i <= 4*1*pv_num_unbroken_1D_joints; i = i + 4)
    {
        fprintf(oneD_joint_vtp_file, "%d ", i); 
    }

    fprintf(oneD_joint_vtp_file, "\n</DataArray>\n");

    /* Done writing polygons */
    fprintf(oneD_joint_vtp_file, "</Polys>\n");

    /* Done writing file */
    fprintf(oneD_joint_vtp_file, "</Piece>\n</PolyData>\n</VTKFile>");

    if(fclose(oneD_joint_vtp_file)) { perror("Could not close 1D joint VTP output file"); exit(1); }
    
    /*
     *
     *
     *
     *
     */

    /* Open the acoustic .vtp file */
    snprintf(tmp_filename, 1024, "%s_acoustic_%d.vtu", pv_file_prefix, timestep);

    /* Open the .vtp file */
    acoustic_vtp_file = fopen(tmp_filename, "w");
    if(acoustic_vtp_file == NULL) { perror("Could not open acoustic VTP output file"); exit(1); }

    fprintf(acoustic_vtp_file, "<?xml version='1.0'?>\n<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>");
    fprintf(acoustic_vtp_file, "<UnstructuredGrid>\n<Piece NumberOfPoints='%d' NumberOfCells='0'>\n", pv_num_acoustic_events);

    /* Write out the data associated with the points */
    fprintf(acoustic_vtp_file, "<PointData Scalars='acoustic_energy acoustic_magnitude acoustic_fracture_energy acoustic_event_mode'>\n");

    /* Write out acoustic energy */
    fprintf(acoustic_vtp_file, "<DataArray type='Float32' Name='acoustic_energy' NumberOfComponents='1' format='ascii'>\n");
  
    for(i = 0; i < pv_num_acoustic_events; ++i)
    {
        fprintf(acoustic_vtp_file, "%.6f ", (double)pv_acoustic_energy[i]);
    }
    fprintf(acoustic_vtp_file, "\n</DataArray>\n");

    /* Write out acoustic magnitude */
    fprintf(acoustic_vtp_file, "<DataArray type='Float32' Name='acoustic_magnitude' NumberOfComponents='1' format='ascii'>\n");
  
    for(i = 0; i < pv_num_acoustic_events; ++i)
    {
        fprintf(acoustic_vtp_file, "%.6f ", (double)pv_acoustic_magnitude[i]);
    }
    fprintf(acoustic_vtp_file, "\n</DataArray>\n");

    /* Write out acoustic fracture energy */
    fprintf(acoustic_vtp_file, "<DataArray type='Float32' Name='acoustic_fracture_energy' NumberOfComponents='1' format='ascii'>\n");
  
    for(i = 0; i < pv_num_acoustic_events; ++i)
    {
        fprintf(acoustic_vtp_file, "%.6f ", (double)pv_acoustic_fracture_energy[i]);
    }
    fprintf(acoustic_vtp_file, "\n</DataArray>\n");

    /* Write out acoustic event mode */
    fprintf(acoustic_vtp_file, "<DataArray type='Float32' Name='acoustic_event_mode' NumberOfComponents='1' format='ascii'>\n");
  
    for(i = 0; i < pv_num_acoustic_events; ++i)
    {
        fprintf(acoustic_vtp_file, "%.6f ", (double)pv_acoustic_event_mode[i]);
    }
    fprintf(acoustic_vtp_file, "\n</DataArray>\n");


    
    /* Done writing point data */
    fprintf(acoustic_vtp_file, "</PointData>\n");

    /* Write out the element data */
    fprintf(acoustic_vtp_file, "<CellData>\n");

    /* Done writing element data */
    fprintf(acoustic_vtp_file, "\n</CellData>\n");

    /* Write out the points */
    fprintf(acoustic_vtp_file, "<Points>\n<DataArray type='Float32' NumberOfComponents='3' format='ascii'>\n");

    for(i = 0; i < pv_num_acoustic_events; ++i)
    {
        fprintf(acoustic_vtp_file, "%.6f %.6f %.6f ", (double)pv_acoustic_coordinate[i].x, (double)pv_acoustic_coordinate[i].y, (double)pv_acoustic_coordinate[i].z); 
    }

    fprintf(acoustic_vtp_file, "\n</DataArray>\n</Points>\n");

    /* Write empty Cell */
    fprintf(acoustic_vtp_file, "<Cells><DataArray type='Int32' Name='connectivity'/><DataArray type='Int32' Name='offsets'/><DataArray type='UInt8' Name='types'/></Cells>\n");

    /* Done writing file */
    fprintf(acoustic_vtp_file, "</Piece>\n</UnstructuredGrid>\n</VTKFile>");

    if(fclose(acoustic_vtp_file)) { perror("Could not close acoustic VTP output file"); exit(1); }
     /*
     *
     *
     *
     *
     */

    /* Open the acoustic .vtp file */
    snprintf(tmp_filename, 1024, "%s_principal_direction_%d.vtu", pv_file_prefix, timestep);

    /* Open the .vtp file */
    principal_direction_vtp_file = fopen(tmp_filename, "w");
    if(principal_direction_vtp_file == NULL) { perror("Could not open acoustic VTP output file"); exit(1); }

    fprintf(principal_direction_vtp_file, "<?xml version='1.0'?>\n<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>");
    fprintf(principal_direction_vtp_file, "<UnstructuredGrid>\n<Piece NumberOfPoints='%d' NumberOfCells='0'>\n", pv_num_elements);
    
    /* Write out the data associated with the points */
    fprintf(principal_direction_vtp_file, "<PointData Vectors='principal_stress_direction_1 principal_stress_direction_2'>\n");
    
    /* Write out principal_stress_direction_1 */
    fprintf(principal_direction_vtp_file, "<DataArray type='Float32' Name='principal_stress_direction_1' NumberOfComponents='3' format='ascii'>\n");
    for(i = 0; i < pv_num_elements; ++i)
    {
        fprintf(principal_direction_vtp_file, "%.6f %.6f %.6f ", (double)pv_principal_stress_direction_1_andrea[i].x, (double)pv_principal_stress_direction_1_andrea[i].y, (double)pv_principal_stress_direction_1_andrea[i].z );
    }
    fprintf(principal_direction_vtp_file, "\n</DataArray>\n");
    
    /* Write out principal_stress_direction_2 */
    fprintf(principal_direction_vtp_file, "<DataArray type='Float32' Name='principal_stress_direction_2' NumberOfComponents='3' format='ascii'>\n");
    for(i = 0; i < pv_num_elements; ++i)
    {
        fprintf(principal_direction_vtp_file, "%.6f %.6f %.6f ", (double)pv_principal_stress_direction_2_andrea[i].x, (double)pv_principal_stress_direction_2_andrea[i].y, (double)pv_principal_stress_direction_2_andrea[i].z );
    }
    fprintf(principal_direction_vtp_file, "\n</DataArray>\n");
       
    /* Done writing point data */
    fprintf(principal_direction_vtp_file, "</PointData>\n");

    /* Write out the element data */
    fprintf(principal_direction_vtp_file, "<CellData>\n");

    /* Done writing element data */
    fprintf(principal_direction_vtp_file, "\n</CellData>\n");
    
    /* Write out the points */
    fprintf(principal_direction_vtp_file, "<Points>\n<DataArray type='Float32' NumberOfComponents='3' format='ascii'>\n");

    for(i = 0; i < pv_num_elements; ++i)
    {
        fprintf(principal_direction_vtp_file, "%.6f %.6f %.6f ", (double)pv_element_centroid_coordinate[i].x, (double)pv_element_centroid_coordinate[i].y, (double)pv_element_centroid_coordinate[i].z); 
    }

    fprintf(principal_direction_vtp_file, "\n</DataArray>\n</Points>\n");
    
    /* Write empty Cell */
    fprintf(principal_direction_vtp_file, "<Cells><DataArray type='Int32' Name='connectivity'/><DataArray type='Int32' Name='offsets'/><DataArray type='UInt8' Name='types'/></Cells>\n");

    /* Done writing file */
    fprintf(principal_direction_vtp_file, "</Piece>\n</UnstructuredGrid>\n</VTKFile>");

    if(fclose(principal_direction_vtp_file)) { perror("Could not close principal direction VTP output file"); exit(1); }
    
    
    /* Write the .pvd file */
}

/*********************PUBLIC********************************************************/
void Yod(
     namep,yd
     )
   CHR *namep; YD yd;
{ YDC ydc=&(yd->ydc);
  YDE yde=&(yd->yde);
  YDI ydi=&(yd->ydi);
  YDN ydn=&(yd->ydn);
  YDO ydo=&(yd->ydo);
  YDPE ydpe=&(yd->ydpe);
  YDPN ydpn=&(yd->ydpn);
  YDPJ ydpj=&(yd->ydpj);
  YDHF ydhf=&(yd->ydhf);
  YDFN ydfn=&(yd->ydfn);
  YDR ydr=&(yd->ydr);
  YDSB ydsb=&(yd->ydsb);
  YDPS ydps=&(yd->ydps);
  
  INT iprop, jprop, ihys, i, inode, nres;
  CHR namef[300];
  CHR cindex[50];
  static INT ncall=0;
  static INT ndefinitions=0;
  static INT ncallBRK=0;
  static INT ncallSFT=0;
  static INT ncallHYFR=0;
  static INT ncallFRAC=0;
  static INT ires=0;    //Counter for the restart data files
  FILE *fout=FILENULL;
  DBL tmp;
  static INT *nodes;

  INT ielem;
  
  //! PARAVIEW VARIABLES
  int inopo;
  int ibc;
  int ip;
  int num_nodes = 0;
  int num_elements;
  int num_joints;
  int num_rebars;
  int num_1D_joints;
  int num_dead_nodes;
  int broken_joint_id = 0;
  int yielded_joint_id = 0;
  int rebar_id = 0;
  int oneD_joint_id = 0;
  int AE_id = 0;
  double timestep_size;
  char *file_prefix;
  pv_vector_t nodal_coord;
  pv_vector_t nodal_coord_0;
  pv_vector_t nodal_coord_1;
  pv_vector_t nodal_coord_2;
  pv_vector_t nodal_coord_3;
  pv_vector_t nodal_vel;
  pv_vector_t nodal_force;
  pv_vector_t nodal_displ;
  pv_vector_t principal_stress_direction_1;
  pv_vector_t principal_stress_direction_2;
  pv_tensor_t element_stress;
  pv_tensor_t element_principal_stress;
  pv_tensor_t element_strain;
  int joint_state;
  int joint_yielding_mode;
  int joint_failure_mode;
  int j,k;
  float L1, L2;
  float q;
  int ipropc;
  int ipropd;
  int boundary_joints = 0;
  int broken_joints = 0;
  int yielded_joints = 0;
  int acoustic_events = 0;
  int unbroken_1D_joints = 0;
  float event_energy;
  float event_magnitude;
  float fracture_energy;
  int mode;
  pv_vector_t AE_coord;
  pv_vector_t element_centroid_coord;

  
  float nodal_fluid_pressure;
  
  double small = 1.0e-10;
  double e1x;
  double e1y;
  double h;
  double s1;
  double s2;
  double o1;
  double o2;
  float sliding;
  float opening;
  float area;
  float length;
  
  float rebar_area;
  float rebar_stress;
  float rebar_force;
  float rebar_strain;
  int rebar_activation_flag;
  
  int i0r,i1r;
  float oneD_joint_normal_stress;
  float oneD_joint_tangential_stress;
  float oneD_joint_normal_strain;
  float oneD_joint_tangential_strain;
  float oneD_joint_state;
  pv_vector_t oneD_joint_force;
  
 
  /* Elastic constants for plane strain transverse isotropy */
  DBL d1peex_pstrain; 
  DBL d1peey_pstrain;
  DBL d1pemx_pstrain;
  DBL d1pemy_pstrain;
  
  /*static FILE *out1=FILENULL;
  if(out1 == FILENULL)
  {	out1=fopen("Yod.txt", "a");
  }*/
  
  DBL voli,volc;
  DBL  B[2][2]; /* left Cauchy-Green strain tensor */
  DBL  D[2][2]; /* rate of deformation (stretching) tensor */
  DBL  E[2][2]; /* strain tensor (small strains) */
  DBL  F[2][2]; /* deformation gradient in global base */
  DBL F0[2][2]; /* initial local base */
  DBL FX[2][2]; /* current local base */
  DBL F0inv[2][2]; /* global base in initial local base */
  DBL FXinv[2][2]; /* global base in current local base */
  DBL  L[2][2]; /* velocity gradient in global base */
  DBL LX[2][2]; /* vel. gradient in current local base = delta x/delta X */
  DBL  T[2][2]; /* Cauchy stress */
  DBL dpemu;
  DBL dpela;
  DBL dpeks;
  
  INT i1usan;
  DBL d1peex;
  DBL d1peey;
  DBL d1pemx;
  DBL d1pemy;
  DBL d1peg;
  DBL dpeem;
  DBL dpenu;
  
  INT output_frequency;
  
  static DBL total_fracture_area = 0.0;
  static DBL total_fracture_length = 0.0;
  
  /* Writing restart data */ 
  /* Initialize the files only during the 1st time step */
  if(ydc->icresf>0)
  { nres=ydc->mcstep/ydc->icresf; //Number of restart data files
    if(nres>0)
    { if(ydc->ncstep==0)
      { i=nres*sizeof(FILE*);
        if(i>0)ydo->f2orsf=(FILE**)MALLOC(i);
        for(i=0;i<nres;i++)
        { ydo->f2orsf[i]=FILENULL;
        }
        nodes=INT1NULL;
        nodes=TalINT1(ydn->nnopst);
        for(ielem=0;ielem<yde->nelemst;ielem++)
        { //printf("st %d: %d \t %d \t %d \n", ielem, yde->i2eltost[0][ielem], yde->i2eltost[1][ielem], yde->i2eltost[2][ielem]);
          for(i=0;i<3;i++)
          { nodes[yde->i2eltost[i][ielem]]=yde->i2elto[i][ielem];            
          }          
        } 
        /*for(inode=0;inode<ydn->nnopst;inode++)
        { printf("inode = %d nodes[inode] = % d \n", inode, nodes[inode]);
        }*/
      }
      //printf("%d \n", ydc->ncstep);
      if(((ydc->ncstep%ydc->icresf)==0) || (ydc->ncstep==(ydc->mcstep-1)))
      { CHRcpy(namef,namep);
        SINTw(cindex,ires,0);
        CHRcat(namef,"res");
        CHRcat(namef,cindex);
        FILE *fp;
        fp = fopen (namef, "w+");

        fprintf(fp, "/* Current restart data stored for \nTime Step: %ld @ Time: %lf */ \n", ydc->ncstep, ydc->dctime);
        fprintf(fp, "\n/YD/YDN/D2NCC 21 %ld %ld \n", ydn->nnodim, ydn->nnopst);
        for(inode=0;inode<ydn->nnopst;inode++)
        { fprintf( fp , "%lf \t %lf \n", ydn->d2ncc[0][nodes[inode]], ydn->d2ncc[1][nodes[inode]] );
        }
        fprintf(fp, "\n/YD/YDN/D2NVC 21 %ld %ld \n", ydn->nnodim, ydn->nnopst);
        for(inode=0;inode<ydn->nnopst;inode++)
        { fprintf( fp , "%lf \t %lf \n", ydn->d2nvc[0][nodes[inode]], ydn->d2nvc[1][nodes[inode]] );
        }
        fclose(fp);
        ires++;
  } } }

  /* Output inputfilename.yhBRK */
  ncallBRK=ncallBRK+1;
  if((ydc->dctime)==(yde->d1etbrk[(yde->nebrk)-1]))
  { if((ydo->foebrk)==FILENULL)
      { CHRcpy(namef,namep);
        CHRcat(namef,"hBRK");
        ydo->foebrk=fopen(namef,"a");
      }
      if((ydo->foebrk)!=FILENULL)
      {
	for(i=0;i<yde->netbrk;i++)
	{ INTw(ydo->foebrk,yde->i1ebrk[(yde->nebrk)-1-i],5); /* ID of the broken joint element */
          CHRwtab(ydo->foebrk);
          //if(yde->i1ebrkf[yde->i1ebrk[(yde->nebrk)-1-i]]==1)
	  //{ CHRw(ydo->foebrk,"3"); }                         /* Yielding in tension and shear */
	  //if(yde->i1ebrkf[yde->i1ebrk[(yde->nebrk)-1-i]]==2)
	  //{ CHRw(ydo->foebrk,"1"); }                         /* Yielding in tension */
	  //if(yde->i1ebrkf[yde->i1ebrk[(yde->nebrk)-1-i]]==3) 
	  //{ CHRw(ydo->foebrk,"2"); }                         /* Yielding in shear */
	  DBLw(ydo->foebrk,yde->d1ebrkf[yde->i1ebrk[(yde->nebrk)-1-i]],10); //! Mode of failure 
          CHRwtab(ydo->foebrk);
          DBLw(ydo->foebrk,ydc->dctime,10); /* time of breakage */
          CHRwtab(ydo->foebrk);
          DBLw(ydo->foebrk,yde->d2ecbrk[0][(yde->nebrk)-1-i],10); /* x coordinate of the breakage (i.e. joint centroid) */
          CHRwtab(ydo->foebrk);
          DBLw(ydo->foebrk,yde->d2ecbrk[1][(yde->nebrk)-1-i],10); /* y coordinate of the breakage (i.e. joint centroid) */
          CHRwtab(ydo->foebrk);
          DBLw(ydo->foebrk,ydn->d2ncc[0][yde->i2elto[0][yde->i1ebrk[(yde->nebrk)-1-i]]],10); /* x coordinate of first joint node */
          CHRwtab(ydo->foebrk);
          DBLw(ydo->foebrk,ydn->d2ncc[1][yde->i2elto[0][yde->i1ebrk[(yde->nebrk)-1-i]]],10); /* y coordinate of first joint node */
          CHRwtab(ydo->foebrk);
          DBLw(ydo->foebrk,ydn->d2ncc[0][yde->i2elto[1][yde->i1ebrk[(yde->nebrk)-1-i]]],10); /* x coordinate of second joint node */
          CHRwtab(ydo->foebrk);
          DBLw(ydo->foebrk,ydn->d2ncc[1][yde->i2elto[1][yde->i1ebrk[(yde->nebrk)-1-i]]],10); /* y coordinate of second joint node */
          CHRwtab(ydo->foebrk);  
          DBLw(ydo->foebrk,yde->d1efe[yde->i1ebrk[(yde->nebrk)-1-i]],10); /* fracture energy released by the breakage */
          CHRwtab(ydo->foebrk);
          DBLw(ydo->foebrk,-(yde->d1efe[yde->i1ebrk[(yde->nebrk)-1-i]])/(yde->d1elbrk[(yde->nebrk)-1-i]),10); /* element fracture energy release rate */
          CHRwtab(ydo->foebrk);
	  DBLw(ydo->foebrk,yde->d1edke[yde->i1ebrk[(yde->nebrk)-1-i]],10); /* kinetic energy released by the breakage */
          CHRwtab(ydo->foebrk);
	  DBLw(ydo->foebrk,yde->d1etmke[yde->i1ebrk[(yde->nebrk)-1-i]],10); /* time associated with maximum kinetic energy */
          CHRwtab(ydo->foebrk);
          if(ydfn->iusefn > 0)
	  { INTw(ydo->foebrk,yde->i1edfnf[yde->i1ebrk[(yde->nebrk)-1-i]],5); } /* Cohesive DFN flag of the broken joint element, 0 = regular joint, 1 = DFN joint */
          CHRwtab(ydo->foebrk);
          DBLw(ydo->foebrk,yde->d1etike[yde->i1ebrk[(yde->nebrk)-1-i]],10); /* initial time of KE monitoring window */
          CHRwcr(ydo->foebrk);
	    
  } } }
  if((ncallBRK>100)||(ydc->ncstep)>=(ydc->mcstep-2))
  { ncallBRK=0;
    if((ydo->foebrk)!=FILENULL)fclose(ydo->foebrk);
      ydo->foebrk=FILENULL;
  }
  /* Output inputfilename.yhSFT */
  ncallSFT=ncallSFT+1;
  if((ydc->dctime)==(yde->d1etsft[(yde->nesft)-1]))
  { if((ydo->foesft)==FILENULL)
      { CHRcpy(namef,namep);
        CHRcat(namef,"hSFT");
        ydo->foesft=fopen(namef,"a");
      }
      if((ydo->foesft)!=FILENULL)
      {
	    for(i=0;i<yde->netsft;i++)
		{
		  INTw(ydo->foesft,yde->i1esft[(yde->nesft)-1-i],5); /* number of the joint element softened */
          CHRwtab(ydo->foesft);
          if(yde->i1esftf[yde->i1esft[(yde->nesft)-1-i]]==1)
		  { CHRw(ydo->foesft,"3"); }                       /* Yielding in tension and shear */
		  if(yde->i1esftf[yde->i1esft[(yde->nesft)-1-i]]==2)
		  { CHRw(ydo->foesft,"1"); }                         /* Yielding in tension */
		  if(yde->i1esftf[yde->i1esft[(yde->nesft)-1-i]]==3) 
		  { CHRw(ydo->foesft,"2"); }                         /* Yielding in shear */
          CHRwtab(ydo->foesft);
          DBLw(ydo->foesft,ydc->dctime,10); /* time of yielding */
          CHRwtab(ydo->foesft);
          DBLw(ydo->foesft,yde->d2ecsft[0][(yde->nesft)-1-i],10); /* x coordinate of yielding (i.e. joint centroid) */
          CHRwtab(ydo->foesft);
          DBLw(ydo->foesft,yde->d2ecsft[1][(yde->nesft)-1-i],10); /* y coordinate of yielding (i.e. joint centroid) */
          CHRwtab(ydo->foesft);
          DBLw(ydo->foesft,ydn->d2ncc[0][yde->i2elto[0][yde->i1esft[(yde->nesft)-1-i]]],10); /* x coordinate of first joint node */
          CHRwtab(ydo->foesft);
          DBLw(ydo->foesft,ydn->d2ncc[1][yde->i2elto[0][yde->i1esft[(yde->nesft)-1-i]]],10); /* y coordinate of first joint node */
          CHRwtab(ydo->foesft);
          DBLw(ydo->foesft,ydn->d2ncc[0][yde->i2elto[1][yde->i1esft[(yde->nesft)-1-i]]],10); /* x coordinate of second joint node */
          CHRwtab(ydo->foesft);
          DBLw(ydo->foesft,ydn->d2ncc[1][yde->i2elto[1][yde->i1esft[(yde->nesft)-1-i]]],10); /* y coordinate of second joint node */
          CHRwcr(ydo->foesft);
  } } }
  if((ncallSFT>100)||(ydc->ncstep)>=(ydc->mcstep-2))
  { ncallSFT=0;
    if((ydo->foesft)!=FILENULL)fclose(ydo->foesft);
      ydo->foesft=FILENULL;
  }
  
  /* Output inputfilename.yhHF */
  ncallHYFR=ncallHYFR+1;
  if((ydhf->iusehf==1)&&(ydhf->ihftyp==2))
  { if((ydo->fohyfr)==FILENULL)
    { CHRcpy(namef,namep);
      CHRcat(namef,"hHF");
      ydo->fohyfr=fopen(namef,"a");
    }
    if((ydo->fohyfr)!=FILENULL)
    { DBLw(ydo->fohyfr,ydc->dctime,10); /* time */
      CHRwtab(ydo->fohyfr);
      fprintf(ydo->fohyfr,"%f",ydhf->flupres);
      //DBLw(ydo->fohyfr,ydhf->flupres,10); /* fluid pressure */
      CHRwtab(ydo->fohyfr);
      fprintf(ydo->fohyfr,"%f",ydhf->flumass);
      //DBLw(ydo->fohyfr,ydhf->flumass,10); /* fluid mass */
      CHRwtab(ydo->fohyfr);
      fprintf(ydo->fohyfr,"%f",ydhf->fluvol);
      //DBLw(ydo->fohyfr,ydhf->flumass,10); /* fluid volume */
      CHRwcr(ydo->fohyfr);
  } } 
  if((ncallHYFR>100)||(ydc->ncstep)>=(ydc->mcstep-2))
  { ncallHYFR=0;
    if((ydo->fohyfr)!=FILENULL)fclose(ydo->fohyfr);
      ydo->fohyfr=FILENULL;
  }
  
  
  
  ncall=ncall+1;
  /* Output history */
  for(ihys=0;ihys<(ydo->nohys);ihys++)
  { if((ydo->d1ohyt[ihys])==(ydc->dctime))
    { if((ydo->f2ohyf[ihys])==FILENULL)
      { CHRcpy(namef,namep);
        SINTw(cindex,ihys,0);
        CHRcat(namef,"h");
        CHRcat(namef,cindex);
        ydo->f2ohyf[ihys]=fopen(namef,"a");
      }
      if((ydo->f2ohyf[ihys])!=FILENULL)
      { tmp=(ydo->d1ohyt[ihys])*(ydo->d1ohyc[ihys]);
        DBLw((ydo->f2ohyf[ihys]),tmp,10);
        CHRwtab(ydo->f2ohyf[ihys]);
        tmp=(ydo->d1ohys[ihys])*(ydo->d1ohyf[ihys]);
        DBLw((ydo->f2ohyf[ihys]),tmp,10);
        CHRwcr(ydo->f2ohyf[ihys]);
  } } }
  if((ncall>100)||((ydc->ncstep)>=(ydc->mcstep-2)))
  { ncall=0;
    for(ihys=0;ihys<(ydo->nohys);ihys++)
    { if((ydo->f2ohyf[ihys])!=FILENULL)fclose(ydo->f2ohyf[ihys]);
      ydo->f2ohyf[ihys]=FILENULL;
  } }
  
  
  /* output animation */
  fout=FILENULL;
  
  
  /* Determine output frequency */
  if((yde->nebrk > ydc->icoutnf) && (ydc->icoutrf > 0)) /* Total number of fractures is greater than threshold value */
  { output_frequency = ydc->icoutrf; }
  else
  { output_frequency = ydc->icoutf; }
  
  //if((ydc->ncstep%ydc->icoutf)==0)
  if((ydc->ncstep%output_frequency)==0)
  { /*CHRcpynoext(namef,namep);
    SINTw(cindex,ydc->icouti,0);
    CHRcat(namef,cindex);
    CHRcat(namef,".ym");
    fout=fopen(namef,"w");
    if(fout!=FILENULL)
    { CHRw(fout,"CODED"); CHRwcr(fout);
      for(iprop=0;iprop<ydpe->nprop;iprop++)
      { if((ydpe->i1ptyp[iprop])==(YTE2TRIELS))
        { Yod2TRIELS(
          yde->nelem ,
          fout       ,
          ydc->dcsizc,ydc->dcsizs,ydc->dcsizv,
          ydpe->d1peks[iprop],ydpe->d1pela[iprop],ydpe->d1pemu[iprop],
          ydpe->d1pero[iprop],
          ydc->icoutp       ,iprop             ,
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],
          yde->i1elpr,yde->i2elto
          );
        }
        else if((ydpe->i1ptyp[iprop])==(YTE2TRISOF))
        { Yod2TRISOF(  /* small strain elastic triangle output */
     /*     yde->nelem ,
          fout       ,
          ydc->dcsizc,ydc->dcsizs,ydc->dcsizv,
          ydpe->d1peks[iprop],ydpe->d1pela[iprop],ydpe->d1pemu[iprop],
          ydpe->d1pero[iprop],
          ydc->icoutp       ,iprop             ,
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nvc[0],
          ydn->d2nvc[1],yde->d2elst[ydpe->i1psde[iprop]]    ,yde->i1elpr,
          yde->i2elto
          );
        }
        else if((ydpe->i1ptyp[iprop])==(YTE2TRIRIG))
        { Yod2TRIRIG(  /* rigid triangle output */
      /*    yde->nelem ,
          fout       ,
          ydc->dcsizc,ydc->dcsizv,
          ydc->icoutp,iprop      ,
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nvc[0],ydn->d2nvc[1],yde->i1elpr,
          yde->i2elto
          );
      } }
      for(jprop=0;jprop<ydpj->npjset;jprop++)
      { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
        { Yod2JOINTSINTACT(  /* 2D joint output */
      /*    yde->nelem ,
          fout       ,
          ydc->dcsizc, ydc->dcsizv,
          ydpj->d1pjft[jprop], ydpj->d1pjgf[jprop],
          ydpj->d1pjpe[jprop],
          ydc->icoutp, (jprop+ydpe->nprop),
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nvc[0],ydn->d2nvc[1],
          yde->d2elst[ydpj->i1psde[jprop]],
          yde->i1elpr,yde->i2elto,yde->d1elfs
          );
          Yod2JOINTSBROKEN(
          yde->nelem ,
          fout       ,
          ydc->dcsizc,ydc->dcsizv,
          ydc->icoutp,(jprop+ydpe->nprop),
          ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nvc[0],ydn->d2nvc[1],
          yde->d2elst[ydpj->i1psde[jprop]],
          yde->i1elpr,yde->i2elto
          );
      } }
      ydc->icouti=ydc->icouti+1;
      fclose(fout);
    } */
  
    //! PARAVIEW OUTPUT
    num_elements = yde->nelemst;
    num_joints = yde->nelem - num_elements;
    file_prefix = namep;
    timestep_size = ydc->dcstec;
    num_rebars = ydr->nrepo/2;
    //num_1D_joints = ydr->nbrjointrb/2;
    num_1D_joints = ydr->nbrjointrb;
    
    //printf("num_1D_joints: %d\n",num_1D_joints);
    
    if(ydc->ncstep==0) //! Initialize files and translate Munjiza IDs only at the first timestep
    { 
    for(ibc=0;ibc < ydpn->npnset;ibc++) //! Loop over nodes that are actually used 
    { for(inopo=0;inopo<ydn->nnopo;inopo++)
      { if((ydn->i1nopr[inopo]==ibc)&&(ydn->d1nmct[inopo]>EPSILON))
        { num_nodes++;
        }
      }
    }
    num_dead_nodes = ydn->nnopo-num_nodes;
    pv_init(num_elements, num_joints, num_rebars, num_1D_joints, num_dead_nodes, num_joints, file_prefix, timestep_size); //! File initialization
      for(iprop=0;iprop<ydpe->nprop;iprop++) //! Munjiza translator for triangles
      { if( (ydpe->i1ptyp[iprop])==(YTE2TRIELS) ||
            (ydpe->i1ptyp[iprop])==(YTE2PLANESTRESS) ||
            (ydpe->i1ptyp[iprop])==(YTE2PLANESTRAIN) ) 
        { for(ielem=0;ielem<yde->nelem;ielem++)
          { if(yde->i1elpr[ielem]==iprop)
            { pv_register_munjiza_element(ielem); } 
          } 
        }
      }
      for(jprop=0;jprop<ydpj->npjset;jprop++) //! Munjiza translator for joints
      { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
        { for(ielem=0;ielem<yde->nelem;ielem++)
          { if(yde->i1elpr[ielem]==(jprop+(ydpe->nprop)))
            { pv_register_munjiza_joint(ielem); }
          }
        }
      }
    }
    //! Triangle nodal coordinates, velocities,forces and displacements
    for(ibc=0;ibc < ydpn->npnset;ibc++) //! Loop over nodes that are actually used 
    { for(inopo=0;inopo<ydn->nnopo;inopo++)
      { if((ydn->i1nopr[inopo]==ibc)&&(ydn->d1nmct[inopo]>EPSILON))
        { nodal_coord.x = ydn->d2ncc[0][inopo];
          nodal_coord.y = ydn->d2ncc[1][inopo];
          nodal_coord.z = 0.0;
          pv_write_nodal_current_coord(inopo, nodal_coord);
       
          nodal_vel.x = ydn->d2nvc[0][inopo];
          nodal_vel.y = ydn->d2nvc[1][inopo];
          nodal_vel.z = 0.0;
          pv_write_nodal_velocity(inopo, nodal_vel);
       
          nodal_force.x = ydn->d2nfc[0][inopo];
          nodal_force.y = ydn->d2nfc[1][inopo];
          nodal_force.z = 0.0;
          pv_write_nodal_force(inopo, nodal_force);
       
          //nodal_displ.x = ydn->d2ncc[0][inopo]-ydn->d2nci[0][inopo];
          //nodal_displ.y = ydn->d2ncc[1][inopo]-ydn->d2nci[1][inopo];
          nodal_displ.x = ydn->d2ncc[0][inopo]-ydn->d2nc0[0][inopo];
          nodal_displ.y = ydn->d2ncc[1][inopo]-ydn->d2nc0[1][inopo];
          nodal_displ.z = 0.0;
          pv_write_nodal_displacement(inopo, nodal_displ);
          
          //fprintf(out1,"inopo: %d\t ydn->i1nowe[inopo]: %d\n",inopo,ydn->i1nowe[inopo]);
          //! Output nodal fluid pressure for hydrofrac
          if(ydhf->iusehf==1)
          { if(ydn->i1nowe[inopo] == 1 || ydn->i1nowe[inopo]==2)
            { if((ydhf->ihftyp==1) || (ydhf->ihftyp==2))
              { nodal_fluid_pressure = ydhf->flupres; }
              else if(ydhf->ihftyp==3)
              { nodal_fluid_pressure = ydn->d1nfp[inopo]; }
            }
            else
            { nodal_fluid_pressure = 0.0; }
          }
          else
          { nodal_fluid_pressure = 0.0; }
          pv_write_nodal_fluid_pressure(inopo, nodal_fluid_pressure); 
        }
      }
    }
    //! Triangle property ID, stress, strain, connectivity table 
    for(iprop=0;iprop<ydpe->nprop;iprop++) //! Loop over triangles
    { if((ydpe->i1ptyp[iprop])==(YTE2TRIELS) ||
        (ydpe->i1ptyp[iprop])==(YTE2PLANESTRESS) ||
        (ydpe->i1ptyp[iprop])==(YTE2PLANESTRAIN) )
      { for(ielem=0;ielem<yde->nelem;ielem++)
        { if(yde->i1elpr[ielem]==iprop)
          { pv_write_element_property_id(ielem, iprop);
            
            dpela=ydpe->d1pela[iprop];
            dpemu=ydpe->d1pemu[iprop];
            dpeks=ydpe->d1peks[iprop];
	    dpeem=ydpe->d1peem[iprop];
	    dpenu=ydpe->d1penu[iprop];
	    
	    /* Transversely isotropic elastic parameters */
	    i1usan=ydpe->i1usan[iprop];
	    d1peex=ydpe->d1peex[iprop];
	    d1peey=ydpe->d1peey[iprop];
	    d1pemx=ydpe->d1pemx[iprop];
	    d1pemy=ydpe->d1pemy[iprop];
	    d1peg=ydpe->d1peg[iprop];
	    
	    
	    
 
            //! Evaluate element strain and stress tensors here
            for(i=1;i<3;i++)
            { F0[0][i-1]=ydn->d2nci[0][(yde->i2elto[i][ielem])]-ydn->d2nci[0][(yde->i2elto[0][ielem])];
              F0[1][i-1]=ydn->d2nci[1][(yde->i2elto[i][ielem])]-ydn->d2nci[1][(yde->i2elto[0][ielem])];
              FX[0][i-1]=ydn->d2ncc[0][(yde->i2elto[i][ielem])]-ydn->d2ncc[0][(yde->i2elto[0][ielem])];
              FX[1][i-1]=ydn->d2ncc[1][(yde->i2elto[i][ielem])]-ydn->d2ncc[1][(yde->i2elto[0][ielem])];
              LX[0][i-1]=ydn->d2nvc[0][(yde->i2elto[i][ielem])]-ydn->d2nvc[0][(yde->i2elto[0][ielem])];
              LX[1][i-1]=ydn->d2nvc[1][(yde->i2elto[i][ielem])]-ydn->d2nvc[1][(yde->i2elto[0][ielem])];
            }
            YMATINV2(F0,F0inv,voli);
            YMATINV2(FX,FXinv,volc);
            for(i=0;i<2;i++)
            { for(j=0;j<2;j++)
              { F[i][j]=R0;
                L[i][j]=R0;
                for(k=0;k<2;k++)
                { F[i][j]=F[i][j]+FX[i][k]*F0inv[k][j];
                  L[i][j]=L[i][j]+LX[i][k]*FXinv[k][j];
                } 
	      } 
	    }
            for(i=0;i<2;i++)
            { for(j=0;j<2;j++)
              { B[i][j]=R0;
                for(k=0;k<2;k++)
                { B[i][j]=B[i][j]+F[i][k]*F[j][k]; // left Cauchy-Green strain 
                }
                D[i][j]=RP5*(L[i][j]+L[j][i]);     // rate of deformation      
                E[i][j]=RP5*B[i][j];               // small strain             
                if(i==j) E[i][j]=E[i][j]-RP5;
	      } 
	    }
	    if(i1usan==1) /* Apply transversely isotropic elastic constitutive law */
            { T[0][0] = (d1peex / (1 - d1pemx * d1pemy)) * (E[0][0] + d1pemy * E[1][1]) + dpeks * D[0][0];
              T[1][1] = (d1peey / (1 - d1pemx * d1pemy)) * (E[1][1] + d1pemx * E[0][0]) + dpeks * D[1][1];
              T[0][1] = (2 * d1peg) * E[0][1] + dpeks * D[0][1];
              T[1][0] = (2 * d1peg) * E[1][0] + dpeks * D[1][0];
            }
            else if(i1usan==2) /* Apply transversely isotropic elastic constitutive law (plane strain) */
            { d1peex_pstrain = d1peex / (1 - d1pemx * d1pemy);
              d1peey_pstrain = d1peey / (1 - d1pemy * d1pemx);
	      d1pemx_pstrain = (d1pemx + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	      d1pemy_pstrain = (d1pemy + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	      T[0][0] = (d1peex_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[0][0] + d1pemy_pstrain * E[1][1]) + dpeks * D[0][0];
              T[1][1] = (d1peey_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[1][1] + d1pemx_pstrain * E[0][0]) + dpeks * D[1][1];
              T[0][1] =  (2 * d1peg) * E[0][1] + dpeks * D[0][1];
              T[1][0] =  (2 * d1peg) * E[1][0] + dpeks * D[1][0];
            }
            else
	    { if (ydpe->i1ptyp[iprop] == YTE2TRIELS)
              { for(i=0;i<2;i++)    
                { for(j=0;j<2;j++)
                { T[i][j]=R2*dpemu*E[i][j]*(voli/volc)+dpeks*D[i][j];
                }
                T[i][i]=T[i][i]+dpela*(volc/voli-voli/volc);
                }
              }
              /* Plane Stress formulation based on E,nu */
              else if (ydpe->i1ptyp[iprop] == YTE2PLANESTRESS)
              {        
                 T[0][0] = (dpeem / (1 - dpenu*dpenu)) * (E[0][0] + dpenu * E[1][1]) + dpeks * D[0][0];
                 T[1][1] = (dpeem / (1 - dpenu*dpenu)) * (E[1][1] + dpenu * E[0][0]) + dpeks * D[1][1];
                 T[0][1] = (dpeem / (1 + dpenu)) * E[0][1] + dpeks * D[0][1];
                 T[1][0] = (dpeem / (1 + dpenu)) * E[1][0] + dpeks * D[1][0];
              }
              /* Plane Strain formulation based on E,nu */
              else if (ydpe->i1ptyp[iprop] == YTE2PLANESTRAIN)
              {
                 T[0][0] = dpeem/((1+dpenu)*(1-2*dpenu)) * ( (1-dpenu)*E[0][0] +  dpenu*E[1][1]) + dpeks * D[0][0];
                 T[1][1] = dpeem/((1+dpenu)*(1-2*dpenu)) * ( (1-dpenu)*E[1][1] +  dpenu*E[0][0]) + dpeks * D[1][1];
                 T[0][1] = (dpeem / (1+dpenu)) * E[0][1] + dpeks * D[0][1];
                 T[1][0] = (dpeem / (1+dpenu)) * E[1][0] + dpeks * D[1][0];
              }
	    }
            
            element_stress.a11=T[0][0];
            element_stress.a12=T[0][1];
            element_stress.a13=0.0;
            element_stress.a21=T[1][0];
            element_stress.a22=T[1][1];
            element_stress.a23=0.0;
            element_stress.a31=0.0;
            element_stress.a32=0.0;
            element_stress.a33=0.0;
	    pv_write_stress(ielem, element_stress);
            	    
            //!Element principal stress tensor here
	    q=0.5*((T[0][0]+T[1][1])+((T[0][0]+T[1][1])/(ABS(T[0][0]+T[1][1])+EPSILON))*
	      SQRT(((T[0][0]-T[1][1])*(T[0][0]-T[1][1]))+4.0*T[0][1]*T[1][0]));

	    //! First eigenvalue
	    //if((T[0][0]+T[1][1]) < EPSILON) //!GPUFEMDEM version
	    if(ABS(q)<EPSILON)
	    { //L1 = -SQRT((T[0][1]*T[1][0]-T[0][0]*T[1][1])); //!GPUFEMDEM version
	      L1 = 0.0;
	    }
	    else 
	    { L1 = q; }
	    //! Second eigenvalue
            //if((T[0][0]+T[1][1]) < EPSILON) //!GPUFEMDEM version
	    if(ABS(q)<EPSILON)
	    { //L2 = SQRT((T[0][1]*T[1][0]-T[0][0]*T[1][1])); //!GPUFEMDEM version
	      L2 = 0.0;
            }
	    else 
	    { L2 = (-T[0][1]*T[1][0]+T[0][0]*T[1][1])/(q); }

	    //! Diagonal matrix with eigenvalues in increasing order
	    if(L1 <= L2) 
	    { element_principal_stress.a11 = L1; 
	      element_principal_stress.a22 = L2; 
	    }
	    else   
	    { element_principal_stress.a11 = L2;
	      element_principal_stress.a22 = L1;
	    }
	    //element_principal_stress.a11=0.0;
	    element_principal_stress.a12=0.0;
	    element_principal_stress.a13=0.0;
	    element_principal_stress.a21=0.0;
	    //element_principal_stress.a22=0.0;
	    element_principal_stress.a23=0.0;
	    element_principal_stress.a31=0.0;
	    element_principal_stress.a32=0.0;
	    element_principal_stress.a33=0.0;
	    pv_write_principal_stress(ielem, element_principal_stress);
            
            //!Element strain
	    element_strain.a11=E[0][0];
            element_strain.a12=E[0][1];
            element_strain.a13=0.0;
            element_strain.a21=E[1][0];
            element_strain.a22=E[1][1];
            element_strain.a23=0.0;
            element_strain.a31=0.0;
            element_strain.a32=0.0;
            element_strain.a33=0.0;
            pv_write_strain(ielem, element_strain);
	    
	    //!Element principal stress directions here	TODO Add additional condition (see GPUFEMDEM)    
	    
	    //principal_stress_direction_1.x = T[0][1]/SQRT((T[0][1]*T[0][1]+(L1-T[0][0])*(L1-T[0][0])));
	    //principal_stress_direction_1.y = (L1-T[0][0])/SQRT((T[0][1]*T[0][1]+(L1-T[0][0])*(L1-T[0][0])));
	    
	    //principal_stress_direction_2.x = T[0][1]/SQRT((T[0][1]*T[0][1]+(L2-T[0][0])*(L2-T[0][0])));
	    //principal_stress_direction_2.y = (L1-T[0][0])/SQRT((T[0][1]*T[0][1]+(L2-T[0][0])*(L2-T[0][0])));
	    	    
	    //principal_stress_direction_1.x = 0.0;
	    //principal_stress_direction_1.y = 0.0;
	    //principal_stress_direction_1.z = 0.0;
	    //principal_stress_direction_2.x = 0.0;
	    //principal_stress_direction_2.y = 0.0;
	    //principal_stress_direction_2.z = 0.0;
	    //pv_write_principal_stress_direction_1(ielem, principal_stress_direction_1);
            //pv_write_principal_stress_direction_2(ielem, principal_stress_direction_2);
	    
	    pv_write_element_connectivity(ielem, yde->i2elto[0][ielem], yde->i2elto[1][ielem], yde->i2elto[2][ielem]);
	    
	    //! Principal direction file
	    
	    if((abs(T[0][1])>EPSILON)&&(abs(T[1][0])>EPSILON))
	    { principal_stress_direction_1.x = T[0][1]/SQRT((T[0][1]*T[0][1]+(L1-T[0][0])*(L1-T[0][0])));
	      principal_stress_direction_1.y = (L1-T[0][0])/SQRT((T[0][1]*T[0][1]+(L1-T[0][0])*(L1-T[0][0])));
	      principal_stress_direction_2.x = (L2-T[1][1])/SQRT((T[0][1]*T[0][1]+(L2-T[1][1])*(L2-T[1][1])));
	      principal_stress_direction_2.y = T[0][1]/SQRT((T[0][1]*T[0][1]+(L2-T[1][1])*(L2-T[1][1])));
	    }
	    else
	    {
	      principal_stress_direction_1.x = 1.0;
	      principal_stress_direction_1.y = 0.0;
	      principal_stress_direction_2.x = 0.0;
	      principal_stress_direction_2.y = 1.0;
	    }
	    element_centroid_coord.x = (ydn->d2ncc[0][(yde->i2elto[0][ielem])]+ydn->d2ncc[0][(yde->i2elto[1][ielem])]+ydn->d2ncc[0][(yde->i2elto[2][ielem])])/3.0;
	    element_centroid_coord.y = (ydn->d2ncc[1][(yde->i2elto[0][ielem])]+ydn->d2ncc[1][(yde->i2elto[1][ielem])]+ydn->d2ncc[1][(yde->i2elto[2][ielem])])/3.0;
	    element_centroid_coord.z = 0.0;
	    pv_write_element_centroid_coordinate(ielem, element_centroid_coord);
	    pv_write_principal_stress_direction_1_andrea(ielem, principal_stress_direction_1);
	    pv_write_principal_stress_direction_2_andrea(ielem, principal_stress_direction_2);
	    
           }
        }
      }
    }
    
    //! NEW joint output
    
    //! Determine the total number of boundary joints
    //! Determine the total number of broken joints (including DFN broken joints)
    //! Determine the total number of yielded joints (not including broken joints)
    

    for(jprop=0;jprop<ydpj->npjset;jprop++) 
    { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
      { for(ielem=0;ielem<yde->nelem;ielem++)
        { ipropc=yde->i1elpr[ielem];
          if(ipropc==jprop+ydpe->nprop)
          { if((yde->i2elto[0][ielem]==yde->i2elto[3][ielem])&&(yde->i2elto[1][ielem]==yde->i2elto[2][ielem]))
            { boundary_joints++;
	      //ipropc=ipropc-YIPROPMAX;
            }
            if ((yde->i1esftf[ielem]>0) && (yde->d1ebrkf[ielem]<0.1))
            { yielded_joints++;
            }
          }
          if((ipropc<0)&&((ipropc+YIPROPMAX)==(jprop+ydpe->nprop)))
          { broken_joints++; 
          }
        }
      }
    }
    
    //! Determine the total number of acoustic events
    for(ielem=0;ielem<yde->nelem;ielem++)
    { if((yde->d1ebrkf[ielem]>0)&&(yde->d1edke[ielem]>0))
      { acoustic_events++; 
    } }
    
    //pv_set_num_broken_joints(boundary_joints+(yde->nebrk));
    pv_set_num_broken_joints(broken_joints+boundary_joints);
    //pv_set_num_yielded_joints(yde->nesft);
    //pv_set_num_yielded_joints(yde->nesft);
    pv_set_num_yielded_joints(yielded_joints);
    //pv_set_num_acoustic_events(yde->nebrk); //! = Broken joints - boundary joints
    pv_set_num_acoustic_events(acoustic_events); //! = Broken joints with acoustic energy greater than zero
    //printf("Broken_joints = %d\n",broken_joints);
    //printf("Boundary_joints = %d\n",boundary_joints);
    //printf("Softened_joints = %d\n", yde->nesft);
    
    total_fracture_area = 0.0; //! Reset area before counting
    total_fracture_length = 0.0; //! Reset length before counting
    
    for(jprop=0;jprop<ydpj->npjset;jprop++) //! Loop over joints
    { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
      { for(ielem=0;ielem<yde->nelem;ielem++)
        { ipropc=yde->i1elpr[ielem];
	  //! Modify the property of boundary joints so they are outputted as broken joints
	  if(ipropc==jprop+ydpe->nprop)
          { if((yde->i2elto[0][ielem]==yde->i2elto[3][ielem])&&(yde->i2elto[1][ielem]==yde->i2elto[2][ielem]))
            { ipropc=ipropc-YIPROPMAX; }
	  }
	  if((ipropc<0)&&((ipropc+YIPROPMAX)==(jprop+ydpe->nprop))) //! Broken joint output
          { 
	    nodal_coord_0.x = ydn->d2ncc[0][yde->i2elto[0][ielem]];
            nodal_coord_0.y = ydn->d2ncc[1][yde->i2elto[0][ielem]];
            nodal_coord_0.z = 0.0;
            
            nodal_coord_1.x = ydn->d2ncc[0][yde->i2elto[1][ielem]];
            nodal_coord_1.y = ydn->d2ncc[1][yde->i2elto[1][ielem]];
            nodal_coord_1.z = 0.0;
            
            nodal_coord_2.x = ydn->d2ncc[0][yde->i2elto[2][ielem]];
            nodal_coord_2.y = ydn->d2ncc[1][yde->i2elto[2][ielem]];
            nodal_coord_2.z = 0.0;
            
            nodal_coord_3.x = ydn->d2ncc[0][yde->i2elto[3][ielem]];
            nodal_coord_3.y = ydn->d2ncc[1][yde->i2elto[3][ielem]];
            nodal_coord_3.z = 0.0;
    
	    pv_write_broken_joint_coordinates(broken_joint_id, nodal_coord_0, nodal_coord_1, nodal_coord_2, nodal_coord_3);
	    
	     //! Joint mode of failure
            /*if (yde->i1ebrkf[ielem]==0)
            { joint_failure_mode = 0; } //! Joint is either a boundary joint or a DFN joint
            else if (yde->i1ebrkf[ielem]==2)
            { joint_failure_mode = 1; } //! Opening mode
            else if (yde->i1ebrkf[ielem]==3)
            { joint_failure_mode = 2; } //! Sliding mode
            else if(yde->i1ebrkf[ielem]==1)
            { joint_failure_mode = 3; } //! Opening + sliding mode
            else if(yde->i1ebrkf[ielem]==4)
            { joint_failure_mode = 4; } //! Excavated joint
            
            pv_write_broken_joint_mode_of_failure(broken_joint_id, joint_failure_mode);*/
	    
	    pv_write_broken_joint_mode_of_failure(broken_joint_id, yde->d1ebrkf[ielem]);
	    
	    //! Joint opening, sliding, and area
	    e1x=RP5*(nodal_coord_1.x+nodal_coord_2.x-nodal_coord_0.x-nodal_coord_3.x);
	    e1y=RP5*(nodal_coord_1.y+nodal_coord_2.y-nodal_coord_0.y-nodal_coord_3.y);
	    h=SQRT(e1x*e1x+e1y*e1y);
            e1x=e1x/(h+small);
            e1y=e1y/(h+small);
	    s1=(nodal_coord_0.y-nodal_coord_3.y)*e1y+(nodal_coord_0.x-nodal_coord_3.x)*e1x;
            s2=(nodal_coord_1.y-nodal_coord_2.y)*e1y+(nodal_coord_1.x-nodal_coord_2.x)*e1x;
            o1=(nodal_coord_0.y-nodal_coord_3.y)*e1x-(nodal_coord_0.x-nodal_coord_3.x)*e1y;
            o2=(nodal_coord_1.y-nodal_coord_2.y)*e1x-(nodal_coord_1.x-nodal_coord_2.x)*e1y;
	    
	    sliding=ABS((s1+s2)/2.0);
	    
	    opening=(o1+o2)/2.0;
	    if(opening<0)
	    { opening = 0; } //! If triangles are penetrating into each other, assume opening is zero (rather than negative)
	    	    
	    area = 0.5 * (((nodal_coord_1.x-nodal_coord_3.x)*(nodal_coord_0.y-nodal_coord_2.y))-((nodal_coord_0.x-nodal_coord_2.x)*(nodal_coord_1.y-nodal_coord_3.y)));
	    if(area < 0)
            { area = 0; } //! If triangles are penetrating into each other, assume area is zero (rather than negative)
	    
	    length = SQRT((nodal_coord_1.x-nodal_coord_0.x)*(nodal_coord_1.x-nodal_coord_0.x)+(nodal_coord_1.y-nodal_coord_0.y)*(nodal_coord_1.y-nodal_coord_0.y));
	    
	    //! If joint belongs to an excavation boundary, set sliding, opening, area and length to zero.
            if((ydpe->i1pexc[yde->i1elpr[yde->i2elnext[0][ielem]]]==1)||(ydpe->i1pexc[yde->i1elpr[yde->i2elnext[1][ielem]]]==1))
            { sliding = 0.0;
              opening = 0.0;
              area = 0.0;
	      length = 0.0;
            } 
            // If joint belongs to an excavation boundary, DFN or model boundary, length is zero
            if(yde->d1ebrkf[ielem]<EPSILON)
	    { length = 0.0; }
                        
            total_fracture_area += area;
	    total_fracture_length += length;

	    pv_write_broken_joint_sliding(broken_joint_id, sliding);
	    pv_write_broken_joint_opening(broken_joint_id, opening);
	    pv_write_broken_joint_area(broken_joint_id, area);
	    
	    broken_joint_id++;
    
	    //! Acoustic Emission Output
	    
	    if((yde->d1ebrkf[ielem]>0)&&(yde->d1edke[ielem]>0)) //! Joint is actually broken (i.e. it's not a boundary joint) (Acoustic Emission)
	    {
	     AE_coord.x = yde->d2ecbrk_NEW[0][ielem]; 
             AE_coord.y = yde->d2ecbrk_NEW[1][ielem]; 
             AE_coord.z = 0.0;
             pv_write_acoustic_coordinate(AE_id, AE_coord);
	     
	     //! AE parameters
	     event_energy = yde->d1edke[ielem];
             pv_write_acoustic_energy(AE_id, event_energy);
	     
	     /*if (event_energy > 0)
	     { event_magnitude = (log10(event_energy*exp(-9)) - 4.8)/1.5; }//! Gutenberg-Richter relationship, NOTE: units: mm, ms, microg
             else
	     { event_magnitude = 10; }*/
	     event_magnitude = (log10(event_energy*exp(-9)) - 4.8)/1.5; //! Gutenberg-Richter relationship, NOTE: units: mm, ms, microg
             pv_write_acoustic_magnitude(AE_id, event_magnitude);
	     
             fracture_energy = yde->d1efe[ielem];
	     pv_write_acoustic_fracture_energy(AE_id, fracture_energy);
	     
	     /*if (yde->i1ebrkf[ielem]==2)
             { mode = 1; } //! Opening mode
             else if (yde->i1ebrkf[ielem]==3)
             { mode = 2; } //! Sliding mode
             else if(yde->i1ebrkf[ielem]==1)
             { mode = 3; } //! Opening + sliding mode*/
             pv_write_acoustic_event_mode(AE_id, yde->d1ebrkf[ielem]);
	     
	     AE_id++;
	    } 
	  }
	  //! TODO Change the following condition so that yielded elements get updated even after they break! Change also pv_set_num_yielded_joints(yde->nesft);
	  if((yde->i1elpr[ielem]==(jprop+(ydpe->nprop)))&&(yde->i1esftf[ielem]>0))//&&(yde->d1ebrkf[ielem]<0.1))/*||(yde->i1ebrkf[ielem]>0)*/ //! Yielded joint output
          { 
	   nodal_coord_0.x = ydn->d2ncc[0][yde->i2elto[0][ielem]];
           nodal_coord_0.y = ydn->d2ncc[1][yde->i2elto[0][ielem]];
           nodal_coord_0.z = 0.0;
           
           nodal_coord_1.x = ydn->d2ncc[0][yde->i2elto[1][ielem]];
           nodal_coord_1.y = ydn->d2ncc[1][yde->i2elto[1][ielem]];
           nodal_coord_1.z = 0.0;
            
           nodal_coord_2.x = ydn->d2ncc[0][yde->i2elto[2][ielem]];
           nodal_coord_2.y = ydn->d2ncc[1][yde->i2elto[2][ielem]];
           nodal_coord_2.z = 0.0;
            
           nodal_coord_3.x = ydn->d2ncc[0][yde->i2elto[3][ielem]];
           nodal_coord_3.y = ydn->d2ncc[1][yde->i2elto[3][ielem]];
           nodal_coord_3.z = 0.0;  
	   
	   pv_write_yielded_joint_coordinates(yielded_joint_id, nodal_coord_0, nodal_coord_1, nodal_coord_2, nodal_coord_3);
	   
	   //! Joint mode of yielding
           if (yde->i1esftf[ielem]==0)
           { joint_yielding_mode = 0; } //! Joint still elastic
           else if (yde->i1esftf[ielem]==2)
           { joint_yielding_mode = 1; } //! Opening mode
           else if(yde->i1esftf[ielem]==3)
           { joint_yielding_mode = 2; } //! Sliding mode
           else if(yde->i1esftf[ielem]==1)
           { joint_yielding_mode = 3; } //! Opening + sliding mode
           
           pv_write_yielded_joint_mode_of_failure(yielded_joint_id, joint_yielding_mode);
	   
	   yielded_joint_id++;
          }  
        }
       }
     }
     
     pv_set_num_rebars(num_rebars); //! Total number of reference points divided by 2
     
     //printf("ydr->nrepo: %d\n",ydr->nrepo);
     //! Loop over rebars
     for(ip=0;ip<(ydr->nrepo);ip=ip+2)
     { //printf("ip: %d\n",ip);
       //printf("rebar_id: %d\n",rebar_id);
       //printf("rebar_id: %d\n",rebar_id);
       nodal_coord_0.x = ydr->d2rccg[0][ip];
       nodal_coord_0.y = ydr->d2rccg[1][ip];
       nodal_coord_0.z = 0.0;
       nodal_coord_1.x = ydr->d2rccg[0][ip+1];
       nodal_coord_1.y = ydr->d2rccg[1][ip+1];
       nodal_coord_1.z = 0.0;
       //printf("nodal_coord_0.x: %f \t nodal_coord_0.y: %f\n",nodal_coord_0.x,nodal_coord_0.y);
       //printf("nodal_coord_1.x: %f \t nodal_coord_1.y: %f\n",nodal_coord_1.x,nodal_coord_1.y);
       pv_write_rebar_coordinates(rebar_id, nodal_coord_0, nodal_coord_1);
       
       rebar_area = ydsb->d1spea[ydr->i1refbar[ip]];
       rebar_stress = ydr->d1rbsig[ip];
       rebar_force = ydr->d1rbfrc[ip];
       rebar_strain = ydr->d1rbstr[ip];
       rebar_activation_flag = ydsb->i1sbac[ydr->i1refbar[ip]];

       pv_write_rebar_area(rebar_id, rebar_area);
       pv_write_rebar_stress(rebar_id, rebar_stress);
       pv_write_rebar_force(rebar_id, rebar_force);
       pv_write_rebar_strain(rebar_id, rebar_strain);
       pv_write_rebar_activation_flag(rebar_id, rebar_activation_flag);
       
       rebar_id++;
     }
     
     
     //! Determine the total number of 1D joints that are not broken
     for (iprop=0;iprop<ydps->nprop;iprop++)
     { for(i=0;ydr->i2relto[0][i]>=0; i++) 
       { if(ydr->i1rprop[ydr->i2relto[0][i]]==iprop && ydr->i1rprop[ydr->i2relto[1][i]]==iprop && ydr->i1myjoint[i]>=R0)
         { unbroken_1D_joints++;
     } } }
     
     //printf("unbroken_1D_joints: %d\n", unbroken_1D_joints);
     
     pv_set_num_unbroken_1D_joints(unbroken_1D_joints); 
     
     //! Loop over 1D joints
     for (iprop=0;iprop<ydps->nprop;iprop++)
     { for(i=0;ydr->i2relto[0][i]>=0; i++) 
       { if(ydr->i1rprop[ydr->i2relto[0][i]]==iprop && ydr->i1rprop[ydr->i2relto[1][i]]==iprop && ydr->i1myjoint[i]>=R0)
         { i0r=ydr->i2relto[0][i];
           i1r=ydr->i2relto[1][i];
         
           nodal_coord_0.x = ydr->d2rccg[0][i0r];
           nodal_coord_0.y = ydr->d2rccg[1][i0r];
           nodal_coord_0.z = 0.0;
           nodal_coord_1.x = ydr->d2rccg[0][i1r];
           nodal_coord_1.y = ydr->d2rccg[1][i1r];
           nodal_coord_1.z = 0.0;

           //oneD_joint_stress = 5e6;
           
           oneD_joint_normal_stress = ydr->d1rjsig[i];
           oneD_joint_tangential_stress = ydr->d1rjtau[i];
           oneD_joint_normal_strain = ydr->d1rjslpnor[i]; //! This is actually a displacement
           oneD_joint_tangential_strain = ydr->d1rjdelnor[i]; //! This is actually a displacement
           oneD_joint_state = ydr->i1rjstnor[i]; //! Valid only if I1SBTY  5, 6 or 7
           oneD_joint_force.x = ydr->d2rjfrv[0][i];
           oneD_joint_force.y = ydr->d2rjfrv[1][i];
           oneD_joint_force.z = 0.0;

           //printf("i: %d\n", i);
           //printf("i0r: %d \t i1r: %d\n", i0r, i1r);
           //printf("nodal_coord_0.x: %f \t nodal_coord_0.y: %f\n",nodal_coord_0.x,nodal_coord_0.y);
           //printf("nodal_coord_1.x: %f \t nodal_coord_1.y: %f\n",nodal_coord_1.x,nodal_coord_1.y);
   
           pv_write_1D_joint_coordinates(oneD_joint_id, nodal_coord_0, nodal_coord_1);
           pv_write_1D_joint_normal_stress(oneD_joint_id, oneD_joint_normal_stress);
           pv_write_1D_joint_tangential_stress(oneD_joint_id, oneD_joint_tangential_stress);
           pv_write_1D_joint_normal_strain(oneD_joint_id, oneD_joint_normal_strain);
           pv_write_1D_joint_tangential_strain(oneD_joint_id, oneD_joint_tangential_strain);
           pv_write_1D_joint_state(oneD_joint_id, oneD_joint_state);
           //pv_write_1D_joint_force(oneD_joint_id, oneD_joint_force); 
         
           oneD_joint_id++;
         }
       }
     }
    
    
    //! OLD joint output
    //! Joint nodal coordinates, joint state, connectivity table
    /*for(jprop=0;jprop<ydpj->npjset;jprop++) //! Loop over joints
    { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
      { for(ielem=0;ielem<yde->nelem;ielem++)
        { if(yde->i1elpr[ielem]==(jprop+(ydpe->nprop))) //! Elastic or yielded joint
          { 
	    nodal_coord.x = ydn->d2ncc[0][yde->i2elto[0][ielem]];
            nodal_coord.y = ydn->d2ncc[1][yde->i2elto[0][ielem]];
            nodal_coord.z = 0.0;
            pv_write_nodal_current_coord(yde->i2elto[0][ielem], nodal_coord);
            
            nodal_coord.x = ydn->d2ncc[0][yde->i2elto[1][ielem]];
            nodal_coord.y = ydn->d2ncc[1][yde->i2elto[1][ielem]];
            nodal_coord.z = 0.0;
            pv_write_nodal_current_coord(yde->i2elto[1][ielem], nodal_coord);
            
            nodal_coord.x = ydn->d2ncc[0][yde->i2elto[2][ielem]];
            nodal_coord.y = ydn->d2ncc[1][yde->i2elto[2][ielem]];
            nodal_coord.z = 0.0;
            pv_write_nodal_current_coord(yde->i2elto[2][ielem], nodal_coord);
            
            nodal_coord.x = ydn->d2ncc[0][yde->i2elto[3][ielem]];
            nodal_coord.y = ydn->d2ncc[1][yde->i2elto[3][ielem]];
            nodal_coord.z = 0.0;
            pv_write_nodal_current_coord(yde->i2elto[3][ielem], nodal_coord);
  
	    //! Joint state
            if (yde->i1esftf[ielem]==0)
            { joint_state = 0; } //! Elastic joint
            else if (yde->i1ebrkf[ielem]>0) 
            { joint_state =  2; } //! Broken joint
            else 
            { joint_state = 1; } //! Softened joint
            pv_write_joint_state(ielem, joint_state);
            
	    //! Joint connectivity
	    pv_write_joint_connectivity(ielem, yde->i2elto[0][ielem], yde->i2elto[1][ielem], yde->i2elto[2][ielem], yde->i2elto[3][ielem]);
	    //!pv_write_joint_connectivity(ielem, yde->i2elto[0][ielem], yde->i2elto[1][ielem], yde->i2elto[1][ielem], yde->i2elto[0][ielem]);
	    
            //! Joint mode of yielding
            //if (yde->i1esftf[ielem]==0)
            //{ joint_yielding_mode = 0; } //! Joint still elastic
            //else if (yde->i1esftf[ielem]==2)
            //{ joint_yielding_mode = 1; } //! Opening mode
            //else if(yde->i1esftf[ielem]==3)
            //{ joint_yielding_mode = 2; } //! Sliding mode
            //else if(yde->i1esftf[ielem]==1)
            //{ joint_yielding_mode = 3; } //! Opening + sliding mode
            //pv_write_joint_mode_of_yielding(ielem, joint_yielding_mode);
            
	  }
         
          for(jprop=0;jprop<ydpj->npjset;jprop++) //! Loop over joints
          { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
            { for(ielem=0;ielem<yde->nelem;ielem++)
              { ipropc=yde->i1elpr[ielem];
                if(ipropc==iprop)
                { if((yde->i2elto[0][ielem]==yde->i2elto[3][ielem])&&(yde->i2elto[1][ielem]==yde->i2elto[2][ielem]))
                  { ipropc=ipropc-YIPROPMAX;
                  } 
	        }
	      }
            }
          }
	  
	  if((ipropc<0)&&((ipropc+YIPROPMAX)==(jprop+ydpe->nprop))) //! Broken joint output
          {
	    //! Joint state
            if (yde->i1esftf[ielem]==0)
            { joint_state = 0; } //! Elastic joint
            else if (yde->i1ebrkf[ielem]>0) 
            { joint_state =  2; } //! Broken joint
            else 
            { joint_state = 1; } //! Softened joint
            pv_write_joint_state(ielem, joint_state);
	    //! Joint mode of failure
            if (yde->i1ebrkf[ielem]==0)
            { joint_failure_mode = 0; } //! Joint is not broken yet
            else if (yde->i1ebrkf[ielem]==2)
            { joint_failure_mode = 1; } //! Opening mode
            else if (yde->i1ebrkf[ielem]==3)
            { joint_failure_mode = 2; } //! Sliding mode
            else if(yde->i1ebrkf[ielem]==1)
            { joint_failure_mode = 3; } //! Opening + sliding mode
            pv_write_joint_mode_of_failure(ielem, joint_failure_mode);
	    //!Joint breakage kinetic energy (i.e. acoustic emission)
	    event_energy = yde->d1edke[ielem];
            pv_write_joint_breakage_kinetic_energy(ielem, event_energy);
	  }
        }
      }
    }*/
    pv_write_timestep(ydc->ncstep); //! timestep output
    
    /* Output inputfilename.yhFRAC */
    ncallFRAC=ncallFRAC+1;
    if((ydo->fofrac)==FILENULL)
    { CHRcpy(namef,namep);
      CHRcat(namef,"hFRAC");
      ydo->fofrac=fopen(namef,"a");
    }
    if((ydo->fofrac)!=FILENULL)
    { fprintf(ydo->fofrac,"%d",ydc->ncstep); /* time step number */
      CHRwtab(ydo->fofrac);
      fprintf(ydo->fofrac,"%.9f",total_fracture_area); /* total fracture area */
      CHRwtab(ydo->fofrac);
      fprintf(ydo->fofrac,"%.9f",total_fracture_length); /* total fracture length */
      CHRwcr(ydo->fofrac);
    } 
    if((ncallFRAC>100)||(ydc->ncstep)>=(ydc->mcstep-2))
    { ncallFRAC=0;
      if((ydo->fofrac)!=FILENULL)fclose(ydo->fofrac);
      ydo->fofrac=FILENULL;
    }
    fflush(ydo->fofrac);
  }
} 


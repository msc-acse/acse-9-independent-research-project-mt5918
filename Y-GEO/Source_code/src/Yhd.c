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
/* File   Ymd.c */
#include "Yproto.h"
#include <assert.h>
#include <math.h>       /* log */

static void YhdINIT(  //! Create i2noid 
            nelem, iprop,
            d1nccx,d1nccy,
            i1elpr,i2elto,
            i1nowe,i2noid,
            i2elnext
            )
INT   nelem; INT   iprop;
DBL *d1nccx; DBL *d1nccy; 
INT *i1elpr; INT **i2elto; 
INT *i1nowe; INT **i2noid;
INT **i2elnext; 

{ INT ielem;
  
  //! Loop over joint elements
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    {
      /*
      printf("Hello, you are on joint element %d\n",ielem);
      printf("x: %f, y: %f\n", d1nccx[i2elto[0][ielem]],d1nccy[i2elto[0][ielem]]);
      printf("x: %f, y: %f\n", d1nccx[i2elto[1][ielem]],d1nccy[i2elto[1][ielem]]);
      printf("x: %f, y: %f\n", d1nccx[i2elto[2][ielem]],d1nccy[i2elto[2][ielem]]);
      printf("x: %f, y: %f\n", d1nccx[i2elto[3][ielem]],d1nccy[i2elto[3][ielem]]);
      */
      
      assert(i2noid[0][i2elto[0][ielem]] == -1);
      assert(i2noid[1][i2elto[1][ielem]] == -1);

      // Check whether the joint is on the exterior
      if(i2elto[0][ielem]==i2elto[3][ielem])
      {
        assert(i2elto[1][ielem]==i2elto[2][ielem]);
        i2noid[0][i2elto[0][ielem]] = i2elto[1][ielem];
        i2noid[1][i2elto[1][ielem]] = i2elto[0][ielem];
      }
      else
      {
        assert(i2elto[1][ielem]!=i2elto[2][ielem]);
        assert(i2noid[0][i2elto[2][ielem]] == -1);
        assert(i2noid[1][i2elto[3][ielem]] == -1);
        i2noid[0][i2elto[0][ielem]] = i2elto[3][ielem];
        i2noid[1][i2elto[1][ielem]] = i2elto[2][ielem];
        i2noid[0][i2elto[2][ielem]] = i2elto[1][ielem];
        i2noid[1][i2elto[3][ielem]] = i2elto[0][ielem];
      }
    }

} }

static void YhdWE(  //! Propagate wet boundary on nodes
            nnopo,ncstep,ibc,
            d1nmct,i1nopr,i1nowe,i2noid
            )
INT nnopo; INT ncstep; INT ibc;
DBL *d1nmct; INT *i1nopr; INT *i1nowe; INT **i2noid;

// { INT inopo; INT jnopo; INT wet;
// 
//   wet = 2 - (ncstep % 2);
// 
//   //! Loop over nodes
//   for(inopo=0;inopo<nnopo;inopo++) 
//   { if((i1nopr[inopo]==ibc)&&(d1nmct[inopo]>EPSILON))
//     { if((i1nowe[inopo]!=wet) && (i1nowe[inopo]>0) && (i1nowe[inopo]<3))
//       {
//         jnopo = inopo;
//         do {
//           //printf("Node idx: %d, wet: %d\n", jnopo,i1nowe[jnopo]);
//           assert(i1nowe[jnopo] != wet);
//           i1nowe[jnopo] = wet;
//           jnopo = i2noid[0][jnopo];
//         } while(jnopo!=inopo);
//       }
//     }
//   }
// }

{ INT inopo; INT jnopo; INT wet; INT old_wet;

  wet = 2 - (ncstep % 2);
  old_wet = 3 - wet;

  //! Loop over nodes
  for(inopo=0;inopo<nnopo;inopo++)
  { if((i1nopr[inopo]==ibc)&&(d1nmct[inopo]>EPSILON))
    { if(i1nowe[inopo]==old_wet)
      { jnopo = inopo;
        while(i1nowe[jnopo]==old_wet || i1nowe[jnopo]==0)
        //while(i1nowe[jnopo]==old_wet)
        { i1nowe[jnopo] = wet;
          jnopo = i2noid[0][jnopo];
        }
} } } }




static void YhdFLUVOL(  //! Calculate volume of fluid cavity 
            nnopo,ibc,
            d1nmct,i1nowe,i2noid,
            d1nccx,d1nccy,
            i1nopr,fluvol,
            dim
            )
INT nnopo; INT ibc; DBL *d1nmct;
INT *i1nowe; INT **i2noid;
DBL *d1nccx; DBL *d1nccy; 
INT *i1nopr; DBL *fluvol;
INT dim;

{ INT inopo; DBL x0; DBL x1; DBL y0; DBL y1; DBL A; DBL l; DBL r0; DBL r1; DBL lx; DBL ly;

  //! Loop over nodes
  for(inopo=0;inopo<nnopo;inopo++) 
  { if((i1nopr[inopo]==ibc)&&(d1nmct[inopo]>EPSILON))
    {
      if(i1nowe[inopo]==1 || i1nowe[inopo]==2)
      {
        if(dim == 2)
        { *fluvol += 0.5 * (d1nccx[i2noid[0][inopo]]*d1nccy[inopo] - d1nccx[inopo]*d1nccy[i2noid[0][inopo]]);
        }
        else if(dim == 3)
        { x0 = d1nccx[inopo]; x1 = d1nccx[i2noid[0][inopo]];
          y0 = d1nccy[inopo]; y1 = d1nccy[i2noid[0][inopo]];
          A = x1*y0 - x0*y1; //2x signed area
          if(fabs(A) > EPSILON)
          { lx = x1 - x0; ly = y1 - y0;
            l  = sqrt(lx*lx + ly*ly);
            r0 = sqrt(x0*x0 + y0*y0);
            r1 = sqrt(x1*x1 + y1*y1);
            *fluvol += A / (3 * l * l * l) * ( A * A * log( (x1*lx + y1*ly + l*r1) / (x0*lx + y0*ly + l*r0) ) + l * (r1 * (x1*lx + y1*ly) - r0 * (x0*lx + y0*ly)));
            assert(*fluvol==*fluvol);
} } } } } }

static void  YhdPUMP( //! Integrate pump equation
             fluvol, flupres, flumass, flurho0, flupres0, flubulk, dcstec,
             ihftyp,dhfflp,dhfflq,hfarow,d2hfaf,dctime,
             d1nfp,gravacc,nnopo, d1nccx, d1nccy, d2wtlev
             )
DBL fluvol; DBL *flupres; DBL *flumass; DBL flurho0; DBL flupres0; DBL flubulk; DBL dcstec;
INT ihftyp; DBL dhfflp; DBL dhfflq; INT hfarow; DBL **d2hfaf; DBL dctime;
DBL *d1nfp; DBL gravacc; INT nnopo; DBL *d1nccx; DBL *d1nccy; DBL **d2wtlev;

{ INT i, inopo;
  DBL water_y_coord;
 
  if(ihftyp==1) //! pressure
  { for(i=1; i<hfarow; i++)
    { if((dctime>=d2hfaf[0][i-1])&&(dctime<=d2hfaf[0][i]))
      {
        *flupres = dhfflp * (d2hfaf[1][i-1] + (d2hfaf[1][i]-d2hfaf[1][i-1]) * (dctime-d2hfaf[0][i-1])/(d2hfaf[0][i]-d2hfaf[0][i-1]));
        return;
    } }
    if(hfarow==0)
    { *flupres = dhfflp;
    }
    else if(dctime<=d2hfaf[0][0])
    { *flupres = dhfflp * d2hfaf[1][0];
    }
    else if(dctime>=d2hfaf[0][hfarow-1])
    { *flupres = dhfflp * d2hfaf[1][hfarow-1];
    }
  }
  else if(ihftyp==2) //! flow rate
  { for(i=1; i<hfarow; i++)
    { if((dctime>=d2hfaf[0][i-1])&&(dctime<=d2hfaf[0][i]))
      {
        *flumass += dhfflq * dcstec * (d2hfaf[1][i-1] + (d2hfaf[1][i]-d2hfaf[1][i-1]) * (dctime-d2hfaf[0][i-1])/(d2hfaf[0][i]-d2hfaf[0][i-1]));
        //return;
    } }
    if(hfarow==0)
    { *flumass += dhfflq * dcstec;
    }
    
    else if(dctime<=d2hfaf[0][0])
    { *flumass += dhfflq * dcstec * d2hfaf[1][0];
    }
    else if(dctime>=d2hfaf[0][hfarow-1])
    { //assert(dctime>=d2hfaf[hfarow][0]);
      *flumass += dhfflq * dcstec * d2hfaf[1][hfarow-1];
    }
    
    *flupres = flupres0 + flubulk * log(*flumass / (fluvol * flurho0));
  }
  else if(ihftyp==3) //! hydrostatic pressure
  { for(inopo=0;inopo<nnopo;inopo++) //! Loop over nodes
    { //! Determine y-coordinate of water level for the given node
      if(d1nccx[inopo]<=d2wtlev[0][0]) // Upstream level
      { water_y_coord = d2wtlev[1][0]; }
      else if((d1nccx[inopo]>d2wtlev[0][0])&&(d1nccx[inopo]<=d2wtlev[0][1])) // Dam: linear interpolation between upstream and downstream level
      { water_y_coord = d2wtlev[1][0] + (d1nccx[inopo]-d2wtlev[0][0]) * ((d2wtlev[1][0]-d2wtlev[1][1])/(d2wtlev[0][0]-d2wtlev[0][1])); }
      else if(d1nccx[inopo]>d2wtlev[0][1]) // Downstream level
      { water_y_coord = d2wtlev[1][1]; }
      
      //! Calculate time-varying nodal pressure
      for(i=1; i<hfarow; i++)
      { if((dctime>=d2hfaf[0][i-1])&&(dctime<=d2hfaf[0][i]))
        { if(d1nccy[inopo] <= water_y_coord) // Apply fluid pressure only to nodes located below the water level
          { d1nfp[inopo] = ((water_y_coord-d1nccy[inopo]) * flurho0 * gravacc) * (d2hfaf[1][i-1] + (d2hfaf[1][i]-d2hfaf[1][i-1]) * (dctime-d2hfaf[0][i-1])/(d2hfaf[0][i]-d2hfaf[0][i-1])); }
          else
          { d1nfp[inopo] = 0.0; }
        } 
      }
      if(hfarow==0)
      { if(d1nccy[inopo] <= water_y_coord) // Apply fluid pressure only to nodes located below the water level
        { d1nfp[inopo] = (water_y_coord-d1nccy[inopo]) * flurho0 * gravacc; }
        else
        { d1nfp[inopo] = 0.0; }
      }
      else if(dctime<=d2hfaf[0][0])
      { if(d1nccy[inopo] <= water_y_coord) // Apply fluid pressure only to nodes located below the water level
        { d1nfp[inopo] = ((water_y_coord-d1nccy[inopo]) * flurho0 * gravacc) * d2hfaf[1][0]; }
        else
        { d1nfp[inopo] = 0.0; }
      }
      else if(dctime>=d2hfaf[0][hfarow-1])
      { if(d1nccy[inopo] <= water_y_coord) // Apply fluid pressure only to nodes located below the water level
        { d1nfp[inopo] = ((water_y_coord-d1nccy[inopo]) * flurho0 * gravacc) * d2hfaf[1][hfarow-1]; }
        else
        { d1nfp[inopo] = 0.0; }
      }
     }
   }
}

static void YhdMASS(  //! Calculate the mass for a given pressure and volume */
            fluvol, flupres, flumass, flurho0, flupres0, flubulk
            )
DBL fluvol; DBL flupres; DBL *flumass; DBL flurho0; DBL flupres0; DBL flubulk;
{
  *flumass = fluvol * flurho0 * exp((flupres - flupres0) / flubulk);
}

static void YhdFLUPR(  //! Update nodal forces due to fluid pressure */
            nnopo,ibc,d1nmct,
            i1nowe,i2noid,
            d1nccx,d1nccy,
            i1nopr,
            d1nfcx,d1nfcy,
            flupres,d1nfp,ihftyp
            )
INT nnopo; INT ibc; DBL *d1nmct;
INT *i1nowe; INT **i2noid;
DBL *d1nccx; DBL *d1nccy; 
INT *i1nopr;
DBL *d1nfcx; DBL *d1nfcy;
DBL flupres; DBL *d1nfp;
INT ihftyp;

{ INT inopo;
 
  //! Loop over nodes
  for(inopo=0;inopo<nnopo;inopo++) 
  { if((i1nopr[inopo]==ibc)&&(d1nmct[inopo]>EPSILON))
    {
      if(i1nowe[inopo]==1 || i1nowe[inopo]==2)
      {
        if((ihftyp==1) || (ihftyp==2)) 
        { d1nfcx[inopo]            += 0.5 * flupres * (d1nccy[inopo] - d1nccy[i2noid[0][inopo]]);
          d1nfcx[i2noid[0][inopo]] += 0.5 * flupres * (d1nccy[inopo] - d1nccy[i2noid[0][inopo]]);

          d1nfcy[inopo]            += 0.5 * flupres * (d1nccx[i2noid[0][inopo]] - d1nccx[inopo]);
          d1nfcy[i2noid[0][inopo]] += 0.5 * flupres * (d1nccx[i2noid[0][inopo]] - d1nccx[inopo]);
        }
        else if(ihftyp==3)
        { d1nfcx[inopo]            += 0.5 * d1nfp[inopo]            * (d1nccy[inopo] - d1nccy[i2noid[0][inopo]]);
          d1nfcx[i2noid[0][inopo]] += 0.5 * d1nfp[i2noid[0][inopo]] * (d1nccy[inopo] - d1nccy[i2noid[0][inopo]]);

          d1nfcy[inopo]            += 0.5 * d1nfp[inopo]            * (d1nccx[i2noid[0][inopo]] - d1nccx[inopo]);
          d1nfcy[i2noid[0][inopo]] += 0.5 * d1nfp[i2noid[0][inopo]] * (d1nccx[i2noid[0][inopo]] - d1nccx[inopo]);
        }

      }
    }
  }
}


/*********************PUBLIC********************************************************/
void Yhd(ydc,  yde,  ydn,  ydo,  ydpe, ydpn, ydpj, ydis, ydhf, ydfn    /***  hydrofrac  ***/
        )
  YDC ydc; YDE yde; YDN ydn; YDO ydo; YDPE ydpe; YDPN ydpn; YDPJ ydpj; YDIS ydis; YDHF ydhf; YDFN ydfn;
{ INT iprop,jprop,i,j,ielem;
  INT ibc, inopo;
  
 if(ydhf->iusehf==1)
 {
  //! Initializing d1nfp */
  if(ydhf->ihftyp==3)
  { if(ydn->d1nfp==DBL1NULL)
    { ydn->d1nfp=TalDBL1(ydn->mnopo);
      for(i=0;i<ydn->mnopo;i++)
       ydn->d1nfp[i]=R0;
    }
  }
  
  //! Create i2noid only at the beginning of the simulation 
  if(ydc->ncstep==0) 
  {
    //! Initialize node connections for joints belonging to a DFN (with "broken" fractures)
    if(ydfn->iusefn == 1)
    { for(ielem=0;ielem<yde->nelem;ielem++)
      { if(yde->i1edfnf[ielem]==1)
        { ydn->i2noid[0][yde->i2elto[0][ielem]]=yde->i2elto[1][ielem];
          ydn->i2noid[1][yde->i2elto[1][ielem]]=yde->i2elto[0][ielem];
          ydn->i2noid[0][yde->i2elto[2][ielem]]=yde->i2elto[3][ielem];
          ydn->i2noid[1][yde->i2elto[3][ielem]]=yde->i2elto[2][ielem];
    } } } 
    //! Loop over joint property sets
    for(jprop=0;jprop<ydpj->npjset;jprop++)
    { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
      {
        YhdINIT( //! Initialize the node connections
        yde->nelem,(jprop+ydpe->nprop),
        ydn->d2ncc[0],ydn->d2ncc[1],
        yde->i1elpr,yde->i2elto,
        ydn->i1nowe,ydn->i2noid,
        yde->i2elnext
        );
  } } }

  //! Loop over boundary condition sets
  for(ibc=0;ibc<ydpn->npnset;ibc++)
  { YhdWE( //! Update i1nowe  
    ydn->nnopo, ydc->ncstep, ibc,
    ydn->d1nmct, ydn->i1nopr, ydn->i1nowe,ydn->i2noid
    );
  }
  
  ydhf->fluvol=0.0;
  //! Loop over boundary condition sets
  for(ibc=0;ibc<ydpn->npnset;ibc++)
  { YhdFLUVOL( //! Calculate volume of fluid
    ydn->nnopo,ibc,
    ydn->d1nmct, ydn->i1nowe,ydn->i2noid,
    ydn->d2ncc[0],ydn->d2ncc[1],
    ydn->i1nopr,&(ydhf->fluvol),
    ydhf->fradim
    );
  } 
  
  if((ydhf->ihftyp==2)&&(ydhf->ihfmsin==0)) 
  {
    //! Initialize the fluid mass based on the pressure and volume only if flow rate is used as hydrofrac input and the mass has not been initialized before
    YhdMASS(
    ydhf->fluvol, ydhf->flupres, &(ydhf->flumass), ydhf->flurho0, ydhf->flupres0, ydhf->flubulk);
    
    ydhf->ihfmsin = 1;
  }
  
  //! Integrate pump equation
  YhdPUMP(
  ydhf->fluvol, &(ydhf->flupres), &(ydhf->flumass), ydhf->flurho0, ydhf->flupres0, ydhf->flubulk,
  ydc->dcstec, ydhf->ihftyp,ydhf->dhfflp,ydhf->dhfflq,ydhf->hfarow,ydhf->d2hfaf,ydc->dctime,
  ydn->d1nfp,ydhf->gravacc,ydn->nnopo,ydn->d2ncc[0],ydn->d2ncc[1],ydhf->d2wtlev
  );
  
  //! Loop over boundary condition sets
  for(ibc=0;ibc<ydpn->npnset;ibc++)
  { YhdFLUPR( //! Update forces due to fluid pressure
    ydn->nnopo,ibc,ydn->d1nmct,
    ydn->i1nowe,ydn->i2noid,
    ydn->d2ncc[0],ydn->d2ncc[1],
    ydn->i1nopr,
    ydn->d2nfc[0],ydn->d2nfc[1],
    ydhf->flupres,ydn->d1nfp,ydhf->ihftyp
    );
  } 
 }
}




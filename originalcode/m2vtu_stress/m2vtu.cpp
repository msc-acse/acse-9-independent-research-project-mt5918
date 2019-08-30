/**********************************************************************/
/** Copyright (C) 2008,                                              **/
/** Queen Mary University of London (QMUL) & Imperial College        **/
/** of Science, Technology and Medicine (ICSTM). All rights reserved.**/
/** Implemented for you by Prof Antonio Munjiza & Dr Jiansheng Xiang **

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

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkZLibDataCompressor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkSystemIncludes.h>
#include "cubicsolver.h"

extern "C" {
#include "frame.h"
#include "Ytypes.h"
}

static void readppm(FILE *fdata, vtkPoints *pts, vtkCellArray *polys, 
             vtkDoubleArray *veloc_vector, vtkDoubleArray *stress_tensor,
			 vtkDoubleArray *sigma1_scalar, vtkDoubleArray *sigma3_scalar,
			 vtkDoubleArray *sigma1_vectoror, vtkDoubleArray *sigma3_vectoror,
			 vtkDoubleArray *sigma_diff_scalar, vtkDoubleArray *sigma_mean_scalar,
			 vtkDoubleArray *cforce_vector,
			 INT npoint, DBL dcsizc, DBL dcsizs, DBL dcsizf, DBL dcsizv)
{
  INT i1num[100],i,j;
  DBL d1num[100];
  DBL x[3],y[3],z[3],v[3],xc,yc,zc,v0[3];
  DBL xg,yg;
  DBL Ncoef;   /* Ncoef is the Normalisation factor */
  DBL **T;
  vtkIdType npts[3];
  CHR c1code[500000];
  DBL sigma1,sigma3,sigma_diff,sigma_mean;
  DBL *v1,*v2,*v3,a1,a2,a3;
  DBL cforce[3];

  T=TalDBL2(3,3);
  v0[0]=R0;
  v0[1]=R0;
  v0[2]=R0;
  j=0;
  v1=TalDBL1(3);
  v2=TalDBL1(3);
  v3=TalDBL1(3);
  a1=R0;
  a2=R0;
  a3=R0;
  cforce[0]=R0;
  cforce[1]=R0;
  cforce[2]=R0;

  if(fdata!=FILENULL)
  { CHRr(fdata,c1code); CHRr(fdata,c1code);
    while((FILEND(fdata)==0))
    { codeCHRtoINT(c1code,i1num);
      codeINTtoDBL(d1num,i1num);
      if(i1num[2]==YTE2TRIELS)
      { for(i=0;i<3;i++)
		{ x[i]=d1num[3+i]*dcsizc;
		  y[i]=d1num[6+i]*dcsizc;
		  z[i]=R0;
		  v[i]=d1num[9+i]*dcsizv;
		  pts->InsertNextPoint(x[i],y[i],z[i]);
		  npts[i]=npoint;
		  npoint++;
		}
        /* Center of mass of element */
        xg=(x[0]+x[1]+x[2])/R3;
        yg=(y[0]+y[1]+y[2])/R3;
        Ncoef=15.441*yg-4.9934;
        v[2]=R0;

        T[0][0]=-d1num[11]*dcsizs; T[0][1]=-d1num[13]*dcsizs; T[0][2]=R0;
        T[1][0]=-d1num[13]*dcsizs; T[1][1]=-d1num[12]*dcsizs; T[1][2]=R0;
		T[2][0]=R0; T[2][1]=R0; T[2][2]=R0;
        
		cforce[0]=d1num[15]*dcsizf; cforce[1]=d1num[16]*dcsizf; cforce[2]=R0;

		sigma1=RP5*(T[0][0]+T[1][1])                 /* SIGMA I    */
                +SQRT((RP5*(T[0][0]-T[1][1]))*(RP5*(T[0][0]-T[1][1]))
                +(T[0][1]*T[0][1]));
		sigma3=RP5*(T[0][0]+T[1][1])                /* SIGMA III   */
                -SQRT((RP5*(T[0][0]-T[1][1]))*(RP5*(T[0][0]-T[1][1]))
                +(T[0][1]*T[0][1]));
        sigma_diff=sigma1-sigma3;           /* Differential stress */
        sigma_mean=(sigma1+sigma3)/R2;
		
		solveSymetricEigenProblem(T,v1,v2,v3,&a1,&a2,&a3);

		for(i=0;i<3;i++)
		{ // veloc_vector->InsertNextTuple3(R0,R0,R0);
		  veloc_vector->InsertNextTuple3(v[0],v[1],v[2]);
		  stress_tensor->InsertNextTuple9(T[0][0],T[0][1],T[0][2],
										  T[1][0],T[1][1],T[1][2],
										  T[2][0],T[2][1],T[2][2]);
		  sigma1_scalar->InsertNextTuple1(sigma1);
		  sigma3_scalar->InsertNextTuple1(sigma3);
		  sigma1_vectoror->InsertNextTuple3(a3*v3[0],a3*v3[1],a3*v3[2]);
		  sigma3_vectoror->InsertNextTuple3(a1*v1[0],a1*v1[1],a1*v1[2]);
		  sigma_diff_scalar->InsertNextTuple1(sigma_diff);
		  sigma_mean_scalar->InsertNextTuple1(sigma_mean);
		  cforce_vector->InsertNextTuple3(cforce[0],cforce[1],cforce[2]);
		}

        polys->InsertNextCell(3,npts);
        xc=(x[0]+x[1]+x[2])/R3;
        yc=(y[0]+y[1]+y[2])/R3;
        zc=R0;
        pts->InsertNextPoint(xc,yc,zc);
        npoint++;

        veloc_vector->InsertNextTuple3(v[0],v[1],v[2]);
        stress_tensor->InsertNextTuple9(T[0][0],T[0][1],T[0][2],
										T[1][0],T[1][1],T[1][2],
										T[2][0],T[2][1],T[2][2]);
		sigma1_scalar->InsertNextTuple1(sigma1);
		sigma3_scalar->InsertNextTuple1(sigma3);
		sigma1_vectoror->InsertNextTuple3(a3*v3[0],a3*v3[1],a3*v3[2]);
		sigma3_vectoror->InsertNextTuple3(a1*v1[0],a1*v1[1],a1*v1[2]);
		sigma_diff_scalar->InsertNextTuple1(sigma_diff);
		sigma_mean_scalar->InsertNextTuple1(sigma_mean);
		cforce_vector->InsertNextTuple3(cforce[0],cforce[1],cforce[2]);
      }
      else if(i1num[2]==YTE2JOINTS)
      {// CHRw(stdout,"drawjoiint");
      }
      else if(i1num[2]==YTEG2RAD5)
      {// CHRw(stdout,"drawG2RP5");
      }
      CHRr(fdata,c1code);
      j=j+1;
} } }

int main(INT argc, char **argv)
{ 
  FILE *fdata=FILENULL;
  FILE *fcoeff=FILENULL;

  CHR name[300];
  CHR c1tmp[100],c1tmp1[100];
  CHR c1name[100];
  CHR c1name1[100];
  CHR c1name_vtk[100];
  DBL dcsizc,dcsizs,dcsizv,dcsizf;
  INT i,npoint,ini,start,finish,interval;

  CHRw(stderr,"Start to convert files\n");
  CHRcpy(c1name1,argv[2]);
  fcoeff=fopen(c1name1,"r");
  if(fcoeff!=FILENULL)
  { CHRr(fcoeff,name);
	while(FILEND(fcoeff)==0) 
	{ 
      if (CHRcmp(name,"/YD/YDC/DCSIZC",14)==0)
	  { DBLr(fcoeff,&dcsizc);
	  }
	  if (CHRcmp(name,"/YD/YDC/DCSIZF",14)==0)
	  { DBLr(fcoeff,&dcsizf);
	  }
	  if (CHRcmp(name,"/YD/YDC/DCSIZS",14)==0)
	  { DBLr(fcoeff,&dcsizs);
	  }
	  if (CHRcmp(name,"/YD/YDC/DCSIZV",14)==0)
	  { DBLr(fcoeff,&dcsizv);
	  }
	  CHRr(fcoeff,name);
	}
  }
  else
  { CHRw(stderr,"Could not open input file - usage -i inputfile");
	CHRwcr(stderr);
	return 0;
  }

  start=atoi(argv[3]); finish=atoi(argv[4]); interval=atoi(argv[5]);
  // start=0; finish=10000; interval=1;
  for(i=start;i<finish;i=i+interval)
  { vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
    
    vtkPoints *pts = vtkPoints::New();
    vtkCellArray *polys = vtkCellArray::New();

    vtkDoubleArray *veloc_vector = vtkDoubleArray::New();
    veloc_vector->SetName("Velocity vector");
    veloc_vector->SetNumberOfComponents(3);

    vtkDoubleArray *stress_tensor = vtkDoubleArray::New();
    stress_tensor->SetName("Stress tensor");
    stress_tensor->SetNumberOfComponents(9);

	vtkDoubleArray *sigma1_scalar = vtkDoubleArray::New();
	sigma1_scalar->SetName("Sigma1");
    sigma1_scalar->SetNumberOfComponents(1);

	vtkDoubleArray *sigma3_scalar = vtkDoubleArray::New();
	sigma3_scalar->SetName("Sigma3");
    sigma3_scalar->SetNumberOfComponents(1);

	vtkDoubleArray *sigma_diff_scalar = vtkDoubleArray::New();
	sigma_diff_scalar->SetName("Differential stress");
    sigma_diff_scalar->SetNumberOfComponents(1);

	vtkDoubleArray *sigma_mean_scalar = vtkDoubleArray::New();
	sigma_mean_scalar->SetName("Mean stress");
    sigma_mean_scalar->SetNumberOfComponents(1);

	vtkDoubleArray *sigma1_vector = vtkDoubleArray::New();
	sigma1_vector->SetName("Sigma1 Vectors");
	sigma1_vector->SetNumberOfComponents(3);

	vtkDoubleArray *sigma3_vector = vtkDoubleArray::New();
	sigma3_vector->SetName("Sigma3 vectors");
	sigma3_vector->SetNumberOfComponents(3);

	vtkDoubleArray *cforce_vector = vtkDoubleArray::New();
	cforce_vector->SetName("Contact force vectors");
	cforce_vector->SetNumberOfComponents(3);

    dataSet->Allocate();

    vtkZLibDataCompressor* myZlibCompressor  = vtkZLibDataCompressor::New();
    myZlibCompressor->SetCompressionLevel(9);

    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();

    npoint=0;

    CHRcpy(c1name,argv[1]);
    SINTw(c1tmp,i,0);
    CHRcat(c1name,c1tmp);
    CHRcpy(c1tmp,c1name);
    CHRcpy(c1tmp1,c1name);
    CHRcat(c1name,".ym");
    fdata=fopen(c1name,"r");
    if(fdata!=FILENULL)
    { readppm(fdata,pts,polys,veloc_vector,stress_tensor,sigma1_scalar,sigma3_scalar,sigma1_vector,
		  sigma3_vector,sigma_diff_scalar,sigma_mean_scalar,cforce_vector,npoint,dcsizc,dcsizs,dcsizf,dcsizv);
      fclose(fdata);
    }
    else break;

    dataSet->SetPoints(pts);
    pts->Delete();

    dataSet->GetPointData()->AddArray(veloc_vector);
    dataSet->GetPointData()->SetActiveAttribute("Velocity vector",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(veloc_vector);
    veloc_vector->Delete();

    dataSet->GetPointData()->AddArray(stress_tensor);
    dataSet->GetPointData()->SetActiveAttribute("Stress tensor",vtkDataSetAttributes::TENSORS);
    dataSet->GetPointData()->SetTensors(stress_tensor);
    stress_tensor->Delete();

	dataSet->GetPointData()->AddArray(sigma1_scalar);
    dataSet->GetPointData()->SetActiveAttribute("Sigma1 Scalar",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(sigma1_scalar);
    sigma1_scalar->Delete();

	dataSet->GetPointData()->AddArray(sigma3_scalar);
    dataSet->GetPointData()->SetActiveAttribute("Sigma3 Scalar",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(sigma3_scalar);
    sigma3_scalar->Delete();

	dataSet->GetPointData()->AddArray(sigma1_vector);
    dataSet->GetPointData()->SetActiveAttribute("Sigma1 Vector",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(sigma1_vector);
    sigma1_vector->Delete();

	dataSet->GetPointData()->AddArray(sigma3_vector);
    dataSet->GetPointData()->SetActiveAttribute("Sigma3 Vector",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(sigma3_vector);
    sigma3_vector->Delete();

	dataSet->GetPointData()->AddArray(sigma_diff_scalar);
    dataSet->GetPointData()->SetActiveAttribute("Differential stress",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(sigma_diff_scalar);
    sigma_diff_scalar->Delete();

	dataSet->GetPointData()->AddArray(sigma_mean_scalar);
    dataSet->GetPointData()->SetActiveAttribute("Mean stress",vtkDataSetAttributes::SCALARS);
    dataSet->GetPointData()->SetScalars(sigma_mean_scalar);
    sigma_mean_scalar->Delete();

	dataSet->GetPointData()->AddArray(cforce_vector);
    dataSet->GetPointData()->SetActiveAttribute("Contact force vector",vtkDataSetAttributes::VECTORS);
    dataSet->GetPointData()->SetVectors(cforce_vector);
    cforce_vector->Delete();

    dataSet->SetCells(VTK_TRIANGLE,polys);
    polys->Delete();

    CHRcpy(c1name_vtk,c1tmp1);
    CHRcat(c1name_vtk,".vtu");
    CHRw(stderr, "Writing vtu xml file ");
    CHRw(stderr, c1name_vtk);
    CHRwcr(stderr);
    writer->SetInput(dataSet);
    writer->SetFileName(c1name_vtk);
    writer->SetCompressor(myZlibCompressor);
    writer->Write();
    dataSet->Delete();
    myZlibCompressor->Delete();
    writer->Delete();
  }
}
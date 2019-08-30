#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkZLibDataCompressor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkSystemIncludes.h>

extern "C" {
#include "frame.h"
#include "Ytypes.h"
}

static void ExtractFracturePattern(FILE *fdata, DBL dcsizc, DBL dcsizs, DBL dcsizf, DBL dcsizv,
	vtkPoints* Points, vtkCellArray* Polygons, 
	vtkDoubleArray *scalar_sheardisp, vtkDoubleArray *scalar_failmode, vtkDoubleArray *scalar_fractype,
	vtkDoubleArray *vector_contforce)
{ /* Define variables */
  CHR c1code[1000];
  INT i1num[100];
  DBL d1num[100];
  INT inode,ifrjt,ibdjt;
  INT nbdjt;					// number of boundary joints
  INT nfrjt;					// number of fracture joints
  DBL **d2bjcx,**d2bjcy;		// boundary joint coordinates
  DBL **d2fjcx,**d2fjcy;		// fracture joint coordinates
  INT npoint;					// number of vtk points
  DBL *d1jsdc;					// joint current shear displacement
  DBL *d1jfmd;					// joint failure mode
  DBL *d1jtyp;					// joint type
  DBL **d2bdcf;					// boundary contact force

  /* Count number of joints */
  nbdjt=0; nfrjt=0;
  CHRr(fdata,c1code);
  CHRr(fdata,c1code);
  while(FILEND(fdata)==0)
  { codeCHRtoINT(c1code,i1num);
	codeINTtoDBL(d1num,i1num);
	if(i1num[2]==YTE2BOUNDS)
	{ nbdjt++;
	}
	else if(i1num[2]==YTE2JOINTS)
	{ nfrjt++;
	}
	CHRr(fdata,c1code);
  }
  rewind(fdata);
  npoint=0;

  /* Extract joint element data */
  d2bjcx=TalDBL2(nbdjt,2);
  d2bjcy=TalDBL2(nbdjt,2);
  d2bdcf=TalDBL2(nbdjt,4);
  d2fjcx=TalDBL2(nfrjt,4);
  d2fjcy=TalDBL2(nfrjt,4);
  d1jsdc=TalDBL1(nfrjt);
  d1jfmd=TalDBL1(nfrjt);
  d1jtyp=TalDBL1(nfrjt);
  ifrjt=0; ibdjt=0;
  if(fdata!=FILENULL)
  { CHRr(fdata,c1code);
	CHRr(fdata,c1code);
	while((FILEND(fdata)==0))
	{ codeCHRtoINT(c1code,i1num);
	  codeINTtoDBL(d1num,i1num);
	  if(i1num[2]==YTE2BOUNDS)
	  { // boundary joints
		d2bjcx[ibdjt][0]=d1num[3]*dcsizc;
		d2bjcx[ibdjt][1]=d1num[4]*dcsizc;
		d2bjcy[ibdjt][0]=d1num[7]*dcsizc;
		d2bjcy[ibdjt][1]=d1num[8]*dcsizc;
		d2bdcf[ibdjt][0]=d1num[11]*dcsizf;
		d2bdcf[ibdjt][1]=d1num[12]*dcsizf;
		ibdjt++;
	  }
	  else if(i1num[2]==YTE2JOINTS)
	  { // fracture joints
		d2fjcx[ifrjt][0]=d1num[3]*dcsizc;
		d2fjcx[ifrjt][1]=d1num[4]*dcsizc;
		d2fjcx[ifrjt][2]=d1num[5]*dcsizc;
		d2fjcx[ifrjt][3]=d1num[6]*dcsizc;
		d2fjcy[ifrjt][0]=d1num[7]*dcsizc;
		d2fjcy[ifrjt][1]=d1num[8]*dcsizc;
		d2fjcy[ifrjt][2]=d1num[9]*dcsizc;
		d2fjcy[ifrjt][3]=d1num[10]*dcsizc;
		d1jfmd[ifrjt]=d1num[11]*R2;
		if(d1jfmd[ifrjt]<0.5)
		{ d1jtyp[ifrjt]=0;
		}
		else
		{ d1jtyp[ifrjt]=1;
		}
		d1jsdc[ifrjt]=DABS(d1num[17]*0.05);
		ifrjt++;
	  }
	  CHRr(fdata,c1code);
  } }

  /* Plot fracture pattern */
  for(ifrjt=0;ifrjt<nfrjt;ifrjt++)
  { // vtk output
	Polygons->InsertNextCell(2);
	for(inode=0;inode<2;inode++)
	{ Points->InsertNextPoint(d2fjcx[ifrjt][inode],d2fjcy[ifrjt][inode],R0);
	  Polygons->InsertCellPoint(npoint); npoint++;
	  scalar_sheardisp->InsertNextTuple1(d1jsdc[ifrjt]);
	  scalar_failmode ->InsertNextTuple1(d1jfmd[ifrjt]);
	  scalar_fractype ->InsertNextTuple1(d1jtyp[ifrjt]);
	  vector_contforce->InsertNextTuple3(R0,R0,R0);
	}
	Polygons->InsertNextCell(2);
	for(inode=0;inode<2;inode++)
	{ Points->InsertNextPoint(d2fjcx[ifrjt][inode+2],d2fjcy[ifrjt][inode+2],R0);
	  Polygons->InsertCellPoint(npoint); npoint++;
	  scalar_sheardisp->InsertNextTuple1(d1jsdc[ifrjt]);
	  scalar_failmode ->InsertNextTuple1(d1jfmd[ifrjt]);
	  scalar_fractype ->InsertNextTuple1(d1jtyp[ifrjt]);
	  vector_contforce->InsertNextTuple3(R0,R0,R0);
  } }

  for(ibdjt=0;ibdjt<nbdjt;ibdjt++)
  { // vtk data
	Polygons->InsertNextCell(2);
	for(inode=0;inode<2;inode++)
	{ Points->InsertNextPoint(d2bjcx[ibdjt][inode],d2bjcy[ibdjt][inode],R0);
	  Polygons->InsertCellPoint(npoint); npoint++;
	  scalar_sheardisp->InsertNextTuple1(R0);
	  scalar_failmode ->InsertNextTuple1(-R1);
	  scalar_fractype ->InsertNextTuple1(-R1);
	  vector_contforce->InsertNextTuple3(d2bdcf[ibdjt][0],d2bdcf[ibdjt][1],R0);
  } }

  FREE(d2bjcx);
  FREE(d2bjcy);
  FREE(d2fjcx);
  FREE(d2fjcy);
  FREE(d1jsdc);
  FREE(d1jfmd);
  FREE(d1jtyp);
  FREE(d2bdcf);
}

int main(int argc, char **argv)
{ /* Define variables */
  FILE *fdata=FILENULL;
  FILE *fcoeff=FILENULL;
  CHR yfname[300];
  CHR prjname[100];
  CHR ymfname[100];
  CHR vtkname[100];
  CHR fileID[100];
  DBL dcsizc,dcsizs,dcsizv,dcsizf;
  INT i,start,finish,interval;

  /* Read the .y file */
  CHRcpy(prjname,argv[2]);
  fcoeff=fopen(prjname,"r");
  if(fcoeff!=FILENULL)
  { CHRr(fcoeff,yfname);
	while(FILEND(fcoeff)==0)
	{ if(CHRcmp(yfname,"/YD/YDC/DCSIZC",14)==0)
	  { DBLr(fcoeff,&dcsizc);
	  }
	  if(CHRcmp(yfname,"/YD/YDC/DCSIZF",14)==0)
	  { DBLr(fcoeff,&dcsizf);
	  }
	  if(CHRcmp(yfname,"/YD/YDC/DCSIZS",14)==0)
	  { DBLr(fcoeff,&dcsizs);
	  }
	  if(CHRcmp(yfname,"/YD/YDC/DCSIZV",14)==0)
	  { DBLr(fcoeff,&dcsizv);
	  }
	  CHRr(fcoeff,yfname);
	}
  }
  else
  { CHRw(stderr,"Could not open input file - usage -i inputfile"); 
	CHRwcr(stderr);   
	return 0;
  }

  /* Process the current .ym file */
  start=atoi(argv[3]);
  finish=atoi(argv[4]);
  interval=atoi(argv[5]);
  for(i=start;i<finish;i=i+interval)
  { // fracture data set
	vtkUnstructuredGrid *dataset = vtkUnstructuredGrid::New();
	dataset->Allocate();

	// fracture point object
	vtkPoints *fracpts = vtkPoints::New();

	// fracture polygon objects
	vtkCellArray *fracpolys = vtkCellArray::New();

	// fracture displacement
	vtkDoubleArray *scalar_sheardisp = vtkDoubleArray::New();
	scalar_sheardisp->SetName("Shear displacement");
	scalar_sheardisp->SetNumberOfComponents(1);

	// fracture failure mode
	vtkDoubleArray *scalar_failmode = vtkDoubleArray::New();
	scalar_failmode->SetName("Failure mode");
	scalar_failmode->SetNumberOfComponents(1);

	// fracture type
	vtkDoubleArray *scalar_fractype = vtkDoubleArray::New();
	scalar_fractype->SetName("Fracture type");
	scalar_fractype->SetNumberOfComponents(1);

	// fracture displacement
	vtkDoubleArray *vector_contforce = vtkDoubleArray::New();
	vector_contforce->SetName("Contact force");
	vector_contforce->SetNumberOfComponents(3);

	// vtk compressor
	vtkZLibDataCompressor* zlibcompressor  = vtkZLibDataCompressor::New();
	zlibcompressor->SetCompressionLevel(9);

	// vtk writer
	vtkXMLUnstructuredGridWriter* datawriter = vtkXMLUnstructuredGridWriter::New();

	// extract fracture pattern
	CHRcpy(ymfname,argv[1]);
	SINTw(fileID,i,0);
	CHRcat(ymfname,fileID);
	CHRcat(ymfname,".ym");
	fdata=fopen(ymfname,"r");
	if(fdata!=FILENULL)
	{ ExtractFracturePattern(fdata,dcsizc,dcsizs,dcsizf,dcsizv,fracpts,fracpolys,scalar_sheardisp,scalar_failmode,scalar_fractype,vector_contforce);
	  fclose(fdata);
	}
	else break;

	// plot fracture pattern
	dataset->SetPoints(fracpts);
	fracpts->Delete();
	dataset->SetCells(VTK_LINE,fracpolys);
	fracpolys->Delete();

	dataset->GetPointData()->AddArray(scalar_sheardisp);
	dataset->GetPointData()->SetActiveAttribute("Shear displacement",vtkDataSetAttributes::SCALARS);
	dataset->GetPointData()->SetScalars(scalar_sheardisp);
	scalar_sheardisp->Delete();

	dataset->GetPointData()->AddArray(scalar_failmode);
	dataset->GetPointData()->SetActiveAttribute("Failure mode",vtkDataSetAttributes::SCALARS);
	dataset->GetPointData()->SetScalars(scalar_failmode);
	scalar_failmode->Delete();

	dataset->GetPointData()->AddArray(scalar_fractype);
	dataset->GetPointData()->SetActiveAttribute("Fracture type",vtkDataSetAttributes::SCALARS);
	dataset->GetPointData()->SetScalars(scalar_fractype);
	scalar_fractype->Delete();

	dataset->GetPointData()->AddArray(vector_contforce);
    dataset->GetPointData()->SetActiveAttribute("Contact force",vtkDataSetAttributes::VECTORS);
    dataset->GetPointData()->SetVectors(vector_contforce);
    vector_contforce->Delete();

	// output the vtu file
	CHRcpy(vtkname,argv[1]);
	CHRcat(vtkname,"_crack");
	SINTw(fileID,i,0);
	CHRcat(vtkname,fileID);
	CHRcat(vtkname,".vtu");
	CHRw(stderr,"Writing vtk xml file ");
	CHRw(stderr, vtkname);
	CHRwcr(stderr);
	datawriter->SetInput(dataset);
	datawriter->SetFileName(vtkname);
	datawriter->Write();
	dataset->Delete();
	datawriter->Delete();
  }
}
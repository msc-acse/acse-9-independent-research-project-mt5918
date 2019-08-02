*realformat "%20.10e"
*intformat "%7i"
*Set Cond Surface-Constraints *nodes *NoRepeat
*Add Cond Point-Constraints *nodes *NoRepeat
*Set Var NConstraints0=CondNumEntities
*set Cond Line-Constraints *nodes *NoRepeat
*Set Var NConstraints1=CondNumEntities
*Set Var NConstraints=NConstraints0+NConstraints1
/YD/YDC/MCSTEP *GenData(1,int)
/YD/YDC/NCSTEP 0
/YD/YDC/DCGRAY *GenData(2)
/YD/YDC/DCSTEC *GenData(3)
/YD/YDC/ICOUTF *GenData(4,int)
/YD/YDC/ICSAVF *GenData(5,int)
/YD/YDC/DCGRST *GenData(6)
/YD/YDC/DCRMPT *GenData(7)
/YD/YDC/DCSIZC *GenData(8)
/YD/YDC/DCSIZF *GenData(9)
/YD/YDC/DCSIZV *GenData(10)
/YD/YDC/DCSIZS *GenData(11)
/YD/YDC/DCSIZD *GenData(12)
/YD/YDC/DCSIZA *GenData(13)
/YD/YDC/DCTIME 0.0
/YD/YDC/ICOUTI 0
*if((strcmp(GenData(16),"64bit")==0))
/YD/YDC/ICOUTP 6
*endif
*if((strcmp(GenData(16),"32bit")==0))
/YD/YDC/ICOUTP 4
*endif
*if((strcmp(GenData(17),"Coulomb")==0))
/YD/YDC/ICFMTY 0
*endif
*if((strcmp(GenData(17),"Barton-Bandis")==0))
/YD/YDC/ICFMTY 1
*endif
*if((strcmp(GenData(18),"roughness")==0))
/YD/YDC/ICIATY 0
*endif
*if((strcmp(GenData(18),"length")==0))
/YD/YDC/ICIATY 1
*endif

*Set elems(triangle)
*Set Var melem=nelem+1
/YD/YDE/MELEM *melem
/YD/YDE/NELEM *nelem
/YD/YDE/MELST 2  
/YD/YDE/NELST 2 
/YD/YDE/MELNO 4  
/YD/YDE/NELNO 3 
/YD/YDE/D2ELST 21 2 0
/YD/YDE/I1ELCF 0
/YD/YDE/I1ELTY 0
/YD/YDE/I1ELPR *nelem 
*loop elems
*ElemsMat 
*end elems

/YD/YDE/I2ELTO 21 3 *nelem 
*loop elems
*ElemsConec
*end elems

/YD/YDI/MICOUP *GenData(14)
/YD/YDI/NICOUP 0 
/YD/YDI/IIECFF -2
/YD/YDI/DIEDI  200.0
/YD/YDI/DIEZON *GenData(15)
/YD/YDI/D1IESL 0
/YD/YDI/I1IECN 0
/YD/YDI/I1IECT 0

*Set Var mnopo=npoin+1
/YD/YDN/MNODIM 3 
/YD/YDN/NNODIM 2
/YD/YDN/MNOPO *mnopo 
/YD/YDN/NNOPO *npoin

/YD/YDN/D2NCC 21 2 *npoin  
*loop nodes
*NodesCoord
*end nodes 

/YD/YDN/D2NCI 21 2 *npoin  
*loop nodes
*NodesCoord 
*end nodes

/YD/YDN/D2NFC 21  0 0 
/YD/YDN/D1NMCT 0
/YD/YDN/D2NVC 21 0 0
/YD/YDN/I1NOBF *npoin 0
/YD/YDN/I1NOPR *npoin 0 
/YD/YDP/MPROP      200
/YD/YDP/NPROP      200
/YD/YDP/D1PEFR     200 0
/YD/YDP/D1PEFT     200 0
/YD/YDP/D1PEGT     200 0
/YD/YDP/D1PEGS     200 0
/YD/YDP/D1PEKS     200 0
/YD/YDP/D1PELA     200 0
/YD/YDP/D1PEMU     200 0
/YD/YDP/D1PEPE     200 0
/YD/YDP/D1PEPC     200 0
/YD/YDP/D1PEPF     200 0
/YD/YDP/D1PBIF     200 0
/YD/YDP/D1PERO     200 0
/YD/YDP/D1PNAP     200 0
/YD/YDP/D1PNAF     200 0
/YD/YDP/D1PNAI     200 0
/YD/YDP/D1PNAT     200 0
/YD/YDP/D1PNAX     200 0
/YD/YDP/D1PNAY     200 0
/YD/YDP/D1PNXX     200 1
/YD/YDP/D1PNXY     200 0
/YD/YDP/D1PNYX     200 0
/YD/YDP/D1PNYY     200 1
/YD/YDP/D1PSEM     200 0
/YD/YDP/D1PICF     200 0
/YD/YDP/D1PCOH     200 0
/YD/YDP/D1PJRC     200 0
/YD/YDP/D1PJCS     200 0
/YD/YDP/D1PJSL     200 0
/YD/YDP/I1PEFR     200 0
/YD/YDP/I1PEJP     200 0
/YD/YDP/I1PEMB     200 0
/YD/YDP/I1PEMN     200 0
/YD/YDP/I1PNFX     200 0
/YD/YDP/I1PNFY     200 0
/YD/YDP/I1PSDE     200 0
/YD/YDP/I1PTYP     200 0

/YD/YDO/MOHYS 1 /YD/YDO/NOHYS 1
/YD/YDO/DOHYP 0.0005
/YD/YDO/D1OHYC 1
1.0
/YD/YDO/D1OHYF 1
1.0 
/YD/YDO/D1OHYS 1
0
/YD/YDO/D1OHYT 1
0
/YD/YDO/D1OHYX 1
0
/YD/YDO/D1OHYY 1
0
/YD/YDO/I1OHYT 1
15

Gen_out_only *GenData(19)

MATERIAL_LIST *nmats 
*loop materials
*MatProp(1,real) *MatProp(2,real) *MatProp(3,real) *MatProp(4,real) *MatProp(5,real) *MatProp(6,real) *MatProp(7,real) *MatProp(8,real) *MatProp(9,real) *MatProp(10,real) *MatProp(11,real) *MatProp(12,real) *MatProp(13,real) *MatProp(14,real) *MatProp(15,real) *MatProp(16,real) *MatProp(17,real) *MatProp(18)
*end materials

*set Cond Element  *elems *CanRepeat
*Set Var NElement=CondNumEntities
Element_List *NElement
*loop elems *OnlyInCond
*ElemsNum *cond(1) *cond(2,int)
*end nodes

*set Cond Surface-Initial-Data  *nodes *NoRepeat
*Add Cond Line-Initial-Data  *nodes *NoRepeat
*Add Cond Point-Initial-Data *nodes *NoRepeat
*Set Var NInitial=CondNumEntities
Initial_Data *NInitial
*loop nodes *OnlyInCond
*NodesNum *cond(1,int) *cond(2,real) *cond(3,real) *cond(4,int) *cond(5,real) *cond(6,real)
*end nodes
*Set Cond Surface-Constraints *nodes *NoRepeat
*Add Cond Line-Constraints *nodes *NoRepeat
*Add Cond Point-Constraints *nodes *NoRepeat
*Set Var NConstraints=CondNumEntities

Constraints_List *NConstraints
*loop nodes *OnlyInCond
*NodesNum *cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6) *cond(7) *cond(8) 
*end nodes

*Set elems(linear)
/YD/YDX/NELEM *nelem 
/YD/YDX/STRIELE 21 2 *nelem 
*loop elems
*ElemsConec
*end elems

$YDOIT
$YSTOP
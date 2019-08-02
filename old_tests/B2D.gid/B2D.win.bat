@ECHO OFF
rem OutputFile: %2\%1.log
del %2\%1.log
cd ..
del 1.tmp
echo "Generating output file ......" >>%2\%1.log
%3\GID_B2D.exe %2\%1.DAT %3\B2D.PAR
if exist 1.tmp goto end
echo "Running job......" >>%2\%1.log
"C:\Program Files (x86)\IC-QMUL\VGeST\Y.exe" %1.y
echo "Converting stress info into VTU files......" >>%2\%1.log
"C:\Program Files (x86)\IC-QMUL\VGeST\m2vtu\m2vtu2D.exe" %1 %1.y
echo "Converting crack info into VTU files......" >>%2\%1.log
"C:\Program Files (x86)\IC-QMUL\VGeST\m2vtu\m2vtu2D_crack.exe" %1 %1.y
echo "Visulizing the results......" >>%2\%1.log
"C:\Program Files (x86)\MayaVi\mayavi" -d %10.vtu -m SurfaceMap -f ExtractTensorComponents -d %1_crack0.vtu -m SurfaceMap
:end
del 1.tmp
del %2\%1.log


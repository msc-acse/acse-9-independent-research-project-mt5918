@ECHO ON

title Generate Y file and run Y2d
rem
rem ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rem Preparations
rem ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rem
for /F "tokens=1 delims=. " %%a in ("%~n0") do (
   set filename=%%a
)
mkdir ..\gid\%filename%
mkdir ..\gid\%filename%\results
cd ..\gid\%filename%\results
rem
echo %filename%
set logfile=..\%filename%.log
echo %logfile%
rem
del %logfile%
del 1.tmp
rem
echo --------------------------------------------------- >> %logfile%
echo Log information >> %logfile%
echo --------------------------------------------------- >> %logfile%
echo Log for >> %logfile%
echo     file: %0 >> %logfile%
echo     1st arg: %1 >> %logfile%
echo     2nd arg: %2 >> %logfile%
rem
echo --------------------------------------------------- >> %logfile%
echo Read userfile >> %logfile%
echo --------------------------------------------------- >> %logfile%
rem
rem ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rem ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rem ~~~~~~~~~~          Set user file          ~~~~~~~~~~~
rem ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rem ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set userfile=../../../user/config.txt >> %logfile%
rem ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rem ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo filename: %filename% >> %logfile%
echo logfile: %logfile% >> %logfile%
set /a n=0
echo n=%n% >> %logfile%
rem
setlocal EnableDelayedExpansion
rem
set N=0
for /f "usebackq tokens=*" %%A in ("%userfile%") do (
  set /a N += 1
  set "Line!N!=%%A"
)
rem
set generate=%Line2%       >> %logfile%
set testfilename=%Line4%   >> %logfile%
set login=%Line6%          >> %logfile%
set hpcdir=%Line8%         >> %logfile%
set localdir=%Line10%      >> %logfile%
set qsubdir=%Line12%       >> %logfile%
set mayavi=%Line14%        >> %logfile%
rem
echo --------------------------------------------------- >> %logfile%
echo Run program >> %logfile%
echo --------------------------------------------------- >> %logfile%
rem
"..\..\B2D.gid\GID_B2D.exe" ..\%filename%.dat ..\%filename%.par >> %logfile%
if exist 1.tmp goto end
"..\..\..\bin\Y.exe" %testfile%
echo "Stress info into VTU files" >> %logfile%
"..\..\..\bin\m2vtu2D.exe" %filename% %testfile% >> %logfile%
echo "Crack info into VTU files" >> %logfile%
"..\..\..\bin\m2vtu2D_crack.exe" %filename% %testfile% >> %logfile%
echo "Visulizing the results" >> %logfile%
"%mayavidir%\mayavi2" -d %filename%0.vtu -m SurfaceMap -f ExtractTensorComponents -d %filename%_crack0.vtu -m SurfaceMap >> %logfile%
echo "---------------------------------------------------" >> %logfile%
echo "Done" >> %logfile%
echo "---------------------------------------------------" >> %logfile%
:end
pause

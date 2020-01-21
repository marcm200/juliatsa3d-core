@echo off
if not exist cube-viewer.exe goto error:
if not exist juliatsa3dcore_d.exe goto error2:
del _in.raw
rem compute from scratch and save as bitcube
juliatsa3Dcore_d func=makin cmd=calc len=8 c=-1,0,0 range=2 twdb=0 revcg=4
rem showing interior cells as pure yellow
cube-viewer file=bc1.ccb cmd=show observer=3 transp=250,250,250,255,255,255
rem using computed data
copy _3d.raw _in.raw /Y
rem compute one refinement level higher and color periodic cycles
juliatsa3Dcore_d func=makin cmd=period len=9 c=-1,0,0 range=2 twdb=0 revcg=4
rem show immediate basins and attraction basins
cube-viewer file=bc1per.ccb cmd=show observer=3 transp=250,250,250,255,255,255
echo.
echo.
echo Note the different vieweing position (termed 3) as compared to the standard internal viewer (termed 1).
echo.
goto end:

:error
echo.
echo Please copy cube-viewer.exe into the current directory and start again.
echo.
echo.
goto end:

:error2
echo.
echo Please copy the desired binary (win32, win64) with (opt) or without (noopt) compiler optimizations
echo into the current directory and rename it to juliatsacore3d_d.exe
echo.
echo.
goto end:

:end
pause

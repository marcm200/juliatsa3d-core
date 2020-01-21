@echo off
if not exist juliatsa3dcore_d.exe goto error:
del _in.raw
juliatsa3Dcore_d func=makin cmd=period len=8 c=-1,0,0 range=2 revcg=4
copy _3d.raw _in.raw /Y
juliatsa3Dcore_d func=makin cmd=period len=9 c=-1,0,0 range=2 revcg=4
copy _3d.raw _in.raw /Y
juliatsa3Dcore_d func=makin cmd=period len=10 c=-1,0,0 range=2 revcg=4
goto end:

:error
echo.
echo Please copy the desired binary (win32, win64) with (opt) or without (noopt) compiler optimizations
echo into the current directory and rename it to juliatsacore3d_d.exe
echo.
echo.

goto end:

:end
pause

echo off
for %%A in (*.png) do (
rem echo %%~nA.pdf %%~nA.pdf
convert %%~nA.png eps2:%%~nA.eps
)
echo on

rem "c:\Program Files (x86)\MiKTeX 2.8\miktex\bin\pdf2ps.exe" %%~nA.pdf %%~nA.eps

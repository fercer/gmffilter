@echo off
IF NOT EXIST build (
echo Creating build folder ...
mkdir .\build
)

IF NOT EXIST bin (
echo Creating bin folder ...
mkdir .\bin
)

SET libname=libgmf.lib
SET executablename=test_gmf.exe
SET required_include_paths=/ID:\Apps\fftw\include /I..\include
SET required_libs=D:\Apps\fftw\lib\fftw3.lib ..\lib\%libname%

cl %required_include_paths% /Fo:build\test_gmf.o src\test_gmf.c %required_libs% /link /out:.\bin\%executablename%

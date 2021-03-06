@echo off
IF NOT EXIST build (
echo Creating build folder ...
mkdir .\build
)

IF NOT EXIST lib (
echo Creating lib folder ...
mkdir .\lib
)

SET libname=libgmf.lib
SET dllname=libgmf.dll
SET required_include_paths=/I%FFTW_INCLUDE_PATH%
SET required_libs=%FFTW_LIBS%
SET macros_definitions=/DBUILDING_GMF_DLL
SET version=release
SET install=false
SET installation_path=


if "%1" == "python" (
    SET required_include_paths=%required_include_paths% /I%PYTHON_3_6_INCLUDE_PATH% /I%NUMPY_CKN_INCLUDE_PATH%
    SET required_libs=%required_libs% %PYTHON_3_6_LIBS% %NUMPY_CKN_LIBS%
    SET macros_definitions=%macros_definitions% /DBUILDING_PYTHON_MODULE
    SET dllname=gmf.pyd
    SET installation_path=%PYTHONPATH%
    if "%2" == "debug" SET version="debug"
    if "%3" == "install" SET install=true
) ELSE (
    if "%1" == "debug" SET version="debug"
    if "%2" == "install" (
        SET install=true
        SET installation_path=%3
    )
)
if "%version%" == "release" SET macros_definitions=%macros_definitions% /DNDEBUG /O2

echo "Required include paths:"
echo %required_include_paths%
cl /LD %macros_definitions% %required_include_paths% /Fo:build\gmf.o include\gmf.c %required_libs% /link /DLL /IMPLIB:.\lib\%libname% /out:.\lib\%dllname%

echo %macros_definitions%
echo %PYTHON_3_6_INCLUDE_PATH%
echo %version%
echo %install%
echo %installation_path%
echo Python path: %PYTHONPATH%

if "%install%" == "true" (
    echo Copy module to: %installation_path%\
    copy .\lib\%dllname% %installation_path%\
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::
::::   MIT License
::::
::::   ParaMonte: plain powerful parallel Monte Carlo library.
::::
::::   Copyright (C) 2012-present, The Computational Data Science Lab
::::
::::   This file is part of the ParaMonte library.
::::
::::   Permission is hereby granted, free of charge, to any person obtaining a 
::::   copy of this software and associated documentation files (the "Software"), 
::::   to deal in the Software without restriction, including without limitation 
::::   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
::::   and/or sell copies of the Software, and to permit persons to whom the 
::::   Software is furnished to do so, subject to the following conditions:
::::
::::   The above copyright notice and this permission notice shall be 
::::   included in all copies or substantial portions of the Software.
::::
::::   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
::::   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
::::   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
::::   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
::::   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
::::   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
::::   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
::::
::::   ACKNOWLEDGMENT
::::
::::   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
::::   As per the ParaMonte library license agreement terms, if you use any parts of 
::::   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
::::   work (education/research/industry/development/...) by citing the ParaMonte 
::::   library as described on this page:
::::
::::       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: Example usage: ..\auxil\bzip.bat --dir .\

@echo off
set ERRORLEVEL=0
cd %~dp0

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: parse arguments
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo.
echo.-- ParaMonte - parsing input arguments...
echo.

:LABEL_parseArgLoop

set FLAG_SUPPORTED=true
set VALUE_SUPPORTED=true
set "DESTINATION_DIR="

if not "%1"=="" (

    echo.-- ParaMonte - processing: %1 %2

    set FLAG=%1
    set VALUE=%2

    set FLAG_SUPPORTED=false
    set VALUE_SUPPORTED=false

    REM --dir

    if "!FLAG!"=="--dir" (
        set FLAG_SUPPORTED=true
        set DESTINATION_DIR=%2
        shift
        if exist "!DESTINATION_DIR!" (
            set VALUE_SUPPORTED=true
            echo. 
            echo. -- ParaMonte - input destination directory exists: !DESTINATION_DIR!
            echo. 
        ) else (
            echo. 
            echo. -- ParaMonte - Fatal Error: input destination directory does not exist: !DESTINATION_DIR!
            echo. 
            cd %~dp0
            set ERRORLEVEL=1
            exit /B 1
        )
    )

    shift
    goto :LABEL_parseArgLoop

)

:LABEL_REPORT_ERR

REM check flag/value support

if "!FLAG_SUPPORTED!"=="true" (
    if "!VALUE_SUPPORTED!" NEQ "true" (
        echo.
        echo.-- ParaMonte - FATAL: The requested input value "!VALUE!" specified 
        echo.-- ParaMonte - FATAL: with the input flag "!FLAG!" is not supported.
        goto LABEL_ERR
    )
) else (
    echo.
    echo.-- ParaMonte - FATAL: The requested input flag "!FLAG!" is not supported.
    goto LABEL_ERR
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: zip subfolders
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cd %~dp0
set "AUXIL_DIR=%~dp0"

if not defined DESTINATION_DIR (
    set DESTINATION_DIR=..\bin
)

call :NORMALIZEPATH "!DESTINATION_DIR!"
if exist "!DESTINATION_DIR!" (
    cd "!DESTINATION_DIR!"
    echo. 
    echo. -- ParaMonte - compressing all subdirectories in the directory: !DESTINATION_DIR!
    echo. 
    for /f "tokens=* usebackq" %%G in (`dir /b /a:d "!DESTINATION_DIR!"`) do (
        if exist "%%~G.zip" (
            echo. -- ParaMonte - WARNING: compressed subdirectory already exists: %%~G.zip
            echo. -- ParaMonte - WARNING: skipping...
        ) else (
            echo. -- ParaMonte - compressing subdirectory: %%~G
            !AUXIL_DIR!\7z.exe a -r -tzip "%%~G.zip" "%%~G"
        )
    )
) else (
    echo. 
    echo. -- ParaMonte - Fatal Error: input destination directory does not exist: !DESTINATION_DIR!
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.
exit /B 0

goto LABEL_EOF

:: subroutines

:getLowerCase
:: Subroutine to convert a variable VALUE to all lower case.
:: The argument for this subroutine is the variable NAME.
FOR %%i IN ("/=/" "+=+" "A=a" "B=b" "C=c" "D=d" "E=e" "F=f" "G=g" "H=h" "I=i" "J=j" "K=k" "L=l" "M=m" "N=n" "O=o" "P=p" "Q=q" "R=r" "S=s" "T=t" "U=u" "V=v" "W=w" "X=x" "Y=y" "Z=z") DO CALL SET "%1=%%%1:%%~i%%"
GOTO:EOF

:getUpperCase
:: Subroutine to convert a variable VALUE to all UPPER CASE.
:: The argument for this subroutine is the variable NAME.
FOR %%i IN ("/=/" "+=+" "a=A" "b=B" "c=C" "d=D" "e=E" "f=F" "g=G" "h=H" "i=I" "j=J" "k=K" "l=L" "m=M" "n=N" "o=O" "p=P" "q=Q" "r=R" "s=S" "t=T" "u=U" "v=V" "w=W" "x=X" "y=Y" "z=Z") DO CALL SET "%1=%%%1:%%~i%%"
GOTO:EOF

:getTitleCase
:: Subroutine to convert a variable VALUE to Title Case.
:: The argument for this subroutine is the variable NAME.
FOR %%i IN ("/=/" "+=+" " a= A" " b= B" " c= C" " d= D" " e= E" " f= F" " g= G" " h= H" " i= I" " j= J" " k= K" " l= L" " m= M" " n= N" " o= O" " p= P" " q= Q" " r= R" " s= S" " t= T" " u= U" " v= V" " w= W" " x= X" " y= Y" " z= Z") DO CALL SET "%1=%%%1:%%~i%%"
GOTO:EOF

:LABEL_ERR

cd %~dp0
set ERRORLEVEL=1
exit /B 1

:NORMALIZEPATH
set DESTINATION_DIR=%~dpfn1
exit /B

:LABEL_EOF

echo.
echo.-- ParaMonte - mission accomplished. 
echo.

exit /B 0

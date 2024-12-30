::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::                                                                                                                            ::::
::::    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ::::
::::                                                                                                                            ::::
::::    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ::::
::::                                                                                                                            ::::
::::    This file is part of the ParaMonte library.                                                                             ::::
::::                                                                                                                            ::::
::::    LICENSE                                                                                                                 ::::
::::                                                                                                                            ::::
::::       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ::::
::::                                                                                                                            ::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: Example usage: ..\auxil\bzip.bat --dir .\

@echo off
set ERRORLEVEL=0
cd %~dp0
set "AUXIL_DIR=%~dp0"
set "paramonte_dir=!AUXIL_DIR!.."

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: parse arguments
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo.
echo. -- ParaMonte - parsing input arguments...
echo.

:LABEL_parseArgLoop

set "DESTINATION_DIR="

if not "%1"=="" (

    echo. -- ParaMonte - processing: %1 %2

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
            cd !paramonte_dir!
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
        echo. -- ParaMonte - FATAL: The requested input value "!VALUE!" specified 
        echo. -- ParaMonte - FATAL: with the input flag "!FLAG!" is not supported.
        goto LABEL_ERR
    )
) else (
    echo.
    echo. -- ParaMonte - FATAL: The requested input flag "!FLAG!" is not supported.
    goto LABEL_ERR
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: zip subfolders
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if not defined DESTINATION_DIR set DESTINATION_DIR=!paramonte_dir!\_bin

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
            !AUXIL_DIR!\7z.exe a -r -tzip "%%~G.zip" "%%~G" || (
                echo.
                echo. -- ParaMonte - FATAL: compression failed for subdirectory: %%~G
                echo. -- ParaMonte - FATAL: gracefully exiting.
                echo.
                cd !paramonte_dir!
                set ERRORLEVEL=1
                exit /B 1
            )
        )
    )
) else (
    echo. 
    echo. -- ParaMonte - Fatal Error: input destination directory does not exist: !DESTINATION_DIR!
    echo. 
    cd !paramonte_dir!
    set ERRORLEVEL=1
    exit /B 1
)
echo.


set "FLAG="
set "VALUE="
set "AUXIL_DIR="
set "DESTINATION_DIR="
cd !paramonte_dir!
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

cd !paramonte_dir!
set ERRORLEVEL=1
exit /B 1

:NORMALIZEPATH
cd !paramonte_dir!
set DESTINATION_DIR=%~dpfn1
exit /B

:LABEL_EOF

echo.
echo. -- ParaMonte - mission accomplished. 
echo.

cd !paramonte_dir!
exit /B 0

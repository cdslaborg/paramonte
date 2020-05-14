####################################################################################################################################
####################################################################################################################################
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   ParaMonte is free software: you can redistribute it and/or modify it
####   under the terms of the GNU Lesser General Public License as published
####   by the Free Software Foundation, version 3 of the License.
####
####   ParaMonte is distributed in the hope that it will be useful,
####   but WITHOUT ANY WARRANTY; without even the implied warranty of
####   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
####   GNU Lesser General Public License for more details.
####
####   You should have received a copy of the GNU Lesser General Public License
####   along with the ParaMonte library. If not, see,
####
####       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
####
####   ACKNOWLEDGMENT
####
####   As per the ParaMonte library license agreement terms,
####   if you use any parts of this library for any purposes,
####   we ask you to acknowledge the ParaMonte library's usage
####   in your work (education/research/industry/development/...)
####   by citing the ParaMonte library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

import os as _os
import sys as _sys
import _paramonte as _pm
import warnings as _warnings

verificationStatusFilePath = _os.path.join( _pm.path.auxil, ".verificationEnabled" )

platform = _sys.platform.lower()
isWin32 = True if platform=="win32" else False
isLinux = True if platform=="linux" else False
isMacOS = True if platform=="darwin" else False

####################################################################################################################################

buildInstructionNote    = "If your platform is non-Windows and is compatible with GNU Compiler Collection (GCC),\n" \
                        + "you can also build the required ParaMonte kernel's shared object files on your system\n" \
                        + "by calling ParaMonte module's build() function from within your Python enviroment, like:\n\n" \
                        + "    import paramonte as pm\n" \
                        + "    pm.build()"

####################################################################################################################################

def verify(reset = True):
    """
    checks (or rechecks) the requirements of the installed ParaMonte library

    Parameters
    ----------
        reset
            boolean whose default value is True. If True, 
            a thorough verification of the existence of the required 
            libraries will performed, as if it is the first ParaMonte 
            module import.

    Returns
    -------
        None

    """

    if not isinstance(reset,bool):
        raise Exception( "\n\nThe input argument reset must be a logical (boolean) value: True or False\n\n")

    # require Python >3.6 for type hints

    _MIN_PYTHON = (3,6)
    if _sys.version_info < _MIN_PYTHON:
        _sys.exit("Python %s.%s or newer is required for ParaMonte. Please install the latest version of Python.\n" % _MIN_PYTHON)

    # ensure messages are printed only for the first-time import

    if reset: writeVerificationStatusFile(verificationEnabled=True)

    with open(verificationStatusFilePath, "r") as verificationStatusFile:
        verificationEnabledStr = verificationStatusFile.readline()

    if verificationEnabledStr=="False":
        verificationEnabled = False
    elif verificationEnabledStr=="True":
        verificationEnabled = True
    else:
        raise Exception ( "\n\nThe internal settings of the ParaMonte library appears to have been messed up\n"
                        + "potentially by the operating system, Python, or some third party applications.\n"
                        + "Please reinstall ParaMonte by typing the following commands \n"
                        + "on a Python-aware command-line interface:\n\n"
                        + "    pip uninstall paramonte\n"
                        + "    pip install paramonte\n\n"
                        )

    if verificationEnabled:

        # ensure 64-bit architecture

        if (_pm.arch=="x64") and (isWin32 or isLinux or isMacOS):

            displayParaMonteBanner()

            # library path

            if not isWin32: setupUnixPath()

            # search for the MPI library

            mpiBinPath = findMPI()
            if mpiBinPath is None:
                _pm.warn( msg   = "The MPI runtime libraries for 64-bit architecture could not be detected on your system. \n"
                                + "The MPI runtime libraries are required for the parallel ParaMonte simulations. \n"
                                + "For Windows and Linux operating systems, you can download and install the Intel MPI runtime \n"
                                + "libraries, free of charge, from Intel website (https://software.intel.com/en-us/mpi-library). \n"
                                + "For macOS (Darwin operating system), you will have to download and install the Open-MPI library \n"
                                + "(https://www.open-mpi.org/).\n\n"
                                + "Alternatively, the ParaMonte library can automatically install these library for you now. \n\n"
                                + "If you don't know how to download and install the correct MPI runtime library version, \n"
                                + "we strongly recommend that you let the ParaMonte library to install this library for you. \n"
                                + "If so, ParaMonte will need access to the World-Wide-Web to download the library \n"
                                + "and will need your administrative permission to install it on your system. \n"
                                + "Therefore, if you have any active firewall on your system such as ZoneAlarm, \n"
                                + "please make sure your firewall allows ParaMonte to access the Internet."
                        , marginTop = 1
                        , marginBot = 1
                        , methodName = _pm.names.paramonte
                        )

                isYes = getUserResponse( msg =  "\n    Do you wish to download and install the MPI runtime library"
                                                "\n    for parallel simulations on your system now (y/n)? " 
                                                )
                if isYes:
                    installMPI()
                    writeVerificationStatusFile(verificationEnabled=True)
                else:
                    _pm.note( msg   = "Skipping the MPI library installation... \n"
                                    + "It is now the user's responsibility to install the \n"
                                    + "required libraries for parallel simulations. \n"
                                    + "If you ever wish to install MPI libraries via ParaMonte again, \n"
                                    + "please try:\n\n"
                                    + "    import paramonte as pm\n"
                                    + "    pm.verify()\n\n"
                                    + "Alternatively, you can also uninstall and reinstall the package \n"
                                    + "via the following commands on a Python-aware command-line interface:\n\n"
                                    + "    pip uninstall paramonte\n"
                                    + "    pip install --user paramonte\n\n"
                                    + "For more information visit:\n\n"
                                    + "    https://www.cdslab.org/pm"
                            , methodName = _pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )
                    writeVerificationStatusFile(verificationEnabled=False)

            else:

                writeVerificationStatusFile(verificationEnabled=False)

            dispFinalMessage()

        else:

            warnForUnsupportedPlatform()
            build()

    return None

####################################################################################################################################

def getUserResponse(msg=""):
    while True:
        answer = input(msg)
        if answer.lower()=="y":
            return True
        elif answer.lower()=="n":
            return False
        else:
            print("Invalid answer. Please enter either y or n, then press enter.")
            continue

####################################################################################################################################

def download(url,filePath):
    import urllib.request
    import shutil
    with urllib.request.urlopen(url) as response, open(filePath, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
    return None

####################################################################################################################################

def warnForUnsupportedPlatform():
    _pm.warn( msg   = "The ParaMonte sampler kernel is currently exclusively available \n"
                    + "on AMD64 (64-bit) architecture for Windows/Linux/Darwin Operating Systems (OS). \n"
                    + "Your system appears to be of a different architecture or OS. As a result, \n"
                    + "the core sampler routines of ParaMonte will not be available on your system. However, \n"
                    + "the generic Python interface of ParaMonte will is available on your system, which can \n"
                    + "be used for post-processing and visualization of the output files from \n"
                    + "already-performed ParaMonte simulations or other similar Monte Carlo simulations. \n"
                    + "There are ongoing efforts, right now as you read this message, to further increase the \n"
                    + "availability of ParaMonte library on a wider-variety of platforms and architectures. \n"
                    + "Stay tuned for updates by visiting, \n\n" 
                    + "    https://www.cdslab.org/paramonte\n\n"
                    + "That said,\n\n"
                    + "if your platform is non-Windows and is compatible with GNU Compiler Collection (GCC),\n\n"
                    + "you can also build the required ParaMonte kernel's shared object files on your system\n"
                    + "by calling ParaMonte module's build() function from within your Python environment."
            , marginTop = 1
            , marginBot = 1
            , methodName = _pm.names.paramonte
            )
    return None

####################################################################################################################################

def getBashrcContents():
    bashrcPath = _os.path.expanduser("~/.bashrc")
    if _os.path.isfile(bashrcPath):
        with open(bashrcPath,"r") as bashrcFile:
            bashrcContents = bashrcFile.read()
    else:
        bashrcContents = ""
        with open(bashrcPath,"w") as bashrcFile:
            pass
    return bashrcContents

####################################################################################################################################

def getBashProfileContents():
    bashProfilePath = _os.path.expanduser("~/.bash_profile")
    bashProfileFileExists = _os.path.isfile(bashProfilePath)
    bashProfileContents = ""
    if bashProfileFileExists:
        with open(bashProfilePath,"r") as bashProfileFile:
            bashProfileContents = bashProfileFile.read()
    if ".bashrc" not in bashProfileContents:
        with open(bashProfilePath,"a+") as bashProfileFile:
            bashProfileFile.write("\n[ -f $HOME/.bashrc ] && . $HOME/.bashrc\n")
    return bashProfileContents

####################################################################################################################################

####################################################################################################################################

def setupUnixPath():

    bashrcContents = getBashrcContents()
    dlibcmd = "export LD_LIBRARY_PATH=" + _pm.path.lib + ":$LD_LIBRARY_PATH"
    if dlibcmd not in bashrcContents:
        _os.system( "chmod 777 ~/.bashrc")
        _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
        _os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte shared library setup >>>' >> ~/.bashrc" )
        _os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc" )
        _os.system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc" )
        _os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
        _os.system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc" )
        _os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte shared library setup <<<' >> ~/.bashrc" )
        _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
        _os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

    localInstallDir = getLocalInstallDir()
    if localInstallDir.root is not None:

        pathcmd = None
        dlibcmd = None
        if localInstallDir.gnu.bin is not None: pathcmd = "export PATH=" + localInstallDir.gnu.bin + ":$PATH"
        if localInstallDir.gnu.lib is not None: dlibcmd = "export LD_LIBRARY_PATH=" + localInstallDir.gnu.lib + ":$LD_LIBRARY_PATH"
        if (pathcmd is not None) or (dlibcmd is not None):
            if (pathcmd not in bashrcContents) or (dlibcmd not in bashrcContents):
                _os.system( "chmod 777 ~/.bashrc")
                _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                _os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte local GNU installation setup >>>' >> ~/.bashrc" )
                if pathcmd is not None:
                    if pathcmd not in bashrcContents:
                        _os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${PATH+x} ]; then' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '    export PATH=.' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '" + pathcmd + "' >>  ~/.bashrc" )
                if dlibcmd is not None:
                    if dlibcmd not in bashrcContents:
                        _os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc" )
            if pathcmd not in bashrcContents or dlibcmd not in bashrcContents:
                _os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte local GNU installation setup <<<' >> ~/.bashrc" )
                _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                _os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

        pathcmd = None
        dlibcmd = None
        if localInstallDir.mpi.bin is not None: pathcmd = "export PATH=" + localInstallDir.mpi.bin + ":$PATH"
        if localInstallDir.mpi.lib is not None: dlibcmd = "export LD_LIBRARY_PATH=" + localInstallDir.mpi.lib + ":$LD_LIBRARY_PATH"
        if (pathcmd is not None) or (dlibcmd is not None):
            if (pathcmd not in bashrcContents) or (dlibcmd not in bashrcContents):
                _os.system( "chmod 777 ~/.bashrc")
                _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                _os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte local MPI installation setup >>>' >> ~/.bashrc" )
                if pathcmd is not None:
                    if pathcmd not in bashrcContents:
                        _os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${PATH+x} ]; then' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '    export PATH=.' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '" + pathcmd + "' >>  ~/.bashrc" )
                if dlibcmd is not None:
                    if dlibcmd not in bashrcContents:
                        _os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc" )
            if pathcmd not in bashrcContents or dlibcmd not in bashrcContents:
                _os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte local MPI installation setup <<<' >> ~/.bashrc" )
                _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                _os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

    return None

####################################################################################################################################

class _Struct:
    pass

def getLocalInstallDir():

    localInstallDir = _Struct()
    localInstallDir.root = None

    localInstallDir.mpi = _Struct()
    localInstallDir.mpi.root = None
    localInstallDir.mpi.bin = None
    localInstallDir.mpi.lib = None

    localInstallDir.gnu = _Struct()
    localInstallDir.gnu.root = None
    localInstallDir.gnu.bin = None
    localInstallDir.gnu.lib = None

    localInstallDir.caf = _Struct()
    localInstallDir.caf.root = None
    localInstallDir.caf.bin = None
    localInstallDir.caf.lib = None


    pmGitRootDir = _os.path.join( _pm.path.root , "paramonte-master" )

    if _os.path.isdir(pmGitRootDir):

        localInstallDir.root = pmGitRootDir

        # mpi

        _ = _os.path.join( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "mpich", "3.2" )
        if _os.path.isdir(_):
            localInstallDir.mpi.root = _
            _ = _os.path.join( localInstallDir.mpi.root, "bin" )
            if _os.path.isdir(_): localInstallDir.mpi.bin = _
            _ = _os.path.join( localInstallDir.mpi.root, "lib" )
            if _os.path.isdir(_): localInstallDir.mpi.lib = _

        # gnu

        _ = _os.path.join( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "gnu", "8.3.0" )
        if _os.path.isdir(_):
            localInstallDir.gnu.root = _
            _ = _os.path.join( localInstallDir.gnu.root, "bin" )
            if _os.path.isdir(_): localInstallDir.gnu.bin = _
            _ = _os.path.join( localInstallDir.gnu.root, "lib64" )
            if _os.path.isdir(_): localInstallDir.gnu.lib = _

        # caf

        _ = _os.path.join( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "opencoarrays", "2.8.0" )
        if _os.path.isdir(_):
            localInstallDir.caf.root = _
            _ = _os.path.join( localInstallDir.caf.root, "bin" )
            if _os.path.isdir(_): localInstallDir.caf.bin = _
            _ = _os.path.join( localInstallDir.caf.root, "lib64" )
            if _os.path.isdir(_): localInstallDir.caf.lib = _

    return localInstallDir

####################################################################################################################################

def findMPI():
    """
    Returns MPI bin directory if it exists, otherwise None.
    """

    if isWin32:

        pathList = _os.environ['PATH'].split(";")
        for path in pathList:
            pathLower = path.lower().replace("\\","")
            if ("mpiintel64bin" in pathLower):
                mpivarsFilePath = _os.path.join( path, "mpivars.bat" )
                if _os.path.exists(mpivarsFilePath):
                    mpivarsCommand = '"' + mpivarsFilePath + '"'
                    _pm.note( msg   = "Intel MPI library for 64-bit architecture detected at: \n\n"
                                    + "    " + path + "\n\n"
                                    + "To perform ParaMonte simulations in parallel on a single node, run the \n"
                                    + "following two commands, in the form and order specified, on a Python-aware, \n"
                                    + "mpiexec-aware command-line interface such as Anaconda 3 command-prompt:\n\n"
                                    + "    " + mpivarsCommand + "\n\n"
                                    + "    mpiexec -localonly -n NUM_PROCESSES python main.py\n\n"
                                    + "where, \n\n"
                                    + "    0.   the first command defines the essential environment variables and, \n"
                                    + "         the second command runs in the simulation in parallel, in which, \n"
                                    + "    1.   you should replace NUM_PROCESSES with the number of processors you \n"
                                    + "         wish to assign to your simulation task, \n"
                                    + "    2.   -localonly indicates a parallel simulation on only a single node (this \n"
                                    + "         flag will obviate the need for MPI library credentials registration). \n"
                                    + "         For more information, visit: \n"
                                    + "         https://software.intel.com/en-us/get-started-with-mpi-for-windows \n"
                                    + "    3.   main.py is the Python file which serves as the entry point to \n"
                                    + "         your simulation, where you call ParaMonte sampler routines. \n\n"
                                    + "Note that the above two commands must be executed on a command-line that recognizes \n"
                                    + "both Python and mpiexec applications, such as the Anaconda command-line interface. \n"
                                    + "For more information, in particular, on how to register to run Hydra services \n"
                                    + "for multi-node simulations on Windows servers, visit:\n\n"
                                    + "https://www.cdslab.org/pm"
                            , marginTop = 1
                            , marginBot = 1
                            , methodName = _pm.names.paramonte
                            )

                    setupFilePath = _os.path.join( _pm.path.root, "setup.bat" )
                    with open(setupFilePath, "w") as setupFile:
                        setupFile.write("@echo off\n")
                        setupFile.write("cd " + path + " && mpivars.bat quiet\n")
                        setupFile.write("cd " + _pm.path.root + "\n")
                        setupFile.write("@echo on\n")

                    return path
                    break

    elif isLinux:

        pathList = _os.environ['PATH'].split(":")
        for path in pathList:
            pathLower = path.lower().replace("/","")
            if ("linuxmpiintel64" in pathLower):
                mpivarsFilePath = _os.path.join( path, "mpivars.sh" )
                if _os.path.exists(mpivarsFilePath):
                    mpivarsCommand = '"' + mpivarsFilePath + '"'
                    _pm.note( msg   = "Intel MPI library for 64-bit architecture detected at: \n\n"
                                    + "    " + path + "\n\n"
                                    + "To perform ParaMonte simulations in parallel on a single node, run the \n"
                                    + "following two commands, in the form and order specified, in a Bash shell, \n\n"
                                    + "    source " + mpivarsCommand + "\n\n"
                                    + "    mpiexec -n NUM_PROCESSES python main.py\n\n"
                                    + "where, \n\n"
                                    + "    0.   the first command defines the essential environment variables and, \n"
                                    + "         the second command runs in the simulation in parallel, in which, \n"
                                    + "    1.   you should replace NUM_PROCESSES with the number of processors you \n"
                                    + "         wish to assign to your simulation task, \n"
                                    + "    2.   main.py is the Python file which serves as the entry point to \n"
                                    + "         your simulation, where you call ParaMonte sampler routines. \n\n"
                                    + "For more information on how to install and use and run parallel ParaMonte \n"
                                    + "simulations on Linux systems, visit:\n\n"
                                    + "https://www.cdslab.org/pm"
                            , marginTop = 1
                            , marginBot = 1
                            , methodName = _pm.names.paramonte
                            )

                    setupFilePath = _os.path.join( _pm.path.root, "setup.sh" )
                    with open(setupFilePath, "w") as setupFile:
                        setupFile.write("source " + mpivarsCommand)

                    return path
                    break

    elif isMacOS:

        import shutil

        gfortranPath = None
        try:
            import subprocess
            gfortranVersion = subprocess.run(args=["gfortran", "--version"],capture_output=True)
            if "GCC 9." in str(gfortranVersion.stdout): gfortranPath = shutil.which("gfortran")
        except:
            pass

        mpiexecPath = None
        try:
            import subprocess
            mpiexecVersion = subprocess.run(args=["mpiexec", "--version"],capture_output=True)
            if "open-mpi" in str(mpiexecVersion.stdout): mpiexecPath = shutil.which("mpiexec")
        except:
            pass

        if (mpiexecPath is not None) and (gfortranPath is not None):
            path = _os.path.dirname(mpiexecPath)
            _pm.note( msg   = "MPI runtime libraries detected at: \n\n"
                    + "    " + path + "\n\n"
                    + "To perform ParaMonte simulations in parallel on a single node, run the \n"
                    + "following command, in the form and order specified, in a Bash shell, \n\n"
                    + "    mpiexec -n NUM_PROCESSES python main.py\n\n"
                    + "where, \n\n"
                    + "    0.   the first command defines the essential environment variables and, \n"
                    + "         the second command runs in the simulation in parallel, in which, \n"
                    + "    1.   you should replace NUM_PROCESSES with the number of processors you \n"
                    + "         wish to assign to your simulation task, \n"
                    + "    2.   main.py is the Python file which serves as the entry point to \n"
                    + "         your simulation, where you call ParaMonte sampler routines. \n\n"
                    + "For more information on how to install and use and run parallel ParaMonte \n"
                    + "simulations on Darwin operating systems, visit:\n\n"
                    + "https://www.cdslab.org/pm"
            , marginTop = 1
            , marginBot = 1
            , methodName = _pm.names.paramonte
            )
            return path

    else:

        LocalInstallDir = getLocalInstallDir()
        if (LocalInstallDir.mpi.bin is not None) and (LocalInstallDir.mpi.lib is not None):
            return LocalInstallDir.mpi.bin

    return None

####################################################################################################################################

def installMPI():

    if isWin32 or isLinux:

        _pm.note( msg = "Downloading the Intel MPI runtime libraries for 64-bit architecture...\n"
                      + "Please make sure your firewall allows access to the Internet."
                , methodName = _pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        if isWin32:
            mpiFileExt = ".exe"
            mpiVersion = "2019.4.245"
            mpiFileName = "w_mpi-rt_p_" + mpiVersion
            downloadList = ["impi.dll","impi.pdb","libfabric.dll"]
            for mpiFileNameExt in downloadList:
                mpiFilePath = _os.path.join( _pm.path.root, mpiFileNameExt )
                download( url = "https://github.com/cdslaborg/paramonte/releases/download/" + _pm.version.kernel.dump() + "/" + mpiFileNameExt
                        , filePath = mpiFilePath
                        )
            
        if isLinux:
            mpiFileExt = ".tgz"
            mpiVersion = "2018.2.199"
            mpiFileName = "l_mpi-rt_" + mpiVersion

        mpiFileNameExt = mpiFileName + mpiFileExt
        mpiFilePath = _os.path.join( _pm.path.root, mpiFileNameExt )
        download( url = "https://github.com/cdslaborg/paramonte/releases/download/" + _pm.version.kernel.dump() + "/" + mpiFileNameExt
                , filePath = mpiFilePath
                )

        _pm.note( msg = "Installing the Intel MPI library for 64-bit architecture...\n"
                      + "file location: " + mpiFilePath
                , methodName = _pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        _pm.warn( msg = "Please do not change the default installation location of the MPI library suggested by the installer.\n"
                      + "If you do change the default path, the onus will be on you to ensure the path to the \n"
                      + "MPI runtime libraries exist in the environmental PATH variable of your session."
                , methodName = _pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        currentDir = _os.getcwd()

        if isWin32:

            err = 0
            err = _os.system(mpiFilePath)
            if err==0:
                writeVerificationStatusFile(verificationEnabled=True)
                _pm.note( msg   = "Intel MPI library installation appears to have succeeded. \n"
                                + "Now close your Python environment and the command-line interface \n"
                                + "and reopen a new fresh (Anaconda) command prompt."
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )
            else:
                _pm.warn( msg = "Intel MPI library installation might have failed. Exit flag: {}.".format(err)
                    , methodName = _pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )

        if isLinux:

            try:

                import tarfile
                tf = tarfile.open(mpiFilePath)
                tf.extractall(path=_pm.path.root)
                mpiExtractDir = _os.path.join( _pm.path.root, mpiFileName )

                _pm.note( msg   = "If needed, use the following serial number when asked by the installer: \n\n"
                                + "    C44K-74BR9CFG\n\n"
                                + "If this is your personal computer, choose \n\n"
                                + "    'install as root'\n\n"
                                + "in the graphical user interface that appears in your session. \n"
                                + "Otherwise, if you are using ParaMonte on a public server, \n"
                                + "for example, on a supercomputer, choose the third option:\n\n"
                                + "   'install as current user'"
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

                mpiInstallScriptPath = _os.path.join( mpiExtractDir, "install_GUI.sh" )
                if not _os.path.exists(mpiInstallScriptPath):
                    _pm.abort   ( msg   = "Internal error occurred.\n"
                                        + "Failed to detect the Intel MPI installation Bash script.\n"
                                        + "Please report this issue at \n\n"
                                        + "    https://github.com/cdslaborg/paramonte/issues\n\n"
                                        + "Visit https://www.cdslab.org/pm for instructions \n"
                                        + "to build ParaMonte object files on your system."
                                , methodName = _pm.names.paramonte
                                , marginTop = 1
                                , marginBot = 1
                                )

            except Exception as e:

                print(str(e))
                _pm.abort   ( msg   = "Unzipping of Intel MPI runtime library tarball failed.\n"
                                    + "Make sure you have tar software installed on your system and try again."
                            , methodName = _pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

            err = _os.system("chmod +x " + mpiInstallScriptPath)
            if err != 0:
                _pm.warn( msg   = "The following action failed:\n\n"
                                + "chmod +x " + mpiInstallScriptPath + "\n\n"
                                + "skipping..."
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

            _os.chdir(mpiExtractDir)

            import subprocess
            try:
                subprocess.check_call( mpiInstallScriptPath, shell = True )
            except Exception as e:
                print(str(e))
                _pm.abort   ( msg   = "Intel MPI runtime libraries installation for \n"
                                    + "64-bit architecture appears to have failed.\n"
                                    + "Please report this error at:\n\n"
                                    + "    https://github.com/cdslaborg/paramonte/issues\n\n"
                                    + "Visit https://www.cdslab.org/pm for more instructions \n"
                                    + "to build and use ParaMonte on your system."
                            , methodName = _pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )


            _pm.note( msg = "Intel MPI runtime libraries installation for \n"
                          + "64-bit architecture appears to have succeeded.\n"
                          + "Searching for the MPI runtime environment setup file..."
                    , methodName = _pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )

            _os.chdir(currentDir)

            setupFilePath = _os.path.join( _pm.path.root, "setup.sh" )

            mpiRootDirNotFound = True
            installationRootDirList = [ "/opt", _pm.path.home ]
            while mpiRootDirNotFound:

                mpiRootDir = []
                mpivarsFilePathDefault = []
                for installationRootDir in installationRootDirList:
                    mpiTrunkDir = _os.path.join( "intel", "compilers_and_libraries_" + mpiVersion, "linux", "mpi", "intel64" )
                    mpiRootDir.append( _os.path.join( installationRootDir, mpiTrunkDir ) )
                    mpivarsFilePathDefault.append( _os.path.join( mpiRootDir[-1] , "bin" , "mpivars.sh" ) )
                    if _os.path.isdir(mpiRootDir[-1]):
                        mpiRootDirNotFound = False
                        break

                if mpiRootDirNotFound:
                    _pm.warn( msg = "Failed to detect the installation root path for Intel MPI runtime libraries \n"
                                  + "for 64-bit architecture on your system. If you specified a different installation \n"
                                  + "root path at the time installation, please copy and paste it below.\n"
                                  + "Note that the installation root path is part of the path that replaces: \n\n"
                                  + "    " + "opt" + "\n\n"
                                  + "in the following path: \n\n"
                                  + "    " + _os.path.join( "opt" , mpiTrunkDir )
                            , methodName = _pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )
                    answer = input  ( "\n    Please type the root path of MPI installation below and press ENTER."
                                    + "\n    If you don't know the root path, simply press ENTER to quit:\n"
                                    )
                    if len(answer.strip())==0:
                        _pm.warn( msg = "Skipping the MPI runtime library environmental path setup..."
                            , methodName = _pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )
                        break
                    else:
                        installationRootDirList = [ answer ]
                        continue

            if mpiRootDirNotFound:

                _pm.warn( msg   = "Failed to find the MPI runtime environment setup file on your system.\n"
                                + "This is highly unusual. Normally, Intel MPI libraries is installed\n"
                                + "in the following directory:\n\n"
                                + "    " + mpiRootDir[0] + "\n\n"
                                + "or,\n\n"
                                + "    " + mpiRootDir[1] + "\n\n"
                                + "If you cannot manually find the Intel MPI installation directory,\n"
                                + "it is likely that the installation might have somehow failed.\n"
                                + "If you do find the installation directory, try to locate the\n"
                                + "'mpivars.sh' file which is normally installed in the following path:\n\n"
                                + "    " + mpivarsFilePathDefault[0] + "\n\n"
                                + "or,\n\n"
                                + "    " + mpivarsFilePathDefault[1] + "\n\n"
                                + "Before attempting to run any parallel ParaMonte simulation, \n"
                                + "make sure you source this file, like the following:\n\n"
                                + "    source " + mpivarsFilePathDefault[0] + "\n\n"
                                + "or,\n\n"
                                + "    source " + mpivarsFilePathDefault[1] + "\n\n"
                                + "where you will have to replace the path in the above with the \n"
                                + "correct path that you find on your system."
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

            else:

                mpiBinDir = _os.path.join( mpiRootDir[-1], "bin" )
                mpiLibDir = _os.path.join( mpiRootDir[-1], "lib" )
                mpivarsFilePath = _os.path.join( mpiBinDir, "mpivars.sh" )
                if _os.path.isfile(mpivarsFilePath):

                    with open(setupFilePath, "w") as setupFile:
                        setupFile.write(mpiBinDir+"\n")
                        setupFile.write(mpiLibDir+"\n")
                        setupFile.write("source " + mpivarsFilePath)

                    _pm.note( msg = "To ensure all MPI routine environmental variables \n"
                                  + "are properly load, source the following Bash script \n"
                                  + "in your Bash environment before calling mpiexec, like:\n\n"
                                  + "    source " + mpivarsFilePath + "\n\n"
                                  + "Alternatively, ParaMonte can also automatically add \n"
                                  + "the required script to your '.bashrc' file, so that \n"
                                  + "all required MPI environmental variables are loaded \n"
                                  + "automatically before any ParaMonte usage from any \n"
                                  + "Bash command line on your system."
                            , methodName = _pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

                    #isYes = getUserResponse ( msg = "\n    MPI runtime variables are essential for parallel ParaMonte"
                    #                                "\n    simulation. Would you like ParaMonte to add the MPI runtime"
                    #                                "\n    environmental variables to your Bash environment (y/n)? " 
                    #                                )
                    #
                    #if isYes:

                    bashrcContents = getBashrcContents()
                    mpivarsFileCommand = "source " + mpivarsFilePath
                    if mpivarsFileCommand not in bashrcContents:
                        _os.system( "chmod 777 ~/.bashrc")
                        _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte MPI runtime library initialization >>>' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '" + mpivarsFileCommand + "' >>  ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte MPI runtime library initialization <<<' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

                        _pm.note( msg = "If you intend to run parallel simulations right now,\n"
                                      + "we highly recommned you to close your current shell environment\n"
                                      + "and open a new Bash shell environment. This is to ensure that all MPI\n"
                                      + "library environmental variables are properly set in your shell environment."
                                , methodName = _pm.names.paramonte
                                , marginTop = 1
                                , marginBot = 1
                                )

                    #else:
                    #    _pm.warn( msg = "skipping...\n"
                    #                  + "It is now your responsibility to ensure that the MPI runtime \n"
                    #                  + "environmental variables in your Bash environment are properly \n"
                    #                  + "set up before attempting to run any parallel ParaMonte simulation. \n"
                    #                  + "You can do so by running the following command in every Bash session:\n\n"
                    #                  + "    source " + mpivarsFilePath + "\n\n"
                    #                  + "Alternatively, ParaMonte can also automatically add \n"
                    #                  + "the required script to your '.bashrc' file, so that \n"
                    #                  + "all required MPI environmental variables are loaded \n"
                    #                  + "automatically before any ParaMonte usage from any \n"
                    #                  + "Bash command line on your system."
                    #            , methodName = _pm.names.paramonte
                    #            , marginTop = 1
                    #            , marginBot = 1
                    #            )

                else:

                    _pm.abort   ( msg   = "ParaMonte was able to detect an MPI library path on your system, however,\n"
                                        + "the MPI installation appears to be corrupted. The required mpivars.sh \n"
                                        + "does not exist:\n\n"
                                        + mpivarsFilePath
                                , methodName = _pm.names.paramonte
                                , marginTop = 1
                                , marginBot = 1
                                )
                    self.Err.abort();

    elif isMacOS:

        _pm.warn( msg   = "To use the ParaMonte kernel routines in parallel on macOS, \n"
                        + "the Open-MPI library will have to be installed on your system. \n"
                        #+ "To ensure full consistency, we recommend building the parallel \n"
                        #+ "object files of ParaMonte library on your system along with Open-MPI.\n\n"
                        + "\n"
                        + "If this installation of the prerequisites is being done from within \n"
                        + "a Jupyter notebook and the installation fails:\n"
                        + "\n"
                        + "    1. quit the Jupyter notebook.\n"
                        + "    2. enter an IPython session on the command-prompt:\n"
                        + "        - On Windows, use Anaconda3 command-prompt.\n"
                        + "        - On Linux / macOS, use the Bash terminal.\n"
                        + "    3. import paramonte as pm\n"
                        + "    4. pm.verify()\n"
                        + "\n"
                        + "Building the ParaMonte library prerequisites on your system..."
                , marginTop = 1
                , marginBot = 1
                , methodName = _pm.names.paramonte
                )
        buildParaMontePrereqsForMac()

    else:

        _pm.warn( msg   = "To use ParaMonte in parallel on this unknown Operating System, \n"
                        + "ParaMonte needs to be built from scratch on your system. \n"
                        + "Building ParaMonte library prerequisites on your system..."
                , marginTop = 1
                , marginBot = 1
                , methodName = _pm.names.paramonte
                )
        build()

####################################################################################################################################

def buildParaMontePrereqsForMac():

    _pm.note( msg = "Checking if Homebrew exists on your system..."
            , methodName = _pm.names.paramonte
            , marginTop = 1
            , marginBot = 1
            )

    import shutil
    import subprocess
    if shutil.which("brew") is None:

        _pm.note( msg = "Failed to detect Homebrew on your system. Installing Homebrew..."
                , methodName = _pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )
        err1 = _os.system('xcode-select --install')
        err2 = _os.system('/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"')
        err3 = _os.system('brew --version')
        if err1 != 0 or err2 != 0 or err3 != 0:
            _pm.abort( msg  = "Failed to install Homebrew on your system.\n"
                            + "Homebrew is required to install and build ParaMonte components and prerequisites.\n"
                            + "Please install Homebrew manually on your system and retry the ParaMonte installation process.\n"
                            + "skipping..."
                    , methodName = _pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )

    # cmake

    cmakeInstallationNeeded = False
    cmakePath = shutil.which("cmake")
    if cmakePath is None:
        cmakeInstallationNeeded = True
        _pm.note( msg           = "cmake installation is missing on your system."
                , methodName    = _pm.names.paramonte
                , marginTop     = 1
                , marginBot     = 1
                )
    else:
        _pm.note( msg           = "cmake installation detected at: " + cmakePath + "\n" + "Checking cmake version..."
                , methodName    = _pm.names.paramonte
                , marginTop     = 0
                , marginBot     = 0
                )
        try:
            cmakeVersion = str(subprocess.run(args=["cmake","--version"],capture_output=True).stdout).split(" ")[2].split("-")[0]
            cmakeVersionList = cmakeVersion.split(".")
            _pm.note( msg           = "current cmake version: " + cmakeVersion
                    , methodName    = _pm.names.paramonte
                    , marginTop     = 0
                    , marginBot     = 0
                    )
            if int(cmakeVersionList[0])>=3 and int(cmakeVersionList[1])>=14:
                _pm.note( msg           = "cmake version is ParaMonte-compatible!"
                        , methodName    = _pm.names.paramonte
                        , marginTop     = 0
                        , marginBot     = 0
                        )
            else:
                cmakeInstallationNeeded = True
                _pm.note( msg           = "cmake version is NOT ParaMonte-compatible."
                        , methodName    = _pm.names.paramonte
                        , marginTop     = 0
                        , marginBot     = 0
                        )
        except:
            cmakeInstallationNeeded = True
            _pm.note( msg           = "Failed to detect the current cmake installation version. skipping..."
                    , methodName    = _pm.names.paramonte
                    , marginTop     = 0
                    , marginBot     = 0
                    )

    if cmakeInstallationNeeded:

        _pm.note( msg = "Installing cmake..."
                , methodName = _pm.names.paramonte
                , marginTop = 0
                , marginBot = 0
                )

        err1 = _os.system("brew install cmake")
        err2 = _os.system("brew link --overwrite cmake")

        if err1 != 0 or err2 != 0:
            _pm.abort   ( msg = "cmake installation or linking failed."
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

        cmakeVersionList = str(subprocess.run(args=["cmake","--version"],capture_output=True).stdout).split(" ")[2].split("-")[0].split(".")
        if int(cmakeVersionList[0])>=3 and int(cmakeVersionList[1])>=14:
            _pm.note( msg           = "cmake installation succeeded."
                    , methodName    = _pm.names.paramonte
                    , marginTop     = 1
                    , marginBot     = 1
                    )
        else:
            _pm.warn( msg   = "Failed to install and link cmake on your system.\n"
                            + "cmake is required to install and build\n"
                            + "ParaMonte components and prerequisites.\n"
                            + "Please install the cmake manually on your system and\n"
                            + "retry the ParaMonte installation process if it fails.\n"
                            + "skipping..."
                    , methodName = _pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )

    # gnu

    _pm.note( msg = "Installing GNU Compiler Collection..."
            , methodName = _pm.names.paramonte
            , marginTop = 1
            , marginBot = 1
            )

    err1 = _os.system("brew install gcc@9")
    err2 = _os.system("brew link gcc@9")

    if err1 != 0 or err2 != 0:
        _pm.warn( msg   = "Failed to install and link GNU Compiler Collection on your system.\n"
                        + "The GNU Compiler Collection is required to install\n"
                        + "and build ParaMonte components and prerequisites.\n"
                        + "Please install the GNU Compiler Collection manually on your\n"
                        + "system and retry the ParaMonte installation process if it fails.\n"
                        + "skipping..."
                , methodName = _pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

    # open-mpi

    _pm.note( msg = "Installing Open-MPI..."
            , methodName = _pm.names.paramonte
            , marginTop = 1
            , marginBot = 1
            )

    err1 = _os.system("brew install open-mpi")
    err2 = _os.system("brew link open-mpi")

    if err1 != 0 or err2 != 0:
        _pm.warn( msg   = "Failed to install and link Open-MPI on your system.\n"
                        + "Open-MPI is required to install and build\n"
                        + "ParaMonte components and prerequisites.\n"
                        + "Please install the Open-MPI manually on your\n"
                        + "system and retry the ParaMonte installation process if it fails.\n"
                        + "skipping..."
                , methodName = _pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )


####################################################################################################################################

def writeVerificationStatusFile(verificationEnabled):
    with open(verificationStatusFilePath, "w") as verificationStatusFile:
        verificationStatusFile.write(str(verificationEnabled))
    return None

####################################################################################################################################

def dispFinalMessage():
    _pm.note( msg   = "To check for the MPI library installation status or display the above messages \n"
                    + "in the future, use the following commands on the Python command-line: \n\n"
                    + "    import paramonte as pm\n"
                    + "    pm.verify()"
            , methodName = _pm.names.paramonte
            , marginTop = 1
            , marginBot = 1
            )
    return None

####################################################################################################################################

def displayParaMonteBanner():
    bannerFilePath = _os.path.join( _pm.path.auxil, ".ParaMonteBanner")
    offset = ( len(_pm.version.interface.dump()) - 5 ) // 2
    print("")
    with open(bannerFilePath,"r") as file:
        for line in file:
            if "Version" in line:
                line = line.replace(" "*offset+"Version 0.0.0","Version "+_pm.version.interface.dump())
            print(line,end="")
    print("")
    return None

####################################################################################################################################

def build(flags=""):

    if isWin32:

        _pm.abort   ( msg   = "ParaMonte library build on Windows Operating Systems (OS) requires \n"
                            + "the installation of the following software on your system:\n\n"
                            + "    -- Microsoft Visual Studio (MSVS) (Community Edition >2017)\n"
                            + "    -- Intel Parallel Studio >2018, which is built on top of MSVS\n\n"
                            + "If you don't have these software already installed on your system, \n"
                            + "please visit the following page for the installation instructions:\n\n"
                            + "    https://www.cdslab.org/pm\n\n"
                            + "Follow the instructions on this website for building ParaMonte on your system."
                    , methodName = _pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )

    else:

        _pm.note( msg   = "You are requesting to build the ParaMonte kernel libraries on your system. \n"
                        + "The kernel library build requires ParaMonte-compatible versions of the following \n"
                        + "compilers and parallelism libraries to be installed on your system: \n\n"
                        + "    GNU compiler collection (GCC >8.3)\n"
                        + "    MPI library (MPICH >3.2) on Linux OS or Open-MPI on Darwin OS\n"
                        + "    OpenCoarrays >2.8\n\n"
                        + "The full installation of these software will require 4 to 5 Gb of free space \n"
                        + "on your system (where the ParaMonte-Python interface is already installed).\n"
                        + "Note that the installation script is in Bash and therefore requires Bash shell.\n"
                        + "An existing recent installation of the GNU Compiler Collection (GCC) on your\n"
                        + "system would be also highly desirable and will significantly cut the build time.\n"
                        + "Also, downloading the prerequisites requires access to the Internet.\n"
                        + "If you have an Internet firewall active on your system, please make sure to\n"
                        + "turn it off before proceeding with the local installation of ParaMonte."
                , methodName = _pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        buildEnabled = getUserResponse  ( msg   = "\n    Do you wish to download and install the ParaMonte library"
                                                + "\n    and its prerequisites on your system now (y/n)? " 
                                        )


        if buildEnabled:

            if isMacOS: buildParaMontePrereqsForMac()

            currentDir = _os.getcwd()

            pmGitTarPath = _os.path.join( _pm.path.root, "master.tar.gz" )
            download( url = "https://github.com/cdslaborg/paramonte/archive/master.tar.gz"
                    , filePath = pmGitTarPath
                    )

            pmGitRootDir = _os.path.join( _pm.path.root, "paramonte-master" )

            try:

                import tarfile
                tf = tarfile.open(pmGitTarPath)
                tf.extractall(path=_pm.path.root) # path=pmGitRootDir)

                pmGitInstallScriptPath = _os.path.join( pmGitRootDir, "install.sh" )
                if not _os.path.exists(pmGitInstallScriptPath):
                    _pm.abort   ( msg   = "Internal error occurred.\n"
                                        + "Failed to detect the ParaMonte installation Bash script.\n"
                                        + "Please report this issue at \n\n"
                                        + "    https://github.com/cdslaborg/paramonte/issues\n\n"
                                        + "Visit https://www.cdslab.org/pm for instructions \n"
                                        + "to build ParaMonte object files on your system."
                                , methodName = _pm.names.paramonte
                                , marginTop = 1
                                , marginBot = 1
                                )

            except Exception as e:

                print(str(e))
                _pm.abort   ( msg   = "Unzipping of the ParaMonte tarball failed.\n"
                                    + "Make sure you have tar software installed on your system and try again."
                            , methodName = _pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

            err = _os.system("chmod +x " + pmGitInstallScriptPath)
            if err != 0:
                _pm.warn( msg   = "The following action failed:\n\n"
                                + "chmod +x " + pmGitInstallScriptPath + "\n\n"
                                + "skipping..."
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

            _os.chdir(pmGitRootDir)

            import subprocess
            try:
                _os.system( "find " + pmGitRootDir + " -type f -iname \"*.sh\" -exec chmod +x {} \;" )
                _os.system( pmGitInstallScriptPath + " --lang python --test_enabled true --exam_enabled false --yes-to-all " + flags )
            except Exception as e:
                print(str(e))
                _pm.abort   ( msg   = "Local installation of ParaMonte failed.\n"
                                    + "Please report this issue at \n\n"
                                    + "    https://github.com/cdslaborg/paramonte/issues"
                            , methodName = _pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

            _os.chdir(currentDir)

            # copy files to module folder

            import glob
            import shutil
            pythonBinDir = _os.path.join( pmGitRootDir , "bin" , "Python" , "paramonte" )
            fileList = glob.glob( _os.path.join( pythonBinDir , "libparamonte_*" ) )

            if len(fileList)==0:

                _pm.abort( msg  = "ParaMonte kernel libraries build and installation appears to have failed. \n"
                                + "You can check this path:\n\n"
                                + pythonBinDir + "\n\n"
                                + "to find out if any shared objects with the prefix 'libparamonte_' have been generated or not.\n"
                                + "Please report this issue at \n\n"
                                + "    https://github.com/cdslaborg/paramonte/issues"
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 2
                        )

            else:

                _pm.note( msg   = "ParaMonte kernel libraries build appears to have succeeded. \n"
                                + "copying the kernel files to the ParaMonte Python module directory..."
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

                for file in fileList:
                    _pm.note( msg   = "file: " + file
                            , methodName = _pm.names.paramonte
                            , marginTop = 0
                            , marginBot = 0
                            )
                    shutil.copy(file, _pm.path.lib)

                _pm.note( msg   = "ParaMonte kernel libraries should be now usable on your system."
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

                setupFilePath = _os.path.join( pmGitRootDir , "build", "prerequisites", "prerequisites", "installations", "opencoarrays", "2.8.0", "setup.sh" )

                if _os.path.exists(setupFilePath):

                    bashrcContents = getBashrcContents()
                    setupFilePathCmd = "source " + setupFilePath
                    if setupFilePathCmd not in bashrcContents:
                        _os.system( "chmod 777 ~/.bashrc")
                        _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte library local installation setup >>>' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '" + setupFilePathCmd + "' >>  ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte library local installation setup <<<' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

                    _pm.warn( msg   = "Whenever you intend to use ParaMonte in the future, before opening your Python session, \n"
                                    + "please execute the following command in your Bash shell to ensure all required paths \n"
                                    + "are properly defined in your environment:\n\n"
                                    + "    " + setupFilePathCmd + " \n\n"
                                    + "mission accomplished."
                            , methodName = _pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

            writeVerificationStatusFile(verificationEnabled=True)

        else:

            _pm.warn( msg   = "Aborting the ParaMonte-for-Python local build on your system."
                    , methodName = _pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )

    return None

####################################################################################################################################

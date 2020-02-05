#**********************************************************************************************************************************
#**********************************************************************************************************************************
#
#  ParaMonte: plain powerful parallel Monte Carlo library.
#
#  Copyright (C) 2012-present, The Computational Data Science Lab
#
#  This file is part of ParaMonte library.
#
#  ParaMonte is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, version 3 of the License.
#
#  ParaMonte is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
#
#**********************************************************************************************************************************
#**********************************************************************************************************************************

import os as _os
import sys as _sys
import _paramonte as _pm
import warnings as _warnings

fileAbsDir = _os.path.dirname(_os.path.abspath(__file__))
verificationStatusFilePath = _os.path.join( fileAbsDir, "verificationEnabled.txt" )

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

    if reset: writeVerificationStatusFile(verified=True)

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

            displayParaMontePythonBanner()

            # search for the MPI library

            mpiBinPath = findMPI()
            if mpiBinPath is None:
                _pm.warn( msg   = "The MPI runtime libraries for 64-bit architecture could not be detected on your system. \n"
                                + "The MPI runtime libraries are required for parallel simulations via ParaMonte library. \n"
                                + "You can download and install the proper versions of an MPI library free of charge \n"
                                + "from either Intel website (https://software.intel.com/en-us/mpi-library) \n"
                                + "for Windows and Linux operating systems) or from MPICH website (www.mpich.org) \n"
                                + "for Mac operating systems. Alternatively, ParaMonte can install these \n"
                                + "library for you automatically now. \n\n"
                                + "If you don't know how to download and install the correct MPI runtime library version, \n"
                                + "we strongly recommend that you let ParaMonte install this library for you now. \n"
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

            writeVerificationStatusFile(verified=False)
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
    _pm.warn( msg   = "The ParaMonte C/Fortran sampler kernel is currently exclusively available \n"
                    + "on AMD64 (64-bit) architecture for Windows/Linux/MacOS Operating Systems (OS). \n"
                    + "Your system appears to be of a different architecture or OS. As a result, \n"
                    + "the core sampler routines of ParaMonte will not be available on your system. However, \n"
                    + "the generic Python interface of ParaMonte will is available on your system, which can \n"
                    + "be used for post-processing and visualization of the output files from \n"
                    + "already-performed ParaMonte simulations or other similar Monte Carlo simulations. \n"
                    + "There are ongoing efforts, right now as you read this message, to further increase the \n"
                    + "availability of ParaMonte library on a wider-variety of platforms and architectures. \n"
                    + "Stay tuned for updates by visiting https://www.cdslab.org/pm\n\n"
                    + "That said,\n\n"
                    + "if your platform is non-Windows and is compatible with GNU Compiler Collection (GCC),\n\n"
                    + "you can also build the required ParaMonte kernel's shared object files on your system\n"
                    + "by calling ParaMonte module's build() function from within your Python enviroment."
            , marginTop = 1
            , marginBot = 1
            , methodName = _pm.names.paramonte
            )
    return None

####################################################################################################################################

class _Struct:
    pass

def getLocalInstallDir():

    localInstallDir = None
    pmGitRootDir = _os.path.join( fileAbsDir , "paramonte-master" )
    if _os.path.isdir(pmGitRootDir):

        localInstallDir = _Struct()
        localInstallDir.root = pmGitRootDir

        # mpi

        localInstallDir.mpi = _Struct()
        localInstallDir.mpi.root = None
        localInstallDir.mpi.bin = None
        localInstallDir.mpi.lib = None
        _ = _os.path.join( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "mpich", "3.2" )
        if _os.path.isdir(_):
            localInstallDir.mpi.root = _
            _ = _os.path.join( localInstallDir.mpi.root, "bin" )
            if _os.path.isdir(_): localInstallDir.mpi.bin = _
            _ = _os.path.join( localInstallDir.mpi.root, "lib" )
            if _os.path.isdir(_): localInstallDir.mpi.lib = _

        # gnu

        localInstallDir.gnu = _Struct()
        localInstallDir.gnu.root = None
        localInstallDir.gnu.bin = None
        localInstallDir.gnu.lib = None
        _ = _os.path.join( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "gnu", "8.3.0" )
        if _os.path.isdir(_):
            localInstallDir.gnu.root = _
            _ = _os.path.join( localInstallDir.gnu.root, "bin" )
            if _os.path.isdir(_): localInstallDir.gnu.bin = _
            _ = _os.path.join( localInstallDir.gnu.root, "lib64" )
            if _os.path.isdir(_): localInstallDir.gnu.lib = _

        # caf

        localInstallDir.caf = _Struct()
        localInstallDir.caf.root = None
        localInstallDir.caf.bin = None
        localInstallDir.caf.lib = None
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
                    mpivarsCommand = '"' + mpivarsFilePath + '" quiet'
                    _pm.note( msg   = "Intel MPI library for 64-bit architecture detected at: \n\n"
                                    + "    " + path + "\n\n"
                                    + "To perform ParaMonte simulations in parallel on a single node, run the \n"
                                    + "following two commands, in the form and order specified, on a Python-aware, \n"
                                    + "mpiexec-aware command-line interface such as Anaconda 3 command-prompt:\n\n"
                                    + "    " + mpivarsCommand + "\n\n"
                                    + "    mpiexec -localonly -np NUM_PROCESSES python main.py\n\n"
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

                    setupFilePath = _os.path.join( fileAbsDir, "setup.bat" )
                    with open(setupFilePath, "w") as setupFile:
                        setupFile.write("call " + mpivarsCommand)

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
                                    + "    mpiexec -np NUM_PROCESSES python main.py\n\n"
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

                    setupFilePath = _os.path.join( fileAbsDir, "setup.sh" )
                    with open(setupFilePath, "w") as setupFile:
                        setupFile.write("source " + mpivarsCommand)

                    return path
                    break

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
            
        if isLinux:
            mpiFileExt = ".tgz"
            mpiVersion = "2018.2.199"
            mpiFileName = "l_mpi-rt_" + mpiVersion

        mpiFileNameExt = mpiFileName + mpiFileExt
        mpiFilePath = _os.path.join( fileAbsDir, mpiFileNameExt )
        download( url = "https://github.com/cdslaborg/paramonte/releases/download/" + _pm.version + "/" + mpiFileNameExt
                , filePath = mpiFilePath
                )

        _pm.note( msg = "Installing the Intel MPI library for 64-bit architecture...\n"
                      + "file location: "  + mpiFilePath
                , methodName = _pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        currentDir = _os.getcwd()

        if isWin32:

            err = 0
            err = _os.system(mpiFilePath)
            if err==0:
                writeVerificationStatusFile(verified=True)
                _pm.note( msg   = "Intel MPI library installation appears to have succeeded. \n"
                                + "A reboot of your system or restart of your Python environment may \n"
                                + "be required before using the parallel features of ParaMonte."
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
                tf.extractall(path=fileAbsDir)
                mpiExtractDir = _os.path.join( fileAbsDir, mpiFileName )

                _pm.note( msg   = "Use the following serial number when asked by the installer: C44K-74BR9CFG"
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

                mpiInstallScriptPath = _os.path.join( mpiExtractDir, "install.sh" )
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

                raise CritError(messages.crit_error_bad_command+" "+str(e))
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
                raise CritError(messages.crit_error_bad_command+" "+str(e))

            _os.chdir(currentDir)

            setupFilePath = _os.path.join( fileAbsDir, "setup.sh" )
            mpiRootDir = _os.path.join( "opt", "intel", "compilers_and_libraries_" + mpiVersion, "linux", "mpi", "intel64" )
            mpiBinDir = _os.path.join( mpiRootDir, "bin" )
            mpiLibDir = _os.path.join( mpiRootDir, "lib" )
            mpivarsFilePath = _os.path.join( mpiBinDir, "mpivars.sh" )
            with open(setupFilePath, "w") as setupFile:
                setupFile.write(mpiBinDir)
                setupFile.write(mpiLibDir)
                setupFile.write("source " + mpivarsFilePath)

    else:

        _pm.warn( msg   = "To use ParaMonte in parallel on Mac Operating Systems, \n"
                        + "ParaMonte needs to build the MPICH MPI library on your system. \n"
                        + "To ensure full consistency we recommend building the parallel \n"
                        + "object files of ParaMonte library on your system along with MPICH.\n\n"
                        + "Requesting a full build of ParaMonte library on your system..."
                , marginTop = 1
                , marginBot = 1
                , methodName = _pm.names.paramonte
                )
        build()

####################################################################################################################################

def writeVerificationStatusFile(verified):
    with open(verificationStatusFilePath, "w") as verificationStatusFile:
        verificationStatusFile.write(str(verified))
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

def displayParaMontePythonBanner():
    bannerFilePath = _os.path.join( fileAbsDir, "ParaMontePythonBanner.txt")
    offset = ( len(_pm.version) - 5 ) // 2
    print("")
    with open(bannerFilePath,"r") as file:
        for line in file:
            if "Version" in line:
                line = line.replace(" "*offset+"Version 1.0.0","Version "+_pm.version)
            print(line,end="")
    print("")
    return None

####################################################################################################################################

def build():

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
                        + "    MPI library (MPICH >3.2)\n"
                        + "    OpenCoarrays >2.8\n\n"
                        + "The installation of these software will require 4 to 5 Gb of free space \n"
                        + "on your system (where the ParaMonte-Python interface is already installed).\n"
                        + "Note that the installation script is in Bash and therefore requires Bash shell.\n"
                        + "Also, downloading the prerequisites requires access to the Internet.\n"
                        + "If you have an Internet firewall active on your system, please make sure to \n"
                        + "turn it off before.\n"
                , methodName = _pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        buildEnabled = getUserResponse( msg =   "\n    Do you wish to download and install the ParaMonte library"
                                                "\n    and its prerequisites on your system now (y/n)? " 
                                                )

        currentDir = _os.getcwd()

        pmGitTarPath = _os.path.join( fileAbsDir, "master.tar.gz" )
        download( url = "https://github.com/cdslaborg/paramonte/archive/master.tar.gz"
                , filePath = pmGitTarPath
                )

        pmGitRootDir = _os.path.join( fileAbsDir, "paramonte-master" )

        try:

            import tarfile
            tf = tarfile.open(pmGitTarPath)
            tf.extractall(path=fileAbsDir) # path=pmGitRootDir)

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

            raise CritError(messages.crit_error_bad_command+" "+str(e))
            _pm.abort   ( msg   = "Unzipping of ParaMonte tarball failed.\n"
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
            subprocess.check_call(  [ pmGitInstallScriptPath
                                    , "--lang python"
                                    , "--test_enabled true"
                                    , "--exam_enabled true"
                                    , "--yes-to-all"
                                    ], shell = True )
        except Exception as e:
            raise CritError(messages.crit_error_bad_command+" "+str(e))

        _os.chdir(currentDir)

        # copy files to module folder

        import glob
        import shutil
        pythonBinDir = _os.path.join( pmGitRootDir , "bin" , "Python" , "paramonte" )
        fileList = glob.glob( _os.path.join( pythonBinDir + "libparamonte_*" ) )

        if len(fileList)>0:
            _pm.note( msg   = "ParaMonte kernel libraries build appears to have succeeded. \n"
                            + "copying the kernel files to the paramonte Python module directory..."
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
                shutil.copy(file, fileAbsDir)
            _pm.note( msg   = "ParaMonte kernel libraries should be now usable on your system."
                    , methodName = _pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )
            setupFilePath = _os.path.join( pmGitRootDir , "build", "prerequisites", "prerequisites", "installations", "opencoarrays", "2.8.0", "setup.sh" )
            if _os.path.exists(setupFilePath):
                _pm.warn( msg   = "Whenever you intend to use ParaMonte in the future, before opening your Python session, \n"
                                + "please execute the following command in your Bash shell to ensure all required paths \n"
                                + "are properly defined in your environment:\n\n"
                                + "source " + setupFilePath + " \n\n"
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )
        else:
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

    return None

####################################################################################################################################

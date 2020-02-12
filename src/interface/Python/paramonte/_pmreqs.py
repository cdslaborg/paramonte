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

pmReleaseVersion = "1.0.0"
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

        if (_pm.arch=="x64") and (isWin32 or isLinux): # or isMacOS):

            displayParaMontePythonBanner()

            # library path

            if not isWin32: setupUnixPath()

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
                    + "on AMD64 (64-bit) architecture for Windows/Linux Operating Systems (OS). \n"
                    + "Your system appears to be of a different architecture or OS. As a result, \n"
                    + "the core sampler routines of ParaMonte will not be available on your system. However, \n"
                    + "the generic Python interface of ParaMonte will is available on your system, which can \n"
                    + "be used for post-processing and visualization of the output files from \n"
                    + "already-performed ParaMonte simulations or other similar Monte Carlo simulations. \n"
                    + "There are ongoing efforts, right now as you read this message, to further increase the \n"
                    + "availability of ParaMonte library on a wider-variety of platforms and architectures. \n"
                    + "If you are a Mac user, please check back by April 1, 2020 for a new release of ParaMonte \n"
                    + "with support for macOS. \n"
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

def getBashrcContents():
    bashrcPath = _os.path.expanduser("~/.bashrc")
    if _os.path.isfile(bashrcPath):
        with open(bashrcPath,"r") as bashrcFile:
            bashFileContents = bashrcFile.read()
    else:
        bashFileContents = ""
        with open(bashrcPath,"w") as bashrcFile:
            pass
    return bashFileContents

####################################################################################################################################

def setupUnixPath():

    bashFileContents = getBashrcContents()
    dlibcmd = "export LD_LIBRARY_PATH=" + fileAbsDir + ":$LD_LIBRARY_PATH"
    if dlibcmd not in bashFileContents:
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
            if (pathcmd not in bashFileContents) or (dlibcmd not in bashFileContents):
                _os.system( "chmod 777 ~/.bashrc")
                _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                _os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte local GNU installation setup >>>' >> ~/.bashrc" )
                if pathcmd is not None:
                    if pathcmd not in bashFileContents:
                        _os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${PATH+x} ]; then' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '    export PATH=.' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '" + pathcmd + "' >>  ~/.bashrc" )
                if dlibcmd is not None:
                    if dlibcmd not in bashFileContents:
                        _os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc" )
            if pathcmd not in bashFileContents or dlibcmd not in bashFileContents:
                _os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte local GNU installation setup <<<' >> ~/.bashrc" )
                _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                _os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

        pathcmd = None
        dlibcmd = None
        if localInstallDir.mpi.bin is not None: pathcmd = "export PATH=" + localInstallDir.mpi.bin + ":$PATH"
        if localInstallDir.mpi.lib is not None: dlibcmd = "export LD_LIBRARY_PATH=" + localInstallDir.mpi.lib + ":$LD_LIBRARY_PATH"
        if (pathcmd is not None) or (dlibcmd is not None):
            if (pathcmd not in bashFileContents) or (dlibcmd not in bashFileContents):
                _os.system( "chmod 777 ~/.bashrc")
                _os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                _os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte local MPI installation setup >>>' >> ~/.bashrc" )
                if pathcmd is not None:
                    if pathcmd not in bashFileContents:
                        _os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${PATH+x} ]; then' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '    export PATH=.' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '" + pathcmd + "' >>  ~/.bashrc" )
                if dlibcmd is not None:
                    if dlibcmd not in bashFileContents:
                        _os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        _os.system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc" )
            if pathcmd not in bashFileContents or dlibcmd not in bashFileContents:
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


    pmGitRootDir = _os.path.join( fileAbsDir , "paramonte-master" )

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
                        setupFile.write("@echo off\n")
                        setupFile.write("cd " + path + " && mpivars.bat quiet\n")
                        setupFile.write("cd " + fileAbsDir + "\n")
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
            downloadList = ["impi.dll","impi.pdb"]
            for mpiFileNameExt in ["impi.dll","impi.pdb","libfabric.dll"]:
                mpiFilePath = _os.path.join( fileAbsDir, mpiFileNameExt )
                download( url = "https://github.com/cdslaborg/paramonte/releases/download/" + pmReleaseVersion + "/" + mpiFileNameExt
                        , filePath = mpiFilePath
                        )
            
        if isLinux:
            mpiFileExt = ".tgz"
            mpiVersion = "2018.2.199"
            mpiFileName = "l_mpi-rt_" + mpiVersion

        mpiFileNameExt = mpiFileName + mpiFileExt
        mpiFilePath = _os.path.join( fileAbsDir, mpiFileNameExt )
        download( url = "https://github.com/cdslaborg/paramonte/releases/download/" + pmReleaseVersion + "/" + mpiFileNameExt
                , filePath = mpiFilePath
                )

        _pm.note( msg = "Installing the Intel MPI library for 64-bit architecture...\n"
                      + "file location: " + mpiFilePath
                , methodName = _pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        _pm.warn( msg = "Please do not change the default installation location MPI suggested by the installer.\n"
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
                tf.extractall(path=fileAbsDir)
                mpiExtractDir = _os.path.join( fileAbsDir, mpiFileName )

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

            from pathlib import Path
            home = str(Path.home())
            setupFilePath = _os.path.join( fileAbsDir, "setup.sh" )

            mpiRootDirNotFound = True
            installationRootDirList = [ "/opt", home ]
            while mpiRootDirNotFound:

                for installationRootDir in installationRootDirList:
                    mpiTrunkDir = _os.path.join( "intel", "compilers_and_libraries_" + mpiVersion, "linux", "mpi", "intel64" )
                    mpiRootDir = _os.path.join( installationRootDir, mpiTrunkDir )
                    if _os.path.isdir(mpiRootDir):
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
                                + "This is highly unusual. Normally, Intel MPI libraries must be installed\n"
                                + "in the following directory:\n\n"
                                + "    " + mpiRootDir + "\n\n"
                                + "If you cannot manually find the Intel MPI installation directory,\n"
                                + "it is likely that the installation might have somehow failed.\n"
                                + "If you do find the installation directory, try to locate the\n"
                                + "'mpivars.sh' file which is normally installed in the following path:\n\n"
                                + "    " + mpivarsFilePath + "\n\n"
                                + "Before attempting to run any parallel ParaMonte simulation, \n"
                                + "make sure you source this file, like the following:\n\n"
                                + "    source " + mpivarsFilePath + "\n\n"
                                + "where you will have to replace the path in the above with the \n"
                                + "correct path that you find on your system."
                        , methodName = _pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

            else:

                mpiBinDir = _os.path.join( mpiRootDir, "bin" )
                mpiLibDir = _os.path.join( mpiRootDir, "lib" )
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

                    bashFileContents = getBashrcContents()
                    mpivarsFileCommand = "source " + mpivarsFilePath
                    if mpivarsFileCommand not in bashFileContents:
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

        _pm.warn( msg   = "To use ParaMonte in parallel on Mac Operating Systems, \n"
                        + "ParaMonte needs to build the MPICH MPI library on your system. \n"
                        + "To ensure full consistency, we recommend building the parallel \n"
                        + "object files of ParaMonte library on your system along with MPICH.\n\n"
                        + "Requesting a full build of ParaMonte library on your system..."
                , marginTop = 1
                , marginBot = 1
                , methodName = _pm.names.paramonte
                )
        build()

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
                        + "    MPI library (MPICH >3.2)\n"
                        + "    OpenCoarrays >2.8\n\n"
                        + "The installation of these software will require 4 to 5 Gb of free space \n"
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

                print(str(e))
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

                    bashFileContents = getBashrcContents()
                    setupFilePathCmd = "source " + setupFilePath
                    if setupFilePathCmd not in bashFileContents:
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

            writeVerificationStatusFile(verificationEnabled=True)

        else:

            _pm.warn( msg   = "Aborting the ParaMonte-for-Python local build on your system."
                    , methodName = _pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )

    return None

####################################################################################################################################

####################################################################################################################################
####################################################################################################################################
####
####   MIT License
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   Permission is hereby granted, free of charge, to any person obtaining a
####   copy of this software and associated documentation files (the "Software"),
####   to deal in the Software without restriction, including without limitation
####   the rights to use, copy, modify, merge, publish, distribute, sublicense,
####   and/or sell copies of the Software, and to permit persons to whom the
####   Software is furnished to do so, subject to the following conditions:
####
####   The above copyright notice and this permission notice shall be
####   included in all copies or substantial portions of the Software.
####
####   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
####   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
####   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
####   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
####   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
####   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
####   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
####
####   ACKNOWLEDGMENT
####
####   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
####   As per the ParaMonte library license agreement terms, if you use any parts of
####   this library for any purposes, kindly acknowledge the use of ParaMonte in your
####   work (education/research/industry/development/...) by citing the ParaMonte
####   library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

import os
import sys
import copy
import typing as tp
import numpy as np
import _paramonte as pm
import warnings

Struct = pm.Struct
newline = pm.newline

verificationStatusFilePath = os.path.join( pm.path.auxil, ".verificationEnabled" )

####################################################################################################################################

buildInstructionNoteWindows = ""
buildInstructionNoteUnix    = ("If your platform is non-Windows and is compatible with the " + newline
                            + "GNU Compiler Collection (GCC), you can also build the required " + newline
                            + "ParaMonte kernel's shared object files on your system by calling " + newline
                            + "the ParaMonte's build() method from within your Python session, " + newline
                            + newline
                            + "    import paramonte as pm" + newline
                            + "    pm.build()"
                            )

####################################################################################################################################
#### verify
####################################################################################################################################

def verify(reset = True):
    """

    checks (or rechecks) the requirements of the installed ParaMonte library.

        **Usage**

            .. code-block:: python

                import paramonte as pm
                pm.verify()

        **Parameters**

            reset

                A boolean whose default value is ``True``. If ``True``,
                a thorough verification of the existence of the required
                libraries will performed, as if it is the first ParaMonte
                module import.

        **Returns**

            None

    """

    if not isinstance(reset,bool):
        raise Exception ( newline
                        + "The input argument reset must be a logical (boolean) value: True or False"
                        + newline
                        )

    #### require Python >3.6 for type hints

    _MIN_PYTHON = (3,6)
    if sys.version_info < _MIN_PYTHON:
        sys.exit("Python %s.%s or newer is required for ParaMonte. Please install the latest version of Python.\n" % _MIN_PYTHON)

    #### ensure messages are printed only for the first-time import

    if reset: writeVerificationStatusFile("True")

    with open(verificationStatusFilePath, "r") as verificationStatusFile:
        verificationEnabledString = verificationStatusFile.readline()

    if verificationEnabledString=="False":
        verificationEnabled = False
    elif verificationEnabledString=="True" or verificationEnabledString=="Testing":
        verificationEnabled = True
    else:
        raise Exception ( newline
                        + "The internal settings of the ParaMonte library appears to have been messed up" + newline
                        + "potentially by the user, the operating system, Python, or other applications." + newline
                        + "Please reinstall ParaMonte by typing the following commands" + newline
                        + "on a Python-aware command-line interface:" + newline
                        + newline
                        + "    pip uninstall paramonte" + newline
                        + "    pip install --user --upgrade paramonte" + newline
                        + newline
                        )

    if verificationEnabled:

        #### ensure 64-bit architecture

        if (pm.platform.arch=="x64") and (pm.platform.isWin32 or pm.platform.isLinux or pm.platform.isMacOS):

            displayParaMonteBanner()

            #### display dependency version message

            displayDependencyVersionMessage()

            #### verify module dependencies

            #### On some systems like TACC, the matplotlib causes segmentation fault that is not controllable in any way.
            #### This is apparently a bug in the older versions of matplotlib. Until it is fully resolved, the following
            #### dependency version check is commented out.

            # verifyDependencyVersion()

            #### library path

            if not pm.platform.isWin32: setupUnixPath()

            #### search for the MPI library

            mpi = findMPI()

            if mpi.install.found and not mpi.path.broken:

                writeVerificationStatusFile("False")

            else:

                if mpi.install.found and mpi.path.broken:
                    msg=( "An MPI library installation appears to exist on your system, however, " + newline
                        + "some components of the library appear to be missing, or the environmental " + newline
                        + "path to the MPI library installation is corrupted. You can inspect the " + newline
                        + "contents of the environmental path variable for potential path " + newline
                        + "corruptions by typing, " + newline
                        + newline
                        + "    import os" + newline
                        + "    os.environ[\"PATH\"]" + newline
                        + newline
                        + "on your Python command line. If you or the ParaMonte library (on your behalf) " + newline
                        + "have already successfully installed an MPI library on your system, " + newline
                        + "you can safely ignore this warning and avoid further reinstallation of the " + newline
                        + "MPI library. Otherwise, you can continue to reinstall the MPI library."
                        )
                else:
                    msg=( "The MPI runtime libraries for 64-bit architecture could not be detected " + newline
                        + "on your system. The MPI runtime libraries are required for the parallel " + newline
                        + "ParaMonte simulations."
                        )

                pm.note ( msg   = msg + newline
                                + "For Windows and Linux operating systems, you can download and install the " + newline
                                + "Intel MPI runtime libraries, free of charge, from the Intel website, " + newline
                                + newline
                                + "    " + pm.website.intel.mpi.home.url + newline
                                + newline
                                + "For macOS (Darwin), you can download and install the Open-MPI library. " + newline
                                + newline
                                + "    " + pm.website.openmpi.home.url + newline
                                + newline
                                + "Alternatively, the ParaMonte library can automatically install these " + newline
                                + "libraries for you now. If you don't know how to download and install the " + newline
                                + "correct MPI runtime library version, we strongly recommend that you let the " + newline
                                + "ParaMonte library to install this library for you. If so, ParaMonte will need " + newline
                                + "access to the World-Wide-Web to download the library and will need your " + newline
                                + "administrative permission to install it on your system. Therefore, if " + newline
                                + "you have any active firewall on your system such as ZoneAlarm, please " + newline
                                + "make sure your firewall allows ParaMonte to access the Internet."
                        , marginTop = 1
                        , marginBot = 1
                        , methodName = pm.names.paramonte
                        )

                if verificationEnabledString=="Testing":

                    writeVerificationStatusFile("False")

                else:

                    isYes = getUserResponse( msg =  "\n    Do you wish to download and install the MPI runtime library"
                                                    "\n    (only needed for parallel simulations) on your system now (y/n)? "
                                            )

                    if isYes:
                        installMPI()
                        writeVerificationStatusFile("Testing")
                    else:
                        pm.note ( msg   = "Skipping the MPI library installation... " + newline
                                        + "It is now the user's responsibility to ensure an MPI runtime library" + newline
                                        + "exists on the system for parallel simulations. " + newline
                                       #+ "If you ever wish to install the MPI libraries via ParaMonte again, " + newline
                                       #+ "you can try: " + newline
                                       #+ newline
                                       #+ "    import paramonte as pm" + newline
                                       #+ "    pm.verify()" + newline
                                       #+ newline
                                        + "For more information visit:" + newline
                                        + newline
                                        + "    " + pm.website.home.url
                        , marginTop = 1
                        , marginBot = 1
                        , methodName = pm.names.paramonte
                        )
                        writeVerificationStatusFile("False")

            dispFinalMessage()

        else:

            warnForUnsupportedPlatform()
            build()

    return None

####################################################################################################################################
#### getUserResponse
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
#### download
####################################################################################################################################

def download(url,filePath):
    import urllib.request
    import shutil
    with urllib.request.urlopen(url) as response, open(filePath, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
    return None

####################################################################################################################################
#### warnForUnsupportedPlatform
####################################################################################################################################

def warnForUnsupportedPlatform():
    pm.warn ( msg   = "The ParaMonte sampler kernel is currently exclusively available on " + newline
                    + "AMD64 (64-bit) architecture for Windows/Linux/Darwin Operating Systems (OS). " + newline
                    + "Your system appears to be of a different architecture or OS. As a result, " + newline
                    + "the core sampler routines of ParaMonte will not be available on your system. " + newline
                    + "However, the generic Python interface of ParaMonte will be available on your " + newline
                    + "system, which can be used for post-processing and visualization of the output " + newline
                    + "files from already-performed ParaMonte simulations or other similar Monte " + newline
                    + "Carlo simulations. There are ongoing efforts, right now as you read this " + newline
                    + "message, to further increase the availability of ParaMonte library on a " + newline
                    + "wider-variety of platforms and architectures. Stay tuned for updates by " + newline
                    + "visiting, " + newline
                    + newline
                    + "    " + pm.website.home.url + newline
                    + newline
                    + "That said, " + newline
                    + newline
                    + "if your platform is non-Windows and is compatible with the GNU Compiler " + newline
                    + "Collection (GCC), you can also build the required ParaMonte kernel's " + newline
                    + "shared object files on your system by calling ParaMonte module's " + newline
                    + "build() function from within your Python environment."
            , marginTop = 1
            , marginBot = 1
            , methodName = pm.names.paramonte
            )
    return None

####################################################################################################################################
#### getBashrcContents
####################################################################################################################################

def getBashrcContents():
    bashrcPath = os.path.expanduser("~/.bashrc")
    if os.path.isfile(bashrcPath):
        with open(bashrcPath,"r") as bashrcFile:
            bashrcContents = bashrcFile.read()
    else:
        bashrcContents = ""
        with open(bashrcPath,"w") as bashrcFile:
            pass
    return bashrcContents

####################################################################################################################################
#### getBashProfileContents
####################################################################################################################################

def getBashProfileContents():
    bashProfilePath = os.path.expanduser("~/.bash_profile")
    bashProfileFileExists = os.path.isfile(bashProfilePath)
    bashProfileContents = ""
    if bashProfileFileExists:
        with open(bashProfilePath,"r") as bashProfileFile:
            bashProfileContents = bashProfileFile.read()
    if ".bashrc" not in bashProfileContents:
        with open(bashProfilePath,"a+") as bashProfileFile:
            bashProfileFile.write("\n[ -f $HOME/.bashrc ] && . $HOME/.bashrc\n")
    return bashProfileContents

####################################################################################################################################
#### setupUnixPath
####################################################################################################################################

def setupUnixPath():

    bashrcContents = getBashrcContents()
    dlibcmd = "export LD_LIBRARY_PATH=" + pm.path.lib + ":$LD_LIBRARY_PATH"
    if dlibcmd not in bashrcContents:
        os.system( "chmod 777 ~/.bashrc")
        os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
        os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte shared library setup >>>' >> ~/.bashrc" )
        os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc" )
        os.system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc" )
        os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
        os.system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc" )
        os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte shared library setup <<<' >> ~/.bashrc" )
        os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
        os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

    localInstallDir = getLocalInstallDir()
    if localInstallDir.root is not None:

        pathcmd = None
        dlibcmd = None
        if localInstallDir.gnu.bin is not None: pathcmd = "export PATH=" + localInstallDir.gnu.bin + ":$PATH"
        if localInstallDir.gnu.lib is not None: dlibcmd = "export LD_LIBRARY_PATH=" + localInstallDir.gnu.lib + ":$LD_LIBRARY_PATH"
        if (pathcmd is not None) or (dlibcmd is not None):
            if (pathcmd not in bashrcContents) or (dlibcmd not in bashrcContents):
                os.system( "chmod 777 ~/.bashrc")
                os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte local GNU installation setup >>>' >> ~/.bashrc" )
                if pathcmd is not None:
                    if pathcmd not in bashrcContents:
                        os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${PATH+x} ]; then' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '    export PATH=.' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '" + pathcmd + "' >>  ~/.bashrc" )
                if dlibcmd is not None:
                    if dlibcmd not in bashrcContents:
                        os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc" )
            if pathcmd not in bashrcContents or dlibcmd not in bashrcContents:
                os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte local GNU installation setup <<<' >> ~/.bashrc" )
                os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

        pathcmd = None
        dlibcmd = None
        if localInstallDir.mpi.bin is not None: pathcmd = "export PATH=" + localInstallDir.mpi.bin + ":$PATH"
        if localInstallDir.mpi.lib is not None: dlibcmd = "export LD_LIBRARY_PATH=" + localInstallDir.mpi.lib + ":$LD_LIBRARY_PATH"
        if (pathcmd is not None) or (dlibcmd is not None):
            if (pathcmd not in bashrcContents) or (dlibcmd not in bashrcContents):
                os.system( "chmod 777 ~/.bashrc")
                os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte local MPI installation setup >>>' >> ~/.bashrc" )
                if pathcmd is not None:
                    if pathcmd not in bashrcContents:
                        os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${PATH+x} ]; then' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '    export PATH=.' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '" + pathcmd + "' >>  ~/.bashrc" )
                if dlibcmd is not None:
                    if dlibcmd not in bashrcContents:
                        os.system( "chmod 777 ~/.bashrc && echo 'if [ -z ${LD_LIBRARY_PATH+x} ]; then' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '    export LD_LIBRARY_PATH=.' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo 'fi' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '" + dlibcmd + "' >>  ~/.bashrc" )
            if pathcmd not in bashrcContents or dlibcmd not in bashrcContents:
                os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte local MPI installation setup <<<' >> ~/.bashrc" )
                os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

    return None

####################################################################################################################################
#### getLocalInstallDir
####################################################################################################################################

def getLocalInstallDir():

    localInstallDir = Struct()
    localInstallDir.root = None

    localInstallDir.mpi = Struct()
    localInstallDir.mpi.root = None
    localInstallDir.mpi.bin = None
    localInstallDir.mpi.lib = None

    localInstallDir.gnu = Struct()
    localInstallDir.gnu.root = None
    localInstallDir.gnu.bin = None
    localInstallDir.gnu.lib = None

    localInstallDir.caf = Struct()
    localInstallDir.caf.root = None
    localInstallDir.caf.bin = None
    localInstallDir.caf.lib = None


    pmGitRootDir = os.path.join( pm.path.root , "paramonte-master" )

    if os.path.isdir(pmGitRootDir):

        localInstallDir.root = pmGitRootDir

        # mpi

        _ = os.path.join( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "mpich", "3.2" )
        if os.path.isdir(_):
            localInstallDir.mpi.root = _
            _ = os.path.join( localInstallDir.mpi.root, "bin" )
            if os.path.isdir(_): localInstallDir.mpi.bin = _
            _ = os.path.join( localInstallDir.mpi.root, "lib" )
            if os.path.isdir(_): localInstallDir.mpi.lib = _

        # gnu

        _ = os.path.join( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "gnu", "8.3.0" )
        if os.path.isdir(_):
            localInstallDir.gnu.root = _
            _ = os.path.join( localInstallDir.gnu.root, "bin" )
            if os.path.isdir(_): localInstallDir.gnu.bin = _
            _ = os.path.join( localInstallDir.gnu.root, "lib64" )
            if os.path.isdir(_): localInstallDir.gnu.lib = _

        # caf

        _ = os.path.join( localInstallDir.root, "build", "prerequisites", "prerequisites", "installations", "opencoarrays", "2.8.0" )
        if os.path.isdir(_):
            localInstallDir.caf.root = _
            _ = os.path.join( localInstallDir.caf.root, "bin" )
            if os.path.isdir(_): localInstallDir.caf.bin = _
            _ = os.path.join( localInstallDir.caf.root, "lib64" )
            if os.path.isdir(_): localInstallDir.caf.lib = _

    return localInstallDir

####################################################################################################################################
#### findMPI
####################################################################################################################################

def findMPI():
    """

    Return a structure containing the paths to
    different components of the MPI library.

    """

    mpi = Struct()
    mpi.path = Struct()
    mpi.install = Struct()
    mpi.install.bin = Struct()
    mpi.install.bin.mpiexec = Struct()
    mpi.install.bin.mpivars = Struct()

    mpi.path.broken = False
    mpi.install.found = False
    mpi.install.bin.found = False
    mpi.install.bin.path = None
    mpi.install.bin.mpiexec.found = False
    mpi.install.bin.mpiexec.path = None
    mpi.install.bin.mpivars.found = False
    mpi.install.bin.mpivars.path = None

    if pm.platform.isWin32:

        pathList = os.environ['PATH'].split(";")
        for thisPath in pathList:

            pathLower = thisPath.lower().replace("\\","")
            if ("mpiintel64bin" in pathLower):

                mpi.install.bin.found = True
                mpi.install.bin.path = thisPath

                mpiexecFilePath = os.path.join(mpi.install.bin.path,"mpiexec.exe")
                if os.path.isfile(mpiexecFilePath):
                    mpi.install.bin.mpiexec.found = True
                    mpi.install.bin.mpiexec.path = mpiexecFilePath

                mpivarsFilePath = os.path.join( thisPath, "mpivars.bat" )
                if os.path.isfile(mpivarsFilePath):

                    mpi.install.bin.mpivars.found = True
                    mpi.install.bin.mpivars.path = mpivarsFilePath

                    mpivarsCommand = '"' + mpivarsFilePath + '"'
                    pm.note ( msg   = "Intel MPI library for 64-bit architecture detected at: " + newline
                                    + newline
                                    + "    " + thisPath + newline
                                    + newline
                                    + "To perform ParaMonte simulations in parallel on a single node, " + newline
                                    + "run the following two commands, in the form and order specified, " + newline
                                    + "on a Python-aware mpiexec-aware command-line interface such as " + newline
                                    + "Anaconda3 Windows command prompt: " + newline
                                    + newline
                                    + "    " + mpivarsCommand + newline
                                    + newline
                                    + "    mpiexec -localonly -n NUM_PROCESSES python main.py" + newline
                                    + newline
                                    + "where, " + newline
                                    + newline
                                    + "    0.   the first command defines the essential environment variables, " + newline
                                    + "         and the second command runs in the simulation in parallel, where, " + newline
                                    + "    1.   you should replace NUM_PROCESSES with the number of processors " + newline
                                    + "         you wish to assign to your simulation task and, " + newline
                                    + "    2.   the flag '-localonly' indicates a parallel simulation on only " + newline
                                    + "         a single node (this flag will obviate the need for the MPI " + newline
                                    + "         library credentials registration). For more information, visit: " + newline
                                    + "         " + pm.website.intel.mpi.windows.url + " " + newline
                                    + "    3.   main.py is the Python file which serves as the entry point to " + newline
                                    + "         your simulation, where you call the ParaMonte sampler routines. " + newline
                                    + newline
                                    + "Note that the above two commands must be executed on a command-line that " + newline
                                    + "recognizes both Python and mpiexec applications, such as the Anaconda " + newline
                                    + "command-line interface. For more information, in particular, on how " + newline
                                    + "to register to run Hydra services for multi-node simulations " + newline
                                    + "on Windows servers, visit: " + newline
                                    + newline
                                    + "    " + pm.website.home.url
                            , marginTop = 1
                            , marginBot = 1
                            , methodName = pm.names.paramonte
                            )

                    setupFilePath = os.path.join( pm.path.auxil, "setup.bat" )
                    with open(setupFilePath, "w") as setupFile:
                        setupFile.write("@echo off\n")
                        setupFile.write("cd " + thisPath + " && mpivars.bat quiet\n")
                        setupFile.write("cd " + pm.path.root + "\n")
                        setupFile.write("@echo on\n")

                    mpi.install.found = mpi.install.bin.found and mpi.install.bin.mpiexec.found and mpi.install.bin.mpivars.found
                    if mpi.install.found: break

    elif pm.platform.isLinux:

        pathList = os.environ['PATH'].split(":")
        for thisPath in pathList:

            pathLower = thisPath.lower().replace("/","")
            if ("linuxmpiintel64" in pathLower):

                mpi.install.bin.found = True
                mpi.install.bin.path = thisPath

                mpiexecFilePath = os.path.join( mpi.install.bin.path, "mpiexec" )
                if os.path.isfile(mpiexecFilePath):
                    mpi.install.bin.mpiexec.found = True
                    mpi.install.bin.mpiexec.path = mpiexecFilePath

                mpivarsFilePath = os.path.join( thisPath, "mpivars.sh" )
                if os.path.exists(mpivarsFilePath):

                    mpi.install.bin.mpivars.found = True
                    mpi.install.bin.mpivars.path = mpivarsFilePath

                    mpivarsCommand = '"' + mpivarsFilePath + '"'
                    pm.note ( msg   = "Intel MPI library for 64-bit architecture detected at: " + newline
                                    + newline
                                    + "    " + thisPath + newline
                                    + newline
                                    + "To perform ParaMonte simulations in parallel on a single node, " + newline
                                    + "run the following two commands, in the form and order specified, " + newline
                                    + "in a Bash shell (terminal), " + newline
                                    + newline
                                    + "    source " + mpivarsCommand + newline
                                    + newline
                                    + "    mpiexec -n NUM_PROCESSES python main.py" + newline
                                    + newline
                                    + "where, " + newline
                                    + newline
                                    + "    0.   the first command defines the essential environment variables" + newline
                                    + "         and the second command runs in the simulation in parallel, where," + newline
                                    + "    1.   you should replace NUM_PROCESSES with the number of processors " + newline
                                    + "         you wish to assign to your simulation task, " + newline
                                    + "    2.   main.py is the Python file which serves as the entry point to " + newline
                                    + "         your simulation, where you call ParaMonte sampler routines. " + newline
                                    + newline
                                    + "For more information on how to install and use and run parallel " + newline
                                    + "ParaMonte simulations on Linux systems, visit: " + newline
                                    + newline
                                    + pm.website.home.url
                            , marginTop = 1
                            , marginBot = 1
                            , methodName = pm.names.paramonte
                            )

                    try:
                        setupFilePath = os.path.join( pm.path.auxil, "setup.sh" )
                        with open(setupFilePath, "w") as setupFile:
                            setupFile.write("source " + mpivarsCommand)
                    except:
                        pm.warn ( msg   = "Failed to create the MPI setup file. " + newline
                                        + "It looks like the ParaMonte library directory is read-only. " + newline
                                        + "This can be potentially problematic. Skipping for now..."
                                , marginTop = 1
                                , marginBot = 1
                                , methodName = pm.names.paramonte
                                )

                    mpi.install.found = mpi.install.bin.found and mpi.install.bin.mpiexec.found and mpi.install.bin.mpivars.found
                    if mpi.install.found: break

    elif pm.platform.isMacOS:

        import shutil

        gfortranPath = None
        try:
            import subprocess
            gfortranVersion = subprocess.run(args=["gfortran", "--version"],capture_output=True)
            if "GCC 10." in str(gfortranVersion.stdout): gfortranPath = shutil.which("gfortran")
        except:
            warnings.warn("Failed to capture the GNU Compiler Collection version...")

        mpi.install.bin.mpiexec.path = None
        try:
            import subprocess
            mpiexecVersion = subprocess.run(args=["mpiexec", "--version"],capture_output=True)
            if "open-mpi" in str(mpiexecVersion.stdout): mpi.install.bin.mpiexec.path = shutil.which("mpiexec")
        except:
            warnings.warn("Failed to capture the mpiexec version...")

        if (mpi.install.bin.mpiexec.path is not None) and (gfortranPath is not None):
            mpi.install.bin.found = True
            mpi.install.bin.mpiexec.found = True
            mpi.install.bin.mpivars.found = True # dummy
            mpi.install.bin.path = os.path.dirname(mpi.install.bin.mpiexec.path)
            pm.note ( msg   = "MPI runtime libraries detected at: " + newline
                            + newline
                            + "    " + mpi.install.bin.path + newline
                            + newline
                            + "To perform ParaMonte simulations in parallel on a single node, " + newline
                            + "run the following command, in the form and order specified, " + newline
                            + "in a Bash shell (terminal), " + newline
                            + newline
                            + "    mpiexec -n NUM_PROCESSES python main.py" + newline
                            + newline
                            + "where, " + newline
                            + newline
                            + "    0.   the first command defines the essential environment variables " + newline
                            + "         and the second command runs in the simulation in parallel, where, " + newline
                            + "    1.   you should replace NUM_PROCESSES with the number of processors " + newline
                            + "         you wish to assign to your simulation task, " + newline
                            + "    2.   main.py is the Python file which serves as the entry point to " + newline
                            + "         your simulation, where you call the ParaMonte sampler routines. " + newline
                            + newline
                            + "For more information on how to install and use and run parallel ParaMonte " + newline
                            + "simulations on the macOS (Darwin) operating systems, visit:" + newline
                            + newline
                            + pm.website.home.url
            , marginTop = 1
            , marginBot = 1
            , methodName = pm.names.paramonte
            )

        mpi.install.found = mpi.install.bin.found and mpi.install.bin.mpiexec.found and mpi.install.bin.mpivars.found

    else:

        LocalInstallDir = getLocalInstallDir()
        if (LocalInstallDir.mpi.bin is not None) and (LocalInstallDir.mpi.lib is not None):

            mpi.install.bin.found = True
            mpi.install.bin.path = LocalInstallDir.mpi.bin

            mpiexecFilePath = os.path.join(mpi.install.bin.path,"mpiexec")
            if os.path.isfile(mpiexecFilePath):
                mpi.install.bin.mpiexec.found = True
                mpi.install.bin.mpiexec.path = mpiexecFilePath

        mpi.install.bin.mpivars.found = mpi.install.bin.found and mpi.install.bin.mpiexec.found # dummy
        mpi.install.found = mpi.install.bin.found and mpi.install.bin.mpiexec.found and mpi.install.bin.mpivars.found

    #### one last try to find the MPI library if not found yet

    if not mpi.install.found:

        mpi.path.broken = True

        if pm.platform.isLinux:

            defaultIntelLinuxMpiPath = getDefaultIntelLinuxMpiPath()
            if defaultIntelLinuxMpiPath.mpiRootDirNotFound:
                return mpi
            else:
                mpi.install.found = True
                pm.warn ( msg   = "The PATH environmental variable of your Bash terminal does not point to " + newline
                                + "any current installation of the Intel MPI runtime libraries on your system, " + newline
                                + "however, ParaMonte has detected a hidden installation of the Intel MPI " + newline
                                + "runtime libraries on your system at, " + newline
                                + newline
                                + "    " + defaultIntelLinuxMpiPath.mpiDefaultRootDirList[-1] + newline
                                + newline
                                + "Include this path to your terminal's PATH environmental variable to ensure " + newline
                                + "the MPI runtime libraries will be properly detected in the future."
                        , marginTop = 1
                        , marginBot = 1
                        , methodName = pm.names.paramonte
                        )
                # mpi.install.bin.path = setupIntelLinuxMpiPath(defaultIntelLinuxMpiPath)

        elif pm.platform.isWin32:

            pm.warn ( msg   = "Failed to detect the Intel MPI library for 64-bit architecture." + newline
                            + "Now searching through the installed applications..." + newline
                            + "This may take some time..."
                    , marginTop = 1
                    , marginBot = 1
                    , methodName = pm.names.paramonte
                    )

            import subprocess
            installedApp = str(subprocess.run(args=["wmic","product","get","Name,","Version"],capture_output=True).stdout)

            if "Intel MPI" in installedApp:
                mpi.install.found = True
                pm.note ( msg = "Possible Intel MPI installation detected:"
                        , marginTop = 0
                        , marginBot = 1
                        , methodName = pm.names.paramonte
                        )
                installedAppList = str(installedApp).replace("\\r","").split("\\n")
                for app in installedAppList:
                    appClean = app.replace(chr(13),"").replace(chr(10),"") # remove cr, nl
                    if "Intel MPI" in appClean:
                        pm.note ( msg = appClean
                                , marginTop = 0
                                , marginBot = 0
                                , methodName = pm.names.paramonte
                                )
    return mpi

####################################################################################################################################
#### getPrereqs
####################################################################################################################################

def getPrereqs(DependencyList = None):

    prereqs = Struct()
    prereqs.mpi = Struct()
    prereqs.mpi.intel = Struct()

    prereqs.list = getDependencyList() if DependencyList is None else DependencyList

    if pm.platform.isLinux:
        intelMpiFilePrefix, intelMpiFileSuffix = "l_mpi-rt_" , ".tgz"
    elif pm.platform.isWin32:
        intelMpiFilePrefix, intelMpiFileSuffix = "w_mpi-rt_p_" , ".exe"
    else:
        return prereqs

    for dependency in prereqs.list:
        fullFilePath = os.path.join( pm.path.lib, dependency )
        if intelMpiFilePrefix in dependency and intelMpiFileSuffix in dependency:
            prereqs.mpi.intel.fullFileName = dependency
            prereqs.mpi.intel.fullFilePath = fullFilePath
            prereqs.mpi.intel.fileName = prereqs.mpi.intel.fullFileName.split(intelMpiFileSuffix)[0]
            prereqs.mpi.intel.version = prereqs.mpi.intel.fileName.split(intelMpiFilePrefix)[1]

    return prereqs

####################################################################################################################################
#### getDefaultIntelLinuxMpiPath
####################################################################################################################################

def getDefaultIntelLinuxMpiPath(prereqs = None):
    if prereqs is None: prereqs = getPrereqs()
    mpiPath = Struct()
    mpiPath.mpiDefaultRootDirList = []
    mpiPath.mpiRootDirNotFound = True
    mpiPath.mpivarsDefaultFilePathList = []
    mpiPath.installationRootDirList = [ "/opt", pm.path.home ]
    mpiPath.mpiTrunkDir = os.path.join("intel", "compilers_and_libraries_" + prereqs.mpi.intel.version, "linux", "mpi", "intel64")
    for installationRootDir in mpiPath.installationRootDirList:
        mpiPath.mpiDefaultRootDirList.append( os.path.join(installationRootDir, mpiPath.mpiTrunkDir) )
        mpiPath.mpivarsDefaultFilePathList.append( os.path.join(mpiPath.mpiDefaultRootDirList[-1],"bin","mpivars.sh") )
        if os.path.isdir(mpiPath.mpiDefaultRootDirList[-1]):
            mpiPath.mpiRootDirNotFound = False
            break
    return mpiPath

####################################################################################################################################
#### getDependencyList
####################################################################################################################################

def getDependencyList():
    fileName = ".dependencies_";
    if pm.platform.isWin32: fileName = fileName + "windows"
    if pm.platform.isMacOS: fileName = fileName + "macos"
    if pm.platform.isLinux: fileName = fileName + "linux"
    with open(os.path.join(pm.path.auxil, fileName), "r") as depFile: lines = depFile.read().splitlines()
    dependencyList = []
    for count,item in enumerate(lines):
        if item[0]!="!": dependencyList.append(item) # remove comment lines
    return dependencyList

####################################################################################################################################
#### installMPI
####################################################################################################################################

def installMPI():

    if pm.platform.isWin32 or pm.platform.isLinux:

        pm.note ( msg = "Downloading the ParaMonte parallel library prerequisites... " + newline
                      + "Please make sure your firewall allows access to the Internet. "
                , methodName = pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        prereqs = getPrereqs()

        for dependency in prereqs.list:

            fullFilePath = os.path.join(pm.path.lib, dependency)
            thisVersion = pm.version.kernel.dump()
            while thisVersion is not None:
                try:
                    download( url = pm.website.github.release.url + "/download/" + thisVersion + "/" + dependency
                            , filePath = fullFilePath
                            )
                    break
                except:
                    thisVersion = getPreviousVersion(thisVersion)

            if thisVersion is None:
                pm.warn ( msg   = "Exhausted all releases of the ParaMonte library in search " + newline
                                + "of the prerequisites, but could not find: " + dependency + newline
                                + "Please report this issue at " + newline
                                + newline
                                + "    " + pm.website.github.issues.url + newline
                                + newline
                                + "In the meantime, visit, " + newline
                                + newline
                                + "    " + pm.website.home.url + newline
                                + newline
                                + "for instructions to manually install the MPI library on your " + newline
                                + "system. Aborting the automatic MPI installation by ParaMonte..."
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )
                return

        pm.note ( msg = "Installing the Intel MPI library for 64-bit architecture... " + newline
                      + "file location: " + prereqs.mpi.intel.fullFilePath
                , methodName = pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        pm.warn ( msg = "Please do not change the default installation location of " + newline
                      + "the MPI library suggested by the installer. If you do change " + newline
                      + "the default path, the onus will be on you to ensure the path " + newline
                      + "to the MPI runtime libraries exist in the environmental PATH " + newline
                      + "variable of your session."
                , methodName = pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        currentDir = os.getcwd()

        if pm.platform.isWin32:

            err = 0
            err = os.system(prereqs.mpi.intel.fullFilePath)
            if err==0:
                #writeVerificationStatusFile("Testing")
                pm.note ( msg   = "Intel MPI library installation appears to have succeeded. " + newline
                                + "Now close your Python environment and the command-line interface " + newline
                                + "and reopen a new fresh (Anaconda) command prompt."
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )
            else:
                pm.warn ( msg = "Intel MPI library installation might have failed. Exit flag: {}.".format(err)
                    , methodName = pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )

        if pm.platform.isLinux:

            try:

                import tarfile
                tf = tarfile.open(prereqs.mpi.intel.fullFilePath)
                tf.extractall(path=pm.path.lib)
                mpiExtractDir = os.path.join(pm.path.lib, prereqs.mpi.intel.fileName)

                pm.note ( msg   = "If this is your personal computer and you have opened your Python " + newline
                                + "session with superuser (sudo) privileges, then you can choose " + newline
                                + newline
                                + "    'install as root'" + newline
                                + newline
                                + "in the graphical user interface that appears in your session. " + newline
                                + "Otherwise, if you are using the ParaMonte library on a public " + newline
                                + "server, for example, on a supercomputer, or you do not have " + newline
                                + "superuser (sudo) privileges on your system, then choose " + newline
                                + "the third option: " + newline
                                + newline
                                + "   'install as current user'"
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

                mpiInstallScriptPath = os.path.join( mpiExtractDir, "install_GUI.sh" )
                if not os.path.exists(mpiInstallScriptPath):
                    pm.abort( msg   = "Internal error occurred." + newline
                                    + "Failed to detect the Intel MPI installation Bash script." + newline
                                    + "Please report this issue at " + newline
                                    + newline
                                    + "    " + pm.website.github.issues.url + newline
                                    + newline
                                    + "Visit " + pm.website.home.url + " for instructions " + newline
                                    + "to build ParaMonte object files on your system."
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

            except Exception as e:

                print(str(e))
                pm.abort( msg   = "Unzipping of Intel MPI runtime library tarball failed." + newline
                                + "Make sure you have tar software installed on your system and try again."
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

            err = os.system("chmod +x " + mpiInstallScriptPath)
            if err != 0:
                pm.warn ( msg   = "The following action failed: " + newline
                                + newline
                                + "    chmod +x " + mpiInstallScriptPath + newline
                                + newline
                                + "skipping..."
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

            os.chdir(mpiExtractDir)

            import subprocess
            try:
                subprocess.check_call( mpiInstallScriptPath, shell = True )
            except Exception as e:
                print(str(e))
                pm.abort   ( msg   = "Intel MPI runtime libraries installation for " + newline
                                    + "64-bit architecture appears to have failed." + newline
                                    + "Please report this error at:" + newline
                                    + newline
                                    + "    " + pm.website.github.issues.url + newline
                                    + newline
                                    + "Visit " + pm.website.home.url + " for more instructions " + newline
                                    + "to build and use the ParaMonte library on your system."
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

            pm.note ( msg = "Intel MPI runtime libraries installation for " + newline
                          + "64-bit architecture appears to have succeeded. " + newline
                          + "Searching for the MPI runtime environment setup file..."
                    , methodName = pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )

            os.chdir(currentDir)

            setupFilePath = os.path.join( pm.path.auxil, "setup.sh" )

            installationRootDirList = [ "/opt", pm.path.home ]
            mpivarsDefaultFilePathList = ["",""]
            mpiRootDir = ["",""]

            mpiRootDirNotFound = True
            while mpiRootDirNotFound:

                mpiRootDir = []
                mpivarsDefaultFilePathList = []
                mpiTrunkDir = os.path.join( "intel", "compilers_and_libraries_" + prereqs.mpi.intel.version, "linux", "mpi", "intel64" )

                for installationRootDir in installationRootDirList:
                    mpiRootDir.append( os.path.join( installationRootDir, mpiTrunkDir ) )
                    mpivarsDefaultFilePathList.append( os.path.join( mpiRootDir[-1] , "bin" , "mpivars.sh" ) )
                    if os.path.isdir(mpiRootDir[-1]):
                        mpiRootDirNotFound = False
                        break

                if mpiRootDirNotFound:
                    pm.warn ( msg = "Failed to detect the installation root path for Intel MPI runtime " + newline
                                  + "libraries for 64-bit architecture on your system. If you specified " + newline
                                  + "a different installation root path at the time of installation, " + newline
                                  + "please copy and paste it below. Note that the installation root " + newline
                                  + "path is part of the path that replaces: " + newline
                                  + newline
                                  + "    " + "opt" + newline
                                  + newline
                                  + "in the following path: " + newline
                                  + newline
                                  + "    " + os.path.join( "opt" , mpiTrunkDir )
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )
                    answer = input  ( "\n    Please type the root path of MPI installation below and press ENTER."
                                    + "\n    If you don't know the root path, simply press ENTER to quit:\n"
                                    )
                    if len(answer.strip())==0:
                        pm.warn ( msg = "Skipping the MPI runtime library environmental path setup..."
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )
                        break
                    else:
                        installationRootDirList = [ answer ]
                        continue

            if mpiRootDirNotFound:

                pm.warn ( msg   = "Failed to find the MPI runtime environment setup file on your system. " + newline
                                + "This is highly unusual. Normally, Intel MPI libraries is installed " + newline
                                + "in the following directory: " + newline
                                + newline
                                + "    " + mpiRootDir[0] + newline
                                + newline
                                + "or," + newline
                                + newline
                                + "    " + mpiRootDir[1] + newline
                                + newline
                                + "If you cannot manually find the Intel MPI installation directory," + newline
                                + "it is likely that the installation might have somehow failed. " + newline
                                + "If you do find the installation directory, try to locate the " + newline
                                + "'mpivars.sh' file which is normally installed in the following path:" + newline
                                + newline
                                + "    " + mpivarsDefaultFilePathList[0] + newline
                                + newline
                                + "or, " + newline
                                + newline
                                + "    " + mpivarsDefaultFilePathList[1] + newline
                                + newline
                                + "Before attempting to run any parallel ParaMonte simulation, " + newline
                                + "make sure you source this file, like the following: " + newline
                                + newline
                                + "    source " + mpivarsDefaultFilePathList[0] + newline
                                + newline
                                + "or, " + newline
                                + newline
                                + "    source " + mpivarsDefaultFilePathList[1] + newline
                                + newline
                                + "where you will have to replace the path in the above with the " + newline
                                + "correct path that you find on your system."
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

            else:

                mpiBinDir = os.path.join( mpiRootDir[-1], "bin" )
                mpiLibDir = os.path.join( mpiRootDir[-1], "lib" )
                mpivarsFilePath = os.path.join( mpiBinDir, "mpivars.sh" )
                if os.path.isfile(mpivarsFilePath):

                    with open(setupFilePath, "w") as setupFile:
                        setupFile.write(mpiBinDir+"\n")
                        setupFile.write(mpiLibDir+"\n")
                        setupFile.write("source " + mpivarsFilePath)

                    pm.note ( msg = "To ensure all MPI routine environmental variables \n"
                                  + "are properly loaded, source the following Bash script \n"
                                  + "in your Bash environment before calling mpiexec, like:\n\n"
                                  + "    source " + mpivarsFilePath + "\n\n"
                                  + "Alternatively, ParaMonte can also automatically add \n"
                                  + "the required script to your '.bashrc' file, so that \n"
                                  + "all required MPI environmental variables are loaded \n"
                                  + "automatically before any ParaMonte usage from any \n"
                                  + "Bash command line on your system."
                            , methodName = pm.names.paramonte
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
                        os.system( "chmod 777 ~/.bashrc")
                        os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte MPI runtime library initialization >>>' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '" + mpivarsFileCommand + "' >>  ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte MPI runtime library initialization <<<' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

                        pm.note ( msg = "If you intend to run parallel simulations right now,\n"
                                      + "we highly recommned you to close your current shell environment\n"
                                      + "and open a new Bash shell environment. This is to ensure that all MPI\n"
                                      + "library environmental variables are properly set in your shell environment."
                                , methodName = pm.names.paramonte
                                , marginTop = 1
                                , marginBot = 1
                                )

                    #else:
                    #    pm.warn ( msg = "skipping...\n"
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
                    #            , methodName = pm.names.paramonte
                    #            , marginTop = 1
                    #            , marginBot = 1
                    #            )

                else:

                    pm.abort   ( msg   = "ParaMonte was able to detect an MPI library path on your system, however,\n"
                                        + "the MPI installation appears to be corrupted. The required mpivars.sh \n"
                                        + "does not exist:\n\n"
                                        + mpivarsFilePath
                                , methodName = pm.names.paramonte
                                , marginTop = 1
                                , marginBot = 1
                                )

    elif pm.platform.isMacOS:

        pm.warn ( msg   = "To use the ParaMonte kernel routines in parallel on macOS, " + newline
                        + "the Open-MPI library will have to be installed on your system. " + newline
                        #+ "To ensure full consistency, we recommend building the parallel " + newline
                        #+ "object files of ParaMonte library on your system along with Open-MPI." + newline
                        + newline
                        + "If this installation of the prerequisites is being done from within " + newline
                        + "a Jupyter notebook and the installation fails:" + newline
                        + newline
                        + "    1. quit the Jupyter notebook." + newline
                        + "    2. enter an IPython session on the command-prompt:" + newline
                        + "        - On Windows, use Anaconda3 command-prompt." + newline
                        + "        - On Linux / macOS, use the Bash terminal." + newline
                        + "    3. import paramonte as pm" + newline
                        + "    4. pm.verify()" + newline
                        + newline
                        + "Building the ParaMonte library prerequisites on your system..."
                , marginTop = 1
                , marginBot = 1
                , methodName = pm.names.paramonte
                )
        _ = buildParaMontePrereqsForMac()

    else:

        pm.warn ( msg   = "To use ParaMonte in parallel on this unknown Operating System, " + newline
                        + "ParaMonte needs to be built from scratch on your system. " + newline
                        + "Building ParaMonte library prerequisites on your system..."
                , marginTop = 1
                , marginBot = 1
                , methodName = pm.names.paramonte
                )
        build()

####################################################################################################################################
#### buildParaMontePrereqsForMac
####################################################################################################################################

def buildParaMontePrereqsForMac():

    pm.note ( msg = "Checking if Homebrew exists on your system..."
            , methodName = pm.names.paramonte
            , marginTop = 1
            , marginBot = 1
            )

    import shutil
    import subprocess
    if shutil.which("brew") is None:

        pm.note ( msg = "Failed to detect Homebrew on your system. Installing Homebrew..."
                , methodName = pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        err1 = os.system('xcode-select --install')
        if err1 != 0 and not os.path.isdir( subprocess.check_output(['xcode-select','-p']).decode('utf-8').replace("\n","").replace(chr(13),"") ):
            pm.warn ( msg = getMacosInstallHelpMsg("xcode-select")
                    , methodName = pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )
            return False

        #err2 = os.system('ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"')
        err2 = os.system('/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"')
        err3 = os.system('brew --version')
        if err2 != 0 or err3 != 0:
            pm.warn ( msg = getMacosInstallHelpMsg("Homebrew")
                    , methodName = pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )
            return False

    #### cmake

    cmakeInstallationNeeded = False
    cmakePath = shutil.which("cmake")
    if cmakePath is None:
        cmakeInstallationNeeded = True
        pm.note ( msg           = "cmake installation is missing on your system."
                , methodName    = pm.names.paramonte
                , marginTop     = 1
                , marginBot     = 1
                )
    else:
        pm.note ( msg           = "cmake installation detected at: " + cmakePath + newline + "Checking cmake version..."
                , methodName    = pm.names.paramonte
                , marginTop     = 0
                , marginBot     = 0
                )
        try:
            cmakeVersion = str(subprocess.run(args=["cmake","--version"],capture_output=True).stdout).split(" ")[2].split("-")[0]
            cmakeVersionList = cmakeVersion.split(".")
            pm.note ( msg           = "current cmake version: " + cmakeVersion
                    , methodName    = pm.names.paramonte
                    , marginTop     = 0
                    , marginBot     = 0
                    )
            if int(cmakeVersionList[0])>=3 and int(cmakeVersionList[1])>=14:
                pm.note ( msg           = "cmake version is ParaMonte-compatible!"
                        , methodName    = pm.names.paramonte
                        , marginTop     = 0
                        , marginBot     = 0
                        )
            else:
                cmakeInstallationNeeded = True
                pm.note ( msg           = "cmake version is NOT ParaMonte-compatible."
                        , methodName    = pm.names.paramonte
                        , marginTop     = 0
                        , marginBot     = 0
                        )
        except:
            cmakeInstallationNeeded = True
            pm.note ( msg           = "Failed to detect the current cmake installation version. skipping..."
                    , methodName    = pm.names.paramonte
                    , marginTop     = 0
                    , marginBot     = 0
                    )

    if cmakeInstallationNeeded:

        pm.note ( msg = "Installing cmake..."
                , methodName = pm.names.paramonte
                , marginTop = 0
                , marginBot = 0
                )

        err1 = os.system("brew install cmake")
        err2 = os.system("brew link --overwrite cmake")

        if err1 != 0 or err2 != 0:
            pm.warn ( msg = "cmake installation or linking failed."
                    , methodName = pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )
            return False

        cmakeVersionList = str(subprocess.run(args=["cmake","--version"],capture_output=True).stdout).split(" ")[2].split("-")[0].split(".")
        if int(cmakeVersionList[0])>=3 and int(cmakeVersionList[1])>=14:
            pm.note ( msg           = "cmake installation succeeded."
                    , methodName    = pm.names.paramonte
                    , marginTop     = 1
                    , marginBot     = 1
                    )
        else:
            pm.warn ( msg = getMacosInstallHelpMsg("cmake")
                    , methodName = pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )
            return False

    #### gnu

    pm.note ( msg = "Installing GNU Compiler Collection..."
            , methodName = pm.names.paramonte
            , marginTop = 1
            , marginBot = 1
            )

    err1 = os.system("brew install gcc@10")
    err2 = os.system("brew link gcc@10")

    if err1 != 0 or err2 != 0:
        pm.warn ( msg = getMacosInstallHelpMsg("GNU Compiler Collection")
                , methodName = pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )
        return False

    #### open-mpi

    pm.note ( msg = "Installing Open-MPI..."
            , methodName = pm.names.paramonte
            , marginTop = 1
            , marginBot = 1
            )

    err1 = os.system("brew install open-mpi")
    err2 = os.system("brew link open-mpi")

    if err1 != 0 or err2 != 0:
        pm.warn ( msg = getMacosInstallHelpMsg("Open-MPI")
                , methodName = pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )
        return False

    return True

####################################################################################################################################
#### getMacosInstallHelpMsg
####################################################################################################################################

def getMacosInstallHelpMsg(app = ""):
    msg = ("Failed to install and link the application " + app + " on your " + newline
        + "system. The " + app + " application is required to install and " + newline
        + "build the ParaMonte components and prerequisites. " + newline
        + newline
        + "If you are performing this installation from within a Jupyter " + newline
        + "Notebook, then, reinstalling from within an ipython environment or " + newline
        + "the python command line (instead of Jupyter Notebook) will likely " + newline
        + "resolve the errors. To do so, open a Bash command line and type, " + newline
        + newline
        + "    ipython || python" + newline
        + newline
        + "Then, inside the (i)python environment, type, " + newline
        + newline
        + "    import paramonte as pm" + newline
        + "    pm.verify()" + newline
        + newline
        + "Otherwise, you can install the application " + app + " manually " + newline
        + "on your system. The " + app + " installation is only a single " + newline
        + "command and takes only a few seconds to install. " + newline
        + "You can get the installation command from this page: " + newline
        + newline
        + "    " + pm.website.home.install.macos.prereqs.cmd.url + newline
        + newline
        + "Once you have manually installed the missing component, retry, " + newline
        + newline
        + "    import paramonte as pm" + newline
        + "    pm.verify()" + newline
        + newline
        + "skipping the installation for now..."
        )
    return msg

####################################################################################################################################
#### writeVerificationStatusFile
####################################################################################################################################

def writeVerificationStatusFile(verificationEnabledString):
    with open(verificationStatusFilePath, "w") as verificationStatusFile:
        verificationStatusFile.write(verificationEnabledString)
    return None

####################################################################################################################################
#### dispFinalMessage
####################################################################################################################################

def dispFinalMessage():
    pm.note ( msg   = "To check for the MPI library installation status or display the above " + newline
                    + "messages in the future, type the following on the Python command-line: " + newline
                    + newline
                    + "    import paramonte as pm" + newline
                    + "    pm.verify()" + newline
                    + newline
                    + "To get started, type the following on the Python command-line," + newline
                    + newline
                    + "    import paramonte as pm" + newline
                    + "    pm.helpme()"
            , methodName = pm.names.paramonte
            , marginTop = 1
            , marginBot = 1
            )
    return None

####################################################################################################################################
#### displayParaMonteBanner
####################################################################################################################################

def displayParaMonteBanner():
    bannerFilePath = os.path.join( pm.path.auxil, ".ParaMonteBanner")
    offset = ( len(pm.version.interface.dump()) - 5 ) // 2
    print("")
    with open(bannerFilePath,"r") as file:
        for line in file:
            if "Version" in line:
                line = line.replace(" "*offset+"Version 0.0.0","Version "+pm.version.interface.dump())
            print(line,end="")
    print("")
    return None

####################################################################################################################################
#### build
####################################################################################################################################

def build(flags=""):
    """

    Builds the ParaMonte library kernel on the user's system from scratch.

        **Parameters**

            flags

                A string containing any of the ParaMonte install script flags.
                If the operating system is Unix-based (e.g., Linux or macOS) then
                the value of ``flags`` must conform to the rules and syntax of
                the flags of the Bash install script of the ParaMonte library
                on GitHub. If the operating system is Windows, then the value
                of ``flags`` must conform to the rules and syntax of the flags
                of the Batch install script of the ParaMonte library on GitHub.
                The default value is an empty string ``""``.

    """

    if pm.platform.isWin32:

        pm.warn ( msg   = "The ParaMonte library build on Windows Operating Systems (OS) " + newline
                        + "requires the installation of the following software on your system: " + newline
                        + newline
                        + "    - Microsoft Visual Studio (MSVS) (Community Edition >2017)" + newline
                        + "    - Intel Parallel Studio >2018, which is built on top of MSVS" + newline
                        + newline
                        + "If you don't have these software already installed on your system, " + newline
                        + "please visit the following page for the installation instructions: " + newline
                        + newline
                        + "    " + pm.website.home.url + newline
                        + newline
                        + "Follow the instructions on this website for building the ParaMonte " + newline
                        + "ParaMonte on your system."
                , methodName = pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

    else:

        pm.note ( msg   = "You are requesting to build the ParaMonte kernel libraries on your system. " + newline
                        + "The kernel library build requires ParaMonte-compatible versions of the " + newline
                        + "following compilers and parallelism libraries installed on your system: " + newline
                        + newline
                        + "    GNU compiler collection (GCC >8.3)" + newline
                        + "    MPI library (MPICH >3.2) on Linux OS or Open-MPI on Darwin OS" + newline
                        + "    OpenCoarrays >2.8" + newline
                        + newline
                        + "The full installation of these software could require 4 to 5 Gb of free " + newline
                        + "space on your system (where the ParaMonte library is already installed)." + newline
                        + "Note that the installation script is in Bash and therefore requires a " + newline
                        + "Bash or Bash-compatible shell. An existing recent installation of the " + newline
                        + "GNU Compiler Collection (GCC) on your system would be also highly " + newline
                        + "desirable and will significantly cut the build time. Also, downloading " + newline
                        + "the prerequisites requires access to the Internet. If you have an " + newline
                        + "Internet firewall active on your system, please make sure to turn " + newline
                        + "it off before proceeding with the local installation of ParaMonte."
                , methodName = pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )

        buildEnabled = getUserResponse  ( msg   = "\n    Do you wish to download and install the ParaMonte library"
                                                + "\n    and its prerequisites on your system now (y/n)? "
                                        )

        if buildEnabled:

            if pm.platform.isMacOS:
                succeeded = buildParaMontePrereqsForMac()
                if not succeeded:
                    pm.warn ( msg   = "The ParaMonte build failed. To get further instructions " + newline
                                    + "to build the ParaMonte library on your macOS, visit, " + newline
                                    + newline
                                    + "    " + pm.website.home.install.macos.url + newline
                                    + newline
                                    + "You can also report this issue at, " + newline
                                    + newline
                                    + "    " + pm.website.github.issues.url + newline
                                    + newline
                                    + "to get direct help. For more information, visit, " + newline
                                    + newline
                                    + "    " + pm.website.home.url
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )
                    return None

            currentDir = os.getcwd()

            pmGitTarPath = os.path.join( pm.path.root, "master.tar.gz" )
            download( url = pm.website.github.archive.master.tar.url
                    , filePath = pmGitTarPath
                    )

            pmGitRootDir = os.path.join( pm.path.root, "paramonte-master" )

            try:

                import tarfile
                tf = tarfile.open(pmGitTarPath)
                tf.extractall(path=pm.path.root) # path=pmGitRootDir)

                pmGitInstallScriptPath = os.path.join( pmGitRootDir, "install.sh" )
                if not os.path.exists(pmGitInstallScriptPath):
                    pm.abort( msg   = "Internal error occurred." + newline
                                    + "Failed to detect the ParaMonte installation Bash script. " + newline
                                    + "Please report this issue at " + newline
                                    + newline
                                    + "    " + pm.website.github.issues.url + newline
                                    + newline
                                    + "Visit, " + " for instructions " + newline
                                    + newline
                                    + "    " + pm.website.home.url
                                    + "to build ParaMonte object files on your system."
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

            except Exception as e:

                print(str(e))
                pm.abort    ( msg   = "Unzipping of the ParaMonte tarball failed.\n"
                                    + "Make sure you have tar software installed on your system and try again."
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

            err = os.system("chmod +x " + pmGitInstallScriptPath)
            if err != 0:
                pm.warn ( msg   = "The following action failed:\n\n"
                                + "chmod +x " + pmGitInstallScriptPath + "\n\n"
                                + "skipping..."
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

            os.chdir(pmGitRootDir)

            import subprocess
            try:
                os.system( "find " + pmGitRootDir + " -type f -iname \"*.sh\" -exec chmod +x {} \;" )
                os.system( pmGitInstallScriptPath + " --lang python --test_enabled true --exam_enabled false --yes-to-all " + flags )
            except Exception as e:
                print(str(e))
                pm.abort   ( msg   = "The Local installation of ParaMonte failed." + newline
                                    + "Please report this issue at " + newline
                                    + newline
                                    + "    " + pm.website.github.issues.url
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

            os.chdir(currentDir)

            # copy files to module folder

            import glob
            import shutil
            pythonBinDir = os.path.join( pmGitRootDir , "bin" , "Python" , "paramonte" )
            fileList = glob.glob( os.path.join( pythonBinDir , "libparamonte_*" ) )

            if len(fileList)==0:

                pm.abort( msg   = "ParaMonte kernel libraries build and installation " + newline
                                + "appears to have failed. You can check this path:" + newline
                                + newline
                                + "    " + pythonBinDir + newline
                                + newline
                                + "to find out if any shared objects with the prefix " + newline
                                + "'libparamonte_' have been generated or not. " + newline
                                + "Please report this issue at " + newline
                                + newline
                                + "    " + pm.website.github.issues.url
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 2
                        )

            else:

                pm.note ( msg   = "ParaMonte kernel libraries build appears to have succeeded. " + newline
                                + "copying the kernel files to the ParaMonte Python module directory..."
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

                for file in fileList:
                    pm.note ( msg = "file: " + file
                            , methodName = pm.names.paramonte
                            , marginTop = 0
                            , marginBot = 0
                            )
                    shutil.copy(file, pm.path.lib)

                pm.note ( msg = "ParaMonte kernel libraries should be now usable on your system."
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

                setupFilePath = os.path.join( pmGitRootDir , "build", "prerequisites", "prerequisites", "installations", "opencoarrays", "2.8.0", "setup.sh" )

                if os.path.exists(setupFilePath):

                    bashrcContents = getBashrcContents()
                    setupFilePathCmd = "source " + setupFilePath
                    if setupFilePathCmd not in bashrcContents:
                        os.system( "chmod 777 ~/.bashrc")
                        os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '# >>> ParaMonte library local installation setup >>>' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '" + setupFilePathCmd + "' >>  ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '# <<< ParaMonte library local installation setup <<<' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && echo '' >> ~/.bashrc" )
                        os.system( "chmod 777 ~/.bashrc && sh ~/.bashrc" )

                    pm.warn ( msg    = "Whenever you intend to use ParaMonte in the future, " + newline
                                    + "before opening your Python session, please execute " + newline
                                    + "the following command in your Bash shell to ensure " + newline
                                    + "all required paths are properly defined in your " + newline
                                    + "environment: " + newline
                                    + newline
                                    + "    " + setupFilePathCmd + newline
                                    + newline
                                    + "mission accomplished."
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

            writeVerificationStatusFile("True")

        else:

            pm.warn ( msg   = "Aborting the ParaMonte-for-Python local build on your system."
                    , methodName = pm.names.paramonte
                    , marginTop = 1
                    , marginBot = 1
                    )

    return None

####################################################################################################################################
#### getVersionTriplet
####################################################################################################################################

def getVersionTriplet(versionDumpString):
    """
    Take an input version string like, "1.1.1" and return an integer triplet list.
    """
    return np.int32(versionDumpString.split("."))

def getPreviousVersion(currentVerionString):
    """
    Take an input version string like, "1.1.1" and return another string representing the version before the input version, like, 1.1.0.
    """
    currentVerionTriplet = getVersionTriplet(currentVerionString)
    previousVerionString = None
    index = 3
    while True:
        index -= 1
        if index<=0:
            break
        else:
            if currentVerionTriplet[index]>0:
                previousVerionTriplet = copy.deepcopy(currentVerionTriplet)
                previousVerionTriplet[index] -= 1
                previousVerionString = ".".join([str(number) for number in previousVerionTriplet])
                break
    return previousVerionString

####################################################################################################################################
#### getDependencyVersion
####################################################################################################################################

dependencyVersionDict = { "numpy": '1.19.2'
                        , "scipy": '1.5.2'
                        , "pandas": '1.1.2'
                        , "seaborn": '0.11.0'
                        , "matplotlib": '3.3.2'
                        }

def getDependencyVersion( pkg : tp.Optional[ str ] = None ):
    """

    Return the minimum required version of the Python library 
    for the successful use of the ParaMonte library visualization
    and post-processing tools.

        **Parameters**

            pkg

                An optional string representing the name of the 
                Python package whose version is being inquired.

        **Returns**

            A string representing the required minimum version 
            of the input ``pkg``. If ``pkg`` is missing or the 
            package dependency does not exist within the ParaMonte 
            library, the dictionary of all dependencies will 
            be returned.

    """
    if pkg is not None:
        try:
            version = dependencyVersionDict[ pkg ]
        except:
            version = dependencyVersionDict 
    else:
        version = dependencyVersionDict
    return version

####################################################################################################################################
#### getDependencyVersion
####################################################################################################################################

def displayDependencyVersionMessage():
    indentedNewLine = newline + "    "
    pm.note( msg    = "The ParaMonte::Kernel samplers have no Python package dependencies " + newline
                    + "beyond numpy. However, the ParaMonte::Python post-processing and " + newline
                    + "visualization tools require the following Python packages, " + newline
                    + indentedNewLine
                    + indentedNewLine.join("{} : {}".format(key, val) for key, val in dependencyVersionDict.items()) + newline
                    + newline
                    + "If you do not intend to use the postprocessing and visualization tools, " + newline
                    + "you can ignore this message. Otherwise, UPDATE THE ABOVE PACKAGES TO " + newline
                    + "THE REQUESTED VERSIONS OR NEWER, SO THAT THE VISUALIZATION TOOLS " + newline
                    + "OF THE ParaMonte::Python LIBRARY FUNCTION PROPERLY."
            , methodName = pm.names.paramonte
            , marginTop = 1
            , marginBot = 1
            )

####################################################################################################################################
#### verifyDependencyVersion
####################################################################################################################################

def verifyDependencyVersion():
    """

    Verify the existence of the required Python packages and 
    their minimum versions on the current system.

        **Parameters**

            None

        **Returns**

            None

    """
    print("")
    for module, version in dependencyVersionDict.items():

        versionIsCompatible = False
        print("checking the ParaMonte::Python dependency on " + module + " ... ", end = "")

        try:

            exec("import " + module)
            installedVersion = eval(module + ".__version__")

            if installedVersion == version:

                versionIsCompatible = True

            else:

                if installedVersion.split(".")[0] != version.split(".")[0]:
                    pm.warn ( msg   = "The current installation version of the " + module + " library on" + newline
                                    + "your system (" + installedVersion + ") is significantly different from " + newline
                                    + "the version (" + version + ") with which the ParaMonte library " + newline
                                    + "has been tested. This could potentially create runtime issues. " + newline
                                    + "Please consider upgrading this library to the most recent " + newline
                                    + "version by typing the following on your command prompt, " + newline
                                    + newline
                                    + "    pip install --user --upgrade " + module + newline
                                    + newline
                                    + "before you begin to use the ParaMonte library. Should "
                                    + "the simulations or the post-processing of the " + newline
                                    + "output files fail, please report it at, " + newline
                                    + newline
                                    + "    " + pm.website.github.issues.url + newline
                                    + newline
                                    + "for a possible solution. skipping for now..."
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

                if installedVersion.split(".")[1] < version.split(".")[1]:
                    pm.warn ( msg   = "The current installation version of the " + module + " library on" + newline
                                    + "your system (" + installedVersion + ") is not the same as the " + newline
                                    + "version (" + version + ") with which the ParaMonte library " + newline
                                    + "has been tested. This may not create any issues, however, " + newline
                                    + "should the simulations or the post-processing of the " + newline
                                    + "output files fail, please upgrade the library via, " + newline
                                    + newline
                                    + "    pip install --user --upgrade " + module + newline
                                    + newline
                                    + "If the error persists, please report it at, " + newline
                                    + newline
                                    + "    " + pm.website.github.issues.url + newline
                                    + newline
                                    + "for a possible solution. skipping for now..."
                            , methodName = pm.names.paramonte
                            , marginTop = 1
                            , marginBot = 1
                            )

        except Exception as e:

            print(str(e))

            if module=="numpy" or module=="pandas":
                pm.abort( msg   = "Failed to import the " + module + " library into your Python session." + newline
                                + "This library is required for the ParaMonte kernel library to perform " + newline
                                + "simulations. Please install the latest version of this " + newline
                                + "library by typing the following on your command prompt: " + newline
                                + newline
                                + "    pip install --user --upgrade " + module + newline
                                + newline
                                + "If the error persists, please report it at, " + newline
                                + newline
                                + "    " + pm.website.github.issues.url + newline
                                + newline
                                + "for a possible solution. skipping for now..."
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )
            else:
                pm.warn ( msg   = "Failed to import the " + module + " library into your Python session." + newline
                                + "This library is required for the post-processing of the ParaMonte " + newline
                                + "simulation output files. Please install the latest version of this " + newline
                                + "library by typing the following on your command prompt: " + newline
                                + newline
                                + "    pip install --user --upgrade " + module + newline
                                + newline
                                + "If the error persists, please report it at, " + newline
                                + newline
                                + "    " + pm.website.github.issues.url + newline
                                + newline
                                + "for a possible solution. skipping for now..."
                        , methodName = pm.names.paramonte
                        , marginTop = 1
                        , marginBot = 1
                        )

        if versionIsCompatible: print("OK")

    return None

####################################################################################################################################
#### verifyDependencyVersion
####################################################################################################################################

def checkForUpdate(package = "paramonte"):
    import subprocess
    import sys
    latestVersion = str(subprocess.run([sys.executable, '-m', 'pip', 'install', '{}==random'.format(package)], capture_output=True, text=True))
    latestVersion = latestVersion[latestVersion.find('(from versions:')+15:]
    latestVersion = latestVersion[:latestVersion.find(')')]
    latestVersion = latestVersion.replace(' ','').split(',')[-1]

    #currentVersion = str(subprocess.run([sys.executable, '-m', 'pip', 'show', '{}'.format(package)], capture_output=True, text=True))
    #currentVersion = currentVersion[currentVersion.find('Version:')+8:]
    #currentVersion = currentVersion[:currentVersion.find('\\n')].replace(' ','')
    currentVersion = pm.version.interface.dump()

    if latestVersion == currentVersion:
        pm.note ( msg   = "You have the latest version of the ParaMonte library. " + newline
                        + "To see the most recent changes to the library, visit, " + newline
                        + newline
                        + "    " + pm.website.home.overview.changes.python.url
        , methodName = pm.names.paramonte
        , marginTop = 1
        , marginBot = 1
        )
    else:

        currentVersionTriplet = currentVersion.split(".")
        latestVersionTriplet = latestVersion.split(".")
        newerVersionAvailable = True
        for current, latest in currentVersionTriplet, latestVersionTriplet:
            if int(current)>latest:
                newerVersionAvailable = False
                return
                break

        if newerVersionAvailable:
            msg = ("A newer version (" + latestVersion + ") of the ParaMonte library appears " + newline
                + "to be available on the PyPI repository. The currently-installed version is: " + currentVersion + newline
                + "You can upgrade to the latest version by typing the following on " + newline
                + "your Bash terminal or Anaconda command prompt: " + newline
                + newline
                + "    pip install --user --upgrade " + package + newline
                + newline
                + "To upgrade from within your Jupyter or IPython session, try, " + newline
                + newline
                + "    !pip install --user --upgrade " + package + newline
                + newline
                + "To see the latest changes to the ParaMonte Python library, visit, " + newline
                + newline
                + "    " + pm.website.home.overview.changes.python.url
                )
        else:
            msg = "Looks like you have a version of ParaMonte that is newer than the PyPI version. Good for you!"

        pm.note ( msg = msg
                , methodName = pm.names.paramonte
                , marginTop = 1
                , marginBot = 1
                )
            
    return None

####################################################################################################################################

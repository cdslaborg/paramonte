!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   WARNING: The Doxygen 2.9 documenter has extreme difficulty parsing the documentation of this module
!   because of the presence of so many Windows-style directory separators (backslash) in the text.<br>
!   As such, any modification to the documentation must be done with extreme care,
!   otherwise, frequent hard-to-debug runtime Doxygen segfaults are likely.<br>

!>  \brief
!>  This module contains classes and procedures for manipulating system file/folder paths.<br>
!>
!>  \details
!>  This module contains constants, routines, and functionalities for performing actions
!>  that are mostly if not entirely independent of the underlying processor runtime shell environment.<br>
!>  The following is a collection of tips and definitions regarding paths in Windows and Unix-like systems.<br>
!>      -#  A **path** is a string of characters used to uniquely identify a location in a directory structure.<br>
!>          It is composed by following the directory tree hierarchy in which components,
!>          separated by a **path separator** character, represent each directory.<br>
!>      -#  The **delimiting character** is most commonly
!>              1.  the (forward) slash (`/`),
!>              2.  the backslash character (`\`), or
!>              3.  the colon (`:`) character,
!>          though some operating systems may use a different separator.<br>
!>      -#  Paths are used extensively in computer science to represent the directory/file relationships common in modern operating systems
!>          and are essential in the construction of **Uniform Resource Locators** (**URLs**).<br>
!>          Resources can be represented by either **absolute** or **relative** paths.<br>
!>      -#  Around 1970, [Unix](https://en.wikipedia.org/wiki/Unix) introduced the slash character (`/`) as its directory separator
!>          and used the dash character `-` as a command line *switch prefix*.<br>
!>      -#  In 1981, the first version of [Microsoft DOS](https://en.wikipedia.org/wiki/DOS) was released.<br>
!>          MS-DOS 1.0 did not support file directories. Also, a major portion of the utility commands packaged with MS-DOS 1.0
!>          came from IBM and their command line syntax used the slash character `/` as a *switch prefix*.<br>
!>          For example, `dir /w` runs the dir command with the wide list format option.<br>
!>          These fundamental difference between Unix and Windows OS still pervade.<br>
!>          When directory support was added to MS-DOS in version 2.0, `/` was kept as the *switch prefix* character for backwards compatibility.<br>
!>          Microsoft chose the backslash character (`\`) as a *directory separator*, which looks similar to the slash character.<br>
!>          However, **modern versions of Windows are slash-agnostic**, allowing a mix of both types of directory separators in a paths.<br>
!>      -#  An **absolute** or **full** **path** points to the same location in a file system, regardless of the current working directory.<br>
!>          To do that, it must include the root directory.<br>
!>      -#  By contrast, a **relative path** starts from some given working directory, avoiding the need to provide the full absolute path.<br>
!>          *A filename can be considered as a relative path based at the current working directory.*<br>
!>
!>  \anchor syspath_pmod_details_unix
!>  On a **Unix-like OS** and most [POSIX](https://en.wikipedia.org/wiki/POSIX)-compliant OS, [⛓](@ref syspath_pmod_details_unix)<br>
!>  <ol>
!>      <li>    The **root directory** is represented by `/`.<br>
!>      <li>    The **directory separator** is represented by `/`.<br>
!>      <li>    The **current directory** is represented by `.`.<br>
!>      <li>    The **parent directory** is represented by `..`.<br>
!>      <li>    The **home directory** is represented by `~`.<br>
!>  </ol>
!>  Examples:<br>
!>  \code{.sh}
!>      /home/user/docs/letter.txt
!>      ./inThisDir
!>      ../../greatGrandParent
!>      ~/.rcinfo
!>  \endcode
!>
!>  \anchor syspath_pmod_details_windows
!>  On **Windows**, [⛓](@ref syspath_pmod_details_windows)<br>
!>      -#  A **standard DOS** or Windows Command Prompt (CMD) path can consist of three components:
!>              -#  A **volume** or **drive letter** followed by the volume separator (:).<br>
!>              -#  A **directory name**. The directory separator character separates subdirectories within the nested directory hierarchy.<br>
!>              -#  An **optional filename**. The directory separator character separates the file path and the filename.<br>
!>
!>          If all three components are present, the path is **absolute**.<br>
!>          If no volume or drive letter is specified and the directory name begins with the directory separator character, the path is **relative** from the root of the current drive.<br>
!>          Otherwise, the path is relative to the current directory.<br>
!>          Although **UNC paths** (beginning with two directory separators as described below) are not recognized by the standard DOS or CMD, they are also **considered as absolute paths**.<br>
!>          The following table shows some possible directory and file paths.<br>
!>          Path                                                        | Description
!>          :-----------------------------------------------------------| :-----------------------------------------------------------------------------
!>          \f$\ms{C:\Documents\Newsletters\Summer2018.pdf}\f$          | An absolute file path from the root of drive `C:`.
!>          \f$\ms{\Program Files\Custom Utilities\StringFinder.exe}\f$ | An absolute path from the root of the current drive.
!>          \f$\ms{2018\January.xlsx}\f$                                | A relative path to a file in a subdirectory of the current directory.
!>          \f$\ms{..\Publications\TravelBrochure.pdf}\f$               | A relative path to a file in a directory starting from the current directory.
!>          \f$\ms{C:\Projects\apilibrary\apilibrary.sln}\f$            | An absolute path to a file from the root of drive `C:`.
!>          \f$\ms{C:Projects\apilibrary\apilibrary.sln}\f$             | A relative path from the current directory of the `C:` drive.
!>      -#  The Windows **path namespaces**.<br>
!>          The Windows operating system has a unified object model that points to all resources, including files.<br>
!>          These object paths are accessible from the console window and are exposed to the Win32 layer through a special folder of symbolic links that legacy DOS and UNC paths are mapped to.<br>
!>          This special folder is accessed via the DOS device path syntax, which is one of:
!>              -#  **Win32 File Namespaces**: `\\?\C:\Test\Foo.txt`<br>
!>                  For file I/O, the `\\?\` prefix to a path string tells the Windows APIs to disable all string parsing and to send the string that follows it straight to the file system.<br>
!>                  For example, if the file system supports large paths and file names, you can exceed the \f$\ms{MAX_PATH}\f$ limits that are otherwise enforced by the Windows APIs.<br>
!>                  Because it turns off automatic expansion of the path string, the `\\?\` prefix also allows the use of `..` and `.` in the path names, which can be useful
!>                  if you are attempting to perform operations on a file with these otherwise reserved relative path specifiers as part of the fully qualified path.<br>
!>              -#  **Win32 Device Namespaces**: `\\.\C:\Test\Foo.txt`<br>
!>                  The `\\.\` prefix will access the Win32 device namespace instead of the Win32 file namespace.<br>
!>                  This is how access to physical disks and volumes is accomplished directly,
!>                  without going through the file system, if the API supports this type of access.<br>
!>                  <br>
!>
!>      -#  A Windows DOS device path consists of the following components:
!>              -#  The device path specifier (`\\.\` or `\\?\`), which identifies the path as a DOS device path.<br>
!>              -#  A symbolic link to the *real* device object (\f$\ms{C:\\}\f$ in the case of a drive name, or \f$\ms{Volume{b75e2c83-0000-0000-0000-602f00000000}}\f$ in the case of a volume GUID).<br>
!>                  The first segment of the DOS device path after the device path specifier identifies the volume or drive.<br>
!>                  (For example, \f$\ms{\\\\?\\C:\\}\f$ and \f$\ms{\\\\.\\BootPartition\\}\f$)<br>
!>                  There is a specific link for UNCs that is called, not surprisingly, UNC.<br>
!>                  For example: \f$\ms{\\\\.\\UNC\\Server\\Share\\Test\\Foo.txt}\f$ or \f$\ms{\\\\?\\UNC\\Server\\Share\\Test\\Foo.txt}\f$
!>                  <br>
!>
!>      -#  A **standard Universal naming convention (UNC) path** which are used to access network resources, have the following format:
!>              -#  A server or host name, which is prefaced by `\\`. The server name can be a NetBIOS machine name or an IP/FQDN address (IPv4 as well as v6 are supported).<br>
!>              -#  A share name, which is separated from the host name by `\`. Together, the server and share name make up the volume.<br>
!>              -#  A directory name. The directory separator character separates subdirectories within the nested directory hierarchy.<br>
!>              -#  An optional filename. The directory separator character separates the file path and the filename.<br>
!>
!>          The following are some examples of UNC paths.<br>
!>          Path                                                        | Description
!>          :-----------------------------------------------------------| :-----------------------------------------------------------------------------
!>          \f[\ms{\\\\system07\\C\$\\}\f]                              | The root directory of the C: drive on \f$\ms{system07}\f$.
!>          \f$\ms{\\\\Server2\Share\Test\Foo.txt}\f$                   | The Foo.txt file in the Test directory of the \f$\ms{\\\\Server2\Share}\f$ volume.
!>
!>          UNC paths must always be fully qualified (absolute). They can include relative directory segments (. and ..),
!>          but these must be part of a fully qualified path. Relative paths can be used only by mapping a UNC path to a drive letter.
!>
!>      -#  The **root directory** in Windows can be represented by either of the following,<br>
!>          \verbatim
                !>\ (relative to current working directory root)
                !><drive_letter>:
                !>\\<server>\<sharename>
                !>\\?\<drive_spec>:
                !>\\?\<server>\<sharename>
                !>\\?\UNC\<server>\<sharename>
                !>\\.\<physical_device>
!>          \endverbatim
!>          where `<something>` is just a placeholder for `something` (The less-than and greater-than symbols have to be dropped).<br>
!>      -#  The **directory separator** is represented by `\`.<br>
!>      -#  The **current directory** is represented by `\` or `/`.<br>
!>      -#  The **parent directory** is represented by `..`.<br>
!>      -#  The **home directory** is represented by `%%userprofile%`.<br>
!>          Note that `%%userprofile%` expands to `%%SystemDrive%\Users\{username}` and is equivalent to the `$HOME` environment variable in Unix/Linux.<br>
!>          Note that `%%homedrive%%%homepath%` also frequently (but not always) refers to **home directory** on Windows.<br>
!>          Note that the `%%homedrive%` environment variable typically expands to the network root profile directory.<br>
!>          Note that the `%%systemdrive%` (typically `c:`) is the partition with the directory `%%systemroot%` (typically <tt>C:\\WINDOWS</tt>).<br>
!>      -#  The `%%OS%` environment variable expands to operating system name on the user workstation.<br>
!>      -#  **Windows path normalization**.<br>
!>          Almost all paths passed to Windows APIs are normalized. During normalization, Windows performs the following steps:
!>              -#  Identifies the path.<br>
!>                  The first step in path normalization is identifying the type of path. Paths fall into one of a few categories:
!>                      -#  They are device paths; that is, they begin with two separators and a question mark or period (\f$\ms{\\\\?}\f$ or \f$\ms{\\\\.}\f$).<br>
!>                      -#  They are UNC paths; that is, they begin with two separators without a question mark or period.<br>
!>                      -#  They are fully qualified DOS paths; that is, they begin with a drive letter, a volume separator, and a component separator (\f$\ms{C:\\}\f$).<br>
!>                      -#  They designate a legacy device (`CON`, `LPT1`).<br>
!>                      -#  They are relative to the root of the current drive; that is, they begin with a single component separator (`\`).<br>
!>                      -#  They are relative to the current directory of a specified drive; that is, they begin with a drive letter, a volume separator, and no component separator (\f$\ms{C:}\f$).<br>
!>                      -#  They are relative to the current directory; that is, they begin with anything else (\f$\ms{temp\\testfile.txt}\f$).<br>
!>
!>                  The type of the path determines whether or not a current directory is applied in some way. It also determines what the "root" of the path is.<br>
!>                  If the path is a legacy DOS device such as `CON`, `COM1`, or `LPT1`, it is converted into a device path by prepending \f$\ms{\\\\.\\}\f$ and returned.<br>
!>                  A path that begins with a legacy device name is always interpreted as a legacy device.<br>
!>                  For example, the DOS device path for \f$\ms{CON.TXT}\f$ is \f$\ms{\\\\.\\CON}\f$, and the DOS device path for \f$\ms{COM1.TXT\\file1.txt}\f$ is \f$\ms{\\\\.\\COM1}\f$.<br>
!>              -#  Applies the current directory to partially qualified (relative) paths.<br>
!>                      -#  If a path isn not fully qualified, Windows applies the current directory to it.<br>
!>                          UNCs and device paths do not have the current directory applied. Neither does a full drive with separator \f$\ms{C:\\}\f$.<br>
!>                      -#  If the path starts with a single component separator, the drive from the current directory is applied.<br>
!>                          For example, if the file path is \f$\ms{\\utilities}\f$ and the current directory is \f$\ms{C:\\temp\\}\f$, normalization produces \f$\ms{C:\\utilities}\f$.<br>
!>                      -#  If the path starts with a drive letter, volume separator, and no component separator, the last current directory set from the command shell for the specified drive is applied.<br>
!>                          If the last current directory was not set, the drive alone is applied.<br>
!>                          For example, if the file path is \f$\ms{D:sources}\f$, the current directory is \f$\ms{C:\\Documents\\}\f$,
!>                          and the last current directory on drive \f$\ms{D:}\f$ was \f$\ms{D:\\sources\\}\f$, the result is \f$\ms{D:\\sources\\sources}\f$.<br>
!>                          These *drive-relative* paths are a common source of program and script logic errors.<br>
!>                          Assuming that a path beginning with a letter and a colon isn not relative is obviously incorrect.<br>
!>              -#  Canonicalizes component and directory separators.<br>
!>                      -#  All forward slashes (`/`) are converted into the standard Windows separator, the backslash (`\`).<br>
!>                          If two adjacent slashes are present, they are collapsed into a single slash.<br>
!>              -#  Evaluates relative directory components (`.` for the current directory and `..` for the parent directory).<br>
!>                  As the path is processed, any components or segments that are composed of a single or a double period (`.` or `..`) are evaluated:<br>
!>                      -#  For a single period, the current segment is removed, since it refers to the current directory.<br>
!>                      -#  For a double period, the current segment and the parent segment are removed, since the double period refers to the parent directory.<br>
!>                      -#  Parent directories are only removed if they are not past the root of the path.<br>
!>                          The root of the path depends on the type of path.<br>
!>                          It is the drive (\f$\ms{C:\\}\f$) for DOS paths, the server/share for UNCs (\f$\ms{\\\\Server\\Share}\f$),
!>                          and the device path prefix for device paths (\f$\ms{\\\\?\\}\f$ or \f$\ms{\\\\.\\}\f$).<br>
!>              -#  Trims certain characters.<br>
!>                  Along with the runs of separators and relative segments removed earlier, some additional characters are removed during normalization:<br>
!>                      -#  If a segment ends in a single period, that period is removed.<br>
!>                          A segment of a single or double period is normalized in the previous step.<br>
!>                          A segment of three or more periods is not normalized and is actually a valid file/directory name.<br>
!>                      -#  If the path does not end in a separator, all trailing periods and spaces (U+0020) are removed.<br>
!>                          If the last segment is simply a single or double period, it falls under the relative components rule above.<br>
!>                          This rule means that one can create a directory name with a trailing space by adding a trailing separator after the space.<br>
!>                          \warning
!>                          You should never create a directory or filename with a trailing space.<br>
!>                          Trailing spaces can make it difficult or impossible to access a directory,
!>                          and applications commonly fail when attempting to handle directories or files whose names include trailing spaces.<br>
!>
!>      -#  **Skipping the Windows path normalization**.<br>
!>          Normally, any path passed to a Windows API is (effectively) passed to the GetFullPathName of Windows API function and normalized.<br>
!>          There is one important exception: a device path that begins with a question mark instead of a period.<br>
!>          Unless the path starts exactly with \f$\ms{\\\\?\\}\f$ (note the use of the canonical backslash), it is normalized.<br>
!>          Why would one want to skip normalization? There are three major reasons:<br>
!>              -#  To get access to paths that are normally unavailable but are legal.<br>
!>                  A file or directory called hidden, for example, is impossible to access in any other way.<br>
!>              -#  To improve performance by skipping normalization if you have already normalized.<br>
!>              -#  On \f$\ms{.NET}\f$ Framework only, to skip the \f$\ms{MAX_PATH}\f$ check for path length to allow for paths that are greater than 259 characters.<br>
!>                  Most APIs allow this, with some exceptions.<br>
!>
!>          Skipping normalization and max path checks is the only difference between the two device path syntaxes (\f$\ms{\\\\?}\f$ and \f$\ms{\\\\.}\f$).<br>
!>          They are otherwise identical. Skipping Windows path normalization can be dangerous, since one can easily create paths that are difficult for other applications to deal with.<br>
!>      -#  **Case-sensitivity and the Windows file system**.<br>
!>          A peculiarity of the Windows file system that non-Windows users and developers find confusing is that path and directory names are case-insensitive.<br>
!>          That is, directory and file names reflect the casing of the strings used when they are created.<br>
!>          If you rename a directory or file to change its case, the directory or file name reflects the case of the string used when you rename it.<br>
!>          However, directory and file name **comparisons are case-insensitive**.<br>
!>          If you search for a file named `"test.txt"`, the \f$\ms{.NET}\f$ file system APIs ignore case in the comparison.<br>
!>          For example, `"Test.txt"`, `"TEST.TXT"`, `"test.TXT"`, and any other combination of uppercase and lowercase letters will match `"test.txt"`.<br>
!>
!>  For more information, see the [Microsoft documentation for path conventions](https://docs.microsoft.com/en-us/dotnet/standard/io/file-path-formats).<br>
!>
!>  \remark
!>  Since the original implementation of this module, other Fortran Path libraries have also come to existence. Notably, there is
!>  this [MIT-licensed C++-based path manipulation Fortran module](https://github.com/scivision/fortran-pathlib/blob/main/API.md).<br>
!>
!>  \remark
!>  Most Unix-like Operating Systems (OS), in particular, Linux and macOS, strive to follow the conventions of the POSIX OS standard.<br>
!>  For many applications, it is necessary to identify the path-separator recognized by the processor shell.<br>
!>  On Unix systems, the task is easy. However, on Windows the path-separator symbol depends on the shell being used.<br>
!>  For example, if the codebase is compiled with \gfortran and used within `Cygwin`, other Unix shell emulators, then only the Unix style path-separator is recognized.<br>
!>  In a Windows Batch terminal, however, both POSIX and Windows path-separators can be used.<br>
!>
!>  \see
!>  [pm_io](@ref pm_io)<br>
!>  [pm_sysInfo](@ref pm_sysInfo)<br>
!>  [pm_sysPath](@ref pm_sysPath)<br>
!>  [pm_sysShell](@ref pm_sysShell)<br>
!>
!>  \test
!>  [test_pm_sysPath](@ref test_pm_sysPath)<br>
!>
!>  \todo
!>  \pmed
!>  The following functionalities need to be added to this module,
!>  <ol>
!>      <li>    Changing file permission: `chmod`
!>  </ol>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Friday 3:09 AM, Dec 8, 2017, Dell Medical School, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>      \cond excluded
#if     WINDOWS_ENABLED
#define MAX_PATH  260_IK
#define MAX_NAME  260_IK
#else
#define MAX_PATH 4096_IK
#define MAX_NAME  255_IK
#endif
#if     FORTRAN_ENABLED
#define BIND2(X)
#else
#define BIND2(X)bind(C, name = X)
#endif
!>      \endcond excluded

module pm_sysPath

    use pm_kind, only: IK, LK, SK
    use pm_container, only: css_type

    implicit none

    character(*, SK), parameter ::  MODULE_NAME = SK_"@pm_sysPath"

    !>  \cond excluded
    integer         , private   ::  i
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    integer(IK)     , parameter ::  MAX_LEN_FILE_NAME = MAX_NAME    !<  \public The scalar `integer` constant of default kind \IK, representing the maximum system-allowed path length (`260` on Windows,  `255` on Unix).<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MAX_LEN_FILE_NAME
#endif

    integer(IK)     , parameter ::  MAX_LEN_FILE_PATH = MAX_PATH    !<  \public The scalar `integer` constant of default kind \IK, representing the maximum system-allowed path length (`260` on Windows, `4096` on Unix).<br>
                                                                    !!          Note that the maximum file length path supported by the Intel compilers is also `4096`.<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: MAX_LEN_FILE_PATH
#endif
    character(*, SK), parameter ::  DIR_SEP_WINDOWS_ALL = SK_"\/"   !<  \public The scalar `character` constant of default kind \SK, containing the directory separator characters recognized by the Windows operating systems.<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: DIR_SEP_WINDOWS_ALL
#endif
    character(*, SK), parameter ::  DIR_SEP_POSIX_ALL = SK_"/"      !<  \public The scalar `character` constant of default kind \SK, containing the directory separator characters recognized by the POSIX-compliant systems.<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: DIR_SEP_POSIX_ALL
#endif
    character(*, SK), parameter ::  DIR_SEP_WINDOWS = SK_"\"        !<  \public The scalar `character` constant of default kind \SK, containing the preferred directory separator character on Windows operating systems.<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: DIR_SEP_WINDOWS
#endif
    character(*, SK), parameter ::  DIR_SEP_POSIX = SK_"/"          !<  \public The scalar `character` constant of default kind \SK, containing the preferred directory separator character on POSIX-compliant systems.<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: DIR_SEP_POSIX
#endif
    character(*, SK), parameter ::  DIR_SEP_ALL = SK_"/\"           !<  \public The scalar `character` constant of default kind \SK, containing the directory separator characters recognized by all platforms supported in this library (POSIX and Windows).<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: DIR_SEP_ALL
#endif

    character(*, SK), parameter ::  PATH_SEP_WINDOWS = SK_";"       !<  \public The scalar `character` constant of default kind \SK, containing the preferred directory separator character on Windows operating systems.<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: PATH_SEP_WINDOWS
#endif
    character(*, SK), parameter ::  PATH_SEP_POSIX = SK_":"         !<  \public The scalar `character` constant of default kind \SK, containing the preferred directory separator character on POSIX-compliant systems.<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: PATH_SEP_POSIX
#endif
    character(*, SK), parameter ::  PATH_SEP_ALL = SK_":;"          !<  \public The scalar `character` constant of default kind \SK, containing the directory separator characters recognized by all platforms supported in this library (POSIX and Windows).<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: PATH_SEP_ALL
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type with no components, whose instantiated objects are
    !>  used to signify the **verbatim** interpretation of paths when passed to various procedures within the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the generic interfaces that take objects of this type as input arguments.<br>
    !>
    !>  \see
    !>  [verbatim](@ref pm_sysPath::verbatim)<br>
    !>  [verbatim_type](@ref pm_sysPath::verbatim_type)<br>
    !>
    !>  \final{verbatim_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:56 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: verbatim_type; end type

    !>  \brief
    !>  The scalar constant of type [verbatim_type](@ref pm_sysPath::verbatim_type) that can be
    !>  used to signify the **verbatim** interpretation of paths when passed to various procedures within the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the generic interfaces that take objects of this type as input arguments.<br>
    !>
    !>  \see
    !>  [verbatim](@ref pm_sysPath::verbatim)<br>
    !>  [verbatim_type](@ref pm_sysPath::verbatim_type)<br>
    !>
    !>  \final{verbatim_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:56 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(verbatim_type), parameter :: verbatim = verbatim_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: verbatim
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `character` constant of default kind \SK, containing the prefix operating systems path-separator symbol (`"\"` in Windows, `"/"` in Unix).<br>
    !>
    !>  \details
    !>  Prefixing Windows paths with this namespace allows Windows APIs to disable all string parsing and to send the string that follows
    !>  it straight to the file system. For example, if the file system supports large paths and file names, you can exceed the `MAX_PATH`
    !>  limits that are otherwise enforced by the Windows APIs. For more information about the normal maximum path limitation, see the
    !>  [Windows API documentation](https://docs.microsoft.com/en-us/windows/win32/fileio/naming-a-file#maximum-path-length-limitation).<br>
    !>
    !>  \warning
    !>  The `\\` at the beginning of the namespace indicates that the path should be passed to the system with minimal modification,
    !>  which means that the path cannot contain forward slashes to represent directory separators, or a period `.` to represent the current
    !>  directory, or double dots `..` to represent the parent directory.<br>
    !>
    !>  \remark
    !>  The presence of `\\` at the beginning of a path implies a path that follows the Universal Naming Convention (**UNC**).<br>
    !>  A UNC name of any format always starts with two backslash characters ("\\") followed by an **absolute** path.<br>
    !>  For example the UNC name for Microsoft Subsystem for Linux (WSL2) is `\\wsl$\`.<br>
    !>
    !>  \remark
    !>  The addition of this namespace prefix allows the use of `.` or `..` in path names (without any special meaning),
    !>  as well as relaxing the 260 character path name limit of Windows OS,
    !>  if the underlying file system supports long paths and file names.<br>
    !>
    !>  \final{WIN32_NAMESPACE_FILE}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:56 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter ::  WIN32_NAMESPACE_FILE = SK_"\\?\"
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: WIN32_NAMESPACE_FILE
#endif

    !>  \brief
    !>  The scalar `character` constant of default kind \SK, containing the operating systems path-separator symbol (`"\"` in Windows, `"/"` in Unix).<br>
    !>
    !>  \details
    !>  Prefixing Windows paths with this namespace allows Windows APIs to disable all string parsing and to send the string that follows
    !>  it straight to the file system. For example, if the file system supports large paths and file names, you can exceed the `MAX_PATH`
    !>  limits that are otherwise enforced by the Windows APIs. For more information about the normal maximum path limitation, see the
    !>  [Windows API documentation](https://docs.microsoft.com/en-us/windows/win32/fileio/naming-a-file#maximum-path-length-limitation).<br>
    !>
    !>  \warning
    !>  The `\\` at the beginning of the namespace indicates that the path should be passed to the system with minimal modification,
    !>  which means that the path cannot contain forward slashes to represent directory separators, or a period `.` to represent the current
    !>  directory, or double dots `..` to represent the parent directory.<br>
    !>
    !>  \remark
    !>  The addition of this namespace prefix allows the use of `.` or `..` in path names,
    !>  as well as relaxing the 260 character path name limit of Windows OS,
    !>  if the underlying file system supports long paths and file names.<br>
    !>
    !>  \final{WIN32_NAMESPACE_DEVICE}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:56 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter ::  WIN32_NAMESPACE_DEVICE = SK_"\\.\"
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: WIN32_NAMESPACE_DEVICE
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `character` of default kind \SK containing the ASCII control characters that need special care in paths.<br>
    !>
    !>  \warning
    !>  The zeroth ASCII character (the null character) is intentionally excluded from the list to avoid certain
    !>  runtime failures that can happen due to the convention of terminating the C strings with the null character.<br>
    !>
    !>  \final{ASCII_CONTROL_STR}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:56 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter ::  ASCII_CONTROL_STR = achar( 1, SK)   // &
                                                        achar( 2, SK)   // &
                                                        achar( 3, SK)   // &
                                                        achar( 4, SK)   // &
                                                        achar( 5, SK)   // &
                                                        achar( 6, SK)   // &
                                                        achar( 7, SK)   // &
                                                        achar( 8, SK)   // &
                                                        achar( 9, SK)   // &
                                                        achar(10, SK)   // &
                                                        achar(11, SK)   // &
                                                        achar(12, SK)   // &
                                                        achar(13, SK)   // &
                                                        achar(14, SK)   // &
                                                        achar(15, SK)   // &
                                                        achar(16, SK)   // &
                                                        achar(17, SK)   // &
                                                        achar(18, SK)   // &
                                                        achar(19, SK)   // &
                                                        achar(21, SK)   // &
                                                        achar(22, SK)   // &
                                                        achar(23, SK)   // &
                                                        achar(24, SK)   // &
                                                        achar(25, SK)   // &
                                                        achar(26, SK)   // &
                                                        achar(27, SK)   // &
                                                        achar(28, SK)   // &
                                                        achar(29, SK)   // &
                                                        achar(30, SK)   // &
                                                        achar(31, SK)   // &
                                                        achar(127,SK)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ASCII_CONTROL_STR
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `character` constant of default kind \SK, containing the Windows reserved characters not allowed in filenames.<br>
    !>  Note that the ASCII control characters (0-31) are also forbidden in Windows paths.<br>
    !>
    !>  \warning
    !>  It is critical for this constant to begin with the character `\`.<br>
    !>
    !>  \final{WINDOWS_RESERVED_STR}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:56 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter :: WINDOWS_RESERVED_STR =  SK_'\/<>:"|?*'//ASCII_CONTROL_STR ! do not change the first character here.<br>
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: WINDOWS_RESERVED_STR
#endif

    !>  \brief
    !>  The vector `character` constant of default kind \SK of `len = 1`
    !>  containing the individual characters in [WINDOWS_RESERVED_STR](@ref pm_sysPath::WINDOWS_RESERVED_STR).<br>
    !>
    !>  \warning
    !>  The first character must always be the Windows directory separator `\`.<br>
    !>
    !>  \final{WINDOWS_RESERVED_CHR}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(1, SK), parameter :: WINDOWS_RESERVED_CHR(*) = [( WINDOWS_RESERVED_STR(i:i), i = 1_IK, len(WINDOWS_RESERVED_STR,IK) )]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: WINDOWS_RESERVED_CHR
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `character` constant of default kind \SK, containing the Windows CMD shell metacharacters.<br>
    !>
    !>  \details
    !>  Metacharacters are characters that have special meaning for Windows CMD.<br>
    !>
    !>  \final{WINDOWS_CMD_METACHAR_STR}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:56 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter :: WINDOWS_CMD_METACHAR_STR = SK_'()%!^"<>&|'
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: WINDOWS_CMD_METACHAR_STR
#endif

    !>  \brief
    !>  The vector `character` constant of default kind \SK of `len = 1`
    !>  containing the individual characters in [WINDOWS_CMD_METACHAR_STR](@ref pm_sysPath::WINDOWS_CMD_METACHAR_STR).<br>
    !>
    !>  \final{WINDOWS_CMD_METACHAR_CHR}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(1, SK), parameter :: WINDOWS_CMD_METACHAR_CHR(*) = [( WINDOWS_CMD_METACHAR_STR(i:i), i = 1_IK, len(WINDOWS_CMD_METACHAR_STR,IK) )]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: WINDOWS_CMD_METACHAR_CHR
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The vector `character` constant of default kind \SK, containing the Windows reserved (forbidden) file names.<br>
    !>  No file or directory name can contain these names alone, possibly mixed with leading or trailing blanks.<br>
    !>
    !>  \final{WINDOWS_RESERVED_DEVICE_NAME}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(4, SK), parameter :: WINDOWS_RESERVED_DEVICE_NAME(*) =[ SK_"CON " &
                                                                    , SK_"PRN " &
                                                                    , SK_"AUX " &
                                                                    , SK_"NUL " &
                                                                    , SK_"COM1" &
                                                                    , SK_"COM2" &
                                                                    , SK_"COM3" &
                                                                    , SK_"COM4" &
                                                                    , SK_"COM5" &
                                                                    , SK_"COM6" &
                                                                    , SK_"COM7" &
                                                                    , SK_"COM8" &
                                                                    , SK_"COM9" &
                                                                    , SK_"LPT1" &
                                                                    , SK_"LPT2" &
                                                                    , SK_"LPT3" &
                                                                    , SK_"LPT4" &
                                                                    , SK_"LPT5" &
                                                                    , SK_"LPT6" &
                                                                    , SK_"LPT7" &
                                                                    , SK_"LPT8" &
                                                                    , SK_"LPT9" &
                                                                    , SK_"con " &
                                                                    , SK_"prn " &
                                                                    , SK_"aux " &
                                                                    , SK_"nul " &
                                                                    , SK_"com1" &
                                                                    , SK_"com2" &
                                                                    , SK_"com3" &
                                                                    , SK_"com4" &
                                                                    , SK_"com5" &
                                                                    , SK_"com6" &
                                                                    , SK_"com7" &
                                                                    , SK_"com8" &
                                                                    , SK_"com9" &
                                                                    , SK_"lpt1" &
                                                                    , SK_"lpt2" &
                                                                    , SK_"lpt3" &
                                                                    , SK_"lpt4" &
                                                                    , SK_"lpt5" &
                                                                    , SK_"lpt6" &
                                                                    , SK_"lpt7" &
                                                                    , SK_"lpt8" &
                                                                    , SK_"lpt9" &
                                                                    ]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: WINDOWS_RESERVED_DEVICE_NAME
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `character` constant of default kind \SK, containing the ASCII characters that
    !>  require escaping via `<ESC>` character (if not double-quoted) within the Unix Shell environment.<br>
    !>
    !>  \details
    !>  The ASCII characters needing escape in POSIX-compliant shells can be obtained via the following Bash script,
    !>  \code{.sh}
    !>
    !>      for i in {0..127} ;do
    !>          printf -v var \\%o $i
    !>          printf -v var $var
    !>          printf -v res "%q" "$var"
    !>          esc=E
    !>          [ "$var" = "$res" ] && esc=-
    !>          printf "%02X %s %-7s\n" $i $esc "$res"
    !>      done |
    !>          column
    !>
    !>  \endcode
    !>
    !>  which should render a table with the first field being hexa-value of byte,
    !>  the second containing `E` if character need to be escaped and third field
    !>  showing escaped presentation of character.<br>
    !>
    !>  \final{POSIX_RESERVED_STR}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter :: POSIX_RESERVED_STR = SK_" " // & ! space character
                                                        SK_"!" // & ! history expansion.<br>
                                                        SK_'"' // & ! shell syntax.<br>
                                                        SK_"#" // & ! comment start when preceded by whitespace; zsh wildcards.<br>
                                                        SK_"$" // & ! shell syntax.<br>
                                                        SK_"&" // & ! shell syntax.<br>
                                                        SK_"'" // & ! shell syntax.<br>
                                                        SK_"(" // & ! even in the middle of a word: ksh extended globs (also available in bash and zsh); zsh wildcards.<br>
                                                        SK_")" // & ! even in the middle of a word: ksh extended globs (also available in bash and zsh); zsh wildcards.<br>
                                                        SK_"*" // & ! sh wildcard.<br>
                                                        SK_"," // & ! only inside brace expansion.<br>
                                                        SK_";" // & ! shell syntax.<br>
                                                        SK_"<" // & ! shell syntax.<br>
                                                        SK_"=" // & ! in zsh, when it is at the beginning of a file name (filename expansion with PATH lookup).<br>
                                                        SK_">" // & ! shell syntax.<br>
                                                        SK_"?" // & ! sh wildcard.<br>
                                                        SK_"[" // & ! sh wildcard.<br>
                                                        SK_"\" // & ! shell syntax.<br>
                                                        SK_"]" // & ! you may get away with leaving it unquoted.<br>
                                                        SK_"^" // & ! history expansion; zsh wildcard.<br>
                                                        SK_"`" // & ! shell syntax.<br>
                                                        SK_"{" // & ! brace expansion.<br>
                                                        SK_"|" // & ! shell syntax.<br>
                                                        SK_"}" // & ! needs to be escaped in zsh, other shells are more lenient when there is no matching opening brace.<br>
                                                        SK_"~" // & ! home directory expansion when at the beginning of a filename; zsh wildcard; safe when it is the last character.<br>
                                                        ASCII_CONTROL_STR
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: POSIX_RESERVED_STR
#endif

    !stupid gfortran (possibly version 8.3) gives error with the above syntax.<br>
    !character(*, SK), parameter :: POSIX_RESERVED_STR = " !"//'"#$&'//"'()*,;<=>?[\]^`{|}~"

    !>  \brief
    !>  The vector `character` constant of default kind \SK of `len = 1` containing the individual characters in [POSIX_RESERVED_STR](@ref pm_sysPath::POSIX_RESERVED_STR).<br>
    !>
    !>  \final{POSIX_RESERVED_CHR}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(1, SK), parameter :: POSIX_RESERVED_CHR(*) = [( POSIX_RESERVED_STR(i:i), i = 1_IK, len(POSIX_RESERVED_STR,IK) )]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: POSIX_RESERVED_CHR
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `character` constant of default kind \SK,
    !>  containing the ASCII characters that require escaping via `ESC` (Escape) character within the Unix Shell environment
    !>  even when the string double-quoted. Note that the history expansion character (`!`) is not included since it is
    !>  meaningless in the POSIX standard.<br>
    !>  See the [GNU Bash manual](https://www.gnu.org/software/bash/manual/html_node/Double-Quotes.html) for more information.<br>
    !>
    !>  \final{POSIX_RESERVED_DQUOTE_STR}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(*, SK), parameter :: POSIX_RESERVED_DQUOTE_STR = SK_"`$\"""//achar(10, SK) ! space character
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: POSIX_RESERVED_DQUOTE_STR
#endif

    !>  \brief
    !>  The vector `character` constant of default kind \SK of `len = 1` containing the individual characters in [POSIX_RESERVED_DQUOTE_STR](@ref pm_sysPath::POSIX_RESERVED_DQUOTE_STR).<br>
    !>
    !>  \final{POSIX_RESERVED_DQUOTE_CHR}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    character(1, SK), parameter :: POSIX_RESERVED_DQUOTE_CHR(*) = [( POSIX_RESERVED_DQUOTE_STR(i:i), i = 1_IK, len(POSIX_RESERVED_DQUOTE_STR,IK) )]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: POSIX_RESERVED_DQUOTE_CHR
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!>  \brief
!!>  This is the [PathPart_type](@ref pm_sysPath::PathPart_type) class providing a convenient
!!>  collection of attributes that are frequently needed for handling and manipulating paths.<br>
!!>
!!>  \example{PathPart_type}
!!>  \include{lineno} example/pm_sysPath/PathPart_type/main.F90
!!>  \compile{PathPart_type}
!!>  \output{PathPart_type}
!!>  \include{lineno} example/pm_sysPath/PathPart_type/main.out.F90
!!>
!!>  \final
!!>
!!>  \author
!!>  Amir Shahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
!type :: PathPart_type
!    character(:, SK), allocatable   :: val      !< The `allocatable` scalar `character` of default kind \SK to contain the path value.<br>
!    character(1, SK)                :: sep      !< The scalar `character` of default kind \SK of length `1` to contain the directory separator (`\` or `/`) recognized by the system Shell.<br>
!    integer(IK)                     :: dir(2)   !< The `integer` vector of length `2` of default kind \IK to contain the starting an the ending indices of the directory segment of the path.<br>
!    integer(IK)                     :: name(2)  !< The `integer` vector of length `2` of default kind \IK to contain the starting an the ending indices of the name of the file, if any exists in the path.<br>
!    integer(IK)                     :: base(2)  !< The `integer` vector of length `2` of default kind \IK to contain the starting an the ending indices of the base of the file name, if any exists in the path.<br>
!    integer(IK)                     :: ext(2)   !< The `integer` vector of length `2` of default kind \IK to contain the starting an the ending indices of the file extension, if any exists in the path (including the dot separator).<br>
!end type
!
!
!
!    type :: pathList_type
!        character(:, SK), allocatable   :: parent   !<  \public     The public `allocatable` scalar `character` of default kind \SK containing the directory path whose contents is being returned.<br>
!        integer(IK)     , private       :: unit     !<  \private    The private scalar `integer` of default kind \IK containing the `unit` connected to the file containing the directory listing results.<br>
!        logical(LK)     , private       :: done     !<  \private    The private scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b> the contents of the directory listing is fully exhausted.<br>
!                                                    !!              This logical flag is used properly close and delete the external files associated with an object of this type upon finalization.<br>
!    !contains
!    !    procedure, pass , private       ::                      isFailedGetNextPathList, isFailedGetNextPathListLen
!    !    generic                         :: isFailedGetNext =>   isFailedGetNextPathList, isFailedGetNextPathListLen
!    !    final                           :: destroyPathList
!    end type
!
!    !>  \cond excluded
!    interface pathList_type
!    module function pathList_typer(path, sort, showdir, showfile, showhidden, failed, errmsg) result(PathList)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: pathList_typer
!#endif
!        use pm_kind, only: SKG => SK
!        character(*,SKG)    , intent(in)    , optional  :: path
!        character(*,SKG)    , intent(in)    , optional  :: sort
!        logical(LK)         , intent(in)    , optional  :: showdir, showfile, showhidden
!        logical(LK)         , intent(out)   , optional  :: failed
!        character(*, SK)    , intent(inout) , optional  :: errmsg
!        character(:,SKG)    , allocatable               :: PathList
!    end function
!    end interface
!    !>  \endcond excluded
!
!    !>  \cond excluded
!    interface
!
!    module function isFailedGetNextPathList(PathList, path, errmsg) result(failed)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetNextPathList
!#endif
!        use pm_kind, only: SKG => SK
!        class(pathList_type), intent(inout)                 :: PathList
!        character(*,SKG)    , intent(out)   , allocatable   :: path
!        character(*, SK)    , intent(inout) , optional      :: errmsg
!        logical(LK)                                         :: failed
!    end function
!
!    module function isFailedGetNextPathListLen(PathList, path, pathlen, errmsg) result(failed)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetNextPathListLen
!#endif
!        use pm_kind, only: SKG => SK
!        class(pathList_type), intent(inout)                 :: PathList
!        character(*,SKG)    , intent(out)                   :: path
!        integer(IK)         , intent(out)                   :: pathlen
!        character(*, SK)    , intent(inout) , optional      :: errmsg
!        logical(LK)                                         :: failed
!    end function
!
!    end interface
!    !>  \endcond excluded
!
!
!    !>  \cond excluded
!    interface
!    module subroutine destroyPathList(PathList)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGetNextPathListLen
!#endif
!        use pm_kind, only: SKG => SK
!        type(pathList_type) , intent(inout)                 :: PathList
!    end subroutine
!    end interface
!    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a `list` of all paths within the specified input `path`.<br>
    !>
    !>  \details
    !>  This generic interface is a convenience wrapper around the more verbose [isFailedList](@ref pm_sysPath::isFailedList).<br>
    !>  Unlike [isFailedList](@ref pm_sysPath::isFailedList) which can gracefully handle the occurrence of runtime error,
    !>  this generic interface halts the program by calling `error stop` if a runtime error occurs.<br>
    !>
    !>  \param[in]      path        :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the directory or file pattern to list.<br>
    !>  \param[in]      sort        :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the sorting methodology used for sorting paths in the listing.<br>
    !>                                  The following values are supported,
    !>                                  <ol>
    !>                                      <li>    `sort = "name"`: Sort all items alphabetically by their names.<br>
    !>                                      <li>    `sort = "tmod"`: Sort all items chronologically by their dates of **last modification** (**oldest first**).<br>
    !>                                      <li>    `sort = "tacc"`: Sort all items chronologically by their dates of **last access** (**oldest first**).<br>
    !>                                      <li>    `sort = "size"`: Sort all items by their size (**smallest first**).<br>
    !>                                      <li>    `sort = "fext"`: Sort all items alphabetically by their file name extension (**empty extensions first**).<br>
    !>                                  </ol>
    !>                                  (**optional**, default = `"name"`)
    !>  \param[in]      showdir     :   The input scalar `logical` of default kind \LK. If `.true.`, the subdirectories (if any exist) in the specified path will be also listed.<br>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[in]      showfile    :   The input scalar `logical` of default kind \LK. If `.true.`, the file (if any exist) in the specified path will be also listed.<br>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[in]      showhidden  :   The input scalar `logical` of default kind \LK. If `.true.`, the hidden files (if any exist) in the specified path will be also listed.<br>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[in]      reversed    :   The input scalar `logical` of default kind \LK. If `.true.`, the sorting of the list will be reversed.<br>
    !>                                  (**optional**, default = `.false.`)
    !>
    !>  \return
    !>  `list`                      :   The output `allocatable` vector of containers of type [css_type](@ref pm_container::css_type)
    !>                                  containing scalar strings of kind \SK, each element of which contains the path to one entry in the directory listing.<br>
    !>
    !>  \interface{ls}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK, IK, LK
    !>      use pm_sysPath, only: ls, css_type
    !>      type(css_type), allocatable :: list(:)
    !>      logical(LK) :: showhidden, showdir, showfile
    !>
    !>      list = ls(path, sort = sort, showdir = showdir, showfile = showfile, showhidden = showhidden, reversed = reversed) ! `list` is a vector of string containers.
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  The procedures of this generic interface can handle characters of any kind (including nongraphical ASCII characters like newline) on Unix systems.<br>
    !>
    !>  \see
    !>  [ls](@ref pm_sysPath::ls)<br>
    !>  [isFailedList](@ref pm_sysPath::isFailedList)<br>
    !>
    !>  \example{ls}
    !>  \include{lineno} example/pm_sysPath/ls/main.F90
    !>  \compilef{ls}
    !>  \output{ls}
    !>  \include{lineno} example/pm_sysPath/ls/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \pmed
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of the procedures of this generic interface for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{ls}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface ls

    module function ls_BSSK(path, sort, showdir, showfile, showhidden, reversed) result(list)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: ls_BSSK
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)                    :: path
        character(*, SK)    , intent(in)    , optional      :: sort
        logical(LK)         , intent(in)    , optional      :: showdir, showfile, showhidden, reversed
        type(css_type)      , allocatable                   :: list(:)
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the result of globing the specified input `pattern`.<br>
    !>
    !>  \details
    !>  This generic interface is a convenience wrapper around the generic interface [isFailedGlob](@ref pm_sysPath::isFailedGlob).<br>
    !>  See the documentation of [isFailedGlob](@ref pm_sysPath::isFailedGlob) for more details on globing and the methodology.<br>
    !>  Unlike [isFailedGlob](@ref pm_sysPath::isFailedGlob) which can gracefully handle the occurrence of runtime error,
    !>  this generic interface halts the program by calling `error stop` if a runtime error occurs.<br>
    !>
    !>  \param[in]      pattern     :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the directory or file pattern to glob.<br>
    !>
    !>  \return
    !>  `list`                      :   The output `allocatable` vector of containers of type [css_type](@ref pm_container::css_type)
    !>                                  containing scalar strings of kind \SK, each element of which contains the path to one entry in the globing.<br>
    !>
    !>  \interface{glob}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: glob, css_type
    !>      type(css_type), allocatable :: list(:)
    !>
    !>      list = glob(pattern) ! `list` is a vector of string containers.
    !>      !
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  The procedures of this generic interface can handle characters of any kind (including nongraphical ASCII characters like newline) on Unix systems.<br>
    !>
    !>  \see
    !>  [ls](@ref pm_sysPath::ls)<br>
    !>  [glob](@ref pm_sysPath::glob)<br>
    !>  [isFailedGlob](@ref pm_sysPath::isFailedGlob)<br>
    !>  [isFailedList](@ref pm_sysPath::isFailedList)<br>
    !>
    !>  \example{glob}
    !>  \include{lineno} example/pm_sysPath/glob/main.F90
    !>  \compilef{glob}
    !>  \output{glob}
    !>  \include{lineno} example/pm_sysPath/glob/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of the procedures of this generic interface for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{glob}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface glob

    module function glob_BSSK(pattern) result(list)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: glob_BSSK
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)                    :: pattern
        type(css_type)      , allocatable                   :: list(:)
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the attempt to globing the specified input `pattern` <b>fails</b>, otherwise, return `.false.`.<br>
    !>
    !>  \details
    !>  **Glob patterns** specify sets of filenames with wildcard characters.<br>
    !>  For example, the Unix Bash shell command `mv *.txt textfiles/` moves all files with names ending in `.txt` from the current directory to the directory `textfiles`.<br>
    !>  Here, `*` is a wildcard standing for any string of characters except `*` and `*.txt` is a **glob pattern**.<br>
    !>  The other common wildcard is the question mark (`?`), which stands for one character.<br>
    !>  For example, `mv ?.txt shorttextfiles/` will move all files named with a single character followed by `.txt` from the current directory to directory `shorttextfiles`,
    !>  while `??.txt` would match all files whose name consists of `2` characters followed by `.txt`.<br>
    !>  In addition to matching filenames, globs are also used widely for matching arbitrary strings (wildcard matching).<br>
    !>  In this capacity a common interface is the popular `fnmatch()`.
    !>
    !>  Origin
    !>  ------
    !>
    !>  The glob command, short for **global**, originates in the earliest versions of Bell Labs Unix.<br>
    !>  The command interpreters of the early versions of Unix (1st through 6th Editions, 1969-1975) relied on a separate program to expand wildcard characters in unquoted arguments to a command: `/etc/glob`.<br>
    !>  That program performed the expansion and supplied the expanded list of file paths to the command for execution.<br>
    !>  Glob was originally written in the B programming language.<br>
    !>  It was the first piece of mainline Unix software to be developed in a high-level programming language.<br>
    !>  Later, this functionality was provided as a C library function, `glob()`, used by programs such as the `shell`.<br>
    !>  It is usually defined based on a function named `fnmatch()`, which tests for whether a string matches a given pattern.<br>
    !>  The program using this function can then iterate through a series of strings (usually filenames) to determine which ones match.<br>
    !>  Both functions are a part of POSIX: the functions defined in POSIX.1 since 2001, and the syntax defined in POSIX.2.<br>
    !>  The idea of defining a separate match function started with `wildmat` (wildcard match), a simple library to match strings against Bourne Shell globs.<br>
    !>  Traditionally, globs do not match hidden files in the form of Unix **dotfiles**; to match them the pattern must explicitly start with `.`.<br>
    !>  For example, `*` matches all visible files while `.*` matches all hidden files.
    !>
    !>  The returned `list` of matching paths can be either,
    !>  <ol>
    !>      <li>    a vector of scalar string containers [css_type](@ref pm_container::css_type).<br>
    !>      <li>    a scalar string container, containing the all directory contents concatenated within a
    !>              single string whose records start and stop elements are returned in the output array `index(1:2,:)`.<br>
    !>              This format is particularly efficient as it involved only one allocation compared to the former containerized storage mode.<br>
    !>  </ol>
    !>
    !>  \warning
    !>  The current implementation of this generic interface recognizes only the two most popular wildcards `?` and `*` for globing.<br>
    !>
    !>  \param[in]      pattern     :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the directory or file pattern to glob.<br>
    !>  \param[out]     list        :   The output object that can be either,
    !>                                  <ol>
    !>                                      <li>    an `allocatable` scalar `character` of default kind \SK containing the path entries in sequence,
    !>                                              such that the next path entry starts immediately after the last character of the previous entry.<br>
    !>                                      <li>    an `allocatable` vector of containers of type [css_type](@ref pm_container::css_type) containing scalar strings of kind \SK,
    !>                                              each element of which contains the path to one entry in the directory listing.<br>
    !>                                  </ol>
    !>                                  **Note that the allocation status or the contents of `list` remain undefined if `failed = .true.` on return.**
    !>  \param[out]     index       :   The output `allocatable` array of type `integer` of default kind \IK of shape `(1:2,1:pathCount) where,
    !>                                  <ol>
    !>                                      <li>    `pathCount` represents the number of paths in the output list, and
    !>                                      <li>    `index(1,:)` is the vector of positions of the starting character of individual paths in `list`, and
    !>                                      <li>    `index(2,:)` is the vector of positions of the ending character of individual paths in `list`.<br>
    !>                                  </ol>
    !>                                  By definition, `index(1,1) = 1` and `index(2,pathCount) = len(list)`.<br>
    !>                                  such that the next path entry starts immediately after the last character of the previous entry.<br>
    !>                                  **Note that the allocation status or the contents of `index` remain undefined if `failed = .true.` on return**.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `list` is a scalar `character` of default kind \SK.)
    !>  \param[inout]   errmsg      :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                                  If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                                  A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                                  (**optional**. Its presence is relevant only if `failed` is also present.)
    !>
    !>
    !>  \return
    !>  `failed`                    :   The output scalar `logical` of default kind \LK.<br>
    !>                                  It is `.true.` <b>if and only if</b> an error occurs while globing `pattern`.<br>
    !>
    !>  \interface{isFailedGlob}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK, IK, LK
    !>      use pm_sysPath, only: isFailedGlob
    !>      character(:, SK), allocatable :: list
    !>      type(css_type), allocatable :: list(:)
    !>      integer(IK), allocatable :: index(:,:)
    !>      logical(LK) :: failed
    !>      logical(LK) :: showdir
    !>      logical(LK) :: showfile
    !>      logical(LK) :: showhidden
    !>
    !>      failed = isFailedGlob(pattern, list, errmsg = errmsg) ! `list` is a vector of string containers.
    !>      failed = isFailedGlob(pattern, list, index, errmsg = errmsg) ! `list` is a scalar string.
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  The interface returning `list` as a scalar string of concatenated paths is significantly faster than a container
    !>  implementation by reducing the number of allocations from an arbitrarily large number to `< 3`, frequently only `1`.<br>
    !>  The runtime performance gain is particularly significant when the returned list contains thousands of paths.<br>
    !>
    !>  \note
    !>  The procedures of this generic interface can handle characters of any kind (including nongraphical ASCII characters like newline) on Unix systems.<br>
    !>
    !>  \see
    !>  [ls](@ref pm_sysPath::ls)<br>
    !>  [glob](@ref pm_sysPath::glob)<br>
    !>  [isFailedGlob](@ref pm_sysPath::isFailedGlob)<br>
    !>  [isFailedList](@ref pm_sysPath::isFailedList)<br>
    !>  [glob (programming)](https://en.wikipedia.org/wiki/Glob_(programming))<br>
    !>
    !>  \example{isFailedGlob}
    !>  \include{lineno} example/pm_sysPath/isFailedGlob/main.F90
    !>  \compilef{isFailedGlob}
    !>  \output{isFailedGlob}
    !>  \include{lineno} example/pm_sysPath/isFailedGlob/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of the procedures of this generic interface for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \todo
    !>  \phigh
    !>  The current implementation of this generic interface recognizes only the two most popular wildcards `?` and `*` for globing.<br>
    !>  This must be expanded to other widely accepted wildcards (e.g., [linux wildcards](https://tldp.org/LDP/GNU-Linux-Tools-Summary/html/x11655.htm))
    !>  or at least the wildcards supported on the respective platforms.<br>
    !>
    !>  \todo
    !>  \phigh
    !>  The current implementation of this generic interface relies on access to supported runtime shells.<br>
    !>  If the runtime shell is unrecognized, the algorithm falls back to `bash` if the shell is detected on the system, otherwise, the procedure fails.<br>
    !>  This should be fixed by changing the approach to rely instead on `scandir()` and `fnmatch()` popular interfaces for cross-platform globing.<br>
    !>
    !>  \final{isFailedGlob}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedGlob

    module function isFailedGlob_SK(pattern, list, index, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGlob_SK
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)                    :: pattern
        character(:,SKG)    , intent(out)   , allocatable   :: list
        integer(IK)         , intent(out)   , allocatable   :: index(:,:)
        character(*, SK)    , intent(inout) , optional      :: errmsg
        logical(LK)                                         :: failed
    end function

    module function isFailedGlob_BSSK(pattern, list, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedGlob_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        character(*,SKG)    , intent(in)                    :: pattern
        type(css_type)      , intent(out)   , allocatable   :: list(:)
        character(*, SK)    , intent(inout) , optional      :: errmsg
        logical(LK)                                         :: failed
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the attempt to fetch the directory listing of the specified input `path` **fails**, otherwise, return `.false.`.<br>
    !>
    !>  \details
    !>  This generic interface partially replicates the common behaviors in the `dir`
    !>  command of Windows terminals and the `ls` command of POSIX-compliant shells.<br>
    !>
    !>  The returned `list` of directory contents can be either,<br>
    !>  <ol>
    !>      <li>    a vector of scalar string containers [css_type](@ref pm_container::css_type).<br>
    !>      <li>    a scalar string container, containing the all directory contents concatenated within a
    !>              single string whose records start and stop elements are returned in the output array `index(1:2,:)`.<br>
    !>              This format is particularly efficient as it involved only one allocation compared to the former containerized storage mode.<br>
    !>  </ol>
    !>
    !>  \param[in]      path        :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the directory or file to list.<br>
    !>                                  Note that the input `path` is interpreted verbatim.<br>
    !>                                  See [isFailedGlob](@ref pm_sysPath::isFailedGlob) for pattern matching (globing), the results of which can be passed as input to this generic interface.<br>
    !>  \param[out]     list        :   The output object that can be either,
    !>                                  <ol>
    !>                                      <li>    an `allocatable` scalar `character` of default kind \SK containing the path entries in sequence,
    !>                                              such that the next path entry starts immediately after the last character of the previous entry.<br>
    !>                                      <li>    an `allocatable` vector of containers of type [css_type](@ref pm_container::css_type) containing scalar strings of kind \SK,
    !>                                              each element of which contains the path to one entry in the directory listing.<br>
    !>                                  </ol>
    !>                                  **Note that the allocation status or the contents of `list` remain undefined if `failed = .true.` on return.**
    !>  \param[out]     index       :   The output `allocatable` array of type `integer` of default kind \IK of shape `(1:2,1:pathCount) where,
    !>                                  <ol>
    !>                                      <li>    `pathCount` represents the number of paths in the output list, and
    !>                                      <li>    `index(1,:)` is the vector of positions of the starting character of individual paths in `list`, and
    !>                                      <li>    `index(2,:)` is the vector of positions of the ending character of individual paths in `list`.<br>
    !>                                  </ol>
    !>                                  By definition, `index(1,1) = 1` and `index(2,pathCount) = len(list)`.<br>
    !>                                  such that the next path entry starts immediately after the last character of the previous entry.<br>
    !>                                  **Note that the allocation status or the contents of `index` remain undefined if `failed = .true.` on return**.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `list` is a scalar `character` of default kind \SK.)
    !>  \param[in]      sort        :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the sorting methodology used for sorting paths in the listing.<br>
    !>                                  The following values are supported,
    !>                                  <ol>
    !>                                      <li>    `sort = "name"`: Sort all items alphabetically by their names.<br>
    !>                                      <li>    `sort = "tmod"`: Sort all items chronologically by their dates of **last modification** (**oldest first**).<br>
    !>                                      <li>    `sort = "tacc"`: Sort all items chronologically by their dates of **last access** (**oldest first**).<br>
    !>                                      <li>    `sort = "size"`: Sort all items by their size (**smallest first**).<br>
    !>                                      <li>    `sort = "fext"`: Sort all items alphabetically by their file name extension (**empty extensions first**).<br>
    !>                                  </ol>
    !>                                  (**optional**, default = `"name"`)
    !>  \param[in]      showdir     :   The input scalar `logical` of default kind \LK. If `.true.`, the subdirectories (if any exist) in the specified path will be also listed.<br>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[in]      showfile    :   The input scalar `logical` of default kind \LK. If `.true.`, the file (if any exist) in the specified path will be also listed.<br>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[in]      showhidden  :   The input scalar `logical` of default kind \LK. If `.true.`, the hidden files (if any exist) in the specified path will be also listed.<br>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[in]      reversed    :   The input scalar `logical` of default kind \LK. If `.true.`, the sorting of the list will be reversed.<br>
    !>                                  (**optional**, default = `.false.`)
    !>  \param[inout]   errmsg      :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                                  If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                                  A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                                  (**optional**. Its presence is relevant only if `failed` is also present.)
    !>
    !>
    !>  \return
    !>  `failed`                    :   The output scalar `logical` of default kind \LK.<br>
    !>                                  It is `.true.` <b>if and only if</b> an error occurs while resolving `path`.<br>
    !>
    !>  \interface{isFailedList}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK, IK, LK
    !>      use pm_sysPath, only: isFailedList
    !>      character(:, SK), allocatable :: string
    !>      type(css_type), allocatable :: list(:)
    !>      integer(IK), allocatable :: index(:,:)
    !>      logical(LK) :: failed
    !>      logical(LK) :: showdir
    !>      logical(LK) :: showfile
    !>      logical(LK) :: showhidden
    !>
    !>      failed = isFailedList(path, list(:), sort = sort, showdir = showdir, showfile = showfile, showhidden = showhidden, reversed = reversed, errmsg = errmsg) ! `list` is a vector of string containers.
    !>      failed = isFailedList(path, string, index, sort = sort, showdir = showdir, showfile = showfile, showhidden = showhidden, reversed = reversed, errmsg = errmsg) ! `string` is a scalar string.
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  The interface returning `list` as a scalar string of concatenated paths is significantly faster than a container
    !>  implementation by reducing the number of allocations from an arbitrarily large number to `< 3`, frequently only `1`.<br>
    !>  The runtime performance gain is particularly significant when the returned list contains thousands of paths.<br>
    !>
    !>  \note
    !>  The procedures of this generic interface can handle characters of any kind (including nongraphical ASCII characters like newline) on Unix systems.<br>
    !>
    !>  \see
    !>  [ls](@ref pm_sysPath::ls)<br>
    !>  [isFailedList](@ref pm_sysPath::isFailedList)<br>
    !>
    !>  \example{isFailedList}
    !>  \include{lineno} example/pm_sysPath/isFailedList/main.F90
    !>  \compilef{isFailedList}
    !>  \output{isFailedList}
    !>  \include{lineno} example/pm_sysPath/isFailedList/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of the procedures of this generic interface for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  Support for Windows powershell must be added.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  Support for input wildcard patterns must be added with the help of [isFailedGlob](@ref pm_sysPath::isFailedGlob).<br>
    !>
    !>  \final{isFailedList}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedList

    module function isFailedList_SK(path, list, index, sort, showdir, showfile, showhidden, reversed, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedList_SK
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)                    :: path
        character(:,SKG)    , intent(out)   , allocatable   :: list
        integer(IK)         , intent(out)   , allocatable   :: index(:,:)
        character(*, SK)    , intent(inout) , optional      :: errmsg
        character(*, SK)    , intent(in)    , optional      :: sort
        logical(LK)         , intent(in)    , optional      :: showdir, showfile, showhidden, reversed
        logical(LK)                                         :: failed
    end function

    module function isFailedList_BSSK(path, list, sort, showdir, showfile, showhidden, reversed, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedList_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        character(*,SKG)    , intent(in)                    :: path
        type(css_type)      , intent(out)   , allocatable   :: list(:)
        character(*, SK)    , intent(inout) , optional      :: errmsg
        character(*, SK)    , intent(in)    , optional      :: sort
        logical(LK)         , intent(in)    , optional      :: showdir, showfile, showhidden, reversed
        logical(LK)                                         :: failed
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a POSIX standard (Unix-like) version of the input POSIX-style-separated `path`
    !>  where all the *relevant* POSIX reserved characters are properly escaped via (`\`).<br>
    !>
    !>  \details
    !>  For more details on the implementation of the procedure, see the documentation of [setPathPosixEscaped](@ref pm_sysPath::setPathPosixEscaped).<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of default kind \SK containing a
    !>                          POSIX-style-separated path whose relevant escape characters shall be escaped.<br>
    !>
    !>  \return
    !>  `pathEscaped`       :   The output `allocatable` scalar of type `character` of default kind \SK containing
    !>                          the modified input path where the relevant escape characters are escaped properly.<br>
    !>
    !>  \interface{getPathPosixEscaped}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathPosixEscaped
    !>      character(:, SK), allocatable :: pathEscaped
    !>
    !>      pathEscaped = getPathPosixEscaped(path)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [setPathPosixEscaped](@ref pm_sysPath::setPathPosixEscaped)<br>
    !>
    !>  \example{getPathPosixEscaped}
    !>  \include{lineno} example/pm_sysPath/getPathPosixEscaped/main.F90
    !>  \compilef{getPathPosixEscaped}
    !>  \output{getPathPosixEscaped}
    !>  \include{lineno} example/pm_sysPath/getPathPosixEscaped/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \pmed The implementation of this procedure can be improved by
    !>  replacing the call to the subroutine version with the actual code.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathPosixEscaped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPathPosixEscaped_D0_SK5(path) result(pathEscaped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathPosixEscaped_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: path
        character(:,SKG)            , allocatable               :: pathEscaped
    end function
#endif

#if SK4_ENABLED
    PURE module function getPathPosixEscaped_D0_SK4(path) result(pathEscaped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathPosixEscaped_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: path
        character(:,SKG)            , allocatable               :: pathEscaped
    end function
#endif

#if SK3_ENABLED
    PURE module function getPathPosixEscaped_D0_SK3(path) result(pathEscaped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathPosixEscaped_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: path
        character(:,SKG)            , allocatable               :: pathEscaped
    end function
#endif

#if SK2_ENABLED
    PURE module function getPathPosixEscaped_D0_SK2(path) result(pathEscaped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathPosixEscaped_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: path
        character(:,SKG)            , allocatable               :: pathEscaped
    end function
#endif

#if SK1_ENABLED
    PURE module function getPathPosixEscaped_D0_SK1(path) result(pathEscaped)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathPosixEscaped_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: path
        character(:,SKG)            , allocatable               :: pathEscaped
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a POSIX standard (Unix-like) version of the input POSIX-style-separated `path`
    !>  where all the *relevant* POSIX reserved characters are properly escaped via (`\`).<br>
    !>
    !>  \details
    !>  The POSIX reserved characters that handled by this procedure are specified in [POSIX_RESERVED_STR](@ref pm_sysPath::POSIX_RESERVED_STR).<br>
    !>  This includes escaping any single or double quotation marks that might be present in the input `path`.<br>
    !>
    !>  \param[out] pathEscaped :   The output scalar `allocatable` `character` of default kind \SK containing the
    !>                              input `path` that is reallocated to contain the input `path` where all instances
    !>                              of POSIX Shell escape character sets ([POSIX_RESERVED_STR](@ref pm_sysPath::POSIX_RESERVED_STR) or
    !>                              [POSIX_RESERVED_DQUOTE_STR](@ref pm_sysPath::POSIX_RESERVED_DQUOTE_STR)) are properly escaped via `\`.<br>
    !>  \param[in]  path        :   The input scalar `character` of default kind \SK containing a POSIX-style string whose characters shall be escaped.<br>
    !>
    !>  \interface{setPathPosixEscaped}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: setPathPosixEscaped
    !>      character(:, SK), allocatable :: pathEscaped, path
    !>
    !>      call setPathPosixEscaped(pathEscaped, path)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathPosixEscaped](@ref pm_sysPath::getPathPosixEscaped)<br>
    !>
    !>  \example{setPathPosixEscaped}
    !>  \include{lineno} example/pm_sysPath/setPathPosixEscaped/main.F90
    !>  \compilef{setPathPosixEscaped}
    !>  \output{setPathPosixEscaped}
    !>  \include{lineno} example/pm_sysPath/setPathPosixEscaped/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \plow This procedure can be converted to a generic interface with an input
    !>  `intent(inout), allocatable :: path` interface for further convenience.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setPathPosixEscaped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPathPosixEscaped_D0_SK5(pathEscaped, path)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathPosixEscaped_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(out)   , allocatable   :: pathEscaped
        character(*,SKG)            , intent(in)                    :: path
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPathPosixEscaped_D0_SK4(pathEscaped, path)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathPosixEscaped_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(out)   , allocatable   :: pathEscaped
        character(*,SKG)            , intent(in)                    :: path
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPathPosixEscaped_D0_SK3(pathEscaped, path)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathPosixEscaped_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(out)   , allocatable   :: pathEscaped
        character(*,SKG)            , intent(in)                    :: path
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPathPosixEscaped_D0_SK2(pathEscaped, path)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathPosixEscaped_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(out)   , allocatable   :: pathEscaped
        character(*,SKG)            , intent(in)                    :: path
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPathPosixEscaped_D0_SK1(pathEscaped, path)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathPosixEscaped_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(out)   , allocatable   :: pathEscaped
        character(*,SKG)            , intent(in)                    :: path
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a POSIX-standard (Unix-like) path from the input Windows-style path.<br>
    !>
    !>  \details
    !>  This procedure calls [setPathPosix](@ref pm_sysPath::setPathPosix) to generate the Windows-style normalized path.<br>
    !>  See the documentation of [setPathPosix](@ref pm_sysPath::setPathPosix) for details of the path normalization.<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of default kind \SK containing
    !>                          the Windows-style path to be converted to the POSIX-style.<br>
    !>  \param[in]  ignore  :   The input scalar `character` of the same kind as `path` of arbitrary length type parameter.<br>
    !>                          If present, any pattern matching the value of `ignore` from left to right will be interpreted verbatim and copied to the output path **as is**.<br>
    !>                          (**optional**, default = `""`)
    !>
    !>  \return
    !>  `pathPosix`         :   The output `allocatable` scalar `character` of default kind \SK
    !>                          containing the input path converted to the POSIX-style.<br>
    !>
    !>  \interface{getPathPosix}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathPosix
    !>      character(:, SK), allocatable :: pathPosix
    !>
    !>      pathPosix = getPathPosix(path, ignore = ignore)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All warnings in the documentation of [setPathPosix](@ref pm_sysPath::setPathPosix)
    !>  are also relevant and apply to this functional procedure.<br>
    !>
    !>  \remark
    !>  The subroutine version [setPathPosix](@ref pm_sysPath::setPathPosix)
    !>  of this functional procedure is expected to be slightly faster than this function.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>
    !>  \example{getPathPosix}
    !>  \include{lineno} example/pm_sysPath/getPathPosix/main.F90
    !>  \compilef{getPathPosix}
    !>  \output{getPathPosix}
    !>  \include{lineno} example/pm_sysPath/getPathPosix/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \plow The performance of this algorithm could be likely improved by
    !>  reimplementing the function to remove the unnecessary copy `pathPosix = path`.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathPosix

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPathPosix_D0_SK5(path, ignore) result(pathPosix)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathPosix_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: path
        character(*,SKG)            , intent(in)    , optional  :: ignore
        character(:,SKG)            , allocatable               :: pathPosix
    end function
#endif

#if SK4_ENABLED
    PURE module function getPathPosix_D0_SK4(path, ignore) result(pathPosix)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathPosix_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: path
        character(*,SKG)            , intent(in)    , optional  :: ignore
        character(:,SKG)            , allocatable               :: pathPosix
    end function
#endif

#if SK3_ENABLED
    PURE module function getPathPosix_D0_SK3(path, ignore) result(pathPosix)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathPosix_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: path
        character(*,SKG)            , intent(in)    , optional  :: ignore
        character(:,SKG)            , allocatable               :: pathPosix
    end function
#endif

#if SK2_ENABLED
    PURE module function getPathPosix_D0_SK2(path, ignore) result(pathPosix)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathPosix_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: path
        character(*,SKG)            , intent(in)    , optional  :: ignore
        character(:,SKG)            , allocatable               :: pathPosix
    end function
#endif

#if SK1_ENABLED
    PURE module function getPathPosix_D0_SK1(path, ignore) result(pathPosix)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathPosix_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: path
        character(*,SKG)            , intent(in)    , optional  :: ignore
        character(:,SKG)            , allocatable               :: pathPosix
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a POSIX-standard (Unix-like) path from the input (Windows-style) path.<br>
    !>
    !>  \details
    !>  See the note below on how an input POSIX-style pathis handled.<br>
    !>  The path nomalization by this procedure is specifically done by<br>
    !>  -#  replacing all instances of the Windows-style path-separator (`\`) in the input path with the POSIX-style path-separator (`/`), and
    !>  -#  collapsing all redundant directory separators such that `\A//B\/`, `///A/B\\`, `\A//B\` all become `/A/B`.<br>
    !>  -#  If the path begins with two directory separators, it is interpreted as a [UNC path](@ref pm_sysPath::getPathHostNameIndex) and the initial hostname portion of the path remains untouched.<br>
    !>  -#  Additionally, if `ignore` is present, then all instances of the pattern `ignore` will be copied to the output `path` **as is**.<br>
    !>      For example,
    !>      -#  if the input path contains a Hex value `\xhh`, then set `ignore = \xhh` or
    !>      -#  if the is already (partially) POSIX-compliant such that `\` should be interprested verbatim then ignore can be set to `\`.<br>
    !>
    !>  \param[inout]   path    :   The input/output `allocatable` scalar `character` of default kind \SK.<br>
    !>                              On input, it contains the Windows-style path to be converted to the POSIX-style.<br>
    !>                              On output, it is reallocated to contain the converted POSIX-style path.<br>
    !>                              On output, all instances of the backslash are replaced with forward-slash.<br>
    !>  \param[in]      ignore  :   The input scalar `character` of the same kind as `path` of arbitrary length type parameter.<br>
    !>                              If present, any pattern matching the value of `ignore` from left to right will be interpreted verbatim and copied to the output path **as is**.<br>
    !>                              (**optional**, default = `""`)
    !>
    !>  \interface{setPathPosix}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: setPathPosix
    !>      character(:, SK), allocatable :: path
    !>
    !>      call setPathPosix(path)
    !>      call setPathPosix(path, ignore = ignore)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Any escape character `\` is interpreted as a Windows-style directory separator by the procedures
    !>  of this generic interface and will be converted to POSIX-style directory separator character `/`.<br>
    !>  Therefore, passing an already POSIX-compliant path to the procedures of this generic interface is fine as long as,
    !>  -#  **no** POSIX special character is escaped via the POSIX escape character `\` in the input `path`, and
    !>  -#  **no** backward slash `\` appears in the `path` unless it serves as a directory separator or `ignore = `"\"`.<br>
    !>
    !>  \warnpure
    !>  The impurity is caused by the dependency on [setInserted](@ref pm_arrayInsert::setInserted).<br>
    !>
    !>  \warning
    !>  Note that no trimming of the input `path` is performed by the procedure.<br>
    !>  Consequently, any leading or trailing whitespace characters will be significant.<br>
    !>
    !>  \remark
    !>  This procedure has the same functionality as its functional equivalent [getPathPosix](@ref pm_sysPath::getPathPosix)
    !>  except for the fact it can be slightly faster than the functional procedure.<br>
    !>
    !>  \note
    !>  The procedure [setPathPosixEscaped](@ref pm_sysPath::setPathPosixEscaped) can be called if any Bash environment special characters as
    !>  specified in [POSIX_RESERVED_STR](@ref pm_sysPath::POSIX_RESERVED_STR) need to be escaped in the output `path` via the escape character `\`.<br>
    !>  Otherwise, [getPathVerbatimPosix](@ref getPathVerbatimPosix) can generate a path whose characters will be interpreted verbatim in POSIX shell environment.<br>
    !>
    !>  \see
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [getPathPosixEscaped](@ref pm_sysPath::getPathPosixEscaped)<br>
    !>  [setPathPosixEscaped](@ref pm_sysPath::setPathPosixEscaped)<br>
    !>
    !>  \example{setPathPosix}
    !>  \include{lineno} example/pm_sysPath/setPathPosix/main.F90
    !>  \compilef{setPathPosix}
    !>  \output{setPathPosix}
    !>  \include{lineno} example/pm_sysPath/setPathPosix/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  An optional argument external procedure `iseq()` should be added to the procedures
    !>  of this generic interface to bring flexibility to searching for `ignore` patterns.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setPathPosix

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setPathPosix_D0_SK5(path, ignore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathPosix_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: ignore
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setPathPosix_D0_SK4(path, ignore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathPosix_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: ignore
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setPathPosix_D0_SK3(path, ignore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathPosix_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: ignore
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setPathPosix_D0_SK2(path, ignore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathPosix_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: ignore
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setPathPosix_D0_SK1(path, ignore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathPosix_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: ignore
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a normalized Windows-style path from the input POSIX-style or Windows-style path,
    !>  such that it conforms to the path conventions of the Windows Command Prompt (**CMD**).<br>
    !>
    !>  \details
    !>  This procedure calls [setPathWindows](@ref pm_sysPath::setPathWindows) to generate the Windows-style normalized path.<br>
    !>  See the documentation of [setPathWindows](@ref pm_sysPath::setPathWindows) for details of the path normalization.<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of kind \SKALL containing
    !>                          the POSIX-style path to be converted to the Windows-style.<br>
    !>  \param[in]  ignore  :   The input scalar `character` of the same kind as `path` of arbitrary length type parameter.<br>
    !>                          If present, any pattern matching the value of `ignore` from left to right will be interpreted verbatim and copied to the output path **as is**.<br>
    !>                          (**optional**, default = `""`)
    !>
    !>  \return
    !>  `pathWindows`       :   The output `allocatable` scalar of type `character` of the same kind as the input `path`
    !>                          containing the modified input path such that it conforms to the Windows Command Prompt style.<br>
    !>
    !>  \interface{getPathWindows}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathWindows
    !>
    !>      pathWindows = getPathWindows(path)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The procedures under this generic interface assume that the input path is valid in both Unix and Windows environments.<br>
    !>  Specifically, the procedures do not make any attempts to remove or modify the appearance of any Windows reserved characters
    !>  ([WINDOWS_RESERVED_STR](@ref pm_sysPath::WINDOWS_RESERVED_STR)) or device names ([WINDOWS_RESERVED_DEVICE_NAME](@ref pm_sysPath::WINDOWS_RESERVED_DEVICE_NAME)).<br>
    !>
    !>  \warning
    !>  Although the procedures under this generic interface can properly handle input paths with trailing whitespace character,
    !>  you should never create a directory or filename with a trailing space.<br>
    !>  Trailing spaces can make it difficult or impossible to access a directory, and applications
    !>  commonly fail when attempting to handle directories or files whose names include trailing spaces.<br>
    !>
    !>  \warning
    !>  Note that **no** trimming of the input `path` is performed by the procedure.<br>
    !>  Consequently, any leading or trailing whitespace characters will be significant.<br>
    !>  Similarly, single and double quotation marks that appear anywhere in path are significant and kept.<br>
    !>
    !>  \pure
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{getPathWindows}
    !>  \include{lineno} example/pm_sysPath/getPathWindows/main.F90
    !>  \compilef{getPathWindows}
    !>  \output{getPathWindows}
    !>  \include{lineno} example/pm_sysPath/getPathWindows/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathWindows

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getPathWindows_D0_SK5(path, ignore) result(pathWindows)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathWindows_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: path
        character(*,SKG)            , intent(in)    , optional  :: ignore
        character(:,SKG)            , allocatable               :: pathWindows
    end function
#endif

#if SK4_ENABLED
    pure module function getPathWindows_D0_SK4(path, ignore) result(pathWindows)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathWindows_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: path
        character(*,SKG)            , intent(in)    , optional  :: ignore
        character(:,SKG)            , allocatable               :: pathWindows
    end function
#endif

#if SK3_ENABLED
    pure module function getPathWindows_D0_SK3(path, ignore) result(pathWindows)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathWindows_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: path
        character(*,SKG)            , intent(in)    , optional  :: ignore
        character(:,SKG)            , allocatable               :: pathWindows
    end function
#endif

#if SK2_ENABLED
    pure module function getPathWindows_D0_SK2(path, ignore) result(pathWindows)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathWindows_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: path
        character(*,SKG)            , intent(in)    , optional  :: ignore
        character(:,SKG)            , allocatable               :: pathWindows
    end function
#endif

#if SK1_ENABLED
    pure module function getPathWindows_D0_SK1(path, ignore) result(pathWindows)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathWindows_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: path
        character(*,SKG)            , intent(in)    , optional  :: ignore
        character(:,SKG)            , allocatable               :: pathWindows
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Normalize the input POSIX-style or Windows-style path to a Windows-style path,
    !>  such that it conforms to the path conventions of the Windows Command Prompt (**CMD**).<br>
    !>
    !>  \details
    !>  The nomalization done by this procedure is only part of the Windows system-level path normalization that is detailed in the documentation of [pm_sysPath](@ref pm_sysPath).<br>
    !>  Here it is specifically done by<br>
    !>  -#  replacing all instances of the POSIX-style path-separator (`/`) in the input path with the Windows-style path-separator (`\`), and
    !>  -#  collapsing all redundant directory separators such that `\A//B\/`, `///A/B\\`, `\A//B\` all become `/A/B`.<br>
    !>  -#  If the path begins with two directory separators, it is interpreted as a [UNC path](@ref pm_sysPath::getPathHostNameIndex) and the initial hostname portion of the path remains untouched.<br>
    !>  -#  Additionally, if `ignore` is present, then all instances of the pattern `ignore` will be copied to the output `path` **as is**.<br>
    !>      For example,
    !>      -#  if the forward slashes `/` should be kept as is, then set `ignore = "/"`.<br>
    !>
    !>  \param[inout]   path    :   The input `allocatable` scalar `character` of kind \SKALL containing the POSIX-style path to be converted to the Windows-style.<br>
    !>                              On output, it contains the modified input path such that it conforms to the Windows-Command-Prompt style.<br>
    !>  \param[in]      ignore  :   The input scalar `character` of the same kind as `path` of arbitrary length type parameter.<br>
    !>                              If present, any pattern matching the value of `ignore` from left to right will be interpreted verbatim and copied to the output path **as is**.<br>
    !>                              (**optional**, default = `""`)
    !>
    !>  \interface{setPathWindows}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: setPathWindows
    !>
    !>      call setPathWindows(path, ignore = ignore)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The procedures under this generic interface assume that the input path is valid in both Unix and Windows environments.<br>
    !>  Specifically, the procedures do not make any attempts to remove or modify the appearance of any Windows reserved characters
    !>  ([WINDOWS_RESERVED_STR](@ref pm_sysPath::WINDOWS_RESERVED_STR)) or device names ([WINDOWS_RESERVED_DEVICE_NAME](@ref pm_sysPath::WINDOWS_RESERVED_DEVICE_NAME)).<br>
    !>
    !>  \warning
    !>  Although the procedures under this generic interface can properly handle input paths with trailing whitespace character, you
    !>  should never create a directory or filename with a trailing space. Trailing spaces can make it difficult or impossible to access
    !>  a directory, and applications commonly fail when attempting to handle directories or files whose names include trailing spaces.<br>
    !>
    !>  \warning
    !>  Note that **no** trimming of the input `path` is performed by the procedure.<br>
    !>  Consequently, any leading or trailing whitespace characters will be significant.<br>
    !>  Similarly, single and double quotation marks that appear anywhere in path are significant and kept.<br>
    !>
    !>  \pure
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{setPathWindows}
    !>  \include{lineno} example/pm_sysPath/setPathWindows/main.F90
    !>  \compilef{setPathWindows}
    !>  \output{setPathWindows}
    !>  \include{lineno} example/pm_sysPath/setPathWindows/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  An optional argument external procedure `iseq()` should be added to the procedures
    !>  of this generic interface to bring flexibility to searching for `ignore` patterns.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setPathWindows

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setPathWindows_D0_SK5(path, ignore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathWindows_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)            , intent(inout) , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: ignore
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setPathWindows_D0_SK4(path, ignore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathWindows_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)            , intent(inout) , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: ignore
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setPathWindows_D0_SK3(path, ignore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathWindows_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)            , intent(inout) , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: ignore
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setPathWindows_D0_SK2(path, ignore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathWindows_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)            , intent(inout) , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: ignore
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setPathWindows_D0_SK1(path, ignore)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathWindows_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)            , intent(inout) , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: ignore
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the index of the last character of the hostname part of the input UNC path.<br>
    !>
    !>  \details
    !>  -#  The Universal Naming Convention (**UNC**) is a standard for identifying servers, printers and other resources in a network.<br>
    !>  -#  According to the US patent 5341499 (by IBM), the UNC was jointly invented by IBM and Microsoft for use in their jointly developed Local Area Network (LAN) software products.<br>
    !>  -#  A UNC path uses **double** slashes or backslashes to precede the name of the computer.<br>
    !>  -#  The host-name portion of a UNC path can consist of either a network name string set by an
    !>      administrator and maintained by a network naming service like DNS or WINS, or by an IP address.<br>
    !>  -#  These hostnames normally refer to either a PC, printer, or other devices.<br>
    !>  -#  The path (disk and directories) within the computer are separated with a single slash or backslash, as usual.<br>
    !>  -#  Note that in under Microsoft DOS/Windows platforms, drive letters (`c:`, `d:`, etc.) are not used in UNC names and if used, the `:` character is dropped.<br>
    !>  -#  The typical structure of a UNC path is the following,
    !>      -#  On Windows systems: `\\hostname\path` or `//hostname/path` has the hostname `\\hostname`.<br>
    !>      -#  On Unix systems: `//hostname/path` has the hostname `//hostname`.<br>
    !>  -#  Note that UNC path names are not universally recognized on all Unix platforms.<br>
    !>      For more information, see the [POSIX standard](https://pubs.opengroup.org/onlinepubs/009695399/xrat/xbd_chap04.html#tag_01_04_11).<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of kind \SKALL of arbitrary length type parameter containing
    !>                          the potentially UNC path whose hostname portion is to be identified and the index of its last character returned.<br>
    !>  \param[in]  dirsep  :   The input scalar `character` of the same kind as `path`, of arbitrary length type parameter
    !>                          containing the directory separators that might be present in the input `path`. For example,
    !>                          -#  On Windows platforms, this can be [DIR_SEP_WINDOWS_ALL](@ref pm_sysPath::DIR_SEP_WINDOWS_ALL).<br>
    !>                          -#  On Unix platforms, this can be [DIR_SEP_POSIX_ALL](@ref pm_sysPath::DIR_SEP_POSIX_ALL).<br>
    !>
    !>  \return
    !>  `index`             :   The output scalar `integer` of default kind \IK containing the
    !>                          index of the last character of hostname portion of the input UNC `path`.<br>
    !>                          If there is no hostname, `index = 0` on output.<br>
    !>
    !>  \interface{getPathHostNameIndex}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathHostNameIndex
    !>      use pm_kind, only: IK
    !>      integer(IK) :: index
    !>
    !>      index = getPathHostNameIndex(path, dirsep)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that **no** trimming of the input `path` is performed by the procedure.<br>
    !>  Consequently, any leading or trailing whitespace characters will be significant.<br>
    !>  Similarly, single and double quotation marks that appear anywhere in path are significant and kept.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{getPathHostNameIndex}
    !>  \include{lineno} example/pm_sysPath/getPathHostNameIndex/main.F90
    !>  \compilef{getPathHostNameIndex}
    !>  \output{getPathHostNameIndex}
    !>  \include{lineno} example/pm_sysPath/getPathHostNameIndex/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathHostNameIndex

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getPathHostNameIndex_SK5(path, dirsep) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathHostNameIndex_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: path, dirsep
        integer(IK)                                             :: index
    end function
#endif

#if SK4_ENABLED
    pure module function getPathHostNameIndex_SK4(path, dirsep) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathHostNameIndex_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: path, dirsep
        integer(IK)                                             :: index
    end function
#endif

#if SK3_ENABLED
    pure module function getPathHostNameIndex_SK3(path, dirsep) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathHostNameIndex_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: path, dirsep
        integer(IK)                                             :: index
    end function
#endif

#if SK2_ENABLED
    pure module function getPathHostNameIndex_SK2(path, dirsep) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathHostNameIndex_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: path, dirsep
        integer(IK)                                             :: index
    end function
#endif

#if SK1_ENABLED
    pure module function getPathHostNameIndex_SK1(path, dirsep) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathHostNameIndex_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: path, dirsep
        integer(IK)                                             :: index
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` is the input Windows-style path **begins** with a drive-letter pattern like `C:`.<br>
    !>
    !>  \details
    !>  A drive letter can appear as the prefix in both absolute (e.g., `C:\` or `C:/`) or relative (e.g., `C:paramonte` or `C:../`) Windows paths.<br>
    !>  When it appears in relative paths, its proper treatment, particularly when joining path segments becomes critical.<br>
    !>  The procedures of this generic interface generate and return `.true.` if the input `path` begins with an Enlish alphabet,
    !>  followed by a colon `:` character. Otherwise, they return `.false.`.<br>
    !>
    !>  \param[in]  path        :   The input scalar `character` of kind \SKALL containing an existing or virtual path.<br>
    !>
    !>  \return
    !>  `pathHasDriveLetter`    :   The output scalar `logical` of default kind \LK that is `.true.` <b>if and only if</b> the input path is an absolute path based on Windows conventions.<br>
    !>
    !>  \interface{hasDriveLetter}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only : LK
    !>      use pm_sysPath, only: hasDriveLetter
    !>      logical(LK) :: pathHasDriveLetter
    !>
    !>      pathHasDriveLetter = hasDriveLetter(path)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPathAbs](@ref pm_sysPath::getPathAbs)<br>
    !>  [isPathAbsPosix](@ref pm_sysPath::isPathAbsPosix)<br>
    !>  [isPathAbsWindows](@ref pm_sysPath::isPathAbsWindows)<br>
    !>
    !>  \example{hasDriveLetter}
    !>  \include{lineno} example/pm_sysPath/hasDriveLetter/main.F90
    !>  \compilef{hasDriveLetter}
    !>  \output{hasDriveLetter}
    !>  \include{lineno} example/pm_sysPath/hasDriveLetter/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final{hasDriveLetter}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface hasDriveLetter

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function hasDriveLetter_SK5(path) result(pathHasDriveLetter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hasDriveLetter_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathHasDriveLetter
    end function
#endif

#if SK4_ENABLED
    pure elemental module function hasDriveLetter_SK4(path) result(pathHasDriveLetter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hasDriveLetter_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathHasDriveLetter
    end function
#endif

#if SK3_ENABLED
    pure elemental module function hasDriveLetter_SK3(path) result(pathHasDriveLetter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hasDriveLetter_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathHasDriveLetter
    end function
#endif

#if SK2_ENABLED
    pure elemental module function hasDriveLetter_SK2(path) result(pathHasDriveLetter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hasDriveLetter_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathHasDriveLetter
    end function
#endif

#if SK1_ENABLED
    pure elemental module function hasDriveLetter_SK1(path) result(pathHasDriveLetter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hasDriveLetter_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathHasDriveLetter
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` is the input path is an absolute Windows path (whether existing or virtual), otherwise return `.false.`.<br>
    !>
    !>  \details
    !>  The definition of an absolute path on Windows is given in the [details section](@ref syspath_pmod_details_windows) of [pm_sysPath](@ref pm_sysPath).<br>
    !>  In general, a Windows path is considered absolute by the procedures of this generic interface if,
    !>  -#  it begins with two directory separators recognized on Windows as stored in [DIR_SEP_WINDOWS_ALL](@ref pm_sysPath::DIR_SEP_WINDOWS_ALL) (i.e., it is a UNC path), or
    !>  -#  it begins with a drive letter, followed by a colon character `:`, followed by at least one of the directory separators recognized on Windows as stored in
    !>      [DIR_SEP_WINDOWS_ALL](@ref pm_sysPath::DIR_SEP_WINDOWS_ALL).<br>
    !>
    !>  A path that begins with one of the directory separators recognized on Windows as stored in [DIR_SEP_WINDOWS_ALL](@ref pm_sysPath::DIR_SEP_WINDOWS_ALL), **is not absolute**.<br>
    !>  Such a path is recognied by Microsoft as a relative path with respect to the current drive.<br>
    !>  This is contrary to the conventions followed by some other languages (e.g., Python).<br>
    !>  This generic interface follows Microsoft conventions for the definition of an absolute path on Windows.<br>
    !>  For example,
    !>  -#  `\paramonte` is a relative path with respect to the current drive of the runtime shell.<br>
    !>  -#  `C:paramonte` is a relative path on drive `C` with respect to the current working directory of the runtime shell.<br>
    !>
    !>  By definition, an absolute Windows path must contain at least `3` characters (or more for UNC-style absolute paths).<br>
    !>  To check whether a path begins with a drive letter, use [hasDriveLetter](@ref pm_sysPath::hasDriveLetter).<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of kind \SKALL containing an existing or virtual path.<br>
    !>                          An empty input path is by definition a relative path (**not** absolute).<br>
    !>
    !>  \return
    !>  `pathIsAbs`         :   The output scalar `logical` of default kind \LK that is `.true.` <b>if and only if</b> the input path is an absolute path based on Windows conventions.<br>
    !>
    !>  \interface{isPathAbsWindows}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only : LK
    !>      use pm_sysPath, only: isPathAbsWindows
    !>      logical(LK) :: pathIsAbs
    !>
    !>      pathIsAbs = isPathAbsWindows(path)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPathAbs](@ref pm_sysPath::getPathAbs)<br>
    !>  [isPathAbsPosix](@ref pm_sysPath::isPathAbsPosix)<br>
    !>  [isPathAbsWindows](@ref pm_sysPath::isPathAbsWindows)<br>
    !>
    !>  \example{isPathAbsWindows}
    !>  \include{lineno} example/pm_sysPath/isPathAbsWindows/main.F90
    !>  \compilef{isPathAbsWindows}
    !>  \output{isPathAbsWindows}
    !>  \include{lineno} example/pm_sysPath/isPathAbsWindows/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final{isPathAbsWindows}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isPathAbsWindows

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isPathAbsWindows_SK5(path) result(pathIsAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isPathAbsWindows_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathIsAbs
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isPathAbsWindows_SK4(path) result(pathIsAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isPathAbsWindows_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathIsAbs
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isPathAbsWindows_SK3(path) result(pathIsAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isPathAbsWindows_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathIsAbs
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isPathAbsWindows_SK2(path) result(pathIsAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isPathAbsWindows_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathIsAbs
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isPathAbsWindows_SK1(path) result(pathIsAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isPathAbsWindows_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathIsAbs
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` is the input path is an absolute POSIX path (whether existing or virtual), otherwise return `.false.`.<br>
    !>
    !>  \details
    !>  A POSIX path is absolute if,
    !>  -#  it begins with a POSIX directory separator as stored in [DIR_SEP_POSIX_ALL](@ref pm_sysPath::DIR_SEP_POSIX_ALL).<br>
    !>
    !>
    !>  \param[in]  path    :   The input scalar `character` of kind \SKALL containing an existing or virtual path.<br>
    !>                          An empty input path is by definition a relative path (**not** absolute).<br>
    !>
    !>  \return
    !>  `pathIsAbs`         :   The output scalar `logical` of default kind \LK that is `.true.` <b>if and only if</b> the input path is an absolute path based on POSIX conventions.<br>
    !>
    !>  \interface{isPathAbsPosix}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only : LK
    !>      use pm_sysPath, only: isPathAbsPosix
    !>      logical(LK) :: pathIsAbs
    !>
    !>      pathIsAbs = isPathAbsPosix(path)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getPathAbs](@ref pm_sysPath::getPathAbs)<br>
    !>  [isPathAbsPosix](@ref pm_sysPath::isPathAbsPosix)<br>
    !>  [isPathAbsWindows](@ref pm_sysPath::isPathAbsWindows)<br>
    !>
    !>  \example{isPathAbsPosix}
    !>  \include{lineno} example/pm_sysPath/isPathAbsPosix/main.F90
    !>  \compilef{isPathAbsPosix}
    !>  \output{isPathAbsPosix}
    !>  \include{lineno} example/pm_sysPath/isPathAbsPosix/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final{isPathAbsPosix}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isPathAbsPosix

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isPathAbsPosix_SK5(path) result(pathIsAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isPathAbsPosix_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathIsAbs
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isPathAbsPosix_SK4(path) result(pathIsAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isPathAbsPosix_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathIsAbs
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isPathAbsPosix_SK3(path) result(pathIsAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isPathAbsPosix_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathIsAbs
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isPathAbsPosix_SK2(path) result(pathIsAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isPathAbsPosix_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathIsAbs
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isPathAbsPosix_SK1(path) result(pathIsAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isPathAbsPosix_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)    , intent(in)                :: path
        logical(LK)                                     :: pathIsAbs
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the absolute path corresponding to the input (potentially *relative*) path.<br>
    !>
    !>  \details
    !>  The procedures of this generic interface take the following steps to generate the absolute path corresponding to the input `path`,
    !>  <ol>
    !>      <li>    If the input `path` is empty, the absolute path to the current directory will be returned.<br>
    !>      <li>    If the input `path` is an absolute path determined via
    !>              [isPathAbsPosix](@ref pm_sysPath::isPathAbsPosix) or [isPathAbsWindows](@ref pm_sysPath::isPathAbsWindows),
    !>              the input `path` is returned as the output absolute path.<br>
    !>      <li>    If the input `path` is relative, then,
    !>              <ol>
    !>                  <li>    If the library is compiled with Intel compilers, the absolute path is determined by a
    !>                          call to `FullPathQQ()` procedure of `ifport` module of the Intel `ifort` compiler.<br>
    !>                  <li>    If the above step fails or if the compiler is non-intel,
    !>                          <ol>
    !>                              <li>    The absolute path is directly inferred by a system call to the relevant
    !>                                      Microsoft or Windows PowerShell command if the shell is available on the platform.<br>
    !>                          </ol>
    !>              </ol>
    !>      <li>    If all of the above approaches fail or do not apply (particularly the case with POSIX-compliant runtime shells),
    !>              the absolute path is constructed by appending the current working directory to the partially qualified input `path`.<br>
    !>  </ol>
    !>
    !>  See the warnings and remarks below.<br>
    !>
    !>  \param[in]      path    :   The input scalar `character` of kind \SK containing the relative path to be resolved to the corresponding absolute path.<br>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK.<br>
    !>                              It is `.true.` <b>if and only if</b> an error occurs while resolving `path`.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                              A length of `2047` characters for `errmsg` is likely enough to capture most error messages in full.<br>
    !>                              (**optional**. It can be present **only if** `failed` is also present.)
    !>
    !>  \return
    !>  `pathAbs`               :   The output `allocatable` scalar `character` of default kind \SK containing the absolute path corresponding to the input (potentially relative) `path`.<br>
    !>                              <ol>
    !>                                  <li>    If the input `path` is relative but the path resolution succeeds, the program sets `pathAbs` to the resolved path and `failed = .false.` upon return.<br>
    !>                                  <li>    If the input `path` is relative but the path resolution fails, the program sets `pathAbs = path` and `failed = .true.` upon return.<br>
    !>                                  <li>    If the input `path` is absolute already, the program sets `pathAbs = path` and `failed = .false.` upon return.<br>
    !>                              </ol>
    !>
    !>  \interface{getPathAbs}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK
    !>      use pm_sysPath, only: getPathAbs
    !>      character(:, SK), allocatable :: pathAbs
    !>
    !>      pathAbs = getPathAbs(path)
    !>      pathAbs = getPathAbs(path, failed)
    !>      pathAbs = getPathAbs(path, failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  The procedures of this generic interface do not make any attempts to resolve any soft links (symlinks) present in the input `path`.<br>
    !>  Consequently, any instances of parent directory pattern (`/../` or `\..\`) that may appear within the input `path` remain untouched
    !>  because resolving such patterns can inadvertently modify symlinks if any are present.<br>
    !>  In most applications, there is no need to resolve soft links or instances of parent directory
    !   patterns because the operating system automatically and gracefully handles such cases.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isPathAbsWindows](@ref pm_sysPath::isPathAbsWindows)<br>
    !>  [isPathAbsPosix](@ref pm_sysPath::isPathAbsPosix)<br>
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{getPathAbs}
    !>  \include{lineno} example/pm_sysPath/getPathAbs/main.F90
    !>  \compilef{getPathAbs}
    !>  \output{getPathAbs}
    !>  \include{lineno} example/pm_sysPath/getPathAbs/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \pvlow
    !>  A subroutine version of this functional interface could be implemented in the future to avoid allocations.<br>
    !>
    !>  \final{getPathAbs}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathAbs

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getPathAbs(path) result(pathAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathAbs
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    :: path
        character(:,SKG), allocatable   :: pathAbs
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getPathAbsFailed(path, failed) result(pathAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathAbsFailed
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    :: path
        logical(LK)     , intent(out)   :: failed
        character(:,SKG), allocatable   :: pathAbs
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getPathAbsFailedMsg(path, failed, errmsg) result(pathAbs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathAbsFailedMsg
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    :: path
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        character(:,SKG), allocatable   :: pathAbs
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **Current Working Directory (CWD)** of the runtime shell.<br>
    !>
    !>  \details
    !>  The name of this interface is intentionally chosen to be different from `getcwd`
    !>  which is an Intel `IFPORT` and GNU Fortran extension procedure name.<br>
    !>
    !>  The implementation of this procedure is compiler-dependent depending on whether the output argument `errmsg` is present or missing.<br>
    !>  <ol>
    !>      <li> If `errmsg` is missing, then the appropriate runtime system shell command `pwd` is called to generate the CWD.<br>
    !>      <li> If `errmsg` is present, then,<br>
    !>      <ol>
    !>          <li>    If the library is built with either GNU or Intel Fortran compilers, then the compiler intrinsic procedure `getcwd()` is called to generate the CWD.<br>
    !>          <li>    If the library is built with other compilers, then the appropriate runtime system shell command `pwd` is called to generate the CWD.<br>
    !>      </ol>
    !>      This behavior is enforced by the Intel and GNU interfaces for `getcwd()` which does not return error message for failures.<br>
    !>  </ol>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK.<br>
    !>                              It is `.true.` <b>if and only if</b> an error occurs while inferring the current working directory.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                              A length of `2047` characters for `errmsg` is likely enough to capture most error messages in full.<br>
    !>                              (**optional**. It can be present **only if** `failed` is also present.)
    !>
    !>  \return
    !>  `dirCurrent`            :   The output `allocatable` scalar `character` of default kind \SK containing the path to the current working directory (including drive letter on windows OS).<br>
    !>                              The default value for `dirCurrent` is `"."` if the procedure fails to infer the full path to the current working directory.<br>
    !>
    !>  \interface{getDirCurrent}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getDirCurrent
    !>      use pm_kind, only: LK
    !>      character(:, SK), allocatable :: dirCurrent
    !>      character(2047, SK) :: errmsg
    !>      logical(LK) :: failed
    !>
    !>      dirCurrent = getDirCurrent()
    !>      dirCurrent = getDirCurrent(failed)
    !>      dirCurrent = getDirCurrent(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [isFailedGetDirTemp](@ref pm_sysShell::isFailedGetDirTemp)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{getDirCurrent}
    !>  \include{lineno} example/pm_sysPath/getDirCurrent/main.F90
    !>  \compilef{getDirCurrent}
    !>  \output{getDirCurrent}
    !>  \include{lineno} example/pm_sysPath/getDirCurrent/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \plow
    !>  A subroutine version of this functional interface could be implemented in
    !>  future to avoid allocations and allow for non-default character kinds.<br>
    !>
    !>  \final{getDirCurrent}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDirCurrent

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getDirCurrent() result(dirCurrent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirCurrent
#endif
        use pm_kind, only: SKG => SK
        character(:,SKG), allocatable   :: dirCurrent
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getDirCurrentFailed(failed) result(dirCurrent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirCurrentFailed
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        character(:,SKG), allocatable   :: dirCurrent
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getDirCurrentFailedMsg(failed, errmsg) result(dirCurrent)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirCurrentFailedMsg
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        character(:,SKG), allocatable   :: dirCurrent
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Home Directory of the current user in the current runtime shell.<br>
    !>
    !>  \details
    !>  The home directory identification is platform-dependent:
    !>  <ol>
    !>      <li>    On **Unix** platforms:
    !>              <ol>
    !>                  <li>    The home directory is retrieved from the envrinment variable `HOME`.<br>
    !>                          If the attempt fails, the home directory will is looked up directly in the
    !>                          password directory.<br>
    !>              </ol>
    !>      <li>    On **Windows** platforms:
    !>              <ol>
    !>                  <li>    The home directory is retrieved from the envrinment variable `USERPROFILE` of if not available,
    !>                          from the concatenation of the contents of `HOMEDRIVE` and `HOMEPATH` envrinment variables.<br>
    !>              </ol>
    !>      <li>    If `user` is specified, the corresponding user home directory will be returned. See examples below.<br>
    !>      <li>    If an error occurs the returned home directory is an empty string.<br>
    !>  </ol>
    !>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK.<br>
    !>                              It is `.true.` <b>if and only if</b> an error occurs while inferring the current home directory.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                              A length of `2047` characters for `errmsg` is likely enough to capture most error messages in full.<br>
    !>                              (**optional**. If missing and an error occurs, then no error will be reported.)
    !>  \param[in]      user    :   The input scalar `character` of default kind \SK containing the username whose home directory is to be returned.<br>
    !>                              (**optional**. If missing, the home directory of the current user is returned.)
    !>
    !>
    !>  \return
    !>  `dirHome`               :   The output `allocatable` scalar `character` of default kind \SK containing the path to the home directory (or empty if an error occurs).<br>
    !>
    !>  \interface{getDirHome}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getDirHome
    !>      use pm_kind, only: LK
    !>      character(:, SK), allocatable :: dirHome
    !>      character(2047, SK) :: errmsg
    !>      logical(LK) :: failed
    !>
    !>      dirHome = getDirHome(user = user, failed = failed, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [isFailedMakeDir](@ref pm_sysPath::isFailedMakeDir)<br>
    !>  [isFailedGetDirTemp](@ref pm_sysShell::isFailedGetDirTemp)<br>
    !>  [isFailedMakeDirTemp](@ref pm_sysPath::isFailedMakeDirTemp)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{getDirHome}
    !>  \include{lineno} example/pm_sysPath/getDirHome/main.F90
    !>  \compilef{getDirHome}
    !>  \output{getDirHome}
    !>  \include{lineno} example/pm_sysPath/getDirHome/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \plow
    !>  A subroutine version of this functional interface could be implemented in
    !>  future to avoid allocations and allow for non-default character kinds.<br>
    !>
    !>  \final{getDirHome}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDirHome

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getDirHome(user, failed, errmsg) result(dirHome)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirHome
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    , optional      :: user
        logical(LK)     , intent(out)   , optional      :: failed
        character(*, SK), intent(inout) , optional      :: errmsg
        character(:,SKG)                , allocatable   :: dirHome
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the input path wherein the initial tilde character (`~`) is expanded to the home directory of current user.<br>
    !>
    !>  \details
    !>  On Unix and Windows, return the argument with an initial component of ~ or ~user replaced by that user’s home directory.<br>
    !>  <ol>
    !>      <li>    On **Unix** platforms:
    !>              <ol>
    !>                  <li>    If the character appearing immediately after `~` is a shell directory separator
    !>                          (i.e., forward slash on Unix platforms, backslash or forward slash on Windows)
    !>                          then `~` is replaced with the current home directory of the user. See examples below.<br>
    !>                  <li>    If the character appearing immediately after `~` is not a shell directory separator
    !>                          (i.e., forward slash on Unix platforms, backslash or forward slash on Windows)
    !>                          then `~` along with the trailing token (until before the first directory separator)
    !>                          is looked up directly in the password directory as `~user`. See examples below.<br>
    !>              </ol>
    !>      <li>    On **Windows** platforms:
    !>              <ol>
    !>                  <li>    If the character appearing immediately after `~` is a shell directory separator
    !>                          (i.e., forward slash on Unix platforms, backslash or forward slash on Windows)
    !>                          then `~` is replaced with the current home directory of the user. See examples below.<br>
    !>                  <li>    If the character appearing immediately after `~` is not a shell directory separator
    !>                          (i.e., forward slash on Unix platforms, backslash or forward slash on Windows)
    !>                          then `~` along with the trailing token (until before the first directory separator)
    !>                          is handled by stripping the last directory component from the user home path.<br>
    !>              </ol>
    !>      <li>    If an error occurs or the first character of `path` is not a `~`, then the original path is returned unchanged.<br>
    !>  </ol>
    !>
    !>  \param[in]      path    :   The input scalar `character` of default kind \SK containing the path whose leading `~` or `~user` is to be expanded.<br>
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK.<br>
    !>                              It is `.true.` <b>if and only if</b> an error occurs while expanding the user home directory.<br>
    !>                              (**optional**. If missing, the program will return the original input `path` unchanged.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                              A length of `2047` characters for `errmsg` is likely enough to capture most error messages in full.<br>
    !>                              (**optional**. It can be present **only if** `failed` is also present.)
    !>
    !>  \return
    !>  `pathExpandedUser`      :   The output `allocatable` scalar `character` of default kind \SK containing
    !>                              the input path where `~` or `~user` is expanded to the specific `user` home directory.<br>
    !>
    !>  \interface{getPathExpandedUser}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathExpandedUser
    !>      use pm_kind, only: LK
    !>      character(:, SK), allocatable :: pathExpandedUser
    !>      character(2047, SK) :: errmsg
    !>      logical(LK) :: failed
    !>
    !>      pathExpandedUser = getPathExpandedUser()
    !>      pathExpandedUser = getPathExpandedUser(failed)
    !>      pathExpandedUser = getPathExpandedUser(failed, errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [isFailedGetDirTemp](@ref pm_sysShell::isFailedGetDirTemp)<br>
    !>  [isFailedMakeDirTemp](@ref pm_sysPath::isFailedMakeDirTemp)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{getPathExpandedUser}
    !>  \include{lineno} example/pm_sysPath/getPathExpandedUser/main.F90
    !>  \compilef{getPathExpandedUser}
    !>  \output{getPathExpandedUser}
    !>  \include{lineno} example/pm_sysPath/getPathExpandedUser/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \plow
    !>  A subroutine version of this functional interface could be implemented in
    !>  future to avoid allocations and allow for non-default character kinds.<br>
    !>
    !>  \final{getPathExpandedUser}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathExpandedUser

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getPathExpandedUser(path) result(pathExpandedUser)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathExpandedUser
#endif
        use pm_kind, only: SKG => SK
        character(*, SK), intent(in)    :: path
        character(:,SKG), allocatable   :: pathExpandedUser
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getPathExpandedUserFailed(path, failed) result(pathExpandedUser)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathExpandedUserFailed
#endif
        use pm_kind, only: SKG => SK
        character(*, SK), intent(in)    :: path
        logical(LK)     , intent(out)   :: failed
        character(:,SKG), allocatable   :: pathExpandedUser
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getPathExpandedUserFailedMsg(path, failed, errmsg) result(pathExpandedUser)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathExpandedUserFailedMsg
#endif
        use pm_kind, only: SKG => SK
        character(*, SK), intent(in)    :: path
        logical(LK)     , intent(out)   :: failed
        character(*, SK), intent(inout) :: errmsg
        character(:,SKG), allocatable   :: pathExpandedUser
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the attempt to set the current working directory to the input `dir` **fails**, otherwise return `.false.` upon success.<br>
    !>
    !>  \details
    !>  The name of this interface is intentionally chosen to be different from `chdir`
    !>  which is an Intel `ifport` and GNU Fortran extension procedure name.<br>
    !>
    !>  \param[in]  dir     :   The input scalar `character` of default kind \SK containing the directory path to which the current working directory should switch.<br>
    !>
    !>  \return
    !>  `failed`            :   The output scalar `logical` of default kind \LK.<br>
    !>                          It is `.true.` <b>if and only if</b> an error occurs while setting the current working directory to the specified input `dir`.<br>
    !>
    !>  \interface{isFailedChangeDir}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: isFailedChangeDir
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>
    !>      failed = isFailedChangeDir(dir)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [isFailedGetDirTemp](@ref pm_sysShell::isFailedGetDirTemp)<br>
    !>  [isFailedMakeDirTemp](@ref pm_sysPath::isFailedMakeDirTemp)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{isFailedChangeDir}
    !>  \include{lineno} example/pm_sysPath/isFailedChangeDir/main.F90
    !>  \compilef{isFailedChangeDir}
    !>  \output{isFailedChangeDir}
    !>  \include{lineno} example/pm_sysPath/isFailedChangeDir/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \plow
    !>  A subroutine version of this functional interface could be implemented in
    !>  future to avoid allocations and allow for non-default character kinds.<br>
    !>
    !>  \final{isFailedChangeDir}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedChangeDir

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function isFailedChangeDir(dir) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedChangeDir
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    :: dir
        logical(LK)                     :: failed
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the attempt to create the requested directory path **fails**, otherwise return `.false.`.<br>
    !>
    !>  \details
    !>  The requested (potentially **nested**) directory is generated by first identifying the processor shell type (CMD, PowerShell, Bash, csh, zsh, ...)
    !>  and then calling the runtime system shell via the Fortran intrinsic subroutine `execute_command_line()` to invoke `mkdir` shell command with appropriate flags.<br>
    !>  Consequently, this procedure should function as expected in Windows CMD, PowerShell, and POSIX-compatible (Unix-like) shells.<br>
    !>
    !>  **Why not use the C library `mkdir()`?**<br>
    !>  While the C interface is nice and simple, it does not create nested directories, which leads to complications.<br>
    !>  In the end, a Fortran-based solution was deemed more appropriate than calling the C libraries on different platforms.<br>
    !>  The procedures under this generic interface open a subprocess to call the shell command `mkdir` to achieve the goal.<br>
    !>  This is likely less efficient, but has the benefit of automatically enabling the creation of the nested directories without
    !>  the need to process and split the input path into multiple directories.<br>
    !>
    !>  \param[in]      path        :   The input scalar `character` of default kind \SK containing a POSIX-style or Windows-style path.<br>
    !>  \param[in]      wait        :   The input scalar `logical` of default kind \LK with the same functionality as the `wait` argument of the Fortran intrinsic `execute_command_line()`.<br>
    !>                                  <ol>
    !>                                      <li>    If `.true.`, the procedure will wait for the `mkdir` action to terminate before returning the control to the Fortran program.<br>
    !>                                      <li>    Otherwise, the process of making the new directory will be done **asynchronously**.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[in]      ntry        :   The input scalar positive `integer` of default kind \IK representing the number of times the requested action should be attempted.<br>
    !>                                  Multiple tries are particularly useful on Windows platforms where applications take the ownership of a
    !>                                  particular resource on the system and may not allow other applications to utilize the resource.<br>
    !>                                  For example, Dropbox is a well-known example of an application that blocks any other applications from making any changes
    !>                                  to a specific path with which it is working, preventing any manipulation of the path by other applications.<br>
    !>                                  (**optional**, default = `1`)
    !>  \param[inout]   errmsg      :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                                  If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                                  If the length of `errmsg` is too short for the output error message, the message tail will be clipped as needed.<br>
    !>                                  A length of `2047` characters for `errmsg` is likely enough to capture most error messages in full.<br>
    !>                                  (**optional**. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `failed`                    :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b> an error occurs while creating the requested directory.<br>
    !>                                  This includes the possibility of the input directory `path` not being detected on the system after making it (only if `wait = .true.`).<br>
    !>
    !>  \interface{isFailedMakeDir}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: isFailedMakeDir
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>
    !>      failed = isFailedMakeDir(path, wait = wait, ntry = ntry, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that this procedure does not verify the conformance of the input
    !>  path to the naming conventions of the current processor shell. For example,
    !>      -#  if the shell is POSIX-conforming, then the path is expected to conform to the POSIX conventions.<br>
    !>      -#  if the shell is Windows CMD or PowerShell, the path is expected to conform to the Windows conventions.<br>
    !>
    !>  \warning
    !>  **The procedures of this generic interface interpret the input path verbatim (as is)**.<br>
    !>  This is done to minimize security vulnerabilities and to properly handle the presence of whitespace or other potentially problematic characters.<br>
    !>
    !>  \warning
    !>  **Illegal characters in paths on Windows**.<br>
    !>  On Windows platforms with CMD or PowerShell runtime shell, this procedure returns `.true.`
    !>  indicating failure if the input `path` contains Windows reserved characters or names that are illegal in paths.<br>
    !>  These include Windows device names [WINDOWS_RESERVED_DEVICE_NAME](@ref pm_sysPath::WINDOWS_RESERVED_DEVICE_NAME)
    !>  or Windows reserved characters [WINDOWS_RESERVED_STR](@ref pm_sysPath::WINDOWS_RESERVED_STR).<br>
    !>
    !>  \warning
    !>  **The input `path` must not exist prior to calling [isFailedMakeDir](@ref pm_sysPath::isFailedMakeDir)**.<br>
    !>  This procedure does not check the existence of the input directory path.<br>
    !>  Consequently, if it already exists, it may lead to a runtime error.<br>
    !>  The existence of the input `path` can be readily checked via [idsir](@ref pm_sysPath::isDir).<br>
    !>
    !>  \warning
    !>  **Asynchronous directory creation**.<br>
    !>  Note that even if the processor supports asynchronous command execution (`wait = .false.`),
    !>  there is no mechanism provided for finding out later whether the command being executed
    !>  asynchronously has terminated or what its exit status `exitstat` was.<br>
    !>
    !>  \warning
    !>  **Nested directory creation when the shell type cannot be identified**.<br>
    !>  If the procedure cannot identify the runtime system shell, then it simply resorts to calling the `mkdir` command that
    !>  is recognized by almost all shells, albeit without adding any flags for creating nested directories.<br>
    !>  In such a scenario, the nested directory generation might fail if the shell requires `mkdir` with specific flags, or other commands.<br>
    !>
    !>  \warning
    !>  The input `ntry` must be a positive (non-zero) `integer`.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \remark
    !>  This procedure automatically generates nested directories if needed.<br>
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [isFailedGetDirTemp](@ref pm_sysShell::isFailedGetDirTemp)<br>
    !>  [isFailedMakeDirTemp](@ref pm_sysPath::isFailedMakeDirTemp)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{isFailedMakeDir}
    !>  \include{lineno} example/pm_sysPath/isFailedMakeDir/main.F90
    !>  \compilef{isFailedMakeDir}
    !>  \output{isFailedMakeDir}
    !>  \include{lineno} example/pm_sysPath/isFailedMakeDir/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \pmed
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  This generic interface should be extended to in the future to support non-default character kinds when they fully supported by the Fortran intrinsic procedures.<br>
    !>
    !>  \final{isFailedMakeDir}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedMakeDir

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function isFailedMakeDir(path, wait, ntry, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedMakeDir
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)            , intent(in)                    :: path
        logical(LK)                 , intent(in)    , optional      :: wait
        integer(IK)                 , intent(in)    , optional      :: ntry
        character(*, SK)            , intent(inout) , optional      :: errmsg
        logical(LK)                                                 :: failed
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the attempt to create a temporary directory **fails**.<br>
    !>  Otherwise return `.false.` and set the output `path` to the newly created temporary directory.<br>
    !>
    !>  \details
    !>  The procedures of this generic interface attempt to create a unique temporary directory wihtin a parent directory that is either,
    !>  -#  specified by the user as the input `parent`, or
    !>  -#  the system-designated temporary directory returned by [isFailedGetDirTemp](@ref pm_sysShell::isFailedGetDirTemp), or
    !>  -#  the current working directory.<br>
    !>
    !>  \param[out]     path        :   The output `allocatable` scalar `character` of default kind \SK containing a POSIX-style
    !>                                  or Windows-style path to the temporary directory the procedure creates.<br>
    !>                                  **If the procedure fails, the allocation status or contents of `path` will be undefined.**
    !>  \param[in]      parent      :   The input scalar `character` of default kind \SK containing the parent directory within
    !>                                  which the temporary directory should be created.<br>
    !>                                  (**optional**, default = either the system-specified temporary directory or the current working directory)
    !>  \param[inout]   errmsg      :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                                  If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                                  If the length of `errmsg` is too short for the output error message, the message tail will be clipped as needed.<br>
    !>                                  A length of `2047` characters for `errmsg` is likely enough to capture most error messages in full.<br>
    !>                                  (**optional**. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `failed`                    :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b> an error occurs while creating the temporary directory.<br>
    !>
    !>  \interface{isFailedMakeDirTemp}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: isFailedMakeDirTemp
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>
    !>      failed = isFailedMakeDirTemp(path, parent = parent, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [isFailedMakeDir](@ref pm_sysPath::isFailedMakeDir)<br>
    !>  [isFailedGetDirTemp](@ref pm_sysShell::isFailedGetDirTemp)<br>
    !>
    !>  \example{isFailedMakeDirTemp}
    !>  \include{lineno} example/pm_sysPath/isFailedMakeDirTemp/main.F90
    !>  \compilef{isFailedMakeDirTemp}
    !>  \output{isFailedMakeDirTemp}
    !>  \include{lineno} example/pm_sysPath/isFailedMakeDirTemp/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \pmed
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  This generic interface should be extended to in the future to support non-default character kinds when they fully supported by the Fortran intrinsic procedures.<br>
    !>
    !>  \final{isFailedMakeDirTemp}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedMakeDirTemp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function isFailedMakeDirTemp(path, parent, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedMakeDirTemp
#endif
        use pm_kind, only: SKG => SK
        character(:,SKG)            , intent(out)   , allocatable   :: path
        character(*,SKG)            , intent(in)    , optional      :: parent
        character(*, SK)            , intent(inout) , optional      :: errmsg
        logical(LK)                                                 :: failed
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the requested copy action **failed**, otherwise return `.false.` indicating success.<br>
    !>
    !>  \details
    !>  The procedures of this generic interface perform the copy action depending on the type of the two input paths:
    !>  -#  If `from` is an **existing file**, then,
    !>      -#  If `to` is an **existing directory**, the source file `from` will be copied (with the same old file base name) into the existing destination folder.<br>
    !>          Example:
    !>          -#  `from = "./file.txt"`, `to = "./nested/existing/dir/"` (Upon return, a new file `file.txt` will exist in the destination folder `to`).<br>
    !>      -#  If `to` is an **existing file**, then the source file will overwrite the existing destination file <b>if and only if</b> the optional argument `forced = .true.` is specified.<br>
    !>          Example:
    !>          -#  `from = "./file.txt"`, `to = "./nested/existing/dir/existing.file"` (Upon return, `existing.file` will have been overwritten only if `forced = .true.`).<br>
    !>      -#  If `to` is a **non-existing path**, then,
    !>          -#  if `to` ends with any of the directory separators recognized by the runtime shell, then `to` is interpreted as a **non-existing destination folder**. If so,
    !>              -#  Firstly, the (potentially nested non-existing) destination directory `to` will be created.<br>
    !>              -#  Secondly, the source file `from` will be copied **into** the newly created destination folder `to`.<br>
    !>                  Example:
    !>                  -#  (Unix/Windows) `from = "./file.txt"`, `to = "./nested/non-existing/path/"` (Upon return, a new file `file.txt` will exist in the destination folder `to`).<br>
    !>                  -#  (Windows) `from = ".\file.txt"`, `to = ".\nested\non-existing\path\"` (Upon return, a new file `file.txt` will exist in the destination folder `to`).<br>
    !>          -#  if `to` does **not** end with any of the directory separators recognized by the runtime shell, then `to` is interpreted as a **non-existing destination file**. If so,
    !>              -#  Firstly, the (potentially nested non-existing) **dirname** of the destination file will be created.<br>
    !>              -#  Secondly, the source file `from` will be copied to the destination folder with a new name that is the base name of `to`.<br>
    !>                  Example:
    !>                  -#  `from = "./file.txt"`, `to = "./nested/non-existing/path"` (A new file named `path` will be created containing the contents of `file.txt`).<br>
    !>  -#  If `from` is an **existing directory**, then,
    !>      -#  If `to` is an **existing directory**, the entire source directory `from` (along with its basename) will be copied **into** the existing destination folder.<br>
    !>          -#  Existing files with the same name as in the source folder `from` will be overwritten only if the optional argument `forced = .true.` is specified.<br>
    !>              Example:
    !>              -#  `from = "./paramonte/examples/"`, `to = "./newdir/"` (A new directory `"./newdir/examples/"` will be created containing the contents of `from`).<br>
    !>      -#  If `to` is an **existing file**, the program will **not** perform the requested copy action and will return `failed = .true.`.<br>
    !>      -#  If `to` is a **non-existing path**, then `to` is interpreted as a non-existing destination directory.<br>
    !>          -#  Firstly, the (potentially nested non-existing) destination directory `to` will be created.<br>
    !>          -#  Secondly, the **contents** of the source directory `from` will be copied **into** the newly-created destination directory `to`.<br>
    !>          -#  This behavior is similar to the behavior of `cp` linux command, except that this procedure automatically creates the **dirname** of `to` if it does not exist.<br>
    !>              Example:
    !>              -#  `from = "./paramonte/examples/"`, `to = "./new/nested/path"` (non-existing). Result: The newly-created directory `to = "./new/nested/path"` will contain the **contents** of `from`.<br>
    !>  -#  If `from` is neither an existing file nor an existing folder, the program will **not** perform the requested copy action and will return `failed = .true.`.<br>
    !>
    !>  **This procedure is capable of copying a file or a directory to nonexistent (potentially nested) directories.**
    !>
    !>  The copy request is done by first identifying the processor shell type (CMD, PowerShell, Bash, csh, zsh, ...)
    !>  and then calling the runtime system shell via the Fortran intrinsic subroutine `execute_command_line()` to invoke `mkdir` shell command with appropriate flags.<br>
    !>  Consequently, this procedure should function as expected in Windows CMD, PowerShell, and POSIX-compatible (Unix-like) shells.<br>
    !>  -#  On POSIX-compliant shells, the commands `cp -arn` and `cp -arf` are used to perform the copy regularly or forcefully (with overwriting), respectively.<br>
    !>  -#  On Windows-based terminals, the commands `xcopy /e /q /s` and `xcopy /e /q /s /Y` are used to perform the copy regularly or forcefully (with overwriting), respectively.<br>
    !>
    !>  \param[in]      from        :   The input scalar `character` of default kind \SK containing a POSIX-style or Windows-style path to a file or directory.<br>
    !>  \param[in]      to          :   The input scalar `character` of default kind \SK containing a POSIX-style or Windows-style path to which the path `from` should be copied.<br>
    !>                                  When the input argument `from` points to a single file, there is ambiguity about how the destination path should be interpreted.<br>
    !>                                  As such, if `to` points to a (potentially non-existing) folder, make sure to end it with a suitable directory separator [getDirSep](@ref pm_sysPath::getDirSep)
    !>                                  before passing it to the procedures under this generic integer.<br>
    !>                                  Otherwise, `to` is assumed to point to a file name.<br>
    !>  \param[in]      recursive   :   The input scalar `logical` of default kind \LK. If `.true.` the copy action will be performed recursively (including all subdirectories).<br>
    !>                                  This input argument is relevant only if the input source path `from` is a directory with subdirectories.<br>
    !>                                  (**optional**, default = `.false.`. **See the warning below regarding the implications of setting `recursive = .true.`**.)
    !>  \param[in]      forced      :   The input scalar `logical` of default kind \LK. If `.true.` the copy action will overwrite the destination `to` if it already exists.<br>
    !>                                  This input argument is relevant only if the input destination path partially or entirely contains the same files in `from`.<br>
    !>                                  (**optional**, default = `.false.`)
    !>  \param[in]      wait        :   The input scalar `logical` of default kind \LK with the same functionality as the `wait` argument of the Fortran intrinsic `execute_command_line()`.<br>
    !>                                  <ol>
    !>                                      <li>    If `.true.`, the procedure will wait for the shell action to terminate before returning the control to the Fortran program.<br>
    !>                                      <li>    Otherwise, the system call will be done **asynchronously** and the success of the action will not be verified.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[in]      ntry        :   The input scalar positive `integer` of default kind \IK representing the number of times the requested action should be attempted.<br>
    !>                                  Multiple tries are particularly useful on Windows platforms where applications take the ownership of a
    !>                                  particular resource on the system and may not allow other applications to utilize the resource.<br>
    !>                                  For example, Dropbox is a well-known example of an application that blocks any other applications from making any changes
    !>                                  to a specific path with which it is working, preventing any manipulation of the path by other applications.<br>
    !>                                  (**optional**, default = `1`)
    !>  \param[inout]   errmsg      :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                                  If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                                  If the length of `errmsg` is too short for the output error message, the message tail will be clipped as needed.<br>
    !>                                  A length of `2047` characters for `errmsg` is likely enough to capture most error messages in full.<br>
    !>                                  (**optional**. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `failed`                    :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b> an error occurs while performing the requested task.<br>
    !>                                  This includes the possibility of the input path `from` not existing before the copy action or
    !>                                  the destination path `to` not existing after the copy action.<br>
    !>
    !>  \interface{isFailedCopy}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK, IK, SK
    !>      use pm_sysPath, only: isFailedCopy
    !>      character(2047, SK) :: errmsg
    !>      logical(LK) :: failed, forced
    !>      logical(LK) :: recursive
    !>      logical(LK) :: wait
    !>      integer(IK) :: ntry
    !>
    !>      failed = isFailedCopy(from, to, recursive = recursive, forced = forced, wait = wait, ntry = ntry, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  **When does the option `recursive = .true.` lead to infinite copy loop?**<br>
    !>  When `to` is a subdirectory of `from`.<br>
    !>  Be warned that setting the input argument `recursive = .true.` when `to` is a subdirectory of `from`
    !>  can lead to a never ending copy cycle which will likely eventually cause the process to crash.<br>
    !>  Use `recursive = .true.` only if `to` is guaranteed to not be a subdirectory of `from`.<br>
    !>
    !>  \warning
    !>  This procedure currently does not copy hidden files consciously. The result is currently platform and usage dependent.<br>
    !>
    !>  \warning
    !>  Note that this procedure does not verify the conformance of the input
    !>  paths to the naming conventions of the current processor shell. For example,
    !>      -#  if the shell is POSIX-conforming, then the path is expected to conform to the POSIX conventions.<br>
    !>      -#  if the shell is Windows CMD or PowerShell, the path is expected to conform to the Windows conventions.<br>
    !>
    !>  \warning
    !>  **This procedure interprets the input paths verbatim**.<br>
    !>  This procedure potentially calls the operating system to create the input destination path directories.<br>
    !>  To minimize security vulnerabilities and to properly handle the presence of whitespace or other potentially problematic
    !>  characters, the input path will be interpreted verbatim by calling [getPathVerbatim](@ref pm_sysPath::getPathVerbatim).<br>
    !>
    !>  \warning
    !>  **Asynchronous copy action**.<br>
    !>  Note that even if the processor supports asynchronous command execution (`wait = .false.`),
    !>  there is no mechanism provided for finding out later whether the command being executed
    !>  asynchronously has terminated or what its exit status `exitstat` was.<br>
    !>
    !>
    !>  \warning
    !>  The input `ntry` must be a positive (non-zero) `integer`.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isFailedCopy](@ref pm_sysPath::isFailedCopy)<br>
    !>  [isFailedMove](@ref pm_sysPath::isFailedMove)<br>
    !>  [isFailedRemove](@ref pm_sysPath::isFailedRemove)<br>
    !>  [isFailedMakeDir](@ref pm_sysPath::isFailedMakeDir)<br>
    !>
    !>  \example{isFailedCopy}
    !>  \include{lineno} example/pm_sysPath/isFailedCopy/main.F90
    !>  \compilef{isFailedCopy}
    !>  \output{isFailedCopy}
    !>  \include{lineno} example/pm_sysPath/isFailedCopy/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \pmed
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  This generic interface should be extended to in the future to support non-default character kinds when they fully supported by the Fortran intrinsic procedures.<br>
    !>
    !>  \final{isFailedCopy}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedCopy
    module function isFailedCopy(from, to, recursive, forced, wait, ntry, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedCopy
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                :: from, to
        logical(LK)     , intent(in)    , optional  :: recursive
        logical(LK)     , intent(in)    , optional  :: forced
        logical(LK)     , intent(in)    , optional  :: wait
        integer(IK)     , intent(in)    , optional  :: ntry
        character(*, SK), intent(inout) , optional  :: errmsg
        logical(LK)                                 :: failed
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the requested move action **failed**, otherwise return `.false.` indicating success.<br>
    !>
    !>  \details
    !>  The procedures of this generic interface perform the move action depending on the type of the two input paths:
    !>  -#  If `from` is an existing file, then,
    !>      -#  If `to` is an existing directory, the source file `from` will be moved (with the same old file base name) into the existing destination folder.<br>
    !>      -#  If `to` is an existing file, then the source file will overwrite the existing destination file <b>if and only if</b> the optional argument `forced = .true.` is specified.<br>
    !>      -#  If `to` is a non-existing path, then `to` is interpreted as a non-existing destination file. If so,
    !>          -#  The (potentially nested non-existing) **dirname** of the destination file will be created.<br>
    !>          -#  The source file `from` will be moved to the destination folder with a new name that is the base name of `to`.<br>
    !>  -#  If `from` is an existing directory, then,
    !>      -#  If `to` is an existing directory, the **contents** of the source directory `from` will be moved **into** the existing destination folder.<br>
    !>          -#  Existing files with the same name as in the source folder `from` will be overwritten only if the optional argument `forced = .true.` is specified.<br>
    !>      -#  If `to` is an existing file, the program will not perform the requested move action and will return `failed = .true.`.<br>
    !>      -#  If `to` is a non-existing path, then `to` is interpreted as a non-existing destination directory.<br>
    !>          -#  Firstly, the (potentially nested) non-existing directory will be generated.<br>
    !>          -#  Secondly, the **contents** of the source directory `from` will be moved **into** the newly-created destination directory `to`.<br>
    !>  -#  If `from` is neither an existing file nor an existing folder, the program will not perform the requested move action and will return `failed = .true.`.<br>
    !>
    !>  **This procedure is capable of moving a file or a directory to nonexistent (potentially nested) directories.**
    !>
    !>  The move request is done by first identifying the processor shell type (CMD, PowerShell, Bash, csh, zsh, ...)
    !>  and then calling the runtime system shell via the Fortran intrinsic subroutine `execute_command_line()` to invoke `mkdir` shell command with appropriate flags.<br>
    !>  Consequently, this procedure should function as expected in Windows CMD, PowerShell, and POSIX-compatible (Unix-like) shells.<br>
    !>  -#  On POSIX-compliant shells, the commands `mv -arn` and `mv -arf` are used to perform the move regularly or forcefully (with overwriting), respectively.<br>
    !>  -#  On Windows-based terminals, the commands `xmove /e /q /s` and `xmove /e /q /s /Y` are used to perform the move regularly or forcefully (with overwriting), respectively.<br>
    !>
    !>  \param[in]      from        :   The input scalar `character` of default kind \SK containing a POSIX-style or Windows-style path to a file or directory.<br>
    !>  \param[in]      to          :   The input scalar `character` of default kind \SK containing a POSIX-style or Windows-style path to which the path `from` should be moved.<br>
    !>                                  When the input argument `from` points to a single file, there is ambiguity about how the destination path should be interpreted.<br>
    !>                                  As such, if `to` points to a (potentially non-existing) folder, make sure to end it with a suitable directory separator [getDirSep](@ref pm_sysPath::getDirSep)
    !>                                  before passing it to the procedures under this generic integer.<br>
    !>                                  Otherwise, `to` is assumed to point to a file name.<br>
    !>  \param[in]      forced      :   The input scalar `logical` of default kind \LK. If `.true.` the move action will overwrite the destination `to` if it already exists.<br>
    !>                                  This input argument is relevant only if the input destination path partially or entirely contains the same files in `from`.<br>
    !>                                  (**optional**, default = `.false.`)
    !>  \param[in]      wait        :   The input scalar `logical` of default kind \LK with the same functionality as the `wait` argument of the Fortran intrinsic `execute_command_line()`.<br>
    !>                                  <ol>
    !>                                      <li>    If `.true.`, the procedure will wait for the shell action to terminate before returning the control to the Fortran program.<br>
    !>                                      <li>    Otherwise, the system call will be done **asynchronously** and the success of the action will not be verified.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[in]      ntry        :   The input scalar positive `integer` of default kind \IK representing the number of times the requested action should be attempted.<br>
    !>                                  Multiple tries are particularly useful on Windows platforms where applications take the ownership of a
    !>                                  particular resource on the system and may not allow other applications to utilize the resource.<br>
    !>                                  For example, Dropbox is a well-known example of an application that blocks any other applications from making any changes
    !>                                  to a specific path with which it is working, preventing any manipulation of the path by other applications.<br>
    !>                                  (**optional**, default = `1`)
    !>  \param[inout]   errmsg      :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                                  If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                                  If the length of `errmsg` is too short for the output error message, the message tail will be clipped as needed.<br>
    !>                                  A length of `2047` characters for `errmsg` is likely enough to capture most error messages in full.<br>
    !>                                  (**optional**. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `failed`                    :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b> an error occurs while performing the requested task.<br>
    !>                                  This includes the possibility of the input path `from` not existing before the move action or
    !>                                  the destination path `to` not existing after the move action.<br>
    !>
    !>  \interface{isFailedMove}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK, IK, SK
    !>      use pm_sysPath, only: isFailedMove
    !>      character(2047, SK) :: errmsg
    !>      logical(LK) :: failed, forced
    !>      logical(LK) :: wait
    !>      integer(IK) :: ntry
    !>
    !>      failed = isFailedMove(from, to, forced = forced, wait = wait, ntry = ntry, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that this procedure does not verify the conformance of the input
    !>  paths to the naming conventions of the current processor shell. For example,
    !>      -#  if the shell is POSIX-conforming, then the path is expected to conform to the POSIX conventions.<br>
    !>      -#  if the shell is Windows CMD or PowerShell, the path is expected to conform to the Windows conventions.<br>
    !>
    !>  \warning
    !>  **This procedure interprets the input paths verbatim**.<br>
    !>  This procedure potentially calls the operating system to create the input destination path directories.<br>
    !>  To minimize security vulnerabilities and to properly handle the presence of whitespace or other potentially problematic
    !>  characters, the input path will be interpreted verbatim by calling [getPathVerbatim](@ref pm_sysPath::getPathVerbatim).<br>
    !>
    !>  \warning
    !>  **Asynchronous move action**.<br>
    !>  Note that even if the processor supports asynchronous command execution (`wait = .false.`),
    !>  there is no mechanism provided for finding out later whether the command being executed
    !>  asynchronously has terminated or what its exit status `exitstat` was.<br>
    !>
    !>
    !>  \warning
    !>  The input `ntry` must be a positive (non-zero) `integer`.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isFailedCopy](@ref pm_sysPath::isFailedCopy)<br>
    !>  [isFailedMove](@ref pm_sysPath::isFailedMove)<br>
    !>  [isFailedRemove](@ref pm_sysPath::isFailedRemove)<br>
    !>  [isFailedMakeDir](@ref pm_sysPath::isFailedMakeDir)<br>
    !>
    !>  \example{isFailedMove}
    !>  \include{lineno} example/pm_sysPath/isFailedMove/main.F90
    !>  \compilef{isFailedMove}
    !>  \output{isFailedMove}
    !>  \include{lineno} example/pm_sysPath/isFailedMove/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{isFailedMove}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedMove
    module function isFailedMove(from, to, forced, wait, ntry, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedMove
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                :: from, to
        logical(LK)     , intent(in)    , optional  :: forced
        logical(LK)     , intent(in)    , optional  :: wait
        integer(IK)     , intent(in)    , optional  :: ntry
        character(*, SK), intent(inout) , optional  :: errmsg
        logical(LK)                                 :: failed
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the requested path removal action from the file system **failed**, otherwise return `.false.` indicating success.<br>
    !>
    !>  \details
    !>  The procedures of this generic interface perform the path removal action whether the input path points to a file or an empty or non-empty directory.<br>
    !>  The removal strategy depends on the nature of the input `path`:<br>
    !>  <ol>
    !>      <li>    If the input `path` points to a single file, then the path is removed from the file system using Fortran intrinsic routines.<br>
    !>      <li>    If the input `path` points to a multiple files via a pattern specification or a directory, then the path is removed via the relevant shell commands.<br>
    !>  </ol>
    !>
    !>  \param[in]      path        :   The input scalar `character` of default kind \SK containing a POSIX-style or Windows-style path to a file or directory.<br>
    !>  \param[in]      recursive   :   The input scalar `logical` of default kind \LK. If `.true.` the removal action will be performed recursively (including all subdirectories).<br>
    !>                                  This input argument is relevant only if the input `path` is a directory potentially with subdirectories.<br>
    !>                                  (**optional**, default = `.false.`)
    !>  \param[in]      forced      :   The input scalar `logical` of default kind \LK. If `.true.` the removal action will be performed without any questions asked (and nonexistent files will be ignored).<br>
    !>                                  Note that setting `force = .true.` is required to ignore errors raised on Windows related to non-existence of the specified path or pattern.<br>
    !>                                  (**optional**, default = `.false.`)
    !>  \param[in]      wait        :   The input scalar `logical` of default kind \LK with the same functionality as the `wait` argument of the Fortran intrinsic `execute_command_line()`.<br>
    !>                                  <ol>
    !>                                      <li>    If `.true.`, the procedure will wait for the shell action to terminate before returning the control to the Fortran program.<br>
    !>                                      <li>    Otherwise, the system call will be done **asynchronously** and the success of the action will not be verified.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = `.true.`)
    !>  \param[in]      ntry        :   The input scalar positive `integer` of default kind \IK representing the number of times the requested action should be attempted.<br>
    !>                                  Multiple tries are particularly useful on Windows platforms where applications take the ownership of a
    !>                                  particular resource on the system and may not allow other applications to utilize the resource.<br>
    !>                                  For example, Dropbox is a well-known example of an application that blocks any other applications from making any changes
    !>                                  to a specific path with which it is working, preventing any manipulation of the path by other applications.<br>
    !>                                  (**optional**, default = `1`)
    !>  \param[inout]   errmsg      :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                                  If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                                  If the length of `errmsg` is too short for the output error message, the message tail will be clipped as needed.<br>
    !>                                  A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                                  (**optional**. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `failed`                    :   The output scalar `logical` of default kind \LK. It is `.true.` <b>if and only if</b> an error occurs while performing the requested task.<br>
    !>                                  This includes the possibility of the input `path` not existing before the removal action.<br>
    !>
    !>  \interface{isFailedRemove}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK, IK, SK
    !>      use pm_sysPath, only: isFailedRemove
    !>      character(2047, SK) :: errmsg
    !>      logical(LK) :: failed, forced
    !>      logical(LK) :: recursive
    !>      logical(LK) :: wait
    !>      integer(IK) :: ntry
    !>
    !>      failed = isFailedRemove(path, recursive = recursive, forced = forced, wait = wait, ntry = ntry, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that this procedure does not verify the conformance of the input
    !>  paths to the naming conventions of the current processor shell. For example,
    !>      -#  if the shell is POSIX-conforming, then the path is expected to conform to the POSIX conventions.<br>
    !>      -#  if the shell is Windows CMD or PowerShell, the path is expected to conform to the Windows conventions.<br>
    !>
    !>  \warning
    !>  **This procedure interprets the input paths verbatim**.<br>
    !>  This procedure potentially calls the operating system to create the input destination path directories.<br>
    !>  To minimize security vulnerabilities and to properly handle the presence of whitespace or other potentially problematic
    !>  characters, the input path will be interpreted verbatim by calling [getPathVerbatim](@ref pm_sysPath::getPathVerbatim).<br>
    !>
    !>  \warning
    !>  **Asynchronous removal action**.<br>
    !>  Note that even if the processor supports asynchronous command execution (`wait = .false.`),
    !>  there is no mechanism provided for finding out later whether the command being executed
    !>  asynchronously has terminated or what its exit status `exitstat` was.<br>
    !>
    !>  \warning
    !>  The input `ntry` must be a positive (non-zero) `integer`.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \recursive
    !>
    !>  \see
    !>  [isFailedCopy](@ref pm_sysPath::isFailedCopy)<br>
    !>  [isFailedMove](@ref pm_sysPath::isFailedMove)<br>
    !>  [isFailedRemove](@ref pm_sysPath::isFailedRemove)<br>
    !>  [isFailedMakeDir](@ref pm_sysPath::isFailedMakeDir)<br>
    !>
    !>  \example{isFailedRemove}
    !>  \include{lineno} example/pm_sysPath/isFailedRemove/main.F90
    !>  \compilef{isFailedRemove}
    !>  \output{isFailedRemove}
    !>  \include{lineno} example/pm_sysPath/isFailedRemove/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{isFailedCopy}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedRemove
    recursive module function isFailedRemove(path, recursive, forced, wait, ntry, errmsg) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedRemove
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                :: path
        logical(LK)     , intent(in)    , optional  :: recursive
        logical(LK)     , intent(in)    , optional  :: forced
        logical(LK)     , intent(in)    , optional  :: wait
        integer(IK)     , intent(in)    , optional  :: ntry
        character(*, SK), intent(inout) , optional  :: errmsg
        logical(LK)                                 :: failed
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the directory separator symbol based on the runtime system shell.<br>
    !>
    !>  \details
    !>  The directory separator returned by this procedure can be used in constructing paths that can be passed to system function
    !>  for example via [isFailedExec](@ref pm_sysShell::isFailedExec) or the Fortran intrinsic `execute_command_line()`.<br>
    !>  -#  The directory separator for Windows-based terminals such as CMD and PowerShell is the backslash character `\`.<br>
    !>  -#  The directory separator for POSIX-compliant terminals such as sh, bash, csh, zsh, ... is the forward slash `/`.<br>
    !>
    !>  \param[in]      mold    :   The input scalar `character` of kind \SKALL of length `1`. The value of `mold` s ignored.<br>
    !>                              However, its kind type parameter is used to to set the kind type parameter of the output `dirsep`.<br>
    !>                              (**optional**, default = `SK_" "` where `SK` refers to the default `character` kind \SK.)
    !>                              <b>if and only if</b> an error occurs while inferring the system shell directory separator.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the system shell directory separator.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `dirsep`                :   The output scalar `character` of the same kind as the input `mold`, of length type parameter `1`,
    !>                              representing the directory separator symbol based on the runtime system shell (`\` or `/`).<br>
    !>                              If the procedure fails and the output argument `failed` is present, it is set to `.true.`.<br>
    !>                              Also, `dirsep` is set to the default value of the `dirsep` component of [shell_type](@ref pm_sysShell::shell_type).<br>
    !>
    !>  \interface{getDirSep}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_sysPath, only: getDirSep
    !>      character(1, SK) :: dirsep
    !>
    !>      dirsep = getDirSep(failed = failed, errmsg = errmsg)
    !>      dirsep = getDirSep(mold, failed = failed, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isDir](@ref pm_sysPath::isDir)<br>
    !>  [isFile](@ref pm_sysPath::isFile)<br>
    !>  [getDirSep](@ref pm_sysPath::getDirSep)<br>
    !>  [getDirSeps](@ref pm_sysPath::getDirSeps)<br>
    !>  [getPathSep](@ref pm_sysPath::getPathSep)<br>
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>  [shellis_type](@ref pm_sysShell::shellis_type)<br>
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>
    !>  \example{getDirSep}
    !>  \include{lineno} example/pm_sysPath/getDirSep/main.F90
    !>  \compilef{getDirSep}
    !>  \output{getDirSep}
    !>  \include{lineno} example/pm_sysPath/getDirSep/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{getDirSep}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDirSep

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getDirSep_SK(failed, errmsg) result(dirsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSep_SK
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: dirsep
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getDirSep_SK5(mold, failed, errmsg) result(dirsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSep_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: dirsep
    end function
#endif

#if SK4_ENABLED
    impure module function getDirSep_SK4(mold, failed, errmsg) result(dirsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSep_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: dirsep
    end function
#endif

#if SK3_ENABLED
    impure module function getDirSep_SK3(mold, failed, errmsg) result(dirsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSep_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: dirsep
    end function
#endif

#if SK2_ENABLED
    impure module function getDirSep_SK2(mold, failed, errmsg) result(dirsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSep_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: dirsep
    end function
#endif

#if SK1_ENABLED
    impure module function getDirSep_SK1(mold, failed, errmsg) result(dirsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSep_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: dirsep
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return all the directory separator symbols supported
    !>  by the current runtime operating system based on the runtime system shell.<br>
    !>
    !>  \details
    !>  The directory separators returned by this generic interface can be used for splitting paths to its components via for example,
    !>  [getDirName](@ref pm_sysPath::getDirName) or the Fortran intrinsic `execute_command_line()`.<br>
    !>  -#  The directory separators for Windows-based terminals such as CMD and PowerShell can be both the backslash character `\` or the forward-slash character `\`.<br>
    !>  -#  The directory separators for POSIX-compliant terminals such as sh, bash, csh, zsh, ... is the forward slash `/`.<br>
    !>
    !>  \param[in]      mold    :   The input scalar `character` of kind \SKALL of length `1`. The value of `mold` s ignored.<br>
    !>                              However, its kind type parameter is used to to set the kind type parameter of the output `dirseps`.<br>
    !>                              (**optional**, default = `SK_" "` where `SK` refers to the default `character` kind \SK.)
    !>                              <b>if and only if</b> an error occurs while inferring the system shell directory separator.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while inferring the system shell directory separator.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `dirseps`               :   The output `allocatable` scalar `character` of the same kind as the input `mold` of arbitrary
    !>                              length type parameter (typically `1` on Unix systems or `2` on Windows), representing
    !>                              the directory separator symbol based on the runtime system shell (`\` or `/`).<br>
    !>                              If the procedure fails and the output argument `failed` is present, it is set to `.true.`.<br>
    !>                              Also, `dirseps` is set to the default value of the `dirseps` component of [shell_type](@ref pm_sysShell::shell_type).<br>
    !>
    !>  \interface{getDirSeps}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_sysPath, only: getDirSeps
    !>      character(:, SK) :: dirseps
    !>      logical(LK) :: failed
    !>
    !>      dirseps = getDirSeps(failed = failed, errmsg = errmsg)
    !>      dirseps = getDirSeps(mold, failed = failed, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isDir](@ref pm_sysPath::isDir)<br>
    !>  [isFile](@ref pm_sysPath::isFile)<br>
    !>  [getDirSep](@ref pm_sysPath::getDirSep)<br>
    !>  [getDirSeps](@ref pm_sysPath::getDirSeps)<br>
    !>  [getPathSep](@ref pm_sysPath::getPathSep)<br>
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>  [shellis_type](@ref pm_sysShell::shellis_type)<br>
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>
    !>  \example{getDirSeps}
    !>  \include{lineno} example/pm_sysPath/getDirSeps/main.F90
    !>  \compilef{getDirSeps}
    !>  \output{getDirSeps}
    !>  \include{lineno} example/pm_sysPath/getDirSeps/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{getDirSeps}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDirSeps

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function getDirSeps_SK(failed, errmsg) result(dirseps)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSeps_SK
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(:,SKG), allocatable               :: dirseps
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getDirSeps_SK5(mold, failed, errmsg) result(dirseps)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSeps_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(:,SKG), allocatable               :: dirseps
    end function
#endif

#if SK4_ENABLED
    impure module function getDirSeps_SK4(mold, failed, errmsg) result(dirseps)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSeps_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(:,SKG), allocatable               :: dirseps
    end function
#endif

#if SK3_ENABLED
    impure module function getDirSeps_SK3(mold, failed, errmsg) result(dirseps)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSeps_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(:,SKG), allocatable               :: dirseps
    end function
#endif

#if SK2_ENABLED
    impure module function getDirSeps_SK2(mold, failed, errmsg) result(dirseps)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSeps_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(:,SKG), allocatable               :: dirseps
    end function
#endif

#if SK1_ENABLED
    impure module function getDirSeps_SK1(mold, failed, errmsg) result(dirseps)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirSeps_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(:,SKG), allocatable               :: dirseps
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the path separator character used by the runtime processor shell
    !>  for separating paths from each other in the `PATH` environmental variable.<br>
    !>
    !>  \details
    !>  The path separator character returned by this generic interface can be used for splitting
    !>  the contents of the processor environmental `PATH` variable to separate directory paths.<br>
    !>  -#  The path separator for Windows-based terminals such as CMD and PowerShell is [PATH_SEP_WINDOWS](@ref pm_sysPath::PATH_SEP_WINDOWS).<br>
    !>  -#  The path separator for POSIX-compliant terminals such as sh, bash, csh, zsh, ... is [PATH_SEP_POSIX](@ref pm_sysPath::PATH_SEP_POSIX).<br>
    !>
    !>  \param[in]      mold    :   The input scalar `character` of kind \SKALL of length `1`. The value of `mold` s ignored.<br>
    !>                              However, its kind type parameter is used to to set the kind type parameter of the output `dirseps`.<br>
    !>                              (**optional**, default = `SK_" "` where `SK` refers to the default `character` kind \SK.)
    !>                              <b>if and only if</b> an error occurs while inferring the system shell directory separator.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                              <b>if and only if</b> an error occurs while the procedure is performing the task.<br>
    !>                              (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `pathsep`               :   The output scalar `character` of default kind \SK of length type parameter of `1`
    !>                              containing the path separator character used by the runtime processor shell for
    !>                              separating paths from each other in the `PATH` environmental variable.<br>
    !>
    !>  \interface{getPathSep}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathSep
    !>      character(1, SK) :: pathsep
    !>
    !>      pathsep = getPathSep(failed = failed, errmsg = errmsg)
    !>      pathsep = getPathSep(mold, failed = failed, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isDir](@ref pm_sysPath::isDir)<br>
    !>  [isFile](@ref pm_sysPath::isFile)<br>
    !>  [getDirSep](@ref pm_sysPath::getDirSep)<br>
    !>  [getDirSeps](@ref pm_sysPath::getDirSeps)<br>
    !>  [getPathSep](@ref pm_sysPath::getPathSep)<br>
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>  [shellis_type](@ref pm_sysShell::shellis_type)<br>
    !>  [shell_type](@ref pm_sysShell::shell_type)<br>
    !>
    !>  \example{getPathSep}
    !>  \include{lineno} example/pm_sysPath/getPathSep/main.F90
    !>  \compilef{getPathSep}
    !>  \output{getPathSep}
    !>  \include{lineno} example/pm_sysPath/getPathSep/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final{getPathSep}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathSep

    impure module function getPathSep_SK(failed, errmsg) result(pathsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathSep_SK
#endif
        use pm_kind, only: SKG => SK
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: pathsep
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    impure module function getPathSep_SK5(mold, failed, errmsg) result(pathsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathSep_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: pathsep
    end function
#endif

#if SK4_ENABLED
    impure module function getPathSep_SK4(mold, failed, errmsg) result(pathsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathSep_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: pathsep
    end function
#endif

#if SK3_ENABLED
    impure module function getPathSep_SK3(mold, failed, errmsg) result(pathsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathSep_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: pathsep
    end function
#endif

#if SK2_ENABLED
    impure module function getPathSep_SK2(mold, failed, errmsg) result(pathsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathSep_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: pathsep
    end function
#endif

#if SK1_ENABLED
    impure module function getPathSep_SK1(mold, failed, errmsg) result(pathsep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathSep_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG), intent(in)                :: mold
        logical(LK)     , intent(out)   , optional  :: failed
        character(*, SK), intent(inout) , optional  :: errmsg
        character(1,SKG)                            :: pathsep
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a unique (directory or file) path name in the specified directory or the default current working directory.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface return a full unique path that is,<br>
    !>  <ol>
    !>      <li>    optionally prefixed with the input user-specified directory `dir`,<br>
    !>      <li>    optionally followed by the user-specified `prefix`,<br>
    !>      <li>    optionally followed by the file extension `ext`, that can be either<br>
    !>              <ol>
    !>                  <li>    a string with pattern `.*`, <b>if a file path</b> is desired as the output, or<br>
    !>                  <li>    a string that ends with an appropriate [path-separator symbol](@ref DIR_SEP_ALL) of the platform, <b>if a directory path</b> is desired as the output.<br>
    !>              </ol>
    !>  </ol>
    !>  The procedures under this generic interface use the current date and time and the process ID to generate new unique non-existing directory or file paths.<br>
    !>  The general format of the output directory path is the following:<br>
    !>  \verbatim
    !>  <dir><prefix><ccyymmdd><sep><hhmmss><sep><sss><sep>pid<sep><pid><ext>
    !>  \endverbatim
    !>  where,
    !>  <ol>
    !>      <li>    `<dir>`         is replaced with the contents of the input argument `dir` or its default value (followed by the path-separator symbol if needed),
    !>      <li>    `<prefix>`      is replaced with the contents of the input argument `prefix` or its default value,
    !>      <li>    `<sep>`         is replaced with the contents of the input argument `sep` or its default value,
    !>      <li>    `<ccyymmdd>`    is replaced with the current century, year, month, and day, each represented by two characters,
    !>      <li>    `<hhmmss>`      is replaced with the current hour `hh`, minute `mm`, second `ss`, in the specified ordered,
    !>      <li>    `<sss>`         is replaced with the current millisecond `mmm`,
    !>      <li>    `<pid>`         is replaced with the user-specified processor ID `pid` which can be the output of [getImageID](@ref pm_parallelism::getImageID).<br>
    !>      <li>    `<ext>`         is replaced with the contents of the input argument `ext` or its default value.<br>
    !>  </ol>
    !>
    !>  If no unique path is found after `1000` attempts, the procedure returns and reports a failure.<br>
    !>
    !>  \param[in]  dir     :   The input scalar `character` of default kind \SK containing the directory where the new path should reside.<br>
    !>                          If `dir` does not end with the proper system shell path-separator symbol, it will be suffixed with the proper symbol.<br>
    !>                          However, it is recommended to pass an input `dir` that is already suffixed with the proper path-separator.<br>
    !>                          The input `dir` is prefixed to the new path name that is to be constructed.<br>
    !>                          (**optional**, default = `""`)
    !>  \param[in]  prefix  :   The input scalar `character` of default kind \SK containing the requested input prefix to the path
    !>                          that appears immediately after the input `dir` in the newly-constructed path.<br>
    !>                          (**optional**, default = `"new"//sep` where `sep` is the other `optional` argument of the interface.)
    !>  \param[in]  sep     :   The input scalar `character` of default kind \SK containing the requested separator character to be inserted
    !>                          between the individual components of the new path to be constructed.<br>
    !>                          This is **different** from the system shell directory separator.<br>
    !>                          (**optional**, default = `"_"`)
    !>  \param[in]  ext     :   The input scalar `character` of default kind \SK containing the requested input file extension or
    !>                          directory path-separator to appear at the end of the path.<br>
    !>                          (**optional**, default = `""`)
    !>  \param[in]  pid     :   The input scalar `integer` of default kind \IK representing the processor ID (e.g., processor MPI rank or Coarray image ID).<br>
    !>                          If present, it will be used in the generated path to make it processor-specific.<br>
    !>                          The processor MPI rank or Coarray image ID can be for example obtained by calling [getImageID](@ref pm_parallelism::getImageID).<br>
    !>                          (**optional**, default = [getImageID()](@ref pm_parallelism::getImageID))
    !>  \param[out] failed  :   The output `logical` of default kind \LK that is `.true.` <b>if and only if</b> the procedure fails to find a new unique (non-existent) path.<br>
    !>                          (**optional**, if missing and the procedure fails, then the program will halt by calling `error stop`.)
    !>
    !>  \return
    !>  `pathNew`           :   The output `allocatable` scalar of type `character` of default kind \SK containing a new unique non-existent path.<br>
    !>
    !>  \interface{getPathNew}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathNew
    !>
    !>      pathNew = getPathNew(dir = dir, prefix = prefix, sep = sep, ext = ext, failed = failed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Avoid the use of Windows reserved (forbidden) characters [WINDOWS_RESERVED_STR](@ref pm_sysPath::WINDOWS_RESERVED_STR) in path names.<br>
    !>  Although such characters work fine on Unix systems, the Windows OS will not be able to gracefully handle such paths.<br>
    !>  In particular, avoid the use of forward or backslash in path names (it is fine to use them as path-separator symbols).<br>
    !>
    !>  \warning
    !>  The input `ext` argument, if it is present and contains a file extension, **must** the proper file extension separator symbol (e.g., `.`).<br>
    !>  The file extension separator (`.`) will **not** be added to the input extension by this procedure.<br>
    !>
    !>  \remark
    !>  On Windows Operating Systems (OS), the choice of the path-separator depends on the processor runtime shell being used.<br>
    !>  The Windows Batch terminal (CMD) can handle paths that are separated via either [DIR_SEP_WINDOWS](@ref pm_sysPath::DIR_SEP_WINDOWS)
    !>  or [DIR_SEP_POSIX](@ref pm_sysPath::DIR_SEP_POSIX) or even both (although paths containing POSIX separators must be quoted).<br>
    !>  However, Unix-like Shells like Git Bash, Cygwin, etc. only recognize [DIR_SEP_POSIX](@ref pm_sysPath::DIR_SEP_POSIX) even on Windows OS.<br>
    !>  Therefore, to avoid the use of wrong path-separator for the input argument `dir` on Windows OS, the following steps are taken:<br>
    !>      1.  If the runtime shell is Windows-based and the input `dir` does not end with a directory separator recognized by Windows
    !>          ([DIR_SEP_WINDOWS_ALL](@ref pm_sysPath::DIR_SEP_WINDOWS_ALL)), then `dir` will be suffixed with a Windows-style path-separator.<br>
    !>      2.  Otherwise, if the input `dir` does not end with a directory separator ([DIR_SEP_POSIX_ALL](@ref pm_sysPath::DIR_SEP_POSIX_ALL))
    !>          recognized by POSIX-style runtime shells (including Fish), then `dir` will be suffixed with a POSIX-style path-separator.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isFailedMakeDir](@ref pm_sysPath::isFailedMakeDir)<br>
    !>  [isFailedMakeDirTemp](@ref pm_sysPath::isFailedMakeDirTemp)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{getPathNew}
    !>  \include{lineno} example/pm_sysPath/getPathNew/main.F90
    !>  \compilef{getPathNew}
    !>  \output{getPathNew}
    !>  \include{lineno} example/pm_sysPath/getPathNew/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \todo
    !>  \phigh
    !>  An optional `suffix` argument can be added in the future. Currently, the path suffix is hard-coded in the procedure, unlike `prefix`.<br>
    !>
    !>  \final{getPathNew}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathNew
    module function getPathNew(dir, prefix, sep, ext, pid, failed) result(pathNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathNew
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)            , intent(in)    , optional  :: dir, prefix, sep, ext
        integer(IK)                 , intent(in)    , optional  :: pid
        logical(LK)                 , intent(out)   , optional  :: failed
        character(:,SKG)            , allocatable               :: pathNew
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a unique (directory or file) path name in the temporary directory of the system or if there is none,
    !>  in the current working directory of the program, optionally with the user-specified file name segments.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface return a full unique path that is,<br>
    !>  <ol>
    !>      <li>    optionally followed by the user-specified `prefix`,<br>
    !>      <li>    optionally followed by the file extension `ext`, that can be either<br>
    !>              <ol>
    !>                  <li>    a string with pattern `.*`, <b>if a file path</b> is desired as the output, or<br>
    !>                  <li>    a string that ends with an appropriate [path-separator symbol](@ref DIR_SEP_ALL) of the platform, <b>if a directory path</b> is desired as the output.<br>
    !>              </ol>
    !>  </ol>
    !>  The procedures under this generic interface use the current date and time and the process ID to generate new unique non-existing directory or file paths.<br>
    !>  The general format of the output directory path is the following:<br>
    !>  \verbatim
    !>  <tempdir><prefix><ccyymmdd><sep><hhmmss><sep><sss><sep>pid<sep><pid><ext>
    !>  \endverbatim
    !>  where,
    !>  <ol>
    !>      <li>    `<tempdir>`     is replaced with the operating system temporary directory or if none is detected, with the current runtime working directory,
    !>      <li>    `<prefix>`      is replaced with the contents of the input argument `prefix` or its default value,
    !>      <li>    `<sep>`         is replaced with the contents of the input argument `sep` or its default value,
    !>      <li>    `<ccyymmdd>`    is replaced with the current century, year, month, and day, each represented by two characters,
    !>      <li>    `<hhmmss>`      is replaced with the current hour `hh`, minute `mm`, second `ss`, in the specified ordered,
    !>      <li>    `<sss>`         is replaced with the current millisecond `mmm`,
    !>      <li>    `<pid>`         is replaced with the user-specified processor ID `pid` which can be the output of [getImageID](@ref pm_parallelism::getImageID).<br>
    !>      <li>    `<ext>`         is replaced with the contents of the input argument `ext` or its default value.<br>
    !>  </ol>
    !>
    !>  If no unique path is found after `1000` attempts, the procedure returns and reports a failure.<br>
    !>
    !>  \param[in]  prefix  :   The input scalar `character` of default kind \SK containing the requested input prefix to the path
    !>                          that appears first in the newly-constructed path.<br>
    !>                          (**optional**, default = `"new"//sep` where `sep` is the other `optional` argument of the interface.)
    !>  \param[in]  sep     :   The input scalar `character` of default kind \SK containing the requested separator character to be inserted
    !>                          between the individual components of the new path to be constructed.<br>
    !>                          This is **different** from the system shell directory separator.<br>
    !>                          (**optional**, default = `"_"`)
    !>  \param[in]  ext     :   The input scalar `character` of default kind \SK containing the requested input file extension or
    !>                          directory path-separator to appear at the end of the path.<br>
    !>                          (**optional**, default = `""`)
    !>  \param[in]  pid     :   The input scalar `integer` of default kind \IK representing the processor ID (e.g., processor MPI rank or Coarray image ID).<br>
    !>                          If present, it will be used in the generated path to make it processor-specific.<br>
    !>                          The processor MPI rank or Coarray image ID can be for example obtained by calling [getImageID](@ref pm_parallelism::getImageID).<br>
    !>                          (**optional**, default = [getImageID()](@ref pm_parallelism::getImageID))
    !>  \param[out] failed  :   The output `logical` of default kind \LK that is `.true.` <b>if and only if</b> the procedure fails to find a new unique (non-existent) path.<br>
    !>                          (**optional**, if missing and the procedure fails, then the program will halt by calling `error stop`.)
    !>
    !>  \return
    !>  `pathTemp`          :   The output `allocatable` scalar of type `character` of default kind \SK,
    !>                          containing a new unique non-existent temporary (directory or file) path,
    !>                          either in the OS temporary directory or in the current working directory.<br>
    !>
    !>  \interface{getPathTemp}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathTemp
    !>
    !>      pathTemp = getPathTemp(prefix = prefix, sep = sep, ext = ext, failed = failed)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Avoid the use of Windows reserved (forbidden) characters [WINDOWS_RESERVED_STR](@ref pm_sysPath::WINDOWS_RESERVED_STR) in path names.<br>
    !>  Although such characters work fine on Unix systems, the Windows OS will not be able to gracefully handle such paths.<br>
    !>  In particular, avoid the use of forward or backslash in path names (it is fine to use them as path-separator symbols).<br>
    !>
    !>  \warning
    !>  The input `ext` argument, if it is present and contains a file extension, **must** the proper file extension separator symbol (e.g., `.`).<br>
    !>  The file extension separator (`.`) will **not** be added to the input extension by this procedure.<br>
    !>
    !>  \remark
    !>  On Windows Operating Systems (OS), the choice of the path-separator depends on the processor runtime shell being used.<br>
    !>  The Windows Batch terminal (CMD) can handle paths that are separated via either [DIR_SEP_WINDOWS](@ref pm_sysPath::DIR_SEP_WINDOWS)
    !>  or [DIR_SEP_POSIX](@ref pm_sysPath::DIR_SEP_POSIX) or even both (although paths containing POSIX separators must be quoted).<br>
    !>  However, Unix-like Shells like Git Bash, Cygwin, etc. only recognize [DIR_SEP_POSIX](@ref pm_sysPath::DIR_SEP_POSIX) even on Windows OS.<br>
    !>  Therefore, to avoid the use of wrong path-separator for the parent directory of the path on Windows OS, the following steps are taken:<br>
    !>      1.  If the runtime shell is Windows-based and the temporary directory does not end with a directory separator recognized by Windows
    !>          ([DIR_SEP_WINDOWS_ALL](@ref pm_sysPath::DIR_SEP_WINDOWS_ALL)), the directory will be suffixed with a Windows-style path-separator.<br>
    !>      2.  Otherwise, if the temporary directory does not end with a directory separator ([DIR_SEP_POSIX_ALL](@ref pm_sysPath::DIR_SEP_POSIX_ALL))
    !>          recognized by POSIX-style runtime shells (including Fish), then the temporary directory will be suffixed with a POSIX-style path-separator.<br>
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathTemp](@ref pm_sysPath::getPathTemp)<br>
    !>  [isFailedMakeDir](@ref pm_sysPath::isFailedMakeDir)<br>
    !>  [isFailedMakeDirTemp](@ref pm_sysPath::isFailedMakeDirTemp)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{getPathTemp}
    !>  \include{lineno} example/pm_sysPath/getPathTemp/main.F90
    !>  \compilef{getPathTemp}
    !>  \output{getPathTemp}
    !>  \include{lineno} example/pm_sysPath/getPathTemp/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \todo
    !>  \phigh
    !>  An optional `suffix` argument can be added in the future. Currently, the path suffix is hard-coded in the procedure, unlike `prefix`.<br>
    !>
    !>  \final{getPathTemp}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathTemp
    module function getPathTemp(prefix, sep, ext, pid, failed) result(pathTemp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathTemp
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)            , intent(in)    , optional  :: prefix, sep, ext
        integer(IK)                 , intent(in)    , optional  :: pid
        logical(LK)                 , intent(out)   , optional  :: failed
        character(:,SKG)            , allocatable               :: pathTemp
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the `path` that results from joining the two input path parts,
    !>  separated by the specified directory separator or the preferred directory separator of the runtime shell.<br>
    !>
    !>  \details
    !>  Each part is inspected for being an absolute path before being joined to previous part.<br>
    !>  If the second path part `p2` is absolute or relative with respect to the current drive,
    !>  then `pathJoined = p2` on output, otherwise `pathJoined = p1//getDirSep()//p2`.<br>
    !>
    !>  \param[in]      p1      :   The input scalar `character` of default kind \SK containing the first part of the path to be constructed.<br>
    !>  \param[in]      p2      :   The input scalar `character` of default kind \SK containing the second part of the path to be constructed.<br>
    !>  \param[out]     failed  :   The output `logical` of default kind \LK that is `.true.` if the procedure fails while inferring the runtime shell type and directory separator.<br>
    !>                              (**optional**. If missing and the procedure fails, then the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. It can be present <b>if and only if</b> `failed` is also present. If missing, no error message will be output.)
    !>
    !>  \return
    !>  `pathJoined`            :   The output `allocatable` scalar of type `character` of default kind \SK containing the joined path constructed from the input parts.<br>
    !>                              If the output argument `failed` is present, its value must be checked before using the output `pathJoined`.<br>
    !>                              In such a case, the output value of `pathJoined` is undefined and not guaranteed to be correct.<br>
    !>
    !>  \interface{getPathJoined}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathJoined
    !>
    !>      pathJoined = getPathJoined(p1, p2)
    !>      pathJoined = getPathJoined(p1, p2, failed, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getDirSep](@ref pm_sysPath::getDirSep)<br>
    !>  [hasDriveLetter](@ref pm_sysPath::hasDriveLetter)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [getPathPosixEscaped](@ref pm_sysPath::getPathPosixEscaped)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>
    !>  \example{getPathJoined}
    !>  \include{lineno} example/pm_sysPath/getPathJoined/main.F90
    !>  \compilef{getPathJoined}
    !>  \output{getPathJoined}
    !>  \include{lineno} example/pm_sysPath/getPathJoined/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{getPathJoined}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathJoined

    impure module function getPathJoined(p1, p2) result(pathJoined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathJoined
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)            , intent(in)                :: p1, p2
        character(:,SKG)            , allocatable               :: pathJoined
    end function

    impure module function getPathJoinedFailed(p1, p2, failed, errmsg) result(pathJoined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathJoinedFailed
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)            , intent(in)                :: p1, p2
        logical(LK)                 , intent(out)               :: failed
        character(*, SK)            , intent(inout) , optional  :: errmsg
        character(:,SKG)            , allocatable               :: pathJoined
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` is the input path is a file (not a directory), otherwise return `.false.`.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface utilize the Fortran intrinsic `inquire()`
    !>  statement to infer the existence of the specified input `path` as a system file.<br>
    !>
    !>  \param[in]      path    :   The input scalar `character` of kind \SKALL containing a POSIX-style or Windows-style path.<br>
    !>  \param[out]     iostat  :   The **optional** output scalar `integer` of default kind \IK.<br>
    !>                              <ol>
    !>                                  <li>    If present and no error occurs, it is set to `0` on output.<br>
    !>                                  <li>    If present and an error occurs, it is set to a positive non-zero value.<br>
    !>                                  <li>    If missing and an error occurs, then the program halts by calling `error stop` followed by the relevant error message.<br>
    !>                              </ol>
    !>                              (**optional**.)
    !>  \param[inout]   iomsg   :   The input/output scalar `character` of default kind \SK containing the error message, if any error occurs.<br>
    !>                              A length type parameter value of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) is generally sufficient for `iomsg` to contain the output error messages.<br>
    !>                              (**optional**. It can be present **only if** the `iostat` output argument is also present.)
    !>
    !>  \return
    !>  `pathIsFile`            :   The output scalar of type `logical` of default kind \LK that is `.true.` if the input path is a file, otherwise `.false.`.<br>
    !>
    !>  \interface{isFile}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_sysPath, only: isFile
    !>      logical(LK) :: pathIsFile
    !>
    !>      pathIsFile = isFile(path)
    !>      pathIsFile = isFile(path, iostat, iomsg = iomsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isOpen](@ref pm_io::isOpen)<br>
    !>  [isDir](@ref pm_sysPath::isDir)<br>
    !>  [isFile](@ref pm_sysPath::isFile)<br>
    !>  [isExtant](@ref pm_sysPath::isExtant)<br>
    !>  [isPreconnected](@ref pm_io::isPreconnected)<br>
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{isFile}
    !>  \include{lineno} example/pm_sysPath/isFile/main.F90
    !>  \compilef{isFile}
    !>  \output{isFile}
    !>  \include{lineno} example/pm_sysPath/isFile/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  This generic interface should be extended in the future to gracefully handle intrinsic `inquire()` exceptions.<br>
    !>
    !>  \final{isFile}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFile

    impure elemental module function isFileDD(path) result(pathIsFile)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFileDD
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    :: path
        logical(LK)                     :: pathIsFile
    end function

    impure elemental module function isFileII(path, iostat, iomsg) result(pathIsFile)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFileII
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                :: path
        logical(LK)                                 :: pathIsFile
        integer(IK)     , intent(out)               :: iostat
        character(*, SK), intent(inout) , optional  :: iomsg
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` is the input path is an extant system directory, otherwise return `.false.`.<br>
    !>
    !>  \details
    !>  The procedures under this generic interface utilize the Fortran intrinsic `inquire()`
    !>  statement to infer the existence of the specified input `path` as a system directory.<br>
    !>
    !>  \param[in]      path    :   The input scalar `character` of kind \SKALL containing a POSIX-style or Windows-style path.<br>
    !>  \param[out]     iostat  :   The **optional** output scalar `integer` of default kind \IK.<br>
    !>                              <ol>
    !>                                  <li>    If present and no error occurs, it is set to `0` on output.<br>
    !>                                  <li>    If present and an error occurs, it is set to a positive non-zero value.<br>
    !>                                  <li>    If missing and an error occurs, then the program halts by calling `error stop` followed by the relevant error message.<br>
    !>                              </ol>
    !>                              (**optional**.)
    !>  \param[inout]   iomsg   :   The input/output scalar `character` of default kind \SK containing the error message, if any error occurs.<br>
    !>                              A length type parameter value of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) is generally sufficient for `iomsg` to contain the output error messages.<br>
    !>                              (**optional**. It can be present **only if** the `iostat` output argument is also present.)
    !>
    !>  \return
    !>  `pathIsDir`             :   The output scalar of type `logical` of default kind \LK.<br>
    !>                              It is `.true.` if the input path is an extant system directory, otherwise `.false.`.<br>
    !>
    !>  \interface{isDir}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_sysPath, only: isDir
    !>      logical(LK) :: pathIsDir
    !>
    !>      pathIsDir = isDir(path)
    !>      pathIsDir = isDir(path, iostat, iomsg = iomsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [isOpen](@ref pm_io::isOpen)<br>
    !>  [isDir](@ref pm_sysPath::isDir)<br>
    !>  [isFile](@ref pm_sysPath::isFile)<br>
    !>  [isExtant](@ref pm_sysPath::isExtant)<br>
    !>  [isPreconnected](@ref pm_io::isPreconnected)<br>
    !>  [getPathNew](@ref pm_sysPath::getPathNew)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathWindows](@ref pm_sysPath::getPathWindows)<br>
    !>  [setPathWindows](@ref pm_sysPath::setPathWindows)<br>
    !>  [getPathHostNameIndex](@ref pm_sysPath::getPathHostNameIndex)<br>
    !>
    !>  \example{isDir}
    !>  \include{lineno} example/pm_sysPath/isDir/main.F90
    !>  \compilef{isDir}
    !>  \output{isDir}
    !>  \include{lineno} example/pm_sysPath/isDir/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  This generic interface should be extended in the future to gracefully handle intrinsic `inquire()` exceptions.<br>
    !>
    !>  \final{isDir}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:09 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isDir

    impure elemental module function isDirDD(path) result(pathIsDir)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDirDD
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                :: path
        logical(LK)                                 :: pathIsDir
    end function

    impure elemental module function isDirII(path, iostat, iomsg) result(pathIsDir)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDirII
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                :: path
        integer(IK)     , intent(out)               :: iostat
        character(*, SK), intent(inout) , optional  :: iomsg
        logical(LK)                                 :: pathIsDir
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input `path` (whether a file or a directory) exists and `.false.` otherwise.<br>
    !>
    !>  \details
    !>  This generic functional interface is a simple convenience wrapper around the
    !>  generic interfaces `isFile()` and `isDir()` and returns the short-circuited result of `isFile() .or. isDir()`.<br>
    !>  As such, its behavior is fully dictated by these generic interfaces.<br>
    !>
    !>  \param[in]      path    :   The input scalar, or array of arbitrary rank, of type `character` of default kind \SK representing the path whose existence is to be checked.<br>
    !>  \param[out]     iostat  :   The **optional** output scalar `integer` of default kind \IK.<br>
    !>                              <ol>
    !>                                  <li>    If present and no error occurs, it is set to `0` on output.<br>
    !>                                  <li>    If present and an error occurs, it is set to a positive non-zero value.<br>
    !>                                  <li>    If missing and an error occurs, then the program halts by calling `error stop` followed by the relevant error message.<br>
    !>                              </ol>
    !>                              (**optional**.)
    !>  \param[inout]   iomsg   :   The input/output scalar `character` of default kind \SK containing the error message, if any error occurs.<br>
    !>                              A length type parameter value of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) is generally sufficient for `iomsg` to contain the output error messages.<br>
    !>                              (**optional**. It can be present **only if** the `iostat` output argument is also present.)
    !>
    !>  \return
    !>  `extant`                :   The output scalar, or array of same rank, shape, and size as the input array-like argument, of type `logical` of default kind \LK.<br>
    !>                              It is `.true.` <b>if and only if</b> the input `path` (whether a file or a directory) exists.<br>
    !>                              Otherwise, it is `.false.`.<br>
    !>
    !>  \interface{isExtant}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_sysPath, only: isExtant
    !>      integer(IK) :: iostat
    !>      logical(LK) :: extant
    !>
    !>      extant = isExtant(path)
    !>      extant = isExtant(path, iostat, iomsg = iomsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isOpen](@ref pm_io::isOpen)<br>
    !>  [isDir](@ref pm_sysPath::isDir)<br>
    !>  [isFile](@ref pm_sysPath::isFile)<br>
    !>  [isExtant](@ref pm_sysPath::isExtant)<br>
    !>  [isPreconnected](@ref pm_io::isPreconnected)<br>
    !>
    !>  \example{isExtant}
    !>  \include{lineno} example/pm_sysPath/isExtant/main.F90
    !>  \compilef{isExtant}
    !>  \output{isExtant}
    !>  \include{lineno} example/pm_sysPath/isExtant/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final{isExtant}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isExtant

    impure elemental module function isExtantDD(path) result(extant)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isExtantDD
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                :: path
        logical(LK)                                 :: extant
    end function

    impure elemental module function isExtantII(path, iostat, iomsg) result(extant)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isExtantII
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)                :: path
        integer(IK)     , intent(out)               :: iostat
        character(*, SK), intent(inout) , optional  :: iomsg
        logical(LK)                                 :: extant
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a path quoted with double quotation marks
    !>  where all instances of double quotes \f$\ms{"}\f$ within the path are removed.<br>
    !>
    !>  \details
    !>  This procedure exists to minimize security vulnerabilities associated with paths that may contain harmful system commands.<br>
    !>  To ensure paths are interpreted as simple strings, they must be quoted with proper quotation marks on both Windows and Unix-like platforms.<br>
    !>  The path returned by this procedure is as close as one can get to a verbatim path in Windows CMD.<br>
    !>  Windows-style paths **cannot contain** double quotation marks. However, their dangling presence can create security vulnerabilities.<br>
    !>  As such, if a path is be used in a call to a CMD runtime shell, it is recommended to pass it to this procedure
    !>  to ensure the path is quoted and all instances of dangling quotation marks \f$\ms{"}\f$ are gracefully escaped.<br>
    !>  See also the warnings below.<br>
    !>
    !>  **When should I quote a Windows CMD path using this generic interface?**<br>
    !>  Quoting CMD paths is crucial for graceful handling of situations where,
    !>  -#  the path is too long (`> 255` characters), or
    !>  -#  the path contains blanks or dangling double-quotation marks, or
    !>  -#  the path mixes forward slash and backward slash as directory separators.<br>
    !>
    !>  \param[in]      path        :   The input scalar `character` of kind \SKALL containing the path to be enclosed with quotation marks.<br>
    !>
    !>  \return
    !>  `pathVerbatim`              :   The output `allocatable` scalar of type `character` of the same kind as `path` containing the modified input `path`
    !>                                  enclosed with double quotation marks where all instances of double quotation marks within the path are removed.<br>
    !>
    !>  \interface{getPathVerbatimCMD}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathVerbatimCMD
    !>
    !>      pathVerbatim = getPathVerbatimCMD(path)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that this procedure does not inspect the input path for the presence of other illegal Windows reserved characters (e.g., [WINDOWS_RESERVED_STR](@ref pm_sysPath::WINDOWS_RESERVED_STR)).<br>
    !>  The presence of illegal characters in a Windows path can lead to runtime errors even if the path is quoted by the procedures of this generic interface.<br>
    !>
    !>  \warning
    !>  The Windows CMD wildcard characters remain interpretable by the CMD shell if any are present in the input path.<br>
    !>  The asterisk `*` and question mark `?` are used as wildcard characters in Windows CMD, as they are in MS-DOS and Windows.<br>
    !>  The **asterisk** matches **any sequence of characters**, whereas the **question mark** matches **any single character**.<br>
    !>
    !>  \warnpure
    !>  The impurity is caused by the call to procedures whose purity depends on the library build configuration.<br>
    !>
    !>  \see
    !>  [getPathVerbatimPowerShell](@ref pm_sysPath::getPathVerbatimPowerShell)<br>
    !>  [getPathVerbatimPosix](@ref pm_sysPath::getPathVerbatimPosix)<br>
    !>  [getPathVerbatimFish](@ref pm_sysPath::getPathVerbatimFish)<br>
    !>  [getPathVerbatimCMD](@ref pm_sysPath::getPathVerbatimCMD)<br>
    !>  [getPathVerbatim](@ref pm_sysPath::getPathVerbatim)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathPosixEscaped](@ref pm_sysPath::getPathPosixEscaped)<br>
    !>  [setPathPosixEscaped](@ref pm_sysPath::setPathPosixEscaped)<br>
    !>
    !>  \example{getPathVerbatimCMD}
    !>  \include{lineno} example/pm_sysPath/getPathVerbatimCMD/main.F90
    !>  \compilef{getPathVerbatimCMD}
    !>  \output{getPathVerbatimCMD}
    !>  \include{lineno} example/pm_sysPath/getPathVerbatimCMD/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final{getPathVerbatimCMD}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathVerbatimCMD

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPathVerbatimCMD_SK5(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimCMD_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK4_ENABLED
    PURE module function getPathVerbatimCMD_SK4(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimCMD_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK3_ENABLED
    PURE module function getPathVerbatimCMD_SK3(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimCMD_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK2_ENABLED
    PURE module function getPathVerbatimCMD_SK2(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimCMD_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK1_ENABLED
    PURE module function getPathVerbatimCMD_SK1(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimCMD_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a path quoted with single quotation marks
    !>  where all instances of single quotes \f$\ms{'}\f$ within the path are replaced with \f$\ms{''}\f$.<br>
    !>
    !>  \details
    !>  This procedure exists to minimize security vulnerabilities associated with paths that may contain harmful system commands.<br>
    !>  To ensure paths are interpreted as simple strings, they must be quoted with proper quotation marks on both Windows and Unix-like platforms.<br>
    !>  The path returned by this procedure is as close as one can get to a verbatim path in Windows PowerShell or PowerShell Core on Unix platforms.<br>
    !>  If a path is be used in a call to a PowerShell runtime shell, it is recommended to pass it to this procedure
    !>  to ensure the path is quoted and all instances of dangling quotation marks \f$\ms{'}\f$ are gracefully escaped.<br>
    !>
    !>  \param[in]      path        :   The input scalar `character` of kind \SKALL containing the path to be enclosed with quotation marks.<br>
    !>
    !>  \return
    !>  `pathVerbatim`              :   The output `allocatable` scalar of type `character` of the same kind as `path` containing the modified input `path`
    !>                                  enclosed with single quotation marks where all instances of single quotation marks within the path are escaped.<br>
    !>
    !>  \interface{getPathVerbatimPowerShell}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathVerbatimPowerShell
    !>
    !>      pathVerbatim = getPathVerbatimPowerShell(path)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Note that this procedure does not inspect the input path for the presence of other illegal Windows reserved characters on Windows platforms (e.g., [WINDOWS_RESERVED_STR](@ref pm_sysPath::WINDOWS_RESERVED_STR)).<br>
    !>  The presence of illegal characters in a Windows path can lead to runtime errors even if the path is quoted by the procedures of this generic interface.<br>
    !>
    !>  \warnpure
    !>  The impurity is caused by the call to procedures whose purity depends on the library build configuration.<br>
    !>
    !>  \see
    !>  [getPathVerbatimPowerShell](@ref pm_sysPath::getPathVerbatimPowerShell)<br>
    !>  [getPathVerbatimPosix](@ref pm_sysPath::getPathVerbatimPosix)<br>
    !>  [getPathVerbatimFish](@ref pm_sysPath::getPathVerbatimFish)<br>
    !>  [getPathVerbatimCMD](@ref pm_sysPath::getPathVerbatimCMD)<br>
    !>  [getPathVerbatim](@ref pm_sysPath::getPathVerbatim)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathPosixEscaped](@ref pm_sysPath::getPathPosixEscaped)<br>
    !>  [setPathPosixEscaped](@ref pm_sysPath::setPathPosixEscaped)<br>
    !>
    !>  \example{getPathVerbatimPowerShell}
    !>  \include{lineno} example/pm_sysPath/getPathVerbatimPowerShell/main.F90
    !>  \compilef{getPathVerbatimPowerShell}
    !>  \output{getPathVerbatimPowerShell}
    !>  \include{lineno} example/pm_sysPath/getPathVerbatimPowerShell/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final{getPathVerbatimPowerShell}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathVerbatimPowerShell

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPathVerbatimPowerShell_SK5(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimPowerShell_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK4_ENABLED
    PURE module function getPathVerbatimPowerShell_SK4(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimPowerShell_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK3_ENABLED
    PURE module function getPathVerbatimPowerShell_SK3(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimPowerShell_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK2_ENABLED
    PURE module function getPathVerbatimPowerShell_SK2(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimPowerShell_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK1_ENABLED
    PURE module function getPathVerbatimPowerShell_SK1(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimPowerShell_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a path quoted with single quotation marks
    !>  where all instances of single quotes \f$\ms{'}\f$ within the path are replaced with \f$\ms{'\''}\f$.<br>
    !>
    !>  \details
    !>  **POSIX-style paths** can contain almost any character including single or double quotation marks.<br>
    !>  To simplify the task and minimize the risks, this procedure quotes all input POSIX style paths with **single** quotation marks.<br>
    !>  Enclosing the path with single quotation marks requires escaping any existing instances of single quotes within the path.<br>
    !>  This is done by replacing all instances of \f$\ms{'}\f$ within the path are replaced with \f$\ms{'\''}\f$.<br>
    !>
    !>  Optionally, if a set of wildcard characters are also specified, all instances of such characters identified
    !>  within the specified `path` are properly excluded from the quotations, such that the resulting output
    !>  path can be used as pattern for globing.<br>
    !>
    !>  \param[in]      path        :   The input scalar `character` of kind \SKALL containing the path to be enclosed with quotation marks.<br>
    !>
    !>  \return
    !>  `pathVerbatim`              :   The output `allocatable` scalar of type `character` of the same kind as `path` containing the modified input `path`
    !>                                  enclosed with single quotation marks where all instances of single quotation marks within the path are escaped.<br>
    !>                                  Wildcard characters, if any are present, are not quoted via single quotation marks.<br>
    !>
    !>  \interface{getPathVerbatimPosix}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathVerbatimPosix
    !>
    !>      pathVerbatim = getPathVerbatimPosix(path)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>  The impurity is caused by the call to procedures whose purity depends on the library build configuration.<br>
    !>
    !>  \see
    !>  [getPathVerbatimPowerShell](@ref pm_sysPath::getPathVerbatimPowerShell)<br>
    !>  [getPathVerbatimPosix](@ref pm_sysPath::getPathVerbatimPosix)<br>
    !>  [getPathVerbatimFish](@ref pm_sysPath::getPathVerbatimFish)<br>
    !>  [getPathVerbatimCMD](@ref pm_sysPath::getPathVerbatimCMD)<br>
    !>  [getPathVerbatim](@ref pm_sysPath::getPathVerbatim)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathPosixEscaped](@ref pm_sysPath::getPathPosixEscaped)<br>
    !>  [setPathPosixEscaped](@ref pm_sysPath::setPathPosixEscaped)<br>
    !>
    !>  \example{getPathVerbatimPosix}
    !>  \include{lineno} example/pm_sysPath/getPathVerbatimPosix/main.F90
    !>  \compilef{getPathVerbatimPosix}
    !>  \output{getPathVerbatimPosix}
    !>  \include{lineno} example/pm_sysPath/getPathVerbatimPosix/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final{getPathVerbatimPosix}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathVerbatimPosix

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPathVerbatimPosix_SK5(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimPosix_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK4_ENABLED
    PURE module function getPathVerbatimPosix_SK4(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimPosix_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK3_ENABLED
    PURE module function getPathVerbatimPosix_SK3(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimPosix_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK2_ENABLED
    PURE module function getPathVerbatimPosix_SK2(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimPosix_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK1_ENABLED
    PURE module function getPathVerbatimPosix_SK1(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimPosix_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a path quoted with single quotation marks where all instances of backslashes within the path are replaced
    !>  with double backslash and all single quotes \f$\ms{'}\f$ are replaced with \f$\ms{\'}\f$.<br>
    !>
    !>  \details
    !>  This procedure exists to minimize security vulnerabilities associated with paths that may contain harmful system commands.<br>
    !>  To ensure paths are interpreted as simple strings in Fish shell, they must be quoted with proper quotation marks on all platforms.<br>
    !>  The path returned by this procedure is as close as one can get to a verbatim path in Fish shell on any platform.<br>
    !>  If a path is be used in a call to a PowerShell runtime shell, it is recommended to pass it to this procedure
    !>  to ensure the path is quoted and all instances of dangling quotation marks \f$\ms{'}\f$are gracefully escaped.<br>
    !>
    !>  Optionally, if a set of wildcard characters are also specified, all instances of such characters identified
    !>  within the specified `path` are properly excluded from the quotations, such that the resulting output
    !>  path can be used as pattern for globing.<br>
    !>
    !>  \param[in]      path        :   The input scalar `character` of kind \SKALL containing the path to be enclosed with quotation marks.<br>
    !>
    !>  \return
    !>  `pathVerbatim`              :   The output `allocatable` scalar of type `character` of the same kind as `path` containing the modified input `path`
    !>                                  enclosed with single quotation marks where all instances of single quotation marks and backslashes within the path are escaped.<br>
    !>                                  Wildcard characters, if any are present, are not quoted via single quotation marks.<br>
    !>
    !>  \interface{getPathVerbatimFish}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathVerbatimFish
    !>
    !>      pathVerbatim = getPathVerbatimFish(path)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>  The impurity is caused by the call to procedures whose purity depends on the library build configuration.<br>
    !>
    !>  \see
    !>  [getPathVerbatimPowerShell](@ref pm_sysPath::getPathVerbatimPowerShell)<br>
    !>  [getPathVerbatimPosix](@ref pm_sysPath::getPathVerbatimPosix)<br>
    !>  [getPathVerbatimFish](@ref pm_sysPath::getPathVerbatimFish)<br>
    !>  [getPathVerbatimCMD](@ref pm_sysPath::getPathVerbatimCMD)<br>
    !>  [getPathVerbatim](@ref pm_sysPath::getPathVerbatim)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathPosixEscaped](@ref pm_sysPath::getPathPosixEscaped)<br>
    !>  [setPathPosixEscaped](@ref pm_sysPath::setPathPosixEscaped)<br>
    !>
    !>  \example{getPathVerbatimFish}
    !>  \include{lineno} example/pm_sysPath/getPathVerbatimFish/main.F90
    !>  \compilef{getPathVerbatimFish}
    !>  \output{getPathVerbatimFish}
    !>  \include{lineno} example/pm_sysPath/getPathVerbatimFish/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \final{getPathVerbatimFish}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathVerbatimFish

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getPathVerbatimFish_SK5(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimFish_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK4_ENABLED
    PURE module function getPathVerbatimFish_SK4(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimFish_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK3_ENABLED
    PURE module function getPathVerbatimFish_SK3(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimFish_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK2_ENABLED
    PURE module function getPathVerbatimFish_SK2(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimFish_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

#if SK1_ENABLED
    PURE module function getPathVerbatimFish_SK1(path) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatimFish_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)    , intent(in)                :: path
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a path quoted with single or double quotation marks
    !>  whose metacharacters are also properly escaped depending on the runtime system shell.<br>
    !>
    !>  \details
    !>  This procedure exists to minimize security vulnerabilities associated with paths that may contain harmful system commands.<br>
    !>  To ensure paths are interpreted as simple strings, they must be quoted with proper quotation marks on both Windows and Unix-like platforms.<br>
    !>  The procedures under this generic interface call one of the following generic interfaces, depending on the specified or determined runtime system shell,<br>
    !>  -#  [getPathVerbatimPowerShell](@ref pm_sysPath::getPathVerbatimPowerShell) in Windows or Unix PowerShell environments.<br>
    !>  -#  [getPathVerbatimPosix](@ref pm_sysPath::getPathVerbatimPosix) in POSIX-compliant shell environments.<br>
    !>  -#  [getPathVerbatimCMD](@ref pm_sysPath::getPathVerbatimCMD) in Windows CMD shell environments.<br>
    !>
    !>  \param[in]      path        :   The input scalar `character` of default kind \SK containing the path to be enclosed with quotation marks.<br>
    !>  \param[out]     failed      :   The output scalar `logical` of default kind \LK that is `.true.`
    !>                                  <b>if and only if</b> an error occurs while inferring the runtime system shell.<br>
    !>                                  (**optional**. If missing and a runtime error occurs, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg      :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                                  If an error occurs, `errmsg` will be set to a descriptive message about the nature of the runtime error.<br>
    !>                                  A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                                  (**optional**. If missing, no error message will be output. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `pathVerbatim`              :   The output `allocatable` scalar of type `character` of the same kind as `path` containing the modified input `path` enclosed
    !>                                  with single or double quotation marks such that it can be treated as a verbatim path in the specified or the runtime shell.<br>
    !>
    !>  \interface{getPathVerbatim}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathVerbatim
    !>      character(255, SK) :: errmsg
    !>      logical(LK) :: failed
    !>
    !>      pathVerbatim = getPathVerbatim(path, failed = failed, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \see
    !>  [getPathVerbatimPowerShell](@ref pm_sysPath::getPathVerbatimPowerShell)<br>
    !>  [getPathVerbatimPosix](@ref pm_sysPath::getPathVerbatimPosix)<br>
    !>  [getPathVerbatimCMD](@ref pm_sysPath::getPathVerbatimCMD)<br>
    !>  [getPathVerbatim](@ref pm_sysPath::getPathVerbatim)<br>
    !>  [getPathPosix](@ref pm_sysPath::getPathPosix)<br>
    !>  [setPathPosix](@ref pm_sysPath::setPathPosix)<br>
    !>  [getPathPosixEscaped](@ref pm_sysPath::getPathPosixEscaped)<br>
    !>  [setPathPosixEscaped](@ref pm_sysPath::setPathPosixEscaped)<br>
    !>
    !>  \example{getPathVerbatim}
    !>  \include{lineno} example/pm_sysPath/getPathVerbatim/main.F90
    !>  \compilef{getPathVerbatim}
    !>  \output{getPathVerbatim}
    !>  \include{lineno} example/pm_sysPath/getPathVerbatim/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{getPathVerbatim}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathVerbatim

    impure module function getPathVerbatim(path, failed, errmsg) result(pathVerbatim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathVerbatim
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)                :: path
        logical(LK)         , intent(out)   , optional  :: failed
        character(*, SK)    , intent(inout) , optional  :: errmsg
        character(:,SKG)    , allocatable               :: pathVerbatim
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the index of the **last** character of the **dirname** (directory part) of the input `path`.<br>
    !>
    !>  \details
    !>  Conventionally, the **dirname** segment of a path corresponds to the part that starts
    !>  from the beginning of the path and ends with the character **before the last directory separator** in the path.<br>
    !>  The behavior is consistent with the Unix `dirname` software the last directory separator from the output `dirname`.<br>
    !>  For example, `getIndexDirName("./paramonte", "/")` yields `1` as per the `dirname` software behavior.<br>
    !>  Also, `getIndexDirName("/", "/")` yields `1`, consistent with the behavior of Unix `dirname`.<br>
    !>  Note that `getIndexDirName("", "/")` yields `0`, meaning that the `dirname` is `.`, again consistent with the behavior of Unix `dirname`.<br>
    !>
    !>  This index can be also readily obtained via the Fortran intrinsic `index()` as `index(path, dirsep, back = .true.)`.<br>
    !>  The start index of the dirname segment of the path is trivial (`1`) and not returned as the output of this procedure.<br>
    !>  For example, if `path = "./paramonte"`, then `path(1:getIndexDirName(path,"/"))` yields `.`.<br>
    !>
    !>  Alternatively, however, one may be interested in extracting the dirname of a path **verbatim**
    !>  while keeping the trailing directory separator character in the output `dirname`.<br>
    !>  This leads to the more intuitive behavior of yielding `1` as the result of `getIndexDirName("/", "/")` (similar to POSIX `dirname` software)
    !>  while yielding `10` as the result of `getIndexDirName("paramonte/", "/")` (unlike the POSIX `dirname` software which returns `0` corresponding to `"."`)
    !>  In order to also allow this interpretation and behavior, the procedures of this generic interface implement this latter style by setting the `style` argument
    !>  of this generic interface to an object of type [verbatim_type](@ref pm_sysPath::verbatim_type) such as the constant [verbatim](@ref pm_sysPath::verbatim).<br>
    !>  For more details on this alternative behavior, see the examples below.<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the path.<br>
    !>  \param[in]  dirsep  :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the runtime directory separator(s).<br>
    !>                          The directory separator can be obtained from either the `dirsep` component of an object of type
    !>                          [shell_type](@ref pm_sysShell::shell_type) or [getDirSep](@ref pm_sysPath::getDirSep).<br>
    !>                          The directory separator is `\` or `/` in Windows-based terminals (e.g., CMD or PowerShell) and `/` in POSIX-compliant shells.<br>
    !>                          When in doubt (for example, in Windows terminals), `dirsep` can be set to multiple characters, for example `dirsep = "/\"`.<br>
    !>                          In such a case, the input `path` will be scanned for the presence of any of the individual characters in `dirsep`.<br>
    !>  \param[in]  style   :   The input scalar that can be,<br>
    !>                          <ol>
    !>                              <li>    the constant [verbatim](@ref pm_sysPath::verbatim) or an object of type [verbatim_type](@ref pm_sysPath::verbatim_type)
    !>                                      implying the use of the ParaMonte style described above in extracting the output `dirname`.<br>
    !>                          </ol>
    !>                          (**optional**. If missing, the default behavior corresponds to that of the `dirname` command on POSIX systems. See examples below.)
    !>
    !>  \return
    !>  `index`             :   The output scalar `integer` of default kind \IK containing the position of the last directory separator in the input `path`.<br>
    !>                          If there is no directory separator in the `path`, the output value is `0`.<br>
    !>
    !>  \interface{getIndexDirName}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getIndexDirName, verbatim
    !>      integer(IK) :: index
    !>
    !>      index = getIndexDirName(path, dirsep)
    !>      index = getIndexDirName(path, dirsep, style)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < len(dirsep)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>  [getExtName](@ref pm_sysPath::getExtName)<br>
    !>  [getFileName](@ref pm_sysPath::getFileName)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getIndexDirName](@ref pm_sysPath::getIndexDirName)<br>
    !>  [getIndexExtName](@ref pm_sysPath::getIndexExtName)<br>
    !>  [getIndexBaseName](@ref pm_sysPath::getIndexBaseName)<br>
    !>
    !>  \example{getIndexDirName}
    !>  \include{lineno} example/pm_sysPath/getIndexDirName/main.F90
    !>  \compilef{getIndexDirName}
    !>  \output{getIndexDirName}
    !>  \include{lineno} example/pm_sysPath/getIndexDirName/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  This procedure should be extended to support non-default character kinds.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  The examples should be extended to cover the optional argument `style`.<br>
    !>
    !>  \final{getIndexDirName}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getIndexDirName

    PURE elemental module function getIndexDirNameDef(path, dirsep) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIndexDirNameDef
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)    :: path
        character(*,SKG)    , intent(in)    :: dirsep
        integer(IK)                         :: index
    end function

    PURE elemental module function getIndexDirNamePM(path, dirsep, style) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIndexDirNamePM
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)    :: path
        character(*,SKG)    , intent(in)    :: dirsep
        type(verbatim_type) , intent(in)    :: style
        integer(IK)                         :: index
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the index of the **first** character of the **basename** part of the input `path`.<br>
    !>
    !>  \details
    !>  By definition, the basename of a path is the segment that appears after the last directory separator in the path.<br>
    !>  As such, the basename of a path can be either empty, a directory, or correspond to a filename (including its extension).<br>
    !>  The end index of the basename segment of the path is trivial (`len(path)`) and not returned as the output of this procedure.<br>
    !>  For example, if `path = "./paramonte"`, then `path(getIndexBaseName(path,"/"):)` yields `paramonte`.<br>
    !>
    !>  See also the details of the documentation of the generic interface [getIndexDirName](@ref pm_sysPath::getIndexDirName) and [getDirName](@ref pm_sysPath::getDirName).<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of default kind \SK containing the path.<br>
    !>  \param[in]  dirsep  :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the runtime directory separator(s).<br>
    !>                          The directory separator can be obtained from either the `dirsep` component of an object of type
    !>                          [shell_type](@ref pm_sysShell::shell_type) or [getDirSep](@ref pm_sysPath::getDirSep).<br>
    !>                          The directory separator is `\` or `/` in Windows-based terminals (e.g., CMD or PowerShell) and `/` in POSIX-compliant shells.<br>
    !>                          When in doubt (for example, in Windows terminals), `dirsep` can be set to multiple characters, for example `dirsep = "/\"`.<br>
    !>                          In such a case, the input `path` will be scanned for the presence of any of the individual characters in `dirsep`.<br>
    !>  \param[in]  style   :   The input scalar that can be,<br>
    !>                          <ol>
    !>                              <li>    the constant [verbatim](@ref pm_sysPath::verbatim) or an object of type [verbatim_type](@ref pm_sysPath::verbatim_type)
    !>                                      implying the use of the ParaMonte style described above in extracting the output `basename`.<br>
    !>                          </ol>
    !>                          (**optional**. If missing, the default behavior corresponds to that of the `basename` command on POSIX systems. See examples below.)
    !>
    !>  \return
    !>  `index`             :   The output scalar `integer` of default kind \IK containing the starting position of the basename of the input `path`.<br>
    !>                          If the input `path` does not contain a basename, the output value for `index` will be `len(path) + 1`.<br>
    !>
    !>  \interface{getIndexBaseName}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getIndexBaseName, verbatim
    !>      integer(IK) :: index
    !>
    !>      index = getIndexBaseName(path, dirsep)
    !>      index = getIndexBaseName(path, dirsep, verbatim)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < len(dirsep)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>  [getExtName](@ref pm_sysPath::getExtName)<br>
    !>  [getFileName](@ref pm_sysPath::getFileName)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getIndexDirName](@ref pm_sysPath::getIndexDirName)<br>
    !>  [getIndexExtName](@ref pm_sysPath::getIndexExtName)<br>
    !>  [getIndexBaseName](@ref pm_sysPath::getIndexBaseName)<br>
    !>
    !>  \example{getIndexBaseName}
    !>  \include{lineno} example/pm_sysPath/getIndexBaseName/main.F90
    !>  \compilef{getIndexBaseName}
    !>  \output{getIndexBaseName}
    !>  \include{lineno} example/pm_sysPath/getIndexBaseName/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  This procedure should be extended to support non-default character kinds.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  The examples should be extended to cover the optional argument `style`.<br>
    !>
    !>  \final{getIndexBaseName}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getIndexBaseName

    PURE elemental module function getIndexBaseNameDef(path, dirsep) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIndexBaseNameDef
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)    :: path
        character(*,SKG)    , intent(in)    :: dirsep
        integer(IK)                         :: index
    end function

    PURE elemental module function getIndexBaseNamePM(path, dirsep, style) result(index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIndexBaseNamePM
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)    :: path
        character(*,SKG)    , intent(in)    :: dirsep
        type(verbatim_type) , intent(in)    :: style
        integer(IK)                         :: index
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the index of the **first** character of the **file extension** of the input `path` (if any extension exists).<br>
    !>
    !>  \details
    !>  By definition, the file extension is the last segment of a path that starts with the last dot `.` character in the path.<br>
    !>  As such, the file extension of a path can be either empty, a directory, or correspond to a true file extension.<br>
    !>  The end index of the file extension segment of the path is trivial (`len(path)`) and not returned as the output of this procedure.<br>
    !>  For example, if `path = "./paramonte.txt"`, then `path(getIndexExtName(path,"/"):)` yields `.txt`.<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of default kind \SK containing the path.<br>
    !>  \param[in]  dirsep  :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the runtime directory separator(s).<br>
    !>                          The directory separator can be obtained from either the `dirsep` component of an object of type
    !>                          [shell_type](@ref pm_sysShell::shell_type) or [getDirSep](@ref pm_sysPath::getDirSep).<br>
    !>                          The directory separator is `\` or `/` in Windows-based terminals (e.g., CMD or PowerShell) and `/` in POSIX-compliant shells.<br>
    !>                          When in doubt (for example, in Windows terminals), `dirsep` can be set to multiple characters, for example `dirsep = "/\"`.<br>
    !>                          In such a case, the input `path` will be scanned for the presence of any of the individual characters in `dirsep`.<br>
    !>
    !>  \return
    !>  `indexExtName`      :   The output scalar `integer` of default kind \IK containing the starting position of the file extension segment of the input `path`.<br>
    !>                          If the input `path` does not contain a file extension, the output value for `indexExtName` will be `len(path) + 1`.<br>
    !>
    !>  \interface{getIndexExtName}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getIndexExtName
    !>      integer(IK) :: indexExtName
    !>
    !>      indexExtName = getIndexExtName(path, dirsep)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < len(dirsep)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>  [getExtName](@ref pm_sysPath::getExtName)<br>
    !>  [getFileName](@ref pm_sysPath::getFileName)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getIndexDirName](@ref pm_sysPath::getIndexDirName)<br>
    !>  [getIndexExtName](@ref pm_sysPath::getIndexExtName)<br>
    !>  [getIndexBaseName](@ref pm_sysPath::getIndexBaseName)<br>
    !>
    !>  \example{getIndexExtName}
    !>  \include{lineno} example/pm_sysPath/getIndexExtName/main.F90
    !>  \compilef{getIndexExtName}
    !>  \output{getIndexExtName}
    !>  \include{lineno} example/pm_sysPath/getIndexExtName/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  This procedure should be extended to support non-default character kinds.<br>
    !>
    !>  \final{getIndexExtName}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getIndexExtName
    PURE elemental module function getIndexExtName(path, dirsep) result(indexExtName)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIndexExtName
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    :: path
        character(*,SKG), intent(in)    :: dirsep
        integer(IK)                     :: indexExtName
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **dirname** (directory name) part of the input `path`.<br>
    !>
    !>  \details
    !>  By definition, the **dirname** segment of a path corresponds to the part that starts
    !>  from the beginning of the path and ends with the character **before the last directory separator** in the path.<br>
    !>  The behavior is consistent with the Unix `dirname` software the last directory separator from the output `dirname`.<br>
    !>  For example, `getDirName("./paramonte", "/")` yields `"."` similar to the `dirname` software behavior.<br>
    !>  Also, `getDirName("/", "/")` yields `"/"`, again similar to the behavior of Unix `dirname`.<br>
    !>
    !>  \note
    !>  Counterintuitively, the UNIX `dirname` command treats paths that end in a trailing slash **not** as a directory.<br>
    !>  The trailing slash represents all files within the directory.<br>
    !>  For example according to the POSIX `dirname` implementation,
    !>  \code{.sh}
    !>      /home/martin/docs/
    !>  \endcode
    !>  represents the listing of all files within the specified directory.<br>
    !>  To point to the directory, one has to add a `.` character after the final directory separator character `/` as in,
    !>  \code{.sh}
    !>      /home/martin/docs/.
    !>  \endcode
    !>  By default, this generic interface follows the behavior implemented in POSIX `dirname`.<br>
    !>  Alternatively, however, one may be interested in extracting the dirname of a path **verbatim**
    !>  while keeping the trailing directory separator character in the output `dirname`.<br>
    !>  For example, one may prefer to treat paths ending with a directory separator characters
    !>  (e.g., [DIR_SEP_POSIX_ALL](@ref pm_sysPath::DIR_SEP_POSIX_ALL) or [DIR_SEP_WINDOWS_ALL](@ref pm_sysPath::DIR_SEP_WINDOWS_ALL))
    !>  as full directory names followed by an empty basename.<br>
    !>  This leads to the more intuitive behavior of yielding `"/"` as the result of `getIndexDirName("/", "/")` (similar to POSIX `dirname` software)
    !>  while yielding `"paramonte"` as the result of `getIndexDirName("paramonte/", "/")` (unlike the POSIX `dirname` software which returns `"."`)
    !>  In order to also allow this interpretation and behavior, the procedures of this generic interface implement this latter style by setting the `style` argument
    !>  of this generic interface to an object of type [verbatim_type](@ref pm_sysPath::verbatim_type) such as the constant [verbatim](@ref pm_sysPath::verbatim).<br>
    !>  For more details on this alternative behavior, see the examples below.<br>
    !>
    !>  See also the details of the documentation of the generic interface [getIndexDirName](@ref pm_sysPath::getIndexDirName).<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the path to split.<br>
    !>  \param[in]  dirsep  :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the runtime directory separator(s).<br>
    !>                          The directory separator can be obtained from either the `dirsep` component of an object of type
    !>                          [shell_type](@ref pm_sysShell::shell_type) or [getDirSep](@ref pm_sysPath::getDirSep).<br>
    !>                          The directory separator is `\` or `/` in Windows-based terminals (e.g., CMD or PowerShell) and `/` in POSIX-compliant shells.<br>
    !>                          When in doubt (for example, in Windows terminals), `dirsep` can be set to multiple characters, for example `dirsep = "/\"`.<br>
    !>                          In such a case, the input `path` will be scanned for the presence of any of the individual characters in `dirsep`.<br>
    !>  \param[in]  style   :   The input scalar that can be,<br>
    !>                          <ol>
    !>                              <li>    the constant [verbatim](@ref pm_sysPath::verbatim) or an object of type [verbatim_type](@ref pm_sysPath::verbatim_type)
    !>                                      implying the use of the ParaMonte style described above in extracting the output `dirname`.<br>
    !>                          </ol>
    !>                          (**optional**. If missing, the default behavior corresponds to that of the `dirname` command on POSIX systems. See examples below.)
    !>
    !>  \return
    !>  `dirname`           :   The output `allocatable` scalar `character` of default kind \SK containing the `dirname` segment of the input `path`.<br>
    !>
    !>  \interface{getDirName}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getDirName, verbatim
    !>      character(:, SK), allocatable :: dirname
    !>
    !>      dirname = getDirName(path, dirsep)
    !>      dirname = getDirName(path, dirsep, style)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < len(dirsep)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  See [getIndexDirName](@ref pm_sysPath::getIndexDirName) for a faster
    !>  method of getting the dirname segment of a path that avoids memory allocation.<br>
    !>
    !>  \see
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>  [getExtName](@ref pm_sysPath::getExtName)<br>
    !>  [getFileName](@ref pm_sysPath::getFileName)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getIndexDirName](@ref pm_sysPath::getIndexDirName)<br>
    !>  [getIndexExtName](@ref pm_sysPath::getIndexExtName)<br>
    !>  [getIndexBaseName](@ref pm_sysPath::getIndexBaseName)<br>
    !>
    !>  \example{getDirName}
    !>  \include{lineno} example/pm_sysPath/getDirName/main.F90
    !>  \compilef{getDirName}
    !>  \output{getDirName}
    !>  \include{lineno} example/pm_sysPath/getDirName/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{getDirName}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getDirName

    PURE module function getDirNameDef(path, dirsep) result(dirname)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirNameDef
#endif
        character(*, SK)    , intent(in)    :: path
        character(*, SK)    , intent(in)    :: dirsep
        character(:, SK)    , allocatable   :: dirname
    end function

    PURE module function getDirNamePM(path, dirsep, style) result(dirname)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirNamePM
#endif
        character(*, SK)    , intent(in)    :: path
        character(*, SK)    , intent(in)    :: dirsep
        type(verbatim_type) , intent(in)    :: style
        character(:, SK)    , allocatable   :: dirname
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **file extension** part of the input `path` (including the leading dot in the file extension).<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the path to split.<br>
    !>  \param[in]  dirsep  :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the runtime directory separator(s).<br>
    !>                          The directory separator can be obtained from either the `dirsep` component of an object of type
    !>                          [shell_type](@ref pm_sysShell::shell_type) or [getDirSep](@ref pm_sysPath::getDirSep).<br>
    !>                          The directory separator is `\` or `/` in Windows-based terminals (e.g., CMD or PowerShell) and `/` in POSIX-compliant shells.<br>
    !>                          When in doubt (for example, in Windows terminals), `dirsep` can be set to multiple characters, for example `dirsep = "/\"`.<br>
    !>                          In such a case, the input `path` will be scanned for the presence of any of the individual characters in `dirsep`.<br>
    !>
    !>  \return
    !>  `extname`           :   The output `allocatable` scalar `character` of default kind \SK containing the file extension segment of the input `path`.<br>
    !>
    !>  \interface{getExtName}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getExtName
    !>      character(:, SK), allocatable :: extname
    !>
    !>      extname = getExtName(path, dirsep)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < len(dirsep)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  See [getIndexExtName](@ref pm_sysPath::getIndexExtName) for a faster
    !>  method of getting the file extension segment of a path that avoids memory allocation.<br>
    !>
    !>  \see
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>  [getExtName](@ref pm_sysPath::getExtName)<br>
    !>  [getFileName](@ref pm_sysPath::getFileName)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getIndexDirName](@ref pm_sysPath::getIndexDirName)<br>
    !>  [getIndexExtName](@ref pm_sysPath::getIndexExtName)<br>
    !>  [getIndexBaseName](@ref pm_sysPath::getIndexBaseName)<br>
    !>
    !>  \example{getExtName}
    !>  \include{lineno} example/pm_sysPath/getExtName/main.F90
    !>  \compilef{getExtName}
    !>  \output{getExtName}
    !>  \include{lineno} example/pm_sysPath/getExtName/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{getExtName}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getExtName
    PURE module function getExtName(path, dirsep) result(extname)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtName
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    :: path
        character(*,SKG), intent(in)    :: dirsep
        character(:,SKG), allocatable   :: extname
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **basename** part of the input `path`
    !>  (i.e., all characters that appear after the last directory separator in the `path`).<br>
    !>
    !>  \details
    !>  See also the details of the documentation of the generic interface [getIndexDirName](@ref pm_sysPath::getIndexDirName) and [getDirName](@ref pm_sysPath::getDirName).<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of default kind \SK of
    !>                          arbitrary length type parameter containing the path to split.<br
    !>  \param[in]  dirsep  :   The input scalar `character` of default kind \SK of
    !>                          arbitrary length type parameter containing the runtime directory separator(s).<br>
    !>                          The directory separator can be obtained from either the `dirsep` component of an object of type
    !>                          [shell_type](@ref pm_sysShell::shell_type) or [getDirSep](@ref pm_sysPath::getDirSep).<br>
    !>                          The directory separator is `\` or `/` in Windows-based terminals (e.g., CMD or PowerShell) and `/` in POSIX-compliant shells.<br>
    !>                          When in doubt (for example, in Windows terminals), `dirsep` can be set to multiple characters, for example `dirsep = "/\"`.<br>
    !>                          In such a case, the input `path` will be scanned for the presence of any of the individual characters in `dirsep`.<br>
    !>  \param[in]  style   :   The input scalar that can be,<br>
    !>                          <ol>
    !>                              <li>    the constant [verbatim](@ref pm_sysPath::verbatim) or an object of type [verbatim_type](@ref pm_sysPath::verbatim_type)
    !>                                      implying the use of the ParaMonte style described above in extracting the output `basename`.<br>
    !>                          </ol>
    !>                          (**optional**. If missing, the default behavior corresponds to that of the `basename` command on POSIX systems. See examples below.)
    !>
    !>  \return
    !>  `basename`          :   The output `allocatable` scalar `character` of default kind \SK containing the basename segment of the input `path`.<br>
    !>
    !>  \interface{getBaseName}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getBaseName, verbatim
    !>      character(:, SK), allocatable :: basename
    !>
    !>      basename = getBaseName(path, dirsep)
    !>      basename = getBaseName(path, dirsep, style)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < len(dirsep)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  See [getIndexBaseName](@ref pm_sysPath::getIndexBaseName) for a faster
    !>  method of getting the `basename` segment of a path that avoids memory allocation.<br>
    !>
    !>  \see
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>  [getExtName](@ref pm_sysPath::getExtName)<br>
    !>  [getFileName](@ref pm_sysPath::getFileName)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getIndexDirName](@ref pm_sysPath::getIndexDirName)<br>
    !>  [getIndexExtName](@ref pm_sysPath::getIndexExtName)<br>
    !>  [getIndexBaseName](@ref pm_sysPath::getIndexBaseName)<br>
    !>
    !>  \example{getBaseName}
    !>  \include{lineno} example/pm_sysPath/getBaseName/main.F90
    !>  \compilef{getBaseName}
    !>  \output{getBaseName}
    !>  \include{lineno} example/pm_sysPath/getBaseName/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{getBaseName}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getBaseName

    PURE module function getBaseNameDef(path, dirsep) result(basename)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBaseNameDef
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)    :: path
        character(*,SKG)    , intent(in)    :: dirsep
        character(:,SKG)    , allocatable   :: basename
    end function

    PURE module function getBaseNamePM(path, dirsep, style) result(basename)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBaseNamePM
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG)    , intent(in)    :: path
        character(*,SKG)    , intent(in)    :: dirsep
        type(verbatim_type) , intent(in)    :: style
        character(:,SKG)    , allocatable   :: basename
    end function

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the filename part of the input `path`.<br>
    !>
    !>  \details
    !>  The **filename** is defined as the segment of the path that starts after the last directory separator in the `path`
    !>  and **before** the start of the file extension of the path (if any file extension exists).<br>
    !>  For example, the filename segment of `path = "./paramonte.txt"` is `paramonte`.<br>
    !>
    !>  \param[in]  path    :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the path to split.<br>
    !>  \param[in]  dirsep  :   The input scalar `character` of default kind \SK of arbitrary length type parameter containing the runtime directory separator(s).<br>
    !>                          The directory separator can be obtained from either the `dirsep` component of an object of type
    !>                          [shell_type](@ref pm_sysShell::shell_type) or [getDirSep](@ref pm_sysPath::getDirSep).<br>
    !>                          The directory separator is `\` or `/` in Windows-based terminals (e.g., CMD or PowerShell) and `/` in POSIX-compliant shells.<br>
    !>                          When in doubt (for example, in Windows terminals), `dirsep` can be set to multiple characters, for example `dirsep = "/\"`.<br>
    !>                          In such a case, the input `path` will be scanned for the presence of any of the individual characters in `dirsep`.<br>
    !>
    !>  \return
    !>  `filename`          :   The output `allocatable` scalar `character` of default kind \SK containing the filename segment of the input `path`.<br>
    !>
    !>  \interface{getFileName}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getFileName
    !>      character(:, SK), allocatable :: filename
    !>
    !>      filename = getFileName(path, dirsep)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < len(dirsep)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getDirName](@ref pm_sysPath::getDirName)<br>
    !>  [getExtName](@ref pm_sysPath::getExtName)<br>
    !>  [getFileName](@ref pm_sysPath::getFileName)<br>
    !>  [getBaseName](@ref pm_sysPath::getBaseName)<br>
    !>  [getIndexDirName](@ref pm_sysPath::getIndexDirName)<br>
    !>  [getIndexExtName](@ref pm_sysPath::getIndexExtName)<br>
    !>  [getIndexBaseName](@ref pm_sysPath::getIndexBaseName)<br>
    !>
    !>  \example{getFileName}
    !>  \include{lineno} example/pm_sysPath/getFileName/main.F90
    !>  \compilef{getFileName}
    !>  \output{getFileName}
    !>  \include{lineno} example/pm_sysPath/getFileName/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{getFileName}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getFileName
    PURE module function getFileName(path, dirsep) result(filename)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFileName
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    :: path
        character(*,SKG), intent(in)    :: dirsep
        character(:,SKG), allocatable   :: filename
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a list of directory paths matching a requested system application.<br>
    !>
    !>  \details
    !>  This is a convenience wrapper for the lower-level generic interface [setPathMatch](@ref pm_sysPath::setPathMatch).<br>
    !>
    !>  \param[in]      key     :   The input scalar `character` of default kind \SK of arbitrary length type parameter
    !>                              containing a list of search keywords separated by the specified input `sep` character.<br>
    !>                              If the input `key` is empty, the entire list of paths will be returned in the output `list`.<br>
    !>                              Similarly, an empty search keyword item in the input `key` evaluates to match all paths.<br>
    !>                              (**optional**. If missing, the entire list of paths will be returned in the output `list`.)
    !>  \param[in]      inc     :   The input scalar `character` of the same type and kind as the input `key` of arbitrary length type parameter
    !>                              containing the (partial) name of a target file, subfolder, or application that should be present in the output directory paths.<br>
    !>                              If the input `inc` is missing or empty, the entire list of paths matching `key` will be returned in the output `list`.<br>
    !>                              (**optional**, default = `""`)
    !>  \param[in]      sep     :   The input scalar non-blank `character` of the type and kind as `key` of length type parameter `1`,
    !>                              containing the character to be used for separating search keywords provided in key.<br>
    !>                              The recommended value for `sep` is the returned value from [getPathSep](@ref pm_sysPath::getPathSep).<br>
    !>                              (**optional**, default = [getPathSep](@ref pm_sysPath::getPathSep))
    !>  \param[out]     failed  :   The output scalar `logical` of default kind \LK that is `.true.` <b>if and only if</b>
    !>                              the procedure accomplishes the task successfully, otherwise, it is `.false.`.<br>
    !>                              (**optional**. If missing and the procedure fails, the program will halt by calling `error stop`.)
    !>  \param[inout]   errmsg  :   The input/output scalar `character` of default kind \SK of arbitrary length type parameter.<br>
    !>                              If present and an error occurs, it is assigned an explanatory message describing the nature of the error that has occurred.<br>
    !>                              A length of [LEN_IOMSG](@ref pm_io::LEN_IOMSG) characters is likely sufficient to capture most error messages in full.<br>
    !>                              (**optional**. Its presence is relevant only if `failed` is also present.)
    !>
    !>  \return
    !>  `list`                  :   The output `allocatable` vector of type [css_type](@ref pm_container::css_type)
    !>                              each element of which contains a system application path found to match the search keywords in the
    !>                              input `key` and to contain an application matching the specified `inc`, **without case-sensitivity**.<br>
    !>                              The output vector will have a size of `0` if no matching application path is found.<br>
    !>
    !>  \interface{getPathMatch}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: getPathMatch
    !>      use pm_container, only: css_type
    !>      type(css_type), allocatable :: list(:)
    !>
    !>      list = getPathMatch(key = key, inc = inc, sep = sep, failed = failed, errmsg = errmsg)
    !>
    !>  \endcode
    !>
    !>  \impure
    !>
    !>  \note
    !>  Dropping both optional input arguments in this generic interface is equivalent to
    !>  retrieving a list of all application paths found in the environmental `PATH` variable of the processor.<br>
    !>
    !>  \see
    !>  [getPathMatch](@ref pm_sysPath::getPathMatch)<br>
    !>  [setPathMatch](@ref pm_sysPath::setPathMatch)<br>
    !>
    !>  \example{getPathMatch}
    !>  \include{lineno} example/pm_sysPath/getPathMatch/main.F90
    !>  \compilef{getPathMatch}
    !>  \output{getPathMatch}
    !>  \include{lineno} example/pm_sysPath/getPathMatch/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{getPathMatch}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getPathMatch
    module function getPathMatch(key, sep, inc, failed, errmsg) result(list)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPathMatch
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        character(*,SKG), intent(in)    , optional  :: key
        character(*,SKG), intent(in)    , optional  :: inc
        character(1,SKG), intent(in)    , optional  :: sep
        logical(LK)     , intent(out)   , optional  :: failed
        character(*,SKG), intent(inout) , optional  :: errmsg
        type(css_type)  , allocatable               :: list(:)
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a list of directory paths matching a requested system application.<br>
    !>
    !>  \details
    !>  This generic interface is meant to provide an easy unified cross-language way of
    !>  detecting applications, particularly, vendor MPI installations on the system.<br>
    !>  The original approach of implementing this functionality separately
    !>  in different programming language environments proved messy.<br>
    !>  This generic interface along with [getPathMatch](@ref pm_sysPath::getPathMatch)
    !>  attempt to reduce multi-language code duplication to only one (Fortran) language.<br>
    !>
    !>  **Search Strategy**<br>
    !>  The search strategy relies on retrieving the output of the command `system_profiler SPApplicationsDataType` in Darwin OS.<br>
    !>  The search strategy relies on retrieving the contents of the processor environmental `PATH` variable in Windows and Linux OS.<br>
    !>  Each identified path will be checked to ensure it is a real path before inclusion in the result.<br>
    !>
    !>  \param[in]      key         :   The input scalar `character` of default kind \SK of arbitrary
    !>                                  non-zero length type parameter containing a list of search keywords
    !>                                  separated by the specified input `sep` character.<br>
    !>                                  An empty `key` (with zero length type parameter) leads to
    !>                                  returning all identified paths in the processor environment.<br>
    !>                                  Similarly, an empty search keyword item in the input `key` evaluates to match all paths.<br>
    !>                                  **The input `key` must be null-terminated when the procedure is called from the C language**.<br>
    !>  \param[in]      inc         :   The input scalar `character` of the same type and kind as the input `key` of arbitrary length type parameter
    !>                                  containing the (partial) name of a target file, subfolder, or application that should be present in the output directory paths.<br>
    !>                                  If the input `inc` is missing or empty, the entire list of paths matching `key` will be returned in the output `list`.<br>
    !>                                  **The input `inc` must be null-terminated when the procedure is called from the C language**.<br>
    !>  \param[in]      sep         :   The input scalar `character` of the type and kind as `key` of length type parameter `1`,
    !>                                  containing the character to be used for separating search keywords provided in key.<br>
    !>                                  The recommended value for `sep` is the returned value from [getPathSep](@ref pm_sysPath::getPathSep).<br>
    !>  \param[out]     list        :   The output scalar `character` of the same type and kind as the input `key`,
    !>                                  containing a list of the identified paths matching the input search keywords **without case-sensitivity**.<br>
    !>                                  The paths, if any more than one are detected, are separated by the specified input `sep` character.<br>
    !>                                  The interpretation of the output contents of `list` depends on the value of the input/output argument `lenList`.<br>
    !>  \param[inout]   lenList     :   The input/output scalar `integer` of default kind \IK.<br>
    !>                                  On input, it must contain the positive-valued length of the slice of `list` buffer to be used for storing the results.<br>
    !>                                  On output, depending on its value, `lenList` can signal different scenarios:<br>
    !>                                  <ol>
    !>                                      <li>    A positive output value larger than its input value implies the slice `list` contains the directory path only partially.<br>
    !>                                              In such a case, user may want to call the routine with a larger input buffer `list` length specified via `lenList` to store the path in full.<br>
    !>                                      <li>    A positive output value smaller than its input value implies the slice `list(1:lenList)` contains the directory path in full.<br>
    !>                                      <li>    A negative output value implies that and error occurred and/or no installed application path matching the specified keywords could be found.<br>
    !>                                              In such a case, the slice `list(1 : abs(lenList))` contains a description of the error or other informative messages.<br>
    !>                                              The interpretation of the value of `abs(lenList)` is subject to the same rules for the positive value of `lenList` laid out above.<br>
    !>                                      <li>    A zero output value implies that no installed application path matching the specified keywords could be found on the system.<br>
    !>                                  </ol>
    !>                                  A value of [MAX_LEN_FILE_PATH](@ref pm_sysPath::MAX_LEN_FILE_PATH)
    !>                                  for `lenList` should be more than sufficient with focused keyword search.<br>
    !>                                  If the input `key` contents are expected to match many paths, a longer buffer `list` may be needed.<br>
    !>                                  In such cases, the length of the environmental `PATH` variable can be used a rule of thumb for input value of `lenList`.<br>
    !>
    !>  \interface{setPathMatch}
    !>  \code{.F90}
    !>
    !>      use pm_sysPath, only: setPathMatch
    !>
    !>      call setPathMatch(key, inc, sep, list, lenList)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < lenList` must hold for the corresponding input arguments.<br>
    !>  The input argument `key` must be null-terminated when the procedure is called from the C language.<br>
    !>
    !>  \impure
    !>
    !>  \note
    !>  This generic interface is primarily intended for access from other programming language environments.<br>
    !>  For a more flexible within-Fortran interface see the generic function interface [getPathMatch](@ref pm_sysPath::getPathMatch).<br>
    !>
    !>  \see
    !>  [getPathMatch](@ref pm_sysPath::getPathMatch)<br>
    !>  [setPathMatch](@ref pm_sysPath::setPathMatch)<br>
    !>
    !>  \example{setPathMatch}
    !>  \include{lineno} example/pm_sysPath/setPathMatch/main.F90
    !>  \compilef{setPathMatch}
    !>  \output{setPathMatch}
    !>  \include{lineno} example/pm_sysPath/setPathMatch/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sysPath](@ref test_pm_sysPath)
    !>
    !>  \todo
    !>  \phigh
    !>  The current Fortran standard 202x does not allow passing characters of non-default kind to the intrinsic Fortran statements and procedures.<br>
    !>  As such, the implementation of this procedure for non-default `character` kinds leads to compile-time kind mismatch errors.<br>
    !>  This procedure should be converted back to a generic interface in the future when non-default character kinds are also fully supported by the intrinsic functions.<br>
    !>
    !>  \final{setPathMatch}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 7:03 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setPathMatch
    module subroutine setPathMatch(key, inc, sep, list, lenList) !BIND2("setPathMatch")
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPathMatch
#endif
        use pm_kind, only: SKG => SK
        character(*,SKG), intent(in)    :: key
        character(*,SKG), intent(in)    :: inc
        character(1,SKG), intent(in)    :: sep
        character(*,SKG), intent(out)   :: list
        integer(IK)     , intent(inout) :: lenList
    end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sysPath ! LCOV_EXCL_LINE

#undef MAX_PATH
#undef MAX_NAME
#undef SEP_PATH
#undef BIND2
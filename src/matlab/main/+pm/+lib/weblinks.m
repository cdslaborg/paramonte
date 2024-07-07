%>  \brief
%>  Return a structure containing tree of weblinks for the
%>  ParaMonte MATLAB library source file and documentation website.
%>
%>  \return
%>  ``tree``    :   The output MATLAB ``struct`` containing the ParaMonte website information.
%>
%>  \interface{weblinks}
%>  \code{.m}
%>
%>      tree = pm.lib.weblinks();
%>
%>  \endcode
%>
%>  \example{weblinks}
%>  \include{lineno} example/lib/weblinks/main.m
%>  \output{weblinks}
%>  \include{lineno} example/lib/weblinks/main.out.m
%>
%>  \final{weblinks}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 7:58 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function tree = weblinks()

    persistent stree;

    if ~isempty(stree)
        tree = stree;
        return
    end

    stree = struct();
    stree.home = struct();
    stree.home.url = "https://www.cdslab.org/paramonte";
    stree.home.install = struct();
    stree_home_install_url = stree.home.url + "/notes/installation";

    %%%% installation Linux

    stree.home.overview = struct();
    stree_home_overview_url = stree.home.url + "/notes/overview";
    stree.home.overview.preface = struct();
    stree.home.overview.changes = struct();
    stree.home.overview.preface.url = stree_home_overview_url + "/preface";
    stree.home.overview.changes.fortran = struct();
    stree.home.overview.changes.python = struct();
    stree.home.overview.changes.matlab = struct();
    stree.home.overview.changes.fortran.url = stree_home_overview_url + "/paramonte-kernel-release-notes";
    stree.home.overview.changes.python.url = stree_home_overview_url + "/paramonte-python-release-notes";
    stree.home.overview.changes.matlab.url = stree_home_overview_url + "/paramonte-matlab-release-notes";

    %%%% installation Linux

    stree.home.install.linux = struct();
    stree.home.install.linux.url = stree_home_install_url + "/linux";

    %%%% installation Windows

    stree.home.install.windows = struct();
    stree.home.install.windows.url = stree_home_install_url + "/windows";

    %%%% installation MATLAB

    stree.home.install.matlab = struct();
    stree.home.install.matlab.url = stree_home_install_url + "/matlab";

    %%%% installation Python

    stree.home.install.python = struct();
    stree.home.install.python.url = stree_home_install_url + "/python";

    %%%% installation macOS

    stree.home.install.macos = struct();
    stree.home.install.macos.url = stree_home_install_url + "/macos";
    stree.home.install.macos.prereqs = struct();
    stree.home.install.macos.prereqs.url = stree.home.install.macos.url + "/#the-compile-time-and-runtime-prerequisites";
    stree.home.install.macos.prereqs.cmd = struct();
    stree.home.install.macos.prereqs.cmd.url = stree.home.install.macos.url + "/#prereqs-install";

    %%%% MATLAB examples

    stree.home.examples = struct();
    stree_home_examples_url = stree.home.url + "/notes/examples";
    stree.home.examples.matlab = struct();
    stree.home.examples.matlab.jupyter = struct();
    stree.home.examples.matlab.postprocess = struct();
    stree.home.examples.matlab.jupyter.url = stree_home_examples_url + "/matlab/jupyter";
    stree.home.examples.matlab.postprocess.url = stree_home_examples_url + "/matlab/postprocess";

    %%%% Python examples

    stree.home.examples = struct();
    stree_home_examples_url = stree.home.url + "/notes/examples";
    stree.home.examples.python = struct();
    stree.home.examples.python.jupyter = struct();
    stree.home.examples.python.postprocess = struct();
    stree.home.examples.python.jupyter.url = stree_home_examples_url + "/python/jupyter";
    stree.home.examples.python.postprocess.url = stree_home_examples_url + "/python/postprocess";

    %%%% Python API

    stree.home.api = struct();
    stree_home_api_url = stree.home.url + "/notes/api";
    stree.home.api.python = struct();
    stree.home.api.python.url = stree_home_api_url + "/python/autoapi/paramonte";

    %%%% ParaDRAM

    stree.home.usage = struct();
    stree_home_usage_url = stree.home.url + "/notes/usage";
    stree.home.usage.paradram = struct();
    stree_home_usage_paradram_url = stree_home_usage_url + "/paradram";
    stree.home.usage.paradram.quickstart = struct();
    stree.home.usage.paradram.quickstart.url = stree_home_usage_paradram_url + "/interface";
    stree.home.usage.paradram.input = struct();
    stree.home.usage.paradram.input.url = stree_home_usage_paradram_url + "/input";
    stree.home.usage.paradram.specifications = struct();
    stree.home.usage.paradram.specifications.url = stree_home_usage_paradram_url + "/specifications";
    stree.home.usage.paradram.restart = struct();
    stree.home.usage.paradram.restart.url = stree_home_usage_paradram_url + "/restart";
    stree.home.usage.paradram.output = struct();
    stree.home.usage.paradram.output.url = stree_home_usage_paradram_url + "/output";

    %%%% GitHub issues

    stree.github = struct();
    stree.github.url = "https://github.com/cdslaborg/paramonte";
    stree.github.issues = struct();
    stree.github.issues.url = "https://github.com/cdslaborg/paramonte/issues";
    stree.github.release = struct();
    stree.github.release.url = stree.github.url + "/releases";
    stree.github.release.latest = struct();
    stree.github.release.latest.url = stree.github.release.url + "/latest";
    stree.github.archive = struct();
    stree_github_archive_url = stree.github.url + "/archive";
    stree.github.archive.main = struct();
    stree.github.archive.main.zip = struct();
    stree.github.archive.main.tar = struct();
    stree.github.archive.main.zip.url = stree_github_archive_url + "/main.zip";
    stree.github.archive.main.tar.url = stree_github_archive_url + "/main.tar.gz";

    %%%% GitHub examples

    stree.github.examples = struct();
    stree.github.examples.url = "https://github.com/cdslaborg/paramontex";

    %%%% Intel MPI

    stree.intel = struct();
    stree.intel.mpi = struct();
    stree.intel.mpi.home = struct();
    stree.intel.mpi.home.url = "https://software.intel.com/en-us/mpi-library";

    %%%% Intel MPI Windows

    stree.intel.mpi.windows = struct();
    stree.intel.mpi.windows.url = "https://software.intel.com/en-us/get-started-with-mpi-for-windows";

    %%%% OpenMPI

    stree.openmpi = struct();
    stree.openmpi.home = struct();
    stree.openmpi.home.url = "https://www.open-mpi.org/";

    tree = stree;

end
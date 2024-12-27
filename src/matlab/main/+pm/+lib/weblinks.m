%>  \brief
%>  Return a structure containing tree of weblinks for the
%>  ParaMonte MATLAB library source file and documentation website.<br>
%>
%>  \return
%>  ``tree``    :   The output MATLAB ``struct`` containing the ParaMonte website information.<br>
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

    %%%% docs

    stree.docs = struct();
    stree.docs.url = "https://www.cdslab.org/paramonte";

    %%%% docs generic

    stree.docs.generic = struct();
    stree.docs.generic.url = stree.docs.url + "/generic/" + pm.lib.version("generic", "major");

    %%%% docs lang

    for lang = ["matlab"]; %["c", "cpp", "fortran", "matlab", "python"];
        stree.docs.(lang) = struct();
        stree.docs.(lang).url = stree.docs.url + "/" + lang + "/" + pm.lib.version(lang, "major");
    end

    %%%% docs generic overview

    stree.docs.generic.overview = struct();
    stree.docs.generic.overview.url = stree.docs.generic.url + "/overview";

    stree.docs.generic.overview.preface = struct();
    stree.docs.generic.overview.preface.url = stree.docs.generic.overview.url + "/preface";

    stree.docs.generic.overview.changes = struct();
    stree.docs.generic.overview.changes.url = stree.docs.generic.overview.url + "/CHANGES.md";

    %%%% docs generic installation

    stree.docs.generic.installation = struct();
    stree.docs.generic.installation.url = stree.docs.generic.url + "/installation";

    %%%% docs generic installation Linux

    stree.docs.generic.installation.linux = struct();
    stree.docs.generic.installation.linux.url = stree.docs.generic.installation.url + "/linux";

    %%%% docs generic installation Windows

    stree.docs.generic.installation.windows = struct();
    stree.docs.generic.installation.windows.url = stree.docs.generic.installation.url + "/windows";

    %%%% docs generic installation MATLAB

    stree.docs.generic.installation.matlab = struct();
    stree.docs.generic.installation.matlab.url = stree.docs.generic.installation.url + "/matlab";

    %%%% docs generic installation Python

    stree.docs.generic.installation.python = struct();
    stree.docs.generic.installation.python.url = stree.docs.generic.installation.url + "/python";

    %%%% docs generic installation macOS

    stree.docs.generic.installation.macos = struct();
    stree.docs.generic.installation.macos.url = stree.docs.generic.installation.url + "/macos";
    stree.docs.generic.installation.macos.prereqs = struct();
    stree.docs.generic.installation.macos.prereqs.url = stree.docs.generic.installation.macos.url + "/#the-compile-time-and-runtime-prerequisites";
    stree.docs.generic.installation.macos.prereqs.cmd = struct();
    stree.docs.generic.installation.macos.prereqs.cmd.url = stree.docs.generic.installation.macos.url + "/#prereqs-install";

    %%%% docs generic MATLAB examples

    stree.docs.generic.examples = struct();
    stree.docs.generic.examples.url = stree.docs.generic.url + "/examples";
    stree.docs.generic.examples.matlab = struct();
    stree.docs.generic.examples.matlab.jupyter = struct();
    stree.docs.generic.examples.matlab.postprocess = struct();
    stree.docs.generic.examples.matlab.jupyter.url = stree.docs.generic.examples.url + "/matlab/jupyter";
    stree.docs.generic.examples.matlab.postprocess.url = stree.docs.generic.examples.url + "/matlab/postprocess";

    %%%% docs generic Python examples

    stree.docs.generic.examples = struct();
    stree.docs.generic.examples.url = stree.docs.generic.url + "/examples";
    stree.docs.generic.examples.python = struct();
    stree.docs.generic.examples.python.jupyter = struct();
    stree.docs.generic.examples.python.postprocess = struct();
    stree.docs.generic.examples.python.jupyter.url = stree.docs.generic.examples.url + "/python/jupyter";
    stree.docs.generic.examples.python.postprocess.url = stree.docs.generic.examples.url + "/python/postprocess";

    %%%% docs generic Python API

    stree.docs.generic.api = struct();
    stree.docs.generic.api.url = stree.docs.generic.url + "/api";
    stree.docs.generic.api.python = struct();
    stree.docs.generic.api.python.url = stree.docs.generic.api.url + "/python/autoapi/paramonte";

    %%%% docs generic usage

    stree.docs.generic.usage = struct();
    stree.docs.generic.usage.url = stree.docs.generic.url + "/usage";

    %%%% docs generic usage sampling

    stree.docs.generic.usage.sampling = struct();
    stree.docs.generic.usage.sampling.url = stree.docs.generic.usage.url + "/sampling";

    %%%% docs generic usage sampling ParaDRAM

    stree.docs.generic.usage.sampling.paradram = struct();
    stree.docs.generic.usage.sampling.paradram.url = stree.docs.generic.usage.sampling.url + "/paradram";

    %%%% docs generic usage sampling ParaDRAM pages

    stree.docs.generic.usage.sampling.paradram.quickstart = struct();
    stree.docs.generic.usage.sampling.paradram.quickstart.url = stree.docs.generic.usage.sampling.paradram.url + "/interface";

    stree.docs.generic.usage.sampling.paradram.input = struct();
    stree.docs.generic.usage.sampling.paradram.input.url = stree.docs.generic.usage.sampling.paradram.url + "/input";

    stree.docs.generic.usage.sampling.paradram.specifications = struct();
    stree.docs.generic.usage.sampling.paradram.specifications.url = stree.docs.generic.usage.sampling.paradram.url + "/specifications";

    stree.docs.generic.usage.sampling.paradram.restart = struct();
    stree.docs.generic.usage.sampling.paradram.restart.url = stree.docs.generic.usage.sampling.paradram.url + "/restart";

    stree.docs.generic.usage.sampling.paradram.output = struct();
    stree.docs.generic.usage.sampling.paradram.output.url = stree.docs.generic.usage.sampling.paradram.url + "/output";

    %%%% GitHub issues

    stree.github = struct();
    stree.github.url = "https://github.com/cdslaborg/paramonte";

    stree.github.issues = struct();
    stree.github.issues.url = "https://github.com/cdslaborg/paramonte/issues";

    %%%% GitHub releases

    stree.github.releases = struct();
    stree.github.releases.url = stree.github.url + "/releases";

    stree.github.releases.latest = struct();
    stree.github.releases.latest.url = stree.github.releases.url + "/latest";

    stree.github.releases.tag = struct();
    stree.github.releases.tag.url = stree.github.releases.url + "/tag";

    stree.github.releases.tag.auxil = struct();
    stree.github.releases.tag.auxil.url = stree.github.releases.tag.url + "/auxil";

    stree.github.releases.download = struct();
    stree.github.releases.download.url = stree.github.releases.url + "/download";

    stree.github.releases.download.auxil = struct();
    stree.github.releases.download.auxil.url = stree.github.releases.download.url + "/auxil";

    %%%% GitHub archive

    stree.github.archive = struct();
    stree.github.archive.url = stree.github.url + "/archive";

    stree.github.archive.main = struct();

    stree.github.archive.main.zip = struct();
    stree.github.archive.main.zip.url = stree.github.archive.url + "/main.zip";

    stree.github.archive.main.tar = struct();
    stree.github.archive.main.tar.url = stree.github.archive.url + "/main.tar.gz";

    %%%% GitHub examples

    stree.github.examples = struct();
    stree.github.examples.url = "https://github.com/cdslaborg/paramontex";

    %%%% external

    stree.external = struct();

    %%%% external Intel

    stree.external.intel = struct();
    stree.external.intel.url = "https://software.intel.com/en-us";

    stree.external.intel.mpi = struct();
    stree.external.intel.mpi.url = stree.external.intel.url + "/mpi-library";

    %%%% external OpenMPI

    stree.external.openmpi = struct();
    stree.external.openmpi.url = "https://www.open-mpi.org/";

    %%%%
    %%%% Copy the tree.
    %%%%

    tree = stree;

end
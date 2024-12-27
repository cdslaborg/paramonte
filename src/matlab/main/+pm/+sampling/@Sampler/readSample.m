%>  \brief
%>  Return a list of objects of superclass [pm.sampling.FileContentsSample](@ref FileContentsSample)
%>  containing the content(s) of the ParaMonte simulation output sample
%>  file(s) whose path(s) match the specified input ``pattern`` or the
%>  simulation specification ``sampler.spec.outputFileName``.<br>
%>
%>  \warning
%>  This method is to be only used for post-processing of the
%>  output sample file(s) of an already finished simulation.<br>
%>  Although possible, this method is **not** meant to be called
%>  by all processes in MPI-parallel simulations.<br>
%>
%>  \param[in]  self    :   The input parent object of class [pm.sampling.Sampler](@ref Sampler)
%>                          which is **implicitly** passed to this dynamic method (not by the user).<br>
%>  \param[in]  pattern :   The input scalar MATLAB string containing the pattern matching
%>                          the desired sample file(s) whose contents is to be read.<br>
%>                          The specified ``pattern`` only needs to partially identify
%>                          the name of the simulation to which the sample file belongs.<br>
%>                          For example, specifying ``"./mydir/mysim"`` as input will
%>                          lead to a search for file(s) beginning with ``"mysim"`` and
%>                          ending with ``"_sample.txt"`` inside the directory ``"./mydir/"``.<br>
%>                          If there are multiple files matching in the input ``pattern``,
%>                          then all such files will be read and returned as elements of a list.<br>
%>                          If the specified pattern is a valid existing URL, the file will be
%>                          downloaded as a temporary file to the local system, its contents
%>                          shall be parsed and the file will be subsequently removed.<br>
%>                          If the input ``pattern`` is empty, then the method will search
%>                          for any possible candidate files with the appropriate suffix
%>                          in the current working directory.<br>
%>                          (**optional**, default = ``sampler.spec.outputFileName`` or ``"./"``)
%>  \param[in]  sep     :   The input MATLAB string containing the field separator used in the file(s).<br>
%>                          (**optional**, default = ``sampler.spec.outputSeparator`` or automatically inferred.)
%>
%>  \return
%>  ``sampleList``      :   The output MATLAB cell array of objects
%>                          of class [pm.sampling.FileContentsSample](@ref FileContentsSample),
%>                          each of which corresponds to the contents
%>                          of a unique sample file.<br>
%>
%>  \interface{readSample}
%>  \code{.m}
%>
%>      sampler = pm.sampling.Sampler();
%>      sampleList = sampler.readSample();
%>      sampleList = sampler.readSample([]);
%>      sampleList = sampler.readSample([], []);
%>      sampleList = sampler.readSample(pattern);
%>      sampleList = sampler.readSample([], sep);
%>      sampleList = sampler.readSample(pattern, sep);
%>
%>  \endcode
%>
%>  \note
%>  See the documentation of the subclasses of [pm.sampling.Sampler](@ref Sampler)
%>  (e.g., [pm.sampling.Paradram](@ref Paradram)) for example usage in action.<br>
%>
%>  \example{readSample}
%>  \code{.m}
%>
%>      sampler.readSample("./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>      sampler.readSample();
%>
%>      sampler.readSample("./out/test_run_", ",");
%>
%>      sampler.spec.outputSeparator = ",";
%>      sampler.readSample("./out/test_run_");
%>
%>      sampler.spec.outputFileName = "./out/test_run_";
%>      sampler.spec.outputSeparator = ",";
%>      sampler.readSample();
%>
%>  \endcode
%>
%>  \final{readSample}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 12:38 AM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function sampleList = readSample(self, pattern, sep)

    if  nargin < 3
        sep = [];
    end

    if  nargin < 2
        pattern = [];
    end

    sampleList = pm.sampling.readSample(self, pattern, sep);

end

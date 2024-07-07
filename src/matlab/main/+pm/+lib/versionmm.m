%>  \brief
%>  Return a scalar MATLAB string containing
%>  the last version before the specified input ParaMonte MATLAB library version number.<br>
%>
%>  \details
%>  The ParaMonte version number follows the
%>  [semantic versioning](https://semver.org/) style: ``major.minor.patch``.<br>
%>  The current (or the input version) is decremented by its patch version.<br>
%>  If the patch version is zero, it is decremented by its minor version.<br>
%>  If the patch version is zero, it is decremented by its major version.<br>
%>
%>  \note
%>  The function name stands for ``version--``.<br>
%>
%>  \param[in]  current :   The input scalar MATLAB string containing the
%>                          current Semantic version of the library in
%>                          triplet format (e.g., ``"1.2.3"``).<br>
%>                          The current ParaMonte library version,
%>                          if needed, can be obtained from
%>                          [pm.lib.version](@ref version).<br>
%>
%>  \return
%>  ``previous``        :   The output scalar MATLAB string containing the
%>                          decremented ParaMonte library version as described above.<br>
%>                          The output version **may not necessarily exist**
%>                          as an actual ParaMonte library version.<br>
%>
%>  \interface{versionmm}
%>  \code{.m}
%>
%>      previous = pm.lib.versionmm(current)
%>
%>  \endcode
%>
%>  \example{versionmm}
%>  \include{lineno} example/lib/versionmm/main.m
%>  \output{versionmm}
%>  \include{lineno} example/lib/versionmm/main.out.m
%>
%>  \final{versionmm}
%>
%>  \author
%>  \JoshuaOsborne, May 21 2024, 8:13 PM, University of Texas at Arlington<br>
%>  \FatemehBagheri, May 20 2024, 1:25 PM, NASA Goddard Space Flight Center (GSFC), Washington, D.C.<br>
%>  \AmirShahmoradi, May 16 2016, 9:03 AM, Oden Institute for Computational Engineering and Sciences (ICES), UT Austin<br>
function previous = versionmm(current)
    currentVerionTriplet = double(strsplit(current, "."));
    previous = [];
    index = 4;
    while true
        index = index - 1;
        if  index <= 0
            return;
        else
            if  currentVerionTriplet(index) > 0
                previousVerionTriplet = currentVerionTriplet;
                previousVerionTriplet(index) = previousVerionTriplet(index) - 1;
                previous = join(string(previousVerionTriplet), ".");
                return;
            end
        end
    end
end
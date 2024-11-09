%>  \brief
%>  Generate and return the relevant Post-Processing Message (ppm) for the current ParaMonte
%>  sampler to be displayed on the MATLAB command line after the sampling is complete.<br>
%>
%>  \details
%>  This is a ``private`` dynamic method of the [pm.sampling.Paradram](@ref Paradram) sampler class.<br>
%>  This method is not meant to be used or accessed by the end users.<br>
%>
%>  \param[in]  self    :   The input parent object of class [pm.sampling.Sampler](@ref Sampler)
%>                          which is **implicitly** passed to this dynamic method (not by the user).<br>
%>
%>  \note
%>  This is an internal method of the class [pm.sampling.Sampler](@ref Sampler).
%>
%>  \final{getppm}
%>
%>  \author
%>  \AmirShahmoradi, September 1, 2012, 12:00 AM, National Institute for Fusion Studies, The University of Texas at Austin%>
function ppm = getppm(self)
ppm = getppm@pm.sampling.Sampler(self) + newline ...
+ "Use the following object method to read the generated basic output chain file and and unroll the contents as a Markov Chain: " + newline ...
+ newline ...
+ pm.io.tab + self.name + ".readChainMarkov() % Return a list of the unrolled contents of the output chain file(s) as Markov Chains." + newline ...
+ newline ...
+ "Beware that the chain unrolling significantly increases the chain size and can be very slow." + newline ...
+ "It can potentially overflow the computer RAM for high-dimensional target density functions." + newline ...
;
end
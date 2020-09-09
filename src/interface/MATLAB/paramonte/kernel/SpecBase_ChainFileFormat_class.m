%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   MIT License
%%%%
%%%%   ParaMonte: plain powerful parallel Monte Carlo library.
%%%%
%%%%   Copyright (C) 2012-present, The Computational Data Science Lab
%%%%
%%%%   This file is part of the ParaMonte library.
%%%%
%%%%   Permission is hereby granted, free of charge, to any person obtaining a 
%%%%   copy of this software and associated documentation files (the "Software"), 
%%%%   to deal in the Software without restriction, including without limitation 
%%%%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%%%%   and/or sell copies of the Software, and to permit persons to whom the 
%%%%   Software is furnished to do so, subject to the following conditions:
%%%%
%%%%   The above copyright notice and this permission notice shall be 
%%%%   included in all copies or substantial portions of the Software.
%%%%
%%%%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%%%%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%%%%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%%%%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
%%%%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
%%%%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%%%%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%
%%%%   ACKNOWLEDGMENT
%%%%
%%%%   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
%%%%   As per the ParaMonte library license agreement terms, if you use any parts of 
%%%%   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
%%%%   work (education/research/industry/development/...) by citing the ParaMonte 
%%%%   library as described on this page:
%%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef  SpecBase_ChainFileFormat_class < handle

    properties (Constant)
        CLASS_NAME                  = "@SpecBase_ChainFileFormat_class"
        MAX_LEN_CHAIN_FILE_FORMAT   = 63
    end

    properties
        isCompact 	                = []
        isVerbose 	                = []
        isBinary 	                = []
        compact                     = ''
        verbose                     = ''
        binary                      = ''
        def                         = ''
        val                         = ''
        desc                        = ''
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function self = SpecBase_ChainFileFormat_class(methodName)
            self.isCompact  = false;
            self.isVerbose  = false;
            self.isBinary   = false;
            self.compact    = "compact";
            self.verbose    = "verbose";
            self.binary     = Constants.FILE_TYPE.binary;
            self.def        = self.compact;

            self.desc       = "chainFileFormat is a string variable that represents the format of the output chain file(s) of " + methodName                                                ...
                            + " simulation. The string value must be enclosed by either single or double quotation marks when provided as input."                                            ...
                            + "Three values are possible:" + newline + newline                                                                                                              ...
                            + "    chainFileFormat = 'compact'" + newline + newline                                                                                                         ...
                            + "            This is the ASCII (text) file format which is human-readable but does not preserve the full accuracy of the "                                    ...
                            + "output values. It is also a significantly slower mode of chain file generation, compared to the binary file format (see below). "                            ...
                            + "If the compact format is specified, each of the repeating MCMC states will be condensed into a single entry (row) in "                                       ...
                            + "the output MCMC chain file. Each entry will be then assigned a sample-weight that is equal to the number of repetitions of "                                 ...
                            + "that state in the MCMC chain. Thus, each row in the output chain file will represent a unique sample from the objective function. "                          ...
                            + "This will lead to a significantly smaller ASCII chain file and faster output size compared to the verbose chain file format (see below)." + newline + newline...
                            + "    chainFileFormat = 'verbose'" + newline + newline                                                                                                         ...
                            + "            This is the ASCII (text) file format which is human-readable but does not preserve the full accuracy of "                                        ...
                            + "the output values. It is also a significantly slower mode of chain file generation, "                                                                        ...
                            + "compared to both compact and binary chain file formats (see above and below). "                                                                              ...
                            + "If the verbose format is specified, all MCMC states will have equal sample-weights of 1 in the output chain file. "                                          ...
                            + "The verbose format can lead to much larger chain file sizes than the compact and binary file formats. "                                                      ...
                            + "This is especially true if the target objective function has a very high-dimensional state space." + newline + newline                                       ...
                            + "    chainFileFormat = '" + self.binary + "'" + newline + newline                                                                                             ...
                            + "            This is the binary file format which is not human-readable, but preserves the exact values in the output "                                       ...
                            + "MCMC chain file. It is also often the fastest mode of chain file generation. If the binary file format is chosen, the chain "                                ...
                            + "will be automatically output in the compact format (but as binary) to ensure the production of the smallest-possible output chain file. "                    ...
                            + "Binary chain files will have the " + Constants.FILE_EXT.binary + " file extensions. Use the binary format if you need full accuracy representation "         ...
                            + "of the output values while having the smallest-size output chain file in the shortest time possible." + newline + newline                                    ...
                            + "The default value is chainFileFormat = '" + self.def + "' as it provides a reasonable trade-off between "                                                    ...
                            + "speed and output file size while generating human-readable chain file contents. Note that the input values are case-insensitive."                            ...
                            ;
        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function set(self, chainFileFormat)

            if isempty(chainFileFormat)
                self.val = strtrim(self.def);
            else
                self.val = strtrim(chainFileFormat);
            end

            if lower(self.val) == lower(self.compact),  self.isCompact  = true; end
            if lower(self.val) == lower(self.verbose),  self.isVerbose  = true; end
            if lower(self.val) == lower(self.binary),   self.isBinary   = true; end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

        function Err = checkForSanity(self, Err, methodName)

            FUNCTION_NAME = "@checkForSanity()";

            if ~(self.isCompact || self.isVerbose || self.isBinary)
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                       ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                        ...
                                + "The input requested chain file format ('" + self.val + "') represented by the variable "     ...
                                + "chainFileFormat cannot be set to anything other than '" + self.compact + "' or '"            ...
                                + self.verbose + "' or '" + self.binary + "'. If you don't know an appropriate value for "      ...
                                + "chainFileFormat, drop it from the input list. " + methodName                                 ...
                                + " will automatically assign an appropriate value to it." + newline + newline                  ...
                                ;
            end

        end

    %*******************************************************************************************************************************
    %*******************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end
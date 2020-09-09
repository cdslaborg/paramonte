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

classdef SpecDRAM_ScaleFactor_class < handle

    properties (Constant)
        CLASS_NAME                  = "@SpecDRAM_ScaleFactor_class"
        MAX_LEN_STRING_SCALE_FACTOR = 127
    end

    properties
        val                         = []
        defVal                      = []
        str                         = []
        defStr                      = []
        desc                        = []
    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    methods (Access = public)

%***********************************************************************************************************************************
%***********************************************************************************************************************************

        function self = SpecDRAM_ScaleFactor_class(nd, methodName)
            self.defStr = "gelman";
            self.defVal = 2.38/sqrt(nd);  % Gelman, Roberts, Gilks (1996): Efficient Metropolis Jumping Rules
            self.desc   = "scaleFactor is a real-valued positive number (which must be given as string), by which the "                                             ...
                        + "covariance matrix of the proposal distribution of " + methodName + " sampler is scaled. Specifically, "                                  ...
                        + "the proposal distribution will be scaled in every direction by the value of scaleFactor. "                                               ...
                        + "It can also be given in units of the string keyword 'gelman' (which is case-INsensitive) after the paper:" + newline + newline           ...
                        + Decoration_class.TAB + "Gelman, Roberts, and Gilks (1996): 'Efficient Metropolis Jumping Rules'." + newline + newline                     ...
                        + "The paper finds that the optimal scaling factor for the covariance matrix of a Multivariate Gaussian proposal "                          ...
                        + "distribution for the Metropolis-Hastings Markov Chain Monte Carlo sampling of a target Multivariate Normal Distribution "                ...
                        + "of dimension ndim is given by:" + newline + newline                                                                                      ...
                        + "    scaleFactor = 2.38/sqrt(ndim)  ,  in the limit of ndim->Infinity." + newline + newline                                               ...
                        + "Multiples of the gelman scale factors are also acceptable as input and can be specified like the following examples:" + newline + newline...
                        + "    scaleFactor = '1'" + newline + newline                                                                                               ...
                        + "            multiplies the ndim-dimensional proposal covariance matrix by 1, essentially no change occurs to "                           ...
                        + "the covariance matrix." + newline + newline                                                                                              ...
                        + "    scaleFactor = ""1""" + newline + newline                                                                                             ...
                        + "            same as the previous example. The double-quotation marks act the same way as single-quotation marks." + newline + newline    ...
                        + "    scaleFactor = '2.5'" + newline + newline                                                                                             ...
                        + "            multiplies the ndim-dimensional proposal covariance matrix by 2.5." + newline + newline                                      ...
                        + "    scaleFactor = '2.5*Gelman'" + newline + newline                                                                                      ...
                        + "            multiplies the ndim-dimensional proposal covariance matrix by 2.5 * 2.38/sqrt(ndim)." + newline + newline                    ...
                        + "    scaleFactor = ""2.5 * gelman""" + newline + newline                                                                                  ...
                        + "            same as the previous example, but with double-quotation marks. space characters are ignored." + newline + newline            ...
                        + "    scaleFactor = ""2.5 * gelman*gelman*2""" + newline + newline                                                                         ...
                        + "            equivalent to gelmanFactor-squared multiplied by 5." + newline + newline                                                     ...
                        + "Note, however, that the result of Gelman et al. paper applies only to multivariate normal proposal distributions, in "                   ...
                        + "the limit of infinite dimensions. Therefore, care must be taken when using Gelman's scaling factor with non-Gaussian "                   ...
                        + "proposals and target objective functions. Currently, only one appearance of the product symbol (*) can be parsed "                       ...
                        + "in the string value of scaleFactor. The presence of other mathematical symbols or multiple appearances of the product "                  ...
                        + "symbol will lead to a simulation crash. Also, note that the prescription of an acceptance range specified by the input "                 ...
                        + "variable 'AccRange' will lead to dynamic modification of the initial input value of scaleFactor throughout sampling, "                   ...
                        + "for adaptiveUpdateCount times. "                                                                                                         ...
                        + "The default scaleFactor string-value is 'gelman' (for all proposals), which is subsequently converted to 2.38/sqrt(ndim)."               ...
                        ;
        end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

        function set(self, scaleFactor)
            if isempty(scaleFactor)
                self.str = self.defStr;
                self.val = self.defVal;
            else
                self.str = strtrim(scaleFactor);
            end
        end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

        % ATTN: This method also assigns the value of ScaleFactor. It MUST be executed by all images.
        function Err = checkForSanity(self, Err, methodName)
            FUNCTION_NAME = "@checkForSanity()";

            % First convert the scaleFactor string to real value:

            String.value = strrep(self.str, " ", "");   % remove the white spaces
            String.value = lower(strtrim(String.value));

            if length(String.value) == 0
                Err.occurred    = true;
                Err.msg         = Err.msg                                                                                   ...
                                + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred. "                                    ...
                                + "The input string value (" + self.str + ") for the variable scaleFactor "                 ...
                                + "is empty. Make sure the input string follows the syntax rules of "                       ...
                                + methodName + " for this variable. Otherwise drop it from the input list. "                ...
                                + methodName + " will automatically assign an appropriate value to it." + newline + newline ...
                                ;
                return
            end

            % Now split the string by "*" to real coefficient and character (gelman) parts for further evaluations

            String.Parts    = split(String.value, "*");
            self.val        = 1;
            for i = 1 : length(String.Parts)
                if lower(String.Parts(i)) == "gelman"
                    self.val = self.val * self.defVal;
                else
                    number   = str2double(String.Parts(i));
                    if isnan(number)
                        Err.occurred    = true;
                        Err.msg         = Err.msg                                                                                                           ...
                                        + self.CLASS_NAME + FUNCTION_NAME + ": Error occurred while reading real number. "                                  ...
                                        + "The input string value for the variable scaleFactor (" + self.str + ") does not appear to follow "               ...
                                        + "the standard syntax rules of " + methodName + " for this variable. '" + String.Parts(i)                          ...
                                        + "' cannot be parsed into any meaningful token. Please correct the input value, or drop it from the input list, "  ...
                                        + "in which case, " + methodName + " will automatically assign an appropriate value to it." + newline + newline     ...
                                        ;
                        return
                    end
                    self.val = self.val * number;
                end
            end

            % Now check if the real value is positive

            if self.val <= 0
                Err.occurred = true;
                Err.msg = Err.msg                                                                                           ...
                        + self.CLASS_NAME + FUNCYION_NAME + ": Error occurred." + newline                                   ...
                        + "The input string value (" + self.str + ") translates to a negative real value: "                 ...
                        + num2str(self.val) + ". "                                                                          ...
                        + "Make sure the input string follows the syntax rules of " + methodName + " for this variable. "   ...
                        + "Otherwise drop it from the input list. " + methodName                                            ...
                        + " will automatically assign an appropriate value to it." + newline + newline                      ...
                        ;
                return
            end
        end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

    end

%***********************************************************************************************************************************
%***********************************************************************************************************************************

end
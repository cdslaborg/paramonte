cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('[loclist, namlist] = pm.str.locname(["library", "paramonte"], "paramonte")')
                [loclist, namlist] = pm.str.locname(["library", "paramonte"], "paramonte")
pm.matlab.show('loclist')
pm.matlab.show( loclist )
pm.matlab.show('namlist')
pm.matlab.show( namlist )
assert(namlist == "paramonte")
assert(loclist == 2)

pm.matlab.show()
pm.matlab.show('[loclist, namlist] = pm.str.locname(["library", "paramonte"], ["paramonte", 1])')
                [loclist, namlist] = pm.str.locname(["library", "paramonte"], ["paramonte", 1])
pm.matlab.show('loclist')
pm.matlab.show( loclist )
pm.matlab.show('namlist')
pm.matlab.show( namlist )
assert(all(namlist == ["paramonte", "library"]))
assert(all(loclist == [2, 1]))
#!/usr/bin/env python3
def f2F(dir):
    """
    Convert all .f90 extensions in the input dir to .F90.
    Note that dir must be defined with respect to the containing folder of f2F.py script file.
    """
    import os
    import glob
    import shutil
    fileList = glob.glob(os.path.join(dir, "**", "*.f90"), recursive = True)
    print(len(fileList))
    for file in fileList:
        split_tup = os.path.splitext(file)
        newfile = split_tup[0] + ".F90"
        print(file + " --> " + newfile)
        shutil.move(file, newfile)
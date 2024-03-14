/*
####################################################################################################################################
####################################################################################################################################
####                                                                                                                            ####
####    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           ####
####                                                                                                                            ####
####    Copyright (C) 2012-present, The Computational Data Science Lab                                                          ####
####                                                                                                                            ####
####    This file is part of the ParaMonte library.                                                                             ####
####                                                                                                                            ####
####    LICENSE                                                                                                                 ####
####                                                                                                                            ####
####       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          ####
####                                                                                                                            ####
####################################################################################################################################
####################################################################################################################################
*/
#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include <stdio.h>

int32_t pm_sys_isdirc(char* path){
    int32_t itis;
    struct stat info;
    if (stat(path, &info) != 0)
        itis = -1; // error occurred, e.g., path does not exist.
    else if (info.st_mode &S_IFDIR)
        itis = 1; // path is directory.
    else
        itis = 0; // path is not directory.
    return itis;
}

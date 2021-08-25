#include <stdarg.h>
#include <ctype.h>
#include "dir_functions.h"

#if defined(_WIN32)
    #include <windows.h>
    void mkdir(char *path, mode_t mode){
        _mkdir(path);
    }
#else
    #include <unistd.h>
    #include <sys/stat.h>
    #include <sys/types.h>
    void CopyFile(char* source, char* target, int sw){
        if(sw)  system(concat("cp -n", source, " ", target, NULL));
        else    system(concat("cp ", source, " ", target, NULL));
    }
#endif

void CopyFile(char* source, char* target){
    CopyFile(source, target, 0);
}

void init_dir(char *folder, char *runnumber){
    init_dir_with_subdir(folder, runnumber, "");
}

void init_dir_with_subdir(char *folder, char *runnumber, char *subdir){
    mkdir(folder, 0777);
    mkdir(concat(folder, "/", runnumber, NULL), 0777);
    mkdir(concat(folder, "/", runnumber, "/", subdir, NULL), 0777);
}
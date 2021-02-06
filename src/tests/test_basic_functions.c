#include <stdio.h>
#include "../basic_functions.h"

void test_concat(){
    printf("CONCAT: %s\n", concat("He","ll","oWor","ld","!",NULL)) ;
    printf("CONCAT: %s\n", concat("He","ll","\3",NULL,"ld","!")) ;
    printf("CONCAT: %s\n", concat("He","ll","oWor\0",NULL,"ld","!")) ;
    printf("CONCAT: %s\n", concat("He",NULL,"ll","","ld","!\3",NULL)) ;
    printf("CONCAT: %s\n", concat(NULL, "He","ll","oWor","!")) ;
    printf("CONCAT: %s\n", concat("He","ll","o!",NULL)) ;
    char* f="RR";
    
    printf("CONCAT: %s\n", concat("He",f,"o!",NULL)) ;
}



#include <stdio.h>
#include <string.h>
#include "taiga_test.h"
#include "tests/test_basic_functions.h"

int test_concat(){
    TAIGA_INIT_TEST("CONCAT");
    TAIGA_ASSERT_COMPARE("HelloWorld!", concat("He","ll","oWor","ld","!",NULL), "concat test 01");
    TAIGA_ASSERT_COMPARE("Hell\3", concat("He","ll","\3",NULL,"ld","!"), "concat test 02");
    TAIGA_ASSERT_COMPARE("HelloWor", concat("He","ll","oWor\0",NULL,"ld","!"), "concat test 03");
    TAIGA_ASSERT_COMPARE("He",  concat("He",NULL,"ll","","ld","!\3",NULL), "concat test 04");
    TAIGA_ASSERT_COMPARE("",  concat(NULL, "He","ll","oWor","!"), "concat test 05");
    TAIGA_ASSERT_COMPARE("Hello!",  concat("He","ll","o!",NULL), "concat test 06");
    char* f="RR";

    TAIGA_ASSERT_COMPARE("HeRRo!", concat("He",f,"o!",NULL), "concat test 07") ;
    return TAIGA_ASSERT_SUMMARY();
}


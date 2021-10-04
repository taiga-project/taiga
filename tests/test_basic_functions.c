#include <stdio.h>
#include <string.h>
#include "taiga_test.h"
#include "test_basic_functions.h"

#include "utils/basic_functions.h"

int test_concat(){
    TAIGA_INIT_TEST("CONCAT");
    TAIGA_ASSERT_COMPARE("HelloWorld!", concat("He","ll","oWor","ld","!",NULL), "concat test 01");
    TAIGA_ASSERT_COMPARE("Hell\3", concat("He","ll","\3",NULL,"ld","!"), "concat test 02");
    TAIGA_ASSERT_COMPARE("HelloWor", concat("He","ll","oWor\0",NULL,"ld","!"), "concat test 03");
    TAIGA_ASSERT_COMPARE("He",  concat("He",NULL,"ll","","ld","!\3",NULL), "concat test 04");
    TAIGA_ASSERT_COMPARE("",  concat(NULL, "He","ll","oWor","!"), "concat test 05");
    TAIGA_ASSERT_COMPARE("Hello!",  concat("He","ll","o!",NULL), "concat test 06");
    char* f="RR";

    TAIGA_ASSERT_COMPARE("HeRRo!", concat("He",f,"o!",NULL), "concat test 07");
    return TAIGA_ASSERT_SUMMARY();
}

int test_interpolate(){
    int N = 8;
    double x[] = {1,2,3,4,5,6,7,8};
    double x2[] = {8,7,6,5,4,3,2,1};
    double y[] = {2,3,5,8,13,21,34,55};
    double x0 = 2.7;
    double y0;

    TAIGA_INIT_TEST("INTERPOLATE");
    y0=linear_interpolate(x, N, y, N, x0);
    TAIGA_ASSERT_EQ(4.4, y0, "increasing order");
    y0=linear_interpolate(x2, N, y, N, x0);
    TAIGA_ASSERT_EQ(24.9, y0, "decreasing order");
    return TAIGA_ASSERT_SUMMARY();
}

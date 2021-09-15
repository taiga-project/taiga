#include <stdio.h>
#include <string.h>
#include "taiga_test.h"

double return_double(double x){
    return x;
}

int main(){

    TAIGA_INIT_TEST_NO_TITLE();
    TAIGA_ASSERT_ALMOST_EQ(0.0,-0.0,"almost equal");
    TAIGA_ASSERT_SUMMARY();

    TAIGA_INIT_TEST("SAMPLE TEST");
    double fifty_six_double = 56.0;
    double fifty_seven_double = 57.0;
    int fifty_six = 56;
    int fifty_seven = 57;
    char* fifty_six_string = "56";
    TAIGA_ASSERT_EQ(fifty_six_double, fifty_six,"double vs int (same)");
    TAIGA_ASSERT_EQ(fifty_six_double, fifty_seven,"double vs int (different)");
    TAIGA_ASSERT_EQ(fifty_six, fifty_six_double,"int vs double (same)");
    TAIGA_ASSERT_EQ(fifty_seven, fifty_six_double,"int vs double (different)");
    TAIGA_ASSERT_EQ(fifty_six_double, return_double(fifty_six),"double vs double function (same)");
    TAIGA_ASSERT_EQ(fifty_six_double, return_double(fifty_seven),"double vs double function (different)");
    TAIGA_ASSERT_EQ(fifty_six, return_double(fifty_six),"int vs double function (same)");
    TAIGA_ASSERT_EQ(fifty_six, return_double(fifty_seven),"int vs double function (different)");
    TAIGA_ASSERT_EQ(fifty_six, fifty_six_string,"int vs string");
    TAIGA_ASSERT_EQ(fifty_six_string, fifty_six,"string vs int");
    TAIGA_ASSERT_EQ("56.2", fifty_six_string,"different strings");
    TAIGA_ASSERT_EQ("Hello", "hello", "string test (different)");
    TAIGA_ASSERT_EQ(fifty_six_string, "56", "string test (same)");
    TAIGA_ASSERT_EQ("56", fifty_six_string, "string test reversed (same)");
    TAIGA_ASSERT_ALMOST_EQ(fifty_six,56.000005,"almost equal");
    TAIGA_ASSERT_ALMOST_EQ(fifty_six_double,56.000006,"almost not equal");
    return TAIGA_ASSERT_SUMMARY();
}

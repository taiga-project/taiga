#ifndef TAIGA_TEST_H
#define TAIGA_TEST_H

#include <stdio.h>
#include <string.h>

#define IS_PASSED 0
#define NUMBER_OF_TESTS number_of_tests
#define NUMBER_OF_FAILS number_of_fails

#define UNCHANGEABLE_VALUE 5
#define CHANGEABLE_VALUE foo

#define PRINT_VALUE(x) printf("%s", #x)

#define TAIGA_TEST_MESSAGE(test_name, status, expected, actual) { \
    ++number_of_tests;                                            \
    if (status == IS_PASSED){                                     \
        printf("\033[0;32m");                                     \
        printf("[   ok   ] ");                                    \
        printf("\033[0m");                                        \
        printf("%s\n", test_name);                                \
    }else{                                                        \
        ++number_of_fails;                                        \
        printf("\033[0;31m");                                     \
        printf("[ failed ] ");                                    \
        printf("\033[0m");                                        \
        printf("%s\n", test_name);                                \
        printf("\texpected: ");                                   \
        PRINT_VALUE(expected);                                    \
        printf("\n\tactual:   ");                                 \
        PRINT_VALUE(actual);                                      \
        printf("\n");                                             \
    }                                                             \
}

#define TAIGA_ASSERT_EQ(expected, actual, test_name) { \
    int is_failed = (expected!=actual);                \
    TAIGA_TEST_MESSAGE(test_name, is_failed, expected, actual); \
}


#define TAIGA_ASSERT_ALMOST_EQ(expected, actual, test_name) { \
    double diff = 1-actual/expected;                          \
    double tolerance = 1e-7;                                  \
    int is_failed = ((diff>-tolerance)&(diff<tolerance))?0:1; \
    TAIGA_TEST_MESSAGE(test_name, is_failed, expected, actual); \
}

#define TAIGA_ASSERT_SUMMARY() ({ \
    printf("==========\n");      \
    if (NUMBER_OF_FAILS == 0) {  \
        printf("\033[0;32m");    \
        printf("[ PASSED ] ");   \
        if (NUMBER_OF_TESTS == 1) { \
            printf("1 test");    \
        } else {                 \
            printf("%d tests", NUMBER_OF_TESTS); \
        }                        \
        printf("\n");            \
    } else {                     \
        printf("\033[0;31m");    \
        printf("[ FAILED ] ");   \
        if (NUMBER_OF_TESTS == 1) { \
            printf("1 test: 1 failed, 0 passed"); \
        } else {                 \
            printf("%d tests: %d failed, %d passed", \
            NUMBER_OF_TESTS, NUMBER_OF_FAILS,     \
            (NUMBER_OF_TESTS-NUMBER_OF_FAILS));   \
        }                        \
        printf("\n");            \
    }                             \
    NUMBER_OF_FAILS;              \
})

#define TAIGA_INIT_TEST() { \
       NUMBER_OF_TESTS = 0;  \
       NUMBER_OF_FAILS = 0;  \
}

int number_of_tests = 0;
int number_of_fails = 0;

#endif //TAIGA_TEST_H

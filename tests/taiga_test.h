#ifndef TAIGA_TEST_H
#define TAIGA_TEST_H

#define PASSED 0
#define PASSED_DETAIL -1
#define FAILED 1
#define TAIGA_ALMOST_ZERO 1e-300

#define NUMBER_OF_TESTS number_of_tests
#define NUMBER_OF_FAILS number_of_fails

#define PRINT_VALUE(x) {                                          \
    const int string_size = 50;                                   \
    char double_x_str[string_size];                               \
    double double_x;                                              \
    snprintf(double_x_str, string_size, "%.16lg", x);             \
    sscanf(double_x_str, "%lf", &double_x);                       \
    if (TAIGA_ALMOST_EQ(double_x, TAIGA_ALMOST_ZERO) == FAILED) { \
        printf("%s", double_x_str);                               \
    } else if (sizeof(x) == sizeof(int)) {                        \
        char int_x_str[string_size];                              \
        snprintf(int_x_str, string_size, "%ld", x);               \
        printf("%d", x);                                          \
    } else if(x) {                                                \
        int printf_len = printf("\"%s\"", x);                     \
    } else {                                                      \
        printf("%lf", x);                                         \
    }                                                             \
}

#define TAIGA_TEST_MESSAGE(test_name, status, expected, actual) { \
    ++number_of_tests;                                            \
    if (status == PASSED) {                                       \
        printf("\033[0;32m");                                     \
        printf("[   ok   ] ");                                    \
        printf("\033[0m");                                        \
        printf("%s\n", test_name);                                \
    } else if (status == PASSED_DETAIL) {                         \
        printf("\033[0;32m");                                     \
        printf("[   ok   ] ");                                    \
        printf("\033[0m");                                        \
        printf("%s\n", test_name);                                \
        printf("\t\t\texpected:   ");                             \
        PRINT_VALUE(expected);                                    \
        printf("\n\t\t\tactual:     ");                           \
        PRINT_VALUE(actual);                                      \
        if (expected != 0) {                                      \
            printf("\n\t\t\tdifference: ");                       \
            PRINT_VALUE(actual-expected);                         \
        }                                                         \
        printf("\n");                                             \
    } else {                                                      \
        ++number_of_fails;                                        \
        printf("\033[0;31m");                                     \
        printf("[ failed ] ");                                    \
        printf("\033[0m");                                        \
        printf("%s\n", test_name);                                \
        printf("\t\t\texpected: ");                               \
        PRINT_VALUE(expected);                                    \
        printf("\n\t\t\tactual:   ");                             \
        PRINT_VALUE(actual);                                      \
        if (expected != 0) {                                      \
            printf("\n\t\t\tdifference: ");                       \
            PRINT_VALUE(actual-expected);                         \
        }                                                         \
        printf("\n");                                             \
    }                                                             \
}

#define TAIGA_ASSERT_EQ(expected, actual, test_name) {          \
    int is_failed = ((expected)!=(actual));                     \
    TAIGA_TEST_MESSAGE(test_name, is_failed, expected, actual); \
}

#define TAIGA_ASSERT_COMPARE(expected, actual, test_name) {     \
    int is_failed = (strcmp(expected, actual)!=0);              \
    TAIGA_TEST_MESSAGE(test_name, is_failed, expected, actual); \
}

#define TAIGA_ALMOST_EQ(diff, tolerance) ({ \
    ((diff>-tolerance)&(diff<tolerance))? PASSED_DETAIL:FAILED; \
})

#define TAIGA_ASSERT_ALMOST_EQ(expected, actual, test_name) {   \
    double diff = 1-actual/expected;                            \
    double tolerance = 1e-7;                                    \
    int is_failed = TAIGA_ALMOST_EQ(diff, tolerance);           \
    TAIGA_TEST_MESSAGE(test_name, is_failed, expected, actual); \
}

#define TAIGA_ASSERT_ALMOST_EQ_MAX_DIFF(expected, actual, tolerance, test_name) { \
    double diff = actual-expected;                                                \
    int is_failed = TAIGA_ALMOST_EQ(diff, tolerance);                             \
    TAIGA_TEST_MESSAGE(test_name, is_failed, expected, actual);                   \
}

#define TAIGA_ASSERT_SUMMARY() ({                    \
    if (NUMBER_OF_TESTS < 0) {                       \
        printf("----------\n");                      \
    }                                                \
    if (NUMBER_OF_FAILS == 0) {                      \
        printf("\033[0;32m");                        \
        printf("[ PASSED ] ");                       \
        if (NUMBER_OF_TESTS == 1) {                  \
            printf("1 test");                        \
        } else {                                     \
            printf("%d tests", NUMBER_OF_TESTS);     \
        }                                            \
        printf("\n");                                \
    } else {                                         \
        printf("\033[0;31m");                        \
        printf("[ FAILED ] ");                       \
        if (NUMBER_OF_TESTS == 1) {                  \
            printf("1 test: 1 failed, 0 passed");    \
        } else {                                     \
            printf("%d tests: %d failed, %d passed", \
            NUMBER_OF_TESTS, NUMBER_OF_FAILS,        \
            (NUMBER_OF_TESTS-NUMBER_OF_FAILS));      \
        }                                            \
        printf("\n");                                \
    }                                                \
    printf("\033[0m");                               \
    printf("==========\n");                          \
    NUMBER_OF_FAILS;                                 \
})

#define TAIGA_INIT_TEST(title) {                        \
       NUMBER_OF_TESTS = 0;                             \
       NUMBER_OF_FAILS = 0;                             \
       printf("++++++++++ %s (%s)\n", title, __FILE__); \
}

#define TAIGA_INIT_TEST_NO_TITLE() { \
    TAIGA_INIT_TEST("");             \
}

int number_of_tests;// = 0;
int number_of_fails;// = 0;

#endif //TAIGA_TEST_H
#include <stdio.h>
#include <string.h>
#include "taiga_test.h"
#include "test_basic_functions.h"
#include "test_solver.h"
#include "test_bspline.h"
#include "test_bezier.h"

int main(){
    return /*test_concat() +
    test_interpolate() +
    test_to_string() +
    test_bspline() +
    test_solver() +*/
    test_bezier();
}

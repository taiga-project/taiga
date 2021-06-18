#include "taiga_test.c"

int main(){
    TAIGA_ASSERT_EQ(56,56,"Árvíztűrő tükörfúrógép");
    TAIGA_ASSERT_EQ("56","56","Árvíztűrő tükörfúrógép");
    TAIGA_ASSERT_ALMOST_EQ(56,56.000005,"Árvíztűrő tükörfúrógép");
    TAIGA_ASSERT_ALMOST_EQ(56,56.000006,"Árvíztűrő tükörfúrógép");
    return TAIGA_ASSERT_SUMMARY();
}

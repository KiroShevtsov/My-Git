#include <gtest/gtest.h>
#include "../Source/Slae.hxx"
TEST(Sweep, SweepTest){
    // std::vector<double> solve = {1, 4, 5};
    
    // Vector x = SweepMethod({0, 0}, {1, 1, 1}, {0, 0}, solve);
    // ASSERT_EQ(x, Vector({1, 4, 5})) << x;

    std::vector<double> solve = {
        115.0 / 297,     
        20.0 / 27,    
        4077.0 / 2673,      
        -1218.0 / 2673
        
    };
    Vector d({5, 9, 7, 2});
    Vector x = SweepMethod({4, 1, 4}, {11, 8, 5, 9}, {1, 1, 3}, d);
    const double tolerance = 1e-3;
    for(std::size_t i = 0; i < d.dim(); ++i){
        EXPECT_NEAR(x[i], solve[i], tolerance);
    }
}

#include <Slae.hxx>
#include <iostream>
#include <map>
int main(){
    try{
        auto s = std::pair<int, int>(0, 0);
        auto f = std::pair<int, int>(0, 1);
        auto g = std::pair<int, int>(1, 1);
        auto h = std::pair<int, int>(2, 1);
        auto k = std::pair<int, int>(2, 2);
        std::map<std::pair<size_t, size_t>, double> exp = {{s, 1}, {f, 2}, {g, 4}, {h, 2}, {k, 6}};
        CsrMatrix m(3, 3, exp);
        CsrMatrix sum = m * m * m;
        std::cout << sum(1, 1);
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
}
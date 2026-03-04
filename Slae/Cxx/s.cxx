#include <Slae.hxx>
#include <iostream>
int main(){
    Matrix m1(2, 2, {1, 2, 3, 4});
    std::cout << m1 + m1;
}
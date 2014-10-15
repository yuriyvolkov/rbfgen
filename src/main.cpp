#include <iostream>
#include <rbf.h>
int main(){

    ainet::RBFnet<unsigned long> net(10,10);

    std::cout << net << std::endl;
    return 0;
}

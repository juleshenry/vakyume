#include <iostream>

int funct(int num) 
{
    std::cout << "Num = " << num << std::endl;
    return 0;
}

extern "C" {
    int my_function(int a)
    {
        return funct(a);
    }
}
// NeutrinoSignificance.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <exception>
#include "../SignificanceLib/SignificanceMC.h"

int main(int argc, char* argv[])
{
    GridPoint observed{ 0,0,0 };
    for (int i = 0; i < argc - 1; i++)
    {
        size_t word;
        std::stringstream s{ argv[i + 1] };
        s >> word;
        observed[i] = word;
    }

    SignificanceMC mc{ observed };
    double p = mc.GetPValue();
    std::cout << "Returned p = " << p << std::endl;

}


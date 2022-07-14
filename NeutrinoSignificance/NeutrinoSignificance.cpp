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
    SignificanceMC::PResult pResult = mc.GetPValue(20, 0.005);
    std::cout << "pValue:  " << pResult.p << std::endl;
    std::cout << "sigmaP = " << pResult.sigmaP << std::endl;
    std::cout << "1 sigma confidence interval: {" << pResult.pLow << " - " << pResult.pHigh << "}" << std::endl;
    std::cout << "Significance (in standard deviations): " << pResult.z << std::endl;
    std::cout << "Significance confidence interval: {" << pResult.zLow << " - " << pResult.zHigh << "}" << std::endl;
}


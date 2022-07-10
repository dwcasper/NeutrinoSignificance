// NeutrinoSignificance.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <exception>
#include "../SignificanceLib/Significance.h"

int main(int argc, char* argv[])
{
    size_t observed{ 0 };
    for (int i = 0; i < argc - 1; i++)
    {
        size_t word;
        std::stringstream s{ argv[i + 1] };
        s >> word;
        if (i == 0 && word > (Significance::grid0 >> Significance::shift0) || i == 1 && word > (Significance::grid1 >> Significance::shift1) || i == 2 && word >> (Significance::grid2 >> Significance::shift2))
            throw std::out_of_range("Input measurement out of allowed range.");
        observed |= (word << (i == 0 ? Significance::shift0 : (i == 1 ? Significance::shift1 : Significance::shift2)));
    }

    std::cout << "Decoded inputs:" << std::endl;
    std::cout << ((observed & Significance::grid2) >> Significance::shift2) << std::endl;
    std::cout << ((observed & Significance::grid1) >> Significance::shift1) << std::endl;
    std::cout << ((observed & Significance::grid0) >> Significance::shift0) << std::endl;

    Significance* significance = new Significance(1.0e-8, observed);

    std::cout << "Tail = " << significance->GetTail() << std::endl;
    std::cout << "q0Observed = " << significance->Get_q0Observed() << std::endl;
    char c;
    std::cout << "Enter a character to terminate" << std::endl;
    std::cin >> c;

    delete significance;

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

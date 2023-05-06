#include <iostream>
#include "func.h"
#include "func.cpp"

int main()
{
    Solution sol = Solution();

    std::string FileName = "C:\\Users\\dell\\Desktop\\workspace\\.vscode\\ComputationalMethod\\lab\\lab5\\point.txt";
    int n = 21;

    std::vector<double> Ms(n);
    sol.Spline(FileName, Ms, 0, 0);


    std::vector<double> Ms2(n);
    sol.Spline(FileName, Ms2, 0, 1);
    return 0;
}
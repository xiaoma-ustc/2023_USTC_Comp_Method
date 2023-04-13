#include"func.h"
#include"func.cpp"

int main()
{
    
    std::vector<std::vector<double>> a{{1.0/9.0, 1.0/8.0, 1.0/7.0, 1.0/6.0, 1.0/5.0},
                                        {1.0/8.0, 1.0/7.0, 1.0/6.0, 1.0/5.0, 1.0/4.0},
                                        {1.0/7.0, 1.0/6.0, 1.0/5.0, 1.0/4.0, 1.0/3.0},
                                        {1.0/6.0, 1.0/5.0, 1.0/4.0, 1.0/3.0, 1.0/2.0},
                                        {1.0/5.0, 1.0/4.0, 1.0/3.0, 1.0/2.0, 1.0/1.0}};
    Solution sol = Solution(5, a, "A1.md");
    sol.MinEigIPM(5, 1e-5);

    std::vector<std::vector<double>> b{{4, -1, 1, 3},
                                        {16, -2, -2, 5},
                                        {16, -3, -1, 7},
                                        {6, -4, 2, 9}};
    Solution sol2 = Solution(4, b, "A2.md");
    sol2.MinEigIPM(4, 1e-5);

    return 0;
}
#ifndef FUNC
#define FUNC

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
class Solution
{
private:
public:

    Solution();
    void Thomas(int n, std::vector<std::vector<double>>& Mn3, std::vector<double>& b, std::vector<double>& x);
    void Spline(std::string FileName, std::vector<double>& Ms, int condition, int flag);

};

#endif
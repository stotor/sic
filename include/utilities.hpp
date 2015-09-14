#ifndef utilities_hpp
#define utilities_hpp

#include <vector>
#include <fstream>

#define PI (4 * atan(1.0))

void write_data(std::vector<double> &array, std::ofstream &file, int n_values);

int mod(int a, int b);

#endif /* utilities_hpp */

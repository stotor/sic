#include <vector>
#include <fstream>
#include "utilities.hpp"

void write_data(std::vector<double> &array, std::ofstream &file, int n_values)
{
  for (int i = 0; i < n_values; i++) {
    file << array[i] << " ";
  }
  file << std::endl;
  return;
}

int mod(int a, int b)
{
  return (a%b + b) % b;
}

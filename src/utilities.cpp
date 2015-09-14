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

void write_data_tavg(std::vector<double> &array, std::vector<double> &array_old, std::ofstream &file, 
		     int n_values)
{
  for (int i = 0; i < n_values; i++) {
    file << (array[i] + array_old[i]) / 2.0 << " ";
  }
  file << std::endl;
  return;
}

void calc_tavg_arry(std::vector<double> &array_tavg, std::vector<double> &array, 
		    std::vector<double> &array_old, 
		    int n_values)
{
  for (int i = 0; i < n_values; i++) {
    array_tavg[i] = (array[i] + array_old[i]) / 2.0;
  }
  return;
}

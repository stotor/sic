#include <vector>
#include <fstream>
#include "utilities.hpp"

double sum_of_squares(std::vector<double> &array, int n_values)
{
  double sum = 0.0;
  for (int i = 0; i < n_values; i++) {
    sum = sum + array[i] * array[i];
  }
  return sum;
}

void save_old_values(std::vector<double> &array, std::vector<double> &array_old,
		     int n_values)
{
  for (int i = 0; i < n_values; i++) {
    array_old[i] = array[i];
  }
  return;
}

void write_data(std::vector<double> &array, std::ofstream &file, int n_values)
{
  for (int i = 0; i < n_values; i++) {
    //    file << array[i];
    file.write((char *) &array[i], sizeof(double));
  }
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
    //file.write((char *) &(array[i] + array_old[i]) / 2.0, sizeof(double));
    file << (array[i] + array_old[i]) / 2.0;
  }
  return;
}

void calc_tavg_array(std::vector<double> &array_tavg, std::vector<double> &array,
		     std::vector<double> &array_old, 
		     int n_values)
{
  for (int i = 0; i < n_values; i++) {
    array_tavg[i] = (array[i] + array_old[i]) / 2.0;
  }
  return;
}

void half_int_to_int(std::vector<double> &field, 
		     std::vector<double> &shifted_field, 
		     int n_g)
{
  for (int i = 0; i < n_g; i++) {
    shifted_field[i] = (field[mod((i-1),n_g)] + field[i]) / 2.0;
  }
  return;
}

#ifndef utilities_hpp
#define utilities_hpp

#include <vector>
#include <fstream>

#include "mpi.h"

#define PI (4 * atan(1.0))

void data_to_file(std::vector<double> data, std::string filename);

void save_old_values_int(std::vector<int> &array, std::vector<int> &array_old,
			 int n_values);

void save_old_values_double(std::vector<double> &array, std::vector<double> &array_old,
			    int n_values);

void write_data(std::vector<double> &array, std::ofstream &file, int n_values);

int mod(int a, int b);

void half_int_to_int(std::vector<double> &field, 
		     std::vector<double> &shifted_field, 
		     int n_g);

double random_double(void);

void sum_array_to_root(double *array, int n_values, MPI_Comm COMM, int my_rank);

#endif /* utilities_hpp */

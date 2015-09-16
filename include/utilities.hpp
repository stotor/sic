#ifndef utilities_hpp
#define utilities_hpp

#include <vector>
#include <fstream>

#define PI (4 * atan(1.0))

void write_data(std::vector<double> &array, std::ofstream &file, int n_values);

int mod(int a, int b);

void write_data_tavg(std::vector<double> &array, std::vector<double> &array_old, std::ofstream &file, 
		     int n_values);

void calc_tavg_arry(std::vector<double> &array_tavg, std::vector<double> &array, 
		    std::vector<double> &array_old, 
		    int n_values);

void half_int_to_int(std::vector<double> &field, 
		     std::vector<double> &shifted_field, 
		     int n_g);

#endif /* utilities_hpp */

#ifndef fields_hpp
#define fields_hpp

#include <vector>
#include <string>
#include <fstream>

class Field {
public:
  Field(int n_g, std::string field_name) : field(n_g), field_old(n_g), field_tavg(n_g)
  {
    this->n_g = n_g;
    output_stream.open(field_name.c_str(), std::ios::binary);
  }
  
  ~Field() {
    output_stream.close();
  }
  
  // Attributes
  int n_g;
  std::vector<double> field, field_old, field_tavg, energy_history;
  std::ofstream output_stream;
  void write_field();
  void write_field_tavg();
  void set_field_to_zero();
  void calculate_energy();
  void calculate_energy_tavg();
};


// Initialize e_x
void initialize_e_x(std::vector<double> &rho, std::vector<double> &e_x, 
		    double dx, int n_g);

void e_x_poisson_solve(std::vector<double> &rho, std::vector<double> &e_x, 
		       double dx, int n_g, std::vector<double> &phi);

void initialize_transverse_em_fields(std::vector<double> &e_y, 
				     std::vector<double> &b_z, int n_g, 
				     double dx, double e_y_1, double b_z_1,
				     int mode);

void advance_b_z(std::vector<double> &b_z, std::vector<double> &e_y, double dt, 
		 double dx, int n_g);

void advance_e_x(std::vector<double> &e_x, std::vector<double> &j_x, double dt, 
		 double dx, int n_g);

void advance_e_y(std::vector<double> &e_y, std::vector<double> &b_z, 
		 std::vector<double> &j_y, double dt, 
		 double dx, int n_g);

double sum_of_squares(std::vector<double> &array, int n_values);

#endif /* fields_hpp */

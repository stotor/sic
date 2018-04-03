#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>

#include "fields.hpp"
#include "utilities.hpp"

void Field::write_energy_history()
{
  data_to_file(energy_history, (field_name+"_ene"));
  return; 
}

void Field::calculate_energy()
{
  energy_history.push_back(0.5 * sum_of_squares(field, n_g));
  return;
}

void Field::write_field()
{
  write_data(field, output_stream, n_g);
  return;
}

void Field::set_field_to_zero()
{
  for (int i = 0; i < n_g; i++) {
    field[i] = 0.0;
  }
  return;
}

void initialize_transverse_em_fields(std::vector<double> &e_y, 
				     std::vector<double> &b_z, int n_g, 
				     double dx, double e_y_1, double b_z_1,
				     int mode)
{
  double k = 2.0 * PI * mode / (n_g * dx);
  for (int i = 0; i < n_g; i++) {
    e_y[i] = e_y_1 * cos(k * i * dx);
    b_z[i] = b_z_1 * cos(k * (i + 0.5) * dx);
  }
  return;
}

void initialize_fields_weibel(std::vector<double> &e_y, 
			      std::vector<double> &b_z,
			      int n_g, 
			      double dx,
			      int mode_max,
			      double amplitude)
{
  double k, phase;
  srand(0);
  for (int i = 0; i < n_g; i++) {
    e_y[i] = 0.0;
    b_z[i] = 0.0;
  }
  for (int mode = 1; mode <= mode_max; mode++) {  
    k = 2.0 * PI * mode / (n_g * dx);
    phase = 2.0 * PI * random_double();
    for (int i = 0; i < n_g; i++) {
      b_z[i] += amplitude * cos(k * (i + 0.5) * dx + phase);
    }
  }
  return;
}

void initialize_e_x(std::vector<double> &rho, std::vector<double> &e_x, 
		    double dx, int n_g)
{
  double e_x_ave = 0.0;
  e_x[0] = 0.0;
  for (int i = 1; i < n_g; i++) {
    e_x[i] = e_x[i-1] + dx * rho[i];
  }

  // Subtract average field
  for (int i = 0; i < n_g; i++) {
    e_x_ave = e_x_ave + e_x[i];
  }
  e_x_ave = e_x_ave / double(n_g);
  for (int i = 0; i < n_g; i++) {
    e_x[i] = e_x[i] - e_x_ave;
  }

  return;
}

void e_x_poisson_solve(std::vector<double> &rho, std::vector<double> &e_x, 
		       double dx, int n_g, std::vector<double> &phi)
{
  // Calculate electric potential
  phi[0] = 0;
  for (int i = 0; i < n_g; i++) {
    phi[0] = phi[0] + (i+1) * rho[i];
  }
  phi[0] = phi[0] / double(n_g);
  phi[1] = rho[0] + 2.0 * phi[0];

  for (int i = 2; i < n_g; i++) {
    phi[i] = rho[i-1] + 2.0 * phi[i-1] - phi[i-2];
  }

  // Calculate e_x from phi
  for (int i = 0; i < (n_g-1); i++) {
    e_x[i] = (-1.0 * dx) * (phi[i+1] - phi[i]);
  }
  e_x[n_g-1] = (-1.0 * dx) * (phi[0] - phi[n_g-1]);
  return;
}

void advance_b_z(std::vector<double> &b_z, std::vector<double> &e_y, double dt, 
		 double dx, int n_g)
{
  for (int i = 0; i < (n_g - 1); i++) {
    b_z[i] = b_z[i] - (dt / dx) * (e_y[i+1] - e_y[i]);
  }
  b_z[n_g-1] = b_z[n_g-1] - (dt / dx) * (e_y[0] - e_y[n_g-1]);
  return;
}

void advance_e_x(std::vector<double> &e_x, std::vector<double> &j_x, double dt, 
		 double dx, int n_g)
{
  for (int i = 0; i < n_g; i++) {
    e_x[i] = e_x[i] - dt * j_x[i];
  }
  return;
}

void advance_e_y(std::vector<double> &e_y, std::vector<double> &b_z, 
		 std::vector<double> &j_y, double dt, 
		 double dx, int n_g)
{
  e_y[0] = e_y[0] - (dt / dx) * (b_z[0] - b_z[n_g-1]) - dt * j_y[0];
  for (int i = 1; i < n_g; i++) {
    e_y[i] = e_y[i] - (dt / dx) * (b_z[i] - b_z[i-1]) - dt * j_y[i];
  }
  return;
}

// void advance_b_z_2d(std::vector<double> &b_z, std::vector<double> &e_y, double dt, 
// 		    double dx, int n_g)
// {
//   for (int j = 0; j < (n_g - 1); j++) {
//     for (int i = 0; i < (n_g - 1); i++) {
//       b_z[j, i] = b_z[j, i] - (dt / dx) * (e_y[j, i+1] - e_y[j, i] - e_x[j+1, i] - e_x[j, i]);
//     }
//     b_z[j, n_g-1] = b_z[j, n_g-1] - (dt / dx) * (e_y[j, 0] - e_y[j, n_g-1] - e_x[j+1, i] - e_x[j, i]);
//   }
//   for (int i = 0; i < (n_g - 1); i++) {
//     b_z[n_g-1, i] = b_z[n_g-1, i] - (dt / dx) * (e_y[n_g-1, i+1] - e_y[n_g-1, i] - e_x[0, i] - e_x[n_g-1, i]);
//   }
//   b_z[n_g-1, n_g-1] = b_z[n_g-1, n_g-1] - (dt / dx) * (e_y[n_g-1, 0] - e_y[n_g-1, n_g-1] - e_x[0, n_g-1] - e_x[n_g-1, n_g-1]);
//   return;
// }

// void advance_e_x_2d(std::vector<double> &e_x, std::vector<double> &j_x, double dt, 
// 		    double dx, int n_g)
// {
//   for (int j = 0; j < n_g; j++) {
//       for (int i = 0; i < n_g; i++) {
// 	  e_x[j, i] = e_x[j, i] - (dt / dx) * (b_z[j+1, i] - b_z[j, i]) - dt * j_x[j, i];
//       }
//   }
//   return;
// }

// void advance_e_y_2d(std::vector<double> &e_y, std::vector<double> &b_z, 
// 		    std::vector<double> &j_y, double dt, 
// 		    double dx, int n_g)
// {
//   //  e_y[0] = e_y[0] - (dt / dx) * (b_z[0] - b_z[n_g-1]) - dt * j_y[0];
//   for (int j = 1; j < n_g; j++) {
//     for (int i = 1; i < n_g; i++) {
//       e_y[j, i] = e_y[j, i] - (dt / dx) * (b_z[j, i] - b_z[j, i-1]) - dt * j_y[j, i];
//     }
//   }
//   return;
// }

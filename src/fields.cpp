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

void add_weibel_perturbation(std::vector<double> &b_z,
			     int n_g,
			     double dx)
{
  double k, phase, pi;
  int mode_max = 64;
  double amplitude = pow(10.0, -6.0);
  pi = 4.0 * atan(1.0);
  srand(0);
  for (int mode = 1; mode <= mode_max; mode++) {
    k = 2.0 * pi * mode / (n_g * dx);
    phase = 2.0 * pi * random_double();
    for (int i = 0; i < n_g; i++) {
      b_z[i] += amplitude * cos(k * (i + 0.5) * dx + phase);
    }
  }
  return;
}

void initialize_transverse_em_fields(std::vector<double> &e_y,
				     std::vector<double> &e_z,
				     std::vector<double> &b_x,
				     std::vector<double> &b_y,
				     std::vector<double> &b_z,
				     int n_g, 
				     double dx,
				     int simulation_type)
{
  if (simulation_type == -1 or simulation_type == 0 or simulation_type == 1) {
    for (int i = 0; i < n_g; i++) {
      e_y[i] = 0.0;
      e_z[i] = 0.0;
    
      b_x[i] = 0.0;
      b_y[i] = 0.0;
      b_z[i] = 0.0;
    }

    if (simulation_type == 1) {
      add_weibel_perturbation(b_z, n_g, dx);
    }
  } else if (simulation_type == 2) {
    for (int i = 0; i < n_g; i++) {
      e_y[i] = 2.0 * 0.40214 * cos(0.44526860656 * i * dx);
      e_z[i] = 0.0;
    
      b_x[i] = 1.7;
      b_y[i] = 2.0 * 0.80428 * sin(0.44526860656 * (i + 0.5) * dx);
      b_z[i] = 0.0;
    }
  }
  return;
}

void initialize_e_x(std::vector<double> &rho, std::vector<double> &e_x, 
		    double dx, int n_g, bool gravity)
{
  double e_x_ave = 0.0;
  e_x[0] = 0.0;
  for (int i = 1; i < n_g; i++) {
    e_x[i] = e_x[i-1] + dx * rho[i];
    e_x_ave = e_x_ave + e_x[i];
  }
  
  e_x_ave = e_x_ave / double(n_g);
  
  for (int i = 0; i < n_g; i++) {
    e_x[i] = e_x[i] - e_x_ave;
  }
  if (gravity) {
    for (int i = 0; i < n_g; i++) {
      e_x[i] = -1.0 * e_x[i];
    }
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

void advance_b_y(std::vector<double> &b_y, std::vector<double> &e_z, double dt, 
		 double dx, int n_g)
{
  for (int i = 0; i < (n_g - 1); i++) {
    b_y[i] = b_y[i] + (dt / dx) * (e_z[i+1] - e_z[i]);
  }
  b_y[n_g-1] = b_y[n_g-1] + (dt / dx) * (e_z[0] - e_z[n_g-1]);
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
		 double dx, int n_g, bool gravity)
{
  if (gravity) {
    for (int i = 0; i < n_g; i++) {
      e_x[i] = e_x[i] + dx * j_x[i];
    }
  } else {
    for (int i = 0; i < n_g; i++) {
      e_x[i] = e_x[i] - dx * j_x[i];
    }
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

void advance_e_z(std::vector<double> &e_z, std::vector<double> &b_y,
		 std::vector<double> &j_z, double dt, 
		 double dx, int n_g)
{
  e_z[0] = e_z[0] + (dt / dx) * (b_y[0] - b_y[n_g-1]) - dt * j_z[0];
  for (int i = 1; i < n_g; i++) {
    e_z[i] = e_z[i] + (dt / dx) * (b_y[i] - b_y[i-1]) - dt * j_z[i];
  }
  return;
}

#include <vector>
#include "fields.hpp"

// Initialize EM fields
void initialize_fields(std::vector<double> e_x, std::vector<double> e_y, std::vector<double> b_z, std::vector<double> j_x,
		       std::vector<double> j_y, int n_g)
{
  // EM wave
  // int mode = 4;
  // for (int i = 0; i < n_g; i++) {
  //   e_x[i] = 0.0;
  //   e_y[i] = cos(double(mode) * 2.0 * PI * (double(i) / double(n_g)));
  //   b_z[i] = cos(double(mode) * 2.0 * PI * (double(i) / double(n_g)));
  //   // BELOW FOR TESTING EM WAVE
  //   j_x[i] = 0.0;
  //   j_y[i] = 0.0;
  // }
  // Zero
  for (int i = 0; i < n_g; i++) {
    e_x[i] = 0.0;
    e_y[i] = 0.0;
    b_z[i] = 0.0;
    j_x[i] = 0.0;
    j_y[i] = 0.0;
  }
}

////////////////////////////////////////////////////////////////
// Fields
void save_old_values(std::vector<double> &array, std::vector<double> &array_old,
		     int n_values)
{
  for (int i = 0; i < n_values; i++) {
    array_old[i] = array[i];
  }
  return;
}

// Initialize e_x
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

// Not currently used
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
/**
 * SIC
 * - 1D electromagnetic code implementing the
 *   simplex-in-cell current deposit
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>

#include "speciesgroup.hpp"
#include "fields.hpp"
#include "utilities.hpp"

int main(int argc, char *argv[])
{
  ///////////////////////////////
  // INITIALIZATION
  // Define simulation parameters
  if (argc != 4) {
    std::cout << "Usage: ./sic <method> <n_ppc> <simulation_type>" << std::endl;
    return 1;
  }

  std::stringstream ss;

  ss << argv[1];
  int method;
  ss >> method;

  ss.str(std::string());
  ss.clear();
  double n_ppc;
  ss << argv[2];
  ss >> n_ppc;

  ss.str(std::string());
  ss.clear();
  int simulation_type; // 0: Electrostatic wave, 1: Electromagnetic wave, 2: Weibel
  ss << argv[3];
  ss >> simulation_type;
  
  int n_species;
  std::vector<double> u_x_drift, u_y_drift;
  double u_x_1, u_y_1, e_y_1, b_z_1;
  int mode = 1;
  int mode_max = 32;
  double amplitude;
  int n_t, n_g;
  double k, dx, dt;

  switch (simulation_type) {
  case 0:
    n_t = 200;
    n_g = 1000;
    dx = 0.5;
    dt = 0.2;
    k = 2.0 * PI * mode / (n_g * dx);
    n_species = 1;
    u_x_drift.push_back(0.0);
    u_y_drift.push_back(0.0);
    u_x_1 = 0.01;
    u_y_1 = 0.0;
    e_y_1 = 0.0;
    b_z_1 = 0.0;
    break;
  case 1:
    n_t = 200;
    n_g = 100;
    dx = 0.05;
    dt = 0.04;
    k = 2.0 * PI * mode / (n_g * dx);

    n_species = 1;
    u_x_drift.push_back(0.0);
    u_y_drift.push_back(0.0);
    u_x_1 = 0.0;
    u_y_1 = 0.0025;
    e_y_1 = u_y_1 * sqrt(1 + k*k);
    b_z_1 = u_y_1 * k;
    break;
  case 2:
    // NOT READY YET
    n_t = 200;
    n_g = 1000;
    dx = 0.05;
    dt = 0.04;
    k = 2.0 * PI * mode / (n_g * dx);

    n_species = 2;
    u_x_drift.push_back(0.0);
    u_x_drift.push_back(0.0);
    u_y_drift.push_back(-0.1);
    u_y_drift.push_back(0.1);
    u_x_1 = 0.0;
    u_y_1 = 0.0;
    e_y_1 = 0.0;
    b_z_1 = 0.0;
    break;
  case 3:
    n_t = 200;
    n_g = 100;
    dx = 0.05;
    dt = 0.04;
    k = 2.0 * PI * mode / (n_g * dx);
    n_species = 1;
    u_x_drift.push_back(0.0);
    u_y_drift.push_back(0.0);
    u_x_1 = 0.1;
    u_y_1 = 0.1;
    e_y_1 = 0.0;
    b_z_1 = 0.0;
    break;
  case 4:
    n_t = 500;
    n_g = 128;
    dx = 0.1;
    dt = 0.09;
    n_species = 2;
    u_x_drift.push_back(0.0);
    u_x_drift.push_back(0.0);
    u_y_drift.push_back(0.5);
    u_y_drift.push_back(-0.5);
    u_x_1 = 0.0;
    u_y_1 = 0.0;
    mode = 0;
    mode_max = 32;
    amplitude = pow(10.0, -6.0);
    break;
  }
  
  if (fmod((n_ppc * double(n_g)), 1.0) != 0.0) {
    std::cout << "Error: n_ppc * n_g is not an integer.";
    return 1;
  }
  
  // Initialize particles 
  SpeciesGroup particles(n_species, method, n_ppc, dt, dx, n_g);
  particles.initialize_species(n_ppc, u_x_drift, u_y_drift, mode, u_x_1, u_y_1);

  // Initialize fields
  Field e_x(n_g, "e_x");
  Field e_y(n_g, "e_y");
  Field b_z(n_g, "b_z");
  Field j_x(n_g, "j_x");
  Field j_y(n_g, "j_y");  
  Field rho(n_g, "rho");
  //particles.deposit_rho(rho.field, n_g);
  //initialize_e_x(rho.field, e_x.field, dx, n_g);

  //  initialize_transverse_em_fields(e_y.field, b_z.field, n_g, dx, e_y_1, b_z_1, mode);
  initialize_fields_weibel(e_x.field, e_y.field, b_z.field, n_g, dx, mode_max, amplitude);

  // Evolve u_x and u_y backwards in time to t = - 1/2 dt
  //  particles.initial_velocity_deceleration(e_x.field, e_y.field, b_z.field);
  //  advance_b_z(b_z.field, e_y.field, (-0.5 * dt), dx, n_g);

  ////////////////////////////////////////
  // MAIN LOOP
  for (int t = 0; t < n_t; t++) {
    // Print integer time diagnostics: e_x, e_y, rho
    e_x.write_field();
    e_y.write_field();
    b_z.write_field();
    particles.deposit_rho(rho.field, n_g);
    rho.write_field();
    e_x.calculate_energy();
    e_y.calculate_energy();
    b_z.calculate_energy();

    // Calculate u_x and u_y, at t+1/2
    particles.save_u_x_old();
    particles.save_u_y_old();
    particles.advance_velocity(e_x.field, e_y.field, b_z.field);

    // Update particle positions to t+1
    particles.save_x_old();
    particles.advance_x();

    // Calculate j_x and j_y at t+1/2
    particles.deposit_j_x(j_x.field);
    particles.deposit_j_y(j_y.field);
    
    // Print diagnostics for fields defined at half-integer times
    j_x.write_field();
    j_y.write_field();

    // Update magnetic field
    advance_b_z(b_z.field, e_y.field, dt/2.0, dx, n_g);
    advance_e_x(e_x.field, j_x.field, dt, dx, n_g);
    advance_e_y(e_y.field, b_z.field, j_y.field, dt, dx, n_g);
    advance_b_z(b_z.field, e_y.field, dt/2.0, dx, n_g);
    
    std::cout << t << std::endl;
  }
  // Write out energy diagnostics
  e_x.write_energy_history();
  e_y.write_energy_history();
  b_z.write_energy_history();
  particles.write_energy_history();

  return 0;
}   
// End of main

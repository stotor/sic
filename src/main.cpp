/**
 * 1DEM
 * - 1D electromagnetic code
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
  if (argc != 3) {
    std::cout << "Usage: ./1DEM <method> <n_ppc>" << std::endl;
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

  int n_t = 20;
  int n_g = 100;
  double dx = 0.5;
  double dt = 0.4;

  int simulation_type = 1; // 0: Electrostatic wave, 1: Electromagnetic wave, 2: Weibel
  int n_species;
  int n_p = n_g * n_ppc;
  double v_x_1, v_y_1, e_y_1, b_z_1;

  int mode = 1;
  double k = 2.0 * PI * mode / (n_g * dx);
  double omega;
  std::vector<double> u_drift_x;
  std::vector<double> u_drift_y;
  switch (simulation_type) {
  case 0:
    n_species = 1;
    v_x_1 = 0.0025;
    v_y_1 = 0.0;
    u_drift_x.push_back(0.0);
    u_drift_y.push_back(0.0);
    e_y_1 = 0.0;
    b_z_1 = 0.0;      
    break;
  case 1:
    n_species = 1;
    v_x_1 = 0.0;
    v_y_1 = 0.0025;
    omega = sqrt(1.0 + k*k);
    e_y_1 = v_y_1 * omega;
    b_z_1 = v_y_1 * k;
    break;
  case 2:
    // NOT READY YET
    n_species = 2;
    v_y_1 = 0.0025;
    break;
  }
  
  // Initialize particles 
  SpeciesGroup particles(n_species, method, n_p, dt, dx, n_g);
  particles.initialize_species(n_g, n_ppc, dx);  //v1, mode, 

  // Initialize fields
  Field e_x(n_g, "e_x");
  Field e_y(n_g, "e_y");
  Field b_z(n_g, "b_z");
  Field j_x(n_g, "j_x");
  Field j_y(n_g, "j_y");  
  Field rho(n_g, "rho");
  particles.deposit_rho(rho.field, n_g);
  initialize_e_x(rho.field, e_x.field, dx, n_g);

  initialize_transverse_em_fields(e_y.field, b_z.field, n_g, dx, e_y_1, b_z_1, mode);

  // Evolve b_z, u_x and u_y backwards in time to t = - 1/2 dt
  particles.initial_velocity_deceleration(e_x.field, e_y.field, b_z.field);
  advance_b_z(b_z.field, e_y.field, (-0.5 * dt), dx, n_g);

  ////////////////////////////////////////
  // MAIN LOOP
  for (int t = 0; t < n_t; t++) {
    // Print integer time diagnostics: e_x, e_y, rho
    e_x.write_field();
    e_y.write_field();
    rho.write_field();
    e_x.calculate_energy();
    e_y.calculate_energy();

    // Update magnetic field
    save_old_values(b_z.field, b_z.field_old, n_g);
    advance_b_z(b_z.field, e_y.field, dt, dx, n_g);
    calc_tavg_array(b_z.field_tavg, b_z.field, b_z.field_old, n_g);
    b_z.calculate_energy_tavg();

    // Calculate u_x and u_y, at t+1/2
    particles.save_u_x_old();
    particles.save_u_y_old();
    particles.advance_velocity(e_x.field, e_y.field, b_z.field_tavg);

    // Update particle positions to t+1
    particles.save_x_old();
    particles.advance_x();

    // Calculate j_x and j_y at t+1/2
    particles.deposit_j_x(j_x.field);
    particles.deposit_j_y(j_y.field);
    
    // Print diagnostics for fields defined at half-integer times
    b_z.write_field();
    j_x.write_field();
    j_y.write_field();

    // Print diagnostics that are averaged to integer times: field and particle energies, phase space

    // Advance e_x and e_y to t+1; calculate rho at t+1 for diagnostic purposes
    advance_e_x(e_x.field, j_x.field, dt, dx, n_g);
    advance_e_y(e_y.field, b_z.field, j_y.field, dt, dx, n_g);
    particles.deposit_rho(rho.field, n_g);

    std::cout << t << std::endl;
  }
  return 0;
}   
// End of main

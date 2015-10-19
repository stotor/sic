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
  std::cout << method << std::endl;

  ss.str(std::string());
  ss.clear();
  double n_ppc;
  ss << argv[2];
  ss >> n_ppc;
  std::cout << n_ppc << std::endl;

  int n_t = 101;
  int n_g = 10;
  double dx = 0.5;
  double dt = 0.2;

  int n_species = 1;

  int n_p = n_g * n_ppc;

  // Initialize fields
  Field e_x(n_g, "e_x");
  Field e_y(n_g, "e_y");
  Field b_z(n_g, "b_z");
  Field j_y(n_g, "j_y"); 
  Field j_x(n_g, "j_x"); 
  Field rho(n_g, "rho");

  // Define species
  SpeciesGroup particles(n_species, method, n_p, dt, dx, n_g);

  initialize_fields(e_x.field, e_y.field, b_z.field, j_x.field, j_y.field, n_g);
  particles.initialize_species(n_g, n_ppc, dx);

  // Calculate e_x at t = 0 from Poisson's equation
  particles.deposit_rho(rho.field, n_g);
  initialize_e_x(rho.field, e_x.field, dx, n_g);

  // Calculate b_z, u_x, and u_y at t = - 1/2
  advance_b_z(b_z.field, e_y.field, (-0.5 * dt), dx, n_g);

  ////////////////////////////////////////
  // MAIN LOOP
  for (int t = 0; t < n_t; t++) {
    // Print diagnostics for e_x, e_y, and x
    e_x.write_field();
    e_y.write_field();

    particles.deposit_rho(rho.field, n_g);

    write_data(rho.field, rho.output_stream, n_g);

    // Save b_z old
    save_old_values(b_z.field, b_z.field_old, n_g);

    // Calculate b_z at t+1/2
    advance_b_z(b_z.field, e_y.field, dt, dx, n_g);

    // Calculate b_z_tavg
    calc_tavg_arry(b_z.field_tavg, b_z.field, b_z.field_old, n_g);

    // Print diagnostics for b_z
    b_z.write_field_tavg();

    particles.save_u_x_old();
    particles.save_u_y_old();

    // Calculated e_x at integer values to eliminate self force
    half_int_to_int(e_x.field, e_x.field_int, n_g);

    // Calculate u_x and u_y, at t+1/2
    particles.advance_velocity(e_x.field_int, e_y.field, b_z.field_tavg);

    // Print diagnostics for u_x and u_y

    // Save j_x_old and j_y_old
    save_old_values(j_x.field, j_x.field_old, n_g);
    save_old_values(j_y.field, j_y.field_old, n_g);

    particles.save_x_old();

    // Calculate x at t+1
    for (int i = 0; i < n_species; i++) {
      particles.species[i].advance_x();
    }

    // Calculate j_x and j_y at t+1/2
    particles.deposit_j_x(j_x.field);
    particles.deposit_j_y(j_y.field);
    write_data_tavg(j_x.field, j_x.field_old, j_x.output_stream, n_g);
    write_data_tavg(j_y.field, j_y.field_old, j_y.output_stream, n_g);

    // Calculate e_x and e_y at t+1
    advance_e_x(e_x.field, j_x.field, dt, dx, n_g);
    advance_e_y(e_y.field, b_z.field, j_y.field, dt, dx, n_g);

    std::cout << t << std::endl;
  }
  
  return 0;

}   
// End of main

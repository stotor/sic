/**
 * sic
 * - The greatest PIC code of all time
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>

#include "mpi.h"

#include "speciesgroup.hpp"
#include "fields.hpp"
#include "utilities.hpp"

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  int num_procs, my_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // Load input parameters for initialization
  if (argc != 4) {
    std::cout << "Usage: ./sic <method> <n_ppc> <simulation_type>" << std::endl;
    return 1;
  }

  bool center_fields = true;
  int interp_order = 1;

  std::stringstream ss;

  ss << argv[1];
  int method;
  ss >> method;

  ss.str(std::string());
  ss.clear();
  long long n_ppc;
  ss << argv[2];
  ss >> n_ppc;

  ss.str(std::string());
  ss.clear();
  int simulation_type;
  ss << argv[3];
  ss >> simulation_type;
  
  int n_species;
  std::vector<double> u_x_drift, u_y_drift;
  double u_x_1, u_y_1, e_y_1, b_z_1;
  int mode, mode_max;
  double amplitude;
  int n_t, n_g;
  double k, dx, dt;

  double refinement_length = 0.1;

  switch (simulation_type) {
  case 0:
    // Electrostatic wave
    mode = 1;
    n_t = 500;
    n_g = 128;
    dx = 0.1;
    dt = 0.09;
    n_species = 1;
    u_x_drift.push_back(0.0);
    u_y_drift.push_back(0.0);
    u_x_1 = 0.0025;
    u_y_1 = 0.0;
    e_y_1 = 0.0;
    b_z_1 = 0.0;
    break;
  case 1:
    // Electromagnetic wave
    mode = 1;
    n_t = 500;
    n_g = 128;
    dx = 0.1;
    dt = 0.09;
    k = 2.0 * PI * mode / (n_g * dx);
    n_species = 1;
    u_x_drift.push_back(0.0);
    u_y_drift.push_back(0.0);
    u_x_1 = 0.0025;
    u_y_1 = 0.0025;
    e_y_1 = u_y_1 * sqrt(1 + k*k);
    b_z_1 = u_y_1 * k;
    break;
  case 2:
    // Two-stream instability
    mode = 10;
    n_t = 500;
    n_g = 128;
    dx = 0.1;
    dt = 0.09;
    n_species = 2;
    u_x_drift.push_back(0.176425);
    u_x_drift.push_back(-0.176425);
    u_y_drift.push_back(0.0);
    u_y_drift.push_back(0.0);
    u_x_1 = 0.00176425;
    u_y_1 = 0.0;
    e_y_1 = 0.0;
    b_z_1 = 0.0;
    break;
  case 3:
    // Weibel
    n_t = 500;
    n_g = 128;
    dx = 0.1;
    dt = 0.09;
    n_species = 2;
    u_x_drift.push_back(0.1);
    u_x_drift.push_back(-0.1);
    u_y_drift.push_back(0.5);
    u_y_drift.push_back(-0.5);
    u_x_1 = 0.0;
    u_y_1 = 0.0;
    mode_max = 64;
    amplitude = pow(10.0, -6.0);
    break; 
  }
  
  if (fmod((n_ppc * double(n_g)), 1.0) != 0.0) {
    std::cout << "Error: n_ppc * n_g is not an integer.";
    return 1;
  }
  
  SpeciesGroup particles(n_species, method, n_ppc, dt, dx, n_g, center_fields, interp_order, num_procs);
  particles.initialize_species(n_ppc, u_x_drift, u_y_drift, mode, u_x_1, u_y_1, my_rank, num_procs);

  Field e_x(n_g, "e_x", my_rank);
  Field e_y(n_g, "e_y", my_rank);
  Field b_z(n_g, "b_z", my_rank);
  Field j_x(n_g, "j_x", my_rank);
  Field j_y(n_g, "j_y", my_rank);
  Field rho(n_g, "rho", my_rank);
  
  particles.deposit_rho(rho.field, n_g, my_rank);
  sum_array_to_root(&rho.field[0], n_g, MPI_COMM_WORLD, my_rank);

  if (my_rank==0) {
    initialize_e_x(rho.field, e_x.field, dx, n_g);
    if (simulation_type == 3) {
      initialize_fields_weibel(e_y.field, b_z.field, n_g, dx, mode_max, amplitude);
    } else {
      initialize_transverse_em_fields(e_y.field, b_z.field, n_g, dx, e_y_1, b_z_1, mode);
    }
  }

  MPI_Bcast(&e_x.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&e_y.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&b_z.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  particles.initial_velocity_deceleration(e_x.field, e_y.field, b_z.field);

  for (int t = 0; t < n_t; t++) {
    particles.deposit_rho(rho.field, n_g, my_rank);
    sum_array_to_root(&rho.field[0], n_g, MPI_COMM_WORLD, my_rank);
    
    if (my_rank==0) {
      e_x.write_field();
      e_y.write_field();
      b_z.write_field();
      rho.write_field();
      e_x.calculate_energy();
      e_y.calculate_energy();
      b_z.calculate_energy();
    }

    particles.advance_velocity(e_x.field, e_y.field, b_z.field);

    particles.save_x_old();
    particles.advance_x();

    particles.deposit_j_x(j_x.field);
    particles.deposit_j_y(j_y.field);

    sum_array_to_root(&j_x.field[0], n_g, MPI_COMM_WORLD, my_rank);
    sum_array_to_root(&j_y.field[0], n_g, MPI_COMM_WORLD, my_rank);
    
    if (my_rank==0) {
      j_x.write_field();
      j_y.write_field();
      advance_e_x(e_x.field, j_x.field, dt, dx, n_g);
      advance_b_z(b_z.field, e_y.field, dt/2.0, dx, n_g);
      advance_e_y(e_y.field, b_z.field, j_y.field, dt, dx, n_g);
      advance_b_z(b_z.field, e_y.field, dt/2.0, dx, n_g);
    }

    MPI_Bcast(&e_x.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&e_y.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b_z.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank==0) {
      std::cout << t << std::endl;
      std::cout << particles.species[0].n_p << std::endl;
    }
    if (refinement_length) {
      particles.refine_segments(refinement_length);
    }
  }
  
  particles.write_energy_history(n_t, my_rank, MPI_COMM_WORLD);
  
  if (my_rank==0) {
    e_x.write_energy_history();
    e_y.write_energy_history();
    b_z.write_energy_history();
  }

  MPI_Finalize();

  return 0;
}

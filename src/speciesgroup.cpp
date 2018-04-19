#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "mpi.h"

#include "speciesgroup.hpp"
#include "particlespecies.hpp"
#include "utilities.hpp"

void SpeciesGroup::write_energy_history(int n_t, int my_rank, MPI_Comm COMM)
{
  for (int i = 0; i < n_species; i++) {
    species[i].write_energy_history(n_t, my_rank, COMM);
  }
  return;
}

void SpeciesGroup::initial_velocity_deceleration(std::vector<double> &e_x_int, 
						 std::vector<double> &e_y,
						 std::vector<double> &b_z)
{
  for (int i = 0; i < n_species; i++) {
    species[i].initial_velocity_deceleration(e_x_int, e_y, b_z);
  }
  return;
}

void SpeciesGroup::deposit_j_x(std::vector<double> &j_x)
{
  for (int i = 0; i < n_g; i++) {
    j_x[i] = 0.0;
  }
 
  for (int i = 0; i < n_species; i++) {
    if (method==0) { 
      species[i].deposit_j_x_ngp(j_x);
    }
    if (method==1) { 
      species[i].deposit_j_x(j_x);
    }
    else if (method==2) {
      species[i].deposit_j_x_segments_zero(j_x);
    }
    else if (method==3) {
      species[i].deposit_j_x_segments_linear(j_x);
    }
  }
  return;
}

void SpeciesGroup::deposit_j_y(std::vector<double> &j_y)
{
  for (int i = 0; i < n_g; i++) {
    j_y[i] = 0.0;
  }
  for (int i = 0; i < n_species; i++) {
    if (method==0) { 
      species[i].deposit_j_y_ngp(j_y);
    }
    if (method==1) { 
      species[i].deposit_j_y(j_y);
    }
    else if (method==2) {
      species[i].deposit_j_y_segments_zero(j_y);
    }
    else if (method==3) {
      species[i].deposit_j_y_segments_linear(j_y);
    }
  }
  return;
}

void SpeciesGroup::save_x_old()
{
  for (int i = 0; i < n_species; i++) {
    save_old_values(species[i].x, species[i].x_old, species[i].n_p);
  }
  return;
}

void SpeciesGroup::advance_velocity(std::vector<double> e_x_int, 
				    std::vector<double> e_y,
				    std::vector<double> b_z)
{
  for (int i = 0; i < n_species; i++) {
    species[i].advance_velocity(e_x_int, e_y, b_z);
  }
  return;
}

void SpeciesGroup::initialize_species(long long n_ppc, 
				      std::vector<double> u_x_drift, 
				      std::vector<double> u_y_drift, 
				      int mode, 
				      double u_x_1, 
				      double u_y_1,
				      int my_rank,
				      int num_procs)
{
  for (int i = 0; i < n_species; i++) {
    species[i].initialize_species(i, n_ppc, u_x_drift[i], u_y_drift[i], mode, u_x_1, u_y_1, my_rank, num_procs);
  }
  return;
}

void SpeciesGroup::deposit_rho(std::vector<double> &rho, int n_g, int my_rank)
{
  // First add neutralizing background density
  if (my_rank==0) { 
    for (int i = 0; i < n_g; i++) {
      rho[i] = 1.0 * n_species;
    }
  } else { 
    for (int i = 0; i < n_g; i++) {
      rho[i] = 0.0;
    }
  }
  
  for (int i = 0; i < n_species; i++) {
    if (method==0) { 
      species[i].deposit_rho_ngp(rho);
    }
    if (method==1) { 
      species[i].deposit_rho(rho);
    }
    else if (method==2) { 
      species[i].deposit_rho_segments_zero(rho);
    }
    else if (method==3) { 
      species[i].deposit_rho_segments_linear(rho);
    }
  }
  return;
}

void SpeciesGroup::advance_x()
{
  for (int i = 0; i < n_species; i++) {
    species[i].advance_x();
  }
  return;
}

void SpeciesGroup::refine_segments(double refinement_length)
{
  for (int i = 0; i < n_species; i++) {
    species[i].refine_segments(refinement_length);
  }
  return;
}

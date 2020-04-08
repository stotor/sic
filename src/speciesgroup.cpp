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

void SpeciesGroup::write_particle_diagnostics(int n_t, int my_rank, MPI_Comm COMM)
{
  for (int i = 0; i < n_species; i++) {
    species[i].write_particle_diagnostics(n_t, my_rank, COMM);
  }
  return;
}

void SpeciesGroup::initial_velocity_deceleration(std::vector<double> &e_x, 
						 std::vector<double> &e_y,
						 std::vector<double> &e_z,
						 std::vector<double> &b_x,
						 std::vector<double> &b_y,
						 std::vector<double> &b_z)
{
  for (int i = 0; i < n_species; i++) {
    species[i].initial_velocity_deceleration(e_x, e_y, e_z, b_x, b_y, b_z);
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
      species[i].deposit_j_x_pic_0(j_x);
    } else if (method==1) { 
      species[i].deposit_j_x_pic_1(j_x);
    }
    else if (method==2) {
      species[i].deposit_j_x_sic_0(j_x);
    }
    else if (method==3) {
      species[i].deposit_j_x_sic_1(j_x);
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
      species[i].deposit_j_y_pic_0(j_y);
    } else if (method==1) { 
      species[i].deposit_j_y_pic_1(j_y);
    }
    else if (method==2) {
      species[i].deposit_j_y_sic_0(j_y);
    }
    else if (method==3) {
      species[i].deposit_j_y_sic_1(j_y);
    }
  }
  return;
}

void SpeciesGroup::deposit_j_z(std::vector<double> &j_z)
{
  for (int i = 0; i < n_g; i++) {
    j_z[i] = 0.0;
  }
  for (int i = 0; i < n_species; i++) {
    if (method==0) { 
      species[i].deposit_j_z_pic_0(j_z);
    } else if (method==1) { 
      species[i].deposit_j_z_pic_1(j_z);
    }
    else if (method==2) {
      species[i].deposit_j_z_sic_0(j_z);
    }
    else if (method==3) {
      species[i].deposit_j_z_sic_1(j_z);
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

void SpeciesGroup::advance_velocity(std::vector<double> &e_x,
				    std::vector<double> &e_y,
				    std::vector<double> &e_z,
				    std::vector<double> &b_x,
				    std::vector<double> &b_y,
				    std::vector<double> &b_z)
{
  for (int i = 0; i < n_species; i++) {
    species[i].advance_velocity(e_x, e_y, e_z, b_x, b_y, b_z);
  }
  return;
}

void SpeciesGroup::initialize_species(long long n_ppc, 
				      int my_rank,
				      int num_procs)
{
  for (int i = 0; i < n_species; i++) {
    species[i].initialize_species(i, n_ppc, my_rank, num_procs);
  }
  return;
}

void SpeciesGroup::deposit_rho(std::vector<double> &rho)
{
  for (int i = 0; i < n_g; i++) {
    rho[i] = 0.0;
  }
  
  for (int i = 0; i < n_species; i++) {
    if (method==0) { 
      species[i].deposit_rho_pic_0(rho);
    } else if (method==1) { 
      species[i].deposit_rho_pic_1(rho);
    }
    else if (method==2) { 
      species[i].deposit_rho_sic_0(rho);
    }
    else if (method==3) { 
      species[i].deposit_rho_sic_1(rho);
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


void SpeciesGroup::communicate_ghost_particles(MPI_Comm COMM)
{
  for (int i = 0; i < n_species; i++) {
    species[i].communicate_ghost_particles(COMM);
  }
  return;
}

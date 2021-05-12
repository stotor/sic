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

void SpeciesGroup::deposit_j_x(std::vector<double> &j_x, int my_rank, MPI_Comm COMM)
{
  for (int i = 0; i < n_g; i++) {
    j_x[i] = 0.0;
  }
  for (int i = 0; i < n_species; i++) {
    switch (species[i].method) {
    case 0 :
      species[i].deposit_j_x_pic_0(j_x);
      break;
    case 1 :
      species[i].deposit_j_x_pic_1(j_x);
      break;
    case 2 :
      species[i].deposit_j_x_pic_2(j_x);      
      break;
    case 3 :
      species[i].deposit_j_x_pic_3(j_x);      
      break;
    case 4 :
      species[i].deposit_j_x_pic_4(j_x);
      break;
    case 5 :
      species[i].deposit_j_x_sic_0(j_x);
      break;
    case 6 :
      species[i].deposit_j_x_sic_1(j_x);
      break;
    case 7 :
      species[i].deposit_j_x_sic_2(j_x);
      break;
    case 8 :
      species[i].deposit_j_x_sic_3(j_x);
      break;      
    case 9 :
      species[i].deposit_j_x_sic_4(j_x);
      break;      
    case 10 :
      species[i].deposit_j_x_sic_higher_order_0(j_x);
      break;
    case 11 :
      species[i].deposit_j_x_sic_higher_order_1(j_x);
      break;
    default:
    std::cout << "Error, selected interpolation order not implemented." << std::endl;
    }
  }
    
  sum_array_to_root(&j_x[0], n_g, COMM, my_rank);
  return;
}

void SpeciesGroup::deposit_j_y(std::vector<double> &j_y, int my_rank, MPI_Comm COMM)
{
  for (int i = 0; i < n_g; i++) {
    j_y[i] = 0.0;
  }
  for (int i = 0; i < n_species; i++) {
    switch (species[i].method) {
    case 0 :
      species[i].deposit_j_y_pic_0(j_y);
      break;
    case 1 :
      species[i].deposit_j_y_pic_1(j_y);
      break;
    case 2 :
      species[i].deposit_j_y_pic_2(j_y);
      break;
    case 3 :
      species[i].deposit_j_y_pic_3(j_y);
      break;
    case 4 :
      species[i].deposit_j_y_pic_4(j_y);
      break;
    case 5 :
      species[i].deposit_j_y_sic_0(j_y);
      break;
    case 6 :
      species[i].deposit_j_y_sic_1(j_y);
      break;
    case 7 :
      species[i].deposit_j_y_sic_2(j_y);
      break;
    case 8 :
      species[i].deposit_j_y_sic_3(j_y);
      break;
    case 9 :
      species[i].deposit_j_y_sic_4(j_y);
    case 10 :
      species[i].deposit_j_y_sic_higher_order_0(j_y);
      break;
    case 11 :
      species[i].deposit_j_y_sic_higher_order_1(j_y);
      break;
    default:
    std::cout << "Error, selected interpolation order not implemented." << std::endl;
    }
  }
  sum_array_to_root(&j_y[0], n_g, COMM, my_rank);
  return;
}

void SpeciesGroup::deposit_j_z(std::vector<double> &j_z, int my_rank, MPI_Comm COMM)
{
  for (int i = 0; i < n_g; i++) {
    j_z[i] = 0.0;
  }
  for (int i = 0; i < n_species; i++) {
    switch (species[i].method) {
    case 0 :
      species[i].deposit_j_z_pic_0(j_z);
      break;
    case 1 :
      species[i].deposit_j_z_pic_1(j_z);
      break;
    case 2 :
      species[i].deposit_j_z_pic_2(j_z);
      break;
    case 3 :
      species[i].deposit_j_z_pic_3(j_z);
      break;
    case 4 :
      species[i].deposit_j_z_pic_4(j_z);
      break;
    case 5 :
      species[i].deposit_j_z_sic_0(j_z);
      break;
    case 6 :
      species[i].deposit_j_z_sic_1(j_z);
      break;
    case 7 :
      species[i].deposit_j_z_sic_2(j_z);
      break;
    case 8 :
      species[i].deposit_j_z_sic_3(j_z);
      break;
    case 9 :
      species[i].deposit_j_z_sic_4(j_z);
      break;
    case 10 :
      species[i].deposit_j_z_sic_higher_order_0(j_z);
      break;
    case 11 :
      species[i].deposit_j_z_sic_higher_order_1(j_z);
      break;
    default:
    std::cout << "Error, selected interpolation order not implemented." << std::endl;
    }
  }
  sum_array_to_root(&j_z[0], n_g, COMM, my_rank);
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

void SpeciesGroup::initialize_species(std::vector<double> n_ppc,
				      int my_rank,
				      int num_procs,
				      std::vector<int> method,
				      int simulation_type)
{
  for (int i = 0; i < n_species; i++) {
    species[i].initialize_species(i, n_ppc[i], my_rank, num_procs, method[i], simulation_type);
  }
  return;
}

void SpeciesGroup::deposit_rho(std::vector<double> &rho, double rho_bg, int my_rank, MPI_Comm COMM)
{
  if (my_rank==0) {
    for (int i = 0; i < n_g; i++) {
      rho[i] = rho_bg;
    }
  } else {
    for (int i = 0; i < n_g; i++) {
      rho[i] = 0.0;
    }
  }  
  for (int i = 0; i < n_species; i++) {
    switch (species[i].method) {
    case 0 :
      species[i].deposit_rho_pic_0(rho);
      break;
    case 1 :
      species[i].deposit_rho_pic_1(rho);
      break;
    case 2 :
      species[i].deposit_rho_pic_2(rho);
      break;
    case 3 :
      species[i].deposit_rho_pic_3(rho);
      break;
    case 4 :
      species[i].deposit_rho_pic_4(rho);
      break;
    case 5 :
      species[i].deposit_rho_sic_0(rho);
      break;
    case 6 :
      species[i].deposit_rho_sic_1(rho);
      break;
    case 7 :
      species[i].deposit_rho_sic_2(rho);
      break;
    case 8 :
      species[i].deposit_rho_sic_3(rho);
      break;
    case 9 :
      species[i].deposit_rho_sic_4(rho);
      break;
    case 10 :
      species[i].deposit_rho_sic_higher_order_0(rho);
      break;
    case 11 :
      species[i].deposit_rho_sic_higher_order_1(rho);
      break;
    default:
    std::cout << "Error, selected interpolation order not implemented." << std::endl;
    }
  }

  sum_array_to_root(&rho[0], n_g, COMM, my_rank);
  return;
}

void SpeciesGroup::write_phase(int t, int my_rank)
{
  for (int i = 0; i < n_species; i++) {
    species[i].write_phase(i, t, my_rank);
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
    if (species[i].method > 4) {
      species[i].communicate_ghost_particles(COMM);
    }
  }
  return;
}

void SpeciesGroup::calculate_segment_density(MPI_Comm COMM)
{
  for (int i = 0; i < n_species; i++) {
    if (species[i].method > 4) {
      species[i].calculate_segment_density(COMM);
    }
  }
  return;
}

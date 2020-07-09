#ifndef speciesgroup_hpp
#define speciesgroup_hpp

#include <vector>
#include <fstream>

#include "mpi.h"

#include "particlespecies.hpp"

class SpeciesGroup {
public:
  SpeciesGroup(int n_species, double dt, double dx, int n_g, bool center_fields, int interp_order) :
    species(n_species, ParticleSpecies(dt, dx, n_g, center_fields, interp_order))
  {
    this->n_species = n_species;
    this->dt = dt;
    this->dx = dx;
    this->n_g = n_g;
  }
  ~SpeciesGroup() {}

  int n_species, n_g;
  double dt, dx;

  // Attributes
  std::vector<ParticleSpecies> species;

  // Methods
  void deposit_rho(std::vector<double> &rho, double rho_bg, int my_rank, MPI_Comm COMM);
  void deposit_j_x(std::vector<double> &j_x, int my_rank, MPI_Comm COMM);
  void deposit_j_y(std::vector<double> &j_y, int my_rank, MPI_Comm COMM);
  void deposit_j_z(std::vector<double> &j_z, int my_rank, MPI_Comm COMM);
  
  void save_x_old();
  void advance_x();
  void refine_segments(double refinement_length);
  void advance_velocity(std::vector<double> &e_x,
			std::vector<double> &e_y,
			std::vector<double> &e_z,
			std::vector<double> &b_x,
			std::vector<double> &b_y,			
			std::vector<double> &b_z);
  void initial_velocity_deceleration(std::vector<double> &e_x,
				     std::vector<double> &e_y,
				     std::vector<double> &e_z,
				     std::vector<double> &b_x,
				     std::vector<double> &b_y,
				     std::vector<double> &b_z);
  void write_particle_diagnostics(int n_t, int my_rank, MPI_Comm COMM);
  void write_phase(int t, int my_rank);  
  void initialize_species(std::vector<double> n_ppc, 
			  int my_rank,
			  int num_procs,
			  std::vector<int> method,
			  int simulation_type);
  void communicate_ghost_particles(MPI_Comm COMM);
  void calculate_segment_density(MPI_Comm COMM);  
};

#endif /* speciesgroup_hpp */

#ifndef speciesgroup_hpp
#define speciesgroup_hpp

#include <vector>
#include <fstream>

#include "mpi.h"

#include "particlespecies.hpp"

class SpeciesGroup {
public:
  SpeciesGroup(int n_species, int method, long long n_ppc, double dt, double dx, int n_g, bool center_fields, int interp_order, int num_procs) :
    species(n_species, ParticleSpecies(method, n_ppc, dt, dx, n_g, center_fields, interp_order, num_procs))
  {
    this->n_species = n_species;
    this->method = method;
    this->dt = dt;
    this->dx = dx;
    this->n_g = n_g;
  }
  ~SpeciesGroup() {}

  int n_species, method, n_g;
  double dt, dx;

  // Attributes
  std::vector<ParticleSpecies> species;

  // Methods
  void deposit_rho(std::vector<double> &rho, int n_g, int my_rank, MPI_Comm COMM);
  void save_x_old();
  void advance_x();
  void refine_segments(double refinement_length);
  void advance_velocity(std::vector<double> e_x_int, std::vector<double> e_y, 
			std::vector<double> b_z);
  void deposit_j_x(std::vector<double> &j_x, MPI_Comm COMM);
  void deposit_j_y(std::vector<double> &j_y, MPI_Comm COMM);
  void initial_velocity_deceleration(std::vector<double> &e_x_int, 
				     std::vector<double> &e_y,
				     std::vector<double> &b_z);
  void write_particle_diagnostics(int n_t, int my_rank, MPI_Comm COMM);
  void initialize_species(long long n_ppc, 
			  std::vector<double> u_x_drift, 
			  std::vector<double> u_y_drift, 
			  int mode, 
			  double u_x_1, 
			  double u_y_1,
			  int my_rank,
			  int num_procs);
  void communicate_ghost_particles(MPI_Comm COMM);
  void u_x_perturbation(double amplitude, int mode_max);
  void initialize_beat_heating(int mode_1, int mode_2,
			       double phase_1, double phase_2,
			       double vel_amp);
};

#endif /* speciesgroup_hpp */

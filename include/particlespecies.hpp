#ifndef particlespecies_hpp
#define particlespecies_hpp

#include <vector>
#include <fstream>
#include <string>

#include "mpi.h"

class ParticleSpecies {
public:
  ParticleSpecies(int method, long long n_ppc, double dt, double dx, int n_g, bool center_fields, int interp_order, int num_procs)
  {
    this->method = method;
    this->n_p = (n_ppc * n_g) / num_procs;
    this->n_ppp = n_p;
    this->dt = dt;
    this->dx = dx;
    this->n_g = n_g;
    this->center_fields = center_fields;
    this->interp_order = interp_order;
    x.resize(n_p);
    u_x.resize(n_p);
    u_y.resize(n_p);
    x_old.resize(n_p);
    charge.resize(n_p);
    rqm = -1.0;
    lagrangian_id.resize(n_p);
    return;
  }

  ~ParticleSpecies()
  {
  }

  // Simulation box parameters
  double dt, dx;
  int n_g;

  // Attributes
  long long n_p, n_ppp;
  int method, interp_order;
  bool center_fields;

  std::vector<double> x, u_x, u_y, x_old, lagrangian_id,
    energy_history, momentum_x_history, momentum_y_history, n_p_history;

  // Charge to mass ratio, and particle charge divided by grid spacing
  double rqm;
  std::vector<double> charge;

  std::string species_name;

  // Methods
  void initialize_species(int species_number, 
			  long long n_ppc,
			  double u_x_drift, 
			  double u_y_drift, 
			  int mode, 
			  double u_x_1, 
			  double u_y_1,
			  int my_rank,
			  int num_procs);
  void advance_x();
  void split_segment_linear(int i);
  void split_segment_lagrange_3(int i);
  void refine_segments(double refinement_length);
  void advance_velocity(std::vector<double> &e_x, std::vector<double> &e_y, 
			std::vector<double> &b_z);
  void initial_velocity_deceleration(std::vector<double> &e_x, 
				     std::vector<double> &e_y,
				     std::vector<double> &b_z);
  void deposit_rho(std::vector<double> &rho);
  void deposit_rho_ngp(std::vector<double> &rho);
  void deposit_rho_segments_zero(std::vector<double> &rho);
  void deposit_rho_segments_linear(std::vector<double> &rho);
  void deposit_j_x(std::vector<double> &j_x);
  void deposit_j_x_ngp(std::vector<double> &j_x);
  void deposit_j_x_segments_zero(std::vector<double> &j_x);
  void deposit_j_x_segments_linear(std::vector<double> &j_x);
  void deposit_j_y(std::vector<double> &j_y);
  void deposit_j_y_ngp(std::vector<double> &j_y);
  void deposit_j_y_segments_zero(std::vector<double> &j_y);
  void deposit_j_y_segments_linear(std::vector<double> &j_y);
  void write_phase(std::ofstream &x_ofstream, std::ofstream &u_x_ofstream, 
		   std::ofstream &u_y_ofstream);
  void write_particle_diagnostics(int n_t, int my_rank, MPI_Comm COMM);
  void communicate_ghost_particles(MPI_Comm COMM);
};

#endif /* particlespecies_hpp */

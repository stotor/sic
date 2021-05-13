#ifndef particlespecies_hpp
#define particlespecies_hpp

#include <vector>
#include <fstream>
#include <string>

#include "mpi.h"

class ParticleSpecies {
public:
  ParticleSpecies(double dt, double dx, int n_g, bool center_fields)
  {
    this->dt = dt;
    this->dx = dx;
    this->n_g = n_g;
    this->center_fields = center_fields;
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

  std::vector<double> x, u_x, u_y, u_z, x_old, lagrangian_id,
    energy_history, momentum_x_history, momentum_y_history, momentum_z_history,
    n_p_history, density, density_old, density_tavg;

  double rqm; // Reciprocal charge to mass ratio
  std::vector<double> charge; // Particle charge divided by grid spacing

  std::string species_name;

  // Methods
  void initialize_species(int species_number,
			  double n_ppc,
			  int my_rank,
			  int num_procs,
			  int method,
			  int simulation_type);
  void advance_x();
  void split_segment_linear(int i);
  void split_segment_lagrange_3(int i);
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
  void deposit_rho_pic_0(std::vector<double> &rho);
  void deposit_rho_pic_1(std::vector<double> &rho);
  void deposit_rho_pic_2(std::vector<double> &rho);
  void deposit_rho_pic_3(std::vector<double> &rho);
  void deposit_rho_pic_4(std::vector<double> &rho);
  void deposit_rho_sic_0(std::vector<double> &rho);
  void deposit_rho_sic_1(std::vector<double> &rho);
  void deposit_rho_sic_2(std::vector<double> &rho);
  void deposit_rho_sic_3(std::vector<double> &rho);
  void deposit_rho_sic_4(std::vector<double> &rho);
  void deposit_rho_sic_higher_order_0(std::vector<double> &rho);
  void deposit_rho_sic_higher_order_1(std::vector<double> &rho);
  void deposit_rho_sic_center(std::vector<double> &rho);
  void deposit_rho_sic_higher_order_center(std::vector<double> &rho);
  
  void deposit_j_x_pic_0(std::vector<double> &j_x);
  void deposit_j_x_pic_1(std::vector<double> &j_x);
  void deposit_j_x_pic_2(std::vector<double> &j_x);
  void deposit_j_x_pic_3(std::vector<double> &j_x);
  void deposit_j_x_pic_4(std::vector<double> &j_x);
  void deposit_j_x_sic_0(std::vector<double> &j_x);
  void deposit_j_x_sic_1(std::vector<double> &j_x);
  void deposit_j_x_sic_2(std::vector<double> &j_x);
  void deposit_j_x_sic_3(std::vector<double> &j_x);
  void deposit_j_x_sic_4(std::vector<double> &j_x);
  void deposit_j_x_sic_higher_order_0(std::vector<double> &j_x);
  void deposit_j_x_sic_higher_order_1(std::vector<double> &j_x);
  void deposit_j_x_sic_center(std::vector<double> &j_x);
  void deposit_j_x_sic_higher_order_center(std::vector<double> &j_x);

  void deposit_j_y_pic_0(std::vector<double> &j_y);
  void deposit_j_y_pic_1(std::vector<double> &j_y);
  void deposit_j_y_pic_2(std::vector<double> &j_y);
  void deposit_j_y_pic_3(std::vector<double> &j_y);
  void deposit_j_y_pic_4(std::vector<double> &j_y);
  void deposit_j_y_sic_0(std::vector<double> &j_y);
  void deposit_j_y_sic_1(std::vector<double> &j_y);
  void deposit_j_y_sic_2(std::vector<double> &j_y);
  void deposit_j_y_sic_3(std::vector<double> &j_y);
  void deposit_j_y_sic_4(std::vector<double> &j_y);
  void deposit_j_y_sic_higher_order_0(std::vector<double> &j_y);
  void deposit_j_y_sic_higher_order_1(std::vector<double> &j_y);
  void deposit_j_y_sic_center(std::vector<double> &j_y);
  void deposit_j_y_sic_higher_order_center(std::vector<double> &j_y);

  void deposit_j_z_pic_0(std::vector<double> &j_z);
  void deposit_j_z_pic_1(std::vector<double> &j_z);
  void deposit_j_z_pic_2(std::vector<double> &j_z);
  void deposit_j_z_pic_3(std::vector<double> &j_z);
  void deposit_j_z_pic_4(std::vector<double> &j_z);
  void deposit_j_z_sic_0(std::vector<double> &j_z);
  void deposit_j_z_sic_1(std::vector<double> &j_z);
  void deposit_j_z_sic_2(std::vector<double> &j_z);
  void deposit_j_z_sic_3(std::vector<double> &j_z);
  void deposit_j_z_sic_4(std::vector<double> &j_z);
  void deposit_j_z_sic_higher_order_0(std::vector<double> &j_z);
  void deposit_j_z_sic_higher_order_1(std::vector<double> &j_z);
  void deposit_j_z_sic_center(std::vector<double> &j_z);
  void deposit_j_z_sic_higher_order_center(std::vector<double> &j_z);

  void write_phase(int species_number, int t, int my_rank);
  void write_particle_diagnostics(int n_t, int my_rank, MPI_Comm COMM);
  void communicate_ghost_particles(MPI_Comm COMM);
  void calculate_segment_density(MPI_Comm COMM);
};

#endif /* particlespecies_hpp */

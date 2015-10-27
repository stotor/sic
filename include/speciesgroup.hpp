#ifndef speciesgroup_hpp
#define speciesgroup_hpp

#include <vector>
#include <fstream>

#include "particlespecies.hpp"

class SpeciesGroup {
public:
  SpeciesGroup(int n_species, int method, int n_ppc, double dt, double dx, int n_g) :
    species(n_species, ParticleSpecies(method, n_ppc, dt, dx, n_g))
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
  void deposit_rho(std::vector<double> &rho, int n_g);
  void save_x_old();
  void advance_x();
  void save_u_x_old();
  void save_u_y_old();
  void advance_velocity(std::vector<double> e_x_int, std::vector<double> e_y, 
			std::vector<double> b_z_tavg);
  void deposit_j_x(std::vector<double> &j_x);
  void deposit_j_y(std::vector<double> &j_y);
  void initial_velocity_deceleration(std::vector<double> &e_x_int, 
				     std::vector<double> &e_y,
				     std::vector<double> &b_z_tavg);
  void write_energy_history();
  void initialize_species(double n_ppc, 
			  std::vector<double> u_x_drift, 
			  std::vector<double> u_y_drift, 
			  int mode, 
			  double u_x_1, 
			  double u_y_1);
};

#endif /* speciesgroup_hpp */

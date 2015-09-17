#ifndef particles_hpp
#define particles_hpp

#include <vector>
#include <fstream>


class ParticleSpecies {
public:
  ParticleSpecies(int n_p);
  ~ParticleSpecies();
  // Simulation box parameters
  double dt, dx;
  int n_g;

  // Attributes
  // Number of particles
  int n_p;

  // Particle phase space attributes
  std::vector<double> x, u_x, u_y, x_old, u_x_old, u_y_old;

  // Charge to mass ratio, and particle charge divided by grid spacing
  std::vector<double> rqm, charge;

  // Methods
  void advance_x();
  void advance_velocity(std::vector<double> &e_x_int, std::vector<double> &e_y, 
			std::vector<double> &b_z_tavg);
  void deposit_rho(std::vector<double> &rho);
  void deposit_rho_segments_zero(std::vector<double> &rho);
  void deposit_rho_segments_linear(std::vector<double> &rho);
  void deposit_j_x(std::vector<double> &j_x);
  void deposit_j_x_segments_zero(std::vector<double> &j_x);
  void deposit_j_x_segments_linear(std::vector<double> &j_x);
  void deposit_j_y(std::vector<double> &j_y);
  void deposit_j_y_segments_zero(std::vector<double> &j_y);
  void write_phase(std::ofstream &x_ofstream, 
		   std::ofstream &u_x_ofstream, std::ofstream &u_y_ofstream);

};

#endif /* particles_hpp */

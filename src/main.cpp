/**
 * 1DEM
 * - 1D electromagnetic code
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "particles.hpp"
#include "fields.hpp"
#include "utilities.hpp"

int main(int argc, char *argv[])
{
  ///////////////////////////////
  // INITIALIZATION

  // Define simulation parameters
  int n_t = 501;
  int n_g = 200;

  double dx = 0.5;
  double dt = 0.2;

  int n_species = 1;
  double n_ppc = 1.0 / 10.0;

  int n_p = n_g * n_ppc;

  int method = 0;

  bool line_segments;
  if (method==0) {
    line_segments = false;
  } else {
    line_segments = true;
  }


  // Initialize
  std::vector<double> e_x(n_g);
  std::vector<double> e_x_int(n_g);
  std::vector<double> e_y(n_g);
  std::vector<double> b_z(n_g); 
  std::vector<double> j_y(n_g); 
  std::vector<double> j_x(n_g); 
  std::vector<double> rho(n_g); 
  std::vector<double> phi(n_g); 
  std::vector<double> b_z_old(n_g);
  std::vector<double> j_y_old(n_g);
  std::vector<double> j_x_old(n_g);
  std::vector<double> b_z_tavg(n_g);

  // Define species
  std::vector<ParticleSpecies> species;
  for (int i = 0; i < n_species; i++) {
    species.push_back(ParticleSpecies(n_p));
    species[i].dt = dt;
    species[i].dx = dx;
    species[i].n_g = n_g;
  }
  
  // Initialize output file streams
  std::ofstream e_x_ofstream, e_y_ofstream, b_z_ofstream, j_x_ofstream, 
    j_y_ofstream, rho_ofstream, phi_ofstream;

  std::ofstream x_ofstream, u_x_ofstream, u_y_ofstream;

  e_x_ofstream.open("e_x");
  e_y_ofstream.open("e_y");
  b_z_ofstream.open("b_z");
  j_x_ofstream.open("j_x");
  j_y_ofstream.open("j_y");
  rho_ofstream.open("rho");
  phi_ofstream.open("phi");

  x_ofstream.open("x");
  u_x_ofstream.open("u_x");
  u_y_ofstream.open("u_y");

  initialize_fields(e_x, e_y, b_z, j_x, j_y, n_g);
  // Particle initialization
  // Electrostatic wave, electromagnetic wave, Weibel

  // Initialize charge, x, u_x, and u_y at t = 0
  // Weibel initialization
  // double u_y_drift[2] = {-0.1, 0.1};
  // for (int i_species = 0; i_species < n_species; i_species++) {
  //   for (int i_particle = 0; i_particle < species[i_species].n_p; i_particle++) {
  //     species[i_species].charge[i_particle] = (-1.0) * double(n_g) / n_p;
  //     species[i_species].rqm[i_particle] = -1.0;
  //     species[i_species].u_x[i_particle] = 0.01 * ((double) rand() / (RAND_MAX));
  //     species[i_species].u_y[i_particle] = u_y_drift[i_species];
  //     species[i_species].x[i_particle] = (double(i_particle) / species[i_species].n_p) * n_g * dx;      
  //   }
  // }

  // Electrostatic wave initialization
  double wave_amplitude = 0.1;
  int wave_mode = 1;
  for (int i_species = 0; i_species < n_species; i_species++) {
    for (int i_particle = 0; i_particle < species[i_species].n_p; i_particle++) {
      species[i_species].relativistic = false;
      species[i_species].charge[i_particle] = (-1.0) * (1.0 / n_ppc);
      species[i_species].rqm[i_particle] = -1.0;
      species[i_species].x[i_particle] = (double(i_particle) / n_p) * n_g * dx + (double(n_g) * dx / double(n_p)) / 2.0;
      species[i_species].u_x[i_particle] = wave_amplitude * cos(2.0 * PI * double(wave_mode) * species[i_species].x[i_particle] / (n_g * dx));
      species[i_species].u_y[i_particle] = 0.0;
    }
  }

  if (line_segments) {
    for (int i = 0; i < n_species; i++) {
      species[i].charge.push_back(0.0);
      species[i].rqm.push_back(species[i].rqm[0]);
      species[i].u_x.push_back(species[i].u_x[0]);
      species[i].u_y.push_back(species[i].u_y[0]);
      species[i].x.push_back(species[i].x[0] + n_g*dx);
      species[i].n_p += 1;
    }
  }

  // Calculate e_x at t = 0 from Poisson's equation
  for (int i = 0; i < n_g; i++) {
    rho[i] = 1.0;
  }
  for (int i = 0; i < n_species; i++) {
    if (method==0) { 
      species[i].deposit_rho(rho);
    }
    else if (method==1) { 
      species[i].deposit_rho_segments_zero(rho);
    }
    else if (method==2) { 
      species[i].deposit_rho_segments_linear(rho);
    }
  }
  initialize_e_x(rho, e_x, dx, n_g);

  // Calculate b_z, u_x, and u_y at t = - 1/2
  advance_b_z(b_z, e_y, (-0.5 * dt), dx, n_g);

  ////////////////////////////////////////
  // MAIN LOOP
  int ndump_phase = 1;
  for (int t = 0; t < n_t; t++) {
    if (t % ndump_phase == 0) {
      species[0].write_phase(x_ofstream, u_x_ofstream, u_y_ofstream);
    }
    // Print diagnostics for e_x, e_y, and x
    write_data(e_x, e_x_ofstream, n_g);
    write_data(e_y, e_y_ofstream, n_g);

    // Deposit the charge from each species
    // First set rho to background charge density
    for (int i = 0; i < n_g; i++) {
      rho[i] = 1.0;
    }
    for (int i = 0; i < n_species; i++) {
      if (method==0) { 
	species[i].deposit_rho(rho);
      }
      else if (method==1) { 
	species[i].deposit_rho_segments_zero(rho);
      }
      else if (method==2) { 
	species[i].deposit_rho_segments_linear(rho);
      }
    }
    write_data(rho, rho_ofstream, n_g);

    // Save b_z_old
    save_old_values(b_z, b_z_old, n_g);

    // Calculate b_z at t+1/2
    advance_b_z(b_z, e_y, dt, dx, n_g);

    // Calculate b_z_tavg
    calc_tavg_arry(b_z_tavg, b_z, b_z_old, n_g);

    // Print diagnostics for b_z
    write_data(b_z_tavg, b_z_ofstream, n_g);

    // Save u_x_old and u_y_old
    for (int i = 0; i < n_species; i++) {
      save_old_values(species[i].u_x, species[i].u_x_old, species[i].n_p);
      save_old_values(species[i].u_y, species[i].u_y_old, species[i].n_p);
    }
    
    // Calculated e_x at integer values to eliminate self force
    half_int_to_int(e_x, e_x_int, n_g);

    // Calculate u_x and u_y, at t+1/2
    for (int i = 0; i < n_species; i++) {
      species[i].advance_velocity(e_x_int, e_y, b_z_tavg);
    }

    // Print diagnostics for u_x and u_y

    // Save j_x_old and j_y_old
    save_old_values(j_x, j_x_old, n_g);
    save_old_values(j_y, j_y_old, n_g);

    // Save x_old
    for (int i = 0; i < n_species; i++) {
      save_old_values(species[i].x, species[i].x_old, species[i].n_p);
    }

    // Calculate x at t+1
    for (int i = 0; i < n_species; i++) {
      species[i].advance_x();
    }

    // Calculate j_x and j_y at t+1/2
    // First set currents to zero
    for (int i = 0; i < n_g; i++) {
      j_x[i] = 0.0;
      j_y[i] = 0.0;
    }
    // Deposit the current from each species
    for (int i = 0; i < n_species; i++) {
      if (method==0) { 
	species[i].deposit_j_x(j_x);
      }
      else if (method==1) {
	species[i].deposit_j_x_segments_zero(j_x);
      }
      else if (method==2) {
	species[i].deposit_j_x_segments_linear(j_x);
      }
      species[i].deposit_j_y(j_y);
    }

    // Print diagnostics for j_x and j_y
    write_data_tavg(j_x, j_x_old, j_x_ofstream, n_g);
    write_data_tavg(j_y, j_y_old, j_y_ofstream, n_g);

    // Calculate e_x and e_y at t+1
    advance_e_x(e_x, j_x, dt, dx, n_g);
    advance_e_y(e_y, b_z, j_y, dt, dx, n_g);

    std::cout << t << std::endl;
  }
  
  ////////////////////////////////////
  // CLEAN UP

  // Close output file streams
  e_x_ofstream.close();
  e_y_ofstream.close();
  b_z_ofstream.close();
  j_y_ofstream.close();
  rho_ofstream.close();
  phi_ofstream.close();

  x_ofstream.close();
  u_x_ofstream.close();
  u_y_ofstream.close();

  return 0;

}   
// End of main

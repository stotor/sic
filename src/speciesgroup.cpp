#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "speciesgroup.hpp"
#include "particlespecies.hpp"
#include "utilities.hpp"

void SpeciesGroup::initial_velocity_deceleration(std::vector<double> &e_x_int, 
						 std::vector<double> &e_y,
						 std::vector<double> &b_z_tavg)
{
  for (int i = 0; i < n_species; i++) {
    species[i].initial_velocity_deceleration(e_x_int, e_y, b_z_tavg);
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
      species[i].deposit_j_x(j_x);
    }
    else if (method==1) {
      species[i].deposit_j_x_segments_zero(j_x);
    }
    else if (method==2) {
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
      species[i].deposit_j_y(j_y);
    }
    else if (method==1) {
      species[i].deposit_j_y_segments_zero(j_y);
    }
    else if (method==2) {
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
				    std::vector<double> b_z_tavg)
{
  for (int i = 0; i < n_species; i++) {
    species[i].advance_velocity(e_x_int, e_y, b_z_tavg);
  }
  return;
}


void SpeciesGroup::save_u_x_old()
{
  for (int i = 0; i < n_species; i++) {
    save_old_values(species[i].u_x, species[i].u_x_old, species[i].n_p);
  }
  return;
}

void SpeciesGroup::save_u_y_old()
{
  for (int i = 0; i < n_species; i++) {
    save_old_values(species[i].u_y, species[i].u_y_old, species[i].n_p);
  }
  return;
}

void SpeciesGroup::initialize_species(int n_g, double n_ppc, double dx)
{
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
  // double wave_amplitude = 0.025;
  // int wave_mode = 1;
  // for (int i_species = 0; i_species < n_species; i_species++) {
  //   for (int i_particle = 0; i_particle < species[i_species].n_p; i_particle++) {
  //     species[i_species].relativistic = true;
  //     species[i_species].charge[i_particle] = (-1.0) * (1.0 / n_ppc);
  //     species[i_species].rqm[i_particle] = -1.0;
  //     species[i_species].x[i_particle] = (double(i_particle) / species[i_species].n_p) * n_g * dx + (double(n_g) * dx / double(species[i_species].n_p)) / 2.0;
  //     species[i_species].u_x[i_particle] = wave_amplitude * cos(2.0 * PI * double(wave_mode) * species[i_species].x[i_particle] / (n_g * dx));
  //     species[i_species].u_y[i_particle] = 0.0;
  //   }
  // }

  // Electromagnetic wave initialization
  double v1 = 0.0025;
  double u1 = v1 / sqrt(1.0-v1*v1);
  int wave_mode = 1;
  for (int i_species = 0; i_species < n_species; i_species++) {
    for (int i_particle = 0; i_particle < species[i_species].n_p; i_particle++) {
      species[i_species].relativistic = true;
      species[i_species].charge[i_particle] = (-1.0) * (1.0 / n_ppc);
      species[i_species].rqm[i_particle] = -1.0;
      species[i_species].x[i_particle] = (double(i_particle) / species[i_species].n_p) * n_g * dx + (double(n_g) * dx / double(species[i_species].n_p)) / 2.0;
      species[i_species].u_x[i_particle] = 0.0;
      species[i_species].u_y[i_particle] = u1 * sin(2.0 * PI * double(wave_mode) * species[i_species].x[i_particle] / (n_g * dx));

    }
  }

  // Add ghost tracer particle if using line segments
  if ((method==1)||(method==2)) {
    for (int i = 0; i < n_species; i++) {
      species[i].charge.push_back(0.0);
      species[i].rqm.push_back(species[i].rqm[0]);
      species[i].u_x.push_back(species[i].u_x[0]);
      species[i].u_y.push_back(species[i].u_y[0]);
      species[i].x.push_back(species[i].x[0] + n_g*dx);
      species[i].n_p += 1;

      species[i].x_old.push_back(0);
      species[i].u_x_old.push_back(0);
      species[i].u_y_old.push_back(0);
    }
  }

  return;
}

void SpeciesGroup::deposit_rho(std::vector<double> &rho, int n_g)
{
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
  return;
}

void SpeciesGroup::advance_x()
{
  for (int i = 0; i < n_species; i++) {
    species[i].advance_x();
  }
  return;
}

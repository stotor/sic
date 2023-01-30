// Make sure all my routines are working for the case where ix[i] and ix[i+1] are separated by more than one.
//  average routine doesn't, but looks like its only called in the correct limits

// Check right maxs in old code
// Big to do is fully understanding right max
// If too high, probably doesn't make a difference just extra expense

//delete the unused ngp function
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "mpi.h"

#include "particlespecies.hpp"
#include "utilities.hpp"

//ix: dx [0, 1)
void find_ngp(int &ix, double &dx, int &ngp, double &delta)
{
  if (dx<0.5) {
    ngp = ix;
    delta = dx;
  } else {
    ngp = ix + 1;
    delta = dx - 1.0;
  }
  return;
}

double a_minus_b(int &ix_a, double &dx_a,
		 int &ix_b, double &dx_b)
{
  return (ix_a-ix_b) + (dx_a-dx_b);
}

void average(int &ix_a, double &dx_a,
	     int &ix_b, double &dx_b,
	     int &ix_c, double &dx_c)
{
  if (ix_a == ix_b) {
    ix_c = ix_a;
    dx_c = (dx_a + dx_b) / 2.0;
  } else if (ix_a > ix_b) {
    if (dx_a >= (1.0-dx_b)) {
      ix_c = ix_a;
      dx_c = (dx_a - (1.0-dx_b)) / 2.0;
    } else {
      ix_c = ix_b;
      dx_c = (dx_b + (1.0+dx_a)) / 2.0;
    }
  } else {
    if (dx_b >= (1.0-dx_a)) {
      ix_c = ix_b;
      dx_c = (dx_b - (1.0-dx_a)) / 2.0;
    } else {
      ix_c = ix_a;
      dx_c = (dx_a + (1.0+dx_b)) / 2.0;
    }
  }
  return;
}

double interpolate_field_0(std::vector<double> &field,
			   int ix,
			   double x,
			   int n_g)
{
  int ngp;
  double interpolated_field, delta;
  find_ngp(ix, x, ngp, delta);
  interpolated_field = field[mod(ngp, n_g)];
  return interpolated_field;
}

double interpolate_field_1(std::vector<double> &field,
			   int ix,
			   double x,
			   int n_g)
{
  double interpolated_field, w1, w2, delta;
  delta = x - 0.5;
  
  w1 = 0.5 - delta;
  w2 = 0.5 + delta;

  interpolated_field = w1 * field[mod(ix,n_g)]
    + w2 * field[mod((ix+1),n_g)];

  return interpolated_field;
}

double interpolate_field_2(std::vector<double> &field,
			   int ix,
			   double x,
			   int n_g)
{
  int ngp;
  double interpolated_field, w1, w2, w3, delta;

  find_ngp(ix, x, ngp, delta);

  w1 = 0.5 * pow((0.5 - delta), 2);
  w2 = 0.75 - delta * delta;
  w3 = 0.5 * pow((0.5 + delta), 2);

  interpolated_field = w1 * field[mod((ngp-1),n_g)]
    + w2 * field[mod(ngp,n_g)]
    + w3 * field[mod((ngp+1),n_g)];

  return interpolated_field;
}

double interpolate_field_3(std::vector<double> &field,
			   int ix,
			   double x,
			   int n_g)
{
  double interpolated_field, w1, w2, w3, w4, delta;
  delta = x - 0.5;

  w1 = -1.0 * pow((-0.5 + delta), 3) / 6.0;
  w2 = (4.0 - 6.0 * pow((0.5 + delta), 2) + 3.0 * pow((0.5 + delta), 3)) / 6.0;
  w3 = (4.0 - 6.0 * pow((0.5 - delta), 2) + 3.0 * pow((0.5 - delta), 3)) / 6.0;
  w4 = pow((0.5 + delta), 3) / 6.0;

  interpolated_field = w1 * field[mod((ix-1),n_g)]
    + w2 * field[mod(ix,n_g)]
    + w3 * field[mod((ix+1),n_g)]
    + w4 * field[mod((ix+2),n_g)];    

  return interpolated_field;
}

double interpolate_field_4(std::vector<double> &field,
			   int ix,
			   double x,
			   int n_g)
{
  int ngp;
  double interpolated_field, w1, w2, w3, w4, w5, delta;

  find_ngp(ix, x, ngp, delta);

  w1 = pow((1.0 - 2.0*delta), 4) / 384.0;
  w2 = (19.0 - 44.0*delta + 24.0*pow(delta, 2) + 16.0*pow(delta,3) - 16.0*pow(delta, 4))/96.0;
  w3 = 0.5989583333333334 - (5.0*pow(delta, 2))/8.0 + pow(delta, 4)/4.0;
  w4 = (19.0 + 44.0*delta + 24.0*pow(delta, 2) - 16.0*pow(delta, 3) - 16*pow(delta,4))/96.0;
  w5 = pow((1.0 + 2.0*delta), 4) / 384.0;
  
  interpolated_field = w1 * field[mod((ngp-2),n_g)]
    + w2 * field[mod((ngp-1),n_g)]
    + w3 * field[mod(ngp,n_g)]
    + w4 * field[mod((ngp+1),n_g)]
    + w5 * field[mod((ngp+2),n_g)];  

  return interpolated_field;
}

double interpolate_field(std::vector<double> &field,
			 int ix,
			 double x,
			 int n_g,
			 bool shift,
			 int order)
{
  double interpolated_field;
  // Shift position only for fields where zero lies at the first half-integer gridpoint 
  if (shift) {
    x = x - 0.5;
    if (x < 0.0) {
      ix = ix - 1;
      x = x + 1.0;
    }
  }

  switch (order) {
  case 0 :
    interpolated_field = interpolate_field_0(field, ix, x, n_g);
    break;
  case 1 :
    interpolated_field = interpolate_field_1(field, ix, x, n_g);
    break;
  case 2 :
    interpolated_field = interpolate_field_2(field, ix, x, n_g);
    break;
  case 3 :
    interpolated_field = interpolate_field_3(field, ix, x, n_g);
    break;
  case 4 :
    interpolated_field = interpolate_field_4(field, ix, x, n_g);
    break;
  default:
    std::cout << "Error, selected interpolation order not implemented." << std::endl;
  }
  
  // Index of left bounding gridpoint
  return interpolated_field;
}

void ParticleSpecies::initialize_species(int species_number, 
					 double n_ppc, 
					 int my_rank,
					 int num_procs,
					 int method,
					 int simulation_type)
{
  this->n_p = round((n_ppc * n_g) / num_procs);
  this->n_ppp = n_p;
  this->interp_order = method % 5;

  long long i_start, i_end;
  i_start = n_p * my_rank;
  i_end = i_start + n_p;

  ix.resize(n_p);  
  x.resize(n_p);
  u_x.resize(n_p);
  u_y.resize(n_p);
  u_z.resize(n_p);
  ix_old.resize(n_p);  
  x_old.resize(n_p);
  charge.resize(n_p);
  lagrangian_id.resize(n_p);
  
  std::stringstream ss;
  ss << "particles_";
  ss << species_number;
  species_name = ss.str();

  double particle_spacing = 1.0 / n_ppc;
  
  // Add ghost tracer particles if using line segments
  if (method>4) {
    for (int i = 0; i < 2; i++) {
      charge.push_back(0.0);
      u_x.push_back(0.0);
      u_y.push_back(0.0);
      u_z.push_back(0.0);      
      ix.push_back(0);
      ix_old.push_back(0);
      x.push_back(0.0);
      x_old.push_back(0.0);
      lagrangian_id.push_back(0.0);
    }
  }

  std::vector<double> density, rqm_vector, u_x_drift, u_y_drift;
  double k = 0.0;  
  double u_x_1 = 0.0;

  if (simulation_type == -1) {
    // Electrostatic wave
    density.push_back(-1.0);
    rqm_vector.push_back(-1.0);
    u_x_drift.push_back(0.0);
    u_y_drift.push_back(0.0);
    k = (2.0 * (4.0*atan(1.0)) / (n_g));
    u_x_1 = 0.0025;
  } else if (simulation_type == 0) {
    // Two-stream instability
    density.push_back(-1.0);
    density.push_back(-1.0);
    rqm_vector.push_back(-1.0);
    rqm_vector.push_back(-1.0);
    u_x_drift.push_back(0.176425);
    u_x_drift.push_back(-0.176425);
    u_y_drift.push_back(0.0);
    u_y_drift.push_back(0.0);
    
    k = 10.0 * (2.0 * (4.0*atan(1.0)) / (n_g));
    u_x_1 = 0.00176425;
    
  } else if (simulation_type == 1) {
    // Weibel
    density.push_back(-1.0);
    density.push_back(-1.0);
    rqm_vector.push_back(-1.0);
    rqm_vector.push_back(-1.0);
    u_x_drift.push_back(0.1);
    u_x_drift.push_back(-0.1);
    u_y_drift.push_back(0.5);
    u_y_drift.push_back(-0.5);
  } else if (simulation_type == 2) {
    // Whistler heating
    density.push_back(-1.0);
    density.push_back(1.0);
    rqm_vector.push_back(-1.0);
    rqm_vector.push_back(1836.0);
    u_x_drift.push_back(0.0);
    u_x_drift.push_back(0.0);
    u_y_drift.push_back(0.0);
    u_y_drift.push_back(0.0);    
  } else if (simulation_type == 3) {
    // Weibel highres
    density.push_back(-1.0);
    density.push_back(-1.0);
    rqm_vector.push_back(-1.0);
    rqm_vector.push_back(-1.0);
    u_x_drift.push_back(0.1);
    u_x_drift.push_back(-0.1);
    u_y_drift.push_back(0.5);
    u_y_drift.push_back(-0.5);

    k = (2.0 * (4.0*atan(1.0)) / (n_g));
    u_x_1 = 0.00176425;
    
  } else if (simulation_type == -2) {
    // Two-stream instability single mode
    density.push_back(-1.0);
    density.push_back(-1.0);
    rqm_vector.push_back(-1.0);
    rqm_vector.push_back(-1.0);
    u_x_drift.push_back(0.176425);
    u_x_drift.push_back(-0.176425);
    u_y_drift.push_back(0.0);
    u_y_drift.push_back(0.0);
    
    k = (2.0 * (4.0*atan(1.0)) / (n_g));
    u_x_1 = 0.00176425;
 }

  this->rqm = rqm_vector[species_number];
  this->method = method;

  for (long long i = i_start; i < i_end; i++) {
    lagrangian_id[i-i_start] = i - i_start;
    charge[i-i_start] = density[species_number] * (1.0 / n_ppc);

    ix[i-i_start] = i / (long long)n_ppc;
    ix_old[i-i_start] = i / (long long)n_ppc;
    
    x[i-i_start] = particle_spacing / 2.0 + particle_spacing *  (i%(long long)n_ppc);
    x_old[i-i_start] = particle_spacing / 2.0 + particle_spacing * (i%(long long)n_ppc);
    
    u_x[i-i_start] = u_x_drift[species_number];    
    u_y[i-i_start] = u_y_drift[species_number];
    u_z[i-i_start] = 0.0;
    
    if (simulation_type == -2 or simulation_type == 0 or simulation_type==-1 or simulation_type==3) {
      u_x[i-i_start] += u_x_1 * sin(k * ix[i-i_start] + k * x[i-i_start]);
    } else if (simulation_type == 2 and species_number == 0) {
      u_z[i-i_start] += (-2.0 * 0.273055) * cos(0.44526860656 * x[i-i_start] * dx);
    }
  }

  return;
}

void ParticleSpecies::advance_velocity(std::vector<double> &e_x, 
				       std::vector<double> &e_y,
				       std::vector<double> &e_z,
				       std::vector<double> &b_x,
				       std::vector<double> &b_y,
				       std::vector<double> &b_z)
{
  std::vector<double> e_x_int(n_g);
  std::vector<double> b_y_int(n_g);  
  std::vector<double> b_z_int(n_g);
  // Center to eliminate self-forces
  if (center_fields) {
    half_int_to_int(e_x, e_x_int, n_g);
    half_int_to_int(b_y, b_y_int, n_g);
    half_int_to_int(b_z, b_z_int, n_g);
  }

  double e_x_particle, e_y_particle, e_z_particle,
    b_x_particle, b_y_particle, b_z_particle, u_temp_x, u_temp_y, u_temp_z,
    gamma, t_norm, s_norm;
  double energy = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  double momentum_z = 0.0;  
  bool shift;
  
  for (int i = 0; i < n_p; i++) {
    b_x_particle = b_x[0]; // Uniform field

    if (center_fields) {
      shift = false;
      e_x_particle = interpolate_field(e_x_int, ix[i], x[i], n_g, shift, interp_order);
      e_y_particle = interpolate_field(e_y, ix[i], x[i], n_g, shift, interp_order);
      e_z_particle = interpolate_field(e_z, ix[i], x[i], n_g, shift, interp_order);
      b_y_particle = interpolate_field(b_y_int, ix[i], x[i], n_g, shift, interp_order);
      b_z_particle = interpolate_field(b_z_int, ix[i], x[i], n_g, shift, interp_order);
    } else {
      shift = true;
      e_x_particle = interpolate_field(e_x, ix[i], x[i], n_g, shift, interp_order);
      b_y_particle = interpolate_field(b_y, ix[i], x[i], n_g, shift, interp_order);
      b_z_particle = interpolate_field(b_z, ix[i], x[i], n_g, shift, interp_order);
      shift = false;
      e_y_particle = interpolate_field(e_y, ix[i], x[i], n_g, shift, interp_order);
      e_z_particle = interpolate_field(e_z, ix[i], x[i], n_g, shift, interp_order);
    }

    u_x[i] = u_x[i] + e_x_particle * (dt / 2.0) / rqm;
    u_y[i] = u_y[i] + e_y_particle * (dt / 2.0) / rqm;
    u_z[i] = u_z[i] + e_z_particle * (dt / 2.0) / rqm;

    gamma = sqrt(1.0 + pow(u_x[i], 2) + pow(u_y[i], 2) + pow(u_z[i], 2));

    // Calculate time centered energy for diagnostic purposes
    energy = energy + (charge[i] * rqm) * (pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0)) / (1.0 + gamma);

    t_norm = 0.5 * dt / (rqm * gamma);

    b_x_particle = b_x_particle * t_norm;
    b_y_particle = b_y_particle * t_norm;
    b_z_particle = b_z_particle * t_norm;

    u_temp_x = u_x[i] + u_y[i] * b_z_particle - u_z[i] * b_y_particle;
    u_temp_y = u_y[i] + u_z[i] * b_x_particle - u_x[i] * b_z_particle;
    u_temp_z = u_z[i] + u_x[i] * b_y_particle - u_y[i] * b_x_particle;
    
    s_norm = 2.0 / (1.0 + pow(b_x_particle, 2) + pow(b_y_particle, 2) + pow(b_z_particle, 2));
    
    b_x_particle = b_x_particle * s_norm;
    b_y_particle = b_y_particle * s_norm;
    b_z_particle = b_z_particle * s_norm;

    u_x[i] = u_x[i] + u_temp_y * b_z_particle - u_temp_z * b_y_particle;
    u_y[i] = u_y[i] + u_temp_z * b_x_particle - u_temp_x * b_z_particle;
    u_z[i] = u_z[i] + u_temp_x * b_y_particle - u_temp_y * b_x_particle;

    u_x[i] = u_x[i] + e_x_particle * (dt / 2.0) / rqm;
    u_y[i] = u_y[i] + e_y_particle * (dt / 2.0) / rqm;
    u_z[i] = u_z[i] + e_z_particle * (dt / 2.0) / rqm;

    momentum_x += charge[i] * rqm * u_x[i];
    momentum_y += charge[i] * rqm * u_y[i];
    momentum_z += charge[i] * rqm * u_z[i];
  }
  
  energy_history.push_back(energy);
  momentum_x_history.push_back(momentum_x);
  momentum_y_history.push_back(momentum_y);
  momentum_z_history.push_back(momentum_z);  
  n_p_history.push_back(n_p);
  return;
}

void ParticleSpecies::initial_velocity_deceleration(std::vector<double> &e_x,
						    std::vector<double> &e_y,
						    std::vector<double> &e_z,
						    std::vector<double> &b_x,
						    std::vector<double> &b_y,
						    std::vector<double> &b_z)
{
  std::vector<double> e_x_int(n_g);
  std::vector<double> b_y_int(n_g);  
  std::vector<double> b_z_int(n_g);
  // Center to eliminate self-forces
  if (center_fields) {
    half_int_to_int(e_x, e_x_int, n_g);
    half_int_to_int(b_y, b_y_int, n_g);
    half_int_to_int(b_z, b_z_int, n_g);
  }

  double e_x_particle, e_y_particle, e_z_particle,
    b_x_particle, b_y_particle, b_z_particle, u_temp_x, u_temp_y, u_temp_z,
    gamma, t_norm, s_norm;
  bool shift;
  //  double energy = 0.0;
  //  double momentum_x = 0.0;
  //  double momentum_y = 0.0;
  //  double momentum_z = 0.0;  
  
  for (int i = 0; i < n_p; i++) {
    b_x_particle = b_x[0]; // Uniform field

    if (center_fields) {
      shift = false;
      e_x_particle = interpolate_field(e_x_int, ix[i], x[i],  n_g, shift, interp_order);
      e_y_particle = interpolate_field(e_y, ix[i], x[i], n_g, shift, interp_order);
      e_z_particle = interpolate_field(e_z, ix[i], x[i], n_g, shift, interp_order);
      b_y_particle = interpolate_field(b_y_int, ix[i], x[i], n_g, shift, interp_order);
      b_z_particle = interpolate_field(b_z_int, ix[i], x[i], n_g, shift, interp_order);
    } else {
      shift = true;
      e_x_particle = interpolate_field(e_x, ix[i], x[i], n_g, shift, interp_order);
      b_y_particle = interpolate_field(b_y, ix[i], x[i], n_g, shift, interp_order);
      b_z_particle = interpolate_field(b_z, ix[i], x[i], n_g, shift, interp_order);
      shift = false;
      e_y_particle = interpolate_field(e_y, ix[i], x[i], n_g, shift, interp_order);
      e_z_particle = interpolate_field(e_z, ix[i], x[i], n_g, shift, interp_order);
    }

    
    gamma = sqrt(1.0 + pow(u_x[i], 2) + pow(u_y[i], 2) + pow(u_z[i], 2));

    // Calculate time centered energy for diagnostic purposes
    //    energy = energy + (charge[i] * rqm) * (pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0)) / (1.0 + gamma);

    t_norm = 0.5 * (-0.5 * dt) / (rqm * gamma);

    b_x_particle = b_x_particle * t_norm;
    b_y_particle = b_y_particle * t_norm;
    b_z_particle = b_z_particle * t_norm;

    u_temp_x = u_x[i] + u_y[i] * b_z_particle - u_z[i] * b_y_particle;
    u_temp_y = u_y[i] + u_z[i] * b_x_particle - u_x[i] * b_z_particle;
    u_temp_z = u_z[i] + u_x[i] * b_y_particle - u_y[i] * b_x_particle;
    
    s_norm = 2.0 / (1.0 + pow(b_x_particle, 2) + pow(b_y_particle, 2) + pow(b_z_particle, 2));
    
    b_x_particle = b_x_particle * s_norm;
    b_y_particle = b_y_particle * s_norm;
    b_z_particle = b_z_particle * s_norm;

    u_x[i] = u_x[i] + u_temp_y * b_z_particle - u_temp_z * b_y_particle;
    u_y[i] = u_y[i] + u_temp_z * b_x_particle - u_temp_x * b_z_particle;
    u_z[i] = u_z[i] + u_temp_x * b_y_particle - u_temp_y * b_x_particle;

    u_x[i] = u_x[i] + e_x_particle * (-1.0 * dt / 2.0) / rqm;
    u_y[i] = u_y[i] + e_y_particle * (-1.0 * dt / 2.0) / rqm;
    u_z[i] = u_z[i] + e_z_particle * (-1.0 * dt / 2.0) / rqm;

    //    momentum_x += charge[i] * rqm * u_x[i];
    //    momentum_y += charge[i] * rqm * u_y[i];
    //    momentum_z += charge[i] * rqm * u_z[i];
  }
  
  //  energy_history.push_back(energy);
  //  momentum_x_history.push_back(momentum_x);
  //  momentum_y_history.push_back(momentum_y);
  //  momentum_z_history.push_back(momentum_z);  
  //  n_p_history.push_back(n_p);

  return;
}

void ParticleSpecies::deposit_rho_pic_0(std::vector<double> &rho)
{
  int ngp, ix_n;
  double x_n, delta;

  for (int i = 0; i < n_p; i++) {
    ix_n = ix[i];
    x_n = x[i];
    find_ngp(ix_n, x_n, ngp, delta);
    rho[mod(ngp,n_g)] += charge[i];
  }
  return;
}

void ParticleSpecies::deposit_rho_pic_1(std::vector<double> &rho)
{
  int ix_n;
  double x_n, w1, w2;

  for (int i = 0; i < n_p; i++) {
    ix_n = ix[i];    
    x_n = x[i];
    
    w1 = 1.0 - x_n;
    w2 = x_n;

    rho[mod(ix_n,n_g)] += charge[i] * w1;
    rho[mod((ix_n+1),n_g)] += charge[i] * w2;
  }
  return;
}

void ParticleSpecies::deposit_rho_pic_2(std::vector<double> &rho)
{
  int ngp, ix_n;
  double x_n, w1, w2, w3, delta;

  for (int i = 0; i < n_p; i++) {
    ix_n = ix[i];
    x_n = x[i];
    find_ngp(ix_n, x_n, ngp, delta);
    
    w1 = 0.5 * pow((0.5 - delta), 2);
    w2 = 0.75 - delta * delta;
    w3 = 0.5 * pow((0.5 + delta), 2);

    rho[mod((ngp-1),n_g)] += charge[i] * w1;
    rho[mod(ngp,n_g)] += charge[i] * w2;
    rho[mod((ngp+1),n_g)] += charge[i] * w3;
  }
  return;
}

void ParticleSpecies::deposit_rho_pic_3(std::vector<double> &rho)
{
  int ix_n;
  double x_n, w1, w2, w3, w4, delta;

  for (int i = 0; i < n_p; i++) {
    ix_n = ix[i];    
    x_n = x[i];
    delta = x_n - 0.5;

    w1 = -1.0 * pow((-0.5 + delta), 3) / 6.0;
    w2 = (4.0 - 6.0 * pow((0.5 + delta), 2) + 3.0 * pow((0.5 + delta), 3)) / 6.0;
    w3 = (4.0 - 6.0 * pow((0.5 - delta), 2) + 3.0 * pow((0.5 - delta), 3)) / 6.0;
    w4 = pow((0.5 + delta), 3) / 6.0;

    rho[mod((ix_n-1),n_g)] += charge[i] * w1;
    rho[mod(ix_n,n_g)] += charge[i] * w2;
    rho[mod((ix_n+1),n_g)] += charge[i] * w3;
    rho[mod((ix_n+2),n_g)] += charge[i] * w4;
  }
  return;
}

void ParticleSpecies::deposit_rho_pic_4(std::vector<double> &rho)
{
  int ngp, ix_n;
  double x_n, w1, w2, w3, w4, w5, delta;

  for (int i = 0; i < n_p; i++) {
    ix_n = ix[i];
    x_n = x[i];
    find_ngp(ix_n, x_n, ngp, delta);

    w1 = pow((1.0 - 2.0*delta), 4) / 384.0;
    w2 = (19.0 - 44.0*delta + 24.0*pow(delta, 2) + 16.0*pow(delta,3) - 16.0*pow(delta, 4))/96.0;
    w3 = (115.0 / 192.0) - (5.0*pow(delta, 2))/8.0 + pow(delta, 4)/4.0;
    w4 = (19.0 + 44.0*delta + 24.0*pow(delta, 2) - 16.0*pow(delta, 3) - 16*pow(delta,4))/96.0;
    w5 = pow((1.0 + 2.0*delta), 4) / 384.0;

    rho[mod((ngp-2),n_g)] += charge[i] * w1;
    rho[mod((ngp-1),n_g)] += charge[i] * w2;
    rho[mod(ngp,n_g)] += charge[i] * w3;
    rho[mod((ngp+1),n_g)] += charge[i] * w4;
    rho[mod((ngp+2),n_g)] += charge[i] * w5;
  }
  return;
}

void ParticleSpecies::deposit_j_x_pic_0(std::vector<double> &j_x)
{
  int ngp_i, ngp_f;  
  double delta_i, delta_f;
  for (int i = 0; i < n_p; i++) {
    find_ngp(ix_old[i], x_old[i], ngp_i, delta_i);
    find_ngp(ix[i], x[i], ngp_f, delta_f);

    if (ngp_i < ngp_f) {
      j_x[mod(ngp_i,n_g)] += charge[i];
    } else if (ngp_i > ngp_f) {
      j_x[mod(ngp_f,n_g)] -= charge[i];
    }
  }
  return;
}

void ParticleSpecies::deposit_j_x_pic_1(std::vector<double> &j_x)
{
  int ix_i, ix_f;
  double x_i, x_f;
  for (int i = 0; i < n_p; i++) {
      ix_i = ix_old[i];
      ix_f = ix[i];    
      x_i = x_old[i];
      x_f = x[i];
      if (ix_i==ix_f) {
	j_x[mod(ix_i,n_g)] += charge[i] * (x_f - x_i);
      } 
      else if (ix_i < ix_f) {
	j_x[mod(ix_i,n_g)] += charge[i] * (1.0 - x_i);
	j_x[mod(ix_f,n_g)] += charge[i] * x_f;
      }
      else {
	j_x[mod(ix_i,n_g)] += charge[i] * (-1.0 * x_i);
	j_x[mod(ix_f,n_g)] += charge[i] * (x_f - 1.0);
      }
  }
  return;
}

void ParticleSpecies::deposit_j_x_pic_2(std::vector<double> &j_x)
{
  double x0, x1, xa, xb, q_norm, delta0, delta1;
  int ngp0, ngp1, shift;
  for (int i = 0; i < n_p; i++) {
    q_norm = charge[i];

    find_ngp(ix_old[i], x_old[i], ngp0, delta0);
    find_ngp(ix[i], x[i], ngp1, delta1);
    shift = ngp1 - ngp0;
    x0 = delta0;
    x1 = delta1 + shift;
    
    if (shift==0) {
      xa = x0;
      xb = x1;
      j_x[mod(ngp0-1,n_g)] += (q_norm * (xa - xb)*(-1.0 + xa + xb)) / 2.0;
      j_x[mod(ngp0,n_g)] += (q_norm * (-(xa*(1.0 + xa)) + xb + xb*xb)) / 2.0;
    } 
    else {
      xa = x0;
      xb = shift * 0.5;      
      j_x[mod(ngp0-1,n_g)] += (q_norm * (xa - xb)*(-1.0 + xa + xb)) / 2.0;
      j_x[mod(ngp0,n_g)] += (q_norm * (-(xa*(1.0 + xa)) + xb + xb*xb)) / 2.0;

      xa = -1.0 * xb;
      xb = x1 - shift;
      j_x[mod(ngp1-1,n_g)] += (q_norm * (xa - xb)*(-1.0 + xa + xb)) / 2.0;
      j_x[mod(ngp1,n_g)] += (q_norm * (-(xa*(1.0 + xa)) + xb + xb*xb)) / 2.0;

    }
  }

  return;
}

void ParticleSpecies::deposit_j_x_pic_3(std::vector<double> &j_x)
{
  double x0, x1, xa, xb, q_norm;
  int i_l_0, i_l_1, shift;
  for (int i = 0; i < n_p; i++) {
    q_norm = charge[i];
    x0 = x_old[i];
    x1 = x[i];
    i_l_0 = ix_old[i];
    i_l_1 = ix[i];
    shift = i_l_1 - i_l_0;
    x0 = x0 - 0.5;
    x1 = x1 - 0.5 + shift;

    if (shift==0) {
      xa = x0;
      xb = x1;
      j_x[mod((i_l_0-1),n_g)] += -1.0 * (q_norm * (pow((-0.5 + xa), 3) - pow((-0.5 + xb), 3))) / 6.0;
      j_x[mod(i_l_0,n_g)] += (q_norm * (-9.0*xa + 4.0*pow(xa, 3) + 9.0*xb - 4.0*pow(xb, 3))) / 12.0;
      j_x[mod((i_l_0+1),n_g)] += (q_norm * (-1.0*(xa*(3.0 + 6.0*xa + 4.0*xa*xa)) + xb*(3.0 + 6.0*xb + 4.0*xb*xb))) / 24.0;
    } 
    else {
      xa = x0;
      xb = shift * 0.5;
      j_x[mod((i_l_0-1),n_g)] += -1.0 * (q_norm * (pow((-0.5 + xa), 3) - pow((-0.5 + xb), 3))) / 6.0;
      j_x[mod(i_l_0,n_g)] += (q_norm * (-9.0*xa + 4.0*pow(xa, 3) + 9.0*xb - 4.0*pow(xb, 3))) / 12.0;
      j_x[mod((i_l_0+1),n_g)] += (q_norm * (-1.0*(xa*(3.0 + 6.0*xa + 4.0*xa*xa)) + xb*(3.0 + 6.0*xb + 4.0*xb*xb))) / 24.0;

      xa = -1.0 * xb;
      xb = x1 - shift;
      j_x[mod((i_l_1-1),n_g)] += -1.0 * (q_norm * (pow((-0.5 + xa), 3) - pow((-0.5 + xb), 3))) / 6.0;
      j_x[mod(i_l_1,n_g)] += (q_norm * (-9.0*xa + 4.0*pow(xa, 3) + 9.0*xb - 4.0*pow(xb, 3))) / 12.0;
      j_x[mod((i_l_1+1),n_g)] += (q_norm * (-1.0*(xa*(3.0 + 6.0*xa + 4.0*xa*xa)) + xb*(3.0 + 6.0*xb + 4.0*xb*xb))) / 24.0;
    }
  }

  return;
}


void ParticleSpecies::deposit_j_x_pic_4(std::vector<double> &j_x)
{
  double x0, x1, xa, xb, q_norm, delta0, delta1;
  int ngp0, ngp1, shift;
  for (int i = 0; i < n_p; i++) {
    q_norm = charge[i];
    find_ngp(ix_old[i], x_old[i], ngp0, delta0);
    find_ngp(ix[i], x[i], ngp1, delta1);    
    shift = ngp1 - ngp0;
    x0 = delta0;
    x1 = delta1 + shift;
    
    if (shift==0) {
      xa = x0;
      xb = x1;
      j_x[mod((ngp0-2),n_g)] += -1.0 * (q_norm*(-1.0 * pow((1.0 - 2.0*xa), 4) + pow((1.0 - 2.0*xb), 4))) / 384.0;
      j_x[mod((ngp0-1),n_g)] += (q_norm*(xa*(-23.0 + xa*(15.0 + 4.0*xa - 6.0*xa*xa)) + xb*(23.0 + xb*(-15.0 - 4.0*xb + 6.0*xb*xb)))) / 48.0;
      j_x[mod(ngp0,n_g)] += (q_norm*(xa*(-23.0 + xa*(-15.0 + 4.0*xa + 6.0*xa*xa)) + xb*(23.0 + xb*(15.0 - 2.0*xb*(2.0 + 3.0*xb))))) / 48.0;
      j_x[mod((ngp0+1),n_g)] += -1.0 * (q_norm*(xa - xb)*(1.0 + xa + xb)*(1.0 + 2.0*xa*(1.0 + xa) + 2.0*xb*(1.0 + xb))) / 48.0;
    } 
    else {
      xa = x0;
      xb = shift * 0.5;
      j_x[mod((ngp0-2),n_g)] += -1.0 * (q_norm*(-1.0 * pow((1.0 - 2.0*xa), 4) + pow((1.0 - 2.0*xb), 4))) / 384.0;
      j_x[mod((ngp0-1),n_g)] += (q_norm*(xa*(-23.0 + xa*(15.0 + 4.0*xa - 6.0*xa*xa)) + xb*(23.0 + xb*(-15.0 - 4.0*xb + 6.0*xb*xb)))) / 48.0;
      j_x[mod(ngp0,n_g)] += (q_norm*(xa*(-23.0 + xa*(-15.0 + 4.0*xa + 6.0*xa*xa)) + xb*(23.0 + xb*(15.0 - 2.0*xb*(2.0 + 3.0*xb))))) / 48.0;
      j_x[mod((ngp0+1),n_g)] += -1.0 * (q_norm*(xa - xb)*(1.0 + xa + xb)*(1.0 + 2.0*xa*(1.0 + xa) + 2.0*xb*(1.0 + xb))) / 48.0;

      xa = -1.0 * xb;
      xb = x1 - shift;
      j_x[mod((ngp1-2),n_g)] += -1.0 * (q_norm*(-1.0 * pow((1.0 - 2.0*xa), 4) + pow((1.0 - 2.0*xb), 4))) / 384.0;
      j_x[mod((ngp1-1),n_g)] += (q_norm*(xa*(-23.0 + xa*(15.0 + 4.0*xa - 6.0*xa*xa)) + xb*(23.0 + xb*(-15.0 - 4.0*xb + 6.0*xb*xb)))) / 48.0;
      j_x[mod(ngp1,n_g)] += (q_norm*(xa*(-23.0 + xa*(-15.0 + 4.0*xa + 6.0*xa*xa)) + xb*(23.0 + xb*(15.0 - 2.0*xb*(2.0 + 3.0*xb))))) / 48.0;
      j_x[mod((ngp1+1),n_g)] += -1.0 * (q_norm*(xa - xb)*(1.0 + xa + xb)*(1.0 + 2.0*xa*(1.0 + xa) + 2.0*xb*(1.0 + xb))) / 48.0;
    }
  }

  return;
}


void ParticleSpecies::deposit_j_y_pic_0(std::vector<double> &j_y)
{
  int ngp, ix_tavg;
  double x_tavg, gamma, j_y_i, delta;

  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_y_i = charge[i] * u_y[i] / gamma;

    average(ix[i], x[i], ix_old[i], x_old[i], ix_tavg, x_tavg);
    find_ngp(ix_tavg, x_tavg, ngp, delta);

    j_y[mod(ngp,n_g)] += j_y_i;
  }
  return;
}

void ParticleSpecies::deposit_j_y_pic_1(std::vector<double> &j_y)
{
  int ix_tavg;
  double x_tavg, w1, w2, gamma, j_y_i;

  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_y_i = charge[i] * u_y[i] / gamma;

    average(ix[i], x[i], ix_old[i], x_old[i], ix_tavg, x_tavg);
    
    w1 = 1.0 - x_tavg;
    w2 = x_tavg;

    j_y[mod(ix_tavg,n_g)] += w1 * j_y_i;
    j_y[mod((ix_tavg+1),n_g)] += w2 * j_y_i;
  }
  return;
}

void ParticleSpecies::deposit_j_y_pic_2(std::vector<double> &j_y)
{
  int ngp, ix_tavg;
  double w1, w2, w3, delta, x_tavg, gamma, j_y_i;

  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_y_i = charge[i] * u_y[i] / gamma;

    average(ix[i], x[i], ix_old[i], x_old[i], ix_tavg, x_tavg);
    find_ngp(ix_tavg, x_tavg, ngp, delta);    
    
    w1 = 0.5 * pow((0.5 - delta), 2);
    w2 = 0.75 - delta * delta;
    w3 = 0.5 * pow((0.5 + delta), 2);

    j_y[mod((ngp-1),n_g)] += w1 * j_y_i;
    j_y[mod(ngp,n_g)] += w2 * j_y_i;
    j_y[mod((ngp+1),n_g)] += w3 * j_y_i;
  }
  return;
}

void ParticleSpecies::deposit_j_y_pic_3(std::vector<double> &j_y)
{
  int ix_tavg;
  double x_tavg, w1, w2, w3, w4, delta, gamma, j_y_i;

  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_y_i = charge[i] * u_y[i] / gamma;
    
    average(ix[i], x[i], ix_old[i], x_old[i], ix_tavg, x_tavg);
    delta = x_tavg - 0.5;

    w1 = -1.0 * pow((-0.5 + delta), 3) / 6.0;
    w2 = (4.0 - 6.0 * pow((0.5 + delta), 2) + 3.0 * pow((0.5 + delta), 3)) / 6.0;
    w3 = (23.0 + 30.0*delta - 12.0*pow(delta, 2) - 24.0*pow(delta,3)) / 48.0;
    w4 = pow((0.5 + delta), 3) / 6.0;

    j_y[mod((ix_tavg-1),n_g)] += w1 * j_y_i;
    j_y[mod(ix_tavg,n_g)] += w2 * j_y_i;
    j_y[mod((ix_tavg+1),n_g)] += w3 * j_y_i;
    j_y[mod((ix_tavg+2),n_g)] += w4 * j_y_i;
  }
  return;
}

void ParticleSpecies::deposit_j_y_pic_4(std::vector<double> &j_y)
{
  int ngp, ix_tavg;
  double w1, w2, w3, w4, w5, delta, x_tavg, gamma, j_y_i;

  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_y_i = charge[i] * u_y[i] / gamma;

    average(ix[i], x[i], ix_old[i], x_old[i], ix_tavg, x_tavg);
    find_ngp(ix_tavg, x_tavg, ngp, delta);    

    w1 = pow((1.0 - 2.0*delta), 4) / 384.0;
    w2 = (19.0 - 44.0*delta + 24.0*pow(delta, 2) + 16.0*pow(delta,3) - 16.0*pow(delta, 4))/96.0;
    w3 = 0.5989583333333334 - (5.0*pow(delta, 2))/8.0 + pow(delta, 4)/4.0;
    w4 = (19.0 + 44.0*delta + 24.0*pow(delta, 2) - 16.0*pow(delta, 3) - 16*pow(delta,4))/96.0;
    w5 = pow((1.0 + 2.0*delta), 4) / 384.0;

    j_y[mod((ngp-2),n_g)] += w1 * j_y_i;
    j_y[mod((ngp-1),n_g)] += w2 * j_y_i;
    j_y[mod(ngp,n_g)] += w3 * j_y_i;
    j_y[mod((ngp+1),n_g)] += w4 * j_y_i;
    j_y[mod((ngp+2),n_g)] += w5 * j_y_i;
  }
  return;
}

void ParticleSpecies::deposit_j_z_pic_0(std::vector<double> &j_z)
{
  int ngp, ix_tavg;
  double x_tavg, gamma, j_z_i, delta;

  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_z_i = charge[i] * u_z[i] / gamma;

    average(ix[i], x[i], ix_old[i], x_old[i], ix_tavg, x_tavg);
    find_ngp(ix_tavg, x_tavg, ngp, delta);
    
    j_z[mod(ngp,n_g)] += j_z_i;
  }
  return;
}

void ParticleSpecies::deposit_j_z_pic_1(std::vector<double> &j_z)
{
  int ix_tavg;
  double x_tavg, w1, w2, gamma, j_z_i;

  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_z_i = charge[i] * u_z[i] / gamma;

    average(ix[i], x[i], ix_old[i], x_old[i], ix_tavg, x_tavg);

    w1 = 1.0 - x_tavg;
    w2 = x_tavg;

    j_z[mod(ix_tavg,n_g)] += w1 * j_z_i;
    j_z[mod((ix_tavg+1),n_g)] += w2 * j_z_i;
  }
  return;
}

void ParticleSpecies::deposit_j_z_pic_2(std::vector<double> &j_z)
{
  int ngp, ix_tavg;
  double w1, w2, w3, delta, x_tavg, gamma, j_z_i;

  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_z_i = charge[i] * u_z[i] / gamma;

    average(ix[i], x[i], ix_old[i], x_old[i], ix_tavg, x_tavg);
    find_ngp(ix_tavg, x_tavg, ngp, delta);

    w1 = 0.5 * pow((0.5 - delta), 2);
    w2 = 0.75 - delta * delta;
    w3 = 0.5 * pow((0.5 + delta), 2);

    j_z[mod((ngp-1),n_g)] += w1 * j_z_i;
    j_z[mod(ngp,n_g)] += w2 * j_z_i;
    j_z[mod((ngp+1),n_g)] += w3 * j_z_i;
  }
  return;
}

void ParticleSpecies::deposit_j_z_pic_3(std::vector<double> &j_z)
{
  int ix_tavg;
  double w1, w2, w3, w4, delta, x_tavg, gamma, j_z_i;

  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_z_i = charge[i] * u_z[i] / gamma;

    average(ix[i], x[i], ix_old[i], x_old[i], ix_tavg, x_tavg);
    delta = x_tavg - 0.5;    

    w1 = -1.0 * pow((-0.5 + delta), 3) / 6.0;
    w2 = (4.0 - 6.0 * pow((0.5 + delta), 2) + 3.0 * pow((0.5 + delta), 3)) / 6.0;
    w3 = (23.0 + 30.0*delta - 12.0*pow(delta, 2) - 24.0*pow(delta,3)) / 48.0;
    w4 = pow((0.5 + delta), 3) / 6.0;

    j_z[mod((ix_tavg-1),n_g)] += w1 * j_z_i;
    j_z[mod(ix_tavg,n_g)] += w2 * j_z_i;
    j_z[mod((ix_tavg+1),n_g)] += w3 * j_z_i;
    j_z[mod((ix_tavg+2),n_g)] += w4 * j_z_i;
  }
  return;
}

void ParticleSpecies::deposit_j_z_pic_4(std::vector<double> &j_z)
{
  int ngp, ix_tavg;
  double w1, w2, w3, w4, w5, delta, x_tavg, gamma, j_z_i;

  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_z_i = charge[i] * u_z[i] / gamma;

    average(ix[i], x[i], ix_old[i], x_old[i], ix_tavg, x_tavg);
    find_ngp(ix_tavg, x_tavg, ngp, delta);

    w1 = pow((1.0 - 2.0*delta), 4) / 384.0;
    w2 = (19.0 - 44.0*delta + 24.0*pow(delta, 2) + 16.0*pow(delta,3) - 16.0*pow(delta, 4))/96.0;
    w3 = 0.5989583333333334 - (5.0*pow(delta, 2))/8.0 + pow(delta, 4)/4.0;
    w4 = (19.0 + 44.0*delta + 24.0*pow(delta, 2) - 16.0*pow(delta, 3) - 16*pow(delta,4))/96.0;
    w5 = pow((1.0 + 2.0*delta), 4) / 384.0;

    j_z[mod((ngp-2),n_g)] += w1 * j_z_i;
    j_z[mod((ngp-1),n_g)] += w2 * j_z_i;
    j_z[mod(ngp,n_g)] += w3 * j_z_i;
    j_z[mod((ngp+1),n_g)] += w4 * j_z_i;
    j_z[mod((ngp+2),n_g)] += w5 * j_z_i;
  }
  return;
}

void arrange_segment(int &ix_a, double &x_a,
		     int &ix_b, double &x_b,
		     int &ix_l, double &x_l,
		     int &ix_r, double &x_r,
		     double &length)
{
  length = a_minus_b(ix_a, x_a, ix_b, x_b);
  if (length<=0.0) {
    length = length * (-1.0);
    ix_l = ix_a;
    x_l = x_a;
    ix_r = ix_b;
    x_r = x_b;
  } else {
    ix_l = ix_b;
    x_l = x_b;
    ix_r = ix_a;
    x_r = x_a;
  }
  return;
}

void arrange_segment_velocity(int &ix_a, double &x_a,
			      int &ix_b, double &x_b,
			      double &v_a, double &v_b,
			      int &ix_l, double &x_l,
			      int &ix_r, double &x_r,
			      double &length,
			      double &v_l, double &v_r)
{
  length = a_minus_b(ix_a, x_a, ix_b, x_b);
  if (length<=0.0) {
    length = length * (-1.0);
    ix_l = ix_a;
    x_l = x_a;
    ix_r = ix_b;
    x_r = x_b;
    v_l = v_a;
    v_r = v_b;
  } else {
    ix_l = ix_b;
    x_l = x_b;
    ix_r = ix_a;
    x_r = x_a;
    v_l = v_b;
    v_r = v_a;
  }
  return;
}

void deposit_rho_segment_0(std::vector<double> &rho,
			   int ix_l,
			   int ix_r,
			   double x_l,
			   double x_r,
			   double length,
			   double charge,
			   int n_g)
{
  double charge_fraction;
  int bound_left, bound_right;

  bound_left = ix_l;
  bound_right = ix_r + 1;

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    if ((x_l == 0.5) && (x_r == 0.5)) {
      rho[mod((bound_right),n_g)] += charge;
    } else if (x_l >= 0.5) {
      rho[mod((bound_right),n_g)] += charge;
    } else if (x_r <= 0.5) {
      rho[mod((bound_left),n_g)] += charge;
    } else {
      rho[mod((bound_left),n_g)] += (0.5 - x_l) / length * charge;
      rho[mod((bound_right),n_g)] += (x_r - 0.5) / length * charge;
    }
  }
  else {
    // Left end
    if (x_l < 0.5) {
      rho[mod(bound_left,n_g)] += charge * (0.5 - x_l) / length;
      rho[mod((bound_left+1),n_g)] += charge * (0.5) / length;
    } else { 
      rho[mod(bound_left+1,n_g)] += charge * (1.0 - x_l) / length;
    }

    // Pieces connecting two gridpoints
    charge_fraction = (1.0) / length;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      rho[mod(cell,n_g)] += charge * charge_fraction * 0.5;
      rho[mod((cell+1),n_g)] += charge * charge_fraction * 0.5;
    }
    
    // Right end
    if (x_r > 0.5) {
      rho[mod(bound_right-1,n_g)] += charge * (0.5) / length;
      rho[mod(bound_right,n_g)] += charge * (x_r - 0.5) / length;
    } else { 
      rho[mod(bound_right-1,n_g)] += charge * x_r / length;
    }
  }
  return;
}

void deposit_charge_to_left_segment_0(std::vector<double> &j_x, // NEED TO CHECK <= stuff here
				      int ix_l, // Special case for zero length segment?
				      int ix_r,
				      double x_l, 
				      double x_r,
				      double length,
				      double charge,
				      int n_g,
				      int right_max)  
{
  double cell_boundary, weight;
  int bound_left, bound_right;

  bound_left = ix_l;

  // Loop over cell boundaries that could be affected
  // Calculate charge initially to the left of the current cell's right boundary
  for (int cell = bound_left; cell < right_max; cell++) {
    cell_boundary = cell + 0.5;
    if ((ix_r+x_r) <= cell_boundary) {
      weight = 1.0;
    } 
    else if ((ix_l+x_l) >= cell_boundary) {
      weight = 0.0;
    }
    else {
      weight = (cell_boundary - ix_l) / length - x_l / length;
    }
    j_x[mod(cell,n_g)] += charge * weight;
  }
  return;
}


void deposit_rho_segment_1(std::vector<double> &rho,
			   int ix_l,
			   int ix_r,
			   double x_l,
			   double x_r,
			   double length,
			   double charge,
			   int n_g)
{
  double midpoint, charge_fraction, weight;
  int bound_left, bound_right;

  bound_left = ix_l;
  bound_right = ix_r + 1;  

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    midpoint = (x_l + x_r) / 2.0;
    charge_fraction = 1.0;
    weight = 1.0 - midpoint;
    rho[mod(bound_left,n_g)] += charge * charge_fraction * weight;
    rho[mod((bound_left+1),n_g)] += charge * charge_fraction * (1.0 - weight);
  }
  else {
    // Left end
    midpoint = (x_l + 1.0) / 2.0;
    charge_fraction = (1.0 - x_l) / length;
    weight = (1.0 - midpoint);
    rho[mod(bound_left,n_g)] += charge * charge_fraction * weight;
    rho[mod((bound_left+1),n_g)] += charge * charge_fraction * (1.0 - weight);

    // Pieces connecting two gridpoints
    charge_fraction = 1.0 / length;

    weight = 0.5;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      rho[mod(cell,n_g)] += charge * charge_fraction * weight;
      rho[mod((cell+1),n_g)] += charge * charge_fraction * (1.0 - weight);
    }
    
    // Right end
    midpoint = x_r / 2.0;
    charge_fraction = x_r / length;
    weight = 1.0 - midpoint;
    rho[mod((bound_right-1),n_g)] += charge * charge_fraction * weight;
    rho[mod(bound_right,n_g)] += charge * charge_fraction * (1.0-weight);
  }
  return;
}

void deposit_charge_to_left_segment_1(std::vector<double> &j_x,
				      int ix_l, 
				      int ix_r,
				      double x_l, 
				      double x_r,
				      double length,
				      double charge,
				      int n_g,
				      int right_max)
{
  double midpoint, charge_fraction, weight;
  int bound_left, bound_right;

  bound_left = ix_l;
  bound_right = ix_r + 1;    

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    midpoint = (x_l + x_r) / 2.0;
    charge_fraction = 1.0;
    weight = 1.0 - midpoint;
    j_x[mod(bound_left,n_g)] += charge * charge_fraction * weight;
    for (int cell = (bound_left+1); cell < right_max; cell++) {
      j_x[mod(cell,n_g)] += charge * charge_fraction;
    }
  }
  else {
    // Left end
    midpoint = (x_l+1.0) / 2.0;
    charge_fraction = (1.0 - x_l) / length;
    weight = (1.0 - midpoint);
    j_x[mod(bound_left,n_g)] += charge * charge_fraction * weight;
    for (int cell = (bound_left+1); cell < right_max; cell++) {
      j_x[mod(cell,n_g)] += charge * charge_fraction;
    }

    // Pieces connecting two gridpoints
    charge_fraction = 1.0 / length;

    weight = 0.5;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      j_x[mod(cell,n_g)] += charge * charge_fraction * weight;
      for (int boundary = (cell+1); boundary < right_max; boundary++) {
	j_x[mod(boundary,n_g)] += charge * charge_fraction;
      }
    }
    
    // Right end
    midpoint = x_r / 2.0;
    charge_fraction = x_r / length;
    weight = 1.0 - midpoint;
    j_x[mod((bound_right-1),n_g)] += charge * charge_fraction * weight;
    for (int cell = (bound_right); cell < right_max; cell++) {
      j_x[mod(cell,n_g)] += charge * charge_fraction;
    }
  }
  return;
}

void deposit_rho_segment_2(std::vector<double> &rho,
			   int ixl,
			   int ixr,
			   double xl,
			   double xr,
			   double length,
			   double charge,
			   int n_g)
{
  double xa, xb, delta, q_norm, delta_left, delta_right;
  int ngp_left, ngp_right, shift;

  q_norm = charge / length;

  find_ngp(ixl, xl, ngp_left, delta_left);
  find_ngp(ixr, xr, ngp_right, delta_right);  

  // If tracers are in one cell
  if (ngp_left == ngp_right) {
    if (length==0.0) {
      delta = delta_left;

      rho[mod((ngp_left-1),n_g)] += charge * 0.5 * pow((0.5 - delta), 2);
      rho[mod(ngp_left,n_g)] += charge * (0.75 - delta * delta);
      rho[mod((ngp_left+1),n_g)] += charge * 0.5 * pow((0.5 + delta), 2);

    } else { 
      xa = delta_left;
      xb = delta_right;

      rho[mod((ngp_left-1),n_g)] += charge * (1.0 / 24.0) * (3.0 + 4.0 * xa * xa - 6.0 * xb + 4.0 * xb * xb + xa * (-6.0 + 4.0 *  xb));
      rho[mod(ngp_left,n_g)] += charge * (1.0 / 12.0) * (9.0 - 4.0 * (xa * xa + xa * xb + xb * xb));
      rho[mod((ngp_left+1),n_g)] += charge * (1.0 / 24.0) * (3.0 + 4.0 * xa * xa + 6.0 * xb + 4.0 *  xb*xb + xa * (6.0 + 4.0 * xb));
    }
  }
  else {
    // Left end
    xa = delta_left;
    xb = 0.5;

    rho[mod((ngp_left-1),n_g)] += q_norm * (1.0/24.0)*(xa*(-3.0+6.0*xa-4.0*xa*xa)+xb*(3.0-6.0*xb+4.0*xb*xb));
    rho[mod(ngp_left,n_g)] += q_norm * (1.0/12.0)*(-9.0*xa+4.0*pow(xa,3)+9.0*xb-4.0*pow(xb,3));
    rho[mod((ngp_left+1),n_g)] += q_norm * (1.0/24.0)*(-1.0*xa*(3.0+6.0*xa+4.0*xa*xa)+xb*(3.0+6.0*xb+4.0*xb*xb));
    
    // Portions fully covering a cell
    for (int cell = (ngp_left+1); cell < ngp_right; cell++) {
      rho[mod((cell-1),n_g)] += q_norm * (1.0 / 6.0);
      rho[mod(cell,n_g)] += q_norm * (2.0 / 3.0);
      rho[mod((cell+1),n_g)] += q_norm * (1.0 / 6.0);
    }
    
    // Right end
    xa = -0.5;
    xb = delta_right;

    rho[mod((ngp_right-1),n_g)] += q_norm * (1.0/24.0)*(xa*(-3.0+6.0*xa-4.0*xa*xa)+xb*(3.0-6.0*xb+4.0*xb*xb));
    rho[mod(ngp_right,n_g)] += q_norm * (1.0/12.0)*(-9.0*xa+4.0*pow(xa,3)+9.0*xb-4.0*pow(xb,3));
    rho[mod((ngp_right+1),n_g)] += q_norm * (1.0/24.0)*(-1.0*xa*(3.0+6.0*xa+4.0*xa*xa)+xb*(3.0+6.0*xb+4.0*xb*xb));

  }
  return;
}

void deposit_charge_to_left_segment_2(std::vector<double> &j_x,
				      int ixl,
				      int ixr,
				      double xl,
				      double xr,
				      double length,
				      double charge,
				      int n_g,
				      int right_max)
{
  double xa, xb, delta, q_norm, j_run, delta_left, delta_right;
  int ngp_left, ngp_right;

  q_norm = charge / length;

  find_ngp(ixl, xl, ngp_left, delta_left);
  find_ngp(ixr, xr, ngp_right, delta_right);

  // If tracers are in one cell
  if (ngp_left == ngp_right) {
    if (length==0.0) {
      delta = delta_left;

      j_run = charge * 0.5 * pow((0.5 - delta), 2);
      j_x[mod((ngp_left-1),n_g)] += j_run;

      j_run += charge * (0.75 - delta * delta);
      j_x[mod(ngp_left,n_g)] += j_run;

      j_run += charge * 0.5 * pow((0.5 + delta), 2);
      j_x[mod((ngp_left+1),n_g)] += j_run;

      for (int i = ngp_left+2; i < right_max; i++) {
	j_x[mod(i,n_g)] += j_run;
      }
      
    } else { 
      xa = delta_left;
      xb = delta_right;

      j_run = charge * (1.0 / 24.0) * (3.0 + 4.0 * xa * xa - 6.0 * xb + 4.0 * xb * xb + xa * (-6.0 + 4.0 *  xb));      
      j_x[mod((ngp_left-1),n_g)] += j_run;

      j_run += charge * (1.0 / 12.0) * (9.0 - 4.0 * (xa * xa + xa * xb + xb * xb));      
      j_x[mod(ngp_left,n_g)] += j_run;

      j_run += charge * (1.0 / 24.0) * (3.0 + 4.0 * xa * xa + 6.0 *  xb + 4.0 *  xb*xb + xa * (6.0 + 4.0 * xb));
      j_x[mod((ngp_left+1),n_g)] += j_run;

      for (int i = ngp_left+2; i < right_max; i++) {
	j_x[mod(i,n_g)] += j_run;
      }

    }
  }
  else {
    // Left end
    xa = delta_left;
    xb = 0.5;
    j_run = q_norm * (1.0/24.0)*(xa*(-3.0+6.0*xa-4.0*xa*xa)+xb*(3.0-6.0*xb+4.0*xb*xb));
    j_x[mod((ngp_left-1),n_g)] += j_run;

    j_run += q_norm * (1.0/12.0)*(-9.0*xa+4.0*pow(xa,3)+9.0*xb-4.0*pow(xb,3));
    j_x[mod(ngp_left,n_g)] += j_run;

    j_run += q_norm * (1.0/24.0)*(-1.0*xa*(3.0+6.0*xa+4.0*xa*xa)+xb*(3.0+6.0*xb+4.0*xb*xb));
    j_x[mod((ngp_left+1),n_g)] += j_run;

    for (int i = ngp_left+2; i < right_max; i++) {
      j_x[mod(i,n_g)] += j_run;
    }

    // Portions fully covering a cell
    for (int cell = (ngp_left+1); cell < ngp_right; cell++) {
      j_run = q_norm * (1.0 / 6.0);
      j_x[mod((cell-1),n_g)] += j_run;

      j_run += q_norm * (2.0 / 3.0);
      j_x[mod(cell,n_g)] += j_run;

      j_run += q_norm * (1.0 / 6.0);
      j_x[mod((cell+1),n_g)] += j_run;

      for (int i = cell+2; i < right_max; i++) {
	j_x[mod(i,n_g)] += j_run;
      }

    }
    
    // Right end
    xa = -0.5;
    xb = delta_right;
    j_run = q_norm * (1.0/24.0)*(xa*(-3.0+6.0*xa-4.0*xa*xa)+xb*(3.0-6.0*xb+4.0*xb*xb));
    j_x[mod((ngp_right-1),n_g)] += j_run;

    j_run += q_norm * (1.0/12.0)*(-9.0*xa+4.0*pow(xa,3)+9.0*xb-4.0*pow(xb,3));
    j_x[mod(ngp_right,n_g)] += j_run;
    
    j_run += q_norm * (1.0/24.0)*(-1.0*xa*(3.0+6.0*xa+4.0*xa*xa)+xb*(3.0+6.0*xb+4.0*xb*xb));
    j_x[mod((ngp_right+1),n_g)] += j_run;

    for (int i = ngp_right+2; i < right_max; i++) {
      j_x[mod(i,n_g)] += j_run;
    }

  }
  return;
}

void deposit_rho_segment_3(std::vector<double> &rho,
			   int ixl,
			   int ixr,
			   double xl,
			   double xr,
			   double length,
			   double charge,
			   int n_g)
{
  double xa, xb, delta, q_norm;
  int bound_left, bound_right;

  q_norm = charge / length;

  bound_left = ixl;
  bound_right = ixr + 1;

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    if (length==0.0) {
      delta = xl - 0.5;
      rho[mod((bound_left-1),n_g)] += charge * (-1.0 * pow((-0.5 + delta), 3) / 6.0);
      rho[mod(bound_left,n_g)] += charge * (4.0 - 6.0 * pow((0.5 + delta), 2) + 3.0 * pow((0.5 + delta), 3)) / 6.0;
      rho[mod((bound_left+1),n_g)] += charge * (23.0 + 30.0*delta - 12.0*pow(delta, 2) - 24.0*pow(delta,3)) / 48.0;
      rho[mod((bound_left+2),n_g)] += charge * pow((0.5 + delta), 3) / 6.0;
    } else { 
      xa = xl - 0.5;
      xb = xr - 0.5;

      rho[mod((bound_left-1),n_g)] += q_norm * (1.0 / 24.0) * (pow((-0.5 + xa), 4) - pow((-0.5 + xb), 4));
      rho[mod(bound_left,n_g)] += q_norm * (1.0/48.0)*(xa*(-23.0 + xa*(15.0 + 4.0*xa - 6.0*xa*xa)) +  xb*(23.0 + xb*(-15.0 - 4.0*xb + 6.0*xb*xb)));
      rho[mod((bound_left+1),n_g)] += q_norm *(1.0/48.0)*(xa*(-23.0+xa*(-15.0+4.0*xa+6.0*xa*xa)) +xb*(23.0+xb*(15.0-2.0*xb*(2.0+3*xb))));
      rho[mod((bound_left+2),n_g)] += q_norm * (1.0 / 24.0) * (-1.0 * pow((0.5 + xa), 4) + pow((0.5 + xb), 4));
    }
  }
  else {
    // Left end
    xa = xl - 0.5;
    xb = 0.5;

    rho[mod((bound_left-1),n_g)] += q_norm * (1.0 / 24.0) * (pow((-0.5 + xa), 4) - pow((-0.5 + xb), 4));
    rho[mod(bound_left,n_g)] += q_norm * (1.0/48.0)*(xa*(-23.0 + xa*(15.0 + 4.0*xa - 6.0*xa*xa)) +  xb*(23.0 + xb*(-15.0 - 4.0*xb + 6.0*xb*xb)));
    rho[mod((bound_left+1),n_g)] += q_norm *(1.0/48.0)*(xa*(-23.0+xa*(-15.0+4.0*xa+6.0*xa*xa)) +xb*(23.0+xb*(15.0-2.0*xb*(2.0+3*xb))));
    rho[mod((bound_left+2),n_g)] += q_norm * (1.0 / 24.0) * (-1.0 * pow((0.5 + xa), 4) + pow((0.5 + xb), 4));
    
    // Portions connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {

      rho[mod((cell-1),n_g)] += q_norm * (1.0 / 24.0);
      rho[mod(cell,n_g)] += q_norm * (11.0 / 24.0);
      rho[mod((cell+1),n_g)] += q_norm * (11.0 / 24.0);
      rho[mod((cell+2),n_g)] += q_norm * (1.0 / 24.0);
    }
    
    // Right end
    xa = -0.5;
    xb = xr - 0.5;

    rho[mod((bound_right-2),n_g)] += q_norm * (1.0 / 24.0) * (pow((-0.5 + xa), 4) - pow((-0.5 + xb), 4));
    rho[mod((bound_right-1),n_g)] += q_norm * (1.0/48.0)*(xa*(-23.0 + xa*(15.0 + 4.0*xa - 6.0*xa*xa)) +  xb*(23.0 + xb*(-15.0 - 4.0*xb + 6.0*xb*xb)));
    rho[mod(bound_right,n_g)] += q_norm *(1.0/48.0)*(xa*(-23.0+xa*(-15.0+4.0*xa+6.0*xa*xa)) +xb*(23.0+xb*(15.0-2.0*xb*(2.0+3*xb))));
    rho[mod((bound_right+1),n_g)] += q_norm * (1.0 / 24.0) * (-1.0 * pow((0.5 + xa), 4) + pow((0.5 + xb), 4));

  }
  return;
}

void deposit_charge_to_left_segment_3(std::vector<double> &j_x,
				      int ixl,
				      int ixr,
				      double xl,
				      double xr,
				      double length,
				      double charge,
				      int n_g,
				      int right_max)
{
  double xa, xb, delta, q_norm, j_run;
  int bound_left, bound_right;

  q_norm = charge / length;

  bound_left = ixl;
  bound_right = ixr + 1;

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    if (length==0.0) {
      delta = xl - 0.5;
      
      j_run = charge * (-1.0 * pow((-0.5 + delta), 3) / 6.0);
      j_x[mod((bound_left-1),n_g)] += j_run;
      
      j_run += charge * (4.0 - 6.0 * pow((0.5 + delta), 2) + 3.0 * pow((0.5 + delta), 3)) / 6.0;
      j_x[mod(bound_left,n_g)] += j_run;
	
      j_run += charge * (23.0 + 30.0*delta - 12.0*pow(delta, 2) - 24.0*pow(delta,3)) / 48.0;
      j_x[mod((bound_left+1),n_g)] += j_run;
	
      j_run += charge * pow((0.5 + delta), 3) / 6.0;
      j_x[mod((bound_left+2),n_g)] += j_run;

      for (int i = bound_left+3; i < right_max; i++) {
	j_x[mod(i,n_g)] += j_run;
      }
	
    } else { 
      xa = xl - 0.5;
      xb = xr - 0.5;

      j_run = q_norm * (1.0 / 24.0) * (pow((-0.5 + xa), 4) - pow((-0.5 + xb), 4));
      j_x[mod((bound_left-1),n_g)] += j_run;
	
      j_run += q_norm * (1.0/48.0)*(xa*(-23.0 + xa*(15.0 + 4.0*xa - 6.0*xa*xa)) + xb*(23.0 + xb*(-15.0 - 4.0*xb + 6.0*xb*xb)));
      j_x[mod(bound_left,n_g)] += j_run;
	
      j_run += q_norm *(1.0/48.0)*(xa*(-23.0+xa*(-15.0+4.0*xa+6.0*xa*xa)) +xb*(23.0+xb*(15.0-2.0*xb*(2.0+3*xb))));
      j_x[mod((bound_left+1),n_g)] += j_run;
      
      j_run += q_norm * (1.0 / 24.0) * (-1.0 * pow((0.5 + xa), 4) + pow((0.5 + xb), 4));
      j_x[mod((bound_left+2),n_g)] += j_run;
      for (int i = bound_left+3; i < right_max; i++) {
	j_x[mod(i,n_g)] += j_run;
      }

    }
  }
  else {
    // Left end
    xa = xl - 0.5;
    xb = 0.5;
    
    j_run = q_norm * (1.0 / 24.0) * (pow((-0.5 + xa), 4) - pow((-0.5 + xb), 4));
    j_x[mod((bound_left-1),n_g)] += j_run;

    j_run += q_norm * (1.0/48.0)*(xa*(-23.0 + xa*(15.0 + 4.0*xa - 6.0*xa*xa)) + xb*(23.0 + xb*(-15.0 - 4.0*xb + 6.0*xb*xb)));
    j_x[mod(bound_left,n_g)] += j_run;

    j_run += q_norm * (1.0/48.0)*(xa*(-23.0+xa*(-15.0+4.0*xa+6.0*xa*xa))+xb*(23.0+xb*(15.0-2.0*xb*(2.0+3*xb))));
    j_x[mod((bound_left+1),n_g)] += j_run;
    
    j_run += q_norm * (1.0 / 24.0) * (-1.0 * pow((0.5 + xa), 4) + pow((0.5 + xb), 4));
    j_x[mod((bound_left+2),n_g)] += j_run;

    for (int i = bound_left+3; i < right_max; i++) {
      j_x[mod(i,n_g)] += j_run;
    }
    
    // Portions connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      
      j_run = q_norm * (1.0 / 24.0);
      j_x[mod((cell-1),n_g)] += j_run;
	
      j_run += q_norm * (11.0 / 24.0);
      j_x[mod(cell,n_g)] += j_run;
	
      j_run += q_norm * (11.0 / 24.0);
      j_x[mod((cell+1),n_g)] += j_run;
	
      j_run += q_norm * (1.0 / 24.0);
      j_x[mod((cell+2),n_g)] += j_run;

      for (int i = cell+3; i < right_max; i++) {
	j_x[mod(i,n_g)] += j_run;
      }

    }
    
    // Right end
    xa = -0.5;
    xb = xr - 0.5;
    
    j_run = q_norm * (1.0 / 24.0) * (pow((-0.5 + xa), 4) - pow((-0.5 + xb), 4));
    j_x[mod((bound_right-2),n_g)] += j_run;

    j_run += q_norm * (1.0/48.0)*(xa*(-23.0 + xa*(15.0 + 4.0*xa - 6.0*xa*xa)) +  xb*(23.0 + xb*(-15.0 - 4.0*xb + 6.0*xb*xb)));
    j_x[mod((bound_right-1),n_g)] += j_run;
    
    j_run += q_norm *(1.0/48.0)*(xa*(-23.0+xa*(-15.0+4.0*xa+6.0*xa*xa)) +xb*(23.0+xb*(15.0-2.0*xb*(2.0+3*xb))));
    j_x[mod(bound_right,n_g)] += j_run;
      
    j_run += q_norm * (1.0 / 24.0) * (-1.0 * pow((0.5 + xa), 4) + pow((0.5 + xb), 4));
    j_x[mod((bound_right+1),n_g)] += j_run;

    for (int i = bound_right+2; i < right_max; i++) {
      j_x[mod(i,n_g)] += j_run;
    }
    
  }
  return;
}

void deposit_rho_segment_4(std::vector<double> &rho,
			   int ixl,
			   int ixr,
			   double xl,
			   double xr,
			   double length,
			   double charge,
			   int n_g)
{
  double xa, xb, delta, q_norm, delta_left, delta_right;
  int ngp_left, ngp_right;

  q_norm = charge / length;

  find_ngp(ixl, xl, ngp_left, delta_left);
  find_ngp(ixr, xr, ngp_right, delta_right);

  // If tracers are in one cell
  if (ngp_left == ngp_right) {
    if (length==0.0) {
      delta = delta_left;

      rho[mod((ngp_left-2),n_g)] += charge * pow((1.0 - 2.0*delta), 4) / 384.0;
      rho[mod((ngp_left-1),n_g)] += charge * (19.0 - 44.0*delta + 24.0*pow(delta, 2) + 16.0*pow(delta,3) - 16.0*pow(delta, 4))/96.0;
      rho[mod(ngp_left,n_g)] += charge * (115.0 / 192.0) - (5.0*pow(delta, 2))/8.0 + pow(delta, 4)/4.0;
      rho[mod((ngp_left+1),n_g)] += charge * (19.0 + 44.0*delta + 24.0*pow(delta, 2) - 16.0*pow(delta, 3) - 16*pow(delta,4))/96.0;
      rho[mod((ngp_left+2),n_g)] += charge * pow((1.0 + 2.0*delta), 4) / 384.0;
	
    } else { 
      xa = delta_left;
      xb = delta_right;

      rho[mod((ngp_left-2),n_g)] += q_norm * (-1.0*pow((-1.0+2.0*xa),5) + pow((-1.0+2.0*xb),5)) / 3840.0;
      rho[mod((ngp_left-1),n_g)] += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(55.0+2.0*xa*(-10.0+xa*(-5.0+4.0*xa))))+xb*(95.0+2.0*xb*(-55.0+2.0*xb*(10.0+(5.0-4.0*xb)*xb))));
      rho[mod(ngp_left,n_g)] += q_norm * (1.0/960.0)*(-575.0*xa+200.0*pow(xa,3)-48.0*pow(xa,5)+575.0*xb-200.0*pow(xb,3)+48.0*pow(xb,5));
      rho[mod((ngp_left+1),n_g)] += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(-55.0+2.0*xa*(-10.0+xa*(5.0+4.0*xa))))+xb*(95.0+2.0*xb*(55.0-2.0*xb*(-10.0+xb*(5.0+4.0*xb)))));
      rho[mod((ngp_left+2),n_g)] += q_norm * (-1.0*pow((1.0+2.0*xa),5) + pow((1.0+2.0*xb),5)) / 3840.0;
    }
  }
  else {
    // Left end
    xa = delta_left;
    xb = 0.5;

    rho[mod((ngp_left-2),n_g)] += q_norm * (-1.0*pow((-1.0+2.0*xa),5) + pow((-1.0+2.0*xb),5)) / 3840.0;
    rho[mod((ngp_left-1),n_g)] += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(55.0+2.0*xa*(-10.0+xa*(-5.0+4.0*xa))))+xb*(95.0+2.0*xb*(-55.0+2.0*xb*(10.0+(5.0-4.0*xb)*xb))));
    rho[mod(ngp_left,n_g)] += q_norm * (1.0/960.0)*(-575.0*xa+200.0*pow(xa,3)-48.0*pow(xa,5)+575.0*xb-200.0*pow(xb,3)+48.0*pow(xb,5));
    rho[mod((ngp_left+1),n_g)] += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(-55.0+2.0*xa*(-10.0+xa*(5.0+4.0*xa))))+xb*(95.0+2.0*xb*(55.0-2.0*xb*(-10.0+xb*(5.0+4.0*xb)))));
    rho[mod((ngp_left+2),n_g)] += q_norm * (-1.0*pow((1.0+2.0*xa),5) + pow((1.0+2.0*xb),5)) / 3840.0;
    
    // Portions fully covering a cell
    for (int cell = (ngp_left+1); cell < ngp_right; cell++) {
      rho[mod((cell-2),n_g)] += q_norm * (1.0 / 120.0);
      rho[mod((cell-1),n_g)] += q_norm * (13.0 / 60.0);
      rho[mod(cell,n_g)] += q_norm * (11.0 / 20.0);
      rho[mod((cell+1),n_g)] += q_norm * (13.0 / 60.0);
      rho[mod((cell+2),n_g)] += q_norm * (1.0 / 120.0);
    }
    
    // Right end
    xa = -0.5;
    xb = delta_right;

    rho[mod((ngp_right-2),n_g)] += q_norm * (-1.0*pow((-1.0+2.0*xa),5) + pow((-1.0+2.0*xb),5)) / 3840.0;
    rho[mod((ngp_right-1),n_g)] += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(55.0+2.0*xa*(-10.0+xa*(-5.0+4.0*xa))))+xb*(95.0+2.0*xb*(-55.0+2.0*xb*(10.0+(5.0-4.0*xb)*xb))));
    rho[mod(ngp_right,n_g)] += q_norm * (1.0/960.0)*(-575.0*xa+200.0*pow(xa,3)-48.0*pow(xa,5)+575.0*xb-200.0*pow(xb,3)+48.0*pow(xb,5));
    rho[mod((ngp_right+1),n_g)] += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(-55.0+2.0*xa*(-10.0+xa*(5.0+4.0*xa))))+xb*(95.0+2.0*xb*(55.0-2.0*xb*(-10.0+xb*(5.0+4.0*xb)))));
    rho[mod((ngp_right+2),n_g)] += q_norm * (-1.0*pow((1.0+2.0*xa),5) + pow((1.0+2.0*xb),5)) / 3840.0;
  }
  return;
}

void deposit_charge_to_left_segment_4(std::vector<double> &j_x,
				      int ixl,
				      int ixr,
				      double xl,
				      double xr,
				      double length,
				      double charge,
				      int n_g,
				      int right_max)
{
  double xa, xb, delta, q_norm, j_run, delta_left, delta_right;
  int ngp_left, ngp_right;

  q_norm = charge / length;

  find_ngp(ixl, xl, ngp_left, delta_left);
  find_ngp(ixr, xr, ngp_right, delta_right);  

  // If tracers are in one cell
  if (ngp_left == ngp_right) {
    if (length==0.0) {
      delta = delta_left;

      j_run = charge * pow((1.0 - 2.0*delta), 4) / 384.0;
      j_x[mod((ngp_left-2),n_g)] += j_run;
      
      j_run += charge * (19.0 - 44.0*delta + 24.0*pow(delta, 2) + 16.0*pow(delta,3) - 16.0*pow(delta, 4))/96.0;
      j_x[mod((ngp_left-1),n_g)] += j_run;
      
      j_run += charge * (115.0 / 192.0) - (5.0*pow(delta, 2))/8.0 + pow(delta, 4)/4.0;
      j_x[mod(ngp_left,n_g)] += j_run;
      
      j_run += charge * (19.0 + 44.0*delta + 24.0*pow(delta, 2) - 16.0*pow(delta, 3) - 16*pow(delta,4))/96.0;
      j_x[mod((ngp_left+1),n_g)] += j_run;
      
      j_run += charge * pow((1.0 + 2.0*delta), 4) / 384.0;
      j_x[mod((ngp_left+2),n_g)] += j_run;
      
      for (int i = ngp_left+3; i < right_max; i++) {
	j_x[mod(i,n_g)] += j_run;
      }
	
    } else { 
      xa = delta_left;
      xb = delta_right;

      j_run = q_norm * (-1.0*pow((-1.0+2.0*xa),5) + pow((-1.0+2.0*xb),5)) / 3840.0;
      j_x[mod((ngp_left-2),n_g)] += j_run;
      
      j_run += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(55.0+2.0*xa*(-10.0+xa*(-5.0+4.0*xa))))+xb*(95.0+2.0*xb*(-55.0+2.0*xb*(10.0+(5.0-4.0*xb)*xb))));
      j_x[mod((ngp_left-1),n_g)] += j_run;
      
      j_run += q_norm * (1.0/960.0)*(-575.0*xa+200.0*pow(xa,3)-48.0*pow(xa,5)+575.0*xb-200.0*pow(xb,3)+48.0*pow(xb,5));
      j_x[mod(ngp_left,n_g)] += j_run;
      
      j_run += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(-55.0+2.0*xa*(-10.0+xa*(5.0+4.0*xa))))+xb*(95.0+2.0*xb*(55.0-2.0*xb*(-10.0+xb*(5.0+4.0*xb)))));
      j_x[mod((ngp_left+1),n_g)] += j_run;
      
      j_run += q_norm * (-1.0*pow((1.0+2.0*xa),5) + pow((1.0+2.0*xb),5)) / 3840.0;
      j_x[mod((ngp_left+2),n_g)] += j_run;
      
      for (int i = ngp_left+3; i < right_max; i++) {
	j_x[mod(i,n_g)] += j_run;
      }
      
    }
  }
  else {
    // Left end
    xa = delta_left;
    xb = 0.5;

    j_run = q_norm * (-1.0*pow((-1.0+2.0*xa),5) + pow((-1.0+2.0*xb),5)) / 3840.0;
    j_x[mod((ngp_left-2),n_g)] += j_run;
    
    j_run += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(55.0+2.0*xa*(-10.0+xa*(-5.0+4.0*xa))))+xb*(95.0+2.0*xb*(-55.0+2.0*xb*(10.0+(5.0-4.0*xb)*xb))));
    j_x[mod((ngp_left-1),n_g)] += j_run;
    
    j_run += q_norm * (1.0/960.0)*(-575.0*xa+200.0*pow(xa,3)-48.0*pow(xa,5)+575.0*xb-200.0*pow(xb,3)+48.0*pow(xb,5));
    j_x[mod(ngp_left,n_g)] += j_run;
    
    j_run += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(-55.0+2.0*xa*(-10.0+xa*(5.0+4.0*xa))))+xb*(95.0+2.0*xb*(55.0-2.0*xb*(-10.0+xb*(5.0+4.0*xb)))));
    j_x[mod((ngp_left+1),n_g)] += j_run;
    
    j_run += q_norm * (-1.0*pow((1.0+2.0*xa),5) + pow((1.0+2.0*xb),5)) / 3840.0;
    j_x[mod((ngp_left+2),n_g)] += j_run;
    
    for (int i = ngp_left+3; i < right_max; i++) {
      j_x[mod(i,n_g)] += j_run;
    }
    
    // Portions fully covering a cell
    for (int cell = (ngp_left+1); cell < ngp_right; cell++) {
      j_run = q_norm * (1.0 / 120.0);
      j_x[mod((cell-2),n_g)] += j_run;
      
      j_run += q_norm * (13.0 / 60.0);
      j_x[mod((cell-1),n_g)] += j_run;
      
      j_run += q_norm * (11.0 / 20.0);
      j_x[mod(cell,n_g)] += j_run;
      
      j_run += q_norm * (13.0 / 60.0);
      j_x[mod((cell+1),n_g)] += j_run;
      
      j_run += q_norm * (1.0 / 120.0);
      j_x[mod((cell+2),n_g)] += j_run;
      
      for (int i = cell+3; i < right_max; i++) {
	j_x[mod(i,n_g)] += j_run;
      }

    }
    
    // Right end
    xa = -0.5;
    xb = delta_right;

    j_run = q_norm * (-1.0*pow((-1.0+2.0*xa),5) + pow((-1.0+2.0*xb),5)) / 3840.0;
    j_x[mod((ngp_right-2),n_g)] += j_run;

    j_run += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(55.0+2.0*xa*(-10.0+xa*(-5.0+4.0*xa))))+xb*(95.0+2.0*xb*(-55.0+2.0*xb*(10.0+(5.0-4.0*xb)*xb))));
    j_x[mod((ngp_right-1),n_g)] += j_run;
    
    j_run += q_norm * (1.0/960.0)*(-575.0*xa+200.0*pow(xa,3)-48.0*pow(xa,5)+575.0*xb-200.0*pow(xb,3)+48.0*pow(xb,5));
    j_x[mod(ngp_right,n_g)] += j_run;
    
    j_run += q_norm * (1.0/480.0)*(xa*(-95.0+2.0*xa*(-55.0+2.0*xa*(-10.0+xa*(5.0+4.0*xa))))+xb*(95.0+2.0*xb*(55.0-2.0*xb*(-10.0+xb*(5.0+4.0*xb)))));
    j_x[mod((ngp_right+1),n_g)] += j_run;

    j_run += q_norm * (-1.0*pow((1.0+2.0*xa),5) + pow((1.0+2.0*xb),5)) / 3840.0;
    j_x[mod((ngp_right+2),n_g)] += j_run;
    
    for (int i = ngp_right+3; i < right_max; i++) {
      j_x[mod(i,n_g)] += j_run;
    }
    
  }
  return;
}

int find_max(int ix_a, int ix_b, int ix_c, int ix_d)
{
  int max;
  max = ix_a;
  if (ix_b > max) {
    max = ix_b;
  }
  if (ix_c > max) {
    max = ix_c;
  }
  if (ix_d > max) {
    max = ix_d;
  }
  return max;
}

void ParticleSpecies::deposit_j_x_sic(std::vector<double> &j_x)
{
  int ix_l, ix_r, ix_l_old, ix_r_old, right_max;
  double x_l, x_r, x_l_old, x_r_old, length, length_old, q;

  for (int i = 0; i < n_p; i++) {
    q = charge[i];
    arrange_segment(ix_old[i], x_old[i],
		    ix_old[i+1], x_old[i+1],
		    ix_l_old, x_l_old,
		    ix_r_old, x_r_old,
		    length_old);
    arrange_segment(ix[i], x[i],
		    ix[i+1], x[i+1],
		    ix_l, x_l,
		    ix_r, x_r,
		    length);
    // Index of right bounding gridpoint (not really(?), more precisely define this)
    // SOME RIGHT MAXES WERE NGP!!!! NEED TO UNDERSTAND EVERYTHING CAREFULLY    
    right_max = find_max(ix[i], ix_old[i], ix[i+1], ix_old[i+1]);
    switch (method) {
    case 5 :
      right_max = right_max + 1;
      deposit_charge_to_left_segment_0(j_x, ix_l_old, ix_r_old, x_l_old, x_r_old, length_old, q, n_g, right_max);
      deposit_charge_to_left_segment_0(j_x, ix_l, ix_r, x_l, x_r, length, (-1.0)*q, n_g, right_max);
      break;
    case 6 :
      right_max = right_max + 1;
      deposit_charge_to_left_segment_1(j_x, ix_l_old, ix_r_old, x_l_old, x_r_old, length_old, q, n_g, right_max);
      deposit_charge_to_left_segment_1(j_x, ix_l, ix_r, x_l, x_r, length, (-1.0)*q, n_g, right_max);
      break;
    case 7 :      
      //      right_max = find_max(ix[i], ix_old[i], ix[i+1], ix_old[i+1]);
      // right_max = fmax(fmax(x_old[i], x_old[i+1]),fmax(x[i],x[i+1]));
      // right_max = get_nearest_gridpoint(right_max) + 2;
      right_max = right_max + 3; // 3 works, 2 doesn't.  But should check NGP above ^
      
      deposit_charge_to_left_segment_2(j_x, ix_l_old, ix_r_old, x_l_old, x_r_old, length_old, q, n_g, right_max);
      deposit_charge_to_left_segment_2(j_x, ix_l, ix_r, x_l, x_r, length, (-1.0)*q, n_g, right_max);
      break;
    case 8 :      
      right_max = right_max + 3;
      deposit_charge_to_left_segment_3(j_x, ix_l_old, ix_r_old, x_l_old, x_r_old, length_old, q, n_g, right_max);
      deposit_charge_to_left_segment_3(j_x, ix_l, ix_r, x_l, x_r, length, (-1.0)*q, n_g, right_max);
      break;
    case 9 :      
      right_max = right_max + 4;
      deposit_charge_to_left_segment_4(j_x, ix_l_old, ix_r_old, x_l_old, x_r_old, length_old, q, n_g, right_max);
      deposit_charge_to_left_segment_4(j_x, ix_l, ix_r, x_l, x_r, length, (-1.0)*q, n_g, right_max);
      break;
    default:
      std::cout << "Error, selected interpolation order not implemented." << std::endl;
    }
  }
  return;
}

void ParticleSpecies::deposit_rho_sic(std::vector<double> &rho)
{
  int ix_l, ix_r;
  double x_l, x_r, length;
  for (int i = 0; i < n_p; i++) {
    arrange_segment(ix[i], x[i], ix[i+1], x[i+1], ix_l, x_l, ix_r, x_r, length);
    switch (method) {
    case 5 :
      deposit_rho_segment_0(rho, ix_l, ix_r, x_l, x_r, length, charge[i], n_g);
      break;
    case 6 :
      deposit_rho_segment_1(rho, ix_l, ix_r, x_l, x_r, length, charge[i], n_g);
      break;
    case 7 :
      deposit_rho_segment_2(rho, ix_l, ix_r, x_l, x_r, length, charge[i], n_g);
      break;
    case 8 :
      deposit_rho_segment_3(rho, ix_l, ix_r, x_l, x_r, length, charge[i], n_g);
      break;
    case 9 :
      deposit_rho_segment_4(rho, ix_l, ix_r, x_l, x_r, length, charge[i], n_g);
      break;
    default:
      std::cout << "Error, selected interpolation order not implemented." << std::endl;
    }
  }
  return;
}

double interpolate_segment_velocity(double v_left, double v_right, 
				    double length, double distance_from_left)
{
  return v_left + (distance_from_left / length) * (v_right - v_left);
}

void deposit_j_t_segment_0(std::vector<double> &j_t,
			   int ixl,
			   int ixr,
			   double xl, 
			   double xr,
			   double vl,
			   double vr,
			   double length,			     
			   double charge,
			   int n_g)
{
  double charge_fraction, avg_velocity, midpoint, distance_from_left;
  int bound_left, bound_right;

  // Shift tracer positions so that 0 is the left bounary of cell zero
  xl = xl + 0.5;
  if (xl >= 1.0) {
    ixl = ixl + 1;
    xl = xl - 1.0;
  }
  
  xr = xr + 0.5;
  if (xr >= 1.0) {
    ixr = ixr + 1;
    xr = xr - 1.0;
  }

  bound_left = ixl;
  bound_right = ixr + 1;

  if (bound_right == (bound_left + 1)) {
    // If tracers are between two gridpoints
    charge_fraction = 1.0;
    avg_velocity = (vl + vr) / 2.0;
    j_t[mod(bound_left,n_g)] += charge * charge_fraction * avg_velocity;
  }
  else {
    // Left end
    charge_fraction = (1.0 - xl) / length;
    midpoint = (1.0 + xl) / 2.0;
    distance_from_left = midpoint - xl;
    avg_velocity = interpolate_segment_velocity(vl, vr, length, distance_from_left);
    j_t[mod(bound_left,n_g)] += charge * charge_fraction * avg_velocity;
    // Pieces connecting two gridpoints
    charge_fraction = 1.0 / length;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      distance_from_left = (cell-ixl) + 0.5 - xl;
      avg_velocity = interpolate_segment_velocity(vl, vr, length, distance_from_left);
      j_t[mod(cell,n_g)] += charge * charge_fraction * avg_velocity;
    }
    // Right end
    charge_fraction = xr / length;
    midpoint = xr / 2.0;
    distance_from_left = (midpoint - xl) + (ixr - ixl);
    avg_velocity = interpolate_segment_velocity(vl, vr, length, distance_from_left);
    j_t[mod((bound_right-1),n_g)] += charge * charge_fraction * avg_velocity;
  }
  return;
}

// Precision of this and the following function should be improved by using relative
// values for xa and xb rather than the larger absolute values
// (MAYBE?) check this all carefully
double j_t_segment_linear_left_gridpoint(double charge, double xl, 
					 double length, double vl, double vr, 
					 double xa, double xb)
{
  double xip1 = ceil(xb);
  double j_t = (1.0/(6.0*length*length))*charge*(xa-xb);
  j_t = j_t * (3.0*length*vl*(xa+xb-2.0*xip1)-(vl-vr)*(2.0*xa*xa+2.0*xb*xb+6.0*xip1*xl-3.0*xb*(xip1+xl)+xa*(2.0*xb-3.0*(xip1+xl))));
  return j_t;
}

double j_t_segment_linear_right_gridpoint(double charge, double xl, 
					  double length, double vl, double vr, 
					  double xa, double xb)
{
  double xi = floor(xa);
  double j_t = (1.0/(6.0*length*length))*charge*(xa-xb);
  j_t = j_t * (-3.0*length*vl*(xa+xb-2.0*xi)+(vl-vr)*(2.0*xa*xa+2.0*xb*xb+6.0*xi*xl-3.0*xb*(xi+xl)+xa*(2.0*xb-3.0*(xi+xl))));
  return j_t;
}

void deposit_j_t_segment_1(std::vector<double> &j_t,
			   int ixl,
			   int ixr,
			   double xl, 
			   double xr,
			   double vl,
			   double vr,
			   double length,			     
			   double charge,
			   int n_g)
{
  double xa, xb;
  int bound_left, bound_right;

  bound_left = ixl;
  bound_right = ixr + 1;

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    if (length==0.0) {
      j_t[mod(bound_left,n_g)] += charge * ((vl+vr)/2.0) * (1.0-xl);
      j_t[mod((bound_left+1),n_g)] += charge * ((vl+vr)/2.0) * xl;
    } else { 
      xa = xl;
      xb = xr;
      j_t[mod(bound_left,n_g)] += j_t_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
      j_t[mod((bound_left+1),n_g)] += j_t_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
    }
  }
  else {
    // Left end
    xa = xl;
    xb = 1.0;
    j_t[mod(bound_left,n_g)] += j_t_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
    j_t[mod((bound_left+1),n_g)] += j_t_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
    // Portions connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      xa = double(cell - ixl);
      xb = double(cell + 1 - ixl);
      j_t[mod(cell,n_g)] += j_t_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
      j_t[mod((cell+1),n_g)] += j_t_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
    }
    // Right end
    xa = double(ixr - ixl);
    xb = xr + double(ixr - ixl);
    j_t[mod((bound_right-1),n_g)] += j_t_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
    j_t[mod(bound_right,n_g)] += j_t_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
  }
  return;
}

double j_t_3_w1(double xa, double xb, double va, double vb)
{
  return (1.0/480.0)*(xa-xb)*(vb*(-5.0+2.0*xa*(5.0+xa*(-5.0+2.0*xa))+20.0*xb+4.0*xa*(-5.0+2.0*xa)*xb+6.0*(-5.0+2.0*xa)*xb*xb+16.0*pow(xb,3))+va*(-5.0+16.0*pow(xa,3)+6.0*xa*xa*(-5.0+2.0*xb)+4.0*xa*(5.0+xb*(-5.0+2.0*xb))+2.0*xb*(5.0+xb*(-5.0+2.0*xb))));
}

double j_t_3_w1_full(double xa, double xb, double va, double vb)
{
  return (vb*(3.0+10.0*xa)-va*(3.0+10.0*xb))/(240.0*(xa-xb));
}

double j_t_3_w2(double xa, double xb, double va, double vb)
{
  return -1.0*(1.0/480.0)*(xa-xb)*(vb*(115.0+2.0*xa*(-5.0+2.0*xa)*(5.0+3.0*xa)-100.0*xb+4.0*xa*(-5.0+6.0*xa)*xb+6.0*(-5.0+6.0*xa)*xb*xb+48.0*pow(xb,3))+va*(115.0+48.0*pow(xa,3)+4.0*xa*(-5.0+2.0*xb)*(5.0+3.0*xb)+2.0*xb*(-5.0+2.0*xb)*(5.0+3.0*xb)+6.0*xa*xa*(-5.0+6.0*xb)));  
}

double j_t_3_w2_full(double xa, double xb, double va, double vb)
{
  return (11.0*(vb+10.0*vb*xa-va*(1.0+10.0*xb)))/(240.0*(xa-xb));
}

double j_t_3_w3(double xa, double xb, double va, double vb)
{
  return (1.0/480.0)*(xa-xb)*(va*(-115.0+48.0*pow(xa,3)+4.0*xa*(5.0+2.0*xb)*(-5.0+3.0*xb)+2.0*xb*(5.0+2.0*xb)*(-5.0+3.0*xb)+6*xa*xa*(5.0+6.0*xb))+vb*(2.0*xa*(5.0+2.0*xa)*(-5.0+3.0*xa)+4.0*xa*(5.0+6.0*xa)*xb+6.0*(5.0+6.0*xa)*xb*xb+48.0*pow(xb,3)-5.0*(23.0+20.0*xb)));
}

double j_t_3_w3_full(double xa, double xb, double va, double vb)
{
  return (11.0*(va+vb*(-1.0+10.0*xa)-10.0*va*xb))/(240.0*(xa-xb));
}

double j_t_3_w4(double xa, double xb, double va, double vb)
{
  return -(1.0/480.0)*(xa-xb)*(vb*(5.0+2.0*xa*(5.0+xa*(5.0+2.0*xa))+20.0*xb+4.0*xa*(5.0+2.0*xa)*xb+6.0*(5.0+2.0*xa)*xb*xb+16.0*pow(xb,3))+va*(5.0+16.0*pow(xa,3)+6.0*xa*xa*(5.0+2.0*xb)+4.0*xa*(5.0+xb*(5.0+2.0*xb))+2.0*xb*(5.0+xb*(5.0+2.0*xb))));
}


double j_t_3_w4_full(double xa, double xb, double va, double vb)
{
  return (vb*(-3.0+10.0*xa)+va*(3.0-10.0*xb))/(240.0*(xa-xb));
}


double j_t_2_w1(double xa, double xb, double va, double vb)
{
  return (-1.0/48.0)*(xa-xb)*(va*(3.0+6.0*xa*xa+4.0*xa*(-2.0+xb)+2.0*(-2.0+xb)*xb)+vb*(3.0+2.0*(-2.0+xa)*xa-8.0*xb+4.0*xa*xb+6.0*xb*xb));
}

double j_t_2_w1_full(double xa, double xb, double va, double vb)
{
  return (vb+4.0*vb*xa-va*(1.0+4.0*xb))/(24.0*(xa-xb));
}

double j_t_2_w2(double xa, double xb, double va, double vb)
{

  return (1.0/24.0)*(xa-xb)*(-9.0*va-9.0*vb+6.0*va*xa*xa+2.0*vb*xa*xa+4.0*(va+vb)*xa*xb+2.0*(va+3.0*vb)*xb*xb);
}

double j_t_2_w2_full(double xa, double xb, double va, double vb)
{
  return (2.0*vb*xa-2.0*va*xb)/(3.0*xa-3.0*xb);
}

double j_t_2_w3(double xa, double xb, double va, double vb)
{
  return (-1.0/48.0)*(xa-xb)*(vb*(3.0+2.0*xa*(2.0+xa)+8.0*xb+4.0*xa*xb+6.0*xb*xb)+va*(3.0+6.0*xa*xa+4.0*xa*(2.0+xb)+2.0*xb*(2.0 + xb)));
}

double j_t_2_w3_full(double xa, double xb, double va, double vb)
{
  return (va-vb+4.0*vb*xa-4.0*va*xb)/(24.0*xa-24.0*xb);
}

double j_t_4_w1(double xa, double xb, double va, double vb)
{
  return (-1.0/11520.0)*(xa-xb)*(vb*(15.0+4.0*xa*(-10.0+xa*(15.0+4.0*(-3.0+xa)*xa))-80.0*xb+8.0*xa*(15.0+4.0*(-3.0+xa)*xa)*xb+12.0*(15.0+4.0*(-3.0+xa)*xa)*xb*xb+64.0*(-3.0+xa)*pow(xb,3)+80.0*pow(xb,4))+va*(15.0+80.0*pow( xa,4)+64.0*pow( xa,3)*(-3.0+xb)+12.0*xa*xa*(15.0+4.0*(-3.0+xb)*xb)+8.0*xa*(-10.0+xb*(15.0+4.0*(-3.0+xb)*xb))+4.0*xb*(-10.0+xb*(15.0+4.0*(-3.0+xb)*xb))));
}

double j_t_4_w2(double xa, double xb, double va, double vb)
{

  return (1.0/2880.0)*(xa-xb)*(-285.0*vb+va*(5.0*(-57.0+44.0*xb)+4.0*(xa*(110.0+xa*(-45.0+4.0*xa*(-6.0+5.0*xa)))+2.0*xa*(-15.0+xa*(-9.0+8.0*xa))*xb+3.0*(-5.0+4.0*(-1.0+xa)*xa)*pow(xb,2)+2.0*(-3.0+4.0*xa)*pow(xb,3)+4.0*pow(xb,4)))+4.0*vb*(4.0*pow(xa,4)+pow(xa,3)*(-6.0+8.0*xb)+3.0*xa*xa*(-5.0+4.0*(-1.0+xb)*xb)+xb*(110.0+xb*(-45.0+4.0*xb*(-6.0+5.0*xb)))+xa*(55.0+2.0*xb*(-15.0+xb*(-9.0+8.0*xb)))));
}

double j_t_4_w3(double xa, double xb, double va, double vb)
{

  return (-1.0/1920.0)*(xa-xb)*(575.0*va+575.0*vb-300.0*va*xa*xa-100.0*vb*xa*xa+80.0*va*pow(xa,4)+16.0*vb*pow(xa,4)+8.0*xa*(-25.0*(va+vb)+4.0*(2.0*va+vb)*xa*xa)*xb+4.0*(-25.0*(va+3.0*vb)+12.0*(va+vb)*xa*xa)*xb*xb+32.0*(va+2.0*vb)*xa*pow(xb,3)+16.0*(va+5.0*vb)*pow(xb,4));
}

double j_t_4_w4(double xa, double xb, double va, double vb)
{

  return (1.0/2880.0)*(xa-xb)*(-5.0*va*(57.0+44.0*xb)-5.0*vb*(57.0+88.0*xb)+4.0*va*(xa*(-110.0+xa*(-45.0+4.0*xa*(6.0+5.0*xa)))+2.0*xa*(-15.0+xa*(9.0+8.0*xa))*xb+3.0*(-5.0+4.0*xa*(1.0+xa))*xb*xb+2.0*(3.0+4.0*xa)*pow(xb,3)+4.0*pow(xb,4))+4.0*vb*(4.0*pow(xa,4)+pow(xa,3)*(6.0+8.0*xb)+3.0*xa*xa*(-5.0+4.0*xb*(1.0+xb))+xb*xb*(-45.0+4.0*xb*(6.0+5.0*xb))+xa*(-55.0+2.0*xb*(-15.0+xb*(9.0+8.0*xb)))));
}

double j_t_4_w5(double xa, double xb, double va, double vb)
{
  return (-1.0/11520.0)*(xa-xb)*(vb*(15.0+4.0*xa*(10.0+xa*(15.0+4.0*xa*(3.0+xa)))+80*xb+8.0*xa*(15.0+4.0*xa*(3.0+xa))*xb+12.0*(15.0+4.0*xa*(3.0+xa))*xb*xb+64.0*(3.0+xa)*pow(xb,3)+80.0*pow(xb,4))+va*(15.0+80.0*pow(xa,4)+64.0*pow(xa,3)*(3.0+xb)+12.0*xa*xa*(15.0+4.0*xb*(3.0+xb))+8.0*xa*(10.0+xb*(15.0+4.0*xb*(3.0+xb)))+4.0*xb*(10.0+xb*(15.0+4.0*xb*(3.0+xb)))));
}

double j_t_4_w1_full(double xa, double xb, double va, double vb)
{
  return (vb+3.0*vb*xa-va*(1.0+3.0*xb))/(360.0*(xa-xb));
}

double j_t_4_w2_full(double xa, double xb, double va, double vb)
{
  return (13.0*(vb+6.0*vb*xa-va*(1.0+6.0*xb)))/(360.0*(xa-xb));
}

double j_t_4_w3_full(double xa, double xb, double va, double vb)
{
  return (11.0*(vb*xa-va*xb))/(20.0*(xa-xb));
}

double j_t_4_w4_full(double xa, double xb, double va, double vb)
{
  return (13.0*(va+vb*(-1.0+6.0*xa)-6.0*va*xb))/(360.0*(xa-xb));
}

double j_t_4_w5_full(double xa, double xb, double va, double vb)
{
  return (va+vb*(-1.0+3.0*xa)-3.0*va*xb)/(360.0*(xa-xb));
}

void deposit_j_t_segment_2(std::vector<double> &j_t,
			   int ixl,
			   int ixr,
			   double xl, 
			   double xr,
			   double vl,
			   double vr,
			   double length,			     
			   double charge,
			   int n_g)
{
  double xa, xb, va, vb, delta, q_norm, delta_left, delta_right;
  int ngp_left, ngp_right;

  q_norm = charge / length;

  find_ngp(ixl, xl, ngp_left, delta_left);
  find_ngp(ixr, xr, ngp_right, delta_right);  

  // If tracers are in one cell
  if (ngp_left == ngp_right) {
    if (length==0.0) {
      delta = delta_left;
      va = (vl + vr) / 2.0;
      j_t[mod((ngp_left-1),n_g)] += charge * va * 0.5 * pow((0.5 - delta), 2);
      j_t[mod(ngp_left,n_g)] += charge * va * 0.75 - delta * delta;
      j_t[mod((ngp_left+1),n_g)] += charge * va * 0.5 * pow((0.5 + delta), 2);
    } else { 
      xa = delta_left;
      xb = delta_right;
      va = vl;
      vb = vr;

      j_t[mod((ngp_left-1),n_g)] += q_norm * j_t_2_w1(xa, xb, va, vb);
      j_t[mod(ngp_left,n_g)] += q_norm * j_t_2_w2(xa, xb, va, vb);
      j_t[mod((ngp_left+1),n_g)] += q_norm * j_t_2_w3(xa, xb, va, vb);
    }
  }
  else {
    // Left end
    xa = delta_left;
    xb = 0.5;
    va = vl;
    vb = interpolate_segment_velocity(vl, vr, length, (xb-xa));

    j_t[mod((ngp_left-1),n_g)] += q_norm * j_t_2_w1(xa, xb, va, vb);
    j_t[mod(ngp_left,n_g)] += q_norm * j_t_2_w2(xa, xb, va, vb);
    j_t[mod((ngp_left+1),n_g)] += q_norm * j_t_2_w3(xa, xb, va, vb);
    
    // Portions fully covering a cell
    for (int cell = (ngp_left+1); cell < ngp_right; cell++) {
      xa = -0.5;
      xb = 0.5;
      va = interpolate_segment_velocity(vl, vr, length, ((cell-ngp_left)-(0.5+delta_left)));
      vb = interpolate_segment_velocity(vl, vr, length, ((cell-ngp_left)+(0.5-delta_left)));

      j_t[mod((cell-1),n_g)] += q_norm * j_t_2_w1(xa, xb, va, vb);
      j_t[mod(cell,n_g)] += q_norm * j_t_2_w2(xa, xb, va, vb);
      j_t[mod((cell+1),n_g)] += q_norm * j_t_2_w3(xa, xb, va, vb);
    }
    
    // Right end
    xa = -0.5;
    xb = delta_right;
    va = interpolate_segment_velocity(vl, vr, length, (ngp_right-ngp_left)-(0.5+delta_left));
    vb = vr;

    j_t[mod((ngp_right-1),n_g)] += q_norm * j_t_2_w1(xa, xb, va, vb);
    j_t[mod(ngp_right,n_g)] += q_norm * j_t_2_w2(xa, xb, va, vb);
    j_t[mod((ngp_right+1),n_g)] += q_norm * j_t_2_w3(xa, xb, va, vb);

  }
  return;
}

void deposit_j_t_segment_3(std::vector<double> &j_t,
			   int ixl,
			   int ixr,
			   double xl, 
			   double xr,
			   double vl,
			   double vr,
			   double length,			     
			   double charge,
			   int n_g)
{
  double xa, xb, va, vb, delta, q_norm;
  int bound_left, bound_right;

  q_norm = charge / length;

  bound_left = ixl;
  bound_right = ixr + 1;

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    if (length==0.0) {
      delta = xl - 0.5; 
      va = (vl + vr) / 2.0;     
      j_t[mod((bound_left-1),n_g)] += charge * va * (-1.0 * pow((-0.5 + delta), 3) / 6.0);
      j_t[mod(bound_left,n_g)] += charge * va * (4.0 - 6.0 * pow((0.5 + delta), 2) + 3.0 * pow((0.5 + delta), 3)) / 6.0;
      j_t[mod((bound_left+1),n_g)] += charge * va * (23.0 + 30.0*delta - 12.0*pow(delta, 2) - 24.0*pow(delta,3)) / 48.0;
      j_t[mod((bound_left+2),n_g)] += charge * va * pow((0.5 + delta), 3) / 6.0;
    } else { 
      xa = xl - 0.5;
      xb = xr - 0.5;
      va = vl;
      vb = vr;      

      j_t[mod((bound_left-1),n_g)] += q_norm * j_t_3_w1(xa, xb, va, vb);
      j_t[mod(bound_left,n_g)] += q_norm * j_t_3_w2(xa, xb, va, vb);
      j_t[mod((bound_left+1),n_g)] += q_norm * j_t_3_w3(xa, xb, va, vb);
      j_t[mod((bound_left+2),n_g)] += q_norm * j_t_3_w4(xa, xb, va, vb);
    }
  }
  else {
    // Left end
    xa = xl - 0.5;
    xb = 0.5;
    va = vl;
    vb = interpolate_segment_velocity(vl, vr, length, (xb-xa));

    j_t[mod((bound_left-1),n_g)] += q_norm * j_t_3_w1(xa, xb, va, vb);
    j_t[mod(bound_left,n_g)] += q_norm * j_t_3_w2(xa, xb, va, vb);
    j_t[mod((bound_left+1),n_g)] += q_norm *j_t_3_w3(xa, xb, va, vb);
    j_t[mod((bound_left+2),n_g)] += q_norm * j_t_3_w4(xa, xb, va, vb);
    
    // Portions connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      xa = -0.5;
      xb = 0.5;
      va = interpolate_segment_velocity(vl, vr, length, (cell-bound_left)-xl);
      vb = interpolate_segment_velocity(vl, vr, length, (cell-bound_left)+(1.0-xl));

      j_t[mod((cell-1),n_g)] += q_norm * j_t_3_w1(xa, xb, va, vb);
      j_t[mod(cell,n_g)] += q_norm * j_t_3_w2(xa, xb, va, vb);
      j_t[mod((cell+1),n_g)] += q_norm * j_t_3_w3(xa, xb, va, vb);
      j_t[mod((cell+2),n_g)] += q_norm * j_t_3_w4(xa, xb, va, vb);
    }
    
    // Right end
    xa = -0.5;
    xb = xr - 0.5;
    va = interpolate_segment_velocity(vl, vr, length, (bound_right-bound_left-1)-xl);
    vb = vr;
    
    j_t[mod((bound_right-2),n_g)] += q_norm * j_t_3_w1(xa, xb, va, vb);
    j_t[mod((bound_right-1),n_g)] += q_norm * j_t_3_w2(xa, xb, va, vb);
    j_t[mod(bound_right,n_g)] += q_norm * j_t_3_w3(xa, xb, va, vb);
    j_t[mod((bound_right+1),n_g)] += q_norm * j_t_3_w4(xa, xb, va, vb);

  }
  return;
}

void deposit_j_t_segment_4(std::vector<double> &j_t,
			   int ixl,
			   int ixr,
			   double xl, 
			   double xr,
			   double vl,
			   double vr,
			   double length,			     
			   double charge,
			   int n_g)
{
  double xa, xb, va, vb, delta, q_norm, delta_left, delta_right;
  int ngp_left, ngp_right;

  q_norm = charge / length;

  find_ngp(ixl, xl, ngp_left, delta_left);
  find_ngp(ixr, xr, ngp_right, delta_right);  

  // If tracers are in one cell
  if (ngp_left == ngp_right) {
    if (length==0.0) {
      delta = delta_left;
      va = (vl + vr) / 2.0;
      j_t[mod((ngp_left-2),n_g)] += charge * va * pow((1.0 - 2.0*delta), 4) / 384.0;
      j_t[mod((ngp_left-1),n_g)] += charge * va * (19.0 - 44.0*delta + 24.0*pow(delta, 2) + 16.0*pow(delta,3) - 16.0*pow(delta, 4))/96.0;
      j_t[mod(ngp_left,n_g)] += charge * va * (115.0 / 192.0) - (5.0*pow(delta, 2))/8.0 + pow(delta, 4)/4.0;
      j_t[mod((ngp_left+1),n_g)] += charge * va * (19.0 + 44.0*delta + 24.0*pow(delta, 2) - 16.0*pow(delta, 3) - 16*pow(delta,4))/96.0;
      j_t[mod((ngp_left+2),n_g)] += charge * va * pow((1.0 + 2.0*delta), 4) / 384.0;
	
    } else { 
      xa = delta_left;
      xb = delta_right;
      va = vl;
      vb = vr;      

      j_t[mod((ngp_left-2),n_g)] += q_norm * j_t_4_w1(xa, xb, va, vb);
      j_t[mod((ngp_left-1),n_g)] += q_norm * j_t_4_w2(xa, xb, va, vb);
      j_t[mod(ngp_left,n_g)] += q_norm * j_t_4_w3(xa, xb, va, vb);
      j_t[mod((ngp_left+1),n_g)] += q_norm * j_t_4_w4(xa, xb, va, vb);
      j_t[mod((ngp_left+2),n_g)] += q_norm * j_t_4_w5(xa, xb, va, vb);
    }
  }
  else {
    // Left end
    xa = delta_left;
    xb = 0.5;
    va = vl;
    vb = interpolate_segment_velocity(vl, vr, length, (xb-xa));    

    j_t[mod((ngp_left-2),n_g)] += q_norm * j_t_4_w1(xa, xb, va, vb);
    j_t[mod((ngp_left-1),n_g)] += q_norm * j_t_4_w2(xa, xb, va, vb);
    j_t[mod(ngp_left,n_g)] += q_norm * j_t_4_w3(xa, xb, va, vb);
    j_t[mod((ngp_left+1),n_g)] += q_norm * j_t_4_w4(xa, xb, va, vb);
    j_t[mod((ngp_left+2),n_g)] += q_norm * j_t_4_w5(xa, xb, va, vb);
    
    // Portions fully covering a cell
    for (int cell = (ngp_left+1); cell < ngp_right; cell++) {
      xa = -0.5;
      xb = 0.5;
      va = interpolate_segment_velocity(vl, vr, length, (cell-ngp_left)-(0.5+xl));
      vb = interpolate_segment_velocity(vl, vr, length, (cell-ngp_left)+(0.5-xl));

      j_t[mod((cell-2),n_g)] += q_norm * j_t_4_w1(xa, xb, va, vb);
      j_t[mod((cell-1),n_g)] += q_norm * j_t_4_w2(xa, xb, va, vb);
      j_t[mod(cell,n_g)] += q_norm * j_t_4_w3(xa, xb, va, vb);
      j_t[mod((cell+1),n_g)] += q_norm * j_t_4_w4(xa, xb, va, vb);
      j_t[mod((cell+2),n_g)] += q_norm * j_t_4_w5(xa, xb, va, vb);
    }
    
    // Right end
    xa = -0.5;
    xb = delta_right;
    va = interpolate_segment_velocity(vl, vr, length, (ngp_right-ngp_left)-(0.5+xl));
    vb = vr;    

    j_t[mod((ngp_right-2),n_g)] += q_norm * j_t_4_w1(xa, xb, va, vb);
    j_t[mod((ngp_right-1),n_g)] += q_norm * j_t_4_w2(xa, xb, va, vb);
    j_t[mod(ngp_right,n_g)] += q_norm * j_t_4_w3(xa, xb, va, vb);
    j_t[mod((ngp_right+1),n_g)] += q_norm * j_t_4_w4(xa, xb, va, vb);
    j_t[mod((ngp_right+2),n_g)] += q_norm * j_t_4_w5(xa, xb, va, vb);
  }
  return;
}

void ParticleSpecies::deposit_j_y_sic(std::vector<double> &j_y)
{
  int ix_a, ix_b, ix_l, ix_r;
  double x_a, x_b, x_l, x_r, v_a, v_b, v_l, v_r, length;
  for (int i = 0; i < n_p; i++) {
    average(ix_old[i], x_old[i], ix[i], x[i], ix_a, x_a);
    average(ix_old[i+1], x_old[i+1], ix[i+1], x[i+1], ix_b, x_b);
    v_a = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));    
    v_b = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    
    arrange_segment_velocity(ix_a, x_a,
			     ix_b, x_b,
			     v_a, v_b,
			     ix_l, x_l,
			     ix_r, x_r,
			     length,
			     v_l, v_r);
    
    switch (method) {
    case 5 :
      deposit_j_t_segment_0(j_y, ix_l, ix_r, x_l, x_r, v_l, v_r, length, charge[i], n_g);
      break;
    case 6 :
      deposit_j_t_segment_1(j_y, ix_l, ix_r, x_l, x_r, v_l, v_r, length, charge[i], n_g);
      break;      
    case 7 :
      deposit_j_t_segment_2(j_y, ix_l, ix_r, x_l, x_r, v_l, v_r, length, charge[i], n_g);
      break;      
    case 8 :
      deposit_j_t_segment_3(j_y, ix_l, ix_r, x_l, x_r, v_l, v_r, length, charge[i], n_g);
      break;      
    case 9 :
      deposit_j_t_segment_4(j_y, ix_l, ix_r, x_l, x_r, v_l, v_r, length, charge[i], n_g);
      break;      
    default:
      std::cout << "Error, selected interpolation order not implemented." << std::endl;
    }

  }
  return;
}

void ParticleSpecies::deposit_j_z_sic(std::vector<double> &j_z)
{
  int ix_a, ix_b, ix_l, ix_r;
  double x_a, x_b, x_l, x_r, v_a, v_b, v_l, v_r, length;
  for (int i = 0; i < n_p; i++) {
    average(ix_old[i], x_old[i], ix[i], x[i], ix_a, x_a);
    average(ix_old[i+1], x_old[i+1], ix[i+1], x[i+1], ix_b, x_b);
    v_a = u_z[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    v_b = u_z[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    arrange_segment_velocity(ix_a, x_a,
			     ix_b, x_b,
			     v_a, v_b,
			     ix_l, x_l,
			     ix_r, x_r,
			     length,
			     v_l, v_r);
    switch (method) {
    case 5 :
      deposit_j_t_segment_0(j_z, ix_l, ix_r, x_l, x_r, v_l, v_r, length, charge[i], n_g);
      break;
    case 6 :
      deposit_j_t_segment_1(j_z, ix_l, ix_r, x_l, x_r, v_l, v_r, length, charge[i], n_g);
      break;
    case 7 :
      deposit_j_t_segment_2(j_z, ix_l, ix_r, x_l, x_r, v_l, v_r, length, charge[i], n_g);
      break;
    case 8 :
      deposit_j_t_segment_3(j_z, ix_l, ix_r, x_l, x_r, v_l, v_r, length, charge[i], n_g);
      break;
    case 9 :
      deposit_j_t_segment_4(j_z, ix_l, ix_r, x_l, x_r, v_l, v_r, length, charge[i], n_g);
      break;
    default:
      std::cout << "Error, selected interpolation order not implemented." << std::endl;
    }
  }
  return;
}

void ParticleSpecies::write_phase(int species_number, int t, int my_rank)
{
  std::string ix_filename;  
  std::string x_filename;
  std::string u_x_filename;
  
  std::stringstream ss;
  ss << "ix_";
  ss << species_number;
  ss << "_";
  ss << t;
  ss << "_";
  ss << my_rank;
  ix_filename = ss.str();

  ss.str(std::string());
  ss.clear();
  ss << "x_";
  ss << species_number;
  ss << "_";
  ss << t;
  ss << "_";
  ss << my_rank;
  x_filename = ss.str();  

  ss.str(std::string());
  ss.clear();
  ss << "u_x_";
  ss << species_number;
  ss << "_";
  ss << t;
  ss << "_";
  ss << my_rank;
  u_x_filename = ss.str();
  
  data_to_file_int(ix, ix_filename);
  data_to_file(x, x_filename);
  data_to_file(u_x, u_x_filename);
  return;
}

void ParticleSpecies::write_particle_diagnostics(int n_t, int my_rank, MPI_Comm COMM)
{
  sum_array_to_root(&energy_history[0], n_t, COMM, my_rank);
  sum_array_to_root(&momentum_x_history[0], n_t, COMM, my_rank);
  sum_array_to_root(&momentum_y_history[0], n_t, COMM, my_rank);
  sum_array_to_root(&momentum_z_history[0], n_t, COMM, my_rank);  
  sum_array_to_root(&n_p_history[0], n_t, COMM, my_rank);

  if (my_rank==0) {
    data_to_file(energy_history, (species_name+"_ene"));
    data_to_file(momentum_x_history, (species_name+"_p_x"));
    data_to_file(momentum_y_history, (species_name+"_p_y"));
    data_to_file(momentum_z_history, (species_name+"_p_z"));
    data_to_file(n_p_history, (species_name+"_n_p"));
  }
  return; 
}

void ParticleSpecies::advance_x()
{
  double gamma;
  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    x[i] = x[i] + (dt/dx) * u_x[i] / gamma;
    if (x[i] >= 1.0) {
      ix[i] = ix[i] + 1;
      x[i] = x[i] - 1.0;
    } else if (x[i] < 0.0) {
      ix[i] = ix[i] - 1;
      x[i] =  x[i] + 1.0;
    }
  } 
  return;
}

void ParticleSpecies::communicate_ghost_particles(MPI_Comm COMM)
{
  MPI_Status status;
  int num_procs, my_rank, dest, source;
  int tag = 0;
  
  MPI_Comm_size(COMM, &num_procs);
  MPI_Comm_rank(COMM, &my_rank);
  dest = mod(my_rank-1, num_procs);
  source = mod(my_rank+1, num_procs);

  MPI_Sendrecv(&ix[0], 2, MPI_INT, dest, tag,
	       &ix[n_p], 2, MPI_INT,
	       source, tag, COMM, &status);
  MPI_Sendrecv(&ix_old[0], 2, MPI_INT, dest, tag,
	       &ix_old[n_p], 2, MPI_INT,
	       source, tag, COMM, &status);  
  MPI_Sendrecv(&x[0], 2, MPI_DOUBLE, dest, tag,
	       &x[n_p], 2, MPI_DOUBLE,
	       source, tag, COMM, &status);
  MPI_Sendrecv(&x_old[0], 2, MPI_DOUBLE, dest, tag,
	       &x_old[n_p], 2, MPI_DOUBLE,
	       source, tag, COMM, &status);
  MPI_Sendrecv(&charge[0], 2, MPI_DOUBLE, dest, tag,
	       &charge[n_p], 2, MPI_DOUBLE,
	       source, tag, COMM, &status);
  MPI_Sendrecv(&u_x[0], 2, MPI_DOUBLE, dest, tag,
	       &u_x[n_p], 2, MPI_DOUBLE,
	       source, tag, COMM, &status);
  MPI_Sendrecv(&u_y[0], 2, MPI_DOUBLE, dest, tag,
	       &u_y[n_p], 2, MPI_DOUBLE,
	       source, tag, COMM, &status);
  MPI_Sendrecv(&u_z[0], 2, MPI_DOUBLE, dest, tag,
	       &u_z[n_p], 2, MPI_DOUBLE,
	       source, tag, COMM, &status);
  MPI_Sendrecv(&lagrangian_id[0], 2, MPI_DOUBLE, dest, tag,
	       &lagrangian_id[n_p], 2, MPI_DOUBLE,
	       source, tag, COMM, &status);
  
  lagrangian_id[n_p] += n_ppp;  
  lagrangian_id[n_p+1] += n_ppp;
  
  if (my_rank==(num_procs-1)) {
    ix[n_p] += n_g;
    ix[n_p+1] += n_g;
    ix_old[n_p] += n_g;
    ix_old[n_p+1] += n_g;
  }
  
  return;
}

double lagrange_3(double *x, double *y, double x_sample)
{
  double y_interpolated = 0.0;
  double p_i;
  for (int i = 0; i < 3; i++) {
    p_i = 1.0;
    for (int j = 0; j < 3; j++) {
      if (j==i) {
	continue;
      } else {
	p_i *= (x_sample - x[j]) / (x[i] - x[j]);
      }
    }
    y_interpolated += y[i] * p_i;
  }
  return y_interpolated;
}

double lagrange_3_position(double *x, int *iy, double *y, double x_sample)
{
  double y_interpolated = 0.0;
  double p_i;
  for (int i = 0; i < 3; i++) {
    p_i = 1.0;
    for (int j = 0; j < 3; j++) {
      if (j==i) {
	continue;
      } else {
	p_i *= (x_sample - x[j]) / (x[i] - x[j]);
      }
    }
    if (i==0) {
      y_interpolated += y[i] * p_i;
    } else {
      y_interpolated += (y[i]+(iy[i]-iy[0])) * p_i;
    }
  }
  return y_interpolated;
}

void ParticleSpecies::split_segment_lagrange_3(int i)
{
  double x_avg, x_interpolated;
  int ix_avg, ix_interpolated;

  double new_id = (lagrangian_id[i+1]+lagrangian_id[i])/2.0;  

  x_interpolated = lagrange_3_position(&lagrangian_id[i], &ix[i], &x[i], new_id);
  ix_interpolated = floor(x_interpolated);

  ix.insert(ix.begin()+i+1, ix_interpolated + ix[i]);
  x.insert(x.begin()+i+1, (x_interpolated-ix_interpolated));

  // Linear interpolation on x_old to be consistent with Gauss's law
  average(ix_old[i], x_old[i], ix_old[i+1], x_old[i+1], ix_avg, x_avg);
  x_old.insert(x_old.begin()+i+1, x_avg);
  ix_old.insert(ix_old.begin()+i+1, ix_avg);

  u_x.insert(u_x.begin()+i+1, lagrange_3(&lagrangian_id[i], &u_x[i], new_id));
  u_y.insert(u_y.begin()+i+1, lagrange_3(&lagrangian_id[i], &u_y[i], new_id));
  u_z.insert(u_z.begin()+i+1, lagrange_3(&lagrangian_id[i], &u_z[i], new_id));
  lagrangian_id.insert(lagrangian_id.begin()+i+1, new_id);
  charge[i] *= 0.5;
  charge.insert(charge.begin()+i+1, charge[i]);

  n_p += 1;
  return;
}

void ParticleSpecies::split_segment_linear(int i)
{
  int ix_avg;
  double x_avg;  
  double new_id = (lagrangian_id[i+1]+lagrangian_id[i])/2.0;

  average(ix[i], x[i], ix[i+1], x[i+1], ix_avg, x_avg);
  x.insert(x.begin()+i+1, x_avg);
  ix.insert(ix.begin()+i+1, ix_avg);
  
  average(ix_old[i], x_old[i], ix_old[i+1], x_old[i+1], ix_avg, x_avg);
  x_old.insert(x_old.begin()+i+1, x_avg);
  ix_old.insert(ix_old.begin()+i+1, ix_avg);
  
  u_x.insert(u_x.begin()+i+1, (u_x[i+1] + u_x[i])/2.0);
  u_y.insert(u_y.begin()+i+1, (u_y[i+1] + u_y[i])/2.0);
  u_z.insert(u_z.begin()+i+1, (u_z[i+1] + u_z[i])/2.0);
  lagrangian_id.insert(lagrangian_id.begin()+i+1, new_id);
  charge[i] *= 0.5;
  charge.insert(charge.begin()+i+1, charge[i]);
  
  n_p += 1;
  return;
}

void ParticleSpecies::refine_segments(double refinement_length)
{
  double length;
  int i = 0;
  double max = 0.0;
  
  while (i < n_p) {
    length = fabs(a_minus_b(ix[i+1], x[i+1], ix[i], x[i]));
    
    if (length > max) {
      max = length;
    }
    if (length > refinement_length) {
      split_segment_lagrange_3(i);
      //split_segment_linear(i);
    } else {
      i+=1;
    }
  }
  std::cout << max << std::endl;
  return;
}

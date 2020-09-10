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

int get_nearest_gridpoint(double x)
{
  int nearest_gridpoint;
  if ((x-floor(x)) <= 0.5) {
    nearest_gridpoint = floor(x);
  } else {
    nearest_gridpoint = ceil(x);
  }
  return nearest_gridpoint;
}

double interpolate_field_integer(std::vector<double> &field, double x, double dx, int n_g)
{
  int i_lower;
  double alpha, beta;
  // Normalize x to cell length
  x = x / dx;
  // index of left bounding gridpoint
  i_lower = floor(x);
  alpha = 1.0 - (x - i_lower);
  beta = 1.0 - alpha;
  return (alpha * field[mod(i_lower,n_g)]) + (beta * field[mod((i_lower+1),n_g)]);
}

double interpolate_field_half_integer(std::vector<double> &field, double x, double dx, 
				      int n_g)
{
  int i_lower;
  double alpha, beta;
  // Shift position so that zero lies at the first half-integer gridpoint
  x = x - dx / 2.0;
  // Normalize x to cell length
  x = x / dx;
  // Index of left bounding gridpoint
  i_lower = floor(x);
  alpha = 1.0 - (x - i_lower);
  beta = 1.0 - alpha;
  return (alpha * field[mod(i_lower,n_g)]) + (beta * field[mod((i_lower+1),n_g)]);
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

  long long i_start, i_end;
  i_start = n_p * my_rank;
  i_end = i_start + n_p;

  x.resize(n_p);
  u_x.resize(n_p);
  u_y.resize(n_p);
  u_z.resize(n_p);  
  x_old.resize(n_p);
  charge.resize(n_p);
  lagrangian_id.resize(n_p);
  
  std::stringstream ss;
  ss << "particles_";
  ss << species_number;
  species_name = ss.str();

  double particle_spacing = dx / n_ppc;
  
  // Add ghost tracer particles if using line segments
  if (method>1) {
    for (int i = 0; i < 2; i++) {
      charge.push_back(0.0);
      u_x.push_back(0.0);
      u_y.push_back(0.0);
      u_z.push_back(0.0);      
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
    k = (2.0 * (4.0*atan(1.0)) / (n_g * dx));
    //    u_x_1 = 0.0025;
    u_x_1 = 0.1;
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
    
    k = 10.0 * (2.0 * (4.0*atan(1.0)) / (n_g * dx));
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
    
    k = (2.0 * (4.0*atan(1.0)) / (n_g * dx));
    u_x_1 = 0.00176425;
 }


  this->rqm = rqm_vector[species_number];
  this->method = method;

  for (long long i = i_start; i < i_end; i++) {
    lagrangian_id[i-i_start] = i - i_start;    
    charge[i-i_start] = density[species_number] * (1.0 / n_ppc);
    x[i-i_start] = ((long double) i) * particle_spacing + particle_spacing / 2.0;
    x_old[i-i_start] = ((long double) i) * particle_spacing + particle_spacing / 2.0;
    
    u_x[i-i_start] = u_x_drift[species_number];    
    u_y[i-i_start] = u_y_drift[species_number];
    u_z[i-i_start] = 0.0;
    
    if (simulation_type == -2 or simulation_type == 0 or simulation_type==-1) {
      u_x[i-i_start] += u_x_1 * sin(k * x[i-i_start]);
    } else if (simulation_type == 2 and species_number == 0) {
      u_z[i-i_start] += (-2.0 * 0.273055) * cos(0.44526860656 * x[i-i_start]);
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
  int ngp_integer, ngp_half_integer;
  
  for (int i = 0; i < n_p; i++) {
    b_x_particle = b_x[0]; // Uniform field
    if (interp_order==0) {
      if (center_fields) {
	ngp_integer = get_nearest_gridpoint(x[i] / dx);
	e_x_particle = e_x_int[mod(ngp_integer, n_g)];
	e_y_particle = e_y[mod(ngp_integer, n_g)];
	e_z_particle = e_z[mod(ngp_integer, n_g)];
	b_y_particle = b_y_int[mod(ngp_integer, n_g)];	
	b_z_particle = b_z_int[mod(ngp_integer, n_g)];
      } else {
	ngp_integer = get_nearest_gridpoint(x[i] / dx);
	ngp_half_integer = get_nearest_gridpoint((x[i]-0.5) / dx);
	e_x_particle = e_x[mod(ngp_half_integer, n_g)];
	e_y_particle = e_y[mod(ngp_integer, n_g)];
	e_z_particle = e_z[mod(ngp_integer, n_g)];	
	b_y_particle = b_y[mod(ngp_half_integer, n_g)];
	b_z_particle = b_z[mod(ngp_half_integer, n_g)];	
      }
    } else {
      if (center_fields) {
	e_x_particle = interpolate_field_integer(e_x_int, x[i], dx, n_g);
	e_y_particle = interpolate_field_integer(e_y, x[i], dx, n_g);
	e_z_particle = interpolate_field_integer(e_z, x[i], dx, n_g);
	b_y_particle = interpolate_field_integer(b_y_int, x[i], dx, n_g);
	b_z_particle = interpolate_field_integer(b_z_int, x[i], dx, n_g);
      } else {
	e_x_particle = interpolate_field_half_integer(e_x, x[i], dx, n_g);
	e_y_particle = interpolate_field_integer(e_y, x[i], dx, n_g);
	e_z_particle = interpolate_field_integer(e_z, x[i], dx, n_g);
	b_y_particle = interpolate_field_half_integer(b_y, x[i], dx, n_g);
	b_z_particle = interpolate_field_half_integer(b_z, x[i], dx, n_g);
      }
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
  //  double energy = 0.0;
  //  double momentum_x = 0.0;
  //  double momentum_y = 0.0;
  //  double momentum_z = 0.0;  
  int ngp_integer, ngp_half_integer;
  
  for (int i = 0; i < n_p; i++) {
    b_x_particle = b_x[0]; // Uniform field
    if (interp_order==0) {
      if (center_fields) {
	ngp_integer = get_nearest_gridpoint(x[i] / dx);
	e_x_particle = e_x_int[mod(ngp_integer, n_g)];
	e_y_particle = e_y[mod(ngp_integer, n_g)];
	e_z_particle = e_z[mod(ngp_integer, n_g)];
	b_y_particle = b_y_int[mod(ngp_integer, n_g)];	
	b_z_particle = b_z_int[mod(ngp_integer, n_g)];
      } else {
	ngp_integer = get_nearest_gridpoint(x[i] / dx);
	ngp_half_integer = get_nearest_gridpoint((x[i]-0.5) / dx);
	e_x_particle = e_x[mod(ngp_half_integer, n_g)];
	e_y_particle = e_y[mod(ngp_integer, n_g)];
	e_z_particle = e_z[mod(ngp_integer, n_g)];	
	b_y_particle = b_y[mod(ngp_half_integer, n_g)];
	b_z_particle = b_z[mod(ngp_half_integer, n_g)];	
      }
    } else {
      if (center_fields) {
	e_x_particle = interpolate_field_integer(e_x_int, x[i], dx, n_g);
	e_y_particle = interpolate_field_integer(e_y, x[i], dx, n_g);
	e_z_particle = interpolate_field_integer(e_z, x[i], dx, n_g);
	b_y_particle = interpolate_field_integer(b_y_int, x[i], dx, n_g);
	b_z_particle = interpolate_field_integer(b_z_int, x[i], dx, n_g);
      } else {
	e_x_particle = interpolate_field_half_integer(e_x, x[i], dx, n_g);
	e_y_particle = interpolate_field_integer(e_y, x[i], dx, n_g);
	e_z_particle = interpolate_field_integer(e_z, x[i], dx, n_g);
	b_y_particle = interpolate_field_half_integer(b_y, x[i], dx, n_g);
	b_z_particle = interpolate_field_half_integer(b_z, x[i], dx, n_g);
      }
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

void deposit_rho_sic_0_segment(std::vector<double> &rho,
			       double x_tracer_a, 
			       double x_tracer_b,
			       double charge,
			       int n_g,
			       double dx)
{
  double left, right, length, charge_fraction;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  length = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    if ((left == bound_left + 0.5) && (right == bound_left + 0.5)) {
      rho[mod((bound_right),n_g)] += charge;
    } else if (left >= (bound_left + 0.5)) {
      rho[mod((bound_right),n_g)] += charge;
    } else if (right <= (bound_left + 0.5)) {
      rho[mod((bound_left),n_g)] += charge;
    } else {
      rho[mod((bound_left),n_g)] += (bound_left + 0.5 - left) / length * charge;
      rho[mod((bound_right),n_g)] += (right - (bound_left + 0.5)) / length * charge;
    }
  }
  else {
    // Left end
    if (left < bound_left + 0.5) {
      rho[mod(bound_left,n_g)] += charge * (bound_left + 0.5 - left) / length;
      rho[mod((bound_left+1),n_g)] += charge * (0.5) / length;
    } else { 
      rho[mod(bound_left+1,n_g)] += charge * (bound_left + 1.0 - left) / length;
    }

    // Pieces connecting two gridpoints
    charge_fraction = (1.0) / length;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      rho[mod(cell,n_g)] += charge * charge_fraction * 0.5;
      rho[mod((cell+1),n_g)] += charge * charge_fraction * 0.5;
    }
    
    // Right end
    if (right > bound_right - 0.5) {
      rho[mod(bound_right-1,n_g)] += charge * (0.5) / length;
      rho[mod(bound_right,n_g)] += charge * (right - (bound_right-0.5)) / length;
    } else { 
      rho[mod(bound_right-1,n_g)] += charge * (right - (bound_right-1.0)) / length;
    }
  }
  return;
}

void deposit_rho_sic_1_segment(std::vector<double> &rho,
			       double x_tracer_a, 
			       double x_tracer_b,
			       double charge,
			       int n_g,
			       double dx)
{
  double left, right, length, midpoint, charge_fraction, weight;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  length = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    midpoint = (left + right) / 2.0;
    charge_fraction = 1.0;
    weight = bound_right - midpoint;
    rho[mod(bound_left,n_g)] += charge * charge_fraction * weight;
    rho[mod((bound_left+1),n_g)] += charge * charge_fraction * (1.0 - weight);
  }
  else {
    // Left end
    midpoint = (left + bound_left + 1) / 2.0;
    charge_fraction = (bound_left + 1 - left) / length;
    weight = (bound_left + 1 - midpoint);
    rho[mod(bound_left,n_g)] += charge * charge_fraction * weight;
    rho[mod((bound_left+1),n_g)] += charge * charge_fraction * (1.0 - weight);

    // Pieces connecting two gridpoints
    charge_fraction = (1.0) / length;

    weight = 0.5;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      rho[mod(cell,n_g)] += charge * charge_fraction * weight;
      rho[mod((cell+1),n_g)] += charge * charge_fraction * (1.0 - weight);
    }
    
    // Right end
    midpoint = (right + bound_right - 1) / 2.0;
    charge_fraction = (right - (bound_right - 1)) / length;
    weight = bound_right - midpoint;
    rho[mod((bound_right-1),n_g)] += charge * charge_fraction * weight;
    rho[mod(bound_right,n_g)] += charge * charge_fraction * (1.0-weight);
  }
  return;
}

void deposit_charge_to_left_segment_linear(std::vector<double> &j_x,
					   double x_tracer_a, 
					   double x_tracer_b,
					   double charge,
					   int n_g,
					   double dx,
					   int right_max)
{
  double left, right, length, midpoint, charge_fraction, weight;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  length = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    midpoint = (left + right) / 2.0;
    charge_fraction = 1.0;
    weight = bound_right - midpoint;
    j_x[mod(bound_left,n_g)] += charge * charge_fraction * weight;
    for (int cell = (bound_left+1); cell < right_max; cell++) {
      j_x[mod(cell,n_g)] += charge * charge_fraction;
    }
  }
  else {
    // Left end
    midpoint = (left + bound_left + 1) / 2.0;
    charge_fraction = (bound_left + 1 - left) / length;
    weight = (bound_left + 1 - midpoint);
    j_x[mod(bound_left,n_g)] += charge * charge_fraction * weight;
    for (int cell = (bound_left+1); cell < right_max; cell++) {
      j_x[mod(cell,n_g)] += charge * charge_fraction;
    }

    // Pieces connecting two gridpoints
    charge_fraction = (1.0) / length;

    weight = 0.5;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      j_x[mod(cell,n_g)] += charge * charge_fraction * weight;
      for (int boundary = (cell+1); boundary < right_max; boundary++) {
	j_x[mod(boundary,n_g)] += charge * charge_fraction;
      }
    }
    
    // Right end
    midpoint = (right + bound_right - 1) / 2.0;
    charge_fraction = (right - (bound_right - 1)) / length;
    weight = bound_right - midpoint;
    j_x[mod((bound_right-1),n_g)] += charge * charge_fraction * weight;
    for (int cell = (bound_right); cell < right_max; cell++) {
      j_x[mod(cell,n_g)] += charge * charge_fraction;
    }
  }
  return;
}


void ParticleSpecies::deposit_j_x_sic_1(std::vector<double> &j_x)
{
  double right_max;
  for (int i = 0; i < n_p; i++) {
    // Index of right bounding gridpoint
    right_max = fmax(fmax(x_old[i], x_old[i+1]),fmax(x[i],x[i+1]));
    right_max = right_max / dx;
    right_max = ceil(right_max);
    
    deposit_charge_to_left_segment_linear(j_x, x_old[i], x_old[i+1], (charge[i] * (dx / dt)), n_g, dx, right_max);
    deposit_charge_to_left_segment_linear(j_x, x[i], x[i+1], (-1.0) * (charge[i] * (dx / dt)), n_g, dx, right_max);
  }
  return;
}

// Higher order routines below

double full_xg(double nl, double nr, double n0, double l0, int xg, double  x0, double xa, double xb)
{
  return  l0*n0*((1.0-l0/3.0-x0+xg)-l0*nr/(3.0*(nl+nr)));
}

double full_xgp1(double nl, double nr, double n0, double l0, int xg, double x0, double xa, double xb)
{
  return  l0*n0*(l0+l0*nr/(nl+nr)+3.0*(x0-xg))/3.0;
}

double general_xg(double nl, double nr, double n0, double l0, int xg, double  x0, double xa, double xb)
{
  return n0*(xa-xb)/(3.0*l0*(nl+nr))*(3*l0*nl*(-2.0+xa+xb-2.0*xg)+(nl-nr)*(3.0*(xa+xb)-2.0*(xa*xa+xa*xb+xb*xb)+3.0*x0*(-2.0+xa+xb-2.0*xg)+3.0*(xa+xb)*xg));
}

double general_xgp1(double nl, double nr, double n0, double l0, int xg, double x0, double xa, double xb)
{
  return n0*(xa-xb)/(3.0*l0*(nl+nr))*(-3.0*l0*nl*(xa+xb-2.0*xg)-(nl-nr)*(-2.0*(xa*xa+xa*xb+xb*xb)+3.0*x0*(xa+xb-2.0*xg)+3.0*(xa+xb)*xg));
}

double cover_xg(double nl, double nr, double n0, double l0, int xg, double x0, double xa, double xb)
{      
  return n0*(3.0*nl*l0+(nl-nr)*(-1.0+3.0*(x0-xg)))/(3*l0*(nl+nr));
}

double cover_xgp1(double nl, double nr, double n0, double l0, int xg, double x0, double xa, double xb)
{      
  return n0*(3.0*nl*l0+(nl-nr)*(-2.0+3.0*(x0-xg)))/(3*l0*(nl+nr));
}

double left_end_xg(double nl, double nr, double n0, double l0, int xg, double x0, double xa, double xb)
{
  return n0*pow(1.0-x0+xg,2)*(nl*3.0*l0+(nr-nl)*(1.0-x0+xg))/(3.0*l0*(nl+nr));
}

double left_end_xgp1(double nl, double nr, double n0, double l0, int xg, double x0, double xa, double xb)
{
  return n0*(1.0-x0+xg)/(3.0*l0*(nl+nr))*((nl-nr)*(x0-xg-1.0)*(2.0+x0-xg)+nl*(3.0*l0*(1.0+x0-xg)));
}

double right_end_xg(double nl, double nr, double n0, double l0, int xg, double x0, double xa, double xb)
{
  return n0*(xg-x0-l0)*(l0*l0*(nl+2.0*nr)+3.0*l0*nl+l0*(2.0*nl+nr)*(-3.0+x0-xg)+(nl-nr)*(-3.0+x0-xg)*(x0-xg))/(3.0*l0*(nl+nr));
}

double right_end_xgp1(double nl, double nr, double n0, double l0, int xg, double x0, double xa, double xb)
{
  return (n0*(l0*(nl+2.0*nr)+(nl-nr)*(x0-xg))*pow(l0+x0-xg,2))/(3.0*l0*(nl+nr));
}

double q_unweighted(double nl, double nr, double n0, double l0, double x0, double xa, double xb)
{
  return (n0*(-2.0*l0*nl-(nl-nr)*(2*x0-xa-xb))*(xa-xb))/(l0*(nl+nr));
}

double full_j_t_xg(double nl, double nr, double vl, double vr, double charge, double l0, int xg, double x0)
{
  return charge*(-l0*(nl*(vl+vr)+nr*(vl+3.0*vr))-2.0*(nl*(2.0*vl+vr)+nr*(vl+2.0*vr))*(-1.0+x0-xg))/(6.0*(nl+nr));
}

double full_j_t_xgp1(double nl, double nr, double vl, double vr, double charge, double l0, int xg, double x0)
{
  return charge*(l0*(nl+nr)*vl+l0*(nl+3.0*nr)*vr+2.0*(nl*(2.0*vl+vr)+nr*(vl+2.0*vr))*(x0-xg))/(6.0*(nl+nr));
}

double cover_j_t_xg(double nl, double nr, double vl, double vr, double n0, double l0, int xg, double x0, double xa, double xb)
{      
  return n0/(6.0*l0*l0*(nl+nr))*(nr*(vr+2.0*vr*(-2.0+3.0*x0-3.0*xg)*(x0-xg)+
	   vl*(-1.0+4.0*x0-6.0*pow(x0-xg,2)-4.0*xg+l0*(2.0-6.0*x0+6.0*xg)))+ 
	   nl*(vl*(1.0+6.0*l0*l0-4.0*x0+4.0*l0*(-1.0+3.0*x0-3.0*xg)+6.0*pow(x0-xg,2)+4.0*xg)+
	   vr*(-1.0+4.0*x0-6.0*pow(x0-xg,2)-4.0*xg+l0*(2.0-6.0*x0+6.0*xg))));
}

double cover_j_t_xgp1(double nl, double nr, double vl, double vr, double n0, double l0, int xg, double x0, double xa, double xb)
{      
  return n0/(6.0*l0*l0*(nl+nr))*(3.0*nr*vr+2.0*nr*vr*(-4.0+3.0*x0-3.0*xg)*(x0-xg)+
	   nl*vl*(3.0+6.0*l0*l0-8.0*x0+4.0*l0*(-2.0+3.0*x0-3.0*xg)+6.0*pow(x0-xg,2)+8.0*xg)+ 
	   nr*vl*(-3.0+8.0*x0-6.0*pow(x0-xg,2)-8.0*xg+l0*(4.0-6.0*x0+6.0*xg))+ 
	   nl*vr*(-3.0+8.0*x0-6.0*pow(x0-xg,2)-8.0*xg+l0*(4.0-6.0*x0+6.0*xg)));
}

double left_end_j_t_xg(double nl, double nr, double vl, double vr, double n0, double l0, int xg, double x0, double xa, double xb)
{
  return n0*pow(1.0-x0+xg,2)*(-1.0*nr*(-1.0+x0-xg)*(vl*(-1.0+2.0*l0+x0-xg)+ 
           vr*(1.0-x0+xg))+nl*(-1.0*vr*(-1.0+x0-xg)*(-1.0+2.0*l0+x0-xg)+
           vl*(6.0*l0*l0+4.0*l0*(-1.0+x0-xg)+pow(1.0-x0+xg,2))))/(6.0*l0*l0*(nl+nr));
}

double left_end_j_t_xgp1(double nl, double nr, double vl, double vr, double n0, double l0, int xg, double x0, double xa, double xb)
{
  return n0/(6.0*l0*l0*(nl+nr))*(nr*(2.0*l0*vl*(2.0+x0-xg)+vl*(-1.0+x0-xg)*(3.0+x0-xg)-vr*(-1.0+x0-xg)*(3.0+x0-xg))*pow(1.0-x0+xg,2)-
           nl*(-1.0+x0-xg)*(-vr*(2.0*l0*(2.0+x0-xg)+(-1.0+x0-xg)*(3.0+x0-xg))*(-1.0+x0-xg)+vl*(6.0*l0*l0*(1.0+x0-xg)+4.0*l0*(-1.0+x0-xg)*(2.0+x0-xg)+(3.0+x0-xg)*pow(1-x0+xg,2))));
}

double right_end_j_t_xg(double nl, double nr, double vl, double vr, double n0, double l0, int xg, double x0, double xa, double xb)
{
  return n0/(6.0*l0*l0*(nl+nr))*(-pow(l0,4)*(nl*(vl+vr)+nr*(vl+3.0*vr))-2.0*pow(l0,3)*(nl*(2.0*vl+vr)+nr*(vl+2.0*vr))*(-1.0+x0-xg)-
				  6.0*l0*l0*nl*vl*(-2.0+x0-xg)*(x0-xg) -2.0*l0*(2.0*nl*vl-nr*vl-nl*vr)*(-3.0+x0-xg)*pow(x0-xg,2) -(nl-nr)*(vl-vr)*(-4.0+x0-xg)*pow(x0-xg,3));
}

double right_end_j_t_xgp1(double nl, double nr, double vl, double vr, double n0, double l0, int xg, double x0, double xa, double xb)
{
  return n0*(l0*l0*(nl*(vl+vr)+nr*(vl+3.0*vr))+2.0*l0*(nl*vl-nr*vr)*(x0-xg)+(nl-nr)*(vl-vr)*pow(x0-xg,2))*pow(l0+x0-xg,2)/(6.0*l0*l0*(nl+nr));
}

double j_t_unweighted(double nl, double nr, double vl, double vr, double n0, double l0, double x0, double xa, double xb)
{
  return n0/(3.0*l0*l0*(nl+nr))*(3.0*l0*(2.0*nl*vl-nr*vl-nl*vr)*(2.0*x0-xa-xb)*(xb-xa)+6.0*l0*l0*nl*vl*(xb-xa)+2.0*(nl-nr)*(vl-vr)*(xb-xa)*(3.0*x0*x0+xa*xa+xa*xb+xb*xb-3.0*x0*(xa+xb)));
}

double j_t_unweighted_full(double nl, double nr, double vl, double vr, double charge)
{
  return (charge*(nl*(2.0*vl+vr)+nr*(vl+2.0*vr)))/(3.0*(nl+nr));
}
		    
void deposit_charge_to_left_segment_higher_order_1(std::vector<double> &j_x,
						   double x_tracer_a, 
						   double x_tracer_b,
						   double charge,
						   int n_g,
						   double dx,
						   int right_max,
						   double nl,
						   double nr)
{
  double left, right, xg, x0, l0, n0;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  x0 = left;
  l0 = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);
    
  if (l0 == 0.0) {

    j_x[mod(bound_left,n_g)] += charge * (bound_right - left);
    for (int cell = (bound_left+1); cell < right_max; cell++) {
      j_x[mod(cell,n_g)] += charge;
    }

  } else {
    
    n0 = charge / l0;

    if (left == x_tracer_b / dx) {
      std::swap(nl, nr);
    }

    if (nl == 0.0 && nr == 0.0) {
      nl = 1.0;
      nr = 1.0;
    } else if (nl == 0.0) {
      nl = 1.0;
      nr = 0.0;
    } else if (nr == 0.0) {
      nl = 0.0;
      nr = 1.0;
    }

    if (bound_right == (bound_left + 1)) {
      xg = bound_left;
      j_x[mod(bound_left,n_g)] += full_xg(nl, nr, n0, l0, xg, x0, x0, x0+l0);
    }
    else {

      // Left end
      xg = bound_left;
      j_x[mod(bound_left,n_g)] += left_end_xg(nl, nr, n0, l0, xg, x0, x0, bound_left+1);
    
      // Pieces connecting two gridpoints
      for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
	j_x[mod(cell,n_g)] += q_unweighted(nl, nr, n0, l0, x0, x0, cell);
	j_x[mod(cell,n_g)] += cover_xg(nl, nr, n0, l0, cell, x0, cell, cell+1);	
      }

      // Right end
      xg = bound_right-1;
    
      j_x[mod(xg,n_g)] += q_unweighted(nl, nr, n0, l0, x0, x0, xg);
      j_x[mod(xg,n_g)] += right_end_xg(nl, nr, n0, l0, xg, x0, xg, x0+l0);
    }
  
    for (int cell = (bound_right); cell < right_max; cell++) {
      j_x[mod(cell,n_g)] += charge;
    }
  }
  
  return;
}

void deposit_charge_to_left_segment_higher_order_0(std::vector<double> &j_x,
						   double x_tracer_a, 
						   double x_tracer_b,
						   double charge,
						   int n_g,
						   double dx,
						   int right_max,
						   double nl,
						   double nr)
{
  double left, right, x0, l0, n0;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx + 0.5;
  right = fmax(x_tracer_a, x_tracer_b) / dx + 0.5;
  x0 = left;
  l0 = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);    

  if (bound_right == (bound_left + 1) or l0==0.0) {
    j_x[mod(bound_left,n_g)] += charge;
  } else {
    n0 = charge / l0;

    if (left == x_tracer_b / dx) {
      std::swap(nl, nr);
    }

    if (nl == 0.0 && nr == 0.0) {
      nl = 1.0;
      nr = 1.0;
    } else if (nl == 0.0) {
      nl = 1.0;
      nr = 0.0;
    } else if (nr == 0.0) {
      nl = 0.0;
      nr = 1.0;
    }
    
    // Left end
    j_x[mod(bound_left,n_g)] += q_unweighted(nl, nr, n0, l0, x0, x0, bound_left+1);
    
    // Pieces connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      j_x[mod(cell,n_g)] += q_unweighted(nl, nr, n0, l0, x0, x0, cell+1);
    }

    // Right end
    j_x[mod((bound_right-1),n_g)] += charge;
  }
  
  for (int cell = (bound_right); cell < right_max; cell++) {
    j_x[mod(cell,n_g)] += charge;
  }
  
  return;
}


void deposit_rho_segment_higher_order_1(std::vector<double> &rho,
					double x_tracer_a, 
					double x_tracer_b,
					double charge,
					int n_g,
					double dx,
					double nl,
					double nr)
{
  double left, right, l0, xg, xa, xb, x0, n0;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  x0 = left;
  l0 = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);

  // If tracers are between two gridpoints
  if (l0 == 0.0) {
    rho[mod(bound_left,n_g)] += charge * (bound_right - left);
    rho[mod(bound_left+1,n_g)] += charge * (1.0 - (bound_right - left));
  } else {
    n0 = charge / l0;

    if (left == x_tracer_b / dx) {
      std::swap(nl, nr);
    }
    
    if (nl == 0.0 && nr == 0.0) {
      nl = 1.0;
      nr = 1.0;
    } else if (nl == 0.0) {
      nl = 1.0;
      nr = 0.0;
    } else if (nr == 0.0) {
      nl = 0.0;
      nr = 1.0;
    }
    
    if (bound_right == (bound_left + 1)) {
      xg = bound_left;
      xa = x0;
      xb = x0+l0;
      rho[mod(bound_left,n_g)] += full_xg(nl, nr, n0, l0, xg, x0, xa, xb);
      rho[mod(bound_left+1,n_g)] += full_xgp1(nl, nr, n0, l0, xg, x0, xa, xb);
    }
    else {
      // Left end
      xg = bound_left;
      xa = x0;
      xb = bound_left+1;
      rho[mod(bound_left,n_g)] += left_end_xg(nl, nr, n0, l0, xg, x0, xa, xb);
      rho[mod(bound_left+1,n_g)] += left_end_xgp1(nl, nr, n0, l0, xg, x0, xa, xb);
      // Pieces connecting two gridpoints
      for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
	xg = cell;
	xa = cell;
	xb = cell+1;
	rho[mod(cell,n_g)] += cover_xg(nl, nr, n0, l0, xg, x0, xa, xb);
	rho[mod(cell+1,n_g)] += cover_xgp1(nl, nr, n0, l0, xg, x0, xa, xb);
      }
      // Right end
      xg = bound_right-1;
      xa = bound_right-1;
      xb = x0+l0;
      rho[mod(bound_right-1,n_g)] += right_end_xg(nl, nr, n0, l0, xg, x0, xa, xb);
      rho[mod(bound_right,n_g)] += right_end_xgp1(nl, nr, n0, l0, xg, x0, xa, xb);
    }
  }
  return;
}

void deposit_rho_segment_higher_order_0(std::vector<double> &rho,
					double x_tracer_a, 
					double x_tracer_b,
					double charge,
					int n_g,
					double dx,
					double nl,
					double nr)
{
  double left, right, l0, xa, xb, x0, n0;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx + 0.5;
  right = fmax(x_tracer_a, x_tracer_b) / dx + 0.5;
  x0 = left;
  l0 = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1) or l0==0.0) {
    rho[mod(bound_left,n_g)] += charge;
  } else {
    n0 = charge / l0;

    if (left == x_tracer_b / dx) {
      std::swap(nl, nr);
    }
    
    if (nl == 0.0 && nr == 0.0) {
      nl = 1.0;
      nr = 1.0;
    } else if (nl == 0.0) {
      nl = 1.0;
      nr = 0.0;
    } else if (nr == 0.0) {
      nl = 0.0;
      nr = 1.0;
    }

    // Left end
    xa = x0;
    xb = bound_left+1;
    rho[mod(bound_left,n_g)] += q_unweighted(nl, nr, n0, l0, x0, xa, xb);

    // Pieces connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      xa = cell;
      xb = cell+1;
      rho[mod(cell,n_g)] += q_unweighted(nl, nr, n0, l0, x0, xa, xb);
    }
    // Right end
    xa = bound_right-1;
    xb = x0+l0;
    rho[mod(bound_right-1,n_g)] += q_unweighted(nl, nr, n0, l0, x0, xa, xb);
  }
  return;
}


void ParticleSpecies::deposit_j_x_sic_higher_order_0(std::vector<double> &j_x)
{
  double right_max;
  
  for (int i = 0; i < n_p; i++) {
    // Index of right bounding gridpoint
    right_max = fmax(fmax(x_old[i],x_old[i+1]),fmax(x[i],x[i+1]));
    right_max = right_max / dx + 0.5; // account for shift in the calculations below
    right_max = ceil(right_max);
    
    deposit_charge_to_left_segment_higher_order_0(j_x, x[i], x[i+1], (-1.0) * (charge[i] * (dx / dt)), n_g, dx, right_max, density[(i+1)-1], density[(i+1)+1]);
    deposit_charge_to_left_segment_higher_order_0(j_x, x_old[i], x_old[i+1], (charge[i] * (dx / dt)), n_g, dx, right_max, density_old[(i+1)-1], density_old[(i+1)+1]);
  }
  return;
}

void ParticleSpecies::deposit_rho_sic_higher_order_0(std::vector<double> &rho)
{
  for (int i = 0; i < n_p; i++) {
    deposit_rho_segment_higher_order_0(rho, x[i], x[i+1], charge[i], n_g, dx, density[(i+1)-1], density[(i+1)+1]);
  }
  return;
}

void ParticleSpecies::deposit_j_x_sic_higher_order_1(std::vector<double> &j_x)
{
  double right_max;
  
  for (int i = 0; i < n_p; i++) {
    // Index of right bounding gridpoint
    right_max = fmax(fmax(x_old[i],x_old[i+1]),fmax(x[i],x[i+1]));
    right_max = right_max / dx;
    right_max = ceil(right_max);
    
    deposit_charge_to_left_segment_higher_order_1(j_x, x[i], x[i+1], (-1.0) * (charge[i] * (dx / dt)), n_g, dx, right_max, density[(i+1)-1], density[(i+1)+1]);
    deposit_charge_to_left_segment_higher_order_1(j_x, x_old[i], x_old[i+1], (charge[i] * (dx / dt)), n_g, dx, right_max, density_old[(i+1)-1], density_old[(i+1)+1]);
  }
  return;
}

void ParticleSpecies::deposit_rho_sic_higher_order_1(std::vector<double> &rho)
{
  for (int i = 0; i < n_p; i++) {
    deposit_rho_segment_higher_order_1(rho, x[i], x[i+1], charge[i], n_g, dx, density[(i+1)-1], density[(i+1)+1]);
  }
  return;
}


// Higher order routines above

double linear_interp(double a, double b, double x0, double l0, double x)
{
  return a + (b - a) * (x - x0) / l0;
}

void deposit_rho_segment_higher_order_center(std::vector<double> &rho,
					     double x_tracer_a, 
					     double x_tracer_b,
					     double charge,
					     int n_g,
					     double dx,
					     double nl,
					     double nr)
{
  double left, right, l0, x0, n0, a, b;
  int bound_left, bound_right, swap;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  x0 = left;
  l0 = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);

  n0 = charge / l0;

  swap = 0;
  if (left == x_tracer_b / dx) {
    swap = 1;
    std::swap(nl, nr);
  }
    
  if (nl == 0.0 && nr == 0.0) {
    nl = 1.0;
    nr = 1.0;
  } else if (nl == 0.0) {
    nl = 1.0;
    nr = 0.0;
  } else if (nr == 0.0) {
    nl = 0.0;
    nr = 1.0;
  }

  a = 2.0 * nl * n0 / (nl + nr);
  b = 2.0 * nr * n0 / (nl + nr);

  if (bound_left==left) {
    if (swap == 0) {
      rho[mod(bound_left,n_g)] += linear_interp(a, b, x0, l0, bound_left);
    }
  }
  for (int cell = (bound_left+1); cell < bound_right; cell++) {
    rho[mod(cell,n_g)] += linear_interp(a, b, x0, l0, cell);
  }
  if (bound_right==right) {
    if (swap == 1) {
      rho[mod(bound_right,n_g)] += linear_interp(a, b, x0, l0, bound_right);
    }
  }

  return;
}

void ParticleSpecies::deposit_rho_sic_higher_order_center(std::vector<double> &rho)
{
  for (int i = 0; i < n_p; i++) {
    deposit_rho_segment_higher_order_center(rho, x[i], x[i+1], charge[i], n_g, dx, density[(i+1)-1], density[(i+1)+1]);
  }
  return;
}

void deposit_j_segment_higher_order_center(std::vector<double> &j,
					   double x_tracer_a, 
					   double x_tracer_b,
					   double charge,
					   int n_g,
					   double dx,
					   double nl,
					   double nr,
					   double vl,
					   double vr)
{
  double left, right, l0, x0, n0, a, b;
  int bound_left, bound_right, swap;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  x0 = left;
  l0 = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);

  n0 = charge / l0;

  swap = 0;
  if (left == x_tracer_b / dx) {
    std::swap(nl, nr);
    std::swap(vl, vr);
    swap = 1;
  }
    
  if (nl == 0.0 && nr == 0.0) {
    nl = 1.0;
    nr = 1.0;
  } else if (nl == 0.0) {
    nl = 1.0;
    nr = 0.0;
  } else if (nr == 0.0) {
    nl = 0.0;
    nr = 1.0;
  }

  a = 2.0 * nl * n0 / (nl + nr);
  b = 2.0 * nr * n0 / (nl + nr);

  if (bound_left==left) {
    if (swap == 0) {
      j[mod(bound_left,n_g)] += linear_interp(a, b, x0, l0, bound_left) * linear_interp(vl, vr, x0, l0, bound_left);
    }
  }
  for (int cell = (bound_left+1); cell < bound_right; cell++) {
    j[mod(cell,n_g)] += linear_interp(a, b, x0, l0, cell) * linear_interp(vl, vr, x0, l0, bound_left);
  }
  if (bound_right==right) {
    if (swap == 1) {
      j[mod(bound_right,n_g)] += linear_interp(a, b, x0, l0, bound_right) * linear_interp(vl, vr, x0, l0, bound_left);
    }
  }

  return;
}


void ParticleSpecies::deposit_j_x_sic_higher_order_center(std::vector<double> &j_x)
{
  double xl, xr, vl, vr, nl, nr;
  for (int i = 0; i < n_p; i++) {
    xl = (x_old[i] + x[i]) / 2.0 + dx / 2.0;
    xr = (x_old[i+1] + x[i+1]) / 2.0 + dx / 2.0;
    vl = u_x[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    vr = u_x[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    nl = density_tavg[(i+1)-1];
    nr = density_tavg[(i+1)+1];
    deposit_j_segment_higher_order_center(j_x, xl, xr, charge[i], n_g, dx, nl, nr, vl, vr);
  }
  return;
}

void ParticleSpecies::deposit_j_y_sic_higher_order_center(std::vector<double> &j_y)
{
  double xl, xr, vl, vr, nl, nr;
  for (int i = 0; i < n_p; i++) {
    xl = (x_old[i] + x[i]) / 2.0;
    xr = (x_old[i+1] + x[i+1]) / 2.0;
    vl = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    vr = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    nl = density_tavg[(i+1)-1];
    nr = density_tavg[(i+1)+1];
    deposit_j_segment_higher_order_center(j_y, xl, xr, charge[i], n_g, dx, nl, nr, vl, vr);
  }
  return;
}

void ParticleSpecies::deposit_j_z_sic_higher_order_center(std::vector<double> &j_z)
{
  double xl, xr, vl, vr, nl, nr;
  for (int i = 0; i < n_p; i++) {
    xl = (x_old[i] + x[i]) / 2.0;
    xr = (x_old[i+1] + x[i+1]) / 2.0;
    vl = u_z[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    vr = u_z[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    nl = density_tavg[(i+1)-1];
    nr = density_tavg[(i+1)+1];
    deposit_j_segment_higher_order_center(j_z, xl, xr, charge[i], n_g, dx, nl, nr, vl, vr);
  }
  return;
}

// Higher order center routines above

void ParticleSpecies::deposit_rho_sic_0(std::vector<double> &rho)
{
  for (int i = 0; i < n_p; i++) {
    deposit_rho_sic_0_segment(rho, x[i], x[i+1], charge[i], n_g, dx);
  }
  return;
}

void ParticleSpecies::deposit_rho_sic_1(std::vector<double> &rho)
{
  for (int i = 0; i < n_p; i++) {    
    deposit_rho_sic_1_segment(rho, x[i], x[i+1], charge[i], n_g, dx);
  }
  return;
}


void ParticleSpecies::deposit_j_x_sic_0(std::vector<double> &j_x)
{
  double left_initial, right_initial, left_final, right_final, length_initial, length_final,
    charge_initial, charge_final, cell_boundary;
  int bound_left, bound_right;
  for (int i = 0; i < n_p; i++) {
    // Initial line segment
    left_initial = fmin(x_old[i], x_old[i+1]) / dx;
    right_initial = fmax(x_old[i], x_old[i+1]) / dx;
    length_initial = right_initial - left_initial;
    // Final line segment
    left_final = fmin(x[i], x[i+1]) / dx;
    right_final = fmax(x[i], x[i+1]) / dx;
    length_final = right_final - left_final;
    // Bounding box
    bound_left = floor(fmin(left_initial,left_final));
    bound_right = ceil(fmax(right_initial,right_final));
    // Loop over cell boundaries that could be affected
    for (int cell = bound_left; cell < bound_right; cell++) {
      // Calculate charge initially to the left of the current cell's right boundary
      cell_boundary = cell + 0.5;
      if (right_initial <= cell_boundary) {
	charge_initial = 1.0;
      } 
      else if (left_initial >= cell_boundary) {
	charge_initial = 0.0;
      }
      else {
	charge_initial = (cell_boundary - left_initial) / length_initial;
      }
      // Calculate charge finally to the left of the current cell's right boundary
      if (right_final <= cell_boundary) {
	charge_final = 1.0;
      } 
      else if (left_final >= cell_boundary) {
	charge_final = 0.0;
      }
      else {
	charge_final = (cell_boundary - left_final) / length_final;
      }
      j_x[mod(cell,n_g)] += charge[i] * (charge_initial - charge_final) * (dx / dt);
    }
  }
  return;
}

void ParticleSpecies::deposit_rho_pic_1(std::vector<double> &rho)
{
  double lower_weight, x_normalized;
  int lower;
  for (int i = 0; i < n_p; i++) {
    x_normalized = x[i] / dx;
    lower = floor(x_normalized);
    lower_weight = double(lower) + 1.0 - x_normalized;
    rho[mod(floor(lower),n_g)] += charge[i] * lower_weight;
    rho[mod(floor(lower+1),n_g)] += charge[i] * (1.0 - lower_weight);
  }
  return;
}

void ParticleSpecies::deposit_rho_pic_0(std::vector<double> &rho)
{
  double x_normalized;
  int nearest_gridpoint;
  for (int i = 0; i < n_p; i++) {
    x_normalized = x[i] / dx;
    nearest_gridpoint = get_nearest_gridpoint(x_normalized);
    rho[mod(nearest_gridpoint,n_g)] += charge[i];
  }
  return;
}

void ParticleSpecies::deposit_j_y_pic_0(std::vector<double> &j_y)
{
  double x_tavg, gamma;
  int nearest_gridpoint;
  for (int i = 0; i < n_p; i++) {
    x_tavg = (x[i] + x_old[i]) / 2.0;
    x_tavg = x_tavg / dx;
    nearest_gridpoint = get_nearest_gridpoint(x_tavg);
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_y[mod(nearest_gridpoint,n_g)] += charge[i] * u_y[i] / gamma;
  }
  return;
}

void ParticleSpecies::deposit_j_z_pic_0(std::vector<double> &j_z)
{
  double x_tavg, gamma;
  int nearest_gridpoint;
  for (int i = 0; i < n_p; i++) {
    x_tavg = (x[i] + x_old[i]) / 2.0;
    x_tavg = x_tavg / dx;
    nearest_gridpoint = get_nearest_gridpoint(x_tavg);
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_z[mod(nearest_gridpoint,n_g)] += charge[i] * u_z[i] / gamma;
  }
  return;
}


double interpolate_segment_velocity(double v_left, double v_right, 
				    double length, double distance_from_left)
{
  return v_left + (distance_from_left / length) * (v_right - v_left);
}

void deposit_j_t_segment_zero(std::vector<double> &j_t, double left, 
			      double right, double v_t_left, 
			      double v_t_right, 
			      double charge, int n_g, double dx)
{
  double length, charge_fraction, avg_velocity, midpoint, distance_from_left;
  int bound_left, bound_right;
  left = left / dx;
  right = right / dx;
  length = right - left;

  // Shift tracer positions so that 0 is the left bounary of cell zero
  left = left + 0.5;
  right = right + 0.5;

  bound_left = floor(left);
  bound_right = ceil(right);

  if (bound_right == (bound_left + 1)) {
    // If tracers are between two gridpoints
    charge_fraction = 1.0;
    avg_velocity = (v_t_left + v_t_right) / 2.0;
    j_t[mod(bound_left,n_g)] += charge * charge_fraction * avg_velocity;
  }
  else {
    // Left end
    charge_fraction = (bound_left + 1.0 - left) / length;
    midpoint = (bound_left + 1.0 + left) / 2.0;
    distance_from_left = left - midpoint;
    avg_velocity = interpolate_segment_velocity(v_t_left, v_t_right, length, distance_from_left);
    j_t[mod(bound_left,n_g)] += charge * charge_fraction * avg_velocity;
    // Pieces connecting two gridpoints
    charge_fraction = (1.0) / length;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      distance_from_left = cell + 0.5 - left;
      avg_velocity = interpolate_segment_velocity(v_t_left, v_t_right, length, distance_from_left);
      j_t[mod(cell,n_g)] += charge * charge_fraction * avg_velocity;
    }
    // Right end
    charge_fraction = (right - (bound_right - 1)) / length;
    midpoint = (right + bound_right) / 2.0;
    distance_from_left = midpoint - left;
    avg_velocity = interpolate_segment_velocity(v_t_left, v_t_right, length, distance_from_left);
    j_t[mod((bound_right-1),n_g)] += charge * charge_fraction * avg_velocity;
  }
  return;
}

void ParticleSpecies::deposit_j_y_sic_0(std::vector<double> &j_y)
{
  double x_i_tavg, x_ip1_tavg, left, right, v_y_left, v_y_right;
  for (int i = 0; i < n_p; i++) {
    x_i_tavg = (x_old[i] + x[i]) / 2.0;
    x_ip1_tavg = (x_old[i+1] + x[i+1]) / 2.0;
    if (x_ip1_tavg >= x_i_tavg) {
      left = x_i_tavg;
      right = x_ip1_tavg;
      v_y_left = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
      v_y_right = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    }
    else {
      left = x_ip1_tavg;
      right = x_i_tavg;
      v_y_left = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
      v_y_right = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    }
    deposit_j_t_segment_zero(j_y, left, right, v_y_left, v_y_right, charge[i], n_g, dx);
  }
  return;
}

void ParticleSpecies::deposit_j_z_sic_0(std::vector<double> &j_z)
{
  double x_i_tavg, x_ip1_tavg, left, right, v_z_left, v_z_right;
  for (int i = 0; i < n_p; i++) {
    x_i_tavg = (x_old[i] + x[i]) / 2.0;
    x_ip1_tavg = (x_old[i+1] + x[i+1]) / 2.0;
    if (x_ip1_tavg >= x_i_tavg) {
      left = x_i_tavg;
      right = x_ip1_tavg;
      v_z_left = u_z[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
      v_z_right = u_z[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    }
    else {
      left = x_ip1_tavg;
      right = x_i_tavg;
      v_z_left = u_z[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
      v_z_right = u_z[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    }
    deposit_j_t_segment_zero(j_z, left, right, v_z_left, v_z_right, charge[i], n_g, dx);
  }
  return;
}

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

void deposit_j_t_segment_linear(std::vector<double> &j_t, double xl, double xr, 
				double vl, double vr, double charge, int n_g,
				double dx)
{
  double length, xa, xb;
  int bound_left, bound_right;
  xl = xl / dx;
  xr = xr / dx;
  length = xr - xl;

  bound_left = floor(xl);
  bound_right = ceil(xr);

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    if (length==0.0) {
      j_t[mod(bound_left,n_g)] += charge * ((vl+vr)/2.0) * (ceil(xr)-xl);
      j_t[mod((bound_left+1),n_g)] += charge * ((vl+vr)/2.0) * (xl-floor(xl));
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
    xb = double(bound_left+1);
    j_t[mod(bound_left,n_g)] += j_t_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
    j_t[mod((bound_left+1),n_g)] += j_t_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
    // Portions connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      xa = double(cell);
      xb = double(cell+1);
      j_t[mod(cell,n_g)] += j_t_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
      j_t[mod((cell+1),n_g)] += j_t_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
    }
    // Right end
    xa = double(bound_right-1);
    xb = xr;
    j_t[mod((bound_right-1),n_g)] += j_t_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
    j_t[mod(bound_right,n_g)] += j_t_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
  }
  return;
}

void ParticleSpecies::deposit_j_y_sic_1(std::vector<double> &j_y)
{
  double x_i_tavg, x_ip1_tavg, left, right, v_y_left, v_y_right;
  for (int i = 0; i < n_p; i++) {
    x_i_tavg = (x_old[i] + x[i]) / 2.0;
    x_ip1_tavg = (x_old[i+1] + x[i+1]) / 2.0;
    if (x_ip1_tavg >= x_i_tavg) {
      left = x_i_tavg;
      right = x_ip1_tavg;
      v_y_left = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
      v_y_right = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    }
    else {
      left = x_ip1_tavg;
      right = x_i_tavg;
      v_y_left = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
      v_y_right = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    }
    deposit_j_t_segment_linear(j_y, left, right, v_y_left, v_y_right, charge[i], n_g, dx);
  }
  return;
}

void ParticleSpecies::deposit_j_z_sic_1(std::vector<double> &j_z)
{
  double x_i_tavg, x_ip1_tavg, left, right, v_z_left, v_z_right;
  for (int i = 0; i < n_p; i++) {
    x_i_tavg = (x_old[i] + x[i]) / 2.0;
    x_ip1_tavg = (x_old[i+1] + x[i+1]) / 2.0;
    if (x_ip1_tavg >= x_i_tavg) {
      left = x_i_tavg;
      right = x_ip1_tavg;
      v_z_left = u_z[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
      v_z_right = u_z[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    }
    else {
      left = x_ip1_tavg;
      right = x_i_tavg;
      v_z_left = u_z[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
      v_z_right = u_z[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    }
    deposit_j_t_segment_linear(j_z, left, right, v_z_left, v_z_right, charge[i], n_g, dx);
  }
  return;
}

void deposit_j_t_segment_higher_order_0(std::vector<double> &j_t,
					double x_tracer_a, 
					double x_tracer_b,
					double charge,
					int n_g,
					double dx,
					double nl,
					double nr,
					double vl,
					double vr)
{
  double left, right, l0, xa, xb, x0, n0;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx + 0.5;
  right = fmax(x_tracer_a, x_tracer_b) / dx + 0.5;
  x0 = left;
  l0 = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);

  if (left == x_tracer_b / dx) {
    std::swap(nl, nr);
    std::swap(vl, vr);
  }
    
  if (nl == 0.0 && nr == 0.0) {
    nl = 1.0;
    nr = 1.0;
  } else if (nl == 0.0) {
    nl = 1.0;
    nr = 0.0;
  } else if (nr == 0.0) {
    nl = 0.0;
    nr = 1.0;
  }
    
  if (bound_right == (bound_left + 1)) {
    // If tracers are between two gridpoints
    xa = x0;
    xb = x0+l0;
    j_t[mod(bound_left,n_g)] += j_t_unweighted_full(nl, nr, vl, vr, charge);
  } else {
    n0 = charge / l0;
    // Left end
    xa = x0;
    xb = bound_left+1;
    j_t[mod(bound_left,n_g)] += j_t_unweighted(nl, nr, vl, vr, n0, l0, x0, xa, xb);
    // Pieces connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      xa = cell;
      xb = cell+1;
      j_t[mod(cell,n_g)] += j_t_unweighted(nl, nr, vl, vr, n0, l0, x0, xa, xb);
    }
    // Right end
    xa = bound_right-1;
    xb = x0+l0;
    j_t[mod(bound_right-1,n_g)] += j_t_unweighted(nl, nr, vl, vr, n0, l0, x0, xa, xb);
  }

  return;
}

void deposit_j_t_segment_higher_order_1(std::vector<double> &j_t,
					double x_tracer_a, 
					double x_tracer_b,
					double charge,
					int n_g,
					double dx,
					double nl,
					double nr,
					double vl,
					double vr)
{
  double left, right, l0, xg, xa, xb, x0, n0;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  x0 = left;
  l0 = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);

  if (left == x_tracer_b / dx) {
    std::swap(nl, nr);
    std::swap(vl, vr);
  }
    
  if (nl == 0.0 && nr == 0.0) {
    nl = 1.0;
    nr = 1.0;
  } else if (nl == 0.0) {
    nl = 1.0;
    nr = 0.0;
  } else if (nr == 0.0) {
    nl = 0.0;
    nr = 1.0;
  }
    
  if (bound_right == (bound_left + 1)) {
    // If tracers are between two gridpoints
    xg = bound_left;
    xa = x0;
    xb = x0+l0;
    j_t[mod(bound_left,n_g)] += full_j_t_xg(nl, nr, vl, vr, charge, l0, xg, x0);
    j_t[mod(bound_left+1,n_g)] += full_j_t_xgp1(nl, nr, vl, vr, charge, l0, xg, x0);
  } else {
    n0 = charge / l0;
    // Left end
    xg = bound_left;
    xa = x0;
    xb = bound_left+1;
    j_t[mod(bound_left,n_g)] += left_end_j_t_xg(nl, nr, vl, vr, n0, l0, xg, x0, xa, xb);
    j_t[mod(bound_left+1,n_g)] += left_end_j_t_xgp1(nl, nr, vl, vr, n0, l0, xg, x0, xa, xb);
    // Pieces connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      xg = cell;
      xa = cell;
      xb = cell+1;
      j_t[mod(cell,n_g)] += cover_j_t_xg(nl, nr, vl, vr, n0, l0, xg, x0, xa, xb);
      j_t[mod(cell+1,n_g)] += cover_j_t_xgp1(nl, nr, vl, vr, n0, l0, xg, x0, xa, xb);
    }
    // Right end
    xg = bound_right-1;
    xa = bound_right-1;
    xb = x0+l0;
    j_t[mod(bound_right-1,n_g)] += right_end_j_t_xg(nl, nr, vl, vr, n0, l0, xg, x0, xa, xb);
    j_t[mod(bound_right,n_g)] += right_end_j_t_xgp1(nl, nr, vl, vr, n0, l0, xg, x0, xa, xb);
  }

  return;
}


void ParticleSpecies::deposit_j_y_sic_higher_order_0(std::vector<double> &j_y)
{
  double xl, xr, vl, vr, nl, nr;
  for (int i = 0; i < n_p; i++) {
    xl = (x_old[i] + x[i]) / 2.0;
    xr = (x_old[i+1] + x[i+1]) / 2.0;
    vl = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    vr = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    nl = density_tavg[(i+1)-1];
    nr = density_tavg[(i+1)+1];
    deposit_j_t_segment_higher_order_0(j_y, xl, xr, charge[i], n_g, dx, nl, nr, vl, vr);
  }
  return;
}

void ParticleSpecies::deposit_j_z_sic_higher_order_0(std::vector<double> &j_z)
{
  double xl, xr, vl, vr, nl, nr;
  for (int i = 0; i < n_p; i++) {
    xl = (x_old[i] + x[i]) / 2.0;
    xr = (x_old[i+1] + x[i+1]) / 2.0;
    vl = u_z[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    vr = u_z[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    nl = density_tavg[(i+1)-1];
    nr = density_tavg[(i+1)+1];
    deposit_j_t_segment_higher_order_0(j_z, xl, xr, charge[i], n_g, dx, nl, nr, vl, vr);
  }
  return;
}

void ParticleSpecies::deposit_j_y_sic_higher_order_1(std::vector<double> &j_y)
{
  double xl, xr, vl, vr, nl, nr;
  for (int i = 0; i < n_p; i++) {
    xl = (x_old[i] + x[i]) / 2.0;
    xr = (x_old[i+1] + x[i+1]) / 2.0;
    vl = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    vr = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    nl = density_tavg[(i+1)-1];
    nr = density_tavg[(i+1)+1];
    deposit_j_t_segment_higher_order_1(j_y, xl, xr, charge[i], n_g, dx, nl, nr, vl, vr);
  }
  return;
}

void ParticleSpecies::deposit_j_z_sic_higher_order_1(std::vector<double> &j_z)
{
  double xl, xr, vl, vr, nl, nr;
  for (int i = 0; i < n_p; i++) {
    xl = (x_old[i] + x[i]) / 2.0;
    xr = (x_old[i+1] + x[i+1]) / 2.0;
    vl = u_z[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    vr = u_z[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0) + pow(u_z[i+1], 2.0));
    nl = density_tavg[(i+1)-1];
    nr = density_tavg[(i+1)+1];
    deposit_j_t_segment_higher_order_1(j_z, xl, xr, charge[i], n_g, dx, nl, nr, vl, vr);
  }
  return;
}

void ParticleSpecies::deposit_j_x_pic_0(std::vector<double> &j_x)
{
  double x_initial, x_final;
  int ngp_initial, ngp_final;
  for (int i = 0; i < n_p; i++) {
    x_initial = x_old[i] / dx;
    x_final = x[i] / dx;
    ngp_initial = get_nearest_gridpoint(x_initial);
    ngp_final = get_nearest_gridpoint(x_final);
    if (ngp_initial < ngp_final) {
      j_x[mod(ngp_initial,n_g)] += charge[i] * (dx / dt);
    } else if (ngp_initial > ngp_final) {
      j_x[mod(ngp_final,n_g)] -= charge[i] * (dx / dt);
    }
  }
  return;
}

void ParticleSpecies::deposit_j_x_pic_1(std::vector<double> &j_x)
{
  double x_initial, x_final;
  for (int i = 0; i < n_p; i++) {
      x_initial = x_old[i] / dx;
      x_final = x[i] / dx;
      if (floor(x_initial)==floor(x_final)) {
	j_x[mod(floor(x_initial),n_g)] += charge[i] * (x_final - x_initial) * (dx / dt);
      } 
      else if (floor(x_initial) < floor(x_final)) {
	j_x[mod(floor(x_initial),n_g)] += charge[i] * (1.0 - (x_initial - floor(x_initial))) * (dx / dt);
	j_x[mod(floor(x_final),n_g)] += charge[i] * (x_final - floor(x_final)) * (dx / dt);
      }
      else {
	j_x[mod(floor(x_initial),n_g)] += charge[i] * (floor(x_initial) - x_initial) * (dx / dt);
	j_x[mod(floor(x_final),n_g)] += charge[i] * ((x_final - floor(x_final)) - 1.0) * (dx / dt);
      }
  }
  return;
}

void ParticleSpecies::deposit_j_y_pic_1(std::vector<double> &j_y)
{
  int i_lower;
  double alpha, beta, x_tavg, gamma;
  for (int i = 0; i < n_p; i++) {
    x_tavg = (x[i] + x_old[i]) / 2.0;
    x_tavg = x_tavg / dx;
    i_lower = floor(x_tavg);
    alpha = 1.0 - (x_tavg - i_lower);
    beta = 1.0 - alpha;

    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_y[mod(i_lower,n_g)] += alpha * charge[i] * u_y[i] / gamma;
    j_y[mod((i_lower+1),n_g)] += beta * charge[i] * u_y[i] / gamma;
  }
  return;
}

void ParticleSpecies::deposit_j_z_pic_1(std::vector<double> &j_z)
{
  int i_lower;
  double alpha, beta, x_tavg, gamma;
  for (int i = 0; i < n_p; i++) {
    x_tavg = (x[i] + x_old[i]) / 2.0;
    x_tavg = x_tavg / dx;
    i_lower = floor(x_tavg);
    alpha = 1.0 - (x_tavg - i_lower);
    beta = 1.0 - alpha;

    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0) + pow(u_z[i], 2.0));
    j_z[mod(i_lower,n_g)] += alpha * charge[i] * u_z[i] / gamma;
    j_z[mod((i_lower+1),n_g)] += beta * charge[i] * u_z[i] / gamma;
  }
  return;
}


double put_in_box(double x, double box_length)
{
  if (x < 0.0) {
    x = x + box_length;
  } 
  else if (x >= box_length) {
    x = x - box_length;
  }
  return x;
}


void ParticleSpecies::write_phase(int species_number, int t, int my_rank)
{
  std::string x_filename;
  std::string u_x_filename;
  
  std::stringstream ss;
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
    x[i] = x[i] + dt * u_x[i] / gamma;
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

void ParticleSpecies::split_segment_lagrange_3(int i)
{
  double new_id = (lagrangian_id[i+1]+lagrangian_id[i])/2.0;
  x.insert(x.begin()+i+1, lagrange_3(&lagrangian_id[i], &x[i], new_id));
  // Linear interpolation on x_old to be consistent with Gauss's law
  x_old.insert(x_old.begin()+i+1, (x_old[i+1] + x_old[i])/2.0);
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
  double new_id = (lagrangian_id[i+1]+lagrangian_id[i])/2.0;
  x.insert(x.begin()+i+1, (x[i+1] + x[i])/2.0);
  x_old.insert(x_old.begin()+i+1, (x_old[i+1] + x_old[i])/2.0);
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
    length = fabs(x[i+1] - x[i]);
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

void ParticleSpecies::communicate_ghost_particles(MPI_Comm COMM)
{
  MPI_Status status;
  int num_procs, my_rank, dest, source;
  int tag = 0;
  
  MPI_Comm_size(COMM, &num_procs);
  MPI_Comm_rank(COMM, &my_rank);
  dest = mod(my_rank-1, num_procs);
  source = mod(my_rank+1, num_procs);
  
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
  
  if (my_rank==(num_procs-1)) {
    x[n_p] += n_g * dx;
    x[n_p+1] += n_g * dx;
    x_old[n_p] += n_g * dx;
    x_old[n_p+1] += n_g * dx;
  }
  
  return;
}

void ParticleSpecies::calculate_segment_density(MPI_Comm COMM)
{
  MPI_Status status;
  int num_procs, my_rank, dest, source;
  int tag = 0;

  MPI_Comm_size(COMM, &num_procs);
  MPI_Comm_rank(COMM, &my_rank);

  density.resize(n_p+2);
  density_old.resize(n_p+2);
  density_tavg.resize(n_p+2);  

  for (int i = 0; i < n_p; i++) {
    density[i+1] = charge[i] / fabs(x[i+1]-x[i]);
    if (fabs(x[i+1]-x[i]) == 0.0 or std::isinf(density[i+1])) {
      density[i+1] = 0.0;
    }
    
    density_old[i+1] = charge[i] / fabs(x_old[i+1]-x_old[i]);
    if (fabs(x_old[i+1]-x_old[i]) == 0.0 or std::isinf(density_old[i+1])) {
      density_old[i+1] = 0.0;
    }

    density_tavg[i+1] = 2.0 * charge[i] / fabs(x[i+1]+x_old[i+1]-x[i]-x_old[i]);
    if (fabs(x[i+1]+x_old[i+1]-x[i]-x_old[i]) == 0.0 or std::isinf(density_tavg[i+1])) {
      density_tavg[i+1] = 0.0;
    }
  }
  
  dest = mod(my_rank+1, num_procs);
  source = mod(my_rank-1, num_procs);

  MPI_Sendrecv(&density[n_p], 1, MPI_DOUBLE, dest, tag,
	       &density[0], 1, MPI_DOUBLE,
	       source, tag, COMM, &status);
  MPI_Sendrecv(&density_old[n_p], 1, MPI_DOUBLE, dest, tag,
	       &density_old[0], 1, MPI_DOUBLE,
	       source, tag, COMM, &status);
  MPI_Sendrecv(&density_tavg[n_p], 1, MPI_DOUBLE, dest, tag,
	       &density_tavg[0], 1, MPI_DOUBLE,
	       source, tag, COMM, &status);

  dest = mod(my_rank-1, num_procs);
  source = mod(my_rank+1, num_procs);

  MPI_Sendrecv(&density[1], 1, MPI_DOUBLE, dest, tag,
	       &density[n_p+1], 1, MPI_DOUBLE,
	       source, tag, COMM, &status);
  MPI_Sendrecv(&density_old[1], 1, MPI_DOUBLE, dest, tag,
	       &density_old[n_p+1], 1, MPI_DOUBLE,
	       source, tag, COMM, &status);  
  MPI_Sendrecv(&density_tavg[1], 1, MPI_DOUBLE, dest, tag,
	       &density_tavg[n_p+1], 1, MPI_DOUBLE,
	       source, tag, COMM, &status);  

  
  return;
}

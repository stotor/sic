#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "mpi.h"

#include "particlespecies.hpp"
#include "utilities.hpp"

////////////////////////////////////////////////////////////////
// Particles

// GRADIENT CALCULATIONS DO NOT YET WORK WITH REFINEMENT

void ParticleSpecies::initialize_species(int species_number, 
					 long long n_ppc, 
					 double u_x_drift, 
					 double u_y_drift, 
					 int mode, 
					 double u_x_1, 
					 double u_y_1,
					 int my_rank,
					 int num_procs)
{
  long long i_start, i_end;
  i_start = n_p * my_rank;
  i_end = i_start + n_p;
  
  std::stringstream ss;
  ss << "particles_";
  ss << species_number;
  species_name = ss.str();

  double k = 2.0 * PI * double(mode) / (n_g * dx);
  double particle_spacing = dx / double(n_ppc);

  // Add ghost tracer particles if using line segments
  if ((method==2)||(method==3)||(method==4)) {
    for (int i = 0; i < 2; i++) {
      charge.push_back(0.0);
      u_x.push_back(0.0);
      u_y.push_back(0.0);
      x.push_back(0.0);
      x_old.push_back(0);
      lagrangian_id.push_back(0.0);
    }
  }
  
  for (long long i = i_start; i < i_end; i++) {
    charge[i-i_start] = (-1.0) * (1.0 / n_ppc);
    x[i-i_start] = ((long double) i) * particle_spacing + particle_spacing / 2.0;
    u_x[i-i_start] = u_x_drift + u_x_1 * sin(k * x[i-i_start]);
    u_y[i-i_start] = u_y_drift + u_y_1 * sin(k * x[i-i_start]);
    lagrangian_id[i-i_start] = i - i_start;
  }

  return;
}

void ParticleSpecies::initialize_beat_heating(int mode_1, int mode_2,
					      double phase_1, double phase_2,
					      double vel_amp)
{
  double k1 = 2.0 * PI * mode_1 / (n_g * dx);
  double k2 = 2.0 * PI * mode_2 / (n_g * dx);

  for (long long i = 0; i < n_p; i++) {
    u_y[i] += vel_amp * sin(k1 * x[i] + phase_1) + (-1.0) * vel_amp * sin(k2 * x[i] + phase_2);
  }
  return;
}


void ParticleSpecies::u_x_perturbation(double amplitude, int mode_max)
{
  double k, phase;
  srand(1);
  for (int mode = 1; mode <= mode_max; mode++) {  
    k = 2.0 * PI * mode / (n_g * dx);
    phase = 2.0 * PI * random_double();
    for (long long i = 0; i < n_p; i++) {
      u_x[i] += amplitude * cos(k * x[i] + phase);
    }
  }
  return;
}


double j_y_segment_linear_left_gridpoint(double charge, double xl, 
					 double length, double vl, double vr, 
					 double xa, double xb)
{
  double xip1 = ceil(xb);
  double j_y = (1.0/(6.0*length*length))*charge*(xa-xb);
  j_y = j_y * (3.0*length*vl*(xa+xb-2.0*xip1)-(vl-vr)*(2.0*xa*xa+2.0*xb*xb+6.0*xip1*xl-3.0*xb*(xip1+xl)+xa*(2.0*xb-3.0*(xip1+xl))));
  return j_y;
}

double j_y_segment_linear_right_gridpoint(double charge, double xl, 
					  double length, double vl, double vr, 
					  double xa, double xb)
{
  double xi = floor(xa);
  double j_y = (1.0/(6.0*length*length))*charge*(xa-xb);
  j_y = j_y * (-3.0*length*vl*(xa+xb-2.0*xi)+(vl-vr)*(2.0*xa*xa+2.0*xb*xb+6.0*xi*xl-3.0*xb*(xi+xl)+xa*(2.0*xb-3.0*(xi+xl))));
  return j_y;
}

void deposit_j_y_segment_linear(std::vector<double> &j_y, double xl, double xr, 
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
      j_y[mod(bound_left,n_g)] += charge * ((vl+vr)/2.0) * (ceil(xr)-xl);
      j_y[mod((bound_left+1),n_g)] += charge * ((vl+vr)/2.0) * (xl-floor(xl));
    } else { 
      xa = xl;
      xb = xr;
      j_y[mod(bound_left,n_g)] += j_y_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
      j_y[mod((bound_left+1),n_g)] += j_y_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
    }
  }
  else {
    // Left end
    xa = xl;
    xb = double(bound_left+1);
    j_y[mod(bound_left,n_g)] += j_y_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
    j_y[mod((bound_left+1),n_g)] += j_y_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
    // Portions connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      xa = double(cell);
      xb = double(cell+1);
      j_y[mod(cell,n_g)] += j_y_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
      j_y[mod((cell+1),n_g)] += j_y_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
    }
    // Right end
    xa = double(bound_right-1);
    xb = xr;
    j_y[mod((bound_right-1),n_g)] += j_y_segment_linear_left_gridpoint(charge, xl, length, vl, vr, xa, xb);
    j_y[mod(bound_right,n_g)] += j_y_segment_linear_right_gridpoint(charge, xl, length, vl, vr, xa, xb);
  }
  return;
}

void ParticleSpecies::deposit_j_y_segments_linear(std::vector<double> &j_y)
{
  double x_i_tavg, x_ip1_tavg, left, right, v_y_left, v_y_right;
  for (int i = 0; i < n_p; i++) {
    x_i_tavg = (x_old[i] + x[i]) / 2.0;
    x_ip1_tavg = (x_old[i+1] + x[i+1]) / 2.0;
    if (x_ip1_tavg >= x_i_tavg) {
      left = x_i_tavg;
      right = x_ip1_tavg;
      v_y_left = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
      v_y_right = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0));
    }
    else {
      left = x_ip1_tavg;
      right = x_i_tavg;
      v_y_left = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0));
      v_y_right = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
    }
    deposit_j_y_segment_linear(j_y, left, right, v_y_left, v_y_right, charge[i], n_g, dx);
  }
  return;
}

void ParticleSpecies::write_particle_diagnostics(int n_t, int my_rank, MPI_Comm COMM)
{
  sum_array_to_root(&energy_history[0], n_t, COMM, my_rank);
  sum_array_to_root(&momentum_x_history[0], n_t, COMM, my_rank);
  sum_array_to_root(&momentum_y_history[0], n_t, COMM, my_rank);
  sum_array_to_root(&n_p_history[0], n_t, COMM, my_rank);

  if (my_rank==0) {
    data_to_file(energy_history, (species_name+"_ene"));
    data_to_file(momentum_x_history, (species_name+"_p_x"));
    data_to_file(momentum_y_history, (species_name+"_p_y"));
    data_to_file(n_p_history, (species_name+"_n_p"));
  }

  return; 
}

void deposit_rho_segment_zero_new(std::vector<double> &rho, double x_tracer_a, 
				double x_tracer_b, double charge, int n_g, double dx)
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

double linear_correction(double gradient, double midpoint,
			 double x_a, double x_b)
{
  return gradient * 0.5 * ((x_b - midpoint) * (x_b - midpoint) - (x_a - midpoint) * (x_a - midpoint));
}

void deposit_rho_single_gradient(std::vector<double> &rho, double x_tracer_a, 
				 double x_tracer_b, double charge, int n_g,
				 double dx, double gradient)
{
  double left, right, length, charge_fraction;
  double midpoint;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  length = right - left;
  bound_left = floor(left);
  bound_right = ceil(right);
  midpoint = (right + left) / 2.0;

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    if ((left == bound_left + 0.5) && (right == bound_left + 0.5)) {
      rho[mod((bound_right),n_g)] += charge;
    } else if (left >= (bound_left + 0.5)) {
      rho[mod((bound_right),n_g)] += charge;
    } else if (right <= (bound_left + 0.5)) {
      rho[mod((bound_left),n_g)] += charge;
    } else {
      rho[mod((bound_left),n_g)] += (bound_left + 0.5 - left) / length * charge + linear_correction(gradient, midpoint, left, (bound_left+0.5))*dx*dx;
      rho[mod((bound_right),n_g)] += (right - (bound_left + 0.5)) / length * charge + linear_correction(gradient, midpoint, (bound_left+0.5), right)*dx*dx;
    }
  }
  else {
    // Left end
    if (left < bound_left + 0.5) {
      rho[mod(bound_left,n_g)] += charge * (bound_left + 0.5 - left) / length + linear_correction(gradient, midpoint, left, (bound_left+0.5))*dx*dx;
      rho[mod((bound_left+1),n_g)] += charge * (0.5) / length + linear_correction(gradient, midpoint, bound_left+0.5, bound_left+1)*dx*dx;
    } else { 
      rho[mod(bound_left+1,n_g)] += charge * (bound_left + 1.0 - left) / length + linear_correction(gradient, midpoint, left, bound_left+1.0)*dx*dx;
    }

    // Pieces connecting two gridpoints
    charge_fraction = (1.0) / length;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      rho[mod(cell,n_g)] += charge * charge_fraction * 0.5 + linear_correction(gradient, midpoint, cell, cell+0.5)*dx*dx;
      rho[mod((cell+1),n_g)] += charge * charge_fraction * 0.5 + linear_correction(gradient, midpoint, cell+0.5, cell+1)*dx*dx;
    }
    
    // Right end
    if (right > bound_right - 0.5) {
      rho[mod(bound_right-1,n_g)] += charge * (0.5) / length + linear_correction(gradient, midpoint, bound_right-1, bound_right-0.5)*dx*dx;
      rho[mod(bound_right,n_g)] += charge * (right - (bound_right-0.5)) / length + linear_correction(gradient, midpoint, bound_right-0.5, right)*dx*dx;
    } else { 
      rho[mod(bound_right-1,n_g)] += charge * (right - (bound_right-1.0)) / length + linear_correction(gradient, midpoint, bound_right-1.0, right)*dx*dx;
    }
  }
  return;
}


void deposit_rho_segment_linear(std::vector<double> &rho, double x_tracer_a, 
				double x_tracer_b, double charge, int n_g, double dx)
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

double interpolate_segment_velocity(double v_left, double v_right, 
				    double length, double distance_from_left)
{
  return v_left + (distance_from_left / length) * (v_right - v_left);
}

void deposit_j_y_segment_zero(std::vector<double> &j_y, double left, 
			      double right, double v_y_left, 
			      double v_y_right, 
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
    avg_velocity = (v_y_left + v_y_right) / 2.0;
    j_y[mod(bound_left,n_g)] += charge * charge_fraction * avg_velocity;
  }
  else {
    // Left end
    charge_fraction = (bound_left + 1.0 - left) / length;
    midpoint = (bound_left + 1.0 + left) / 2.0;
    distance_from_left = left - midpoint;
    avg_velocity = interpolate_segment_velocity(v_y_left, v_y_right, length, distance_from_left);
    j_y[mod(bound_left,n_g)] += charge * charge_fraction * avg_velocity;
    // Pieces connecting two gridpoints
    charge_fraction = (1.0) / length;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      distance_from_left = cell + 0.5 - left;
      avg_velocity = interpolate_segment_velocity(v_y_left, v_y_right, length, distance_from_left);
      j_y[mod(cell,n_g)] += charge * charge_fraction * avg_velocity;
    }
    // Right end
    charge_fraction = (right - (bound_right - 1)) / length;
    midpoint = (right + bound_right) / 2.0;
    distance_from_left = midpoint - left;
    avg_velocity = interpolate_segment_velocity(v_y_left, v_y_right, length, distance_from_left);
    j_y[mod((bound_right-1),n_g)] += charge * charge_fraction * avg_velocity;
  }
  return;
}

double j_y_grad_portion(double rho_ave, double v_ave, double rho_grad,
			    double v_grad, double x_ave, double x_a,
			    double x_b)
{
  return rho_ave * v_ave * (x_b - x_a)
    + 0.5 * (pow((x_b-x_ave), 2) - pow((x_a-x_ave), 2)) * (rho_ave * v_grad + rho_grad * v_ave)
    + (1.0 / 3.0) *(pow((x_b-x_ave), 3) - pow((x_a-x_ave), 3)) * rho_grad * v_grad;
}

void deposit_j_y_segment_gradient(std::vector<double> &j_y,
				  double left, double right,
				  double v_y_left, double v_y_right, 
				  double charge, int n_g, double dx,
				  double gradient)
{
  double length, x_ave, v_ave, rho_ave, v_grad, rho_grad;
  int bound_left, bound_right;
  left = left / dx;
  right = right / dx;
  length = right - left;

  bound_left = floor(left);
  bound_right = ceil(right);
  x_ave = (right + left) / 2.0;
  v_ave = (v_y_left + v_y_right) / 2.0;
  rho_ave = (charge / length) * dx;
  v_grad = ((v_y_right - v_y_left) / length) * dx;
  rho_grad = gradient * dx;

  // If tracers are between two gridpoints
  if (bound_right == (bound_left + 1)) {
    if ((left == bound_left + 0.5) && (right == bound_left + 0.5)) {
      j_y[mod((bound_right),n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, left, right);
    } else if (left >= (bound_left + 0.5)) {
      j_y[mod((bound_right),n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, left, right);
    } else if (right <= (bound_left + 0.5)) {
      j_y[mod((bound_left),n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, left, right);
    } else {
      j_y[mod((bound_left),n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, left, bound_left+0.5);
      j_y[mod((bound_right),n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, bound_left+0.5, right);
    }
  }
  else {
    // Left end
    if (left < bound_left + 0.5) {
      j_y[mod(bound_left,n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, left, bound_left+0.5);
      j_y[mod((bound_left+1),n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, bound_left+0.5, bound_left+1.0);
    } else { 
      j_y[mod(bound_left+1,n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, left, bound_left+1.0);
    }

    // Pieces connecting two gridpoints
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      j_y[mod(cell,n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, cell, cell+0.5);
      j_y[mod((cell+1),n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, cell+0.5, cell+1.0);
    }
    
    // Right end
    if (right > bound_right - 0.5) {
      j_y[mod(bound_right-1,n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, bound_right-1, bound_right-0.5);
      j_y[mod(bound_right,n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, bound_right-0.5, right);
    } else { 
      j_y[mod(bound_right-1,n_g)] += j_y_grad_portion(rho_ave, v_ave, rho_grad, v_grad, x_ave, bound_right-1.0, right);
    }
  }

  return;
}

void ParticleSpecies::deposit_j_y_segments_gradient(std::vector<double> &j_y, MPI_Comm COMM)
{
  double x_i_tavg, x_ip1_tavg, left, right, v_y_left, v_y_right;

  MPI_Status status;
  int num_procs, my_rank, dest, source;
  int tag = 0;
  MPI_Comm_size(COMM, &num_procs);
  MPI_Comm_rank(COMM, &my_rank);

  std::vector<double> x_tavg, segment_density, grad_tavg;
  x_tavg.resize(n_p+1);
  segment_density.resize(n_p+1);
  grad_tavg.resize(n_p+1);

  for (int i = 0; i < (n_p+1); i++) {
    x_tavg[i] = (x[i]+x_old[i]) / 2.0;
  }

  for (int i = 0; i < (n_p+1); i++) {
    segment_density[i] = charge[i] / fabs(x_tavg[i+1]-x_tavg[i]);
  }

  for (int i = 1; i < n_p; i++) {
    grad_tavg[i] = 0.5 * (segment_density[i+1] - segment_density[i-1]) / (x_tavg[i+1]-x_tavg[i]);
  }

  dest = mod(my_rank+1, num_procs);
  source = mod(my_rank-1, num_procs);

  MPI_Sendrecv(&grad_tavg[n_p-1], 1, MPI_DOUBLE, dest, tag,
	       &grad_tavg[0], 1, MPI_DOUBLE,
	       source, tag, COMM, &status);
  
  for (int i = 0; i < n_p; i++) {
    x_i_tavg = x_tavg[i];
    x_ip1_tavg = x_tavg[i+1];
    if (x_ip1_tavg >= x_i_tavg) {
      left = x_i_tavg;
      right = x_ip1_tavg;
      v_y_left = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
      v_y_right = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0));
    }
    else {
      left = x_ip1_tavg;
      right = x_i_tavg;
      v_y_left = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0));
      v_y_right = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
    }
    deposit_j_y_segment_gradient(j_y, left, right, v_y_left, v_y_right, charge[i], n_g, dx, grad_tavg[i]);
  }
  return;
}

void ParticleSpecies::deposit_j_y_segments_zero(std::vector<double> &j_y)
{
  double x_i_tavg, x_ip1_tavg, left, right, v_y_left, v_y_right;
  for (int i = 0; i < n_p; i++) {
    x_i_tavg = (x_old[i] + x[i]) / 2.0;
    x_ip1_tavg = (x_old[i+1] + x[i+1]) / 2.0;
    if (x_ip1_tavg >= x_i_tavg) {
      left = x_i_tavg;
      right = x_ip1_tavg;
      v_y_left = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
      v_y_right = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0));
    }
    else {
      left = x_ip1_tavg;
      right = x_i_tavg;
      v_y_left = u_y[i+1] / sqrt(1.0 + pow(u_x[i+1], 2.0) + pow(u_y[i+1], 2.0));
      v_y_right = u_y[i] / sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
    }
    deposit_j_y_segment_zero(j_y, left, right, v_y_left, v_y_right, charge[i], n_g, dx);
  }
  return;
}


void deposit_rho_segment_zero(std::vector<double> &rho, double x_tracer_a, 
			      double x_tracer_b, double charge, int n_g, double dx)
{
  double left, right, length, charge_fraction;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  length = right - left;

  // Shift tracer positions so that 0 is the left bounary of cell zero

  bound_left = floor(left);
  bound_right = ceil(right);

  if (bound_right == (bound_left + 1)) {
    // If tracers are between two gridpoints
    charge_fraction = 1.0;
    rho[mod(bound_left,n_g)] += charge * charge_fraction;
  }
  else {

    // Left end
    charge_fraction = (bound_left + 1.0 - left) / length;
    rho[mod(bound_left,n_g)] += charge * charge_fraction;

    // Pieces connecting two gridpoints
    charge_fraction = (1.0) / length;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      rho[mod(cell,n_g)] += charge * charge_fraction;
    }
    
    // Right end
    charge_fraction = (right - (bound_right - 1)) / length;
    rho[mod((bound_right-1),n_g)] += charge * charge_fraction;
  }
  return;
}


void deposit_charge_to_left_segment_linear(std::vector<double> &j_x, double x_tracer_a, 
					   double x_tracer_b, double charge, int n_g, double dx, int right_max)
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


void ParticleSpecies::deposit_j_x_segments_linear(std::vector<double> &j_x)
{
  double right_max;
  for (int i = 0; i < n_p; i++) {
    // Index of right bounding gridpoint
    right_max = fmax(fmax(x_old[i], x_old[i+1]),fmax(x[i],x[i+1]));
    right_max = right_max / dx;
    right_max = ceil(right_max);
    
    // Initial line segment
    deposit_charge_to_left_segment_linear(j_x, x_old[i], x_old[i+1], (charge[i] * (dx / dt)), n_g, dx, right_max);
    deposit_charge_to_left_segment_linear(j_x, x[i], x[i+1], (-1.0) * (charge[i] * (dx / dt)), n_g, dx, right_max);
  }
  return;
}

void ParticleSpecies::deposit_rho_segments_zero(std::vector<double> &rho)
{
  for (int i = 0; i < n_p; i++) {    
    //deposit_rho_segment_zero(rho, x[i], x[i+1], charge[i], n_g, dx);
    deposit_rho_segment_zero_new(rho, x[i], x[i+1], charge[i], n_g, dx);
  }
  return;
}

void ParticleSpecies::deposit_rho_segments_gradient(std::vector<double> &rho,
						    MPI_Comm COMM)
{
  MPI_Status status;
  int num_procs, my_rank, dest, source;
  int tag = 0;
  MPI_Comm_size(COMM, &num_procs);
  MPI_Comm_rank(COMM, &my_rank);

  std::vector<double> segment_density;
  segment_density.resize(n_p+1);

  for (int i = 0; i < (n_p+1); i++) {
    segment_density[i] = charge[i] / fabs(x[i+1]-x[i]);
  }

  for (int i = 1; i < n_p; i++) {
    gradient_old[i] = 0.5 * (segment_density[i+1] - segment_density[i-1]) / (x[i+1]-x[i]);
  }

  dest = mod(my_rank+1, num_procs);
  source = mod(my_rank-1, num_procs);

  MPI_Sendrecv(&gradient_old[n_p-1], 1, MPI_DOUBLE, dest, tag,
	       &gradient_old[0], 1, MPI_DOUBLE,
	       source, tag, COMM, &status);
  
  for (int i = 0; i < n_p; i++) {
    deposit_rho_single_gradient(rho, x[i], x[i+1], charge[i], n_g, dx,
				gradient_old[i]);
  }
  return;
}

void ParticleSpecies::deposit_rho_segments_linear(std::vector<double> &rho)
{
  for (int i = 0; i < n_p; i++) {    
    deposit_rho_segment_linear(rho, x[i], x[i+1], charge[i], n_g, dx);
  }
  return;
}


void ParticleSpecies::deposit_j_x_segments_zero(std::vector<double> &j_x)
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

void ParticleSpecies::deposit_j_x_segments_gradient(std::vector<double> &j_x, MPI_Comm COMM)
{
  double left_initial, right_initial, left_final, right_final, length_initial, length_final,
    charge_initial, charge_final, cell_boundary, midpoint_initial, midpoint_final;
  int bound_left, bound_right;

  MPI_Status status;
  int num_procs, my_rank, dest, source;
  int tag = 0;
  MPI_Comm_size(COMM, &num_procs);
  MPI_Comm_rank(COMM, &my_rank);

  std::vector<double> segment_density;
  segment_density.resize(n_p+1);

  for (int i = 0; i < (n_p+1); i++) {
    segment_density[i] = charge[i] / fabs(x[i+1]-x[i]);
  }

  for (int i = 1; i < n_p; i++) {
    gradient[i] = 0.5 * (segment_density[i+1] - segment_density[i-1]) / (x[i+1]-x[i]);
  }

  dest = mod(my_rank+1, num_procs);
  source = mod(my_rank-1, num_procs);

  MPI_Sendrecv(&gradient[n_p-1], 1, MPI_DOUBLE, dest, tag,
	       &gradient[0], 1, MPI_DOUBLE,
	       source, tag, COMM, &status);
  
  for (int i = 0; i < n_p; i++) {
    // Initial line segment
    left_initial = fmin(x_old[i], x_old[i+1]) / dx;
    right_initial = fmax(x_old[i], x_old[i+1]) / dx;
    length_initial = right_initial - left_initial;
    midpoint_initial = (left_initial + right_initial) / 2.0;
    // Final line segment
    left_final = fmin(x[i], x[i+1]) / dx;
    right_final = fmax(x[i], x[i+1]) / dx;
    length_final = right_final - left_final;
    midpoint_final = (left_final + right_final) / 2.0;
    // Bounding box
    bound_left = floor(fmin(left_initial,left_final));
    bound_right = ceil(fmax(right_initial,right_final));
    // Loop over cell boundaries that could be affected
    for (int cell = bound_left; cell < bound_right; cell++) {
      // Calculate charge initially to the left of the current cell's right boundary
      cell_boundary = cell + 0.5;
      if (right_initial <= cell_boundary) {
	charge_initial = charge[i];
      } 
      else if (left_initial >= cell_boundary) {
	charge_initial = 0.0;
      }
      else {
	charge_initial = charge[i] * (cell_boundary - left_initial) / length_initial
	  + linear_correction(gradient_old[i], midpoint_initial, left_initial, cell_boundary)*dx*dx;
      }
      // Calculate charge finally to the left of the current cell's right boundary
      if (right_final <= cell_boundary) {
	charge_final = charge[i];
      } 
      else if (left_final >= cell_boundary) {
	charge_final = 0.0;
      }
      else {
	charge_final = charge[i] * (cell_boundary - left_final) / length_final
	  + linear_correction(gradient[i], midpoint_final, left_final, cell_boundary)*dx*dx;
      }
      j_x[mod(cell,n_g)] += (charge_initial - charge_final) * (dx / dt);
    }
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


void ParticleSpecies::write_phase(std::ofstream &x_ofstream, 
				  std::ofstream &u_x_ofstream, std::ofstream &u_y_ofstream)
{
  std::vector<double> x_in_box(n_p);
  for (int i = 0; i < n_p; i++) {
    x_in_box[i] = put_in_box(x[i], n_g*dx);
  }

  write_data(x_in_box, x_ofstream, n_p);
  write_data(u_x, u_x_ofstream, n_p);
  write_data(u_y, u_y_ofstream, n_p);
  return;
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

void ParticleSpecies::deposit_rho(std::vector<double> &rho)
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

void ParticleSpecies::deposit_rho_ngp(std::vector<double> &rho)
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

void ParticleSpecies::deposit_j_y_ngp(std::vector<double> &j_y)
{
  double x_tavg, gamma;
  int nearest_gridpoint;
  for (int i = 0; i < n_p; i++) {
    x_tavg = (x[i] + x_old[i]) / 2.0;
    x_tavg = x_tavg / dx;
    nearest_gridpoint = get_nearest_gridpoint(x_tavg);
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
    j_y[mod(nearest_gridpoint,n_g)] += charge[i] * u_y[i] / gamma;
  }
  return;
}

void ParticleSpecies::deposit_j_x_ngp(std::vector<double> &j_x)
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

void ParticleSpecies::deposit_j_x(std::vector<double> &j_x)
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

void ParticleSpecies::deposit_j_y(std::vector<double> &j_y)
{
  int i_lower;
  double alpha, beta, x_tavg, gamma;
  for (int i = 0; i < n_p; i++) {
    x_tavg = (x[i] + x_old[i]) / 2.0;
    x_tavg = x_tavg / dx;
    i_lower = floor(x_tavg);
    alpha = 1.0 - (x_tavg - i_lower);
    beta = 1.0 - alpha;

    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
    j_y[mod(i_lower,n_g)] += alpha * charge[i] * u_y[i] / gamma;
    j_y[mod((i_lower+1),n_g)] += beta * charge[i] * u_y[i] / gamma;
  }
  return;
}

void ParticleSpecies::advance_velocity(std::vector<double> &e_x, 
				       std::vector<double> &e_y, 
				       std::vector<double> &b_z)
{
  // Shift e_x and b_z integer values to eliminate self forces
  std::vector<double> e_x_int(n_g);
  std::vector<double> b_z_int(n_g);
  if (center_fields) {
    half_int_to_int(e_x, e_x_int, n_g);
    half_int_to_int(b_z, b_z_int, n_g);
  }

  double e_x_particle, e_y_particle, b_z_particle, gamma_centered, t, s;
  double energy = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  int ngp_integer, ngp_half_integer;
  for (int i = 0; i < n_p; i++) {
    // Weight fields to the particle location
    if (interp_order==0) {
      if (center_fields) {
	ngp_integer = get_nearest_gridpoint(x[i] / dx);
	e_x_particle = e_x_int[mod(ngp_integer, n_g)];
	e_y_particle = e_y[mod(ngp_integer, n_g)];
	b_z_particle = b_z_int[mod(ngp_integer, n_g)];
      } else {
	ngp_integer = get_nearest_gridpoint(x[i] / dx);
	ngp_half_integer = get_nearest_gridpoint((x[i]-0.5) / dx);
	e_x_particle = e_x[mod(ngp_half_integer, n_g)];
	e_y_particle = e_y[mod(ngp_integer, n_g)];
	b_z_particle = b_z[mod(ngp_half_integer, n_g)];
      }
    } else {
      if (center_fields) {
	e_x_particle = interpolate_field_integer(e_x_int, x[i], dx, n_g);
	e_y_particle = interpolate_field_integer(e_y, x[i], dx, n_g);
	b_z_particle = interpolate_field_integer(b_z_int, x[i], dx, n_g);
      } else {
	e_x_particle = interpolate_field_half_integer(e_x, x[i], dx, n_g);
	e_y_particle = interpolate_field_integer(e_y, x[i], dx, n_g);
	b_z_particle = interpolate_field_half_integer(b_z, x[i], dx, n_g);
      }
    }

    // Relativistic Boris push
    // First half electric impulse
    //Need to gefine tem and gamma
    u_x[i] = u_x[i] + rqm * e_x_particle * (dt / 2.0);
    u_y[i] = u_y[i] + rqm * e_y_particle * (dt / 2.0);
    u_z[i] = u_z[i] + rqm * e_z_particle * (dt / 2.0);    

    tem = 0.5 * dt / rqm;
    gamma =  sqrt(1.0 + pow(u_x[i], 2) + pow(u_y[i], 2) + pow(u_z[i], 2));

    b_x_particle = b_x_particle * tem / gamma;
    b_y_particle = b_y_particle[1] * tem / gamma;
    b_z_particle = b_z_particle[2] * tem / gamma;
    
    u_temp[0] = u[0] + u[1] * b_p[2] - u[2] * b_p[1];
    u_temp[1] = u[1] + u[2] * b_p[0] - u[0] * b_p[2];
    u_temp[2] = u[2] + u[0] * b_p[1] - u[1] * b_p[0];
    
    otsq = 2.0 / (1.0 + pow(b_p[0], 2) + pow(b_p[1], 2) + pow(b_p[2], 2));
    
    b_p[0] = b_p[0] * otsq;
    b_p[1] = b_p[1] * otsq;
    b_p[2] = b_p[2] * otsq;
    
    u[0] = u[0] + u_temp[1] * b_p[2] - u_temp[2] * b_p[1];
    u[1] = u[1] + u_temp[2] * b_p[0] - u_temp[0] * b_p[2];
    u[2] = u[2] + u_temp[0] * b_p[1] - u_temp[1] * b_p[0];

    // Check the energy calculation
    //    gamma_centered = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
    //    energy = energy + (charge[i] / rqm) * (pow(u_x[i], 2.0) + pow(u_y[i], 2.0)) / (1.0 + gamma_centered);
    
    // Second half electric impulse
    u_x[i] = u_x[i] + rqm * e_x_particle * (dt / 2.0);
    u_y[i] = u_y[i] + rqm * e_y_particle * (dt / 2.0);
    u_z[i] = u_z[i] + rqm * e_z_particle * (dt / 2.0);    

    momentum_x += fabs(charge[i]) * u_x[i];
    momentum_y += fabs(charge[i]) * u_y[i];
    momentum_z += fabs(charge[i]) * u_z[i];
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
						    std::vector<double> &b_z)
{
  int ngp_integer, ngp_half_integer;
  // Shift e_x and b_z integer values to eliminate self forces
  std::vector<double> e_x_int(n_g);
  std::vector<double> b_z_int(n_g);
  if (center_fields) {
    half_int_to_int(e_x, e_x_int, n_g);
    half_int_to_int(b_z, b_z_int, n_g);
  }
  
  double e_x_particle, e_y_particle, b_z_particle, gamma_centered, t, s;
  for (int i = 0; i < n_p; i++) {
    // Weight fields to the particle location
    if (interp_order==0) {
      if (center_fields) {
	ngp_integer = get_nearest_gridpoint(x[i] / dx);
	e_x_particle = e_x_int[mod(ngp_integer, n_g)];
	e_y_particle = e_y[mod(ngp_integer, n_g)];
	b_z_particle = b_z_int[mod(ngp_integer, n_g)];
      } else {
	ngp_integer = get_nearest_gridpoint(x[i] / dx);
	ngp_half_integer = get_nearest_gridpoint((x[i]-0.5) / dx);
	e_x_particle = e_x[mod(ngp_half_integer, n_g)];
	e_y_particle = e_y[mod(ngp_integer, n_g)];
	b_z_particle = b_z[mod(ngp_half_integer, n_g)];
      }
    } else {
      if (center_fields) {
	e_x_particle = interpolate_field_integer(e_x_int, x[i], dx, n_g);
	e_y_particle = interpolate_field_integer(e_y, x[i], dx, n_g);
	b_z_particle = interpolate_field_integer(b_z_int, x[i], dx, n_g);
      } else {
	e_x_particle = interpolate_field_half_integer(e_x, x[i], dx, n_g);
	e_y_particle = interpolate_field_integer(e_y, x[i], dx, n_g);
	b_z_particle = interpolate_field_half_integer(b_z, x[i], dx, n_g);
      }
    }
    
    // Update velocities
    gamma_centered = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
    
    // Velocity rotation from magnetic field
    t = rqm * b_z_particle * (-1.0 * dt / 2.0) / (2.0 * gamma_centered);
    s = (2.0 * t) / (1.0 + pow(t, 2.0));
    
    u_x[i] = u_x[i] + u_y[i] * t;
    u_y[i] = u_y[i] - u_x[i] * s;
    u_x[i] = u_x[i] + u_y[i] * t;
    
    // Half electric impulse
    u_x[i] = u_x[i] + rqm * e_x_particle * (-1.0 * dt / 2.0);
    u_y[i] = u_y[i] + rqm * e_y_particle * (-1.0 * dt / 2.0);
  }
  return;
}

void ParticleSpecies::advance_x()
{
  double gamma;
  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
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
  
  while (i < n_p) {
    length = fabs(x[i+1] - x[i]);
    if (length > refinement_length) {
      split_segment_lagrange_3(i);
      //split_segment_linear(i);
    } else {
      i+=1;
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

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

void ParticleSpecies::initialize_species(int species_number, 
					 double n_ppc, 
					 double u_x_drift, 
					 double u_y_drift, 
					 int mode, 
					 double u_x_1, 
					 double u_y_1,
					 int my_rank,
					 int num_procs)
{
  int i_start, i_end, n_ppp;
  n_ppp = (n_g * n_ppc) / num_procs;
  i_start = n_ppp * my_rank;
  i_end = i_start + n_ppp;
  
  std::stringstream ss;
  ss << "particles_";
  ss << species_number;
  species_name = ss.str();

  double k = 2.0 * PI * double(mode) / (n_g * dx);
  double particle_spacing = double(n_g) * dx / double(n_p);

  // Add ghost tracer particle if using line segments
  if ((method==2)||(method==3)) {
    n_p += 1;
    i_end += 1;

    charge.push_back(0.0);
    rqm.push_back(0.0);
    u_x.push_back(0.0);
    u_y.push_back(0.0);
    x.push_back(0.0);
    x_old.push_back(0);
  }

  for (int i = i_start; i < i_end; i++) {
    charge[i-i_start] = (-1.0) * (1.0 / n_ppc);
    rqm[i-i_start] = -1.0;
    x[i-i_start] = double(i) * particle_spacing + particle_spacing / 2.0;
    u_x[i-i_start] = u_x_drift + u_x_1 * sin(k * x[i]);
    u_y[i-i_start] = u_y_drift + u_y_1 * sin(k * x[i]);
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
  for (int i = 0; i < (n_p-1); i++) {
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

void ParticleSpecies::write_energy_history(int n_t, int my_rank, MPI_Comm COMM)
{
  sum_array_to_root(&energy_history[0], n_t, COMM, my_rank);
  sum_array_to_root(&momentum_x_history[0], n_t, COMM, my_rank);
  sum_array_to_root(&momentum_y_history[0], n_t, COMM, my_rank);

  if (my_rank==0) {
    data_to_file(energy_history, (species_name+"_ene"));
    data_to_file(momentum_x_history, (species_name+"_p_x"));
    data_to_file(momentum_y_history, (species_name+"_p_y"));
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

void ParticleSpecies::deposit_j_y_segments_zero(std::vector<double> &j_y)
{
  double x_i_tavg, x_ip1_tavg, left, right, v_y_left, v_y_right;
  for (int i = 0; i < (n_p-1); i++) {
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
  for (int i = 0; i < (n_p-1); i++) {
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
  for (int i = 0; i < (n_p-1); i++) {    
    //deposit_rho_segment_zero(rho, x[i], x[i+1], charge[i], n_g, dx);
    deposit_rho_segment_zero_new(rho, x[i], x[i+1], charge[i], n_g, dx);
  }
  return;
}

void ParticleSpecies::deposit_rho_segments_linear(std::vector<double> &rho)
{
  for (int i = 0; i < (n_p-1); i++) {    
    deposit_rho_segment_linear(rho, x[i], x[i+1], charge[i], n_g, dx);
  }
  return;
}


void ParticleSpecies::deposit_j_x_segments_zero(std::vector<double> &j_x)
{
  double left_initial, right_initial, left_final, right_final, length_initial, length_final,
    charge_initial, charge_final, cell_boundary;
  int bound_left, bound_right;
  for (int i = 0; i < (n_p-1); i++) {
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
    u_x[i] = u_x[i] + rqm[i] * e_x_particle * (dt / 2.0);
    u_y[i] = u_y[i] + rqm[i] * e_y_particle * (dt / 2.0);
    
    gamma_centered = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
    energy = energy + (charge[i] / rqm[i]) * (pow(u_x[i], 2.0) + pow(u_y[i], 2.0)) / (1.0 + gamma_centered);
    
    // Velocity rotation from magnetic field
    t = rqm[i] * b_z_particle * dt / (2.0 * gamma_centered);
    s = (2.0 * t) / (1.0 + pow(t, 2.0));
    
    u_x[i] = u_x[i] + u_y[i] * t;
    u_y[i] = u_y[i] - u_x[i] * s;
    u_x[i] = u_x[i] + u_y[i] * t;
    
    // Second half electric impulse
    u_x[i] = u_x[i] + rqm[i] * e_x_particle * (dt / 2.0);
    u_y[i] = u_y[i] + rqm[i] * e_y_particle * (dt / 2.0);

    momentum_x += fabs(charge[i]) * u_x[i];
    momentum_y += fabs(charge[i]) * u_y[i];
  }
  energy_history.push_back(energy);
  momentum_x_history.push_back(momentum_x);
  momentum_y_history.push_back(momentum_y);
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
    t = rqm[i] * b_z_particle * (-1.0 * dt / 2.0) / (2.0 * gamma_centered);
    s = (2.0 * t) / (1.0 + pow(t, 2.0));
    
    u_x[i] = u_x[i] + u_y[i] * t;
    u_y[i] = u_y[i] - u_x[i] * s;
    u_x[i] = u_x[i] + u_y[i] * t;
    
    // Half electric impulse
    u_x[i] = u_x[i] + rqm[i] * e_x_particle * (-1.0 * dt / 2.0);
    u_y[i] = u_y[i] + rqm[i] * e_y_particle * (-1.0 * dt / 2.0);
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

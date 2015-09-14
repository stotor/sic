#include <vector>
#include <cmath>
#include "particles.hpp"
#include "utilities.hpp"

////////////////////////////////////////////////////////////////
// Particles

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

void deposit_rho_segment_zero(std::vector<double> &rho, double x_tracer_a, 
			      double x_tracer_b, double charge, int n_g, double dx)
{
  double left, right, length, charge_fraction;
  int bound_left, bound_right;
  left = fmin(x_tracer_a, x_tracer_b) / dx;
  right = fmax(x_tracer_a, x_tracer_b) / dx;
  length = right - left;

  // Shift tracer positions so that 0 is the left bounary of cell zero
  left = left + 0.5;
  right = right + 0.5;

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
    deposit_rho_segment_zero(rho, x[i], x[i+1], charge[i], n_g, dx);
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

ParticleSpecies::~ParticleSpecies()
{
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

ParticleSpecies::ParticleSpecies(int n_p)
{
  this->n_p = n_p;
  x.resize(n_p);
  u_x.resize(n_p);
  u_y.resize(n_p);
  x_old.resize(n_p);
  u_x_old.resize(n_p);
  u_y_old.resize(n_p);
  charge.resize(n_p);
  rqm.resize(n_p);

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

void boris_push(double &u_x, double &u_y, double rqm, double e_x, double e_y,
		double b_z, double dt)
{  
  // Variables necessary for intermediate steps of the calculation
  double gamma_centered, t, s;

  // First half electric impulse
  u_x = u_x + rqm * e_x * (dt / 2.0);
  u_y = u_y + rqm * e_y * (dt / 2.0);

  gamma_centered = sqrt(1.0 + pow(u_x, 2.0) + pow(u_y, 2.0));

  // Velocity rotation from magnetic field
  t = rqm * b_z * dt / (2 * gamma_centered);
  s = (2.0 * t) / (1.0 + pow(t, 2.0));

  u_x = u_x + u_y * t;
  u_y = u_y - u_x * s;
  u_x = u_x + u_y * t;
  
  // Second half electric impulse
  u_x = u_x + rqm * e_x * (dt / 2.0);
  u_y = u_y + rqm * e_y * (dt / 2.0);

  return;
}

void ParticleSpecies::advance_velocity(std::vector<double> &e_x, 
				       std::vector<double> &e_y, 
				       std::vector<double> &b_z_tavg)
{
  double e_x_particle, e_y_particle, b_z_particle;

  for (int i = 0; i < n_p; i++) {
    // Weight fields to the particle location
    e_x_particle = interpolate_field_half_integer(e_x, x[i], dx, n_g);
    e_y_particle = interpolate_field_integer(e_y, x[i], dx, n_g);
    b_z_particle = interpolate_field_half_integer(b_z_tavg, x[i], dx, n_g);
    // Update velocities
    boris_push(u_x[i], u_y[i], rqm[i], e_x_particle, e_y_particle, 
	       b_z_particle, dt);
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

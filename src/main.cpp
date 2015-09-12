/**
 * 1DEM
 * - 1D electromagnetic code
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define PI (4 * atan(1.0))

// Initialize EM fields
void initialize_fields(std::vector<double> e_x, std::vector<double> e_y, std::vector<double> b_z, std::vector<double> j_x,
		       std::vector<double> j_y, int n_g)
{
  // EM wave
  // int mode = 4;
  // for (int i = 0; i < n_g; i++) {
  //   e_x[i] = 0.0;
  //   e_y[i] = cos(double(mode) * 2.0 * PI * (double(i) / double(n_g)));
  //   b_z[i] = cos(double(mode) * 2.0 * PI * (double(i) / double(n_g)));
  //   // BELOW FOR TESTING EM WAVE
  //   j_x[i] = 0.0;
  //   j_y[i] = 0.0;
  // }
  // Zero
  for (int i = 0; i < n_g; i++) {
    e_x[i] = 0.0;
    e_y[i] = 0.0;
    b_z[i] = 0.0;
    j_x[i] = 0.0;
    j_y[i] = 0.0;
  }
};


void write_data(std::vector<double> &array, std::ofstream &file, int n_values)
{
  for (int i = 0; i < n_values; i++) {
    file << array[i] << " ";
  }
  file << std::endl;
  return;
};

int mod(int a, int b)
{
  return (a%b + b) % b;
};

////////////////////////////////////////////////////////////////
// Particles

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
  void advance_velocity(std::vector<double> &e_x, std::vector<double> &e_y, 
			std::vector<double> &b_z_tavg);
  void deposit_j_x(std::vector<double> &j_x);
  void deposit_j_x_segments_zero(std::vector<double> &j_x);
  void deposit_j_x_segments_linear(std::vector<double> &j_x);
  void deposit_j_y(std::vector<double> &j_y);
  void write_phase(std::ofstream &x_ofstream, 
		   std::ofstream &u_x_ofstream, std::ofstream &u_y_ofstream);

};


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
  }
  else {
    // Left end
    midpoint = (left + bound_left + 1) / 2.0;
    charge_fraction = (bound_left + 1 - left) / length;
    weight = (bound_left + 1 - midpoint);
    for (int cell = bound_left; cell < bound_right; cell++) {
      j_x[mod(cell,n_g)] += charge * charge_fraction * weight;

    }

    // Pieces connecting two gridpoints
    charge_fraction = (1.0) / length;
    weight = 0.5;
    for (int cell = (bound_left+1); cell < (bound_right-1); cell++) {
      for (int boundary = cell; boundary < bound_right; boundary++) {
	j_x[mod(boundary,n_g)] += charge * charge_fraction * weight;
      }
    }
    
    // Right end
    midpoint = (right + bound_right - 1) / 2.0;
    charge_fraction = (right - (bound_right - 1)) / length;
    weight = bound_right - midpoint;
    j_x[mod((bound_right-1),n_g)] += charge * charge_fraction * weight;
  }
  return;
};

void ParticleSpecies::deposit_j_x_segments_linear(std::vector<double> &j_x)
{
  for (int i = 0; i < (n_p-1); i++) {
    // Initial line segment
    deposit_charge_to_left_segment_linear(j_x, x_old[i], x_old[i+1], (charge[i] * (dx / dt)), n_g, dx);
    deposit_charge_to_left_segment_linear(j_x, x[i], x[i+1], (-1.0) * (charge[i] * (dx / dt)), n_g, dx);
  }
  return;
};

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
};

ParticleSpecies::~ParticleSpecies()
{
};

double put_in_box(double x, double box_length)
{
  if (x < 0.0) {
    x = x + box_length;
  } 
  else if (x >= box_length) {
    x = x - box_length;
  }
  return x;
};


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
};

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
};


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
};

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
};

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
};

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
};

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
};

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
};


void ParticleSpecies::advance_x()
{
  double gamma;
  for (int i = 0; i < n_p; i++) {
    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
    x[i] = x[i] + dt * u_x[i] / gamma;
  }
  return;
};

////////////////////////////////////////////////////////////////
// Fields
void save_old_values(std::vector<double> &array, std::vector<double> &array_old,
		     int n_values)
{
  for (int i = 0; i < n_values; i++) {
    array_old[i] = array[i];
  }
  return;
};

void advance_b_z(std::vector<double> &b_z, std::vector<double> &e_y, double dt, 
		 double dx, int n_g)
{
  for (int i = 0; i < (n_g - 1); i++) {
    b_z[i] = b_z[i] - (dt / dx) * (e_y[i+1] - e_y[i]);
  }
  b_z[n_g-1] = b_z[n_g-1] - (dt / dx) * (e_y[0] - e_y[n_g-1]);
  return;
};

void advance_e_x(std::vector<double> &e_x, std::vector<double> &j_x, double dt, 
		 double dx, int n_g)
{
  for (int i = 0; i < n_g; i++) {
    e_x[i] = e_x[i] - dt * j_x[i];
  }
  return;
};

void advance_e_y(std::vector<double> &e_y, std::vector<double> &b_z, 
		 std::vector<double> &j_y, double dt, 
		 double dx, int n_g)
{
  e_y[0] = e_y[0] - (dt / dx) * (b_z[0] - b_z[n_g-1]) - dt * j_y[0];
  for (int i = 1; i < n_g; i++) {
    e_y[i] = e_y[i] - (dt / dx) * (b_z[i] - b_z[i-1]) - dt * j_y[i];
  }
  return;
};

void write_data_tavg(std::vector<double> &array, std::vector<double> &array_old, std::ofstream &file, 
		     int n_values)
{
  for (int i = 0; i < n_values; i++) {
    file << (array[i] + array_old[i]) / 2.0 << " ";
  }
  file << std::endl;
  return;
};

void calc_tavg_arry(std::vector<double> &array_tavg, std::vector<double> &array, 
		    std::vector<double> &array_old, 
		    int n_values)
{
  for (int i = 0; i < n_values; i++) {
    array_tavg[i] = (array[i] + array_old[i]) / 2.0;
  }
  return;
};

int main(int argc, char *argv[])
{
  ///////////////////////////////
  // INITIALIZATION

  // Define simulation parameters
  int n_t = 1000;
  int n_g = 200;

  double dx = 0.5;
  double dt = 0.45;

  int n_species = 1;
  int n_p = 20000;

  bool line_segments = true;

  // Initialize
  std::vector<double> e_x(n_g);
  std::vector<double> e_y(n_g);
  std::vector<double> b_z(n_g); 
  std::vector<double> j_y(n_g); 
  std::vector<double> j_x(n_g); 
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
    j_y_ofstream;

  std::ofstream x_ofstream, u_x_ofstream, u_y_ofstream;

  e_x_ofstream.open("e_x");
  e_y_ofstream.open("e_y");
  b_z_ofstream.open("b_z");
  j_x_ofstream.open("j_x");
  j_y_ofstream.open("j_y");

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
  double wave_amplitude = 0.01;
  int wave_mode = 4;
  for (int i_species = 0; i_species < n_species; i_species++) {
    for (int i_particle = 0; i_particle < species[i_species].n_p; i_particle++) {
      species[i_species].charge[i_particle] = (-1.0) * double(n_g) / n_p;
      species[i_species].rqm[i_particle] = -1.0;
      species[i_species].u_x[i_particle] = wave_amplitude * cos(2.0 * PI * double(wave_mode) * (double(i_particle) / species[i_species].n_p));      
      species[i_species].u_y[i_particle] = 0.0;
      species[i_species].x[i_particle] = (double(i_particle) / species[i_species].n_p) * n_g * dx;      
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

  // Calculate b_z, u_x, and u_y at t = - 1/2
  advance_b_z(b_z, e_y, (-0.5 * dt), dx, n_g);

  ////////////////////////////////////////
  // MAIN LOOP
  int ndump_phase = 0;
  for (int t = 0; t < n_t; t++) {
    // if (t % ndump_phase == 0) {
    //   species[0].write_phase(x_ofstream, u_x_ofstream, u_y_ofstream);
    // }
    // Print diagnostics for e_x, e_y, and x
    write_data(e_x, e_x_ofstream, n_g);
    write_data(e_y, e_y_ofstream, n_g);

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

    // Calculate u_x and u_y, at t+1/2
    for (int i = 0; i < n_species; i++) {
      species[i].advance_velocity(e_x, e_y, b_z_tavg);
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
      //species[i].deposit_j_x_segments_zero(j_x);
      species[i].deposit_j_x_segments_linear(j_x);
      //species[i].deposit_j_x(j_x);
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

  x_ofstream.close();
  u_x_ofstream.close();
  u_y_ofstream.close();


  return 0;

}   
// End of main

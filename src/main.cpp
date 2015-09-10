/**
 * 1DEM
 * - 1D electromagnetic code
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define PI (4 * atan(1.0))

////////////////////////////////////////////////////////////////
// Particles

class ParticleSpecies {
public:
  ParticleSpecies(int n_p);
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
  void deposit_j_y(std::vector<double> &j_y);
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
  i_lower = int(x);
  alpha = 1.0 - (x - i_lower);
  beta = 1.0 - alpha;
  return (alpha * field[i_lower%n_g]) + (beta * field[(i_lower+1)%n_g]);
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
  i_lower = int(x);
  alpha = 1.0 - (x - i_lower);
  beta = 1.0 - alpha;
  return (alpha * field[i_lower%n_g]) + (beta * field[(i_lower+1)%n_g]);
};

void ParticleSpecies::deposit_j_x(std::vector<double> &j_x)
{
  double x_a, x_b;
  for (int i = 0; i < n_p; i++) {
      x_a = x_old[i] / dx;
      x_b = x[i] / dx;
      if (int(x_a)==int(x_b)) {
	j_x[int(x_a)%n_g] += charge[i] * (x_b - x_a) * (dx / dt);
      } 
      else if (int(x_a) < int(x_b)) {
	j_x[int(x_a)%n_g] += charge[i] * (1.0 - (x_a - int(x_a))) * (dx / dt);
	j_x[int(x_b)%n_g] += charge[i] * (x_b - int(x_b)) * (dx / dt);
      }
      else {
	j_x[int(x_a)%n_g] += charge[i] * (int(x_a) - x_a) * (dx / dt);
	j_x[int(x_b)%n_g] += charge[i] * (1.0 - (x_b - int(x_b))) * (dx / dt);
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
    i_lower = int(x_tavg);
    alpha = 1.0 - (x_tavg - i_lower);
    beta = 1.0 - alpha;

    gamma = sqrt(1.0 + pow(u_x[i], 2.0) + pow(u_y[i], 2.0));
    j_y[i_lower%n_g] += alpha * charge[i] * u_y[i] / gamma;
    j_y[(i_lower+1)%n_g] += beta * charge[i] * u_y[i] / gamma;
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

void write_data(std::vector<double> &array, std::ofstream &file, int n_values)
{
  for (int i = 0; i < n_values; i++) {
    file << array[i] << " ";
  }
  file << std::endl;
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
  int n_g = 100;

  double dx = 0.5;
  double dt = 0.5;

  int n_species = 2;
  int n_p = 10000;

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
  }
  
  // Initialize output file streams
  std::ofstream e_x_ofstream, e_y_ofstream, b_z_ofstream, j_x_ofstream, 
    j_y_ofstream;
  e_x_ofstream.open("e_x");
  e_y_ofstream.open("e_y");
  b_z_ofstream.open("b_z");
  j_x_ofstream.open("j_x");
  j_y_ofstream.open("j_y");

  // Initialize EM fields
  // EM wave
  int mode = 4;
  for (int i = 0; i < n_g; i++) {
    e_x[i] = 0.0;
    e_y[i] = cos(double(mode) * 2.0 * PI * (double(i) / double(n_g)));
    b_z[i] = cos(double(mode) * 2.0 * PI * (double(i) / double(n_g)));
    // BELOW FOR TESTING EM WAVE
    j_x[i] = 0.0;
    j_y[i] = 0.0;
  }
  // Zero
  // for (int i = 0; i < n_g; i++) {
  //   e_x[i] = 0.0;
  //   e_y[i] = 0.0;
  //   b_z[i] = 0.0;
  //   j_x[i] = 0.0;
  //   j_y[i] = 0.0;
  // }

  // Initialize species numerical parameters
  for (int i = 0; i < n_species; i++) {
    species[i].dt = dt;
    species[i].dx = dx;
    species[i].n_g = n_g;
    species[i].n_p = n_p;
  }
  std::cout << "Initialized" << std::endl;
  // Initialize charge, x, u_x, and u_y at t = 0
  double u_y_drift[2] = {-0.0, 0.0};
  for (int i_species = 0; i_species < n_species; i_species++) {
    for (int i_particle = 0; i_particle < species[i_species].n_p; i_particle++) {
      species[i_species].charge[i_particle] = 0.0;//double(n_g) / n_p;
      species[i_species].rqm[i_particle] = 0.0; //-1.0;
      species[i_species].u_x[i_particle] = 0.1 * cos(2.0 * PI * double(i_particle) / double(n_p));
      species[i_species].u_y[i_particle] = u_y_drift[i_species];
      species[i_species].x[i_particle] = (double(i_particle) / species[i_species].n_p) * n_g * dx;
      
    }
  }

  // Calculate b_z, u_x, and u_y at t = - 1/2
  advance_b_z(b_z, e_y, (-0.5 * dt), dx, n_g);

  ////////////////////////////////////////
  // MAIN LOOP

  for (int t = 0; t < n_t; t++) {
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
      species[i].deposit_j_x(j_x);
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

  return 0;

}   
// End of main

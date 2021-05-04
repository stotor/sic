/*
           .                                                      ,;           
          ;Wt                        t                   i      f#i            
         f#EEj            ..       : ED.                LE    .E#t             
       .E#f E#,          ,W,     .Et E#K:              L#E   i#W,   :KW,      L
      iWW;  E#t         t##,    ,W#t E##W;            G#W.  L#D.     ,#W:   ,KG
     L##LffiE#t        L###,   j###t E#E##t          D#K. :K#Wfff;    ;#W. jWi 
    tLLG##L E#t      .E#j##,  G#fE#t E#ti##f        E#K.  i##WLLLLt    i#KED.  
      ,W#i  E#t     ;WW; ##,:K#i E#t E#t ;##D.    .E#E.    .E#L         L#W.   
     j#E.   E#t    j#E.  ##f#W,  E#t E#ELLE##K:  .K#E        f#E:     .GKj#K.  
   .D#j     E#t  .D#L    ###K:   E#t E#L;;;;;;, .K#D          ,WW;   iWf  i#K. 
  ,WK,      E#t :K#t     ##D.    E#t E#t       .W#G            .D#; LK:    t#E 
  EG.       E#t ...      #G      ..  E#t      :W##########Wt     tt i       tDj
  ,         ,;.          j                    :,,,,,,,,,,,,,.                  
                L.                      
            t   EW:        ,ft          
            Ej  E##;       t#E          
            E#, E###t      t#E          
            E#t E#fE#f     t#E          
  .......   E#t E#t D#G    t#E .......  
  GEEEEEEf. E#t E#t  f#E.  t#E GEEEEEEf.
            E#t E#t   t#K: t#E          
            E#t E#t    ;#W,t#E          
            E#t E#t     :K#D#E          
            E#t E#t      .E##E          
            E#t ..         G#E          
            ,;.             fE          
                             ,                            
        .,         ,;                              
       ,Wt       f#i            i              i   
      i#D.     .E#t            LE             LE   
     f#f      i#W,            L#E            L#E   
   .D#i      L#D.            G#W.           G#W.   
  :KW,     :K#Wfff;         D#K.           D#K.    
  t#f      i##WLLLLt       E#K.           E#K.     
   ;#G      .E#L         .E#E.          .E#E.      
    :KE.      f#E:      .K#E           .K#E        
     .DW:      ,WW;    .K#D           .K#D         
       L#,      .D#;  .W#G           .W#G          
        jt        tt :W##########Wt :W##########Wt 
                     :,,,,,,,,,,,,,.:,,,,,,,,,,,,,.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>

#include "mpi.h"

#include "speciesgroup.hpp"
#include "fields.hpp"
#include "utilities.hpp"

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  int num_procs, my_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (argc != 3 && argc != 4) {
    std::cout << "Usage: ./sic <method> <n_ppc> [<refinement length>]" << std::endl;
    return 1;
  }

  bool center_fields = true;

  std::stringstream ss;

  std::vector<int> method;
  method.push_back(-1);
  method.push_back(-1);
  ss << argv[1];
  ss >> method[1];
  method[0] = method[1];
  //  method[0] = 1;
  int interp_order = method[0] % 4;
  
  ss.str(std::string());
  ss.clear();
  std::vector<double> n_ppc;
  n_ppc.push_back(-1);
  n_ppc.push_back(-1);  
  ss << argv[2];
  ss >> n_ppc[1];
  n_ppc[0] = n_ppc[1];
  //n_ppc[0] = 131072;

  double refinement_length;
  if (argc == 4) {
    ss.str(std::string());
    ss.clear();
    ss << argv[3];
    ss >> refinement_length;
  } else {
    refinement_length = 0.0;
  }

  int n_species, n_t, n_g;
  double dx, dt, rho_bg;
  bool gravity = false;

  int simulation_type = 1;

  if (simulation_type==0 or simulation_type==1) {
    // Weibel and two-stream parameters
    n_t = 500;
    n_g = 128;
    dx = 0.1;
    dt = 0.09;
    n_species = 2;
    rho_bg = 2.0;
  }
  else if (simulation_type==2) {
    // Whistler heating parameters
    //    n_t = 5000;
    n_t = 2500;
    n_g = 1000/2;
    dx = 0.014111*2;
    dt = 0.01411*2;
    n_species = 2;
    rho_bg = 0.0;
  }
  else if (simulation_type==-1) {
    // Plasma wave parameters
    n_t = 500;
    n_g = 128;
    dx = 0.1;
    dt = 0.09;
    n_species = 1;
    rho_bg = 1.0;
    gravity = false;
  }
  else if (simulation_type==-2) {
    // Two-stream one mode
    n_t = 2000*4;
    n_g = 128*4;
    dx = 0.01/4;
    dt = 0.009/4;
    n_species = 2;
    rho_bg = 2.0;
  }

  SpeciesGroup particles(n_species, dt, dx, n_g, center_fields, interp_order);
  particles.initialize_species(n_ppc, my_rank, num_procs, method, simulation_type);
  
  particles.communicate_ghost_particles(MPI_COMM_WORLD);
  particles.calculate_segment_density(MPI_COMM_WORLD);  

  Field e_x(n_g, "e_x", my_rank);
  Field e_y(n_g, "e_y", my_rank);
  Field e_z(n_g, "e_z", my_rank);
  
  Field b_x(n_g, "b_x", my_rank);
  Field b_y(n_g, "b_y", my_rank);  
  Field b_z(n_g, "b_z", my_rank);
  
  Field j_x(n_g, "j_x", my_rank);
  Field j_y(n_g, "j_y", my_rank);
  Field j_z(n_g, "j_z", my_rank);
  
  Field rho(n_g, "rho", my_rank);
  
  particles.deposit_rho(rho.field, rho_bg, my_rank, MPI_COMM_WORLD);

  if (my_rank==0) {
    initialize_e_x(rho.field, e_x.field, dx, n_g, gravity);
    initialize_transverse_em_fields(e_y.field, e_z.field, b_x.field, b_y.field, b_z.field, n_g, dx, simulation_type);
  }

  MPI_Bcast(&e_x.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&e_y.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&e_z.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&b_x.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&b_y.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
  MPI_Bcast(&b_z.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  particles.initial_velocity_deceleration(e_x.field, e_y.field, e_z.field, b_x.field, b_y.field, b_z.field);

  for (int t = 0; t < n_t; t++) {
    particles.deposit_rho(rho.field, rho_bg, my_rank, MPI_COMM_WORLD);
    //    if (t%10==0) {
    //      particles.write_phase(t, my_rank);
    //    }

    if (my_rank==0) {
      e_x.write_field();
      e_y.write_field();
      e_z.write_field();
      b_x.write_field();
      b_y.write_field();      
      b_z.write_field();
      rho.write_field();
      
      e_x.calculate_energy();
      e_y.calculate_energy();
      e_z.calculate_energy();
      b_x.calculate_energy();
      b_y.calculate_energy();      
      b_z.calculate_energy();
    }

    particles.advance_velocity(e_x.field, e_y.field, e_z.field, b_x.field, b_y.field, b_z.field);

    particles.save_x_old();
    particles.advance_x();

    particles.communicate_ghost_particles(MPI_COMM_WORLD);
    particles.calculate_segment_density(MPI_COMM_WORLD);
    if (refinement_length) {
      particles.refine_segments(refinement_length);
    }

    particles.deposit_j_x(j_x.field, my_rank, MPI_COMM_WORLD);
    particles.deposit_j_y(j_y.field, my_rank, MPI_COMM_WORLD);
    particles.deposit_j_z(j_z.field, my_rank, MPI_COMM_WORLD);

    if (my_rank==0) {
      j_x.write_field();
      j_y.write_field();
      j_z.write_field();
      
      advance_e_x(e_x.field, j_x.field, dt, dx, n_g, gravity);
      
      advance_b_z(b_z.field, e_y.field, dt/2.0, dx, n_g);
      advance_e_y(e_y.field, b_z.field, j_y.field, dt, dx, n_g);
      advance_b_z(b_z.field, e_y.field, dt/2.0, dx, n_g);
      
      advance_b_y(b_y.field, e_z.field, dt/2.0, dx, n_g);
      advance_e_z(e_z.field, b_y.field, j_z.field, dt, dx, n_g);
      advance_b_y(b_y.field, e_z.field, dt/2.0, dx, n_g);
    }

    MPI_Bcast(&e_x.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&e_y.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&e_z.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(&b_x.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b_y.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b_z.field[0], n_g, MPI_DOUBLE, 0, MPI_COMM_WORLD);    

    if (my_rank==0) {
      std::cout << t << std::endl;
    }
  }
  
  particles.write_particle_diagnostics(n_t, my_rank, MPI_COMM_WORLD);
  
  if (my_rank==0) {
    e_x.write_energy_history();
    e_y.write_energy_history();
    e_z.write_energy_history();
    b_x.write_energy_history();
    b_y.write_energy_history();    
    b_z.write_energy_history();
  }

  MPI_Finalize();

  return 0;
}

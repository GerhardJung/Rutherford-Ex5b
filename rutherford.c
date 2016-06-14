/******************************************************************************************************************************************
 * rutherford.c:
 * - calculate rutherford formula using simulations (electrostatic and hard-core) and theory for gold atoms
 * - analyse trajectories
 * - determine radius from simulations
 ******************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>	//sqrt(x)
#include <unistd.h>	//getopt

#define M_PI 3.14159265358979323846


//read in variables
double impact_start;
double impact_step;
double impact_stop;
double energy_start;
double energy_step;
double energy_stop;
char *outFile;
char *trajFile;
int traj_flag;

//global data (constants etc.)
int impact_steps;
int energy_steps;

//atomic units
double fm2au = 0.000018903;
double mev2au = 36749.32;

//constants
int Z1 = 79;		//atomic number of gold
int Z2 = 2;		//atomic number of helium
double c = 137;		//speed of light (in a.u.)
double mp = 1836;	//proton mass (in a.u.)
double m0 = 4.0*mp;	//approx. helium mass
double radius = 12*fm2au;	//approx. gold radius

//simulation
double x0 = -10000*fm2au;
double dt = 0.0000001;


//functions
void calculate_theory (double *result);
void calculate_simulation (double *result,double *closest_approach);
void output (double *theory, double *simulation, double *approach);

//auxiliary functions
double calc_gamma (double energy);
double calc_beta (double energy);
double calc_speed (double energy);
static double Acot (double x);
void simulate (double impact, double energy, double &angle_result, double &closest_approach);
void single_timestep (double &x, double &y, double &vx, double &vy, int &collision);
int continue_simulation (double x, double old_vx, double old_vy, double vx, double vy);

// how to use program
static void show_usage(char* argv[]){
  std::cout << "Usage: " << argv[0] 
	    << " impact_parameter start [fm]"
	    << " impact_parameter step [fm]"
	    << " impact_parameter stop [fm]"
	    << " alpha_energy start  [MeV]"
	    << " alpha_energy step [MeV]"
	    << " alpha_energy stop [MeV]"
	    << " output File"
	    << std::endl;
}

// read in parameters from command line
void read_command(int argc, char* argv[])
{
  if(argc < 8){ 
    show_usage(argv);
    exit(1);
  }
  impact_start = atof(argv[1])*fm2au;
  impact_step = atof(argv[2])*fm2au;
  impact_stop = atof(argv[3])*fm2au;
  energy_start = atof(argv[4])*mev2au;
  energy_step = atof(argv[5])*mev2au;
  energy_stop = atof(argv[6])*mev2au;
  outFile = argv[7];
  
  int opt;
  while ((opt = getopt(argc, argv, "t:")) != -1) {
    switch (opt) {
      case 't': 
	traj_flag = 1;
	trajFile = argv[optind-1];
	break;
    }
  }
					

  std::cout << "Parameters:"  
	    << "\n impact_parameter start [fm] " << impact_start/fm2au
	    << "\n impact_parameter step [fm] " << impact_step/fm2au
	    << "\n impact_parameter stop [fm] " << impact_stop/fm2au
	    << "\n energy start [MeV] " << energy_start/mev2au
	    << "\n energy step [MeV] " << energy_step/mev2au
	    << "\n energy stop [MeV] " << energy_stop/mev2au
	    << "\n output file " << outFile
	    << std::endl;
	    
  if (traj_flag) {
    std::cout << " trajectory file " << trajFile
    << std::endl;
    
  }
}

//main
int main(int argc, char* argv[]){
  printf("\nSTARTING Rutherford-Analyser\n");
  printf("READING Parameter\n");
  read_command(argc, argv);
  
  impact_steps = int((impact_stop - impact_start)/impact_step+0.5) + 1;
  energy_steps = int((energy_stop - energy_start)/energy_step+0.5) + 1;
  
  double *theory_angles = new double[impact_steps*energy_steps];
  double *simulation_angles = new double[impact_steps*energy_steps];
  double *closest_approach = new double[impact_steps*energy_steps];
  
  //initialize traj file
  if (traj_flag) {
    FILE * output;
    output = fopen(trajFile, "w");
    
    fprintf(output,"atom 0 name helium radius 0.001\n");
    fprintf(output,"atom 1 name gold radius 0.005\n");
    
    fclose(output);
  }
  
  //calculate theoretical and simulation results
  printf("CACLCULATING theory\n");
  calculate_theory(theory_angles);
  printf("CACLCULATING simulation\n");
  calculate_simulation(simulation_angles,closest_approach);
  
  //analyse the results
  printf("ANALYSE results\n");
  //TODO
  
  
  //output the results
  printf("PRINT results\n\n");
  output(theory_angles, simulation_angles, closest_approach);
  
  
}

//calculate theory
void calculate_theory(double *result){
  int i,e;
  for (i=0; i<impact_steps; i++){
    for (e=0; e<energy_steps; e++){
      double impact = impact_start+i*impact_step;
      double energy = energy_start+e*energy_step;
      double cotth2 = impact*m0*calc_gamma(energy)*calc_speed(energy)*calc_speed(energy)/Z1/Z2;
      result[i*energy_steps+e]=2.0*Acot(cotth2);
    }
  }
}

//calculate simulation
void calculate_simulation(double *result,double *closest_approach){
  int i,e;
  for (i=0; i<impact_steps; i++){
    for (e=0; e<energy_steps; e++){
      double impact = impact_start+i*impact_step;
      double energy = energy_start+e*energy_step;
      double angle_result = 0;
      double approach_result = 0;
      simulate (impact, energy, angle_result, approach_result);
      result[i*energy_steps+e]=angle_result;
      closest_approach[i*energy_steps+e]=approach_result;
    }
  }
}

//output the results
void output (double *theory, double *simulation, double *approach){
  FILE * output;
  int i,e;
  output = fopen(outFile, "w");
  
  fprintf(output,"Rutherford-Analyser\n\n");
  fprintf(output,"impact [fm] energy [MeV] angle theory [Rad] angle simulation [Rad] closest appr. sim. [fm]\n");
  for (e=0; e<energy_steps; e++){
    for (i=0; i<impact_steps; i++){
      double impact = impact_start+i*impact_step;
      double energy = energy_start+e*energy_step;
      fprintf(output,"%f %f %f %f %f\n",impact/fm2au,energy/mev2au,theory[i*energy_steps+e],simulation[i*energy_steps+e],approach[i*energy_steps+e]/fm2au);
    }
    fprintf(output,"\n");
  }
  
  fclose(output);
}

// performing the actual simulation
void simulate (double impact, double energy, double &angle_result, double &closest_approach){
  //initial conditions of gold atom: (0,0). Helium: (-500 fm, impact).
  double x = x0;
  double y = impact;
  closest_approach = sqrt(x0*x0 + impact*impact);
  double vx = calc_speed(energy); 
  double vy = 0; 
  double old_vx, old_vy;
  
  int collision = 0;
  
  FILE * output;
  if (traj_flag) {
    output = fopen(trajFile, "a");
  }
  
  do {
    // new timestep
    old_vx = vx;
    old_vy = vy;
    // perform simulation
    single_timestep (x, y, vx, vy, collision);
    if (sqrt(x*x + y*y) < closest_approach) closest_approach = sqrt(x*x + y*y);
    if (traj_flag) {
      fprintf(output,"timestep\n");
      fprintf(output,"%f %f 0\n",x,y);
      fprintf(output,"0 0 0\n\n");
    }
  } while (continue_simulation(x, old_vx, old_vy, vx, vy));
  
  if (traj_flag) {
    fclose(output);
  }
  
  // finished simulation
  double costh = vx/sqrt(vx*vx+vy*vy);
  angle_result = acos(costh);
}

// perform one single timestep
void single_timestep (double &x, double &y, double &vx, double &vy, int &collision){
  //printf("before: x=%f, y=%f, vx=%f, vy=%f\n",x,y,vx,vy);
  
  double gamma = 1/sqrt(1-(vx*vx+vy*vy)/(c*c));
  // velocity verlet integrator
  // calculate acceleration
  double r = sqrt(x*x+y*y);
  double f_pre = Z1*Z2/(r*r*r);
  double ax = f_pre*x/(gamma*m0);
  double ay = f_pre*y/(gamma*m0);
  
  x += vx*dt + 0.5*ax*dt*dt;
  y += vy*dt + 0.5*ay*dt*dt;
  
  // test for collison
  if (x*x+y*y<radius*radius && !collision){
    //printf("collison at (%f,%f) with (%f,%f)\n",x,y,vx,vy);
    // calculate scattering angle
    double my_cos = acos(-(x*vx+y*vy)/sqrt(x*x+y*y)/sqrt(vx*vx+vy*vy));
    //printf("angle = %f\n",my_cos);
    // calculate vector to add on old velocity
    double c = sqrt(2*(vx*vx+vy*vy)*(1-cos(2*my_cos)));
    //printf("c = %f\n",c);
    vx = -vx + c/sqrt(x*x+y*y)*y;
    vy = -vy - c/sqrt(x*x+y*y)*x;
    
    collision = 1;
    //printf("after collison at (%f,%f) with (%f,%f)\n",x,y,vx,vy);
    
  } else {
    // recalculate acceleration
    r = sqrt(x*x+y*y);
    f_pre = Z1*Z2/(r*r*r);
    double ax_dt = f_pre*x/(gamma*m0);
    double ay_dt = f_pre*y/(gamma*m0);
  
    vx += (ax+ax_dt)*0.5*dt;
    vy += (ay+ay_dt)*0.5*dt;
  
  }
  
  //printf("after: x=%f, y=%f, vx=%f, vy=%f\n",x,y,vx,vy);
  
}

// test whtether simulation should be continued
int continue_simulation (double x, double old_vx, double old_vy, double vx, double vy){
  // dont stop before approaching the core
  if (x < 0 && x > x0) return 1;
  if (x < x0) return 0;
  
  double ratio = fabs(vx/vy - old_vx/old_vy);
  if (ratio > dt) return 1;
  else return 0;
}

//calculate the reletivistic gamma
double calc_gamma (double energy){
  return (1+energy/(m0*c*c));
}

//calculate the reletivistic beta
double calc_beta (double energy){
  return sqrt(1.0-1.0/(calc_gamma(energy)*calc_gamma(energy)));
}

//calculate the speed of the particle
double calc_speed (double energy){
  return calc_beta(energy)*c;
}

static double Acot(double x)
{
    return x == 0 ? M_PI/2.0 : atan(1/x);
}
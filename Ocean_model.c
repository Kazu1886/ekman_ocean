/*
Project 2: A Simple Ocean Model
This code models a set of equations for a 2-layer model of Ekman heat transport defined by F. Codron 2012: Ekman Heat Transport for Slab Oceans
The model is designed to take input from the Met Office Unified Model (UM) in .dat file form. As a result, this model uses the same spatial 
horizontal resolution as the UM (90 latitude grid points and 144 longitude points, using spacing of 2.0 and 2.5 degrees respectively). The code
solves the heat transport using a finite difference method, and a sparse linear algebra GSL solver.
Model can be ran on an Ubuntu terminal using the following set of commands:
	gcc -Wall -I/GSL_PATH/gsl/include -c pr2_ocean_transport_v7.c
	gcc -L/GSL_PATH/gsl/lib pr2_ocean_transport_v7.o -lgsl -lgslcblas -lm
	./a.out [-v {NO_TRANSPORT, DIFUSION_ONLY, DIFUSION_EKMAN}] 
where the options for -v are integers defined in enum Transport.
Created by Jake K. Eager
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

#define C_D 0.0013 // drag coefficient, dimensionless
#define C_V 4000.0 // specific heat capacity of water J/kg/K
#define EMISSIVITY 0.985 // effective blackbody emissivity of ocean water
#define HEAT_TRANSFER_COEFFICIENT 3000. /* heat transfer coefficient for water W/m2/K 50-3000*/
#define RHO_AIR 1.22 // air density kg/m3
#define RHO_WATER 1027 // km/m3
#define R_PLANET (1.12384*6.371e6) // radius of planet in m (as fraction of Earth radius)
#define EPSILON (0.00001) // 0.00001 pretty good need to find appropriate value
#define OMEGA (6.501e-6) // angular frequency planet, about its central axis
#define SIGMA (5.67e-8) // stefan boltzman constant
#define H_S 50.0 // thickness of surface layer in m
#define H_D 150.0 // thickness of deep layer below
#define D (25000.0) // horizontal diffusion coefficient m2/s
// #define D (25000.0*50.0/(H_S+H_D)) // horizontal diffusion coefficient m2/s MIGHT NEED THIS FIX, UNSURE...
#define N_LATS 90 /* number of latitude points */
#define N_LONS 144 /* number of longitude points */
#define LAT_MIN -89.0
#define DELTA_LAT 2.0
#define LON_MIN 1.25
#define DELTA_LON 2.5
#define HOURS_PER_DAY 24.
#define MINUTES_PER_HOUR 60.
#define SECONDS_PER_MINUTE 60.
#define DELTA_T (1.0*MINUTES_PER_HOUR*SECONDS_PER_MINUTE) // time step of 1 hours
#define TIME_STEPS (15000*2) // for time step of 12 hours, 20,000 time steps needed for 10,000 day run
#define TIME_OUTPUT_FREQ (25.0*HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE) // print update frequency in seconds, but first term is number of days
#define DATA_OUTPUT_FREQ (100.0*HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE) // output frequency in seconds, but first term is number of days
#define TIME_OUTPUT_TOL (1e-6) // days
#define T_OFFSET 0.0 // IF you want to restart run, specify start time here, so it doesn't save over stuff, and change input temp files to output ones
#define TOL 10
#define MAX_ITER (1e-6)
// input data files
#define MAX_FILE_LINE_SIZE 4000
#define SW_FLUX_NET_IN_DATA "input_data/ProCb/surface_net_downward_shortwave_flux.dat"
#define LW_FLUX_IN_DATA "input_data/ProCb/surface_downwelling_longwave_flux_in_air.dat"
#define LATENT_FLUX_DATA "input_data/ProCb/surface_upward_latent_heat_flux.dat"
#define SENSIBLE_FLUX_DATA "input_data/ProCb/surface_upward_sensible_heat_flux.dat"
// #define INITIAL_SURFACE_TEMP_DATA "output_data/ProCb/T_surf_14000_days.dat"
// #define INITIAL_DEEP_TEMP_DATA "output_data/ProCb/T_deep_14000_days.dat"
#define INITIAL_SURFACE_TEMP_DATA "input_data/ProCb/surface_temperature.dat"
#define INITIAL_DEEP_TEMP_DATA "input_data/ProCb/surface_temperature.dat"
#define U_WIND_DATA "input_data/ProCb/x_wind.dat"
#define V_WIND_DATA "input_data/ProCb/y_wind.dat"
#define LATS_FILE "input_data/lats.dat"
#define LONS_FILE "input_data/lons.dat"
// output data files
#define MAX_FNAME_CHAR 100
#define PYTHON_EXE "python"
#define PYTHON_SCRIPT 		"plotting_scripts/zonal_mean_temp_new.py"
#define OUTPUT_UPWARD_LW_SURFACE_FLUX "output_data/ProCb/upward_lw_surface_flux_"
#define OUTPUT_SURFACE_TEMP_DATA "output_data/ProCb/T_surf_"
#define OUTPUT_DEEP_TEMP_DATA "output_data/ProCb/T_deep_"
#define OUTPUT_DATA_EXT "_days.dat"


typedef enum coords { THETA, PHI, N_COORDS } Coords;
typedef enum depths { SURFACE, DEEP, N_DEPTHS } Depths;
typedef enum transport {NO_TRANSPORT, DIFUSION_ONLY, DIFUSION_EKMAN} Transport;
typedef enum errors { MALLOC_ERR, FREE_ERR, FOPEN_ERR, DEPTH_ERR, DERIVATIVE_ERR, READ_DATA_ERR, SOLVER_ERR, VERSION_ERR } Errors;

typedef struct grid_vals { // values of latitude and longitude
	double lat[N_LATS];
	double lon[N_LONS];
} Grid_vals;

static int select_version(int argc, char **argv);
static void time_stepper(int n_times, int n_lats, int n_lons, int version);

int main(int argc, char **argv)
{
	int version = select_version(argc, argv);
	time_stepper(TIME_STEPS,N_LATS,N_LONS,version);

	return 0;
}

/* Selects the whether to solve using no horizontal transport, diffusion only or diffusion and Ekman transport. */
static int select_version(int argc, char **argv) 
{
	extern char *optarg;
	int c, err = 0; 
	// char *case_number=NULL;
	int version=DIFUSION_EKMAN; // default case is the step function test case
	static char usage[] = "usage: %s [-v version_number] [name ...]\n"; 

	while ((c = getopt(argc, argv, "v:")) != -1)
		switch (c) {
		case 'v':
			version=atoi(optarg);
			if(version!=NO_TRANSPORT && version!=DIFUSION_EKMAN && version!=DIFUSION_ONLY){
				fprintf(stderr, "Version selected not valid, please enter [-v {%i, %i, %i}]\n", NO_TRANSPORT, DIFUSION_ONLY, DIFUSION_EKMAN);
				exit(VERSION_ERR);
			}
			break;
		case '?':
			err = 1;
			break;
		}
	if (err) {
		fprintf(stderr, usage, argv[0]);
		exit(1);
	}
	return version;
}

/* Checks memory allocation has been successful */
static void *xmalloc(size_t n)
{ 
    void *p = malloc(n);
    if(p == NULL)   {
        fprintf(stderr,"Out of memory.\n");
        exit(MALLOC_ERR);
    }
    return p;
}

/* checks pointer being freed is not NULL */
static void xfree(void *p)
{
	if(p == NULL)   {
        fprintf(stderr,"Pointer is NULL, so cannot free.\n");
        exit(FREE_ERR);
    }
	free(p);
}

/* allocates a 2d array of pointers */
static double **create_2d_pointer(int n, int m)
{
	double **p = (double **)xmalloc(n*sizeof(double*));
	for(int i=0;i<n;i++){
		p[i]=(double *)xmalloc(m*sizeof(double));
	}
	return p;
}

/* frees a 2d array of pointers */
static void destroy_2d_pointer(double **p, int n)
{
	for(int i=0;i<n;i++){
		xfree(p[i]);
	}
	xfree(p);
}

/* allocates a 3d array of pointers */
static double ***create_3d_pointer(int n, int m, int l)
{
	double ***p = (double ***)xmalloc(n*sizeof(double**));
	for(int i=0;i<n;i++){
		p[i]=(double **)xmalloc(m*sizeof(double*));
		for(int j=0;j<m;j++){
			p[i][j]=(double *)xmalloc(l*sizeof(double));
		}
	}
	return p;
}

/* frees a 3d array of pointers */
static void destroy_3d_pointer(double ***p, int n, int m)
{
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			xfree(p[i][j]);
		}
		xfree(p[i]);
	}
	xfree(p);
}

/* checks file opened correctly */
static void *xfopen(char *filename, char *command)
{ 
	FILE *outfile;
	outfile = fopen(filename,command);
	if(outfile==NULL){
		fprintf(stderr, "%s could not be opened.\n", filename);
		exit(FOPEN_ERR);
	}
	return outfile;
}

static void construct_fname(char* fname, size_t fname_length, const char* base, const char* time, const char* ext)
{ // construct filename from base, time of run and appropriate extention
    strcpy(fname, base);
    strcat(fname,time);
    strcat(fname,ext);
}

/* reads in a .dat files called filename, that contains data with n_lats rows and n_lons collumns, ie UM output data */
static void read_input_data(double **var, char *filename, int n_lats, int n_lons)
{
	char line[MAX_FILE_LINE_SIZE];
	char *data = line;
	FILE *input;
	input = xfopen(filename, "r" );
	int offset;
	for(int i=0;i<n_lats;i++){
		fgets(line, MAX_FILE_LINE_SIZE, input);
		for(int j=0;j<n_lons;j++){
			int found = sscanf(data, " %le%n", &var[i][j], &offset);
			if(found!=1){
				fprintf(stderr, "%s does not contain enough longitudinal points.\n", filename);
				exit(READ_DATA_ERR);
			}
			data+=offset;
		}
		data=line;
	}
  	fclose(input);
}

/* writes out put data to filename on a grid of n_lats by n_lons. */
static void write_data(double **var, char *fbasename, char *time_str, Grid_vals *grid, int n_lats, int n_lons)
{
	char *fname = (char*) xmalloc((MAX_FNAME_CHAR)*sizeof(char));
	construct_fname(fname, MAX_FNAME_CHAR, fbasename, time_str, OUTPUT_DATA_EXT);
	FILE *output;
	output = xfopen(fname, "w" );
	for(int i=0;i<n_lats;i++){
		for(int j=0;j<n_lons;j++){
			fprintf(output, "%lg\t", var[i][j]);
		}
		fprintf(output, "\n");
	}
	fclose(output);
	xfree(fname);
}

static void calc_lw_out_flux(double **F_lw_out, double **T_surf, Grid_vals *grid, int n_lats, int n_lons)
{
	for(int i=0;i<n_lats;i++){
		for(int j=0;j<n_lons;j++){
			F_lw_out[i][j] = EMISSIVITY*SIGMA*pow(T_surf[i][j],4);
		}
	}
}

static void process_output_data(double ***T, Grid_vals *grid, int n_lats, int n_lons, double time)
{ // processes output data
    double days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE);
    days=round(days);
    printf("Data outputted at %lg\n", days);
    int day_int = days+TIME_OUTPUT_TOL;
    char time_str[10];
    sprintf(time_str, "%d", day_int);
    // char *fname = (char*) xmalloc((MAX_FNAME_CHAR)*sizeof(char));
    double **F_lw_out = create_2d_pointer(n_lats,n_lons);
    calc_lw_out_flux(F_lw_out,T[SURFACE],grid,n_lats,n_lons);
    // construct_fname(fname, MAX_FNAME_CHAR, OUTPUT_UPWARD_LW_SURFACE_FLUX, time_str, OUTPUT_DATA_EXT);
    write_data(F_lw_out, OUTPUT_UPWARD_LW_SURFACE_FLUX, time_str, grid, n_lats, n_lons);
    destroy_2d_pointer(F_lw_out,n_lats);
    // construct_fname(fname, MAX_FNAME_CHAR, OUTPUT_SURFACE_TEMP_DATA, time_str, OUTPUT_DATA_EXT);
    write_data(T[SURFACE], OUTPUT_SURFACE_TEMP_DATA, time_str, grid, n_lats, n_lons);
    // construct_fname(fname, MAX_FNAME_CHAR, OUTPUT_DEEP_TEMP_DATA, time_str, OUTPUT_DATA_EXT);
    write_data(T[DEEP], OUTPUT_DEEP_TEMP_DATA, time_str, grid, n_lats, n_lons);
    // xfree(fname);
}

/* convert angle in degrees to radians for evaluating trigonometric functions in calculus */
static double deg_to_rad(double theta)
{
	return theta*M_PI/180.0;
}

/* calculates the first spatial derivative of A with respect to a (which may be theta or phi), on a spherical surface. Used mainly to 
calculate derivative of mass transport components M_theta and M_phi. Boundary conditions are used for a spherical geometry, using A[i][0]
to be adjacent to A[i][n_lons-1], and also A[0][j] is adjacent to the opposite point in longitude (j+-n_lats/2) at the same latitude. */
static double dA_da(double **A, int i, int j, int a, int n_lats, int n_lons)
{
	double d_theta=deg_to_rad(DELTA_LAT);
	double d_phi=deg_to_rad(DELTA_LON);
	if(a==THETA){
		if(i==n_lats-1){
			int j_temp = j+n_lons/2;
			if(j_temp>=n_lons){
				j_temp -=n_lons;
			}
			return (A[i][j_temp]-A[i-1][j])/(2*d_theta);
		}
		else if(i==0){
			int j_temp = j+n_lons/2;
			if(j_temp>=n_lons){
				j_temp -=n_lons;
			}
			return (A[i+1][j]-A[i][j_temp])/(2*d_theta);
		}
		else{
			return (A[i+1][j]-A[i-1][j])/(2*d_theta);
		}
	}
	else if(a==PHI){
		if(j==n_lons-1){
			return (A[i][0]-A[i][j-1])/(2*d_phi);
		}
		else if(j==0){
			return (A[i][j+1]-A[i][n_lons-1])/(2*d_phi);
		}
		else{
			return (A[i][j+1]-A[i][j-1])/(2*d_phi);
		}
	}
	else{
		fprintf(stderr, "Specified derivative is not respect to either theta or phi\n");
		exit(DERIVATIVE_ERR);
	}
}

/* returns height of layer given index h, from enum heights. */
static double height_of_slab(int h)
{
	if(h==SURFACE){
		return H_S;
	}
	else if(h==DEEP){
		return H_D;
	}
	else{
		fprintf(stderr, "Index for depth is out of range.\n");
		exit(DEPTH_ERR);
	}
}

/* Constructs a grid with latitude and longitude points on the sphere, starting at LAT_MIN and LON_MIN, increasing in increments of DELTA_LAT 
and DELTA_LON */
static void construct_grid(Grid_vals *grid, int n_lats, int n_lons, char *lats_file, char *lons_file)
{
	FILE *lats_out;
	lats_out = xfopen(lats_file, "w" );
	for(int i=0;i<n_lats;i++){
		grid->lat[i]=LAT_MIN+(i*DELTA_LAT);
		fprintf(lats_out, "%i\t%lg\n", i, grid->lat[i]);
	}
	fclose(lats_out);

	FILE *lons_out;
	lons_out = xfopen(lons_file, "w" );
	for(int j=0;j<n_lons;j++){
		grid->lon[j]=LON_MIN+(j*DELTA_LON);
		fprintf(lons_out, "%i\t%lg\n", j, grid->lon[j]);
	}
	fclose(lons_out);
}

/* Calculates the surface stress on the ocean surface as a result of the atmospheric winds near the surface. 
Reads in velocity data from .dat files */
static void calculate_surface_stress(double ***tau, Grid_vals *grid, int n_lats, int n_lons)
{
	double ***v = create_3d_pointer(N_COORDS,N_LATS,N_LONS);
	read_input_data(v[THETA],V_WIND_DATA,n_lats,n_lons);
	read_input_data(v[PHI],U_WIND_DATA,n_lats,n_lons);
	for(int coord=THETA;coord<N_COORDS;coord++){
		for(int i=0;i<n_lats;i++){
			for(int j=0;j<n_lons;j++){
				tau[coord][i][j]=C_D*RHO_AIR*v[coord][i][j]*fabs(v[coord][i][j]); // retains sign of velocity
			}
		}
	}
	destroy_3d_pointer(v,N_COORDS,N_LATS);
}

/* Returns the coriolis parameter at an latitude theta */
static double calculate_coriolis_parameter(double theta)
{ // f=2*Omega*sin(theta)
	return 2*OMEGA*sin(deg_to_rad(theta));
}

/* calculates the mass flux due to Ekman transport */
static void calculate_mass_flux(double ***M, Grid_vals *grid, int n_lats, int n_lons)
{
	double ***tau = create_3d_pointer(N_COORDS,N_LATS,N_LONS);

	calculate_surface_stress(tau, grid, n_lats, n_lons);
	for(int i=0;i<n_lats;i++){
		double f = calculate_coriolis_parameter(grid->lat[i]);
		for(int j=0;j<n_lons;j++){
			M[THETA][i][j]=(EPSILON*tau[THETA][i][j]+f*tau[PHI][i][j])/(pow(EPSILON,2)+pow(f,2));
			M[PHI][i][j]=(EPSILON*tau[PHI][i][j]-f*tau[THETA][i][j])/(pow(EPSILON,2)+pow(f,2));
		}
	}
	destroy_3d_pointer(tau,N_COORDS,N_LATS);
}

/* Plots the script called filename in python */
static void plot_data(char *filename)
{ 
	char command[PATH_MAX];
	snprintf(command, sizeof(command), "%s %s", PYTHON_EXE, filename);
	system( command );
}

/* Returns the horizontal divergence of quanitity M. determines whether there is an up or downward flow between the ocean layers */ 
static double calculate_div_M(double ***M, Grid_vals *grid, int i, int j, int n_lats, int n_lons)
{
	double theta=deg_to_rad(grid->lat[i]);
	double dM_theta_dtheta = dA_da(M[THETA],i,j,THETA,n_lats,n_lons);
	double dM_phi_dphi = dA_da(M[PHI],i,j,PHI,n_lats,n_lons);
	double div_M=1./R_PLANET*(dM_theta_dtheta-tan(theta)*M[THETA][i][j]+1./cos(theta)*dM_phi_dphi);
	return div_M;
}

/* Calculates the index for the matrix used in the solver, before being compressed */
static int calculate_matrix_index(int height, int lat, int lon, int n_lats, int n_lons)
{
	return height*n_lats*n_lons+lat*n_lons+lon;
}

/* At maximum or minimum latitutde, the adjacent latitude point above or bellow, respectively, is located at the same latitude, 
but at the new longitudinal point (from lon to new_lon) */
static int calculate_new_lon(int lon, int n_lons)
{
	int new_lon = lon+n_lons/2;
	if(new_lon>=n_lons){
		new_lon -= n_lons;
	}
	return new_lon;
}

/* calculates and sets the matrix elements for a 2 layer ocean with 144 longitude points and 90 latitude points.
version determines whether the matrix is calculated for a no horizontal transport model, diffusion only model or the full Ekman model. */
static void calculate_matrix(gsl_spmatrix *A, Grid_vals *grid, int n_lats, int n_lons, int version)
{
	double ***M=create_3d_pointer(N_COORDS,N_LATS,N_LONS);
	calculate_mass_flux(M, grid, n_lats, n_lons);
	double d_theta=deg_to_rad(DELTA_LAT);
	double d_phi=deg_to_rad(DELTA_LON);
	double s_dfsn = DELTA_T*D/pow(R_PLANET,2);
	
	for(int height=SURFACE; height<N_DEPTHS; height++){
		// if(height == SURFACE){
		// 	s_dfsn = DELTA_T*D*H_S/pow(R_PLANET,2);
		// }
		// else{
		// 	s_dfsn = DELTA_T*D*H_D/pow(R_PLANET,2);
		// }
		double thickness = height_of_slab(height);
		double s_ekmn = DELTA_T/(RHO_WATER*thickness);
		for(int lat=0; lat<n_lats; lat++){
			double theta = deg_to_rad(grid->lat[lat]);
			for(int lon=0; lon<n_lons; lon++){
				double dM_theta_dtheta = dA_da(M[THETA],lat,lon,THETA,n_lats,n_lons);
				double dM_phi_dphi = dA_da(M[PHI],lat,lon,PHI,n_lats,n_lons);
				double M_theta = M[THETA][lat][lon];
				double M_phi = M[PHI][lat][lon];
				double div_M = calculate_div_M(M,grid,lat,lon,n_lats,n_lons);
				double Aij=0.0;
				int i,j;
				i = calculate_matrix_index(height,lat,lon,n_lats,n_lons);
				/* central point */
				j = i;
				if(version == NO_TRANSPORT){
					Aij=1.;
				}
				else{
					Aij=1.+2.*s_dfsn*(1./pow(d_theta,2)+1./pow(cos(theta)*d_phi,2));
					if(version==DIFUSION_EKMAN){
						// Aij+= s_ekmn/R_PLANET*(dM_theta_dtheta-M_theta*tan(theta)+1./cos(theta)*dM_phi_dphi);
						if(height==SURFACE){
							Aij+= s_ekmn/R_PLANET*(dM_theta_dtheta-M_theta*tan(theta)+1./cos(theta)*dM_phi_dphi);
							if(div_M<0.0){
								Aij+=(+div_M);
							}
						}
						else if(height==DEEP){
							Aij+= -s_ekmn/R_PLANET*(dM_theta_dtheta-M_theta*tan(theta)+1./cos(theta)*dM_phi_dphi);
							if(div_M>0.0){
								Aij+=(-div_M);
							}
						}
					}
				}
				gsl_spmatrix_set(A, i, j, Aij);
				/* Neigbouring points */
				if(version != NO_TRANSPORT){
					/* neighbouring latitude - below */
					Aij=-s_dfsn*(1./pow(d_theta,2)+tan(theta)/(2*d_theta));
					if(version==DIFUSION_EKMAN){
						if(height==SURFACE){
							Aij+= s_ekmn/R_PLANET*(-M_theta/(2.*d_theta));
						}
						else{
							Aij+= -s_ekmn/R_PLANET*(-M_theta/(2.*d_theta));	
						}
					}
					if(lat==0){
						int new_lon = calculate_new_lon(lon,n_lons);
						j = calculate_matrix_index(height,lat,new_lon,n_lats,n_lons);
					}
					else{
						j = calculate_matrix_index(height,lat-1,lon,n_lats,n_lons);
					}
					gsl_spmatrix_set(A, i, j, Aij);
					/* neighbouring latitude - above */
					Aij=-s_dfsn*(1./pow(d_theta,2)-tan(theta)/(2*d_theta));
					if(version==DIFUSION_EKMAN){
						if(height==SURFACE){
							Aij+= s_ekmn/R_PLANET*(M_theta/(2.*d_theta));
						}
						else{
							Aij+= -s_ekmn/R_PLANET*(M_theta/(2.*d_theta));
						}
					}
					if(lat==n_lats-1){
						int new_lon = calculate_new_lon(lon,n_lons);
						j = calculate_matrix_index(height,lat,new_lon,n_lats,n_lons);
					}
					else{
						j = calculate_matrix_index(height,lat+1,lon,n_lats,n_lons);
					}
					gsl_spmatrix_set(A, i, j, Aij);
					/* neighbouring longitude - left */
					Aij=-s_dfsn*1./pow(cos(theta)*d_phi,2);
					if(version==DIFUSION_EKMAN){
						if(height==SURFACE){
							Aij+= s_ekmn/R_PLANET*(-M_phi/(2.*d_phi*cos(theta)));
						}
						else{
							Aij+= -s_ekmn/R_PLANET*(-M_phi/(2.*d_phi*cos(theta)));
						}
					}
					if(lon==0){
						int new_lon = n_lons-1;
						j = calculate_matrix_index(height,lat,new_lon,n_lats,n_lons);
					}
					else{
						j = calculate_matrix_index(height,lat,lon-1,n_lats,n_lons);
					}
					gsl_spmatrix_set(A, i, j, Aij);
					/* neighbouring longitude - right */
					Aij=-s_dfsn*1./pow(cos(theta)*d_phi,2);
					if(version==DIFUSION_EKMAN){
						if(height==SURFACE){
							Aij+= s_ekmn/R_PLANET*(M_phi/(2.*d_phi*cos(theta)));
						}
						else{
							Aij+= -s_ekmn/R_PLANET*(M_phi/(2.*d_phi*cos(theta)));
						}
					}
					if(lon==n_lons-1){
						int new_lon = 0;
						j = calculate_matrix_index(height,lat,new_lon,n_lats,n_lons);
					}
					else{
						j = calculate_matrix_index(height,lat,lon+1,n_lats,n_lons);
					}
					gsl_spmatrix_set(A, i, j, Aij);
					/* interaction between deep and surface layer */
					if(version==DIFUSION_EKMAN){
						if(height==SURFACE && div_M>0.0){
							Aij=(+div_M);
							j = calculate_matrix_index(DEEP,lat,lon,n_lats,n_lons);
							gsl_spmatrix_set(A, i, j, Aij);
						}
						else if(height==DEEP && div_M<0.0){
							Aij=(-div_M);
							j = calculate_matrix_index(SURFACE,lat,lon,n_lats,n_lons);
							gsl_spmatrix_set(A, i, j, Aij);
						}
					}
				}
			}
		}
	}
	destroy_3d_pointer(M,N_COORDS,N_LATS);
}

/* Reads in atmospheric fluxes, net stellar sw flux (subrtracts reflected sw radiation from the incident) 
and incoming lw flux from atmospheric emission */
static void calculate_F_a(double **F_a, int n_lats, int n_lons)
{
	double **F_net_sw_in = create_2d_pointer(n_lats,n_lons);
	double **F_lw_in = create_2d_pointer(n_lats,n_lons);
	double **F_latent_in = create_2d_pointer(n_lats,n_lons);
	double **F_sensible_in = create_2d_pointer(n_lats,n_lons);
	read_input_data(F_net_sw_in,SW_FLUX_NET_IN_DATA,n_lats,n_lons);
	read_input_data(F_lw_in,LW_FLUX_IN_DATA,n_lats,n_lons);
	read_input_data(F_latent_in,LATENT_FLUX_DATA,n_lats,n_lons);
	read_input_data(F_sensible_in,SENSIBLE_FLUX_DATA,n_lats,n_lons);
	for(int i=0;i<n_lats;i++){
		for(int j=0;j<n_lons;j++){
			F_a[i][j] = (F_net_sw_in[i][j] + F_lw_in[i][j]) - F_latent_in[i][j] - F_sensible_in[i][j];
		}
	}
	destroy_2d_pointer(F_net_sw_in,n_lats);
	destroy_2d_pointer(F_lw_in,n_lats);
	destroy_2d_pointer(F_latent_in,n_lats);
	destroy_2d_pointer(F_sensible_in,n_lats);
}

/* Checks if convetion condition is met, then calculates the associated flux. */
static void calculate_F_c(double **F_c, double ***T, int n_lats, int n_lons)
{
	for(int i=0;i<n_lats;i++){
		for(int j=0;j<n_lons;j++){
			if(T[DEEP][i][j]>T[SURFACE][i][j]){ /* convection condition */
				F_c[i][j] = HEAT_TRANSFER_COEFFICIENT*(T[DEEP][i][j]-T[SURFACE][i][j]);
			}
			else{
				F_c[i][j] = 0.0;
			}
		}
	}
}

/* calculated and sets the vector b in matrix Ax=b. EMISSIVITY*SIGMA*T^4 is a blackbody emission from one of the layers, 
to ensure an energy exchange between the layers. 1.0/(RHO_WATER*C_V*H_S) prefactor converts a flux in W/m2 to K/s */
static void calculate_vector_b(double ***T, gsl_vector *b, int n_lats, int n_lons)
{ // calculates b from matrix equation Ax=b
	double **F_a = create_2d_pointer(n_lats,n_lons);
	calculate_F_a(F_a,n_lats,n_lons);
	double **F_c = create_2d_pointer(n_lats,n_lons);
	calculate_F_c(F_c,T,n_lats,n_lons);
	double b_hij = 0.0;
	for(int h=0;h<N_DEPTHS;h++){
		for(int i=0;i<n_lats;i++){
			for(int j=0;j<n_lons;j++){
				double depth = height_of_slab(h);
				if(h==SURFACE){
					b_hij = 1.0/(RHO_WATER*C_V*depth)*(F_c[i][j]+F_a[i][j]-EMISSIVITY*SIGMA*pow(T[SURFACE][i][j],4))*DELTA_T + T[h][i][j];
				}
				else{
					b_hij = 1.0/(RHO_WATER*C_V*depth)*(-F_c[i][j])*DELTA_T + T[h][i][j];
				}
				gsl_vector_set(b, (h*n_lats*n_lons+(i*n_lons+j)), b_hij);
			}
		}
	}
	destroy_2d_pointer(F_a,n_lats);
	destroy_2d_pointer(F_c,n_lats);
}

/* Uses GSL sparse linear algebra to solve Ax=b */
static void calculate_new_T(double ***T, Grid_vals *grid, int n_lats, int n_lons, int version)
{
	int n = N_DEPTHS*n_lats*n_lons;
	gsl_spmatrix *A = gsl_spmatrix_alloc(n ,n);
	gsl_spmatrix *C;
	gsl_vector *b = gsl_vector_alloc(n);
	gsl_vector *x = gsl_vector_alloc(n);

	calculate_matrix(A,grid,n_lats,n_lons,version);
	calculate_vector_b(T,b,n_lats,n_lons);
	C = gsl_spmatrix_ccs(A);

	/* now solve the system with the GMRES iterative solver */
    const double tol = TOL;  /* solution relative tolerance */
    const size_t max_iter = MAX_ITER; /* maximum iterations */
    const gsl_splinalg_itersolve_type *S = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(S, n, 0);
    size_t iter = 0;
    double residual;
    int status;

    /* initial guess x = 0 */
    gsl_vector_set_zero(x);

    /* solve the system A x = b */
    do{
		status = gsl_splinalg_itersolve_iterate(C, b, tol, x, work);

		if (status != GSL_SUCCESS && iter==max_iter-1){ /* Solver exits if solution diverges */
			residual = gsl_splinalg_itersolve_normr(work);
			fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);
		  	fprintf(stderr, "Diverged\n");
		  	exit(SOLVER_ERR);
		}
	} 
	while (status == GSL_CONTINUE && ++iter < max_iter);

    /* output solution */
    for(int h=0;h<N_DEPTHS;h++){ /* sets the new temperature from the gsl vector */
	    for(int i=0;i<n_lats;i++){ 
			for(int j=0;j<n_lons;j++){
				T[h][i][j]=gsl_vector_get(x,h*n_lats*n_lons+(i*n_lons+j));
			}
		}
	}

    gsl_splinalg_itersolve_free(work);

	gsl_spmatrix_free(A);
	gsl_spmatrix_free(C);
	gsl_vector_free(b);
	gsl_vector_free(x);
}

/* Function that steps through selected time. First creates the horizontal grid, then read in temperatures, 
and steps through n_times time steps. Writes the final Temperature data to files, and plots a selected plotting script. */
static void time_stepper(int n_times, int n_lats, int n_lons, int version)
{
	Grid_vals *grid = xmalloc(sizeof(Grid_vals));
	construct_grid(grid,n_lats,n_lons,LATS_FILE,LONS_FILE);
	double ***T=create_3d_pointer(N_DEPTHS,N_LATS,N_LONS);
	read_input_data(T[SURFACE],INITIAL_SURFACE_TEMP_DATA,n_lats,n_lons);
	read_input_data(T[DEEP],INITIAL_DEEP_TEMP_DATA,n_lats,n_lons);
	double time = 0.;
	double days = 0.;
	for(int t=0;t<n_times;t++){
		time = ((double)t+1.)*DELTA_T;
		calculate_new_T(T,grid,n_lats,n_lons,version);
		if (fmod(time,TIME_OUTPUT_FREQ)<TIME_OUTPUT_TOL){
			days = T_OFFSET+time/(HOURS_PER_DAY*MINUTES_PER_HOUR*SECONDS_PER_MINUTE);
			printf("Days passed = %lg\n", days);
		}
		if (fmod(time,DATA_OUTPUT_FREQ)<TIME_OUTPUT_TOL){
			process_output_data(T, grid, n_lats, n_lons, time);
		}
	}
	process_output_data(T, grid, n_lats, n_lons, time);
	plot_data(PYTHON_SCRIPT);
	destroy_3d_pointer(T,N_DEPTHS,N_LATS);
	xfree(grid);
}

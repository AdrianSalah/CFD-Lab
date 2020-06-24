#include "helper.hpp"
#include "visual.hpp"
#include "init.hpp"
#include "sor.hpp"
#include "uvp.hpp"
#include <cstdio>
#include <iostream>
#include <sstream>
#include "boundary_val.hpp"
#include "Timer.h"
#include "parameters.h"

//define scenarios using macros
//#define SCENARIO_NAME "lid_driven_cavity"
//#define SCENARIO_DAT_FILE "../parameters/lid_driven_cavity.dat"
//#define SCENARIO_PGM_FILE "../geometry/lid_driven_cavity.pgm"

int main(int argn, char** args) {

    //select scenario
    if (argn == 1)
        scenarioSpec = 1;

    else if (argn == 2) {

        scenarioSpec = atoi(args[1]);
    }

    else {std::cout << "more arguments given" << std::endl;
        exit(EXIT_FAILURE);}

    //set scenario manually
    scenarioSpec = 5;

    // set paths for SCENARIO_NAME, SCENARIO_DAT_FILE and SCENARIO_PGM_FILE
    set_paths(scenarioSpec, SCENARIO_NAME, SCENARIO_DAT_FILE, SCENARIO_PGM_FILE);


    //check if directory "output" exists, if not creates directory "output"
    check_dir_exists(SCENARIO_NAME);

    FILE *parameterFile;
    FILE *geometryFile;

    const char *input_parameter_file_path = SCENARIO_DAT_FILE.c_str();
    const char *input_geometry_file_path = SCENARIO_PGM_FILE.c_str();

    parameterFile = fopen(input_parameter_file_path, "r");
    geometryFile = fopen(input_parameter_file_path, "r");
    
    //check whether plane_shear.dat exists or not
    if (parameterFile == NULL) {
        printf("Error opening parameter-file");
        exit(EXIT_FAILURE);
    } 
    else if (geometryFile == NULL) {
        printf("Error opening geometry-file");
        exit(EXIT_FAILURE);
    }
    else {
        std::string parameterFile{input_parameter_file_path}; //relative path to dat file
        //ready parameters from dat file and assign values to initalized parameters
        read_parameters(parameterFile, Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, imax, jmax, alpha, omg,
            tau, itermax, eps, dt_value, TI, T_h, T_c, Pr, beta, v_inflow, u_inflow, kappa, heat_flux, CI,
            C_inject, Pr_diffusion, SD_coeff, stoichiometric_coeff, homogeneous_reaction_coef, absorption_coeff, heat_capacity,
            reaction_rate_constant_factor, activation_energy_forward, activation_energy_reverse, activation_energy_catalyst,
            vacant_centers_defficiency_coeff, reaction_heat_effect_Q, is_product);
    }

    cell_array = read_pgm(input_geometry_file_path, *imax, *jmax);

    // Set up grid
    Grid grid(*imax, *jmax, BOUNDARY_SIZE, *PI, *UI, *VI, *TI, CI[A], CI[B], CI[C], CI[D]);

    if (!assert_problem_solvability(cell_array, grid)) {
        printf("PGM file is not solvable");
        exit(EXIT_FAILURE);
    }
 
    //for output to vtk-file
    VTKHelper vtkOutput;
    
    for (int j = grid.jmaxb() - 1; j >= 0; j--){
        for (int i = 0; i < grid.imaxb(); i++){
            //assign cell type
            if (cell_array[i][j] == 0){grid.cell(i,j)._cellType = NOSLIP;}
            else if(cell_array[i][j] == 4){grid.cell(i,j)._cellType = FLUID;}
            else if(cell_array[i][j] == 3){grid.cell(i,j)._cellType = OUTFLOW;}
            else if(cell_array[i][j] == 2){grid.cell(i,j)._cellType = INFLOW;}
            else if(cell_array[i][j] == 1){grid.cell(i,j)._cellType = FREESLIP;}
            else if(cell_array[i][j] == 5){grid.cell(i,j)._cellType = CATALYST;}
            else{
                printf("Error: wrong grey value in geometry-file ");
                exit(EXIT_FAILURE);
            }
        }
    }

    assign_ptr_nbcells(grid);

    // Initializing variables
    double time = 0;                        // time
    int timesteps_total = 0;                // # of iterations for main loop
    int current_timestep_iteration;         // # of iterations for SOR
    double visualization_time_accumulator = 0.0;        // signals when it's time to visualize within the main loop
    int count_failed_SOR = 0;               //# of failed SOR iterations

    //initialize matrices U, V, P, T, C
    matrix<double> U, V, P, T, C;
    init_uvptc(*imax, *jmax, U, V, P, T, *UI, *VI, *PI, *TI, CI, grid);

    //initialize matrices F, G and RS
    matrix<double> F, G, RS;

    //assign initial values FI, GI and RSI on the hole domain for F, G and RS
    init_fgrs(*imax, *jmax, F, G, RS, 0, 0, 0, grid);

    vtkOutput.printVTKFile(grid, *dx, *dy, SCENARIO_NAME, SCENARIO_NAME, timesteps_total);

    // Initialize timer to measure performance
    Timer runtime;

    while (time < *t_end) {
        //here we set time steps manually
        calculate_dt(*Re, *Pr, Pr_diffusion, *tau, dt, *dx, *dy, *imax, *jmax, grid);

        boundaryvalues(*imax, *jmax, grid, *v_inflow, *u_inflow, F, G, *dx, *dy, *beta, *dt, *GX, *GY);

        spec_boundary_val(*imax, *jmax, grid, *v_inflow, *u_inflow, *T_h, *T_c, C_inject, *dx, *dy, *kappa, *heat_flux, *beta, *dt, *GX, *GY, scenarioSpec, time, *t_end);
        
        for (int it = 0; it < LAST; it++)
            calculate_concentration(*Re, Pr_diffusion,
                *alpha, *dt, *dx, *dy, *imax, *jmax, grid, static_cast<ID>(it));

        calculate_chem_kinetics(*dt, *dx, *dy, *imax, *jmax, grid,
            is_product[0],
            is_product[1],
            is_product[2],
            is_product[3],
            stoichiometric_coeff[0],
            stoichiometric_coeff[1],
            stoichiometric_coeff[2],
            stoichiometric_coeff[3],
            homogeneous_reaction_coef[0],
            homogeneous_reaction_coef[1],
            homogeneous_reaction_coef[2],
            homogeneous_reaction_coef[3],
            absorption_coeff[0],
            absorption_coeff[1],
            absorption_coeff[2],
            absorption_coeff[3],
            heat_capacity[0],
            heat_capacity[1],
            heat_capacity[2],
            heat_capacity[3],
            *reaction_rate_constant_factor,
            *activation_energy_forward,
            *activation_energy_reverse,
            *activation_energy_catalyst,
            *SD_coeff,
            *vacant_centers_defficiency_coeff,
            *reaction_heat_effect_Q);

        calculate_temp(*Re, *Pr, *alpha, *dt, *dx, *dy, *imax, *jmax, grid);

        calculate_fg(*Re, *beta, *GX, *GY, *alpha, *dt, *dx, *dy, *imax, *jmax, grid, F, G);

        calculate_rs(*dt, *dx, *dy, *imax, *jmax, F, G, RS, grid);
        

        //reset current number of iterations for SOR
        current_timestep_iteration = 0;

        //reset residual before new SOR iteration
        *res = INFINITY;

        while ((*res > *eps) && (current_timestep_iteration <= *itermax)) {
            sor(*omg, *dx, *dy, *imax, *jmax, grid, RS, res);
            current_timestep_iteration++;
        }
        //count number of failed SOR iterations
        if(*res > *eps){
            //print warning message after failed SOR iteration
            //std::cout << "Warning: current #SOR iterations: " << current_timestep_iteration <<  " exceeded max #SOR iterations: " << *itermax << "!" << std::endl;
            count_failed_SOR++;
        }

        calculate_uv(*dt, *dx, *dy, *imax, *jmax, grid, F, G);
        visualization_time_accumulator += * dt;
        timesteps_total++;
        time += *dt;
        
        // Visualize u v p
        if (visualization_time_accumulator >= *dt_value) {

            //write_vtkFile(SCENARIO_NAME, timesteps_total, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P, T);
            vtkOutput.printVTKFile(grid, *dx, *dy, SCENARIO_NAME, SCENARIO_NAME, timesteps_total);
            solutionProgress(time, *t_end); // Print out total progress with respect to the simulation timerange
            visualization_time_accumulator -= *dt_value;
        }
    }


    //write_vtkFile(SCENARIO_NAME, timesteps_total, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P, T);
    vtkOutput.printVTKFile(grid, *dx, *dy, SCENARIO_NAME, SCENARIO_NAME, timesteps_total);


    // print Temperature in final timestep
    /*
    std::cout << "T temperature" << std::endl;
    for (int i = 0; i < grid.imaxb() / 2; ++i) {
        for (int j = 0; j < grid.jmaxb() / 2; ++j)
            std::cout << T[i][j] << " ";
        std::cout << std::endl;
    }
     */
          
    // Print out the total time required for the solution
    runtime.printTimer();

    //Print total number of timesteps and number of failed SOR iterations
    std::cout << "#total of timesteps: " << timesteps_total << " #failed SOR iterations: " << count_failed_SOR << std::endl;

    //close input file
    fclose(parameterFile);
    fclose(geometryFile);
   
    // Free dynamically allocated memory
    delete[] cell_array;
    delete Re;
    delete UI;
    delete VI;
    delete PI;
    delete GX;
    delete GY;
    delete t_end;
    delete xlength;
    delete ylength;
    delete dt;
    delete dx;
    delete dy;
    delete imax;
    delete jmax;
    delete alpha;
    delete omg;
    delete tau;
    delete itermax;
    delete eps;
    delete dt_value;
    delete TI;
    delete T_h;
    delete T_c;
    delete[] C_inject;
    delete[] CI;
    delete Pr;
    delete[] Pr_diffusion;
    delete res;
    delete beta;
    delete v_inflow;
    delete u_inflow;
    delete kappa;
    delete heat_flux;
    delete SD_coeff;
    delete[] stoichiometric_coeff;
    delete[] homogeneous_reaction_coef;
    delete[] absorption_coeff;
    delete[] heat_capacity;
    delete[] is_product;
    delete reaction_rate_constant_factor;
    delete activation_energy_forward ;
    delete activation_energy_reverse ;
    delete activation_energy_catalyst ;
    delete vacant_centers_defficiency_coeff ;
    delete reaction_heat_effect_Q ;
    return 0;
}

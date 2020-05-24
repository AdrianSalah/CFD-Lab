#include "helper.hpp"
#include "visual.hpp"
#include "init.hpp"
#include "sor.hpp"
#include "uvp.hpp"
#include <cstdio>
#include <iostream>
#include "boundary_val.hpp"
#include "Timer.h"
#include <assert.h>
/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed. Use the predefined matrix<typename> type and give initial values in the constructor.
 * - perform the main loop
 * - at the end: destroy any memory allocated and print some useful statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two-dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop, the following steps are required (for some of the 
 * operations, a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */

int main(int argn, char** args) {

    //initialize all relevant parameters
    double* Re = new double;                /* reynolds number   */
    double* UI = new double;                /* velocity x-direction */
    double* VI = new double;                /* velocity y-direction */
    double* PI = new double;                /* pressure */
    double* GX = new double;                /* gravitation x-direction */
    double* GY = new double;                /* gravitation y-direction */
    double* t_end = new double;             /* end time */
    double* xlength = new double;           /* length of the domain x-dir.*/
    double* ylength = new double;           /* length of the domain y-dir.*/
    double* dt = new double;                /* time step */
    double* dx = new double;                /* length of a cell x-dir. */
    double* dy = new double;                /* length of a cell y-dir. */
    int* imax = new int;                    /* number of cells x-direction*/
    int* jmax = new int;                    /* number of cells y-direction*/
    double* alpha = new double;             /* uppwind differencing factor*/
    double* omg = new double;               /* relaxation factor */
    double* tau = new double;               /* safety factor for time step*/
    int* itermax = new int;                 /* max. number of iterations  */
    double* eps = new double;               /* accuracy bound for pressure*/
    double* dt_value = new double;          /* time for output */
    double* TI = new double;                /* Initial Temperature*/
    double* T_h = new double;               /* Temperature of hot wall*/
    double* T_c = new double;               /* Temperature of cold wall*/
    double* PR = new double;                      /* Prandlt Number*/
    double* res = new double;               /* residual for SOR*/
    double* beta= new double;               /* beta for fg calculation*/
    int **cell_array = new int *;             /* array of geometry */
    //check if directory "output" exists, if not creates directory "output"
    check_dir_exists("output");


    FILE *parameterFile;
    FILE *geometryFile;

    const char *input_parameter_file_path = "../cavity100.dat";
    const char *input_geometry_file_path = "../geometry.pgm";

    parameterFile = fopen(input_parameter_file_path, "r");
    geometryFile = fopen(input_parameter_file_path, "r");

    //check whether cavity100.dat exists or not
    if (parameterFile == NULL) {
        printf("Error opening parameter-file");
        exit(EXIT_FAILURE);
    } else if (geometryFile == NULL) {
        printf("Error opening geometry-file");
        exit(EXIT_FAILURE);
    } else {
        std::string parameterFile{input_parameter_file_path}; //relative path to cavity100.dat file
        //ready parameters from cavity100.dat file and assign values to initalized parameters
        read_parameters(parameterFile, Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, imax, jmax, alpha, omg,
                        tau, itermax, eps, dt_value, TI, T_h, T_c, PR, beta);
    }

    cell_array = read_pgm(input_geometry_file_path);

    if (!assert_problem_solvability(cell_array, *imax, *jmax)) {
        printf("PGM file is not solvable");
        exit(EXIT_FAILURE);
    }

    //for output to vtk-file
    VTKHelper vtkOutput;
    
    // Set up grid
    Grid grid(*imax, *jmax, 1, *PI, *UI, *VI, *TI);

    //TO DO: check wheather imax and jmax same as grid size in geometry file
    for (int j = *jmax-1; j >= 0; j--){
        for (int i = 0; i < *imax; i++){
            //assign cell type
            if (cell_array[i][j] == 0){grid.cell(i,j)._cellType = NOSLIP;}
            else if(cell_array[i][j] == 4){grid.cell(i,j)._cellType = FLUID;}
            else if(cell_array[i][j] == 3){grid.cell(i,j)._cellType = OUTFLOW;}
            else if(cell_array[i][j] == 2){grid.cell(i,j)._cellType = INFLOW;}
            else if(cell_array[i][j] == 1){grid.cell(i,j)._cellType = FREESLIP;}
            else{
                printf("Error: wrong grey value in geometry-file ");
                exit(EXIT_FAILURE);
            }
            //for debugging
            std::cout << cell_array[i][j] << " ";
        }
        //for debugging
        std::cout << std::endl;
    }

    //store pointers to neighbours for inner cells
    // and check if the cell is a FLUID-cell. Increment by 1, if true.
    for(int j = 1; j < *jmax-1; j++) {
        for (int i = 1; i < *imax - 1; i++) {
            grid.cell(i, j)._nbNorth = &grid.cell(i, j + 1);
            grid.cell(i, j)._nbEast = &grid.cell(i + 1, j);
            grid.cell(i, j)._nbWest = &grid.cell(i - 1, j);
            grid.cell(i, j)._nbSouth = &grid.cell(i, j - 1);
            
            // Can it be moved to grid cpp?
            if (grid.cell(i, j)._cellType == FLUID)
                grid.increment_fluid_cells();
        }
    }

    // Checking immediately if the number of FLUID cells is not zero
    assert(grid.get_fluid_cells_quantity() > 0 && "There are no FLUID cells in the domain");

    //neighbour edges
    //bottom left
    grid.cell(0,0)._nbNorth = &grid.cell(0,1);
    grid.cell(0,0)._nbEast = &grid.cell(1,0);
    //top right
    grid.cell(*imax-1, *jmax-1)._nbSouth = &grid.cell(*imax-1,*jmax-2);
    grid.cell(*imax-1, *jmax-1)._nbWest = &grid.cell(*imax-2,*jmax-1);
    //top left
    grid.cell(0,*jmax-1)._nbEast = &grid.cell(1, *jmax-1);
    grid.cell(0, *jmax-1)._nbSouth = &grid.cell(0, *jmax-2);
    //bottom right
    grid.cell(*imax-1, 0)._nbNorth = &grid.cell(*imax-1, 1);
    grid.cell(*imax-1, 0)._nbWest = &grid.cell(*imax-2, 0);

    for(int i = 1; i < *imax-1; i++){
        //bottom
        grid.cell(i, 0)._nbNorth = &grid.cell(i, 1);
        grid.cell(i, 0)._nbEast = &grid.cell(i + 1, 0);
        grid.cell(i, 0)._nbWest = &grid.cell(i-1, 0);
        //top
        grid.cell(i, *jmax-1)._nbSouth  = &grid.cell(i,*jmax-2);
        grid.cell(i, *jmax-1)._nbEast = &grid.cell(i+1, *jmax-1);
        grid.cell(i, *jmax-1)._nbWest = &grid.cell(i-1, *jmax-1);
    }
    for(int j = 1; j <  * jmax-1; j++){
        //left
        grid.cell(0, j)._nbNorth = &grid.cell(0, j+1);
        grid.cell(0, j)._nbSouth = &grid.cell(0, j-1);
        grid.cell(0, j)._nbEast = &grid.cell(1, j);
        //right
        grid.cell(*imax-1, j)._nbNorth = &grid.cell(*imax-1, j+1);
        grid.cell(*imax-1, j)._nbSouth = &grid.cell(*imax-1, j-1);
        grid.cell(*imax-1, j)._nbWest = &grid.cell(*imax-2, j);
    }

    Cell testCell = grid.cell(1,11);
    //displays positions of neighbouring cells of test-cell (i=1, j=11)
    std::cout << testCell._cellType << std::endl;
    std::cout << "W " << testCell._nbWest->_cellType << std::endl;
    std::cout << "E " << testCell._nbEast->_cellType << std::endl;
    std::cout << "N " << testCell._nbNorth->_cellType << std::endl;
    std::cout << "S " << testCell._nbSouth->_cellType << std::endl;


    // Initializing variables
    double time = 0;                        // time
    int timesteps_total = 0;                // # of iterations for main loop
    int current_timestep_iteration;         // # of iterations for SOR
    double visualization_time_accumulator = 0.0;        // signals when it's time to visualize within the main loop
    int count_failed_SOR = 0;               //# of failed SOR iterations

    matrix<double> U, V, P, T;
    init_uvpt(*imax, *jmax, U, V, P, T, *UI, *VI, *PI, *TI, grid);
    //initialize matrices F, G and RS
    matrix<double> F, G, RS;

    //assign initial values FI, GI and RSI on the hole domain for F, G and RS
    init_fgrs(*imax, *jmax, F, G, RS, 0, 0, 0);

    //initialize matrices U, V, P

    // Initialize timer to measure performance
    Timer runtime;

    while (time < *t_end) {
        //here we set time steps manually
        calculate_dt(*Re, *PR, *tau, dt, *dx, *dy, *imax, *jmax, grid);
        boundaryvalues(*imax, *jmax, grid, F, G);
        calculate_fg(*Re, *GX, *GY, *alpha, *dt, *dx, *dy, *imax, *jmax, grid, F, G);
        calculate_rs(*dt, *dx, *dy, *imax, *jmax, F, G, RS);

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
            grid.velocity(U, velocity_type::U);
            grid.velocity(V, velocity_type::V);
            grid.pressure(P);
            //write_vtkFile("cavityData", timesteps_total, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P);
            vtkOutput.printVTKFile(grid, *dx, *dy, "testcase", "output", timesteps_total);
            solutionProgress(time, *t_end); // Print out total progress with respect to the simulation timerange
            visualization_time_accumulator -= *dt_value;
        }
    }

    grid.velocity(U, velocity_type::U);
    grid.velocity(V, velocity_type::V);
    grid.pressure(P);


    //write_vtkFile("cavityData", timesteps_total, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P);
    vtkOutput.printVTKFile(grid, *dx, *dy, "testcase", "output", timesteps_total);


    // Print out the total time required for the solution
    runtime.printTimer();

    //Print total number of timesteps and number of failed SOR iterations
    std::cout << "#total of timesteps: " << timesteps_total << " #failed SOR iterations: " << count_failed_SOR << std::endl;


    //close input file
    fclose(parameterFile);
    fclose(geometryFile);

    
    // Free dynamically allocated memory
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
    delete PR;
    delete res;
    delete beta;

    return 0;
}

#ifndef __PARAMETERS_H_
#define __PARAMETERS_H_

#define BOUNDARY_SIZE 1

int scenarioSpec;
std::string SCENARIO_NAME;
std::string SCENARIO_DAT_FILE;
std::string SCENARIO_PGM_FILE;

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
double* CI = new double[LAST];                /* Initial Concentration */
double* C_inject = new double[LAST];       /* Concentration of the injected chemical */
double* Pr = new double;                /* Prandlt Number*/
double* Pr_diffusion = new double[LAST];      /* Prandlt Number for chemical diffusion*/
double* res = new double;               /* residual for SOR*/
double* beta = new double;               /* beta for fg calculation*/
double* v_inflow = new double;          /* boundary value for inflow BC */
double* u_inflow = new double;          /* boundary value for inflow BC */
double* kappa = new double;             /* thermal conductivity */
double* heat_flux = new double;         /* heat flux */
double* SD_coeff = new double;         /* surface development coefficient */
double* stoichiometric_coeff = new double[LAST];
double* homogeneous_reaction_coef = new double[LAST];
double* absorption_coeff = new double[LAST];
double* heat_capacity = new double[LAST];
bool* is_product = new bool[LAST];
double* reaction_rate_constant_factor = new double;
double* activation_energy_forward = new double;
double* activation_energy_reverse = new double;
double* activation_energy_catalyst = new double;
double* vacant_centers_defficiency_coeff = new double;
double* reaction_heat_effect_Q = new double;
int** cell_array = new int*;           /* array of geometry */

#endif
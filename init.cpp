#include "init.hpp"
#include <fstream>

int read_parameters(std::string szFileName,       /* name of the file */
	double* Re,                /* reynolds number   */
	double* UI,                /* velocity x-direction */
	double* VI,                /* velocity y-direction */
	double* PI,                /* pressure */
	double* GX,                /* gravitation x-direction */
	double* GY,                /* gravitation y-direction */
	double* t_end,             /* end time */
	double* xlength,           /* length of the domain x-dir.*/
	double* ylength,           /* length of the domain y-dir.*/
	double* dt,                /* time step */
	double* dx,                /* length of a cell x-dir. */
	double* dy,                /* length of a cell y-dir. */
	int* imax,                /* number of cells x-direction*/
	int* jmax,                /* number of cells y-direction*/
	double* alpha,             /* uppwind differencing factor*/
	double* omg,               /* relaxation factor */
	double* tau,               /* safety factor for time step*/
	int* itermax,             /* max. number of iterations  */
	double* eps,                 /* accuracy bound for pressure*/
	double* dt_value,            /* time for output */
	double* TI,                   /* Initial temperature */
	double* T_h,                  /* hot wall temperature */
	double* T_c,                  /* cold wall temperature */
	double* PR,                      /*Prandlt Number*/
	double* beta,
	double* v_inflow,          /* inflow y velocity */
	double* u_inflow,          /* inflow x velocity */
	double* kappa,				/* thermal conductivity */
	double* heat_flux,		 /* heat flux */
	double* CI,             /* initial concentration */
	double* C_inject,
	double* Pr_diffusion,            /* Prandlt Number for chemical diffusion*/
	double* SD_coeff,                  /* surface development coefficient */
	double* stoichiometric_coeff,
	double* homogeneous_reaction_coef,
	double* absorption_coeff,
	double* heat_capacity,
	double* reaction_rate_constant_factor,
	double* activation_energy_forward,
	double* activation_energy_reverse,
	double* activation_energy_catalyst,
    double* vacant_centers_defficiency_coeff,
    double* reaction_heat_effect_Q,
	bool* is_product
                     )                     
{
    // Reading Parameters
	std::ifstream file(szFileName);
	if (!file.is_open()) return -1;
	std::string var;
	while (!file.eof() && file.good()) {
		file >> var;
		if (var[0] == '#') {     /* ignore comment line*/
			file.ignore(MAX_LINE_LENGTH, '\n');
		}
		else {
			if (var == "xlength")  file >> *xlength;
			if (var == "ylength")  file >> *ylength;
			if (var == "Re")       file >> *Re;
			if (var == "t_end")    file >> *t_end;
			if (var == "dt")       file >> *dt;
			if (var == "omg")      file >> *omg;
			if (var == "eps")      file >> *eps;
			if (var == "tau")      file >> *tau;
			if (var == "alpha")    file >> *alpha;
			if (var == "dt_value") file >> *dt_value;
			if (var == "UI")       file >> *UI;
			if (var == "VI")       file >> *VI;
			if (var == "GX")       file >> *GX;
			if (var == "GY")       file >> *GY;
			if (var == "PI")       file >> *PI;
			if (var == "itermax")  file >> *itermax;
			if (var == "imax")     file >> *imax;
			if (var == "jmax")     file >> *jmax;
			if (var == "TI")       file >> *TI;
			if (var == "T_h")      file >> *T_h;
			if (var == "T_c")      file >> *T_c;
			if (var == "PR")       file >> *PR;
			if (var == "beta")     file >> *beta;
			if (var == "u_inflow") file >> *u_inflow;
			if (var == "v_inflow") file >> *v_inflow;
			if (var == "kappa")    file >> *kappa;
			if (var == "heat_flux")file >> *heat_flux;
			if (var == "SD_coeff")		  file >> *SD_coeff;

			if (var == "CI") {
				file >> CI[0];
				file >> CI[1];
				file >> CI[2];
				file >> CI[3];
			}
			if (var == "C_inject") {
				file >> C_inject[0];
				file >> C_inject[1];
				file >> C_inject[2];
				file >> C_inject[3];
			}
			if (var == "Pr_diffusion") {
				file >> Pr_diffusion[0];
				file >> Pr_diffusion[1];
				file >> Pr_diffusion[2];
				file >> Pr_diffusion[3];
			}

			if (var == "stoichiometric_coeff") { 
				file >> stoichiometric_coeff[0];
				file >> stoichiometric_coeff[1];
				file >> stoichiometric_coeff[2];
				file >> stoichiometric_coeff[3];
			}
			if (var == "homogeneous_reaction_coef") {
				file >> homogeneous_reaction_coef[0];
				file >> homogeneous_reaction_coef[1];
				file >> homogeneous_reaction_coef[2];
				file >> homogeneous_reaction_coef[3];
			}
			if (var == "absorption_coeff") {
				file >> absorption_coeff[0];
				file >> absorption_coeff[1];
				file >> absorption_coeff[2];
				file >> absorption_coeff[3];
			}
			if (var == "heat_capacity") {
				file >> heat_capacity[0];
				file >> heat_capacity[1];
				file >> heat_capacity[2];
				file >> heat_capacity[3];
			}
			if (var == "is_product") {
				file >> is_product[0];
				file >> is_product[1];
				file >> is_product[2];
				file >> is_product[3];
			}
			if (var == "activation_energy_forward")				file >> *activation_energy_forward;
			if (var == "reaction_rate_constant_factor")		    file >> *reaction_rate_constant_factor;
			if (var == "activation_energy_reverse")				file >> *activation_energy_reverse;
			if (var == "activation_energy_catalyst")		    file >> *activation_energy_catalyst;
			if (var == "vacant_centers_defficiency_coeff")		file >> *vacant_centers_defficiency_coeff;
			if (var == "reaction_heat_effect_Q")				file >> *reaction_heat_effect_Q;

		}
	}
	*dx = *xlength / (double)(*imax);
	*dy = *ylength / (double)(*jmax);

	if (!file.good() && !file.eof()) return -1;

	return 1;
}


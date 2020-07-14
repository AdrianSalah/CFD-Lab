# CFD Lab Project of Group C  
## Chemical transport and catalytic chemical reactions

Our code simulates chemical transport as well as catalytic chemical reactions of four chemical components \\(A\\), \\(B\\), \\(C\\) and \\(D\\). We consider following reactions:
\\[ A + B \rightleftharpoons C  \\]
i.e. the Haber process: 
\\[ N_{2} + 3H_{2} \rightleftharpoons 2NH_{3}  \\]

Two additional scenarios "**catalyst_reactor**" and "**catalyst_reactor2**" where introduced to show the capability of the code. The parameterfile was extended with following parameters which can be set for each component individually
  
  Parametername | Explanation
  ------------- | :-------------:
  is_product	| 1: reaction product, 0: no reaction product, -1: not considered
  SD_coeff      | the surface development coefficient
  CI            | Initial concentrations for all components
  stoichiometric_coeff | Stoichiometric coefficients for the chemical reaction
  C_inject | Injected concentration on inflow
  Pr_diffusion | the prandlt number for chemical transport
  homogeneous_reaction_coef | Reaction coefficients for the homogenous reaction 
  absorption_coeff | Absorption coefficients for each component
  heat_capacity | heat capacity for each component
  
as well as thermodynamic parameters:
  
  Parametername | Explanation
  ------------- | -------------
  reaction_rate_constant_factor	| Constant factor in the rate equation
  activation_energy_forward      | Activation energy for the forward reaction
  activation_energy_reverse      | Activation energy for the reverse reaction
  activation_energy_catalyst | Activation energy for the catalyst reaction
  vacant_centers_defficiency_coeff | ??
  reaction_heat_effect_Q | Released heat by exothermic reaction

 You can find the parameterfile e.g. "***catalyst_reactor.dat***" and "***catalyst_reactor2.dat***" for the specific scenario in the subdirectory "***parameters***".
 In the subdirectory "***geometry***" the predefined pgm-files "***catalyst_reactor.pgm***" and "***catalyst_reactor2.pgm***" can be found. Both geometries contain an inflow boundary condition on the bottom, an outflow boundary condition on the top as well as adiabatic no-slip walls. An additional celltype was added to the PGM-file to model the catalyst. The first scenario contains vertical catalyst blocks inside the reactor with a coarser grid. The second one horizontal catalyst blocks attached to the walls reaching inside the catalyst reactor with a finer grid.

**The Algorithm**
<br>
begin <br>
Read parameters also initial concentrations CI and thermodynamic parameters, new celltype for catalysator <br>
Set t := 0, n := 0 <br>
read pgm file and check it <br>
Assign initial values to u, v, p, T and CA, CB, CC, CD <br>
  while t < tend do <br>
 	select δt with new stability criterion for C <br>
	set boundary values for u, v, P, T and C for catalyst reactor <br>
    compute new concentrations of the 4 components because of chemical transport <br>
	compute new concentrations of components an update temperature depending on reaction type <br>
 	compute new temperature <br>
	compute Fn and Gn according to the modified versions <br>
	compute the right-hand side of the PPE  <br>
 	implicit step: solve the PPE using SOR -> pn+1 <br>
 	explicit step: update un, vn -> un+1, vn+1  <br>
 	t := t + t <br>
 	n := n + 1 <br>
    output of u, v, P, T and CA, CB, CC, CD values and geometry for visualization periodically <br>
  end <br>
output of u, v, P, T and CA, CB, CC, CD values and geometry for visualization <br>
end <br>
<br>
**How to run the code for the project scenarios**

Please compile code and then type: `./sim scenarionumber`  in the command line. To run different scenarios please change `scenarionumber` for an integer that refers to a certain scenario:
* **`scenarionumber = 8`  Cylinder**
* **`scenarionumber = 9`  catalyst reactor**
* **`scenarionumber = 10` catalyst reactor 2**
* **`scenarionumber = 11` catalyst reactor 3**

The scenarions from worksheet 1 and 2 can be executed as well using scenarionumbers above:

Scenarioname  | Scenarionumber
  ------------- | -------------
  Lid driven cavity  | **1**
  Plane shear flow  | **2**
  Karman vortex street  | **3**
  Flow over step  | **4**
  Natrual convection  | **5**
  Fluid trap  | **6**
  Rayleigh benard convection  | **7**


**Output for visualization**

A folder with the respective scenario name is created to store the VTK-files for visualization. 
If no scenario is specified it will run the code for the lid-driven cavity. 
If you pass a parameter `scenarionumber` which is not equal to 1-11 the programm will terminate with an error.
A check is done if the provided geometry file is solvable. Please note one obstacle cell (marked as 0 in the pgm-file) is only allowed to have at maximum two adjacent fluid neighbour cells (marked as 4 in the pgm-file).
The scenario name and the progress of the computation is displayed in the terminal. 
If you encounter strange looking visualisation results please change between channels in Paraview.

**Representation of the cells in the PGM-file**

Cell-Type| Number in PGM-File
  -------------  | -------------
      Fluid Cell | **4**
  No-slip Cell   | **0**
  Free-slip Cell | **1**
      Inflow Cell| **3**
     Outflow Cell| **2**
  Catalyst Cell  | **5**


**This repository contains:**

* directory "***geometry***" which contains PGM-files which are used to read in the geometry for the different scencarios
* directory "***parameters***" which contains DAT-files which are used to read in the parameters for the different scencarios
* an executable "***sim***"
* the headers
* the files with the respective method implementations 
* example files of created visualisations in directory "***visualization_examples***"


## Software Requirements

* VTK 7 or higher
* GCC 9 (optional) 

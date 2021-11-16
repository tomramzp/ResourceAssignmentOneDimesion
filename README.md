# Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads

The following code allows replicating the results from the publication:
 [Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads](https://arxiv.org/abs/2109.09385) 
 Tomás Ramírez, Carlos Mosquera, Nader Alagha
 arXiv, 2021

# SOFTWARE REQUIREMENTS

- MATLAB 2019b or newer version is required. The following MATLAB toolboxes are required:
	- Global Optimization toolbox
	- Parallel Computing Toolbox
	- Optimization Toolbox
	- Statistics and Machine Learning Toolbox

- Additionally, two external packages are also required, MOSEK and CVX: 

	- MOSEK ApS is a software for nonlinear convex optimization. The software requires a user license which can be obtained free of charge for researchers, students or other academic members. More information about the license can be found at https://www.mosek.com/products/academic-licenses/.  MOSEK is employed through CVX ( that is detailed below). It should be reminded to follow the steps in order as described in the guide. In particular,  Step 4  is required to detect the MOSEK licence with CVX, after placing the MOSEK license in the correct folder. 

	- CVX is a powerful and efficient solver for convex optimization problems. A user guide and installation steps can be found in http://web.cvxr.com/cvx/doc/mosek.html
       	  to install both MOSEK and CVX in MATLAB.


# SCRIPTS

		

- "System_Simulations.m": The script simulates different resource assignment strategies for a one-dimensional satellite scenario with a two-colour scheme. A parent folder is created for the storing of the numerical results. This folder is labelled with the time and date of the script execution. Subfolders are created for the considered                       traffic profiles.

 
				INPUT:

				At the beginning of the script, multiple settings can be configured:

					- Section "Monte-Carlo parameters":  Number of Monte-Carlos simulations for the numerical results. For example, the following line of code defines 500 Monter-Carlos simulations:

					  	Nsims_v= 500;

					  If we want to store the data in small batches, for the same total of Monte-Carlo, we can modify the previous line of code to:
				         
					  	Nsims_v= ones(1,100)*5;

					  In that case, 100 batch simulations are stored with 5 Monte-Carlo simulations per batch.

					- "Section  "Resource allocation strategies": The script allows the simulation of the resource strategies from "Flexible User Mapping for Radio Resource Assignment
					  in Advanced Satellite Payloads". Solutions can be enabled or disabled with a boolean variable:

						en_POW: Enables ( with value equal to 1) the optimization of flexible power allocation to cope with the traffic demand. 
						en_BW: Enables ( with value equal to 1) the optimization of flexible bandwidth allocation to cope with the traffic demand.
						en_MA: Enables ( with value equal to 1) the optimization of flexible beam-user mapping (with fixed resources) to cope with the traffic demand. 
						en_BW_MAP: Enables( with value equal to 1) the joint optimization of bandwidth and beam-user mapping to cope with the traffic demand. 
						en_BW_POW: Enables ( with value equal to 1) the joint optimization of bandwidth and power to cope with the traffic demand. 

					- Section "One Dimension satellite scenario": The script allows the configuration of :

					    K : Number of beams ( it must be even)
						R : Beam radius with the Bessel modelling
						beam_cross_roll_off : Beam roll-off.

					- Section "Traffic distribution": Multiple traffic distributions can be defined. Dirichlet distribution models the traffic distribution across beams.

					  By default, the three traffic scenarios from "Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads" are defined. The following 
					  pseudocode allows the introduction of a scenario:
					   
 					  	alpha_v{k}= Vector with the alpha values of the Dirichlet distribution
 					  	label_alpha{k}= Scenario Label
 					  	mkdir(label_alpha{k}), % Create folder

				        Each defined scenario must be associated with a different scenario identifier "k" within the script. 		Also, a scenario label is defined to store the numerical results
					  and to generate later the figures with the performances of the techniques. 

					- Section "System parameters": The script allows the configuration of:

						M : Number of carriers per color ( Half the bandwidth with Two-Color scheme)
						total_band: Total available  bandwidth.
						Req_user_ref: Requested user traffic ( equal to all the users). Normalized to the available  total bandwidth. 



					- Section "Satellite parameters": The script allows the configuration of different parameters for the link budget ( Power budget, transmiter antenna gain, etc).



		        	OUTPUT: 

				A parent folder is created for the storing of the numerical results. The numerical results are stored in the corresponding subfolders of the defined scenarios. Furthermore,  different
                		information about the status of the simulation process is displayed through the command window. 

				Please, ignore the message "NOTE: custom settings have been set for this solver" that is displayed due to the custom settings of MOSEK within CVX.


- "PlotResults.m":  The script displays the average performance in the console window and plots the cumulative distribution functions of different performance metrics for a selected scenario.


				INPUT: 
				With the execution of this script, the user can select a scenario subfolder with the stored numerical results.
					
                		If the "System_Simulations.m" is executed with the SAME INPUT CONFIGURATION in different PC or time instants, the obtained numerical results are stored in different parent 
                		folders with same subfolders labels.  Only in that case where the SAME INPUT CONFIGURATION has been used in the system simulations, all the data for each scenario 
                		can be put together in a same scenario subfolder for a better averaging  of the performances.
                
                		PLEASE, DO NOT MIX NUMERICAL RESULTS FROM DIFFERENT SCENARIOS IN THE SAME FOLDER or NUMERICAL  RESULTS FROM THE SAME SCENARIO ( traffic distribution) BUT FOR DIFFERENT INPUT 
                		PARAMETERS.

				OUTPUT:
			        Average performances and cumulative distribution functions (CDF) are displayed. The following performance metrics are consided:

				- Normalized Quadratic Unmet capacity (NQU)
				- Normalized Unmet Capacity (NU)
				- Total offered rate [Mbps]
				- User offferd rate [Mbps]
				- Minimum user rate [Mbps]

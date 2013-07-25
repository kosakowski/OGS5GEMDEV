
#include "ca_hydration.h"
#include <math.h>
#include <cmath>

#ifndef max
  #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
  #define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#define COMP_MOL_MASS_WATER 0.018016 //TN
#define COMP_MOL_MASS_N2   0.028013 //TN

ca_hydration::ca_hydration(double T_solid, 
	                       double T_gas,  
						   double p_gas,
						   double w_water, 
						   double rho_s_initial,
						   double phi_S,
						   double delta_t,
						   std::string system)
:R(8.314510) // , n_col(3)
{
	x = Eigen::VectorXd(1);
	update_param( T_solid, T_gas, p_gas, w_water, rho_s_initial, phi_S, delta_t, system);
	ca_hydration::rho_s_0 = rho_s_initial;
	if (system.compare("CaOH2") == 0){ //Definition auch in void CSolidProperties::SetSolidReactiveSystemProperties()
		rho_low = 1655.98;
		rho_up = 2200.02;
		reaction_enthalpy = -1.12e+05; //in J/mol; negative for exothermic composition reaction
		reaction_entropy = -143.5; //in J/mol K
	}
	
}


ca_hydration::~ca_hydration(void)
{
}

void ca_hydration::update_param(double T_solid, 
	                            double T_gas,  
								double p_gas,
								double w_water, 
								double rho_s_initial,
								double phi_S,
								double delta_t,
								std::string system)
{
	ca_hydration::T_s     = T_solid;
	ca_hydration::T       = T_gas; 
	ca_hydration::p_gas   = p_gas; // should be in unit bar
	ca_hydration::w_h2o   = w_water; 
	ca_hydration::rho_s   = rho_s_initial; 
	x(0)                  = rho_s_initial;
	ca_hydration::phi_solid = phi_S;
	ca_hydration::dt      = delta_t;
	ca_hydration::reaction_system = system;
}

void ca_hydration::calculate_qR()
{
	// step 1, calculate X_D and X_H
	X_D = (rho_s - rho_up)/(rho_low - rho_up) ;
	qR = 0.0;
	if ( X_D < 0.00001 )
	{
		X_D = max(X_D, 0.00001);
		qR = 0.0; 
		return;
	}
	else if (X_D > 0.99999)
	{
		X_D = min(X_D, 0.99999);
		qR = 0.0;
		return;
	}

	X_H = 1.0 - X_D;

	//Convert mass fraction into mole fraction
	double mol_frac_vapor;
	mol_frac_vapor = COMP_MOL_MASS_N2*ca_hydration::w_h2o/(COMP_MOL_MASS_N2*ca_hydration::w_h2o + COMP_MOL_MASS_WATER*(1.0-ca_hydration::w_h2o)); //calculate mole fraction of vapor from mass fraction - TN
	
	////// step 2, calculate p_eq
	T_s = ca_hydration::T_s;
	p_w_g = mol_frac_vapor * ca_hydration::p_gas; //TN
	if (p_w_g < 1e-3) //TN - avoid illdefined log 
	{
		qR = 0.0;
		return;
	}

	// using the p_eq to calculate the T_eq - Clausius-Clapeyron
	T_eq = (reaction_enthalpy/R) / ((reaction_entropy/R) + log(p_w_g)); // unit of p in bar

	// step 3, calculate dX/dt
	if ( T_s < T_eq ) // hydration
	//if ( p_w_g > p_eq ) // hydration
	{
		/* this is from Schaube
		if ( (T_eq-T) >= 50.0 )
			//dXdt = 13945.0 * exp(-89486.0/R/T) * pow(p_w_g/p_eq - 1.0,0.83) * 3.0 * (1-X_H) * pow(-1.0*log(1.0-X_H),0.666); 
			dXdt = 13945.0 * exp(-89486.0/R/T) * pow(p_w_g/p_eq - 1.0,0.83) * 3.0 * (1-X_H);
		else
			// dXdt = 1.0004e-34 * exp(53.332e3/T) * pow(p_w_g/1.0e5, 6.0) * 2.0 * pow(1.0-X_H, 0.5); 
			// dXdt = 1.0e-308;// 1.0004e-34 * exp(5.3332e4/T) * pow(p_w_g, 6.0) * 2.0 * pow(1.0-X_H, 0.5); 
			dXdt = 1.0004e-34 * exp(5.3332e4/T) * pow(p_w_g, 6.0) * 2.0 * pow(1.0-X_H, 0.5); 
		*/
		// this is from P. Schmidt
		dXdt = -1.0*(1.0-X_H) * (T_s - T_eq) / T_eq;
		//dXdt = 13945.0 * exp(-89486.0/R/T) * pow(p_w_g/p_eq - 1.0,0.83) * 3.0 * (1-X_H) * pow(-1.0*log(1.0-X_H),0.666); 

		// step 4, calculate qR
	}
	else // dehydration
	{
		//dXdt = 1.2627e10 * exp( -1.6407e5/R/T )*pow(1.0-p_w_g/p_eq,2.9)*2.0*pow(1.0 - X_D, 0.5); 
		//dXdt = 13945.0 * exp(-89486.0/R/T) * pow(p_w_g/p_eq - 1.0,0.83) * 3.0 * (1-X_H) * pow(-1.0*log(1.0-X_H),0.666); //TN - test continuity
		//dXdt = 0.0;
	//	// step 4, calculate qR
	//	
		dXdt = -1.0* (1.0-X_D) * (T_s - T_eq) / T_eq;
	}

	qR = (rho_up - rho_low) * dXdt; //TN - reaction rate continuous around T_eq
	//TN - scale qR with mass fraction (smoothes simulation)
	if (qR > 0.0)
		qR *= ca_hydration::w_h2o;

	double k_R;

	if (qR < 0.0)
		k_R = 0.05;
	else
		k_R = 0.2;

	qR *= k_R;

	//qR *= log(ca_hydration::phi_solid)/log(0.2); //TN slows down reaction in less porous media and speads it up with increasing porosity
}

void ca_hydration::set_rho_s(double new_rho_s)
{
	rho_s = new_rho_s; 
}

double ca_hydration::get_qR()
{
	return qR; 
}

void ca_hydration::get_x(Eigen::VectorXd& output_x)
{
	output_x = x; 
}


void ca_hydration::eval(double t, Eigen::VectorXd &y, Eigen::VectorXd &dydx)
{
	assert( y.size() == dydx.size() );

	this->set_rho_s( y(0) ); 
	this->calculate_qR();
	dydx(0) = this->get_qR();  

}
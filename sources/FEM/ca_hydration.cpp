
#include "ca_hydration.h"
#include <math.h>
#include <cmath>

#define SIMPLE_KINETICS //wenn definiert, dann einfache Kinetik, sonst Schaube

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
						   FiniteElement::SolidReactiveSystem system)
:R(8.314510),p_eq(1.0) // , n_col(3)
{
	x = Eigen::VectorXd(1);
	update_param( T_solid, T_gas, p_gas, w_water, rho_s_initial, phi_S, delta_t, system);
	ca_hydration::rho_s_0 = rho_s_initial;
	if (system == FiniteElement::CaOH2){ //Definition auch in void CSolidProperties::SetSolidReactiveSystemProperties()
		rho_low = 1656.0;
		rho_up = 2200.0;
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
								FiniteElement::SolidReactiveSystem system)
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
	const double tol_l = 1.0e-4;
	const double tol_u = 1.0 - tol_l;
	const double tol_rho = 0.1;
	// step 1, calculate X_D and X_H
	X_D = (rho_s - rho_up - tol_rho)/(rho_low - rho_up - 2.0*tol_rho) ;
	qR = 0.0;

	X_D = (X_D < 0.5) ? max(tol_l,X_D) : min(X_D,tol_u); //constrain to interval [tol_l;tol_u]

	X_H = 1.0 - X_D;

	//Convert mass fraction into mole fraction
	double mol_frac_vapor;
	mol_frac_vapor = COMP_MOL_MASS_N2*ca_hydration::w_h2o/(COMP_MOL_MASS_N2*ca_hydration::w_h2o + COMP_MOL_MASS_WATER*(1.0-ca_hydration::w_h2o)); //calculate mole fraction of vapor from mass fraction - TN
	
	////// step 2, calculate equilibrium
	T_s = ca_hydration::T_s;
	p_w_g = max(mol_frac_vapor * ca_hydration::p_gas, 1.0e-3); //TN - avoid illdefined log

	// using the p_eq to calculate the T_eq - Clausius-Clapeyron
	T_eq = (reaction_enthalpy/R) / ((reaction_entropy/R) + log(p_w_g)); // unit of p in bar
	//Alternative: Use T_s as T_eq and calculate p_eq - for Schaube kinetics
	p_eq = exp((reaction_enthalpy/R)/T_s - (reaction_entropy/R));
	

	// step 3, calculate dX/dt
#ifdef SIMPLE_KINETICS
	if ( T_s < T_eq ) // hydration - simple model
#else
	if ( p_w_g > p_eq ) // hydration - Schaube model
#endif
	{
		//X_H = max(tol_l,X_H); //lower tolerance to avoid oscillations at onset of hydration reaction. Set here so that no residual reaction rate occurs at end of hydration.
#ifdef SIMPLE_KINETICS // this is from P. Schmidt
		dXdt = -1.0*(1.0-X_H) * (T_s - T_eq) / T_eq * 0.2 * ca_hydration::w_h2o;
#else //this is from Schaube
		if (X_H == tol_u || rho_s == rho_up)
			dXdt = 0.0;
		else if ( (T_eq-T_s) >= 50.0)
			dXdt = 13945.0 * exp(-89486.0/R/T_s) * pow(p_w_g/p_eq - 1.0,0.83) * 3.0 * (1.0-X_H) * pow(-1.0*log(1.0-X_H),0.666);
		else
			dXdt = 1.0004e-34 * exp(5.3332e4/T_s) * pow(p_w_g, 6.0) * (1.0-X_H);
#endif
	}
	else // dehydration
	{
		//X_D = max(tol_l,X_D); //lower tolerance to avoid oscillations at onset of dehydration reaction. Set here so that no residual reaction rate occurs at end of dehydration.
#ifdef SIMPLE_KINETICS // this is from P. Schmidt
		dXdt = -1.0* (1.0-X_D) * (T_s - T_eq) / T_eq * 0.05;
#else
		if (X_D == tol_u || rho_s == rho_low)
			dXdt = 0.0;
		else if (X_D < 0.2)
			dXdt = -1.9425e12 * exp( -1.8788e5/R/T_s )*pow(1.0-p_w_g/p_eq,3.0)*(1.0 - X_D);
		else
			dXdt = -8.9588e9 * exp( -1.6262e5/R/T_s )*pow(1.0-p_w_g/p_eq,3.0)*2.0*pow(1.0 - X_D, 0.5);
#endif
	}

	qR = (rho_up - rho_low) * dXdt;

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
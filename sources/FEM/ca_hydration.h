#ifndef _CA_HYDRATION_H
#define _CA_HYDRATION_H

#include "Eigen/Eigen"

class ca_hydration
{
public:
	ca_hydration(double T_solid, 
	             double T_gas,  
				 double p_gas,
				 double w_water, 
				 double rho_s_initial,
				 double phi_S,
				 double delta_t,
				 std::string system);
	~ca_hydration(void);

	void update_param(double T_solid, 
		              double T_gas, 
					  double p_gas,
					  double w_water, 
					  double rho_s_initial,
					  double phi_S,
					  double delta_t,
					  std::string system); 
	
	void calculate_qR();
	double get_qR(); 
	void get_x(Eigen::VectorXd& output_x);
	void set_rho_s(double new_rho_s);
	void eval(double t, Eigen::VectorXd &y, Eigen::VectorXd &dydx);

private: 
	//const double rho_caoh2, rho_cao; // density for Ca(OH)2 and CaO
	const double R;  // [J/mol/K]
	double rho_s;    // solid phase density
	double rho_s_0;  // initial solid phase density
	double phi_solid; //solid volume fraction
	double p_gas;    // gas phase pressure; 
	double p_w_g;    // pressure of H2O on gas phase; 
	double p_eq;     // equilibrium pressure; // unit in bar
	double T_eq;     // equilibrium temperature; 
	double T_s;      // solid phase temperature; 
	double T;        // gas phase temperature; 
	double dXdt;     // rate of (de)hydration; 
	double qR;       // rate of solid density change; 
	double w_h2o;    // mass fraction of water in gas phase; 
	double X_D;      // mass fraction of dehydration (CaO) in the solid phase; 
	double X_H;      // mass fraction of hydration in the solid phase; 
	double dt;       // time step size;
	double rho_low; //lower density limit
	double rho_up; //upper density limit
	double reaction_enthalpy;
	double reaction_entropy;
	std::string reaction_system;

	Eigen::VectorXd x; 
	size_t i;        // iterator

};

#endif


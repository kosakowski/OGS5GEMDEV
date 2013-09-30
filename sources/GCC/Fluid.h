//using namespace std; 
class Fluid
{
private:
public:
	Fluid(void);
	~Fluid(void);


/* Data */
	static double TT, PP;

	//The Properties of Gases and Liquids, Fifth Edition   Bruce E. Poling, John M. Prausnitz, John P. O’Connell
	typedef struct
	{
		std::string name;
		std::string formula;
		double Tc; //critical temperature, K
		double Vc; //critical volume, cm3/mol
		double Pc; //critical pressure, bar
		double w;  //omega, acentric factor (See Chap. 2)
		double u;  //mu, dipole moment, D(debye)
		double k;  //association factor Eq. (9-4.11) 
		double M;  //molecular weight, g/mol
		double x;  //mole fraction
	}compound_properties; // eq 9-5.24 to eq 9-5.44
	static std::vector<compound_properties> Component;



/* Methods */
	static double co2_viscosity (double rho, double T);
	static void initial(void);
	static double viscosity_Chung(double T, std::vector<std::string> component_formula, std::vector<double> mole_amount); //low pressure condition
	static double viscosity_Chung_p(double T, double P, double rho, std::vector<std::string> component_formula, std::vector<double> mole_amount); //high pressure condition
	static double viscosity_Chung_CO2(double T, double P);
	static double viscosity_Chung_H2O(double T, double P);
	static double viscosity_TR(double T, double P, std::vector<std::string> component_formula, std::vector<double> mole_amount);

	static void fluid_properties_calc(double T, double P, double mCO2, double &rho_g, double &rho_l, double &solu_CO2, double &eta_g, double &eta_l);
	static void fluid_properties_calc_CH4(double T, double P, double mCH4, double &rho_g, double &rho_l, double &solu_CH4, double &eta_g, double &eta_l);

	static void fluid_properties_calc_H2(double T, double P, double mH2,  double &rho_g, double &rho_l, double &solu_H2, double &eta_g, double &eta_l);

	static double Gex_RT(double T, double P, double m);

	static double ion_limiting_conductivity(double T, double P, std::string ion);

	static double viscosity_XCl(double T, double P, double m);//to test the method

	static double viscosity_CO2_XCl(double T, double P, double mCO2, double mLiCl, double mNaCl, double mKCl);	



	static double Jones_Dole_viscosity(double T, double P, double mCO2, double mLiCl, double mNaCl, double mKCl);

	static double molarity_to_molality_NaCl(double T, double P, double c);

	static double molefraction_to_molality_NaCl(double x);

	static double molality_to_molarity_CO2(double T, double P, double mNaCl, double mCO2);

	static int viscosity_interface();
	static void entrance(void);
};
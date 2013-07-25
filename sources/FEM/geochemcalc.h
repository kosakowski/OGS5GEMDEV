/*
   #ifndef DLLHEADER_H_INCLUDED
   #define DLLHEADER_H_INCLUDED


   extern __declspec (dllimport) double IF97_density(double, double);
   extern __declspec (dllimport) double IF97_viscosity(double, double);
   extern __declspec (dllimport) double IF97_dielectric(double, double);
   extern __declspec (dllimport) double IF97_Psat(double);

   extern __declspec (dllimport) double VLE_density_CO2(double, double);
   extern __declspec (dllimport) double VLE_density_CH4(double, double);
   extern __declspec (dllimport) double VLE_density_H2O(double, double);
   extern __declspec (dllimport) double VLE_Psat_H2O(double);
   extern __declspec (dllimport) double VLE_solubility_CO2(double, double, double);
   extern __declspec (dllimport) double VLE_solubility_CH4(double, double, double);
   extern __declspec (dllimport) double VLE_fraction_H2O(double, double, double);
   extern __declspec (dllimport) double VLE_solubilityNEW_CO2(double T, double P, double mNaCl);
   extern __declspec (dllimport) double VLE_pressure_CO2(double T, double D);

   extern __declspec (dllimport) double Density_CO2brine(double T, double P, double mNaCl, double mCO2);
   extern __declspec (dllimport) double Density_viscosity(double T, double P, double m, int f);
    //T (K), P (bar), m (mol/kg)  flag 0-LiCl, 1-NaCl, 2-KCl


   extern __declspec (dllimport) double HKF_HKFcalcw(double T, double P, int ghs);
   extern __declspec (dllimport) double HKF_HKFcalc(double T, double P, int ghs, std::string name, int sub, int type);

   extern __declspec (dllimport) int HKF_OGS_loadparam(std::string species_name0, int &type, double &charge, double param[][4]);
   extern __declspec (dllimport) int HKF_OGS_calc(double T, double P, double &G, double &H, double &S, int type, double charge, double param[][4]);

   #endif
 */

double IF97_density(double, double){return 0.0; }
double IF97_viscosity(double, double){return 0.0; }
double IF97_dielectric(double, double){return 0.0; }
double IF97_Psat(double){return 0.0; }

double VLE_density_CO2(double, double){return 0.0; }
double VLE_density_CH4(double, double){return 0.0; }
double VLE_density_H2O(double, double){return 0.0; }
double VLE_Psat_H2O(double){return 0.0; }
double VLE_solubility_CO2(double, double, double){return 0.0; }
double VLE_solubility_CH4(double, double, double){return 0.0; }
double VLE_fraction_H2O(double, double, double){return 0.0; }
double VLE_solubilityNEW_CO2(double T, double P, double mNaCl)
{
	(void)T;
	(void)P;
	(void)mNaCl;
	return 0.0;
}
double VLE_pressure_CO2(double T, double D)
{
	(void)T;
	(void)D;
	return 0.0;
}
double Density_CO2brine(double T, double P, double mNaCl, double mCO2)
{
	(void)T;
	(void)P;
	(void)mNaCl;
	(void)mCO2;
	return 0.0;
}
double Density_viscosity(double T, double P, double m, int f)
{
	(void)T;
	(void)P;
	(void)m;
	(void)f;
	return 0.0;
}
//T (K), P (bar), m (mol/kg)  flag 0-LiCl, 1-NaCl, 2-KCl

double HKF_HKFcalcw(double T, double P, int ghs)
{
	(void)T;
	(void)P;
	(void)ghs;
	return 0.0;
}
double HKF_HKFcalc(double T, double P, int ghs, std::string name, int sub, int type)
{
	(void)T;
	(void)P;
	(void)ghs;
	(void)name;
	(void)sub;
	(void)type;
	return 0.0;
}

int HKF_OGS_loadparam(std::string species_name0, int &type, double &charge, double param[][4])
{
	(void)species_name0;
	(void)type;
	(void)charge;
	(void)param;
	return 0;
}
int HKF_OGS_calc(double T,
                 double P,
                 double &G,
                 double &H,
                 double &S,
                 int type,
                 double charge,
                 double param[][4])
{
	(void)T;
	(void)P;
	(void)G;
	(void)H;
	(void)S;
	(void)type;
	(void)charge;
	(void)param;
	return 0;
}

/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

//-------------------------------------
// rf_REACT_GEM.h
// Haibing Shao 23.05.07
// haibing.shao@ufz.de
// GEM Reaction Package
// based on the PSI node-GEM source code
// using the node-GEM code from Paul Sherrer Institute (PSI)
//-------------------------------------
#ifndef RF_REACT_GEM_H
#define RF_REACT_GEM_H

#ifdef GEM_REACT

#include <time.h>
#include <cmath>
#include <string>
#include <boost/thread.hpp>
//#include <thread.hpp>
#include "rf_pcs.h"
// #include "rfmat_cp.h"
#include "GEM/node.h"
#include "rf_mfp_new.h"

#if defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#include "petscksp.h"
#include "petscmat.h"
typedef Mat PETSc_Mat;
typedef Vec PETSc_Vec;
#endif
/**
 * class REACT_GEM for coupling OGS with GEMS
 */
class REACT_GEM
{
private:
    string ipm_input_file_path;
    string dbr_input_file_path;
    string dbr_bc_input_file_path;
    string dch_input_file_path;
    string init_input_file_path;
    string init_input_file_sig;

public:
    REACT_GEM(void);
    ~REACT_GEM(void);	
    
    //FEM
	CFiniteElementStd* Fem_Ele_Std;

    /// Instance of TNode class
    /// HS: 06.2007 Only set one TNode here, repeatedly use its resources for GEM calculation
//		TNode* m_Node;

    // DATABR structure for exchange with GEMIPM
//      DATACH* dCH;                                //pointer to DATACH
//      DATABR* dBR;                                //pointer to DATABR

    // Read function for gem input file
    ios::pos_type Read(ifstream *gem_file);

    // Number of ICs, DCs, Phases and Phases-solutions kept in the node or elem;
    long nIC, nDC, nPH, nPS;

    // Data structure for each node to carry the chemical information (real FMT problems consider many nodes)
    // Names are consistent with the DataBridge structure (also see "\GEM\databr.h")
    long *m_NodeHandle, *m_NodeStatusCH, *m_IterDone, *m_IterDoneCumulative, *m_IterDoneIndex;
    // this is for porosity calculated on volume of solids
    double *m_porosity;
    /// this is used for kinetic law number 4 and other lasws depending on initial amounts of stuff
    double  *m_porosity_initial, *m_volumes_initial;

    /// this we need for porosity coupling to groundwater flow & multiphase flow (_pre: old solution, _pi: previous picard iteration)
    double *m_fluid_volume, *m_gas_volume, *m_fluid_volume_pts, *m_gas_volume_pts,  *m_fluid_volume_pi, *m_gas_volume_pi;

    /// take the fluid density from GEMS for density driven flow
    double *m_fluid_density;
    
    /// indexes, which one in the xDC vector is water, oxygen or hydrogen
    int idx_water, idx_hydrogen, idx_oxygen, idx_carbon, idx_co2g;

    double *m_T, *m_P, *m_Vs, *m_Ms,
    *m_Gs, *m_Hs, *m_IC, *m_pH, *m_pe,
    *m_Eh;

    double *m_xDC, *m_gam, *m_xPH, *m_aPH, *m_vPS, *m_mPS, *m_bPS,
    *m_xPA, *m_dul, *m_dll, *m_bIC, *m_bIC_dummy, *m_rMB, *m_uIC, *m_bSP, *m_phase_volume;

    double  *m_porosity_Elem;

    /// data for transport of IC
    double *m_soluteB, *m_soluteB_corr;
    /// data for charge
    double *m_chargeB, *m_chargeB_pre;
    
    // previous time step DC values (_pre: old solution, _pi: previous picard iteration)
    double *m_xDC_pts, *m_xDC_pi;                          // previous time step Concentration;
    double *m_xDC_MT_delta;                     // delta C from Mass Transport;
    double *m_xDC_Chem_delta, *m_bIC_Chem_delta;                   // delta C from Chemistry;
    // previous time step IC values (_pre: old solution, _pi: previous picard iteration)
    double *m_soluteB_pts, *m_bIC_pts, *m_soluteB_pi, *m_bIC_pi;
    // previous time step kinetic values (_pre: old solution, _pi: previous picard iteration)
    double *dmdt_pts, *m_dul_pts, *m_dll_pts, *dmdt_pi, *m_dul_pi, *m_dll_pi;                              // old kinetically controlled rates
    double *mol_phase_pts, *omega_phase_pts, *omega_components_pts, *mol_phase_pi, *omega_phase_pi, *omega_components_pi; // old saturation indices and mole amounts in phase
    
    
    double *m_excess_water;                     // excess water in m3/s for each node;
    double *m_excess_gas;                       // excess gas in m3/s for each node;
    double *m_saturation;
    double *m_Node_Volume;                      // Volume around the node;

    // this we need for kinetics
    double *mol_phase, *mol_phase_initial, *omega_phase,*omega_components, *rate_components;

    double *dmdt;                               // kinetically controlled rates
    int CalcLimits ( long in, double deltat, TNode* m_Node);
    int CalcLimitsInitial ( long in, TNode* m_Node);
    int SetConstraints ( long in, TNode* m_Node);
    int SetLimitsInitial ( long in, TNode* m_Node);    
    int CalcLimitsSolidSolution ( long in, long ii,int flag_equilibration ,TNode* m_Node);
    double CementHydrationKinetics (double oldvalue,  long kin_phasenr, double dt ); // calculates cement hydration kinetics according ot Lothenbach and Wieland (2008=
    long *m_boundary;                            //holds marker for boundary nodes
    double max_kinetic_timestep;               // variable used for limiting time step
    double MaxChangeKineticsPhase(long in, int ii, int i,TNode *m_Node); // maximum allowed change to 10% concentrations
    bool calc_limits;                          // do this for first Picard iteration
    int CheckConstraints ( long in,TNode* m_Node);
    double MaxDtKinetics (long int,  TNode* m_Node); // function to calculate max timestep in case of kinetics for a node 
    
    int SolveChemistry(long in, TNode* m_Node); // extra time loop for chemistry & kinetics
    
    // CRFProcess *m_pcs;                          // pointer to the PCS Class.
    CRFProcess *m_flow_pcs;                     // pointer to the flow PCS.


    /** Initialization of the GEM TNode Class
    *  return: 0-ok;
    *          1-loading init file failure;
    *          3-dch problem;
    *           4-dbr problem;
            */
    short Init_Nodes(string Project_path);      // Initialization of the GEM TNode Class
    //  return: 0-ok;
    //          1-loading init file failure;
    //          3-dch problem;
    //          4-dbr problem;

    short Init_RUN(string Project_path);                           // Run the node-GEM
    //  return: 0-ok;5-GEM does not converge
    short Run_MainLoop();

    //  return: 0-ok;5-GEM does not converge

    string Get_Init_File_Path(void);
    string Get_DCH_File_Path(void);
    string Get_IPM_File_Path(void);
    string Get_DBR_File_Path(void);

    int Set_Init_File_Path(string m_path);
    int Set_DCH_FILE_PATH(string m_path);
    int Set_IPM_FILE_PATH(string m_path);
    int Set_DBR_FILE_PATH(string m_path);

    bool Load_Init_File(string m_Project_path, TNode* m_Node );
    long* mp_nodeTypes;

    //---flags------
    int initialized_flag;                       //0 - not initialized 1 - initialized
    int flag_iterative_scheme;                  // 0 - sequential non-iterative scheme;
    // 1 - standard iterative scheme;
    int max_gems_iteration_loop;  // maximum number of iterations
    double  iteration_eps;  // leave iteration if this is achieved
    bool flag_concentrations_with_water; // 1 use water itself also in concentrations for H and O ...0: water is substracted from O and H concentrations (only dissolved species are transported)
    bool flag_accept_bad_gem; /// accept bad GEM solutions 
    bool flag_gem_sia; /// allow warmstart option to speedup calculatuons (better convergence of GEM kernel)
    
    int flag_gas_diffusion;                      // simple implemention of gas diffusion for benchmark
    double *m_co2, *m_co2_transport;    
    short GetCO2Value_MT ( long node_Index, int timelevel, double* m_co);
    short SetCO2Value_MT ( long node_Index, int timelevel,  double co2value);
    int heatflag;                               //0-initialized and not heat transport;1-heat_transport;
    int flowflag;                               //0-initialized;1-GROUNDWATER_FLOW;2-LIQUID_FLOW;3-RICHARDS_FLOW;4-FLOW;
    int flag_porosity_change;                   //0-porosity change not coupled into transport; 1=coupled;
    int flag_coupling_hydrology;                //0-without coupling; 1=with coupling;
    int flag_calculate_boundary_nodes;          // set to zero to avoid precipitation/dissolution (porosity change) at boundary nodes
    int gem_pressure_flag;                      //shall we give a constant user defined pressure to gems?
    int flag_transport_b;                       //1: transport only dissolved components of b vector; 0: transport full speciation
    int flag_hayekit;                           // flag to force behaviour for hayekit benchmark...fluid phase is scaled by scaling water only and not solutes
    long m_max_failed_nodes; ///maximum number of failed nodes
    int flag_disable_gems;             //disable gems calculations in main loop ..not for initialization!
    bool *m_calculate_gems;     // size: nNodes: =1 (default) =0 (switch off calculation for specific nodes)
    //--------------

    long nNodes;                                // number of all nodes;
    long nElems;                                // number of all elements;
    int GetHeatFlag_MT(void);
    int GetFlowType_MT(void);
    long GetNodeNumber_MT(void);
    long GetElemNumber_MT(void);
    void GetFluidProperty_MT(void);

    short GetInitialReactInfoFromMassTransport(int timelevel);
    short GetReactInfoFromMassTransport(int timelevel);
    short SetReactInfoBackMassTransport(int timelevel);
    short SetReactInfoBackMassTransportPicardIteration(int timelevel);
    void GetReactInfoFromGEM(long in,  TNode* m_Node );
    void SetReactInfoBackGEM(long in,  TNode* m_Node );
    // necessary for reload with gems
    int WriteReloadGem();
    int ReadReloadGem();

    double GetTempValue_MT(long node_Index, int timelevel);
    double GetPressureValue_MT(long node_Index, int timelevel);
    short GetDCValue_MT(long node_Index, int timelevel, double* m_DC);
    short GetBValue_MT ( long node_i, int timelevel, double *m_soluteB);
    short GetSoComponentValue_MT(long node_Index, int timelevel, double* m_Phase,  TNode* m_Node );
    double GetDCValueSpecies_MT ( long node_Index, int timelevel, int iDc );
    short SetTempValue_MT(long node_Index, int timelevel, double temp);
    short SetPressureValue_MT(long node_Index, int timelevel, double pressure);
//    short SetDCValue_MT(long node_Index, int timelevel, double* m_DC);
    short SetBValue_MT(long node_Index, int timelevel, double* m_soluteB); //, double* m_chargeB);

    int IsThisPointBCIfYesStoreValue ( long index, CRFProcess* m_pcs, double& value );/// taken from rf_REACT_BRNS

    /// repeated_fail used to indicated failed nodes duting Picard iteration..bad style to use it like this..but most convenient 
    double repeated_fail;
    /// Copy current values into previous time step values
    void CopyCurXDCPre ( void );
    void UpdateXDCChemDelta ( void );
    void UpdatebICChemDelta ( void );
    void ResetbICChemDelta ( void );
    void CopyCurBPre ( void );
    double CalcSoluteBDelta ();
    double CalcMaxSoluteBDelta ();
    double CalcSoluteBDeltaNode ( long in );
    double minimum_kinetic_time_step; // for refinement of kinetics
    void StoredSolutionAllPicardIteration ( );
    void StoreOldSolutionAll ( );
    void RestoreOldSolutionAll ( );
    void RestoreSolutionAllPicardIteration ( );
    void StoreOldConcentrationsAll();
    void RestoreOldSolutionNode ( long in );
    short RestoreOldConcMasstransport_MT ( long node_Index);
    void CopyCurKineticsPre ( void );
    // struff related to charge
    double CalculateCharge (long in); //for a given node number in
    void CalculateChargeFromGEMS ( long in, TNode* t_Node);
    /// this is only for porosity interpolation to elemens
    void ConvPorosityNodeValue2Elem(int i_timestep);
    int CalcPorosity(long in,  TNode* m_Node );

    double min_possible_porosity, max_possible_porosity;
    void ScaleVolume_Water(long in);

    // Set porosity in Mass Transport
    int SetPorosityValue_MT(long ele_Index, double m_porosity_Elem, int i_timestep);
    int SetSourceSink_MT(long in, double time_step_size /*in sec*/);

    // pass fluid density back
    double FluidDensity(long elem, int gaussnode, CFiniteElementStd* fem); // if calling routine asks for density at element
    
    // find which one in xDC vector is water
    int FindWater_xDC(TNode* m_Node);
    int Findhydrogen_bIC ( TNode* m_Node);
    int Findoxygen_bIC ( TNode* m_Node);
    int Findcarbon_bIC ( TNode* m_Node);
    int Findco2g_bIC ( TNode* m_Node);

    //kg44 11/2008 for kinetics
    int CalcReactionRate ( long node,  TNode* m_Node );
    double SurfaceAreaPh ( long kin_phasenr,long in ,  TNode* m_Node , double surf_area);

    // concentration related
    int ConcentrationToMass (long l /*idx of node*/,int i_timestep );
    int MassToConcentration (long l /*idx of node*/,int i_timestep,  TNode* m_Node );

    // Unit conversion for pressures
    double Pressure_Pa_2_Bar(double Pre_in_Pa);
    double Pressure_Bar_2_Pa(double Pre_in_Bar);
    double Pressure_M_2_Bar(double Pre_in_M, double flu_density );
    double Pressure_Bar_2_M(double Pre_in_Bar, double flu_density );
    double Pressure_M_2_Pa ( double Pre_in_M, double flu_density );
    double Pressure_Pa_2_M ( double Pre_in_Pa, double flu_density );
    // Calculate the volume of the nodes;
    // given argument is the index of one particular node;
    double GetNodeAdjacentVolume(long Idx_Node);

    // GEMS mass scaling parameter
    double gem_mass_scale;
    // GEM temperature (without coupling to temperature)
    double m_gem_temperature;
    // GEM pressure (needed for Richards flow)
    double m_gem_pressure;

// the next two are always defined, such that it also works in serial version    
    int myrank;
    int mysize;

    // MPI implementation
#if defined(USE_MPI_GEMS) 
    void CleanMPIBuffer(void);
    void CopyToMPIBuffer(long in);
    void GetGEMResult_MPI(void);
#endif


    double GetNodePorosityValue( long node_Index);
    double GetNodePorosityValueInitial( long node_Index);    
    double GetNodeFluidDensityValue( long node_Index);

    // Name lists from DCH file!
    // const long int
    //  MaxICN =      6,      // IC name length
    //  MaxDCN =      16,     // DC name length
    //  MaxPHN =      16;     // PH name length
    char (*m_ICNL)[MaxICN]; // List of IC names in the system, [nIC]  of MaxICN length
    char (*m_DCNL)[MaxDCN]; // List of DC names in the system, [nDC] of MaxDCN length
    char (*m_PHNL)[MaxPHN]; // List of Phase names  [nPH]  of MaxPHN length
#if defined(USE_PETSC)    
    PetscScalar *gem_glob_buff, *gem_glob_x1 , *gem_glob_x0;
    void WriteVTKGEMValuesPETSC(PetscViewer viewer);
    // for synchronizing data
    void SynchronizeData(PetscScalar *data);
    long GetGlobalNodeNumber_MT ( void); 
    long GetLocalNodeNumber_MT(void);
    long loc_NodesNumber_Linear, NodesNumber_Linear,glob_NodesNumber_Linear;
#endif
    void WriteVTKGEMValues(fstream &vtk_file);    
    // timer 
    double GetTimeOfDayDouble();
    
    typedef struct // definition for conditions to switch off GEMS
    {
      int condition_type; //kind of model how to process GEMS nodes initially
      int nB; //number of transport process we would like to control (starting at zero)
      double l_amount; // which lower amount 
      double u_amount; // which upper amount
    }Calculate_GEMS;
    vector<Calculate_GEMS> m_calc;
    
    
    int flag_loose_kinetics; // flag used in calclimits to give more freedom to dul/dll passed to gems
    
    typedef struct
    {
        //kg44 25.11.2008 kinetics...for coupling with GEMS
        //
        string phase_name;
        int phase_number;
        int dc_counter;
        int kinetic_model;                       // only 1 = GEMS implemented right now
        int n_activities;                        // number of species for activities
        string active_species[10];               // name for species ...maximum 10 names
        /**	this vector holds the kinetic material parameters
        *      0,1,2  double E_acid,E_neutral,E_base; // activation energies
        *      3-5  double k_acid, k_neutral,k_base; // dissolution/precipitation rate constants
        *      6-11  double p1,q1,p2,q2,p2,q2; // exponents for omega
        *      12,13, 14  double n_1, n_2, n_3; // exponents for acidic neutral and base cases for species one
        *      append for each species another set of n_1, n_2, n_3 (up to 10 sets -> up to ten species)
        **/
        double kinetic_parameters[41];
        int surface_model;                       // currently only 1 implemented
        double surface_area[10];
        int ss_endmembers; // special model for solid solutions...only read for kinetic model == 5
        double *ss_scaling; // special model for solid solutions...only read for kinetic model == 5
        string endmember_species[10]; // species for cation exchange with endmembers( maximum  10 endmembers) same order as endmembers
    } Kinetic_GEMS;
    vector<Kinetic_GEMS> m_kin;

    typedef struct
    {
        /** for kinetically controlled phases we can set initial constraints
	 *  first model: the initial amount for a component is controled by one elemnt of b vector 
	 * condition applies of b vector is between l_amount and u_amount
	 * 
	 **/
	string component_name; // must be a dependent component name
	int n_comp; // corresponding number in component list!
	int condition_type;
	int nB; //number of b element we use for control (starting at zero) 
	double l_amount; // which lower amount 
        double u_amount; // which upper amount
	double ll_constraint;
	double ul_constraint;
    } Constraints_GEMS; 
    vector<Constraints_GEMS> m_constraints;

    
    // here we define the variables we need for the threads
    unsigned int gem_nThread, gem_nbar; // number of threads, number of threads + 1 (number of threads + master that have to cross a barrier)
    boost::barrier* gem_barrier_start; // start barrier for calculations
    boost::barrier* gem_barrier_finish; //stop barrier when calculations are finished
    boost::thread* gemThread; // the gems worker threads
    boost::mutex rwmutex; // used to lock during write or read operations..example gem init or write dbrs or cout!
    boost::mutex getnode_mutex; // used to lock during getnodeindex
    void gems_worker(int tid, string tinit_path); 
};


#define GEM_FILE_EXTENSION ".gem"
/** This is the function for reading the OGS-GEM specific parameters.
*  Todo: if .gems not exist, decouple gems
*/
extern bool GEMRead(string base_file_name, REACT_GEM *m_GEM_p);
#endif
#endif

// For compilers that don't support precompilation, include "wx/wx.h"
#include "wx/wxprec.h"
#include "wx/filedlg.h"
#include "wx/textctrl.h"
#include "wx/textfile.h"
//#include "wxwin32x32.xpm"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

class BadConversion : public std::runtime_error
{
	public:
		BadConversion ( const std::string& s )
				: std::runtime_error ( s )
		{ }
};

inline std::string stringify ( int x )
{
	std::ostringstream o;
	if ( ! ( o << x ) )
		throw BadConversion ( "stringify(double)" );
	return o.str();
}


#include "gems2ogs.h"

IMPLEMENT_APP ( HelloWorldApp )

// This is executed upon startup, like 'main()' in non-wxWidgets programs.
bool HelloWorldApp::OnInit()
{
	// wxFrame *frame = new wxFrame((wxFrame*) NULL, -1, _T("Hello wxWidgets World"));
	// frame->CreateStatusBar();
	// frame->SetStatusText(_T("Hello World"));

	// Create the main application window
	MyFrame *frame = new MyFrame ( wxT ( "GEM 2 OpenGeoSys Convertor" ) );

	// Show it
	frame->Show ( true );

	// set the GEM initializatio flag to false.
	frame->flag_status = 0;

	frame->water_volume = 0.0;
	frame->idx_water    = -1 ;
	frame->nDC          = 0  ;

	// Start the event loop
	return true;
}

BEGIN_EVENT_TABLE ( MyFrame, wxFrame )
	EVT_MENU ( wxID_ABOUT,   MyFrame::OnAbout )
	EVT_MENU ( wxID_EXIT,    MyFrame::OnQuit )
	EVT_MENU ( wxID_OPEN,    MyFrame::OnOpenDCH )
	EVT_MENU ( ID_Con_IC,    MyFrame::OnConIC )
	EVT_MENU ( ID_Con_BC,    MyFrame::OnConBC )
	EVT_MENU ( ID_Con_MCP,   MyFrame::OnConMCP )
	EVT_MENU ( ID_Con_OUT,   MyFrame::OnConOUT )
	EVT_MENU ( ID_Con_ICB,    MyFrame::OnConICB )
	EVT_MENU ( ID_Con_BCB,    MyFrame::OnConBCB )
	EVT_MENU ( ID_Con_BCBW,    MyFrame::OnConBCBW )
	EVT_MENU ( ID_Con_MCPB,   MyFrame::OnConMCPB )
	EVT_MENU ( ID_Con_OUTB,   MyFrame::OnConOUTB )
	EVT_MENU ( ID_Con_DLUB,   MyFrame::OnConDLUB )

	EVT_MENU ( ID_Con_Set,   MyFrame::OnConSet )
	EVT_MENU ( wxID_SAVE,    MyFrame::OnSave )
	EVT_MENU ( wxID_SAVEAS,  MyFrame::OnSaveAs )
END_EVENT_TABLE()

void MyFrame::OnAbout ( wxCommandEvent& event )
{
	wxString msg;
	msg.Printf ( wxT ( "GEM to GeoSys Convertor\nAuthor: haibing.shao@ufz.de and georg.kosakowski@psi.ch" ) );

	wxMessageBox ( msg, wxT ( "About" ), wxOK | wxICON_INFORMATION, this );
}

void MyFrame::OnQuit ( wxCommandEvent& event )
{
	// Destroy the TNode data structure
	delete m_Node;

	// Destroy the frame
	Close();
}

void MyFrame::OnOpenDCH ( wxCommandEvent &event )
{
	int k,idummy;

	wxString caption = wxT ( "Choose a GEM initialization file" );
	wxString wildcard =
	    wxT ( "GEM initialization files (*.lst)|*.lst|Any files (*.*)|*.*" );
	wxString defaultDir = wxT ( "" );
	wxString defaultFilename = wxEmptyString;

	wxFileDialog dialog ( this, caption, defaultDir, defaultFilename, wildcard, wxFD_OPEN );

	if ( dialog.ShowModal() == wxID_OK )
	{
		wxString path = dialog.GetPath();
                char *dbrfiles_lst_name = 0;

		// Now initialize the TNode structure with proper input file.
		long* mp_nodeTypes;
		mp_nodeTypes = new long;
		* ( mp_nodeTypes ) = 0;

		if ( !m_Node )
		{
		  //	m_Node->~TNode();
			m_Node = new TNode(); 
		}

		// Initializing GEM
		//m_Node->GEM_init ( path.mb_str() , dbrfiles_lst_name, mp_nodeTypes , false );
		m_Node->GEM_init (path.mb_str());

		// Getting direct access to DataCH structure in GEMIPM2K memory
		dCH = m_Node->pCSD();

		// Getting direct access to work node DATABR structure which
		// exchanges data between GEMIPM and FMT parts
		dBR = m_Node->pCNode();

		nDC = dCH->nDC ;

		dBR->NodeStatusCH = NEED_GEM_AIA;

		// run GEM
		idummy = m_Node->GEM_run ( false );
		if ( ! ( idummy == OK_GEM_AIA || idummy == OK_GEM_SIA ) )
		{
			dBR->NodeStatusCH = NEED_GEM_AIA;
			idummy = m_Node->GEM_run (false );
			if ( ! ( idummy == OK_GEM_AIA || idummy == OK_GEM_SIA ) )
			{
				// problem initializing GEM
				wxString msg;
				msg.Printf ( wxT ( "Problem occur while initializing GEM input files!!!" ) );
				wxMessageBox ( msg, wxT ( "Error" ), wxOK | wxICON_INFORMATION, this );
				textCtrl->AppendText ( msg );
			}

		}
		// now gems is ok and we can start

		water_volume=0.0;
                water_volume = m_Node->Ph_Volume ( 0 ); // liquid is always first phase
		water_volume /= dBR->Vs; // now it is relative to 1 m^3 which is equivalent to water filled porosity

		if ( ( water_volume <= 0.0 ) || ( water_volume >1.000001 ) ) // allow for some numerical inaccuracy in volume ...divison by Vs might be not completely accurate
		{
			// problem initializing GEM
			wxString msg;
			msg.Printf ( wxT ( "Problem with water volume! volume= %e"),water_volume  );
			wxMessageBox ( msg, wxT ( "Error" ), wxOK | wxICON_INFORMATION, this );
			textCtrl->AppendText ( msg );
		}
		else // GEM is initialized
		{
			// message handling
			wxString msg;
			msg.Printf ( wxT ( "GEM input files read successfully" ) );
			wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );

			// switch on the initialized flag
			flag_status = 1;

			// show the information of loaded GEM structure in the text control
			msg.Printf ( wxT ( "-----------------------------------\n" ) );
			textCtrl->AppendText ( msg );
			msg.Printf ( wxT ( "---Infomation of the GEM project---\n" ) );
			textCtrl->AppendText ( msg );
			msg.Printf ( wxT ( "Number of ICs:      %d\n" ), dCH->nIC );
			textCtrl->AppendText ( msg );
			msg.Printf ( wxT ( "Number of DCs:      %d\n" ), dCH->nDC );
			textCtrl->AppendText ( msg );
			msg.Printf ( wxT ( "Number of Phases:  %d\n" ),  dCH->nPH );
			textCtrl->AppendText ( msg );
			msg.Printf ( wxT ( "-----------------------------------\n" ) );
			textCtrl->AppendText ( msg );
			m_Node->GEM_write_dbr ( "dbr_for_debug.txt" );
			m_Node->GEM_print_ipm ( "ipm_for_debug.txt" );
		}
	}
}

void MyFrame::OnConIC ( wxCommandEvent &event )
{
	int i;
	double concentration;
	concentration = 0.0;

	if ( flag_status ) // if initialized
	{
		textCtrl->Clear();
		string temp_name,strnum; // used to get the current name of DC;

		for ( i=0 ; i < dCH->nDC ; i++ )
		{
			textCtrl->AppendText ( IC_Head );
			textCtrl->AppendText ( IC_PCS_Head );
			textCtrl->AppendText ( IC_PCS );
			textCtrl->AppendText ( IC_Name_Head );

			temp_name=stringify ( i+1 );
			temp_name.append ( "-" );
			strnum = dCH->DCNL[i]; // get component name;
			temp_name.append ( strnum );
			IC_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );

			textCtrl->AppendText ( IC_Name );
			textCtrl->AppendText ( IC_Geo_Head );
			textCtrl->AppendText ( IC_Geo );
			textCtrl->AppendText ( IC_Value_Head );
			textCtrl->AppendText ( IC_Value_Typ );

			if ( ( dCH->ccDC[i] == 'S' ) || ( dCH->ccDC[i] == 'E' ) || ( dCH->ccDC[i] == 'T' ) )
			{
				// for disolved species do not convert to concentrations...This will be done later, after calling GEMS Init. Otherwise we have to pass water_volume for each node to Geosys
				//			    concentration = dBR->xDC[i] / (water_volume) ;
				concentration = dBR->xDC[i] ;
				//			    concentration = dBR->xDC[i] / (water_volume) ;
				IC_Value.Printf ( wxT ( "%1.14e\n" ), concentration / dBR->Vs ); // scale to a volume of 1m^3
			}
			else
			{
				// for other species, keep it in moles and normalize volume of reactive subsystem
				IC_Value.Printf ( wxT ( "%1.14e\n" ), ( dBR->xDC[i] / dBR->Vs ) ); //scale to a volume of 1 m^3
			}

			textCtrl->AppendText ( IC_Value );
			textCtrl->AppendText ( IC_End );
		}

		flag_status = 2;
	}
	else// if not initialized
	{
		wxString msg;
		msg.Printf ( wxT ( "GEM not initialized. Load GEM data first!" ) );
		wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );
	}
}

void MyFrame::OnConMCP ( wxCommandEvent &event )
{
	int i;

	if ( flag_status ) // if initialized
	{
		textCtrl->Clear();
		string temp_name,strnum; // used to get the current name of DC;

		for ( i=0 ; i < dCH->nDC ; i++ )
		{

			MCP_Head.Printf ( wxT ( "#COMPONENT_PROPERTIES ; comp%i \n" ), i+1 );
			textCtrl->AppendText ( MCP_Head );
			textCtrl->AppendText ( MCP_Name_Head );


			temp_name=stringify ( i+1 );
			temp_name.append ( "-" );
			strnum = dCH->DCNL[i]; // get component name;
			temp_name.append ( strnum );
			MCP_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );

			textCtrl->AppendText ( MCP_Name );
			textCtrl->AppendText ( MCP_Mobile_Head );

			if ( ( dCH->ccDC[i] == 'S' ) || ( dCH->ccDC[i] == 'E' ) || ( dCH->ccDC[i] == 'T' ) )
			{
				// for disolved species, mobile and with an diffusion number
				textCtrl->AppendText ( MCP_Mobile_Yes );
				textCtrl->AppendText ( MCP_Dif_Head );
				textCtrl->AppendText ( MCP_Dif_Yes );
			}
			else
			{
				// for other species, keep it inmobile and without diffusion
				textCtrl->AppendText ( MCP_Mobile_No );
				textCtrl->AppendText ( MCP_Dif_Head );
				textCtrl->AppendText ( MCP_Dif_No );
			}


			textCtrl->AppendText ( BC_End );
		}
// add pH, pe, Eh and NodePorosity as default
		MCP_Head.Printf ( wxT ( "#COMPONENT_PROPERTIES ; comp%i \n" ), 1+ dCH->nDC );
		textCtrl->AppendText ( MCP_Head );
		textCtrl->AppendText ( MCP_Name_Head );
		temp_name = "pH"; // get component name;
		MCP_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );
		textCtrl->AppendText ( MCP_Name );
		textCtrl->AppendText ( MCP_Mobile_Head );
		// for other species, keep it inmobile and without diffusion
		textCtrl->AppendText ( MCP_Mobile_No );
		textCtrl->AppendText ( MCP_Dif_Head );
		textCtrl->AppendText ( MCP_Dif_No );
		MCP_Head.Printf ( wxT ( "#COMPONENT_PROPERTIES ; comp%i \n" ),2+ dCH->nDC );
		textCtrl->AppendText ( MCP_Head );
		textCtrl->AppendText ( MCP_Name_Head );
		temp_name = "pe"; // get component name;
		MCP_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );
		textCtrl->AppendText ( MCP_Name );
		textCtrl->AppendText ( MCP_Mobile_Head );
		// for other species, keep it inmobile and without diffusion
		textCtrl->AppendText ( MCP_Mobile_No );
		textCtrl->AppendText ( MCP_Dif_Head );
		textCtrl->AppendText ( MCP_Dif_No );
		MCP_Head.Printf ( wxT ( "#COMPONENT_PROPERTIES ; comp%i \n" ), 3+ dCH->nDC );
		textCtrl->AppendText ( MCP_Head );
		textCtrl->AppendText ( MCP_Name_Head );
		temp_name = "Eh"; // get component name;
		MCP_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );
		textCtrl->AppendText ( MCP_Name );
		textCtrl->AppendText ( MCP_Mobile_Head );
		// for other species, keep it inmobile and without diffusion
		textCtrl->AppendText ( MCP_Mobile_No );
		textCtrl->AppendText ( MCP_Dif_Head );
		textCtrl->AppendText ( MCP_Dif_No );
		MCP_Head.Printf ( wxT ( "#COMPONENT_PROPERTIES ; comp%i \n" ),4+  dCH->nDC );
		textCtrl->AppendText ( MCP_Head );
		textCtrl->AppendText ( MCP_Name_Head );
		temp_name = "NodePorosity"; // get component name;
		MCP_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );
		textCtrl->AppendText ( MCP_Name );
		textCtrl->AppendText ( MCP_Mobile_Head );
		// for other species, keep it inmobile and without diffusion
		textCtrl->AppendText ( MCP_Mobile_No );
		textCtrl->AppendText ( MCP_Dif_Head );
		textCtrl->AppendText ( MCP_Dif_No );

		// change flag to MCP converted
		flag_status = 4;
	}
	else// if not initialized
	{
		wxString msg;
		msg.Printf ( wxT ( "GEM not initialized. Load GEM data first!" ) );
		wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );
	}


}

void MyFrame::OnConBC ( wxCommandEvent &event )
{
	int i;
	double concentration;
	concentration = 0.0;

	if ( flag_status ) // if initialized
	{
		textCtrl->Clear();
		string temp_name,strnum; // used to get the current name of DC;

		for ( i=0 ; i < dCH->nDC ; i++ )
		{
			textCtrl->AppendText ( BC_Head );
			textCtrl->AppendText ( BC_PCS_Head );
			textCtrl->AppendText ( BC_PCS );
			textCtrl->AppendText ( BC_Name_Head );

			temp_name=stringify ( i+1 );
			temp_name.append ( "-" );
			strnum = dCH->DCNL[i]; // get component name;
			temp_name.append ( strnum );
			BC_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );

			textCtrl->AppendText ( BC_Name );
			textCtrl->AppendText ( BC_Geo_Head );
			textCtrl->AppendText ( BC_Geo );
			textCtrl->AppendText ( BC_Value_Head );
			textCtrl->AppendText ( BC_Value_Typ );

			if ( ( dCH->ccDC[i] == 'S' ) || ( dCH->ccDC[i] == 'E' ) || ( dCH->ccDC[i] == 'T' ) )
			{
				// for disolved species, convert to concentration values in mol/L water
				// for BCs we need concentrations, as BCs are read at the beginning at overwritten each timestep
				// in the current implementation of geosys-gems the boundary conditions are not changed (expect for changes defined in the geosys imput files)
				// small error in porosity during init_gems routine should not be a problem!
				concentration = dBR->xDC[i] / dBR->Vs; //scale to a volume of 1 m^3
				concentration /= ( water_volume ) ; // now it is concentration
				BC_Value.Printf ( wxT ( "%1.14e\n" ), concentration ); // no further scaling with system volume required
			}
			else
			{
				// for other species, keep it in moles and normalize volume of reactive subsystem
				BC_Value.Printf ( wxT ( "%1.14e\n" ), ( dBR->xDC[i] / dBR->Vs ) ); //scale to a volume of 1 m^3
			}
			textCtrl->AppendText ( BC_Value );
			textCtrl->AppendText ( BC_End );
		}

		// change flag to BC converted
		flag_status = 3;
	}
	else// if not initialized
	{
		wxString msg;
		msg.Printf ( wxT ( "GEM not initialized. Load GEM data first!" ) );
		wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );
	}
}

void MyFrame::OnConOUT ( wxCommandEvent &event )
{
	int i;

	if ( flag_status ) // if initialized
	{
		textCtrl->Clear();
		string temp_name,strnum; // used to get the current name of DC;

		textCtrl->AppendText ( OUT_Head );

		for ( i=0 ; i < dCH->nDC ; i++ )
		{

			temp_name=stringify ( i+1 );
			temp_name.append ( "-" );
			strnum = dCH->DCNL[i]; // get component name;
			temp_name.append ( strnum );
			OUT_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );

			textCtrl->AppendText ( OUT_Name );

		}

		textCtrl->AppendText ( OUT_Tail );

		// change flag to MCP converted
		flag_status = 5;
	}
	else// if not initialized
	{
		wxString msg;
		msg.Printf ( wxT ( "GEM not initialized. Load GEM data first!" ) );
		wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );
	}


}

void MyFrame::OnConDLUB ( wxCommandEvent &event )
{
  //example output for constraint
  //
  // $CONSTRAINT_GEMS     ;Aluminate
  //  Periclase 78.4 78.45 1 0 100.0 260.0 ; 3.0e4 3.1e4           ; dll, dul, 1, number of phase for criteria, lower amount, upper amount of phase
  //	
        int i;
	if ( flag_status ) // if initialized
	{
		textCtrl->Clear();
		string temp_name,strnum; // used to get the current name of DC;

		textCtrl->AppendText ( DLU_Head );

		for ( i=0 ; i < dCH->nDC ; i++ )
		{
		  if ( (dBR->dll[i] > 0.0 )|| ((dBR->dul[i]-dBR->dll[i]) <=1e5) ) // do this if obviously some constraints are set
		    {    
                      temp_name=" ";
		      //			temp_name=stringify ( i+1 );
		      //			temp_name.append ( "-" );
			strnum = dCH->DCNL[i]; // get component name;
			temp_name.append ( strnum );
			DLU_Name.Printf ( wxT ( "  %S " ), temp_name.c_str() );
			//IC_Value.Printf ( wxT ( "%1.14e\n" ), ( dBR->bIC[i] / dBR->Vs ) ); //scale to a volume of 1 m^3
			textCtrl->AppendText ( DLU_Name );
                        DLU_Number.Printf (wxT (" %1.14e %1.14e 1 0 %1.14e %1.14e" ),(dBR->dll[i]/dBR->Vs),(dBR->dul[i]/dBR->Vs), ( 0.99 * dBR->bIC[0] / dBR->Vs ), ( 1.01 * dBR->bIC[0] / dBR->Vs )) ;
 		        textCtrl->AppendText ( DLU_Number );
   		        textCtrl->AppendText ( DLU_Tail );

		    }

		}

		// change flag to MCP converted
		flag_status = 5;
	}
	else// if not initialized
	{
		wxString msg;
		msg.Printf ( wxT ( "GEM not initialized. Load GEM data first!" ) );
		wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );
	}


}

void MyFrame::OnConICB ( wxCommandEvent &event )
{
	int i;
	double concentration;
	concentration = 0.0;

	if ( flag_status ) // if initialized
	{
		textCtrl->Clear();
		string temp_name,strnum; // used to get the current name of DC;

		for ( i=0 ; i < dCH->nIC ; i++ )
		{
			textCtrl->AppendText ( IC_Head );
			textCtrl->AppendText ( IC_PCS_Head );
			textCtrl->AppendText ( IC_PCS );
			textCtrl->AppendText ( IC_Name_Head );

			temp_name=stringify ( i+1 );
			temp_name.append ( "-" );
			strnum = dCH->ICNL[i]; // get component name;
			temp_name.append ( strnum );
			IC_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );

			textCtrl->AppendText ( IC_Name );
			textCtrl->AppendText ( IC_Geo_Head );
			textCtrl->AppendText ( IC_Geo );
			textCtrl->AppendText ( IC_Value_Head );
			textCtrl->AppendText ( IC_Value_Typ );

			// here we pass simply the scaled b-vector ...otherwise we need double amount of species initially (solutes and rest)

			IC_Value.Printf ( wxT ( "%1.14e\n" ), ( dBR->bIC[i] / dBR->Vs ) ); //scale to a volume of 1 m^3

			textCtrl->AppendText ( IC_Value );
			textCtrl->AppendText ( IC_End );
		}

		flag_status = 2;
	}
	else// if not initialized
	{
		wxString msg;
		msg.Printf ( wxT ( "GEM not initialized. Load GEM data first!" ) );
		wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );
	}
}

void MyFrame::OnConMCPB ( wxCommandEvent &event )
{
	int i;

	if ( flag_status ) // if initialized
	{
		textCtrl->Clear();
		string temp_name,strnum; // used to get the current name of DC;

		for ( i=0 ; i < dCH->nIC ; i++ )
		{

			MCP_Head.Printf ( wxT ( "#COMPONENT_PROPERTIES ; comp%i \n" ), i+1 );
			textCtrl->AppendText ( MCP_Head );
			textCtrl->AppendText ( MCP_Name_Head );


			temp_name=stringify ( i+1 );
			temp_name.append ( "-" );
			strnum = dCH->ICNL[i]; // get component name;
			temp_name.append ( strnum );
			MCP_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );

			textCtrl->AppendText ( MCP_Name );
			textCtrl->AppendText ( MCP_Mobile_Head );

			// for disolved species, mobile and with an diffusion number
			textCtrl->AppendText ( MCP_Mobile_Yes );
			textCtrl->AppendText ( MCP_Dif_Head );
			textCtrl->AppendText ( MCP_Dif_Yes );


			textCtrl->AppendText ( BC_End );
		}

		// change flag to MCP converted
		flag_status = 4;
	}
	else// if not initialized
	{
		wxString msg;
		msg.Printf ( wxT ( "GEM not initialized. Load GEM data first!" ) );
		wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );
	}


}

void MyFrame::OnConBCB ( wxCommandEvent &event )
{
	int i;
	double concentration;
	concentration = 0.0;

	if ( flag_status ) // if initialized
	{
		textCtrl->Clear();
		string temp_name,strnum; // used to get the current name of DC;

		for ( i=0 ; i < dCH->nIC ; i++ )
		{
			textCtrl->AppendText ( BC_Head );
			textCtrl->AppendText ( BC_PCS_Head );
			textCtrl->AppendText ( BC_PCS );
			textCtrl->AppendText ( BC_Name_Head );

			temp_name=stringify ( i+1 );
			temp_name.append ( "-" );
			strnum = dCH->ICNL[i]; // get component name;
			temp_name.append ( strnum );
			BC_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );

			textCtrl->AppendText ( BC_Name );
			textCtrl->AppendText ( BC_Geo_Head );
			textCtrl->AppendText ( BC_Geo );
			textCtrl->AppendText ( BC_Value_Head );
			textCtrl->AppendText ( BC_Value_Typ );

			// for disolved species, convert to concentration values in mol/L water
			// for BCs we need concentrations, as BCs are read at the beginning at overwritten each timestep
			// in the current implementation of geosys-gems the boundary conditions are not changed (expect for changes defined in the geosys imput files)
			      // small error in porosity during init_gems routine should not be a problem!
			     if ( dCH->ccIC[i] == IC_HYDROGEN )
			     {
				concentration = ( dBR->bPS[i]- ( 2*dBR-> xPA[0] ) ) / dBR->Vs; //scale to a volume of 1 m^3
			     }
			     else if ( dCH->ccIC[i] == IC_OXYGEN )
			     {
				concentration = ( dBR->bPS[i] -dBR-> xPA[0] ) / dBR->Vs; //scale to a volume of 1 m^3
			     }
			     else
			     {
				concentration = dBR->bPS[i] / dBR->Vs; //scale to a volume of 1 m^3
			     }
			concentration /= ( water_volume ) ; // now it is concentration
			BC_Value.Printf ( wxT ( "%1.14e\n" ), concentration ); // no further scaling with system volume required

			textCtrl->AppendText ( BC_Value );
			textCtrl->AppendText ( BC_End );
		}

		// change flag to BC converted
		flag_status = 3;
	}
	else// if not initialized
	{
		wxString msg;
		msg.Printf ( wxT ( "GEM not initialized. Load GEM data first!" ) );
		wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );
	}
}
void MyFrame::OnConBCBW ( wxCommandEvent &event )
{
	int i;
	double concentration;
	concentration = 0.0;

	if ( flag_status ) // if initialized
	{
		textCtrl->Clear();
		string temp_name,strnum; // used to get the current name of DC;

		for ( i=0 ; i < dCH->nIC ; i++ )
		{
			textCtrl->AppendText ( BC_Head );
			textCtrl->AppendText ( BC_PCS_Head );
			textCtrl->AppendText ( BC_PCS );
			textCtrl->AppendText ( BC_Name_Head );

			temp_name=stringify ( i+1 );
			temp_name.append ( "-" );
			strnum = dCH->ICNL[i]; // get component name;
			temp_name.append ( strnum );
			BC_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );

			textCtrl->AppendText ( BC_Name );
			textCtrl->AppendText ( BC_Geo_Head );
			textCtrl->AppendText ( BC_Geo );
			textCtrl->AppendText ( BC_Value_Head );
			textCtrl->AppendText ( BC_Value_Typ );

			// for disolved species, convert to concentration values in mol/L water
			// for BCs we need concentrations, as BCs are read at the beginning at overwritten each timestep
			// in the current implementation of geosys-gems the boundary conditions are not changed (expect for changes defined in the geosys imput files)
			concentration = dBR->bPS[i] / dBR->Vs; //scale to a volume of 1 m^3
			concentration /= ( water_volume ) ; // now it is concentration
			BC_Value.Printf ( wxT ( "%1.14e\n" ), concentration ); // no further scaling with system volume required

			textCtrl->AppendText ( BC_Value );
			textCtrl->AppendText ( BC_End );
		}

		// change flag to BC converted
		flag_status = 3;
	}
	else// if not initialized
	{
		wxString msg;
		msg.Printf ( wxT ( "GEM not initialized. Load GEM data first!" ) );
		wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );
	}
}

void MyFrame::OnConOUTB ( wxCommandEvent &event )
{
	int i;

	if ( flag_status ) // if initialized
	{
		textCtrl->Clear();
		string temp_name,strnum; // used to get the current name of DC;

		textCtrl->AppendText ( OUT_Head );

		for ( i=0 ; i < dCH->nIC ; i++ )
		{

			temp_name=stringify ( i+1 );
			temp_name.append ( "-" );
			strnum = dCH->ICNL[i]; // get component name;
			temp_name.append ( strnum );
			OUT_Name.Printf ( wxT ( "  %S \n" ), temp_name.c_str() );

			textCtrl->AppendText ( OUT_Name );

		}

		textCtrl->AppendText ( OUT_Tail );

		// change flag to MCP converted
		flag_status = 5;
	}
	else// if not initialized
	{
		wxString msg;
		msg.Printf ( wxT ( "GEM not initialized. Load GEM data first!" ) );
		wxMessageBox ( msg, wxT ( "GEM Initialization" ), wxOK | wxICON_INFORMATION, this );
	}


}



void MyFrame::OnConSet ( wxCommandEvent &event )
{


}

void MyFrame::OnSave ( wxCommandEvent &event )
{
	if ( flag_status > 1 )
	{
		wxString caption = wxT ( "Please give the path of output file" );
		wxString wildcard;
		if ( flag_status == 2 )
			wildcard = wxT ( "Geosys IC files (*.ic)|*.ic|Any files (*.*)|*.*" );
		else if ( flag_status == 3 )
			wildcard = wxT ( "Geosys BC files (*.bc)|*.bc|Any files (*.*)|*.*" );

		wxString defaultDir = wxT ( "c:\\" );
		wxString defaultFilename = wxT ( "Output" );

		wxFileDialog dialog ( this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE );

		if ( dialog.ShowModal() == wxID_OK )
		{
			wxString path = dialog.GetPath();

			wxTextFile file;

			if ( file.Create ( path ) )
			{
				int i;
				for ( i=0 ; i < textCtrl->GetNumberOfLines() +1 ; i++ )
				{
					file.AddLine ( textCtrl->GetLineText ( i ) );
				}
				file.Write();


			}
			else
			{
				if ( file.Open ( path ) ) // when the file already exist
				{
					int i;
					file.Clear();
					for ( i=0 ; i < textCtrl->GetNumberOfLines() +1 ; i++ )
					{
						file.AddLine ( textCtrl->GetLineText ( i ) );
					}
					file.Write();

					wxString msg;
					msg.Printf ( wxT ( "File saved successfully. " ) );
					wxMessageBox ( msg, wxT ( "Save File" ), wxOK | wxICON_INFORMATION, this );
				}
				else
				{
					wxString msg;
					msg.Printf ( wxT ( "Can NOT save file!!!" ) );
					wxMessageBox ( msg, wxT ( "Error" ), wxOK | wxICON_INFORMATION, this );
				}
			}
		}
	}
	else
	{
		wxString msg;
		msg.Printf ( wxT ( "Please first load GEM project file!!!" ) );
		wxMessageBox ( msg, wxT ( "Error" ), wxOK | wxICON_INFORMATION, this );
	}
}

void MyFrame::OnSaveAs ( wxCommandEvent &event )
{
	if ( flag_status > 1 )
	{
		wxString caption = wxT ( "Please give the path of output file" );
		wxString wildcard;
		wildcard = wxT ( "Any files (*.*)|*.*" );

		wxString defaultDir = wxT ( "c:\\" );
		wxString defaultFilename = wxT ( "Output" );

		wxFileDialog dialog ( this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE );

		if ( dialog.ShowModal() == wxID_OK )
		{
			wxString path = dialog.GetPath();

			wxTextFile file;

			if ( file.Create ( path ) )
			{
				int i;
				for ( i=0 ; i < textCtrl->GetNumberOfLines() +1 ; i++ )
				{
					file.AddLine ( textCtrl->GetLineText ( i ) );
				}
				file.Write();


			}
			else
			{
				if ( file.Open ( path ) ) // when the file already exist
				{
					int i;
					file.Clear();
					for ( i=0 ; i < textCtrl->GetNumberOfLines() +1 ; i++ )
					{
						file.AddLine ( textCtrl->GetLineText ( i ) );
					}
					file.Write();

					wxString msg;
					msg.Printf ( wxT ( "File saved successfully. " ) );
					wxMessageBox ( msg, wxT ( "Save File" ), wxOK | wxICON_INFORMATION, this );
				}
				else
				{
					wxString msg;
					msg.Printf ( wxT ( "Can NOT save file!!!" ) );
					wxMessageBox ( msg, wxT ( "Error" ), wxOK | wxICON_INFORMATION, this );
				}
			}
		}
	}
	else
	{
		wxString msg;
		msg.Printf ( wxT ( "Please first load GEM project file!!!" ) );
		wxMessageBox ( msg, wxT ( "Error" ), wxOK | wxICON_INFORMATION, this );
	}

}

MyFrame::MyFrame ( const wxString& title )
		:wxFrame ( NULL, wxID_ANY, title )
{
	// Set the frame icon
  //	SetIcon ( wxIcon ( wxwin32x32_xpm ) );

	// The "About" item should be in the help menu
	wxMenu *helpMenu = new wxMenu;
	helpMenu->Append ( wxID_ABOUT, wxT ( "&About...\tF1" ),
	                   wxT ( "Show about dialog" ) );

	// Create a file menu
	wxMenu *fileMenu = new wxMenu;
	fileMenu->Append ( wxID_OPEN, wxT ( "&Open GEM File\tAlt-O" ),
	                   wxT ( "Open GEM Initialization File" ) );
	fileMenu->Append ( wxID_SAVE, wxT ( "Save..." ),
	                   wxT ( "Save current result into IC or BC files" ) );
	fileMenu->Append ( wxID_SAVEAS, wxT ( "Save as..." ),
	                   wxT ( "Save current result into a diffent file" ) );
	fileMenu->Append ( wxID_EXIT, wxT ( "E&xit\tAlt-X" ),
	                   wxT ( "Quit this program" ) );



	// Create a convert menu
	//	wxMenu *convertMenu = new wxMenu;
//	convertMenu->Append ( ID_Con_IC, wxT ( "&IC format\tAlt-I" ),
//	                      wxT ( "Convert to IC format" ) );
//	convertMenu->Append ( ID_Con_BC, wxT ( "&BC format\tAlt-B" ),
//	                      wxT ( "Convert to BC format" ) );
//	convertMenu->Append ( ID_Con_MCP, wxT ( "&MCP format\tAlt-M" ),
//	                      wxT ( "Convert to MCP format" ) );
//	convertMenu->Append ( ID_Con_OUT, wxT ( "&OUT format\tAlt-M" ),
//	                      wxT ( "Convert to OUT format" ) );
	// convertMenu->Append(ID_Con_Set, wxT("Conv &Setting...\tAlt-S"),
	// 				 wxT("Convert Format Settings"));
	// Create a secondconvert menu
	wxMenu *convertBMenu = new wxMenu;
	convertBMenu->Append ( ID_Con_ICB, wxT ( "&IC format\tAlt-I" ),
	                       wxT ( "Convert to IC format" ) );
	convertBMenu->Append ( ID_Con_BCB, wxT ( "&BC format\tAlt-B" ),
	                       wxT ( "Convert to BC format" ) );
	convertBMenu->Append ( ID_Con_BCBW, wxT ( "&BC format with water\tAlt-W" ),
	                       wxT ( "BC format with water in concentrations" ) );
	convertBMenu->Append ( ID_Con_MCPB, wxT ( "&MCP format\tAlt-M" ),
	                       wxT ( "Convert to MCP format" ) );
	convertBMenu->Append ( ID_Con_OUTB, wxT ( "&OUT format\tAlt-O" ),
	                       wxT ( "Convert to OUT format" ) );
	convertBMenu->Append ( ID_Con_DLUB, wxT ( "&constraints format\tAlt-L" ),
	                       wxT ( "Convert to constraints format" ) );


	// Now append the freshly created menu to the menu bar...
	wxMenuBar *menuBar = new wxMenuBar();
	menuBar->Append ( fileMenu, wxT ( "&File" ) );
	//	menuBar->Append ( convertMenu, wxT ( "&Convert to Bvector with water" ) );
	menuBar->Append ( convertBMenu, wxT ( "&Convert to Bvector" ) );
	menuBar->Append ( helpMenu, wxT ( "&Help" ) );

	// ... and attach this menu bar to the frame
	SetMenuBar ( menuBar );

	// Creat the txtCtrl
	textCtrl = new wxTextCtrl ( this, ID_TEXTCTRL,
	                            wxEmptyString, wxDefaultPosition, wxSize ( 240, 100 ),
	                            wxTE_MULTILINE );

	// Creat a status bar just for fun
	CreateStatusBar ( 2 );
	SetStatusText ( wxT ( "Welcome to GEM to GeoSys Convertor!" ) );

	// Initialize the TNode class
		m_Node = new TNode();

	// Set the default values of the output strings
	BC_Head      = wxT ( "#BOUNDARY_CONDITION\n" );
	IC_Head      = wxT ( "#INITIAL_CONDITION\n" );
	BC_Geo       = wxT ( "  POINT POINT0\n" );
	IC_Geo       = wxT ( "  DOMAIN  \n" );
	BC_Geo_Head  = wxT ( " $GEO_TYPE\n" );
	IC_Geo_Head  = wxT ( " $GEO_TYPE\n" );
	BC_End       = wxT ( " \n" );
	IC_End       = wxT ( " \n" );
	BC_PCS		 = wxT ( "  MASS_TRANSPORT\n" );
	IC_PCS		 = wxT ( "  MASS_TRANSPORT\n" );
	BC_PCS_Head  = wxT ( " $PCS_TYPE\n" );
	IC_PCS_Head  = wxT ( " $PCS_TYPE\n" );
	BC_Name      = wxT ( " \n" );
	IC_Name      = wxT ( " \n" );
	BC_Name_Head = wxT ( " $PRIMARY_VARIABLE\n" );
	IC_Name_Head = wxT ( " $PRIMARY_VARIABLE\n" );
	BC_Value     = wxT ( " 0.000000000000000000e+00\n" );
	IC_Value     = wxT ( " 0.000000000000000000e+00\n" );
	BC_Value_Typ = wxT ( "  CONSTANT  " );
	IC_Value_Typ = wxT ( "  CONSTANT  " );
	BC_Value_Head = wxT ( " $DIS_TYPE\n" );
	IC_Value_Head = wxT ( " $DIS_TYPE\n" );

	MCP_Head      = wxT ( "#COMPONENT_PROPERTIES\n" );
	MCP_Name_Head = wxT ( " $NAME  \n" );
	MCP_Name      = wxT ( "  NO NAME  \n" );
	MCP_Mobile_Head=wxT ( " $MOBILE\n" );
	MCP_Mobile_Yes= wxT ( "  1;  MOBILE-Flag: 0=imMOBILEe, 1=MOBILEe/transported\n" );
	MCP_Mobile_No = wxT ( "  0;  MOBILE-Flag: 0=imMOBILEe, 1=MOBILEe/transported\n" );
	MCP_Dif_Head  = wxT ( " $DIFFUSION\n" );
	MCP_Dif_Yes   = wxT ( "  1  1.0e-10\n" );
	MCP_Dif_No    = wxT ( "  0  \n" );

	OUT_Head=   wxT ( "#OUTPUT\n $NOD_VALUES\n" );
	OUT_Tail=   wxT ( "  pH\n  pe\n  Eh\n  NodePorosity\n $ELE_VALUES\n  POROSITY\n $GEO_TYPE\n  DOMAIN\n $DAT_TYPE\n  VTK\n $TIM_TYPE\n   STEPS 316\n#STOP\n" );

        DLU_Head=   wxT (" $CONSTRAINT_GEMS\n");
        DLU_Tail=   wxT ("   ; dll, dul, 1, number of phase for criteria, lower amount, upper amount of phase\n");
}



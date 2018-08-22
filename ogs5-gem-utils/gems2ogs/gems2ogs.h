#ifndef INCLUDED_HELLOWORLDAPP_H
#define INCLUDED_HELLOWORLDAPP_H

#include "node.h"

//

 
// The HelloWorldApp class. This class shows a window
// containing a statusbar with the text "Hello World"
class HelloWorldApp : public wxApp
{
public:
	virtual bool OnInit();
};

class MyFrame : public wxFrame
{
public: 
	// Constructor
	MyFrame(const wxString& title);

	// Event handlers
	void OnQuit(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	void OnOpenDCH(wxCommandEvent& event);
	void OnConIC(wxCommandEvent& event);
	void OnConBC(wxCommandEvent& event);
        void OnConMCP(wxCommandEvent& event);
        void OnConOUT(wxCommandEvent& event);
	void OnConICB(wxCommandEvent& event);
	void OnConBCB(wxCommandEvent& event);
	void OnConBCBW(wxCommandEvent& event);
        void OnConMCPB(wxCommandEvent& event);
        void OnConOUTB(wxCommandEvent& event);
        void OnConDLUB(wxCommandEvent& event);

	void OnConSet(wxCommandEvent& event);
	void OnSave(wxCommandEvent& event);
	void OnSaveAs(wxCommandEvent& event);

	// Ctrl of the text shown in main frame
	wxTextCtrl* textCtrl;
	
	// pointer to the Node-GEMS structure
    TNode* m_Node; // Instance of TNode class
    DATACH* dCH;   // pointer to DATACH
    DATABR* dBR;   // pointer to DATABR

	// Strings used to generate BC and IC files
	wxString BC_Head, IC_Head, BC_Geo,      IC_Geo,  BC_Geo_Head,   IC_Geo_Head, 
			 BC_End,  IC_End,  BC_PCS,      IC_PCS,  BC_PCS_Head,   IC_PCS_Head,
			 BC_Name, IC_Name, BC_Name_Head,IC_Name_Head, BC_Value, IC_Value,
			 BC_Value_Typ, IC_Value_Typ,    IC_Value_Head, BC_Value_Head;	
	wxString MCP_Head, MCP_Name_Head, MCP_Name, MCP_Mobile_Head,
             MCP_Mobile_Yes, MCP_Mobile_No, MCP_Dif_Head, 
			 MCP_Dif_Yes, MCP_Dif_No;
	wxString OUT_Head, OUT_Tail, OUT_Name;
	wxString DLU_Head, DLU_Name, DLU_Tail, DLU_Number;
	// flag. 1-GEM initialized; 0-not initialized; 2- IC converted; 3- BC converted;
	int flag_status;

	// system information
	double water_volume;
	int    idx_water;
	int    nDC; 

private:
	// This calss handles events
	DECLARE_EVENT_TABLE();

	// declaration of the IDs
	enum
	{
		ID_Con_IC   = 801,
		ID_Con_BC   = 802,
		ID_Con_Set  = 803,
                ID_Con_MCP  = 804,
                ID_Con_OUT  = 805,
		ID_Con_ICB   = 806,
		ID_Con_BCB   = 807,
                ID_Con_MCPB  = 808,
                ID_Con_OUTB  = 809,
		ID_Con_BCBW   = 810,
		ID_Con_DLUB   = 811,
		ID_TEXTCTRL = 821,
	};
};

DECLARE_APP(HelloWorldApp)
 
#endif // INCLUDED_HELLOWORLDAPP_H 

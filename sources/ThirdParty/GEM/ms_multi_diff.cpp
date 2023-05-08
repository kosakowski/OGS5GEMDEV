//-------------------------------------------------------------------
// $Id$
//
/// \file ms_multi_diff.cpp
/// Implementation of coping IPM internal structure
//
// Copyright (c) 2017-2020 S.Dmytriyeva, D.Kulik
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#include <fstream>
#include "node.h"
#include "num_methods.h"
#include "v_service.h"
#ifdef IPMGEMPLUGIN
#include "ms_multi.h"
#else
#include "m_param.h"
#include "stepwise.h"
#endif


#ifdef IPMGEMPLUGIN

SPP_SETTING pa_ = {
    " GEMS-GUI v.3.1 r.2150 (rc) " " GEMS3K v.3.1 r.690 (rc) ",
    {    // Typical default set (03.04.2012) new PSSC( logSI ) & uDD()
         2,  /* PC */  2,     /* PD */   -5,   /* PRD */
         1,  /* PSM  */ 130,  /* DP */   1,   /* DW */
         0, /* DT */     30000,   /* PLLG */   1,  /* PE */  7000, /* IIM */
         1000., /* DG */   1e-13,  /* DHB */  1e-20,  /* DS */
         1e-6,  /* DK */  0.01,  /* DF */  0.01,  /* DFM */
         1e-5,  /* DFYw */  1e-5,  /* DFYaq */    1e-5,  /* DFYid */
         1e-5,  /* DFYr,*/  1e-5,  /* DFYh,*/   1e-5,  /* DFYc,*/
         1e-6, /* DFYs, */  1e-17,  /* DB */   1.,   /* AG */
         0.,   /* DGC */   1.0,   /* GAR */  1000., /* GAH */
         1e-3, /* GAS */   12.05,  /* DNS */   1e-13,  /* XwMin, */
         1e-13,  /* ScMin, */  1e-33, /* DcMin, */   1e-20, /* PhMin, */
         1e-5,  /* ICmin */   1e-10,  /* EPS */   1e-3,  /* IEPS */
         1e-10,  /* DKIN  */ 0,  /* tprn */
    },
}; // SPP_SETTING

#endif

const BASE_PARAM& pa_p_= pa_.p;

#ifdef IPMGEMPLUGIN
TMulti::TMulti( TNode* na_ )
{
    pa_standalone.reset( new BASE_PARAM() );
    *pa_standalone = pa_p_;

    pmp = &pm;
    node1 = na_; // parent
    sizeN = 0;
    AA = nullptr;
    BB = nullptr;
    arrL = nullptr;
    arrAN = nullptr;

    U_mean = nullptr;
    U_M2 = nullptr;
    U_CVo = nullptr;
    U_CV = nullptr;
    ICNud = nullptr;
    pm.errorCode[0] ='\0';
    pm.errorBuf[0] ='\0';

    sizeFIs = 0;
    phSolMod = nullptr;
    sizeFIa = 0;
    phSorpMod = nullptr;
    sizeFI = 0;
    phKinMet = nullptr;

    pmp->Guns = nullptr;
    pmp->Vuns = nullptr;
    pmp->tpp_G = nullptr;
    pmp->tpp_S = nullptr;
    pmp->tpp_Vm = nullptr;
    set_def();
}
#endif

bool TMulti::testTSyst() const
{
#ifndef IPMGEMPLUGIN
    return ( TSyst::sm->GetSY()->PYOF != S_OFF);
#else
    return true;
#endif
}

void TMulti::get_PAalp_PSigm( char& PAalp, char& PSigm)
{
#ifndef IPMGEMPLUGIN
    PAalp =  TSyst::sm->GetSY()->PAalp;
    PSigm =  TSyst::sm->GetSY()->PSigm;
#else
    PAalp = PAalp_;
    PSigm = PSigm_;
#endif
}

void TMulti::STEP_POINT( const char* str)
{
#ifndef IPMGEMPLUGIN
    STEP_POINT1(str);
#else
    ipm_logger->debug( std::string("STEP_POINT ")+str);
#endif
}


void TMulti::alloc_IPx( long int LsIPxSum )
{
#ifdef IPMGEMPLUGIN
    if( pm.IPx ) delete[] pm.IPx;
    pm.IPx = new long int[ LsIPxSum];
#else
    pm.IPx = (long int *)aObj[ o_wi_ipxpm ]->Alloc(LsIPxSum, 1, L_);
#endif
}

void TMulti::alloc_PMc( long int LsModSum )
{
#ifdef IPMGEMPLUGIN
    if( pm.PMc ) delete[] pm.PMc;
    pm.PMc = new double[LsModSum];
#else
    pm.PMc = (double *)aObj[ o_wi_pmc]->Alloc( LsModSum, 1, D_);
#endif
}

void TMulti::alloc_DMc( long int LsMdcSum )
{
#ifdef IPMGEMPLUGIN
    if( pm.DMc ) delete[] pm.DMc;
    pm.DMc = new double[LsMdcSum];
#else
    pm.DMc = (double *)aObj[ o_wi_dmc]->Alloc( LsMdcSum, 1, D_ );
#endif
}

void TMulti::alloc_MoiSN( long int LsMsnSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.MoiSN) delete[] pm.MoiSN;
    pm.MoiSN = new double[LsMsnSum];
#else
    pm.MoiSN = (double *)aObj[ o_wi_moisn]->Alloc( LsMsnSum, 1, D_ );
#endif
}

void TMulti::alloc_SitFr( long int LsSitSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.SitFr) delete[] pm.SitFr;
    pm.SitFr = new double[LsSitSum];
#else
    pm.SitFr  = (double *)aObj[ o_wo_sitfr ]->Alloc( LsSitSum, 1, D_ );
#endif
}

void TMulti::alloc_DQFc( long int DQFcSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.DQFc) delete[] pm.DQFc;
    pm.DQFc = new double[DQFcSum];
#else
    pm.DQFc = (double *)aObj[o_wi_dqfc]->Alloc( DQFcSum, 1, D_ );
#endif
}

void TMulti::alloc_PhLin( long int PhLinSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.PhLin) delete[] pm.PhLin;
    pm.PhLin = new long int[PhLinSum][2];
#else
    pm.PhLin = (long int (*)[2])aObj[ o_wi_phlin]->Alloc( PhLinSum, 2, L_ );
#endif
}

void TMulti::alloc_lPhc( long int lPhcSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.lPhc) delete[] pm.lPhc;
    pm.lPhc = new double[lPhcSum];
#else
    pm.lPhc = (double*)aObj[o_wi_lphc]->Alloc( lPhcSum, 1, D_ );
#endif
}

void TMulti::alloc_xSMd( long int xSMdSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.xSMd) delete[] pm.xSMd;
    pm.xSMd = new long int[xSMdSum];
#else
    pm.xSMd = (long int*)aObj[ o_wi_xsmd]->Alloc( xSMdSum, 1, L_ );
#endif
}

void TMulti::alloc_IsoPc( long int IsoPcSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.IsoPc) delete[] pm.IsoPc;
    pm.IsoPc = new double[IsoPcSum];
#else
    pm.IsoPc = (double*)aObj[ o_wi_isopc]->Alloc( IsoPcSum, 1, D_ );
#endif
}

void TMulti::alloc_IsoSc( long int IsoScSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.IsoSc) delete[] pm.IsoSc;
    pm.IsoSc = new double[IsoScSum];
#else
    pm.IsoSc = (double*)aObj[ o_wi_isosc]->Alloc( IsoScSum, 1, D_ );
#endif
}

void TMulti::alloc_IsoCt( long int IsoCtSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.IsoCt) delete[] pm.IsoCt;
    pm.IsoCt = new char[IsoCtSum];
#else
    pm.IsoCt = (char*)aObj[ o_wi_isoct]->Alloc( IsoCtSum, 1, A_ );
#endif
}

void TMulti::alloc_EImc( long int EImcSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.EImc) delete[] pm.EImc;
    pm.EImc = new double[EImcSum];
#else
    pm.EImc = (double*)aObj[ o_wi_eimc]->Alloc( EImcSum, 1, D_ );
#endif
}

void TMulti::alloc_mCDc( long int mCDcSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.mCDc) delete[] pm.mCDc;
    pm.mCDc = new double[mCDcSum];
#else
    pm.mCDc = (double*)aObj[ o_wi_mcdc]->Alloc( mCDcSum, 1, D_ );
#endif
}

void TMulti::alloc_xSKrC( long int xSKrCSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.xSKrC) delete[] pm.xSKrC;
    pm.xSKrC = new long int[xSKrCSum];
#else
    pm.xSKrC = (long int*)aObj[ o_wi_jcrdc]->Alloc( xSKrCSum, 1, L_ );
#endif
}

void TMulti::alloc_ocPRkC( long int ocPRkC_feSArC_Sum )
{
#ifdef IPMGEMPLUGIN
    if(pm.ocPRkC) delete[] pm.ocPRkC;
    pm.ocPRkC = new long int[ocPRkC_feSArC_Sum][2];
#else
    pm.ocPRkC = (long int(*)[2])aObj[ o_wi_ocprkc]->Alloc( ocPRkC_feSArC_Sum, 2, L_ );
#endif
}

void TMulti::alloc_feSArC( long int ocPRkC_feSArC_Sum )
{
#ifdef IPMGEMPLUGIN
    if(pm.feSArC) delete[] pm.feSArC;
    pm.feSArC = new double[ocPRkC_feSArC_Sum];
#else
    pm.feSArC = (double*)aObj[ o_wi_fsac]->Alloc( ocPRkC_feSArC_Sum, 1, D_ );
#endif
}

void TMulti::alloc_rpConC( long int rpConCSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.rpConC) delete[] pm.rpConC;
    pm.rpConC = new double[rpConCSum];
#else
    pm.rpConC = (double*)aObj[ o_wi_krpc]->Alloc( rpConCSum, 1, D_ );
#endif
}

void TMulti::alloc_apConC( long int apConCSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.apConC) delete[] pm.apConC;
    pm.apConC = new double[apConCSum];
#else
    pm.apConC = (double*)aObj[ o_wi_apconc]->Alloc( apConCSum, 1, D_ );
#endif
}

void TMulti::alloc_AscpC( long int AscpCSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.AscpC) delete[] pm.AscpC;
    pm.AscpC = new double[AscpCSum];
#else
    pm.AscpC = (double*)aObj[ o_wi_ascpc]->Alloc( AscpCSum, 1, D_ );
#endif
}

void TMulti::alloc_UMpcC( long int UMpcSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.UMpcC) delete[] pm.UMpcC;
    pm.UMpcC = new double[UMpcSum];
#else
    pm.UMpcC = (double*)aObj[ o_wi_umpc]->Alloc( UMpcSum, 1, D_ );
#endif
}

void TMulti::alloc_xICuC( long int xICuCSum )
{
#ifdef IPMGEMPLUGIN
    if(pm.xICuC) delete[] pm.xICuC;
    pm.xICuC = new long int[xICuCSum];
#else
    pm.xICuC = (long int *)aObj[o_wi_xicuc ]->Alloc( xICuCSum, 1, L_ );
#endif

}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculating demo partial pressures of gases (works only in GEMS-PSI)
//
void TMulti::GasParcP()
{
#ifndef IPMGEMPLUGIN
    long int k,  i, jj=0;
    long int jb, je, j;

    if( !pm.PG )
        return;

    char (*SMbuf)[MAXDCNAME] =
            (char (*)[MAXDCNAME])aObj[ o_w_tprn]->Alloc( pm.PG, 1, MAXDCNAME );
    pm.Fug = (double *)aObj[ o_wd_fug]->Alloc( pm.PG, 1, D_ );
    pm.Fug_l = (double *)aObj[ o_wd_fugl]->Alloc( pm.PG, 1, D_ );
    pm.Ppg_l = (double *)aObj[ o_wd_ppgl]->Alloc( pm.PG, 1, D_ );

    for( k=0, je=0; k<pm.FIs; k++ ) // phase
    {
        jb = je;
        je = jb+pm.L1[k];
        if( pm.PHC[k] == PH_GASMIX || pm.PHC[k] == PH_PLASMA
                || pm.PHC[k] == PH_FLUID )
        {
            for( j=jb; j<je; j++,jj++ )
            {  // fixed 02.03.98 DK

                copyValues(SMbuf[jj], pm.SM[j], MAXDCNAME );
                pm.Fug_l[jj] = -(pm.G0[j] + pm.fDQF[j]);
                if( pm.P > 1e-9 )
                    pm.Fug_l[jj] += log(pm.P);
                for( i=0; i<pm.N; i++ )
                    pm.Fug_l[jj] += *(pm.A+j*pm.N+i) * pm.U[i];
                if( pm.Fug_l[jj] > -37. && pm.Fug_l[jj] < 16. )
                    pm.Fug[jj] = exp( pm.Fug_l[jj] );
                else  pm.Fug[jj] = 0.0;
                // Partial pressure
                pm.Ppg_l[jj] = pm.Fug_l[jj] - pm.lnGam[j];
                pm.Fug_l[jj] *= .43429448;
                pm.Ppg_l[jj] *= .43429448;
            }
            // break;
        }
    }
#endif
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Linking DOD for executing Phase mixing model scripts
void TMulti::pm_GC_ods_link( long int k, long int jb, long int jpb, long int jdb, long int ipb )
{
#ifndef IPMGEMPLUGIN
    ErrorIf( k < 0 || k >= pm.FIs , "CalculateActivityCoefficients():", "Invalid link: k=0||>FIs" );
    aObj[ o_nsmod]->SetPtr( pm.sMod[k] );
    aObj[ o_nncp]->SetPtr( pm.LsMod+k*3 );
    aObj[ o_nncd]->SetPtr( pm.LsMdc+k*3 );
    aObj[ o_ndc]->SetPtr(  pm.L1+k );
    aObj[ o_nez]->SetPtr( pm.EZ+jb );
    aObj[o_nez]->SetN(  pm.L1[k]);
    aObj[ o_npcv]->SetPtr( pm.PMc+jpb );
    aObj[o_npcv]->SetDim( pm.LsMod[k*3], pm.LsMod[k*3+2]);
    //  Object for indexation of interaction parameters
    aObj[ o_nu]->SetPtr( pm.IPx+ipb ); // added 07.12.2006  KD
    aObj[o_nu]->SetDim( pm.LsMod[k*3], pm.LsMod[k*3+1]);
    //
    aObj[ o_ndcm]->SetPtr( pm.DMc+jdb );
    aObj[o_ndcm]->SetDim( pm.L1[k], pm.LsMdc[k*3] );
    aObj[ o_nmvol]->SetPtr( pm.Vol+jb );
    aObj[o_nmvol]->SetN( pm.L1[k]);
    aObj[ o_nppar]->SetPtr(pm.G0+jb );  // changed 10.12.2008 by DK
    aObj[o_nppar]->SetN(  pm.L1[k]);
    //    aObj[ o_ngtn]->SetPtr( pm.G0+jb );
    aObj[ o_ngtn]->SetPtr( pm.fDQF+jb );     // changed 05.12.2006 by DK
    aObj[o_ngtn]->SetN( pm.L1[k] );
    aObj[ o_ngam]->SetPtr( pm.Gamma+jb ); // Gamma calculated
    aObj[o_ngam]->SetN( pm.L1[k] );
    aObj[ o_nlngam]->SetPtr( pm.lnGam+jb ); // ln Gamma calculated
    aObj[o_nlngam]->SetN( pm.L1[k]);
    aObj[ o_nas]->SetPtr(  pm.A+pm.N*jb );
    aObj[o_nas]->SetDim(  pm.L1[k], pm.N );
    aObj[ o_nxa]->SetPtr(  pm.XF+k );
    aObj[ o_nxaa]->SetPtr(  pm.XFA+k );
    if( pm.FIat > 0 )
    {
        aObj[ o_nxast]->SetPtr( pm.XFTS[k] );
        aObj[ o_nxcec]->SetPtr( pm.MASDT[k] );
    }
    else
    {
        aObj[ o_nxast]->SetPtr( 0 );
        aObj[ o_nxcec]->SetPtr( 0 );
    }
    //
    aObj[ o_nbmol]->SetPtr( pm.FVOL+k );  // phase volume
    aObj[ o_nxx]->SetPtr(  pm.X+jb );
    aObj[o_nxx]->SetN( pm.L1[k]);
    aObj[ o_nwx]->SetPtr(  pm.Wx+jb );
    aObj[o_nwx]->SetN( pm.L1[k]);
    aObj[ o_nmju]->SetPtr( pm.Fx+jb );
    aObj[o_nmju]->SetN( pm.L1[k]);
    aObj[ o_nqp]->SetPtr( pm.Qp+k*QPSIZE );
    aObj[ o_nqd]->SetPtr( pm.Qd+k*QDSIZE );   // Fixed 7.12.04 by KD

    // phase excess properties
    aObj[o_ngte]->SetPtr( &pm.GPh[k][0] );
    aObj[o_nhte]->SetPtr( &pm.HPh[k][0] );
    aObj[o_nste]->SetPtr( &pm.SPh[k][0] );
    aObj[o_nvte]->SetPtr( &pm.VPh[k][0] );
    aObj[o_ncpte]->SetPtr( &pm.CPh[k][0] );
    aObj[o_nate]->SetPtr( &pm.APh[k][0] );
    aObj[o_nute]->SetPtr( &pm.UPh[k][0] );
#endif
}


/// Output to "ipmlog.txt" file Warnings
long int TMulti::testMulti()
{
#ifdef IPMGEMPLUGIN
    if( pm.MK || pm.PZ )
    {
        if( base_param()->PSM >= 2 )
        {
            TNode::ipmlog_file->warn(" {} : {}:{}", char_array_to_string(pm.stkey, EQ_RKLEN), pm.errorCode, pm.errorBuf);
        }
        return 1L;
    }
    return 0L	;
#else

    if( pm.MK || pm.PZ )
    {
        if( base_param()->PSM >= 2 )
        {
            TNode::ipmlog_file->warn(" {} : {}:{}", char_array_to_string(pm.stkey, EQ_RKLEN), pm.errorCode, pm.errorBuf);
        }
        if( showMss )
        {
            addErrorMessage(" \nContinue?");
            switch( vfQuestion3(0, pm.errorCode, pm.errorBuf,
                                "&Yes", "&No", "&Yes to All" ))
            {
            case VF3_3:
                showMss=0l;
            case VF3_1:
                break;
            case VF3_2:
                Error(pmp->errorCode, pmp->errorBuf);
            }
        }
        return 1L;
    }
    return 0L;
#endif
}

bool TMulti::calculateActivityCoefficients_scripts( long int LinkMode, long int k, long int jb,
                                                    long int jpb, long int jdb, long int ipb, double pmpXFk)
{
#ifndef IPMGEMPLUGIN
    // This part running Phase math scripts is not used in standalone GEMS3K
    // Link DOD and set sizes of work arrays
    pm_GC_ods_link( k, jb, jpb, jdb, ipb );
    pm.is=0;
    pm.js=0;
    pm.next=1;
    char* sMod = pm.sMod[k];
    const BASE_PARAM *pa_p = base_param();

    switch( LinkMode )
    { // check the calculation mode
    case LINK_TP_MODE: // running TP-dependent scripts
        if(( sMod[SPHAS_DEP] == SM_TPDEP || sMod[SPHAS_DEP] == SM_UXDEP ) && qEp[k]->nEquat() )
        {	// Changed on 26.02.2008 to try TW DQF scripts - DK
            qEp[k]->CalcEquat();
        }
        if((sMod[DCOMP_DEP] == SM_TPDEP || sMod[DCOMP_DEP] == SM_UXDEP) && qEd[k]->nEquat() )
        {
            switch( sMod[DCE_LINK] )
            {
            case SM_PUBLIC:  // one script for all species
                for( pm.js=0, pm.is=0; pm.js<pm.L1[k]; pm.js++ )
                    qEd[k]->CalcEquat();
                break;
            case SM_PRIVATE_: // separate group of equations per species
                qEd[k]->CalcEquat();
                break;
            }
        }
        break;

    case LINK_PP_MODE: // Mode of calculation of integral solution phase properties
        switch( pm.PHC[k] )
        {
        case PH_AQUEL:
        case PH_LIQUID:
        case PH_SINCOND:
        case PH_SINDIS:
        case PH_HCARBL:
        case PH_SIMELT:
        case PH_IONEX:
        case PH_ADSORPT:
        case PH_GASMIX:
        case PH_PLASMA:
        case PH_FLUID:  // How to pull this stuff out of the script (pointers to integral phase properties added)
            // SolModExcessProp( k, sMod[SPHAS_TYP] ); // extracting integral phase properties
            // SolModIdealProp( jb, k, sMod[SPHAS_TYP] );
            // SolModStandProp( jb, k, sMod[SPHAS_TYP] );
            // SolModDarkenProp( jb, k, sMod[SPHAS_TYP] );
            break;
        default:
            break;
        }
        break;

    case LINK_UX_MODE:  // the model is dependent on current concentrations on IPM iteration
        switch( pm.PHC[k] )
        {  //
        case PH_AQUEL:
            if(!(pmpXFk > pm.DSM && pm.X[pm.LO] > pm.XwMinM && pm.IC > pa_p->ICmin ))
                return false;
            break;
        case PH_GASMIX:
        case PH_PLASMA:
        case PH_FLUID:
            if( !(pmpXFk > pm.DSM && pm.XF[k] > pa_p->PhMin))
                return false;
            break;
        case PH_LIQUID:
        case PH_SIMELT:
        case PH_SINCOND:
        case PH_SINDIS:
        case PH_IONEX:
        case PH_ADSORPT:
        case PH_HCARBL:  // solid and liquid mixtures
            if( !(pmpXFk > pm.DSM) )
                return false;
            SolModActCoeff( k, sMod[SPHAS_TYP] );  // Added to introduce multi-site ideal term 29.11.2010
            break;
        case PH_POLYEL:  // PoissonBoltzmann( q, jb, je, k ); break;
        case PH_SORPTION: // electrostatic potenials from Gouy-Chapman eqn
            if( !(pm.PHC[0] == PH_AQUEL && pmpXFk > pm.DSM
                  && (pm.XFA[0] > pm.XwMinM && pm.XF[0] > pm.DSM )))
                return false;
            break;
        default:
            return false;
        } // end switch

        if( sMod[SPHAS_DEP] == SM_UXDEP && qEp[k]->nEquat() )
            // Equations for the whole phase
            qEp[k]->CalcEquat();
        if( sMod[DCOMP_DEP] == SM_UXDEP && qEd[k]->nEquat() )
        {  // Equations for species
            switch( sMod[DCE_LINK] )
            {
            case SM_PUBLIC:  // one script for all species
                for( pm.js=0, pm.is=0; pm.js<pm.L1[k]; pm.js++ )
                    qEd[k]->CalcEquat();
                break;
            case SM_PRIVATE_:  // separate group of equations for each species
                qEd[k]->CalcEquat();
                break;
            }
        }
        break;
    default:
        Error("CalculateActivityCoefficients()","Invalid LinkMode for a scripted solution model");
    } // end switch

#endif
    return true;
}

// Before Calculations
// Calculation by IPM (preparing for calculation, unpacking data) in GUI part
void TMulti::initalizeGEM_IPM_Data_GUI()
{
#ifndef IPMGEMPLUGIN
    // for GEMIPM unpackDataBr( bool uPrimalSol );
    // to define quantities

    bool newInterval = false;

    //   MultiKeyInit( key ); //into PMtest

    ipm_logger->trace(" pm.pBAL =  {}", pm.pBAL);
    if( !pm.pBAL )
        newInterval = true;    // to rebuild lookup arrays

    if( pm.pBAL < 2  || pm.pTPD < 2 )
    {
        SystemToLookup();
    }
    if( pm.pBAL < 2  )
    {
        // Allocating list of phases currently present in non-zero quantities
        MultiSystemInit( );
    }

    // Allocating list of phases currently present in non-zero quantities
    if( !pm.SFs )
        pm.SFs = (char (*)[MAXPHNAME+MAXSYMB])aObj[ o_wd_sfs]->Alloc(
                    pm.FI, 1, MAXPHNAME+MAXSYMB );

    // no old solution => must be simplex
    if( pm.pESU == 0 )
        pm.pNP = 0;

    // build new TNode
    if( !node1 )
    {
        node1 = new TNode( pmp );
        newInterval = true;
    }
    else if( !node1->TestTPGrid(pm.Tai, pm.Pai ))
        newInterval = true;

    if( newInterval )
    {   // build/rebuild internal lookup arrays
        node1->MakeNodeStructures(window(), true, pm.Tai, pm.Pai );
    }

    ipm_logger->trace("newInterval = {}   pm.pTPD =  {}", newInterval, pm.pTPD);
    // New: TKinMet stuff
    if( pm.pKMM <= 0 )
    {
        KinMetModLoad();  // Call point to loading parameters for kinetic models
        pm.pKMM = 1;
    }
#endif
}

void TMulti::multiConstInit_PN()
{
#ifndef IPMGEMPLUGIN
    pm.PZ = 0; // IPM default
    //  pm.FitVar[0] = pa->aqPar[0]; // setting T,P dependent b_gamma parameters
    //  pm.pH = pm.Eh = pm.pe = 0.0;
#else
    pm.PZ = base_param()->DW;  // in IPM
    //  pm.FitVar[0] = 0.0640000030398369;
#endif
}

void TMulti::GEM_IPM_Init_gui1()
{
#ifndef IPMGEMPLUGIN
    if( pm.pIPN <=0 )  // mixing models finalized in any case (AIA or SIA)
    {
        // not done if these models are already present in MULTI !
        pm.PD = TProfil::pm->pa.p.PD;
        SolModLoad();   // Call point to loading scripts and parameters for mixing models
    }
    pm.pIPN = 1;
#endif
}

void TMulti::GEM_IPM_Init_gui2()
{
#ifndef IPMGEMPLUGIN
    // dynamic work arrays - loading initial data
    for(int k=0; k<pm.FI; k++ )
    {
        pm.XFs[k] = pm.YF[k];
        pm.Falps[k] = pm.Falp[k];
        memcpy( pm.SFs[k], pm.SF[k], MAXPHNAME+MAXSYMB );
    }
#endif
}

/// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
//
void TMulti::DC_LoadThermodynamicData(TNode* aNa ) // formerly CompG0Load()
{
    double TK, PPa;

#ifndef IPMGEMPLUGIN
    TNode* na;
    if( aNa )
    {
        na = aNa;// for reading GEMIPM files task
        TK =  na->cTK();
        PPa = na->cP();
    }
    else
    {
        na = node1;
        TK =  pm.TC+C_to_K;
        PPa = pm.P*bar_to_Pa;
    }
#else
    TNode* na = node1;
    TK =  na->cTK();
    PPa = na->cP();
#endif

#ifndef IPMGEMPLUGIN
    if( !aNa )
    {
        double T = TK-C_to_K;
        TMTparm::sm->GetTP()->curT=T;
        TMTparm::sm->GetTP()->curP=PPa/bar_to_Pa;
    }
#endif

    // try generate thermodynamic data from ThermoEngine
    if( !na->load_all_thermodynamic_from_thermo( TK, PPa ))
    {
        load_all_thermodynamic_from_grid(na, TK, PPa );
    }
    pm.pTPD = 2;
}

/// Load Thermodynamic Data from DATACH to MULTI using Lagrangian Interpolator
void TMulti::load_all_thermodynamic_from_grid(TNode* aNa, double TK, double PPa )
{
    long int j, jj, k, xTP, jb, je=0;
    double Go, Gg=0., Ge=0., Vv, h0=0., S0 = 0., Cp0= 0., a0 = 0., u0 = 0.;
    double P = PPa/bar_to_Pa;
    DATACH  *dCH = aNa->pCSD();

//    ipm_logger->info("Calc Lookup T: {}  P: {}", TK, PPa);                 Temporarily disabled 23.Jan.2022
    if( dCH->nTp <1 || dCH->nPp <1 || aNa->check_TP( TK, PPa ) == false )
    {
        Error("load_all_thermodynamic_from_grid: ",
               std::string(" Temperature ")+std::to_string(TK)+" or pressure "+
               std::to_string(PPa)+" out of range, or no T/D data are provided" );
        return;
    }

    pm.T = pm.Tc = TK;
    pm.TC = pm.TCc = TK-C_to_K;
    pm.Pc = P;
    if( P < 1e-5 )
    { // Pressure at saturated H2O vapour at given temperature
        long int xT = aNa->check_grid_T(TK);
        if(xT>= 0)
            P = dCH->Psat[xT]/bar_to_Pa;
        else
            P =  LagranInterp( &PPa, dCH->TKval, dCH->Psat, PPa, TK, dCH->nTp, 1,6 )/bar_to_Pa;
    }
    pm.P = P;
    pm.RT = R_CONSTANT * pm.Tc;
    pm.FRT = F_CONSTANT/pm.RT;
    pm.lnP = log( P );

    xTP = aNa->check_grid_TP( TK, PPa );

    for( k=0; k<5; k++ )
    {
        jj =  k * aNa->gridTP();
        if( xTP >= 0 )
        {
            pm.denW[k] = dCH->denW[jj+xTP]/1e3;
            pm.epsW[k] = dCH->epsW[jj+xTP];
            pm.denWg[k] = dCH->denWg[jj+xTP]/1e3;
            pm.epsWg[k] = dCH->epsWg[jj+xTP];
        }
        else
        {
            pm.denW[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->denW+jj,
                                       PPa, TK, dCH->nTp, dCH->nPp,6 )/1e3;// from test denW enough
            pm.epsW[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->epsW+jj,
                                       PPa, TK, dCH->nTp, dCH->nPp,5 );// from test epsW enough
            pm.denWg[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->denWg+jj,
                                        PPa, TK, dCH->nTp, dCH->nPp,5 )/1e3;
            pm.epsWg[k] = LagranInterp( dCH->Pval, dCH->TKval, dCH->epsWg+jj,
                                        PPa, TK, dCH->nTp, dCH->nPp,5 );
        }
    }

#ifdef  USE_THERMO_LOG
    std::fstream f_log("thermodynamic-log-lookup.csv", std::ios::out/*|std::ios::app*/ );
    f_log << "\nCalc ThermoEngine;T;" << TK << ";P;" << PPa << "\n";
    f_log << "denW";
    for( jj=0; jj<5; jj++)
       f_log << ";" << floating_point_to_string(pm.denW[jj]);
    f_log << "\nepsW";
    for( jj=0; jj<5; jj++)
       f_log << ";" << floating_point_to_string(pm.epsW[jj]);
    f_log << "\ndenWg";
    for( jj=0; jj<5; jj++)
       f_log << ";" << floating_point_to_string(pm.denWg[jj]);
    f_log << "\nepsWg";
    for( jj=0; jj<5; jj++)
       f_log << ";" << floating_point_to_string(pm.epsWg[jj]);
#endif
    long int xVol =  getXvolume();

    for( k=0; k<pm.FI; k++ )
    {
        jb = je;
        je += pm.L1[k];
        // load t/d data from DC - to be extended for DCH->H0, DCH->S0, DCH->Cp0, DCH->DD
        // depending on the presence of these arrays in DATACH and Multi structures
        for( j=jb; j<je; j++ )
        {
            jj =  j * aNa->gridTP();
            if( xTP >= 0 )
            {
                Go = dCH->G0[ jj+xTP];
                Vv = dCH->V0[ jj+xTP]*1e5;
                if( dCH->S0 ) S0 = dCH->S0[ jj+xTP];
                if( dCH->H0 ) h0 = dCH->H0[ jj+xTP];
                if( dCH->Cp0 ) Cp0 = dCH->Cp0[ jj+xTP];
                if( dCH->A0 ) a0 = dCH->A0[ jj+xTP];
                if( dCH->U0 ) h0 = dCH->U0[ jj+xTP];
            }
            else
            {
                Go = LagranInterp( dCH->Pval, dCH->TKval, dCH->G0+jj,
                                   PPa, TK, dCH->nTp, dCH->nPp, 6 ); // from test G0[Ca+2] enough
                Vv = LagranInterp( dCH->Pval, dCH->TKval, dCH->V0+jj,
                                   PPa, TK, dCH->nTp, dCH->nPp, 5 )*1e5;
                if( dCH->S0 ) S0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->S0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp, 4 ); // from test S0[Ca+2] enough
                if( dCH->H0 ) h0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->H0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp,5 );
                if( dCH->Cp0 ) Cp0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->Cp0+jj,
                                                    PPa, TK, dCH->nTp, dCH->nPp, 3 ); // from test Cp0[Ca+2] not more
                if( dCH->A0 ) a0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->A0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp,5 );
                if( dCH->U0 ) u0 =  LagranInterp( dCH->Pval, dCH->TKval, dCH->U0+jj,
                                                  PPa, TK, dCH->nTp, dCH->nPp,5 );
            }
#ifndef IPMGEMPLUGIN
            if( TSyst::sm->GetSY()->Guns )  // This is used mainly in UnSpace calculations
                Gg = TSyst::sm->GetSY()->Guns[pm.muj[j]];    // User-set increment to G0 from project system
            // SDGEX     if( syp->GEX && syp->PGEX != S_OFF )   // User-set increment to G0 from project system
            //            Ge = syp->GEX[pm.muj[j]];     //now Ge is integrated into pm.G0 (since 07.03.2008) DK
#else
            if( pm.tpp_G )
                pm.tpp_G[j] = Go;
            if( pm.Guns )
                Gg = pm.Guns[j];
            else
                Gg = 0.;

            Ge = 0.;
#endif
            pm.G0[j] = ConvertGj_toUniformStandardState( Go+Gg+Ge, j, k ); // formerly Cj_init_calc()
            // Inside this function, pm.YOF[k] can be added!

#ifndef IPMGEMPLUGIN
            if( TMTparm::sm->GetTP()->PtvVm != S_ON )
                pm.Vol[j] = 0.;
            else
#endif
                switch( pm.PV )
                { // put molar volumes of components into A matrix or into the vector of molar volumes
                // to be checked!
                case VOL_CONSTR:
#ifndef IPMGEMPLUGIN
                    if( TSyst::sm->GetSY()->Vuns )
                        Vv += TSyst::sm->GetSY()->Vuns[j];
#else
                    if( pm.Vuns )
                        Vv += pm.Vuns[j];
#endif
                    if( xVol >= 0. )
                        pm.A[j*pm.N+xVol] = Vv;
                    // [[fallthrough]];
                case VOL_CALC:
                case VOL_UNDEF:
#ifndef IPMGEMPLUGIN
                    if( TSyst::sm->GetSY()->Vuns )
                        Vv += TSyst::sm->GetSY()->Vuns[j];
#else
                    if( pm.tpp_Vm )
                        pm.tpp_Vm[j] = Vv;
                    if( pm.Vuns )
                        Vv += pm.Vuns[j];
#endif
                    pm.Vol[j] = Vv  * 10.;
                    break;
                }
            if( pm.S0 ) pm.S0[j] = S0;
            if( pm.H0 ) pm.H0[j] = h0;
            if( pm.Cp0 ) pm.Cp0[j] = Cp0;
            if( pm.A0 ) pm.A0[j] = a0;
            if( pm.U0 ) pm.U0[j] = u0;

#ifdef  USE_THERMO_LOG
            f_log << "\n" << std::string(dCH->DCNL[j], 0, MaxDCN) << ";" << floating_point_to_string(Go)
                   << ";" << floating_point_to_string(pm.G0[j])
                   << ";" << floating_point_to_string(pm.Vol[j]);
            if( dCH->S0 ) f_log << ";" << floating_point_to_string(pm.S0[j]);
            if( dCH->H0 ) f_log << ";" << floating_point_to_string(pm.H0[j]);
            if( dCH->Cp0 ) f_log << ";" << floating_point_to_string(pm.Cp0[j]);
            if( dCH->A0 ) f_log << ";" << floating_point_to_string(pm.A0[j]);
            if( dCH->U0 ) f_log << ";" << floating_point_to_string(pm.U0[j]);
#endif
        }  // j
    } // k
}


//--------------------- End of ms_multi_diff.cpp ---------------------------

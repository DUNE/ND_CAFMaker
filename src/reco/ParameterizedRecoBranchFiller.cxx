/// \file ParameterizedBranchFiller.cxx
///
/// Fill reco branches using parameterized reconstruction.
///

#include "ParameterizedRecoBranchFiller.h"

#include "TF1.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "dumpTree.h"
#include "Params.h"


namespace cafmaker
{
  const double mmu = 0.1056583745;

  /// angular resolution function
  // LAr driven smearing, maybe we want to change for gas?
  const TF1 tsmear( "tsmear", "0.162 + 3.407*pow(x,-1.) + 3.129*pow(x,-0.5)", 0., 999.9 );


  void ParameterizedRecoBranchFiller::_FillRecoBranches(caf::StandardRecord & sr, const cafmaker::dumpTree & dt, const params & par) const
  {
    //--------------------------------------------------------------------------
    // Parameterized reconstruction
    //--------------------------------------------------------------------------
    if( !par.IsGasTPC ) {
      // Loop over final-state particles
      double longest_mip = 0.;
      double longest_mip_KE = 0.;
      int longest_mip_charge = 0;
      sr.reco_lepton_pdg = 0;
      int electrons = 0;
      double electron_energy = 0.;
      int reco_electron_pdg = 0;
      for( int i = 0; i < dt.nFS; ++i ) {
        int pdg = dt.fsPdg[i];
        double p = sqrt(dt.fsPx[i]*dt.fsPx[i] + dt.fsPy[i]*dt.fsPy[i] + dt.fsPz[i]*dt.fsPz[i]);
        double KE = dt.fsE[i] - sqrt(dt.fsE[i]*dt.fsE[i] - p*p);

        if( (abs(pdg) == 13 || abs(pdg) == 211) && dt.fsTrkLen[i] > longest_mip ) {
          longest_mip = dt.fsTrkLen[i];
          longest_mip_KE = KE;
          sr.reco_lepton_pdg = pdg;
          if( pdg == 13 || pdg == -211 ) longest_mip_charge = -1;
          else longest_mip_charge = 1;
        }

        // pi0 as nu_e
        if( pdg == 111 ) {
          TVector3 g1, g2;
          TLorentzVector pi0( dt.fsPx[i], dt.fsPy[i], dt.fsPz[i], dt.fsE[i] );
          decayPi0( pi0, g1, g2 );
          double g1conv = fRando->Exp( 14. ); // conversion distance
          bool compton = (fRando->Rndm() < 0.15); // dE/dX misID probability for photon
          // if energetic gamma converts in first wire, and other gamma is either too soft or too colinear
          if( g1conv < 2.0 && compton && (g2.Mag() < 50. || g1.Angle(g2) < 0.01) ) electrons++;
          electron_energy = g1.Mag();
          reco_electron_pdg = 111;
        }
      }

      // True CC reconstruction
      if( abs(dt.lepPdg) == 11 ) { // true nu_e
        recoElectron( sr, par );
        electrons++;
        reco_electron_pdg = dt.lepPdg;
      } else if( abs(dt.lepPdg) == 13 ) { // true nu_mu
        if     ( dt.muGArLen > 50. ) recoMuonTracker( sr, par ); // gas TPC match
        else if( dt.muonReco == 1 ) recoMuonLAr( sr, par ); // LAr-contained muon, this might get updated to NC...
        else if( dt.muonReco == 3 && dt.muECalLen > 5. ) recoMuonECAL( sr, par ); // ECAL-stopper
        else { // exiting but poorly-reconstructed muon
          sr.Elep_reco = longest_mip * 0.0022;
          sr.reco_q = 0;
          sr.reco_numu = 1; sr.reco_nue = 0; sr.reco_nc = 0;
          sr.muon_contained = 0; sr.muon_tracker = 0; sr.muon_ecal = 0; sr.muon_exit = 1;

          double true_tx = 1000.*atan(sr.LepMomX / sr.LepMomZ);
          double true_ty = 1000.*atan(sr.LepMomY / sr.LepMomZ);
          double evalTsmear = tsmear.Eval(sr.Elep_reco - mmu);
          if( evalTsmear < 0. ) evalTsmear = 0.;
          double reco_tx = true_tx + fRando->Gaus(0., evalTsmear/sqrt(2.));
          double reco_ty = true_ty + fRando->Gaus(0., evalTsmear/sqrt(2.));
          sr.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );
        }
      } else { // NC -- set PID variables, will get updated later if fake CC
        sr.Elep_reco = 0.;
        sr.reco_q = 0;
        sr.reco_numu = 0; sr.reco_nue = 0; sr.reco_nc = 1;
        sr.muon_contained = 0; sr.muon_tracker = 0; sr.muon_ecal = 0; sr.muon_exit = 0;
      }

      // CC/NC confusion
      if( electrons == 1 && dt.muonReco <= 1 ) { // NC or numuCC reco as nueCC
        sr.Elep_reco = electron_energy*0.001;
        sr.reco_q = 0;
        sr.reco_numu = 0; sr.reco_nue = 1; sr.reco_nc = 0;
        sr.muon_contained = 0; sr.muon_tracker = 0; sr.muon_ecal = 0; sr.muon_exit = 0;
        sr.reco_lepton_pdg = reco_electron_pdg;
      } else if( dt.muonReco <= 1 && !(abs(dt.lepPdg) == 11 && sr.Elep_reco > 0.) && (longest_mip < par.CC_trk_length || longest_mip_KE/longest_mip > 3.) ) {
        // reco as NC
        sr.Elep_reco = 0.;
        sr.reco_q = 0;
        sr.reco_numu = 0; sr.reco_nue = 0; sr.reco_nc = 1;
        sr.muon_contained = 0; sr.muon_tracker = 0; sr.muon_ecal = 0; sr.muon_exit = 0;
        sr.reco_lepton_pdg = 0;
      } else if( (abs(dt.lepPdg) == 12 || abs(dt.lepPdg) == 14) && longest_mip > par.CC_trk_length && longest_mip_KE/longest_mip < 3. ) { // true NC reco as CC numu
        sr.Elep_reco = longest_mip_KE*0.001 + mmu;
        if( par.fhc ) sr.reco_q = -1;
        else {
          double michel = fRando->Rndm();
          if( longest_mip_charge == 1 && michel < par.michelEff ) sr.reco_q = 1; // correct mu+
          else if( michel < par.michelEff*0.25 ) sr.reco_q = 1; // incorrect mu-
          else sr.reco_q = -1; // no reco Michel
        }
        sr.reco_numu = 1; sr.reco_nue = 0; sr.reco_nc = 0;
        sr.muon_contained = 1; sr.muon_tracker = 0; sr.muon_ecal = 0; sr.muon_exit = 0;
      }

      // Hadronic energy calorimetrically
      sr.Ev_reco = sr.Elep_reco + dt.hadTot*0.001;
      sr.Ehad_veto = dt.hadCollar;
      sr.eRecoP = dt.hadP*0.001;
      sr.eRecoN = dt.hadN*0.001;
      sr.eRecoPip = dt.hadPip*0.001;
      sr.eRecoPim = dt.hadPim*0.001;
      sr.eRecoPi0 = dt.hadPi0*0.001;
      sr.eRecoOther = dt.hadOther*0.001;

      sr.pileup_energy = 0.;
      if( fRando->Rndm() < par.pileup_frac ) sr.pileup_energy = fRando->Rndm() * par.pileup_max;
      sr.Ev_reco += sr.pileup_energy;
    } else {
      // gas TPC: FS particle loop look for long enough tracks and smear momenta
      sr.Ev_reco = 0.;
      sr.nFSP = dt.nFS;
      for( int i = 0; i < dt.nFS; ++i ) {
        double ptrue = 0.001*sqrt(dt.fsPx[i]*dt.fsPx[i] + dt.fsPy[i]*dt.fsPy[i] + dt.fsPz[i]*dt.fsPz[i]);
        double mass = 0.001*sqrt(dt.fsE[i]*dt.fsE[i] - dt.fsPx[i]*dt.fsPx[i] - dt.fsPy[i]*dt.fsPy[i] - dt.fsPz[i]*dt.fsPz[i]);
        sr.pdg[i] = dt.fsPdg[i];
        sr.ptrue[i] = ptrue;
        sr.trkLen[i] = dt.fsTrkLen[i];
        sr.trkLenPerp[i] = dt.fsTrkLenPerp[i];
        // track length cut 6cm according to T Junk
        if( dt.fsTrkLen[i] > 0. && dt.fsPdg[i] != 2112 ) { // basically select charged particles; somehow neutrons ocasionally get nonzero track length
          double pT = 0.001*sqrt(dt.fsPy[i]*dt.fsPy[i] + dt.fsPz[i]*dt.fsPz[i]); // transverse to B field, in GeV
          double nHits = dt.fsTrkLen[i] / par.gastpc_padPitch; // doesn't matter if not integer as only used in eq
          // Gluckstern formula, sigmapT/pT, with sigmaX and L in meters
          double fracSig_meas = sqrt(720./(nHits+4)) * (0.01*par.gastpc_padPitch/sqrt(12.)) * pT / (0.3 * par.gastpc_B * 0.0001 * dt.fsTrkLenPerp[i]*dt.fsTrkLenPerp[i]);
          // multiple scattering term
          double fracSig_MCS = 0.052 / (par.gastpc_B * sqrt(par.gastpc_X0*dt.fsTrkLenPerp[i]*0.0001));

          double sigmaP = ptrue * sqrt( fracSig_meas*fracSig_meas + fracSig_MCS*fracSig_MCS );
          double preco = fRando->Gaus( ptrue, sigmaP );
          double ereco = sqrt( preco*preco + mass*mass ) - mass; // kinetic energy
          if( abs(dt.fsPdg[i]) == 211 ) ereco += mass; // add pion mass
          else if( dt.fsPdg[i] == 2212 && preco > 1.5 ) ereco += 0.1395; // mistake pion mass for high-energy proton
          sr.partEvReco[i] = ereco;

          // threshold cut
          if( dt.fsTrkLen[i] > par.gastpc_len ) {
            sr.Ev_reco += ereco;
            if( dt.fsPdg[i] == 211 || (dt.fsPdg[i] == 2212 && preco > 1.5) ) sr.gastpc_pi_pl_mult++;
            else if( dt.fsPdg[i] == -211 ) sr.gastpc_pi_min_mult++;
          }

          if( (dt.fsPdg[i] == 13 || dt.fsPdg[i] == -13) && dt.fsTrkLen[i] > 100. ) { // muon, don't really care about nu_e CC for now
            sr.partEvReco[i] += mass;
            sr.Elep_reco = sqrt(preco*preco + mass*mass);
            // angle reconstruction
            double true_tx = 1000.*atan(sr.LepMomX / sr.LepMomZ);
            double true_ty = 1000.*atan(sr.LepMomY / sr.LepMomZ);
            double evalTsmear = tsmear.Eval(sr.Elep_reco - mmu);
            if( evalTsmear < 0. ) evalTsmear = 0.;
            double reco_tx = true_tx + fRando->Gaus(0., evalTsmear/sqrt(2.));
            double reco_ty = true_ty + fRando->Gaus(0., evalTsmear/sqrt(2.));
            sr.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );
            // assume perfect charge reconstruction
            sr.reco_q = (dt.fsPdg[i] > 0 ? -1 : 1);
            sr.reco_numu = 1; sr.reco_nue = 0; sr.reco_nc = 0;
            sr.muon_tracker = 1;
          }
        } else if( dt.fsPdg[i] == 111 || dt.fsPdg[i] == 22 ) {
          double ereco = 0.001 * fRando->Gaus( dt.fsE[i], 0.1*dt.fsE[i] );
          sr.partEvReco[i] = ereco;
          sr.Ev_reco += ereco;
        }
      }
    }

  }


  /// Fill reco variables for muon reconstructed in magnetized tracker
  void ParameterizedRecoBranchFiller::recoMuonTracker(caf::StandardRecord &sr, const params &par ) const
  {
    // smear momentum by resolution
    double p = sqrt(sr.LepE*sr.LepE - mmu*mmu);
    double reco_p = fRando->Gaus( p, p*par.trk_muRes );
    sr.Elep_reco = sqrt(reco_p*reco_p + mmu*mmu);

    double true_tx = 1000.*atan(sr.LepMomX / sr.LepMomZ);
    double true_ty = 1000.*atan(sr.LepMomY / sr.LepMomZ);
    double evalTsmear = tsmear.Eval(sr.Elep_reco - mmu);
    if( evalTsmear < 0. ) evalTsmear = 0.;
    double reco_tx = true_tx + fRando->Gaus(0., evalTsmear/sqrt(2.));
    double reco_ty = true_ty + fRando->Gaus(0., evalTsmear/sqrt(2.));
    sr.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );

    // assume perfect charge reconstruction
    sr.reco_q = (sr.LepPDG > 0 ? -1 : 1);

    // assume always muon for tracker-matched
    sr.reco_numu = 1; sr.reco_nue = 0; sr.reco_nc = 0;
    sr.muon_contained = 0; sr.muon_tracker = 1; sr.muon_ecal = 0; sr.muon_exit = 0;
    sr.Ev_reco = sr.Elep_reco;
  }

  /// Fill reco muon variables for muon contained in LAr
  void ParameterizedRecoBranchFiller::recoMuonLAr(caf::StandardRecord &sr, const params &par ) const
  {
    // range-based, smear kinetic energy
    double ke = sr.LepE - mmu;
    double reco_ke = fRando->Gaus( ke, ke*par.LAr_muRes );
    sr.Elep_reco = reco_ke + mmu;

    double true_tx = 1000.*atan(sr.LepMomX / sr.LepMomZ);
    double true_ty = 1000.*atan(sr.LepMomY / sr.LepMomZ);
    double evalTsmear = tsmear.Eval(sr.Elep_reco - mmu);
    if( evalTsmear < 0. ) evalTsmear = 0.;
    double reco_tx = true_tx + fRando->Gaus(0., evalTsmear/sqrt(2.));
    double reco_ty = true_ty + fRando->Gaus(0., evalTsmear/sqrt(2.));
    sr.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );


    // assume negative for FHC, require Michel for RHC
    if( par.fhc ) sr.reco_q = -1;
    else {
      double michel = fRando->Rndm();
      if( sr.LepPDG == -13 && michel < par.michelEff ) sr.reco_q = 1; // correct mu+
      else if( sr.LepPDG == 13 && michel < par.michelEff*0.25 ) sr.reco_q = 1; // incorrect mu-
      else sr.reco_q = -1; // no reco Michel
    }

    sr.reco_numu = 1; sr.reco_nue = 0; sr.reco_nc = 0;
    sr.muon_contained = 1; sr.muon_tracker = 0; sr.muon_ecal = 0; sr.muon_exit = 0;
    sr.Ev_reco = sr.Elep_reco;
  }

  /// Fill reco variables for muon reconstructed in magnetized tracker
  void ParameterizedRecoBranchFiller::recoMuonECAL(caf::StandardRecord & sr, const params &par ) const
  {
    // range-based KE
    double ke = sr.LepE - mmu;
    double reco_ke = fRando->Gaus( ke, ke*par.ECAL_muRes );
    sr.Elep_reco = reco_ke + mmu;

    double true_tx = 1000.*atan(sr.LepMomX / sr.LepMomZ);
    double true_ty = 1000.*atan(sr.LepMomY / sr.LepMomZ);
    double evalTsmear = tsmear.Eval(sr.Elep_reco - mmu);
    if( evalTsmear < 0. ) evalTsmear = 0.;
    double reco_tx = true_tx + fRando->Gaus(0., evalTsmear/sqrt(2.));
    double reco_ty = true_ty + fRando->Gaus(0., evalTsmear/sqrt(2.));
    sr.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );

    // assume perfect charge reconstruction -- these are fairly soft and should curve a lot in short distance
    sr.reco_q = (sr.LepPDG > 0 ? -1 : 1);

    // assume always muon for ecal-matched
    sr.reco_numu = 1; sr.reco_nue = 0; sr.reco_nc = 0;
    sr.muon_contained = 0; sr.muon_tracker = 0; sr.muon_ecal = 1; sr.muon_exit = 0;
    sr.Ev_reco = sr.Elep_reco;
  }

  /// Fill reco variables for true electron
  void ParameterizedRecoBranchFiller::recoElectron(caf::StandardRecord & sr, const params &par ) const
  {
    sr.reco_q = 0; // never know charge
    sr.reco_numu = 0;
    sr.muon_contained = 0; sr.muon_tracker = 1; sr.muon_ecal = 0; sr.muon_exit = 0;

    // fake efficiency...threshold of 300 MeV, eff rising to 100% by 700 MeV
    if( fRando->Rndm() > (sr.LepE-0.3)*2.5 ) { // reco as NC
      sr.Elep_reco = 0.;
      sr.reco_nue = 0; sr.reco_nc = 1;
      sr.Ev_reco = sr.LepE; // include electron energy in Ev anyway, since it won't show up in reco hadronic energy
    } else { // reco as CC
      sr.Elep_reco = fRando->Gaus( sr.LepE, sr.LepE*(par.em_const + par.em_sqrtE/sqrt(sr.LepE)) );
      sr.reco_nue = 1; sr.reco_nc = 0;
      sr.Ev_reco = sr.Elep_reco;
    }

    double true_tx = 1000.*atan(sr.LepMomX / sr.LepMomZ);
    double true_ty = 1000.*atan(sr.LepMomY / sr.LepMomZ);
    double evalTsmear = 3. + tsmear.Eval(sr.Elep_reco - mmu);
    if( evalTsmear < 0. ) evalTsmear = 0.;
    double reco_tx = true_tx + fRando->Gaus(0., evalTsmear/sqrt(2.));
    double reco_ty = true_ty + fRando->Gaus(0., evalTsmear/sqrt(2.));
    sr.theta_reco = 0.001*sqrt( reco_tx*reco_tx + reco_ty*reco_ty );

  }

  void ParameterizedRecoBranchFiller::decayPi0( const TLorentzVector & pi0, TVector3 &gamma1, TVector3 &gamma2 ) const
  {
    double e = pi0.E();
    double mp = 134.9766; // pi0 mass

    double beta = sqrt( 1. - (mp*mp)/(e*e) ); // velocity of pi0
    double theta = 3.1416*fRando->Rndm(); // theta of gamma1 w.r.t. pi0 direction
    double phi = 2.*3.1416*fRando->Rndm(); // phi of gamma1 w.r.t. pi0 direction

    double p = mp/2.; // photon momentum in pi0 rest frame
    TLorentzVector g1( 0., 0., p, p ); // pre-rotation photon 1
    TLorentzVector g2( 0., 0., -p, p ); // pre-rotation photon 2 is opposite

    // rotate to the fRandom decay axis in pi0 rest frame. choice of rotation about x instead of y is arbitrary
    g1.RotateX( theta );
    g2.RotateX( theta );
    g1.RotateZ( phi );
    g2.RotateZ( phi );

    // boost to lab frame with pi0 velocity. pi0 direction is z axis for this
    g1.Boost( 0., 0., beta );
    g2.Boost( 0., 0., beta );

    // make gamma1 the more energetic one
    if( g1.E() > g2.E() ) {
      gamma1 = g1.Vect();
      gamma2 = g2.Vect();
    } else {
      gamma1 = g2.Vect();
      gamma2 = g1.Vect();
    }

    // rotate from frame where pi0 is z' direction into neutrino frame
    TVector3 pi0dir = pi0.Vect().Unit(); // actually w.r.t. neutrino direction
    gamma1.RotateUz( pi0dir );
    gamma2.RotateUz( pi0dir );
  }

}
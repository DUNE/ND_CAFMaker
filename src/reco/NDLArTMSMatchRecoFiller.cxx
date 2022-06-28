#include "NDLArTMSMatchRecoFiller.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRNDTrackAssn.h"
#include "MLNDLArRecoBranchFiller.h"
#include "TMSRecoBranchFiller.h"


namespace cafmaker
{
  NDLArTMSMatchRecoFiller::NDLArTMSMatchRecoFiller()
  {
    // nothing to do
    SetConfigured(true);
    name = "LArTMSMatcher";
  }

  void NDLArTMSMatchRecoFiller::_FillRecoBranches(std::size_t evtIdx,
                                                  caf::StandardRecord &sr,
                                                  const cafmaker::dumpTree &dt,
                                                  const cafmaker::Params &par) const
  {
    // match tracks using the info that should have been filled by the ND-LAr and TMS reco fillers
    //MatchTracks(sr);
    unsigned int n_TMS_tracks = sr.nd.tms.ntracks;  //# tracks in TMS in evtIdx
    unsigned int n_LAr_tracks = sr.nd.lar.ntracks;  //# tracks in LAr in evtIdx

    sr.nd.trkmatch.resize(n_LAr_tracks);

    for(unsigned int ilar = 0; ilar < n_LAr_tracks; ++ilar )
    {
      x1_lar = 10.*sr.nd.lar.tracks[ilar].start.x;   //mm
      z1_lar = 10.*sr.nd.lar.tracks[ilar].start.z;   //mm
      x2_lar = 10.*sr.nd.lar.tracks[ilar].end.x;     //mm
      z2_lar = 10.*sr.nd.lar.tracks[ilar].end.z;     //mm

      if(z2_lar > 9000.0)  //check: LAr track exits through the end
      {
        for(unsigned int itms = 0; itms < n_TMS_tracks; ++itms )
        {
          x1_tms = sr.nd.tms.tracks[itms].start.x;   //mm
          z1_tms = sr.nd.tms.tracks[itms].start.z;   //mm
          x2_tms = sr.nd.tms.tracks[itms].end.x;     //mm
          z2_tms = sr.nd.tms.tracks[itms].end.z;     //mm
          
          if(z1_tms >= 11000.0 && z1_tms <12000.0) //check: TMS track starts in first few planes.
          {
            std::cout<<"ND_tracks = "<<x1_lar<<"\t"<<z1_lar<<"\t"<<x2_lar<<"\t"<<z2_lar<<std::endl;
            std::cout<<"TMS_tracks = "<<x1_tms<<"\t"<<z1_tms<<"\t"<<x2_tms<<"\t"<<z2_tms<<std::endl;

            //Draw a straight line, extrapolate it to Z = 1000.0 cm and calculate slope
            TF1 *LAr = DrawLines(x1_lar, z1_lar, x2_lar, z2_lar, z1_lar, z2_lar);
            x_LAr = LAr->Eval(10000.);
            slope_LAr = (x2_lar - x1_lar)/(z2_lar - z1_lar);

            //Draw a straight line, extrapolate it to Z = 1000.0 cm and calculate slope
            TF1 *TMS_ex = DrawLines(x1_tms, z1_tms, x2_tms, z2_tms, z1_tms, z2_tms);
            x_TMS = TMS_ex->Eval(10000.);
            slope_TMS = (x2_tms - x1_tms)/(z2_tms - z1_tms);

            slope_TMS_LAr = abs((slope_TMS - slope_LAr)/(1.0 + slope_TMS * slope_LAr));  //tan(theta) = |(m1 - m1)/(1 + m1*m2)|
            ang = TMath::ATan(slope_TMS_LAr);  //theta = tan^-1()
            residual = x_LAr - x_TMS;   //distance between two tracks
            Costheta = cos(ang);        //angle between two tracks
          }
        }
      }
        if(abs(residual) < 400. && Costheta > 0.95){  //matching 
            sr.nd.LArID = evtIdx;
            sr.nd.TMSID = evtIdx;
            sr.nd.Residual = residual;
            sr.nd.cosTheta = Costheta;

            caf::SRNDTrackAssn match;
            match.larid = evtIdx;
            match.tmsid = evtIdx;
            match.transdispl = residual;
            match.angdispl = Costheta;

            sr.nd.trkmatch.emplace_back(match);
            std::cout<<"match track increment = "<<sr.nd.ntrkmatch<<std::endl;
            std::cout<<"match:larid:tmsid = "<<evtIdx<<std::endl;
            std::cout<<"match:residual = "<<residual<<std::endl;
            std::cout<<"match:Costheta = "<<Costheta<<std::endl;
      }
    }
  }
}

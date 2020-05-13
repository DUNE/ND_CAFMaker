#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array

import pyGeoEff

def loop( events, tgeo, tout ):

    offset = [ 0., 5.5, 411. ]
    collarLo = [ -320., -120., 30. ]
    collarHi = [ 320., 120., 470. ]

    # Initialize geometric efficiency module.
    geoEff = pyGeoEff.geoEff(args.seed)
    # Multiple of 64 doesn't waste bits
    geoEff.setNthrows(4992)
    # Use neutrino decay position, rather than fixed neutrino direction as symmetry axis
    geoEff.setUseFixedBeamDir(False)
    # Decay position in detector coordinates. Rough estimate from neutrino direction in mcc11v4 -- might want to update. In cm.
    geoEff.setDecayPos(-args.offaxis*100, -5155, -55400)
    # 30 cm veto
    geoEff.setVetoSizes([30])
    # 20, 30 and 40 MeV threshold
    geoEff.setVetoEnergyThresholds([20, 30, 40])
    # Active detector dimensions
    geoEff.setActiveX(collarLo[0]-30, collarHi[0]+30)
    geoEff.setActiveY(collarLo[1]-30, collarHi[1]+30)
    geoEff.setActiveZ(collarLo[2]-30, collarHi[2]+30)
    # Range for translation throws. Use full active volume but fix X.
    geoEff.setRangeX(-1, -1)
    geoEff.setRandomizeX(False)
    geoEff.setRangeY(collarLo[1]-30, collarHi[1]+30)
    geoEff.setRangeZ(collarLo[2]-30, collarHi[2]+30)
    # Set offset between MC coordinate system and volumes defined above.
    geoEff.setOffsetX(offset[0])
    geoEff.setOffsetY(offset[1])
    geoEff.setOffsetZ(offset[2])

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = events.GetEntries()

    print "Starting loop over %d entries" % N
    ient = 0 # This is unnecessary
    iwritten = 0
    for ient in range(N):

        if ient % 100 == 0:
            print "Event %d of %d..." % (ient,N)
        events.GetEntry(ient)

        for ivtx,vertex in enumerate(event.Primaries):

            ## initialize output variables
            t_ievt[0] = ient;
            t_vtx[0]=0.0; t_vtx[1]=0.0; t_vtx[2]=0.0;
            t_p3lep[0]=0.0; t_p3lep[1]=0.0; t_p3lep[2]=0.0;
            t_lepDeath[0]=0.0; t_lepDeath[1]=0.0; t_lepDeath[2]=0.0;
            t_lepPdg[0] = 0
            t_lepKE[0] = 0.
            t_muonExitPt[0] = 0.0; t_muonExitPt[1] = 0.0; t_muonExitPt[2] = 0.0; 
            t_muonExitMom[0] = 0.0; t_muonExitMom[1] = 0.0; t_muonExitMom[2] = 0.0; 
            t_muonReco[0] = -1;
            t_muon_endVolName.replace(0, ROOT.std.string.npos, "")
            t_muGArLen[0]=0.0;
            t_hadTot[0] = 0.
            t_hadP[0] = 0.
            t_hadN[0] = 0.
            t_hadPip[0] = 0.
            t_hadPim[0] = 0.
            t_hadPi0[0] = 0.
            t_hadOther[0] = 0.
            t_hadCollar[0] = 0.
            t_nFS[0] = 0
            ## done

            # now ID numucc
            reaction=vertex.Reaction

            # set the vertex location for output
            for i in range(3): 
                t_vtx[i] = vertex.Position[i] / 10. - offset[i] # cm

            # fiducial vertex pre-cut
            #if abs(t_vtx[0]) > 310. or abs(t_vtx[1]) > 110. or t_vtx[2] < 40. or t_vtx[2] > 360.:
            #    continue
            

            geoEff.setVertex(vertex.Position[0] / 10., vertex.Position[1] / 10.,vertex.Position[2] / 10.)
            
            # Renew throws every 100th event written to the output file.
            if (iwritten % 100) == 0 :
                geoEff.throwTransforms()
                t_geoEffThrowsY.clear()
                for i in geoEff.getCurrentThrowTranslationsY() :
                    t_geoEffThrowsY.push_back(i)
                t_geoEffThrowsZ.clear()
                for i in geoEff.getCurrentThrowTranslationsZ() :
                    t_geoEffThrowsZ.push_back(i)
                t_geoEffThrowsPhi.clear()
                for i in geoEff.getCurrentThrowRotations() :
                    t_geoEffThrowsPhi.push_back(i)
                tGeoEfficiencyThrowsOut.Fill()
                
            ileptraj = -1
            nfsp = 0
            nHadrons = 0
            # get the lepton kinematics from the edepsim file
            fsParticleIdx = {}
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.Momentum[3]
                p = (particle.Momentum[0]**2 + particle.Momentum[1]**2 + particle.Momentum[2]**2)**0.5
                m = (e**2 - p**2)**0.5
                t_fsPdg[nfsp] = particle.PDGCode
                t_fsPx[nfsp] = particle.Momentum[0]
                t_fsPy[nfsp] = particle.Momentum[1]
                t_fsPz[nfsp] = particle.Momentum[2]
                t_fsE[nfsp] = e
                fsParticleIdx[particle.TrackId] = nfsp
                nfsp += 1
                pdg = particle.PDGCode
                if abs(pdg) in [11,12,13,14]:
                    ileptraj = particle.TrackId
                    t_lepPdg[0] = pdg
                    # set the muon momentum for output
                    for i in range(3): t_p3lep[i] = particle.Momentum[i]
                    t_lepKE[0] = e - m

            assert ileptraj != -1, "There isn't a lepton??"
            t_nFS[0] = nfsp

            # If there is a muon, determine how to reconstruct its momentum and charge
            if abs(t_lepPdg[0]) == 13:
                leptraj = event.Trajectories[ileptraj]
                for p in leptraj.Points:
                    pt = p.Position
                    node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                    volName = node.GetName()
                    active = False
                    if "LAr" in volName or "PixelPlane" in volName or "sPlane" in volName: # in active volume, update exit points
                        t_muonExitPt[0] = pt.X() / 10. - offset[0]
                        t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                        t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                        t_muonExitMom[0] = p.Momentum.x()
                        t_muonExitMom[1] = p.Momentum.y()
                        t_muonExitMom[2] = p.Momentum.z()
                    else:
                        t_muonExitPt[0] = pt.X() / 10. - offset[0]
                        t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                        t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                        break

                endpt = leptraj.Points[-1].Position

                node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

                t_lepDeath[0] = endpt.X()/10. - offset[0]
                t_lepDeath[1] = endpt.Y()/10. - offset[1]
                t_lepDeath[2] = endpt.Z()/10. - offset[2]

                endVolName = node.GetName()
                t_muon_endVolName.replace(0, ROOT.std.string.npos, endVolName)
                if "LArActive" in endVolName: t_muonReco[0] = 1 # contained
                elif "ECal" in endVolName: t_muonReco[0] = 2 # ECAL stopper
                else: t_muonReco[0] = 0 # endpoint not in active material, but might still be reconstructed by curvature if GAr length > 0

                # look for muon hits in the gas TPC
                hits = []
                for key in event.SegmentDetectors:
                    if key.first in ["TPC_Drift1", "TPC_Drift2"]:
                        hits += key.second

                tot_length = 0.0
                for hit in hits:
                    if hit.PrimaryId == ileptraj: # hit is due to the muon
                        # TG4HitSegment::TrackLength includes all delta-rays, which spiral in gas TPC and give ridiculously long tracks
                        hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                        hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                        tot_length += (hStop-hStart).Mag()

                # look for muon hits in the ECAL
                ehits = []
                for key in event.SegmentDetectors:
                    if key.first in ["BarrelECal_vol", "EndcapECal_vol"]:
                        ehits += key.second

                etot_length = 0.0
                for hit in ehits:
                    if hit.PrimaryId == ileptraj: # hit is due to the muon
                        # TG4HitSegment::TrackLength includes all delta-rays, which spiral in gas TPC and give ridiculously long tracks
                        hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                        hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                        etot_length += (hStop-hStart).Mag()
                
                t_muGArLen[0] = tot_length
                t_muECalLen[0] = etot_length

            # hadronic containment -- find hits in ArgonCube
            hits = []
            for key in event.SegmentDetectors:
                if key.first == "ArgonCube":
                    hits += key.second

            # Truth-matching energy -- make dictionary of trajectory --> primary pdg
            traj_to_pdg = {}
            # For pi0s, stop at the photons to determine the visible energy due to gamma1/gamma2
            tid_to_gamma = {}
            gamma_tids = []
            for traj in event.Trajectories:
                mom = traj.ParentId
                tid = traj.TrackId
                if event.Trajectories[mom].PDGCode == 111 and event.Trajectories[tid].PDGCode == 22 and event.Trajectories[mom].ParentId == -1:
                    gamma_tids.append(tid)
                while mom != -1:
                    tid = mom
                    mom = event.Trajectories[mom].ParentId
                    if mom in gamma_tids:
                        tid_to_gamma[tid] = mom
                traj_to_pdg[traj] = event.Trajectories[tid].PDGCode

            collar_energy = 0.
            total_energy = 0.
            
            geoEff_EDepPosition = []
            geoEff_EDepEnergy = []

            track_length = [0. for i in range(nfsp)]
            dEdX = [[] for i in range(nfsp)]
            this_step = [[0.,0.] for i in range(nfsp)]
            end_point = [None for i in range(nfsp)]
            int_energy = [0. for i in range(nfsp)]
            trk_calo = [0. for i in range(nfsp)]
            gamma_energy = {}
            for g in gamma_tids:
                gamma_energy[g] = 0.
            for hit in hits:
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )

                # Don't use edep-sim's PrimaryId, which thinks you want to associate absoltely everything with the primary
                # Instead, get the actual contributors (usually only one) and take the biggest
                tid = hit.Contrib[0]

                traj = event.Trajectories[tid]
                if traj.ParentId == -1: # primary particle
                    idx = fsParticleIdx[hit.PrimaryId]
                    trk_calo[idx] += hit.EnergyDeposit
                    end_point[idx] = hStop
                    dx = (hStop-hStart).Mag()
                    track_length[idx] += dx
                    this_step[idx][1] += dx
                    this_step[idx][0] += hit.EnergyDeposit
                    if this_step[idx][1] > 0.5:
                        dEdX[idx].append( (this_step[idx][0], this_step[idx][1], hStart) ) # MeV/cm
                        this_step[idx] = [0., 0.]
                else: # non-primary energy
                    for k,ep in enumerate(end_point):
                        if ep is None: continue
                        if (hStart-ep).Mag() < 10.:
                            int_energy[k] += hit.EnergyDeposit

                if tid in tid_to_gamma:
                    gamma_energy[tid_to_gamma[tid]] += hit.EnergyDeposit

                if hit.PrimaryId != ileptraj: # here we do want to associate stuff to the lepton
                    hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                    total_energy += hit.EnergyDeposit
                    # check if hit is in collar region
                    if hStart.x() < collarLo[0] or hStart.x() > collarHi[0] or hStart.y() < collarLo[1] or hStart.y() > collarHi[1] or hStart.z() < collarLo[2] or hStart.z() > collarHi[2]:
                        collar_energy += hit.EnergyDeposit
                    
                    # Set up arrays for geometric efficiency
                    for dim in range(3) :
                        geoEff_EDepPosition.append((hit.Start[dim] + hit.Stop[dim])/2./10.)
                    geoEff_EDepEnergy.append(hit.EnergyDeposit)

                    # Determine primary particle
                    pdg = traj_to_pdg[traj]
                    if pdg in [11, -11, 13, -13]: continue # lepton
                    elif pdg == 2212: t_hadP[0] += hit.EnergyDeposit
                    elif pdg == 2112: t_hadN[0] += hit.EnergyDeposit
                    elif pdg == 211: t_hadPip[0] += hit.EnergyDeposit
                    elif pdg == -211: t_hadPim[0] += hit.EnergyDeposit
                    elif pdg == 111: t_hadPi0[0] += hit.EnergyDeposit
                    else: t_hadOther[0] += hit.EnergyDeposit

            t_hadTot[0] = total_energy
            t_hadCollar[0] = collar_energy

            geoEff.setHitSegEdeps(geoEff_EDepEnergy)
            geoEff.setHitSegPoss(geoEff_EDepPosition)

            geoEffThrowResultsList = geoEff.getHadronContainmentThrows()
            
            t_geoEffThrowResults.clear()
            for i in range(len(geoEffThrowResultsList)) :
                iVec = ROOT.std.vector('vector< uint64_t >')()
                for j in range (len(geoEffThrowResultsList[i])) :
                    jVec = ROOT.std.vector('uint64_t')()
                    for k in range(len(geoEffThrowResultsList[i][j])) :
                        jVec.push_back(geoEffThrowResultsList[i][j][k])
                    iVec.push_back(jVec)
                t_geoEffThrowResults.push_back(iVec)
                                       
            for i in range(nfsp):
                t_fsTrkLen[i] = track_length[i]
                # Average of anything in the last 3cm of track
                t_fsTrkEnddEdX[i] = 0.
                t_fsTrkFrontdEdX[i] = 0.
                totlen = 0.
                j = len(dEdX[i])-1
                while j >= 0:
                    totlen += dEdX[i][j][1]
                    t_fsTrkEnddEdX[i] += dEdX[i][j][0]
                    j -= 1
                    if totlen > 3.: break
                if totlen > 0.: t_fsTrkEnddEdX[i] /= totlen
              
                totlen = 0.
                for k in range(j+1):
                    totlen += dEdX[i][k][1]
                    t_fsTrkFrontdEdX[i] += dEdX[i][k][0]
                if totlen > 0.: t_fsTrkFrontdEdX[i] /= totlen

                t_fsTrkEndpointBall[i] = int_energy[i]
                t_fsTrkCalo[i] = trk_calo[i]

                t_fsGamma1[i] = 0.
                t_fsGamma2[i] = 0.

            # Photons
            for t in gamma_tids:
                mom = event.Trajectories[t].ParentId
                if t_fsGamma1[mom] == 0.: t_fsGamma1[mom] = gamma_energy[t]
                elif t_fsGamma2[mom] == 0.: t_fsGamma2[mom] = gamma_energy[t]
                else:
                    print "Pi0 has more than two photons wtf"


            tout.Fill()
            iwritten += 1
        ient += 1 # This is unnecessary
        

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--infile', help='Input file name', default="edep.root")
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--offaxis', help='Off-axis position in metres', default=0., type = "float")
    parser.add_option('--seed', help='Seed for geometric efficiency throws', default=0, type = "int")

    (args, dummy) = parser.parse_args()

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tout = ROOT.TTree( "tree","tree" )
    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')
    t_Ev = array('f', [0.])
    tout.Branch('Ev',t_Ev,'Ev/F')
    t_p3lep = array('f',3*[0.0])
    tout.Branch('p3lep',t_p3lep,'p3lep[3]/F')
    t_vtx = array('f',3*[0.0])
    tout.Branch('vtx',t_vtx,'vtx[3]/F')
    t_lepDeath = array('f',3*[0.0])
    tout.Branch('lepDeath',t_lepDeath,'lepDeath[3]/F')
    t_lepPdg = array('i',[0])
    tout.Branch('lepPdg',t_lepPdg,'lepPdg/I')
    t_lepKE = array('f',[0])
    tout.Branch('lepKE',t_lepKE,'lepKE/F')
    t_muonExitPt = array('f',3*[0.0])
    tout.Branch('muonExitPt',t_muonExitPt,'muonExitPt[3]/F')
    t_muonExitMom = array('f',3*[0.0])
    tout.Branch('muonExitMom',t_muonExitMom,'muonExitMom[3]/F')
    t_muonReco = array('i',[0])
    tout.Branch('muonReco',t_muonReco,'muonReco/I')
    t_muon_endVolName = ROOT.std.string()
    tout.Branch('muon_endVolName', t_muon_endVolName)
    t_muGArLen = array('f',[0])
    tout.Branch('muGArLen',t_muGArLen,'muGArLen/F')
    t_muECalLen = array('f',[0])
    tout.Branch('muECalLen',t_muECalLen,'muECalLen/F')
    t_hadTot = array('f', [0.] )
    tout.Branch('hadTot', t_hadTot, 'hadTot/F' )
    t_hadP = array('f', [0.] )
    tout.Branch('hadP', t_hadP, 'hadP/F' )
    t_hadN = array('f', [0.] )
    tout.Branch('hadN', t_hadN, 'hadN/F' )
    t_hadPip = array('f', [0.] )
    tout.Branch('hadPip', t_hadPip, 'hadPip/F' )
    t_hadPim = array('f', [0.] )
    tout.Branch('hadPim', t_hadPim, 'hadPim/F' )
    t_hadPi0 = array('f', [0.] )
    tout.Branch('hadPi0', t_hadPi0, 'hadPi0/F' )
    t_hadOther = array('f', [0.] )
    tout.Branch('hadOther', t_hadOther, 'hadOther/F' )
    t_hadCollar = array('f', [0.] )
    tout.Branch('hadCollar', t_hadCollar, 'hadCollar/F' )
    t_nFS = array('i',[0])
    tout.Branch('nFS',t_nFS,'nFS/I')
    t_fsPdg = array('i',100*[0])
    tout.Branch('fsPdg',t_fsPdg,'fsPdg[nFS]/I')
    t_fsPx = array('f',100*[0.])
    tout.Branch('fsPx',t_fsPx,'fsPx[nFS]/F')
    t_fsPy = array('f',100*[0.])
    tout.Branch('fsPy',t_fsPy,'fsPy[nFS]/F')
    t_fsPz = array('f',100*[0.])
    tout.Branch('fsPz',t_fsPz,'fsPz[nFS]/F')
    t_fsE = array('f',100*[0.])
    tout.Branch('fsE',t_fsE,'fsE[nFS]/F')
    t_fsTrkLen = array('f',100*[0.])
    tout.Branch('fsTrkLen',t_fsTrkLen,'fsTrkLen[nFS]/F')
    t_fsTrkFrontdEdX = array('f',100*[0.])
    tout.Branch('fsTrkFrontdEdX',t_fsTrkFrontdEdX,'fsTrkFrontdEdX[nFS]/F')
    t_fsTrkEnddEdX = array('f',100*[0.])
    tout.Branch('fsTrkEnddEdX',t_fsTrkEnddEdX,'fsTrkEnddEdX[nFS]/F')
    t_fsTrkEndpointBall = array('f', 100*[0.])
    tout.Branch('fsTrkEndpointBall',t_fsTrkEndpointBall,'fsTrkEndpointBall[nFS]/F')
    t_fsGamma1 = array('f',100*[0.])
    tout.Branch('fsGamma1',t_fsGamma1,'fsGamma1[nFS]/F')
    t_fsGamma2 = array('f',100*[0.])
    tout.Branch('fsGamma2',t_fsGamma2,'fsGamma2[nFS]/F')
    t_fsTrkCalo = array('f',100*[0.])
    tout.Branch('fsTrkCalo',t_fsTrkCalo,'fsTrkCalo[nFS]/F')

    # Geometric efficiency stuff
    t_geoEffThrowResults = ROOT.std.vector('std::vector< std::vector < uint64_t > >')()
    tout.Branch('geoEffThrowResults', t_geoEffThrowResults)

    # Separate TTree to store translations and rotations of throws
    tGeoEfficiencyThrowsOut = ROOT.TTree( "geoEffThrows","geoEffThrows")
    t_geoEffThrowsY = ROOT.std.vector('float')()
    tGeoEfficiencyThrowsOut.Branch("geoEffThrowsY", t_geoEffThrowsY)
    t_geoEffThrowsZ = ROOT.std.vector('float')()
    tGeoEfficiencyThrowsOut.Branch("geoEffThrowsZ", t_geoEffThrowsZ)
    t_geoEffThrowsPhi = ROOT.std.vector('float')()
    tGeoEfficiencyThrowsOut.Branch("geoEffThrowsPhi", t_geoEffThrowsPhi)

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )
    #dspt = ROOT.TChain( "DetSimPassThru/gRooTracker", "other thing" )

    tf = ROOT.TFile( args.infile )
    tf.MakeProject("EDepSimEvents","*","RECREATE++")

    events = tf.Get( "EDepSimEvents" )
    tgeo = tf.Get("EDepSimGeometry")
    loop( events, tgeo, tout )

    fout.cd()
    tout.Write()
    tGeoEfficiencyThrowsOut.Write()
    
    




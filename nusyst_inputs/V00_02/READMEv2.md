# Purpose

  Inputs for the neutrino interaction systematics provided by V00_02 of the
  nusystematics package. This release is an RC for the TDR LBL oscillation
  analysis.

# Content
  
  Two ISystProvider_tools in use require ROOT inputs.

  MKSinglePiTemplate_tool: Uses EnuQ2W templates for NEUT MK SPP/GENIE R-S SPP
    these are included in 3 root files, one for each of Enu = 1,2,3 GeV-centred 
    templates.

  EbLepMomShift_tool: Uses EnuFSLepCTheta templates to apply final state 
    charged lepton momentum shifts to approximate the effect of changing the
    binding energy in L-S CCQE.

  FSILikeEAvailShift_tool: Uses ERecoil/q0,q0,q3 templates to shift the
  Ar40 calorimetric response factor (ERecoil/q0) without changing the total
  xsec in (q0,q3).

# Longevity

  These inputs are only meant for use with the V00_02 release, the fully 
  validated TDR release will be tagged later and the inputs will be distributed
  via a versioned data ups product.

///
/// \author  L.Di Noto for getting SAND reco evt
///// \date    Apr. 2022

#ifndef ND_CAFMAKER_SANDEVT_H
#define ND_CAFMAKER_SANDEVT_H

#include "fwd.h"
#include "struct.h"
// fixme: this will need to be put back to the actual response_helper type when DIRT-II finishes model recommendations
#include <string>

class SANDEvt{
	private:
	SANDEvt();
	static SANDEvt* me;
        event* fEvent;

	public:
 	virtual ~SANDEvt();
    	static SANDEvt* Get();

	void SetSANDEvt(event * evt){ fEvent=evt;}
	double GetSANDEvtEnureco();
	bool GetSANDEvtLepton(particle &Lepton);
	

};

#endif //ND_CAFMAKER_SANDEVT

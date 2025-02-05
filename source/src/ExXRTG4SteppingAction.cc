// ExXRTPSteppingAction.hh
//
// Giuseppe Vacanti (cosine Science & Computing BV)
// October 14, 2009
//
// This file is part of the XRTG4 Geant4 extension.

//     The XRTG4 Geant4 extension is free software: you can
//     redistribute it and/or modify it under the terms of the GNU
//     General Public License as published by the Free Software
//     Foundation, either version 2 of the License, or (at your
//     option) any later version.

//     XRTG4 Geant4 extension is distributed in the hope that it will
//     be useful, but WITHOUT ANY WARRANTY; without even the implied
//     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//     PURPOSE.  See the GNU General Public License for more details.

//     You should have received a copy of the GNU General Public
//     License along with XRTG4 Geant4 extension.  If not, see
//     <http://www.gnu.org/licenses/>.
//

//#include "Analysis.hh" // for AnalysisManager
#include "G4AnalysisManager.hh"

#include "ExXRTG4SteppingAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4XrayUserEventInformation.hh"
#include "G4XrayUserTrackInformation.hh"

ExXRTG4SteppingAction::ExXRTG4SteppingAction() : G4UserSteppingAction() {}
ExXRTG4SteppingAction::~ExXRTG4SteppingAction() {}
void ExXRTG4SteppingAction::UserSteppingAction(const G4Step* aStep) {

  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  const G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  //G4int eventID = analysisManager->GetCurrentEvent()->GetEventID();

  const G4Track * aTrack = aStep->GetTrack();
  //G4int currentStepNum = aTrack->GetCurrentStepNumber();
  
  G4XrayUserTrackInformation * track_info = dynamic_cast<G4XrayUserTrackInformation *>(aTrack->GetUserInformation());
  const std::vector<G4XrayUserTrackInformation::interaction_data_t> & data = track_info->GetTrackInformation();

  /*
  std::string filename = "test.txt";
  std::ofstream writing_file;
  writing_file.open(filename, std::ios::app);
  */
  const G4DynamicParticle *photon = aTrack->GetDynamicParticle();
  const G4double photonEnergy = photon->GetKineticEnergy();
  
  if(aTrack->GetTrackStatus()!=fAlive){ /// case of track terminated
    G4cout << G4endl;
    G4cout << "Stepping action: "
	   << "\n event ID=" << eventID
	   << "\n track ID=" << aTrack->GetTrackID()
	   << "\n track status=" << aTrack->GetTrackStatus()
	   << "\n track kinEnergy=" << aTrack->GetKineticEnergy()
	   << "\n track position=" << aTrack->GetPosition()
	   << "\n PreStepPoint=" << aStep->GetPreStepPoint()->GetPosition()
	   << "\n PostStepPoint=" << aStep->GetPostStepPoint()->GetPosition()
	   << "\n volume=" << aTrack->GetVolume()->GetName()
	   << "\n track info records so far=" << data.size()
	   << G4endl;

    size_t data_size = data.size();
    for(size_t k = 0; k<data.size(); k++){
      const G4XrayUserTrackInformation::interaction_data_t & idata = data[k];
      analysisManager->FillNtupleIColumn(0, eventID);
      analysisManager->FillNtupleIColumn(1, data_size);
      analysisManager->FillNtupleIColumn(2, k);
      analysisManager->FillNtupleDColumn(3, idata.position[0]);
      analysisManager->FillNtupleDColumn(4, idata.position[1]);
      analysisManager->FillNtupleDColumn(5, idata.position[2]);
      analysisManager->FillNtupleDColumn(6, idata.incidentDirection[0]);
      analysisManager->FillNtupleDColumn(7, idata.incidentDirection[1]);
      analysisManager->FillNtupleDColumn(8, idata.incidentDirection[2]);
      analysisManager->FillNtupleDColumn(9, photonEnergy);
      analysisManager->AddNtupleRow();

      /*
      writing_file << eventID << " " << data_size << " " << k << " " 
		   << idata.position[0] << " " << idata.position[1] << " " << idata.position[2] << " "
		   << idata.incidentDirection[0] << " " <<idata.incidentDirection[1] << " " <<idata.incidentDirection[2]<< " "
		   <<  std::endl;
      */

    }

  }
}


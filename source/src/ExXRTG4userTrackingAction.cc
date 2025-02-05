// ExXRTG4userTrackingAction.cc
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

#include "ExXRTG4UserTrackingAction.hh"

#include "G4XrayUserTrackInformation.hh"
#include "G4Track.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"

ExXRTG4UserTrackingAction::ExXRTG4UserTrackingAction() : G4UserTrackingAction() {}

ExXRTG4UserTrackingAction::~ExXRTG4UserTrackingAction() {}

void ExXRTG4UserTrackingAction::PreUserTrackingAction(const G4Track * aTrack) {
  G4TrackingManager * tman = G4EventManager::GetEventManager()->GetTrackingManager();
  G4XrayUserTrackInformation * uti = new G4XrayUserTrackInformation;
  G4XrayUserTrackInformation::interaction_data_t track_info_data;
  track_info_data.status = G4XrayUserTrackInformation::just_started;
  track_info_data.position = aTrack->GetVertexPosition();
  track_info_data.incidentDirection = aTrack->GetVertexMomentumDirection();
  track_info_data.reflectedDirection = aTrack->GetMomentumDirection();
  track_info_data.intensity = aTrack->GetWeight();
  track_info_data.volume = aTrack->GetVolume();
  uti->AddInformation(track_info_data);
  tman->SetUserTrackInformation(uti);   
 
}

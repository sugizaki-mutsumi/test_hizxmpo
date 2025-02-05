// ExXRTG4userTrackingAction.hh
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

#ifndef EXXRTG4USERTRACKINGACTION_HH
#define EXXRTG4USERTRACKINGACTION_HH 1

#include "G4UserTrackingAction.hh"

class ExXRTG4UserTrackingAction : public G4UserTrackingAction {
public:
  ExXRTG4UserTrackingAction();
  ~ExXRTG4UserTrackingAction();
  void PreUserTrackingAction(const G4Track *);

};


#endif // EXXRTG4USERTRACKINGACTION_HH

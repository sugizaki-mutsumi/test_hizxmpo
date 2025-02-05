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

#ifndef EXXRTG4STEPPINGACTION_HH
#define EXXRTG4STEPPINGACTION_HH 1

#include "G4UserSteppingAction.hh"

class G4Step;

class ExXRTG4SteppingAction : public G4UserSteppingAction {
public:
  ExXRTG4SteppingAction();
  ~ExXRTG4SteppingAction();
  void UserSteppingAction(const G4Step*);
};

#endif // EXXRTG4STEPPINGACTION_HH

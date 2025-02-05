// ExXRTG4PhysicsList.hh
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


#include "ExXRTG4PhysicsList.hh"

#include "G4XrayGrazingAngleScattering.hh"

#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// Added by MS
//#include "G4DecayPhysics.hh"
//#include "G4RadioactiveDecayPhysics.hh"
//#include "G4EmStandardPhysics.hh"

ExXRTG4PhysicsList::ExXRTG4PhysicsList() : G4VModularPhysicsList() {

  // EM physics
  //RegisterPhysics(new G4DecayPhysics());
  //RegisterPhysics(new G4RadioactiveDecayPhysics());
  //RegisterPhysics(new G4EmStandardPhysics());
}

ExXRTG4PhysicsList::~ExXRTG4PhysicsList() {}

void ExXRTG4PhysicsList::ConstructParticle() {

  G4Gamma::GammaDefinition();
}

void ExXRTG4PhysicsList::ConstructProcess() {
  //Em physics
  //G4VPhysicsConstructor* emList = new G4EmStandardPhysics();
  //emList->ConstructProcess();

  
  AddTransportation();
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
  G4ProcessManager* pmanager = particle->GetProcessManager();
  G4String particleName = particle->GetParticleName();
    
  const G4int xray_verbosity = 1;
  const G4double xray_intensity_threshold = 1e-5; // if xray_mode = intesity, photons below this threshold are dropped.
  G4XrayGrazingAngleScattering::photon_mode_t xray_mode = G4XrayGrazingAngleScattering::particle;
  pmanager->AddDiscreteProcess(new G4XrayGrazingAngleScattering(xray_verbosity, xray_intensity_threshold, xray_mode));
  
}

void ExXRTG4PhysicsList::SetCuts() {}

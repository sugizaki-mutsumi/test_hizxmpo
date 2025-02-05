// xrtg4.cc
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


// This is an example showing how to make use of the X-ray grazing
// angle scattering process.
#include "G4SystemOfUnits.hh"

/// User action including primary generator
#include "UserActionInitialization.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

/*
// comment out
#include "G4VisManager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
*/

#include "G4VisExecutive.hh"
// add
#include "G4UIExecutive.hh"


/// Geometry
#include "G4XraySpecularReflectingSurface.hh"
#include "ExXRTG4DetectorConstruction.hh"


/// Physics lists
#include "ExXRTG4PhysicsList.hh"

//#include "G4PhysListFactory.hh"
//#include "G4VModularPhysicsList.hh"

//#include "MyPhysicsList.hh"
//#include "G4EmStandardPhysics.hh"


/// Primary generator and User action
/*
#include "ExXRTG4PrimaryGeneratorAction.hh"
#include "ExXRTG4SteppingAction.hh"
#include "ExXRTG4UserTrackingAction.hh"
*/


int main(int argc, char * argv[]) {

 
  G4RunManager* runManager = new G4RunManager;
  
  // micro-roughness reduces the reflectivity;


  //const G4double micro_roughness = 2.5*nm;
  const G4double micro_roughness = 0.55*nm;

  // true sets the reflectivity to be 1 always: useful to debug the
  // geometry, also inconjunction with a large angle.

//   const G4bool  force_unit_reflectivity = true;
//   const G4double angle = 30*deg;

  const G4bool  force_unit_reflectivity = false;

  /*
  const G4double angle = 0.5*deg;
  // photon energy
  const G4double energy = 1*keV;
  */
  
  // Setup mandatory user initialization: Geometry
  G4XraySpecularReflectingSurface * xray_surface_property =
    new G4XraySpecularReflectingSurface("TargetSiliconSurface", micro_roughness, force_unit_reflectivity);
  runManager->SetUserInitialization(new ExXRTG4DetectorConstruction(xray_surface_property));

  // Setup mandatory user initialization: Physics list
  //runManager->SetUserInitialization(new ExXRTG4PhysicsList);


  /// Test MyPhysicsList
  //runManager->SetUserInitialization(new MyPhysicsList);

  
  /// Use factgory
  //G4PhysListFactory factory;
  //G4VModularPhysicsList* physicsList = factory.GetReferencePhysList( "FTFP_BERT_LIV" );
  //runManager->SetUserInitialization( physicsList );
  
  G4VModularPhysicsList* xrtg4physicsList = new ExXRTG4PhysicsList;
  //xrtg4physicsList->RegisterPhysics(new G4EmStandardPhysics());
  runManager->SetUserInitialization( xrtg4physicsList );

  //physicsList->RegisterPhysics(xrtg4physicsList);


  //physicsList->RegisterPhysics( new ExXRTG4PhysicsList );
  //runManager->SetUserInitialization( physicsList );
  
  // Set up user initialization: User Actions
  runManager->SetUserInitialization( new UserActionInitialization );

  /*
  runManager->SetUserAction(new ExXRTG4PrimaryGeneratorAction(energy, angle));
  runManager->SetUserAction(new ExXRTG4SteppingAction);
  runManager->SetUserAction(new ExXRTG4UserTrackingAction);
  */

  /// Initialize G4 Kernel
  runManager->Initialize();



  /// Geant4 Tutorial code

  // Create visualization environment
  auto visManager = new G4VisExecutive;
  visManager->Initialize();
   
  // Get UI (User Interface) manager
  G4UImanager* uiManager = G4UImanager::GetUIpointer();

  // Switch batch or interactive mode
  if ( argc == 1 ) {  // Interactive mode - no command argument
    auto uiExec = new G4UIExecutive( argc, argv );
    uiManager->ApplyCommand( "/control/execute GlobalSetup.mac" ); // default macro
    uiExec->SessionStart();
    delete uiExec;
  } else {            // Batch mode - 1st command argument is a macro-filename
    G4String macroName = argv[1];
    uiManager->ApplyCommand( "/control/execute " + macroName );
  }
  
  // Job termination
  //delete uiExec;
  delete visManager;
  delete runManager;
  
  return 0;
}

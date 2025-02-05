#include "RunAction.hh"
//#include "Analysis.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

RunAction::RunAction()
  : G4UserRunAction()
{
    // Get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    //G4cout << "Using " << analysisManager->GetType() << G4endl;
    analysisManager->SetActivation(true);

    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFileName("mpo_test");
    
    // Creating EventHit Tree
    //G4int id;
    analysisManager->CreateNtuple("tEventHits","EventHits");
    analysisManager->CreateNtupleIColumn("EventID");     // column Id = 0
    analysisManager->CreateNtupleIColumn("NumHitData");  // column Id = 1
    analysisManager->CreateNtupleIColumn("HitID");       // column Id = 2
    analysisManager->CreateNtupleDColumn("posX");        // column Id = 3
    analysisManager->CreateNtupleDColumn("posY");        // column Id = 4
    analysisManager->CreateNtupleDColumn("posZ");        // column Id = 5
    analysisManager->CreateNtupleDColumn("dirX");        // column Id = 6
    analysisManager->CreateNtupleDColumn("dirY");        // column Id = 7
    analysisManager->CreateNtupleDColumn("dirZ");        // column Id = 8
    analysisManager->CreateNtupleDColumn("energy");      // column Id = 9
    
}

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  analysisManager->OpenFile();
  //G4cout << "File " << fileName << " Open" <<G4endl;
}

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // save the file and close
  analysisManager->Write();
  analysisManager->CloseFile();
}

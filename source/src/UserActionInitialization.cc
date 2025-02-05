//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// UserActionInitialization.cc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "UserActionInitialization.hh"
#include "PrimaryGenerator.hh"
#include "RunAction.hh"

/// add for ExXRTG4
#include "ExXRTG4SteppingAction.hh"
#include "ExXRTG4UserTrackingAction.hh"


//------------------------------------------------------------------------------
  UserActionInitialization::UserActionInitialization()
  : G4VUserActionInitialization()
//------------------------------------------------------------------------------
{}

//------------------------------------------------------------------------------
  UserActionInitialization::~UserActionInitialization()
//------------------------------------------------------------------------------
{}

//------------------------------------------------------------------------------
  void UserActionInitialization::Build() const
//------------------------------------------------------------------------------
{
  SetUserAction( new PrimaryGenerator() );

  // SetUserAction
  SetUserAction(new ExXRTG4SteppingAction );
  SetUserAction(new ExXRTG4UserTrackingAction );

  SetUserAction( new RunAction() );
}

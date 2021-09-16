/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"

//.........................................................................

ActionInitialization::ActionInitialization(DetectorConstruction* detector)
 : G4VUserActionInitialization(),
   fDetector(detector)
{}

ActionInitialization::~ActionInitialization()
{}

//.........................................................................

void ActionInitialization::BuildForMaster() const
{
  RunAction* runAction = new RunAction(fDetector, 0);
  SetUserAction(runAction);
}

void ActionInitialization::Build() const
{
  PrimaryGeneratorAction* primary = new PrimaryGeneratorAction();
  SetUserAction(primary);

  RunAction* runAction = new RunAction(fDetector, primary );
  SetUserAction(runAction);
  
  EventAction* event = new EventAction(primary);
  SetUserAction(event);  
  
  TrackingAction* trackingAction = new TrackingAction(fDetector,event);
  SetUserAction(trackingAction);
  
  SteppingAction* steppingAction = new SteppingAction(fDetector, event);
  SetUserAction(steppingAction);
}  

G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const
{
  return new SteppingVerbose();
}  

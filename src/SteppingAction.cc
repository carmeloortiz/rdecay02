//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4OpticalPhoton.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event)
: G4UserSteppingAction(), fDetector(det), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());    
  G4int  event = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //which volume ?
  //
  G4LogicalVolume* lVolume = aStep->GetPreStepPoint()->GetTouchableHandle()
                             ->GetVolume()->GetLogicalVolume();
  G4int iVol = 0;
  G4Track* track = aStep->GetTrack();
  G4double tID = track->GetTrackID();



  if (lVolume == fDetector->GetLogicTarget())   iVol = 1;
  if (lVolume == fDetector->GetLogicDetector()) iVol = 2;


  // count processes
  // 
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  const G4VProcess* tprocess   = endPoint->GetProcessDefinedStep();
  const G4VProcess* sprocess   = aStep->GetPreStepPoint()->GetProcessDefinedStep();
  run->CountProcesses(tprocess, iVol);


  G4TouchableHandle touch = endPoint->GetTouchableHandle();
  G4VPhysicalVolume* eVolume = touch->GetVolume();
  G4String eVname("null");
  if (eVolume)
    {
      eVname = eVolume->GetName();  
      //            std::cout << "SteppingAction eVname: " << eVname << std::endl;
      if (lVolume == fDetector->GetLogicSiPM() || eVolume->GetLogicalVolume() == fDetector->GetLogicSiPM() || 
	  eVname.find("SiPM")!=std::string::npos || (lVolume->GetName()).find("SiPM")!=std::string::npos)     iVol = 3;
    }

  // energy deposit
  //
  G4double edepStep = aStep->GetTotalEnergyDeposit();

  tprocess = aStep->GetPostStepPoint()->GetProcessDefinedStep();

  const G4ParticleDefinition* particle = track->GetParticleDefinition();  
  G4int pID       = particle->GetPDGEncoding();

  if (edepStep <= 0. && (pID!=0 && pID!=-22) ) return; // the deception version of G4 uses -22 for optical photons; my Mac's uses 0.
  G4double time   = aStep->GetPreStepPoint()->GetGlobalTime();
  G4double weight = aStep->GetPreStepPoint()->GetWeight();   

  if (iVol!=3)
    fEventAction->AddEdep(iVol, edepStep, time, weight);

  //fill ntuple id = 2
  G4int id = 4;   
  const G4double length = aStep->GetStepLength();
  const G4ThreeVector pos(aStep->GetPreStepPoint()->GetPosition());
  const G4ThreeVector tpos(aStep->GetPostStepPoint()->GetPosition());



  std::string startp("null");
  std::string endp("null");


  analysisManager->FillNtupleDColumn(id,0, edepStep);
  analysisManager->FillNtupleDColumn(id,1, time/s);
  analysisManager->FillNtupleDColumn(id,2, weight);
  analysisManager->FillNtupleDColumn(id,3, pos[0]/mm);
  analysisManager->FillNtupleDColumn(id,4, pos[1]/mm);
  analysisManager->FillNtupleDColumn(id,5, pos[2]/mm);
  analysisManager->FillNtupleDColumn(id,6, length/mm);
  analysisManager->FillNtupleDColumn(id,7, event);
  analysisManager->FillNtupleDColumn(id,8, pID);
  analysisManager->FillNtupleDColumn(id,9, tID);

  analysisManager->FillNtupleSColumn(id,10, track->GetVolume()->GetLogicalVolume()->GetName());
  if (eVolume)
    analysisManager->FillNtupleDColumn(id,11, eVolume->GetCopyNo());
  else
    analysisManager->FillNtupleDColumn(id,11, 0);
  analysisManager->FillNtupleSColumn(id,12, track->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName());
  analysisManager->FillNtupleDColumn(id,13, track->GetCurrentStepNumber());
  if (sprocess)
      startp = sprocess->GetProcessName();
  if (tprocess)
    endp = tprocess->GetProcessName();
  analysisManager->FillNtupleSColumn(id,14, startp);
  analysisManager->FillNtupleSColumn(id,15, endp);
  analysisManager->FillNtupleDColumn(id,16, tpos[0]/mm);
  analysisManager->FillNtupleDColumn(id,17, tpos[1]/mm);
  analysisManager->FillNtupleDColumn(id,18, tpos[2]/mm);
  analysisManager->FillNtupleSColumn(id,19,eVname);

  if (eVname=="SiPM") {
	  fEventAction->AddEdep(3, 1.0, time, weight);	  
  }
  analysisManager->AddNtupleRow(id);      



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

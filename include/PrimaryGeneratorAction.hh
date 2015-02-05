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
//      ----------------------------------------------------------------
//                          AFRODITE (iThemba Labs)
//      ----------------------------------------------------------------
//
//      Github repository: https://www.github.com/KevinCWLi/AFRODITE
//
//      Main Author:    K.C.W. Li
//      Edits: Steve Peterson & Vijitha Ramanathan
//
//      email: likevincw@gmail.com
//

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

/// The primary generator action class with particle gum.
///
/// It defines a single particle hitting the afrodite target
/// perpendicular to the its surface. The particle characteristics
/// can be changed via the macros defined in PrimaryGeneratorAction
/// or using the G4 build-in commands of G4ParticleGun class.

class PrimaryGeneratorMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    // Sets the default event characteristics
    void          SetDefaultKinematics();

    // Generates the primary event via the ParticleGun method.
    void          GeneratePrimaries(G4Event* event);
    
    // Set methods
    void           SetRandomFlag(G4bool value);

    // Set the verbosity level.
    inline void    SetVerbosity(G4int a)      {verbosityLevel = a;};

    // Methods to change the parameters of primary particle generation interactively
    void           SetStartingParticleType(G4String);
    

private:
    // G4 particle gun
    G4ParticleGun*               fParticleGun; 
    // Initialize messenger to accept macros commands    
    PrimaryGeneratorMessenger*   gunMessenger;
    
    // Verbosity
    G4int                        verbosityLevel;
    
    // parameters of primary particle generation
    G4String                     particle_type;
    G4ParticleDefinition*        particle_definition;
    G4double                     mean_particle_energy;
    G4double                     particle_energy;
    G4ThreeVector                starting_particle_position;
    G4ThreeVector                particle_position;
    G4ThreeVector                starting_particle_momentum_direction;
    G4ThreeVector                particle_momentum_direction;


    G4double    mx;
    G4double    my;
    G4double    mz;
    
    G4double    theta;
    G4double    phi;
    
    G4double    beamEnergy;
    G4double    energy;
    G4double    recoilExcitationEnergy;
    G4double    alphaSeperationEnergy;
    G4double    protonSeperationEnergy;
    G4double    daughterExcitationEnergy;
    G4int       alphaDecayMode;
    
    G4ThreeVector ejectileDirection;
    G4ThreeVector recoilDirection;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif



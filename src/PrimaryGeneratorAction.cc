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

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"

#include "G4IonTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
{
  // Set default particle attributes
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
  SetDefaultKinematics();

  // Create a messenger for this class
  //gunMessenger = new PrimaryGeneratorMessenger(this);

  // Set verbosityLevel
  verbosityLevel = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  //delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetDefaultKinematics()
{

  // Define primary particles: protons, He3, neutron, mu-, gamma, e-, alpha, geantino
  particle_type = "gamma";

  // Define the parameters for the energy of primary particles
  G4double default_mean_particle_energy = 6.0*MeV;
  mean_particle_energy = default_mean_particle_energy;

  // Define the parameters of the initial position (using previous default, for now)
  G4ThreeVector default_particle_position = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
  starting_particle_position = default_particle_position;

  // Define the parameters of the momentum of primary particles
  G4ThreeVector default_particle_momentum_direction = G4ThreeVector(1.0, 0.0, 0.0);
  starting_particle_momentum_direction = default_particle_momentum_direction;

  // Set basic default values
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particle_type);
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(mean_particle_energy);
  fParticleGun->SetParticlePosition(starting_particle_position);
  fParticleGun->SetParticleMomentumDirection(starting_particle_momentum_direction);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    
    // Set the starting primary particles
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(particle_type);
    particle_definition = particle;
    fParticleGun->SetParticleDefinition(particle);

    ///////////////////////////////////////////////////////////////
    // Particle Energy Characteristics
    ///////////////////////////////////////////////////////////////

    particle_energy = mean_particle_energy;
    
    /*
     G4double InitialEnergy = G4RandFlat::shoot( 0., 20.);
     //G4double InitialEnergy = G4RandGamma::shoot( 0.5, 1.);
    */


    ///////////////////////////////////////////////////////////////
    // Particle Position Characteristics
    ///////////////////////////////////////////////////////////////

    particle_position = starting_particle_position;


    ///////////////////////////////////////////////////////////////
    // Particle Momentum/Direction Characteristics
    ///////////////////////////////////////////////////////////////

    // ISOTROPIC
    //   Exploiting the spherical symmetry of normal distributions

    /*
     mx = G4RandGauss::shoot( 0, 1.);
     my = G4RandGauss::shoot( 0, 1.);
     mz = G4RandGauss::shoot( 0, 1.);
    */

    // ISOTROPIC
    //   Alternate method
    G4double cosTheta = 2 * G4UniformRand() - 1.0;
    G4double phi      = CLHEP::twopi * G4UniformRand();
    G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    mx = sinTheta * std::cos(phi);
    my = sinTheta * std::sin(phi);
    mz = cosTheta;
    // Set particle momentum direction
    particle_momentum_direction = G4ThreeVector(mx, my, mz).unit();


    ///////////////////////////////////////////////////////////////
    // Particle Time Distribution
    ///////////////////////////////////////////////////////////////
    
    /*
     initialclocktime = G4RandGauss::shoot( 0, 1.5);
     fParticleGun->SetParticleTime(initialclocktime*ns);
    */
    

    ///////////////////////////////////////////////////////////////
    // Radioactive Decay-enabled Particles
    ///////////////////////////////////////////////////////////////
    
    /*
     //G4int Z = 3, A = 11;
     //G4int Z = 6, A = 21;
     //G4int Z = 10, A = 18;
     //G4int Z = 11, A = 22;
     G4int Z = 27, A = 60;
     //G4int Z = 63, A = 152;
     
     G4double ionCharge   = 0.*eplus;
     G4double excitEnergy = 0.*MeV;
     
     //G4ParticleDefinition* ion = G4ParticleTable::GetParticleTable()->GetIon(Z,A,excitEnergy);
     G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
     //ion->SetPDGLifeTime(1*ns);
     fParticleGun->SetParticleDefinition(ion);
     fParticleGun->SetParticleCharge(ionCharge);
    */
    
    /*
     ////////    4He, +1 charge
     G4int Z = 2, A = 4;
     
     G4double ionCharge   = 1.*eplus;
     G4double excitEnergy = 0.*MeV;
     
     //G4ParticleDefinition* ion = G4ParticleTable::GetParticleTable()->GetIon(Z,A,excitEnergy);
     G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
     //ion->SetPDGLifeTime(1*ns);
     fParticleGun->SetParticleDefinition(ion);
     fParticleGun->SetParticleCharge(ionCharge);
    */


    ///////////////////////////////////////////////////
    // 16O(a,a')
    ///////////////////////////////////////////////////
    /*
     G4double test = G4UniformRand();
     
     //if(test<=0.5) beamEnergy = G4RandGauss::shoot( 200, 0.5); // MeV
     //if(test>0.5) beamEnergy = G4RandGauss::shoot( 180, 0.5); // MeV
     
     if(test>0.0 && test<=0.1) beamEnergy = G4RandGauss::shoot( 215, 0.5); // MeV
     if(test>0.1 && test<=0.2) beamEnergy = G4RandGauss::shoot( 210, 0.5); // MeV
     if(test>0.2 && test<=0.3) beamEnergy = G4RandGauss::shoot( 200, 0.5); // MeV
     if(test>0.3 && test<=0.4) beamEnergy = G4RandGauss::shoot( 190, 0.5); // MeV
     if(test>0.4 && test<=0.5) beamEnergy = G4RandGauss::shoot( 180, 0.5); // MeV
     if(test>0.5 && test<=0.6) beamEnergy = G4RandGauss::shoot( 185, 0.5); // MeV
     if(test>0.6 && test<=0.7) beamEnergy = G4RandGauss::shoot( 170, 0.5); // MeV
     if(test>0.7 && test<=0.8) beamEnergy = G4RandGauss::shoot( 160, 0.5); // MeV
     if(test>0.8 && test<=0.9) beamEnergy = G4RandGauss::shoot( 150, 0.5); // MeV
     if(test>0.9 && test<=1.0) beamEnergy = G4RandGauss::shoot( 140, 0.5); // MeV
     
     //beamEnergy = G4RandGauss::shoot( 200, 0.5); // MeV
     //beamEnergy = 200; // MeV
     
     recoilExcitationEnergy = 15.097;  // MeV
     alphaSeperationEnergy = 7.16192;    // MeV
     protonSeperationEnergy = 12.1274;  // MeV
     //  Elastically scattered outgoing alpha
     mx = G4RandGauss::shoot( 0, 0.00001);
     my = G4RandGauss::shoot( 0, .05);
     //mx = G4RandGauss::shoot( 0, 0.00001);
     //my = G4RandGauss::shoot( 0, 0.00001);
     //mx = 0;
     //my = 0;
     
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, 10.));
     energy = beamEnergy - recoilExcitationEnergy; // MeV
     fParticleGun->SetParticleEnergy(energy*MeV);
     fParticleGun->GeneratePrimaryVertex(anEvent);
     
     
     //////////////////////////////////////////////////////////////
     //                      16O -> 4He + 12C
     //////////////////////////////////////////////////////////////
     //  Alpha 0 isotropic decay:    15.097 MeV state of 16O -> to Ground state of 12C
     //  Alpha 1 anisotropic decay:  15.097 MeV state of 16O -> to 4.43891 MeV state of 12C
     //  Assuming recoil nucleus takes no kinetic energy from inelastic scatter
     
     G4double pAlphaDecayMode = G4UniformRand();
     
     if(pAlphaDecayMode<=0.5) alphaDecayMode = 0;
     if(pAlphaDecayMode>0.5) alphaDecayMode = 1;
     
     //alphaDecayMode = 1;
     
     ////    ALPHA 1
     if(alphaDecayMode == 0) daughterExcitationEnergy = 0; //    MeV
     
     ////    Alpha 0
     if(alphaDecayMode == 1) daughterExcitationEnergy = 4.43891; //    MeV
     
     mx = G4RandFlat::shoot( -1., 1.);
     my = G4RandFlat::shoot( -1., 1.);
     mz = G4RandFlat::shoot( -1., 1.);
     
     fParticleGun->SetParticleEnergy( (recoilExcitationEnergy - alphaSeperationEnergy - daughterExcitationEnergy)*(3/4)*MeV);
     fParticleGun->SetParticleMomentumDirection(ejectileDirection);
     fParticleGun->GeneratePrimaryVertex(anEvent);
    */
    
    
    ////    Recoil - 12C
    /*
     fParticleGun->SetParticleEnergy( (recoilExcitationEnergy - alphaSeperationEnergy - daughterExcitationEnergy)*(1/4)*MeV);
     
     recoilDirection = -ejectileDirection;
     fParticleGun->SetParticleMomentumDirection(recoilDirection);
    */
    
    /*
     Z = 6, A = 12;
     
     G4ParticleDefinition* recoil = G4IonTable::GetIonTable()->GetIon(Z, A, 4.43891*MeV);
     //G4ParticleDefinition* recoil = G4ParticleTable::GetParticleTable()->GetIon(Z, A, 0.*keV);
     fParticleGun->SetParticleCharge(0.*eplus);
     //recoil->SetPDGLifeTime(1*ns);
     //recoil->SetPDGStable(true);
     fParticleGun->SetParticleDefinition(recoil);
     
     fParticleGun->SetParticleEnergy( (recoilExcitationEnergy - alphaSeperationEnergy - daughterExcitationEnergy)*(1/4)*MeV);
     
     recoilDirection = -ejectileDirection;
     fParticleGun->SetParticleMomentumDirection(recoilDirection);
    */
    
    //fParticleGun->GeneratePrimaryVertex(anEvent);
    //G4double Lifetime = fParticleGun->GetParticleDefinition()->GetPDGLifeTime();
    //G4double Lifetime = recoil->GetPDGLifeTime();
    //G4cout << "Here is the alphaSeperationEnergy    "<< alphaSeperationEnergy << G4endl;
    
    
    
    ////     Gamma Decay - from Recoil Nucleus
    /*
     if(alphaDecayMode == 1)
     {
     G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
     fParticleGun->SetParticleDefinition(particleDefinition);
     
     fParticleGun->SetParticleEnergy(daughterExcitationEnergy*MeV);
     
     mx = G4RandFlat::shoot( -1., 1.);
     my = G4RandFlat::shoot( -1., 1.);
     mz = G4RandFlat::shoot( -1., 1.);
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
     
     fParticleGun->GeneratePrimaryVertex(anEvent);
     }
    */

    ////////////////////////////////////////////////////////
    //          FINALIZING PARTICLE SETTINGS 
    ////////////////////////////////////////////////////////

    // Set final beam energy
    fParticleGun->SetParticleEnergy(particle_energy);

    // Set final particle position
    fParticleGun->SetParticlePosition(particle_position);

    // Set final particle momentum direction
    fParticleGun->SetParticleMomentumDirection(particle_momentum_direction);

    // Generate a primary particle
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    ////////////////////////////////////////////////////////
    //          PARTICLE INFORMATION OUTPUT
    ////////////////////////////////////////////////////////

    // Print out starting particle details
    if (verbosityLevel > 0) {
      if (anEvent->GetEventID() == 0 || verbosityLevel > 1) {
	G4cout << G4endl;
        G4cout << "----- Primary Generator Source Details for Event: " << anEvent->GetEventID()
	       << " -----" << G4endl;
        G4cout << "  Particle name: " << particle_definition->GetParticleName() << G4endl;
        G4cout << "         Energy: " << particle_energy/MeV << " MeV" << G4endl;
        G4cout << "       Position: " << particle_position/mm << " mm" << G4endl;
        G4cout << "      Direction: " << particle_momentum_direction/mrad << " mrad" << G4endl;
	G4cout << G4endl;
      }
    }
    
    /*
     G4double Lifetime = fParticleGun->GetParticleDefinition()->GetPDGLifeTime();
     G4double Spin = fParticleGun->GetParticleDefinition()->GetPDGSpin();
     G4double Isospin = fParticleGun->GetParticleDefinition()->GetPDGIsospin();
     G4double Parity = fParticleGun->GetParticleDefinition()->GetPDGiGParity();
     G4String ParticleName = fParticleGun->GetParticleDefinition()->GetParticleName();

     G4cout << "Here is the Lifetime    "<< Lifetime << G4endl;
     G4cout << "Here is the Spin    "<< Spin << G4endl;
     G4cout << "Here is the Isospin    "<< Isospin << G4endl;
     G4cout << "Here is the Parity    "<< Parity << G4endl;
     G4cout << "Here is the ParticleName    "<< ParticleName << G4endl;
    */

}

// ********************************************************************
void PrimaryGeneratorAction::SetStartingParticleType (G4String str)
{ particle_type = str; }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


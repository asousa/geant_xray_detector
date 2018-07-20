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
// $Id: DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  // G4Material* env_mat = nist->FindOrBuildMaterial("G4_");

    // Material: Vacuum
  G4Material* env_mat = new G4Material("Vacuum",
              1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
              kStateGas, 2.73*kelvin, 3.0E-18*pascal );
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  // G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* world_mat = new G4Material("Vacuum",
              1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
              kStateGas, 2.73*kelvin, 3.0E-18*pascal );
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
      

  G4double detector_dimX = 10/2.*mm;
  G4double detector_dimY = 2.2/2.*mm;
  G4double detector_dimZ = 5./2.*mm;
  G4double baffle_thickness = 0.5/2.*mm;
  G4double baffle_height    = 20/2.*mm;

  G4int n_detectors_right = 2;
  G4int n_detectors_left  = 1;

  

  // // Shape 1
  // G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_W");
  // G4ThreeVector pos1 = G4ThreeVector(0, 0, -3*cm);
        
  // // Conical section shape       
  // G4double shape1_rmina =  0.*cm, shape1_rmaxa = 4.*cm;
  // G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  // G4double shape1_hz = 0.5*mm;
  // G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  // G4Cons* solidShape1 =    
  //   new G4Cons("window", 
  //   shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
  //   shape1_phimin, shape1_phimax);
                      
  // G4LogicalVolume* windowShape =                         
  //   new G4LogicalVolume(solidShape1,         //its solid
  //                       shape1_mat,          //its material
  //                       "window");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos1,                    //at position
  //                   windowShape,             //its logical volume
  //                   "window",                //its name
  //                   logicEnv,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  
  // ----------------------------------------------------------------
  // Detector elements:
  // ----------------------------------------------------------------
  // G4Material* detector_mat = nist->FindOrBuildMaterial("G4_Si");

  G4Element* Cd = new G4Element("Cadmium","Cd",48., 112.41*g/mole);
  G4Element* Zn = new G4Element("Zinc","Zn", 30., 65.38*g/mole);
  G4Element* Te = new G4Element("Tellurium","Te", 52., 127.60*g/mole);
  G4Material* CZT = new G4Material("CZT", 5.8*g/cm3, 3);
  CZT->AddElement(Cd, 48*perCent);
  CZT->AddElement(Zn, 02*perCent);
  CZT->AddElement(Te, 50*perCent);



  G4ThreeVector detector_pos  = G4ThreeVector(0, 0, 0);
  // // G4ThreeVector detector_pos2 = G4ThreeVector(0,  2*detector_dimY + 2*baffle_thickness, 0);
  // G4ThreeVector detector_pos3 = G4ThreeVector(0, -2*detector_dimY - 2*baffle_thickness, 0);

  std::ostringstream detname;

  G4VSolid* detector_solid = new G4Box("detector", 
                   detector_dimX, detector_dimY, detector_dimZ);

  for (int i=-n_detectors_left; i<= n_detectors_right; i++){
    detname.str("");
    detname << "detector" << i;
    
    detector_pos  = G4ThreeVector(0, 2*i*(detector_dimY + baffle_thickness), 0);

    G4LogicalVolume* detector =                         
    new G4LogicalVolume(detector_solid,      //its solid
                        CZT,        //its material
                        detname.str());         //its name
               
    new G4PVPlacement(0,                       //no rotation
                    detector_pos,                    //at position
                    detector,                //its logical volume
                    detname.str(),                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  }







  // G4VSolid* detector_solid = new G4Box("detector", 
  //                  detector_dimX, detector_dimY, detector_dimZ);

  // G4LogicalVolume* detector =                         
  //   new G4LogicalVolume(detector_solid,      //its solid
  //                       detector_mat,        //its material
  //                       "detector");         //its name
               
  //   new G4PVPlacement(0,                       //no rotation
  //                   detector_pos,                    //at position
  //                   detector,                //its logical volume
  //                   "detector",                //its name
  //                   logicEnv,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking
           

  // // side detectors:

  //   G4ThreeVector detector_sidepos = G4ThreeVector(0, 2*(detector_dimY + baffle_thickness), 0);
  //   G4LogicalVolume* detector_side =                         
  //     new G4LogicalVolume(detector_solid,      //its solid
  //                           detector_mat,        //its material
  //                           "detector_right_1");         //its name
      
  //   new G4PVPlacement(0,                       //no rotation
  //                       detector_sidepos,                    //at position
  //                       detector_side,                //its logical volume
  //                       "detector_right_1",                //its name
  //                       logicEnv,                //its mother  volume
  //                       false,                   //no boolean operation
  //                       0,                       //copy number
  //                       checkOverlaps);          //overlaps checking    }              

  //   detector_sidepos = G4ThreeVector(0,  -2*(detector_dimY + baffle_thickness), 0);

  //   detector_side =                         
  //       new G4LogicalVolume(detector_solid,      //its solid
  //                           detector_mat,        //its material
  //                           "detector_left_1");         //its name
      
  //       new G4PVPlacement(0,                       //no rotation
  //                       detector_sidepos,                    //at position
  //                       detector_side,                //its logical volume
  //                       "detector_left_1",                //its name
  //                       logicEnv,                //its mother  volume
  //                       false,                   //no boolean operation
  //                       0,                       //copy number
  //                       checkOverlaps);          //overlaps checking    }


  //   detector_sidepos = G4ThreeVector(0,  -4*(detector_dimY + baffle_thickness), 0);

  //   detector_side =                         
  //       new G4LogicalVolume(detector_solid,      //its solid
  //                           detector_mat,        //its material
  //                           "detector_left_2");         //its name
      
  //       new G4PVPlacement(0,                       //no rotation
  //                       detector_sidepos,                    //at position
  //                       detector_side,                //its logical volume
  //                       "detector_left_2",                //its name
  //                       logicEnv,                //its mother  volume
  //                       false,                   //no boolean operation
  //                       0,                       //copy number
  //                       checkOverlaps);          //overlaps checking    }




  // ----------------------------------------------------------------
  // Baffles:
  // ----------------------------------------------------------------

  G4Material* baffle_material = nist->FindOrBuildMaterial("G4_W");
  G4VSolid*   baffle_solid = new G4Box("baffle", detector_dimX, baffle_thickness,  baffle_height + detector_dimZ);
  // G4ThreeVector baffle1_pos = G4ThreeVector(0,  detector_dimY + baffle_thickness,  baffle_height);
  // G4ThreeVector baffle2_pos = G4ThreeVector(0, -detector_dimY - baffle_thickness,  baffle_height);

  // G4LogicalVolume* baffle1 =                         
  // new G4LogicalVolume(baffle_solid,      //its solid
  //                       baffle_material,        //its material
  //                       "baffle1");         //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   baffle1_pos,             //at position
  //                   baffle1,                 //its logical volume
  //                   "baffle1",               //its name
  //                   logicEnv,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking



  // G4LogicalVolume* baffle2 =                         
  // new G4LogicalVolume(baffle_solid,      //its solid
  //                       baffle_material,        //its material
  //                       "baffle2");         //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   baffle2_pos,             //at position
  //                   baffle2,                 //its logical volume
  //                   "baffle2",               //its name
  //                   logicEnv,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking
  

  G4ThreeVector baffle_pos;
  std::ostringstream baffle_name;
  for (int i=-1*(n_detectors_left + 1); i <= n_detectors_right; i++) {
    baffle_name.str("");
    baffle_name << "baffle" << i;
    baffle_pos = G4ThreeVector(0, (2*i + 1)*(detector_dimY + baffle_thickness),  baffle_height);
    G4LogicalVolume* baffle =                         
    new G4LogicalVolume(baffle_solid,      //its solid
                          baffle_material,        //its material
                          baffle_name.str());         //its name
                 
    new G4PVPlacement(0,                       //no rotation
                      baffle_pos,             //at position
                      baffle,                 //its logical volume
                      baffle_name.str(),               //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

  }
  // ----------------------------------------------------------------

  // Set scoring volume
  // fScoringVolume = detector;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//
// Modidfed by Mutsumi S. 
//

// ExXRTG4DetectorConstruction.hh
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
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include <random>
#include "G4SystemOfUnits.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4AssemblyVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "ExXRTG4DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4XraySpecularReflectingSurface.hh"



ExXRTG4DetectorConstruction::ExXRTG4DetectorConstruction(G4XraySpecularReflectingSurface * property)
  : G4VUserDetectorConstruction(), xray_surface_property(property) {}

ExXRTG4DetectorConstruction::~ExXRTG4DetectorConstruction() {}

G4VPhysicalVolume* ExXRTG4DetectorConstruction::Construct() {

//G4NistManager* nist = G4NistManager::Instance();


/*
  Reflaction index data 
  text data format
  Energy[eV]  beta  delta
  
  get data from https://henke.lbl.gov/optical_constants/getdb2.html
  and swap 2nd and 3rd rows.

  sample files
  si.dat 
  ir.dat
  al.dat
*/

  G4Material *vacuum = new G4Material("Vacuum", 1.0 , 1.01*g/mole, 1.0E-25*g/cm3, kStateGas, 2.73*kelvin, 3.0E-18*pascal );
  //G4Material* silicon = nist->FindOrBuildMaterial("G4_AIR");
  G4Material *silicon = new G4Material("Si", 14, 28.0855*g/mole, 2.329*g/cm3);
  G4Material *pt = new G4Material("pt", 78, 195.08*g/mole, 21.45*g/cm3);
  G4Material *al = new G4Material("al", 13., 26.98*g/mole,2.70*g/cm3);
  G4Material *ir = new G4Material("ir", 77., 192.22*g/mole,22.42*g/cm3);
  

  
  //---------world volume-------------
  G4Box* solidWorld = new G4Box("World", 1.5*m, 1.5*m, 1.5*m);  // full size  3 x 3 x 3 m
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, vacuum, "logicWorld"); 
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "physWorld", 0, false, 0, 0);
  


  //---------------- MPO mirror geometry ------------
  /// EP
  /// curvature radius 750 mm
  /// focal length = 375 mm
  /// HiZ
  /// curvature radius 600 mm
  /// focal length = 300 mm
  const G4double frad = 600.*mm; 
  const G4double flen = frad/2.0;
  const G4ThreeVector vfrz = G4ThreeVector(0,0,frad);


  /// pore outer size = 26.0 um (half size 13.0 um, margin 13.0-25.0/2 = 0.5 um)
  /// pore inner size = 20.0 um (half size 10.0 um)
  /// pore height = 1.2 mm  = 1200 um
  /// pore XY pitch = 20.0 + 5.0 = 25.0 um 

  const G4double poreHoleWid   = 20.0*um;
  const G4double poreWallWid   = 5.0*um;
  const G4double poreDepth     = 1.2*mm;
  const G4double poreStepAngle = atan((poreHoleWid+poreWallWid)/2./frad)*2.;

  //const G4double dxout = (poreHoleWid+poreWallWid)*0.5;
  const G4double dxout = 26.0*um *0.5;
  const G4double dxin  = poreHoleWid*0.5;
  const G4double dzth  = poreDepth*0.5; 
  
  G4double dx1, dx2, dy1, dy2, dz;

  dz = dzth;
  dx2 = dxout*(frad+dz)/frad;
  dx1 = dx2*(frad-2*dz)/frad;
  dy2 = dx2;
  dy1 = dx1;
  G4VSolid* poreOuterTrd = new G4Trd("poreOuterTrd", dx1, dx2, dy1, dy2, dz);
  
  dz  = dzth*1.1;  /// pore hole height margin 10 %
  dx2 = dxin*(frad+dz)/frad;
  dx1 = dx2*(frad-2*dz)/frad;
  dy2 = dx2;
  dy1 = dx1;
  G4VSolid* poreInnerTrd = new G4Trd("poreInnderTrd", dx1, dx2, dy1, dy2, dz);

  G4VSolid* poreSolid = new G4SubtractionSolid("poreSolid", poreOuterTrd, poreInnerTrd, 0, G4ThreeVector(0,0,0));

  G4LogicalVolume* poreLogic= new G4LogicalVolume(poreSolid, ir, "poreLogic");
  new G4LogicalSkinSurface("TargetSurface", poreLogic, xray_surface_property);


  G4Colour colorBrown(0.7, 0.4, 0.1);
  G4VisAttributes* copperVisAttributes = new G4VisAttributes(colorBrown);
  poreLogic ->SetVisAttributes(copperVisAttributes);
  
  /// MPO size 25.um * 1600 = 40. mm (40. mm x 40. mm) 
  const G4int nxpole = 1600;
  const G4int nypole = 1600;
  
  G4AssemblyVolume* assemblyMPO = new G4AssemblyVolume();

  G4RotationMatrix Rx, Ry, Rm;
  G4ThreeVector Ta;
  G4Transform3D transform3d;
  G4double xangle, yangle;
  
  for (G4int iy=0; iy<nypole; iy++){
    yangle = poreStepAngle*(iy-nypole/2+0.5);
    Ry = G4RotationMatrix(); /// initialize
    Ry.rotateY(yangle);

    for (G4int ix=0; ix<nxpole; ix++){
      xangle = poreStepAngle*(ix-nxpole/2+0.5);
      Rx = G4RotationMatrix(); /// initialize
      Rx.rotateX(xangle);
      
      Rm = Ry*Rx;
      Ta = Rm*vfrz;
      transform3d = G4Transform3D(Rm, Ta);

      assemblyMPO->AddPlacedVolume(poreLogic, transform3d);
    }
  }
  

  /// place frame
  const G4double frameDep = 10*mm;
  const G4double frameWid = 3*mm;
  const G4double frameWidAngle = atan(frameWid/frad);
  const G4double mpoWidAngle  = poreStepAngle*nxpole;
  const G4double mpoStepAngle = mpoWidAngle+frameWidAngle;

  
  const G4double frmin = frad-frameDep/2;
  const G4double frmax = frad+frameDep/2;

  G4double fsphi, fdphi, fstheta, fdtheta;
  const G4int nxmpo = 3;
  const G4int nympo = 3;
  //const G4int nympo = 1;

  G4VSolid *frameXbar[nympo+1], *frameYbar[nxmpo+1];
  G4LogicalVolume *frameXbarLogic[nympo+1], *frameYbarLogic[nxmpo+1];


  Rm = G4RotationMatrix();  /// initialize
  Rm.rotateY(-90*deg);
  transform3d = G4Transform3D(Rm, G4ThreeVector(0,0,0));

  fsphi   = -0.5*nxmpo*mpoWidAngle - 0.5*(nxmpo+1)*frameWidAngle;
  fdphi   =    nxmpo*mpoWidAngle + (nxmpo+1)*frameWidAngle;
  fdtheta = frameWidAngle;
  for(G4int iy=0; iy<nympo+1; iy++){
    fstheta = 90*deg + (-0.5*nympo+iy)*mpoWidAngle+(-0.5*(nympo+1)+iy)*frameWidAngle;
    frameXbar[iy] = new G4Sphere("frameSphere", frmin, frmax, fsphi, fdphi, fstheta, fdtheta);
    frameXbarLogic[iy] = new G4LogicalVolume(frameXbar[iy], al, "frameLogic", 0,0,0);
    new G4LogicalSkinSurface("frameSurface", frameXbarLogic[iy], xray_surface_property);
    new G4PVPlacement(transform3d, frameXbarLogic[iy], "framePhysi", logicWorld, false, 0,0);
  }

  fstheta = 90*deg -0.5*nympo*mpoWidAngle - 0.5*(nympo+1)*frameWidAngle;
  fdtheta =    nympo*mpoWidAngle + (nympo+1)*frameWidAngle;
  fdphi   = frameWidAngle;
  for(G4int ix=0; ix<nxmpo+1; ix++){
    fsphi   = (-0.5*nxmpo+ix)*mpoWidAngle+(-0.5*(nxmpo+1)+ix)*frameWidAngle;
    frameYbar[ix] = new G4Sphere("frameSphere", frmin, frmax, fsphi, fdphi, fstheta, fdtheta);
    frameYbarLogic[ix] = new G4LogicalVolume(frameYbar[ix], al, "frameLogic", 0,0,0);
    new G4LogicalSkinSurface("frameSurface", frameYbarLogic[ix], xray_surface_property);
    new G4PVPlacement(transform3d, frameYbarLogic[ix], "framePhysi", logicWorld, false, 0,0);
  }


  /// Place Wall
  const G4double halfphi   = (nxmpo*mpoWidAngle + (nxmpo+1)*frameWidAngle)/2.;
  const G4double halftheta = (nympo*mpoWidAngle + (nympo+1)*frameWidAngle)/2.;
  const G4double walldepth = 3*mm;
  const G4double wallTopZ = frad * cos(halftheta);
  const G4double wallBotZ = frad/2.0;
  G4double trnX, trnY, trnZ;


  /// YZ wall

  /// Trapezoid
  /*
  dx1 = walldepth/2.0;
  dx2 = walldepth/2.0;
  dy1 = wallBotZ*tan(halfphi);
  dy2 = wallTopZ*tan(halfphi);
  dz  = (wallTopZ - wallBotZ)/2.;
  G4VSolid *wallYZsolid = new G4Trd("WallYZsolid", dx1, dx2, dy1, dy2, dz);
  G4LogicalVolume *wallYZlogic = new G4LogicalVolume(wallYZsolid, al, "WallYZlogic");
  new G4LogicalSkinSurface("WallSurface", wallYZlogic, xray_surface_property);

  trnX = frad*tan(halftheta) + walldepth/2.0;
  trnY = 0.0;
  trnZ = wallBotZ +  (wallTopZ - wallBotZ)/2.;
  transform3d = G4Transform3D(G4RotationMatrix(), G4ThreeVector(trnX, trnY, trnZ));
  new G4PVPlacement(transform3d, wallYZlogic, "WallYZphysi+", logicWorld, false, 0,0);
  transform3d = G4Transform3D(G4RotationMatrix(), G4ThreeVector(-trnX, trnY, trnZ));
  new G4PVPlacement(transform3d, wallYZlogic, "WallYZphysi-", logicWorld, false, 0,0);
  */
  
  dx1 = walldepth/2.0;
  //dy1 = wallBotZ*tan(halfphi);
  dy1 = wallTopZ*tan(halfphi);
  dz  = (wallTopZ - wallBotZ)/2.;
  G4VSolid *wallYZsolid = new G4Box("WallYZsolid", dx1, dy1, dz);
  G4LogicalVolume *wallYZlogic = new G4LogicalVolume(wallYZsolid, al, "WallYZlogic");
  new G4LogicalSkinSurface("WallSurface", wallYZlogic, xray_surface_property);

  trnX = frad*tan(halftheta) + walldepth/2.0 - 1.0*mm; // Margin 1.0mm
  trnY = 0.0;
  trnZ = wallBotZ +  (wallTopZ - wallBotZ)/2.;
  transform3d = G4Transform3D(G4RotationMatrix(), G4ThreeVector(trnX, trnY, trnZ));
  new G4PVPlacement(transform3d, wallYZlogic, "WallYZphysi+", logicWorld, false, 0,0);
  transform3d = G4Transform3D(G4RotationMatrix(), G4ThreeVector(-trnX, trnY, trnZ));
  new G4PVPlacement(transform3d, wallYZlogic, "WallYZphysi-", logicWorld, false, 0,0);

  
  /// XZ wall
  dx1 = wallTopZ*tan(halftheta);
  dy1 = walldepth/2.0;
  dz  = (wallTopZ - wallBotZ)/2.;
  G4VSolid *wallXZsolid = new G4Box("WallXZsolid", dx1, dy1, dz);
  G4LogicalVolume *wallXZlogic = new G4LogicalVolume(wallXZsolid, al, "WallXZlogic");
  new G4LogicalSkinSurface("WallSurface", wallXZlogic, xray_surface_property);

  trnX = 0.0;
  trnY = frad*tan(halfphi) + walldepth/2.0 - 1.0*mm; // Margin 1.0mm
  trnZ = wallBotZ +  (wallTopZ - wallBotZ)/2.;
  transform3d = G4Transform3D(G4RotationMatrix(), G4ThreeVector(trnX, trnY, trnZ));
  new G4PVPlacement(transform3d, wallXZlogic, "WallYZphysi+", logicWorld, false, 0,0);
  transform3d = G4Transform3D(G4RotationMatrix(), G4ThreeVector(trnX, -trnY, trnZ));
  new G4PVPlacement(transform3d, wallXZlogic, "WallYZphysi-", logicWorld, false, 0,0);


  
  /// Place MPO
  for(G4int iy=0; iy<nympo; iy++){
    Ry = G4RotationMatrix();  /// initialize
    Ry.rotateY( (iy-(nympo-1)*0.5)*mpoStepAngle );

    for(G4int ix=0; ix<nxmpo; ix++){
      Rx = G4RotationMatrix(); /// initialize
      Rx.rotateX( (ix-(nxmpo-1)*0.5)*mpoStepAngle );
      //Rx.rotateX(ix*mpoStepAngle);

      Rm = Ry*Rx;
      transform3d = G4Transform3D(Rm, G4ThreeVector(0,0,0));
      
      /// Place MPO
      assemblyMPO->MakeImprint(logicWorld, transform3d);
    }
  }


  //Imager box
  G4double imagerthink = 2.0*mm;
  G4Box* imagerBox = new G4Box("imagerBox", 20.*cm, 20.*cm, imagerthink/2.0);   // full size 40 cm x 40 cm x 2 mm 
  G4LogicalVolume* imagerLogicVolume = new G4LogicalVolume(imagerBox, silicon, "imagerLogicVolume", 0,0,0);
  //new G4LogicalSkinSurface("ImagerSurface", imagerLogicVolume, xray_surface_property);

  new G4PVPlacement(0, G4ThreeVector(0, 0, flen-imagerthink/2.0), imagerLogicVolume, "imager", logicWorld, false, 0,0);
  
  return physWorld;

}

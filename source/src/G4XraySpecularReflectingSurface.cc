//  Giuseppe Vacanti (cosine Science & Computing BV)
//  May 15, 2008
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


#include <G4RotationMatrix.hh>

#include <G4XraySpecularReflectingSurface.hh>
#include <G4XrayGrazingAngleScattering.hh>

#include<sstream>

using namespace std;

G4XraySpecularReflectingSurface::G4XraySpecularReflectingSurface(const G4String& name, G4double roughness, G4bool force_unit_reflectivity)
  : G4XraySurfaceProperty(name), m_roughness(roughness), m_force_unit_reflectivity(force_unit_reflectivity) {
  if(m_roughness < 0){
    ostringstream msg;
    msg << "roughness=" << roughness << " must be non negative.";
    G4Exception("4XraySpecularReflect",
		"0",
		JustWarning,
		msg.str().c_str());
  }
}


G4XraySpecularReflectingSurface::G4XraySpecularReflectingSurface(G4double roughness, G4bool force_unit_reflectivity)
  : G4XraySurfaceProperty("G4XraySpecularReflectingSurface"), m_roughness(roughness), m_force_unit_reflectivity(force_unit_reflectivity) {
  if(m_roughness < 0){
    ostringstream msg;
    msg << "roughness=" << roughness << " must be non negative.";
    G4Exception("G4XraySpecularReflectingSurface",
		"0",
		JustWarning,
		msg.str().c_str());
  }
}

double G4XraySpecularReflectingSurface::roughness() const { return m_roughness; }


 G4double G4XraySpecularReflectingSurface::reflect(G4double photonEnergy,
			 G4complex // preRefractionIndex
						   ,
			 G4complex postRefractionIndex,
			 G4ThreeVector & surfaceNormal,
			 G4ThreeVector & photonDirection,
			 const G4ThreeVector & // photonPosition
						   ) const {
  if(photonEnergy <= 0){
    G4Exception("G4XraySpecularReflectingSurface::reflect",
		"0",
		EventMustBeAborted,
		"photon energy is <= 0.");
  }
  
  const G4double cosa = photonDirection * surfaceNormal;
  const G4double grazinga = acos(cosa) - M_PI_2;  
  G4double reflectivity = 0;
  if(m_force_unit_reflectivity)
    reflectivity = 1;
  else {
    G4XrayGrazingAngleScattering::ReflectivityFromVacuum(grazinga, postRefractionIndex, reflectivity);
    
    const G4double specular_fraction = G4XrayGrazingAngleScattering::SpecularFraction(photonEnergy,
										      grazinga, m_roughness,postRefractionIndex);
    reflectivity *= specular_fraction;
  }
  specular_reflect(surfaceNormal, photonDirection);
  return reflectivity;
}

void G4XraySpecularReflectingSurface::specular_reflect(const G4ThreeVector & surfaceNormal,
						       G4ThreeVector & photonDirection) {
  const G4double cosa = photonDirection * surfaceNormal;
  const G4double rotationa = 2*acos(cosa) - M_PI;
  const G4ThreeVector rotax = photonDirection.cross(surfaceNormal);
  const G4RotationMatrix rotmat(rotax, rotationa);
  photonDirection = rotmat*photonDirection;
}




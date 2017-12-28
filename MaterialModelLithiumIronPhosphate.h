// -*- C++ -*-
#ifndef MATERIAL_MODEL_CHEMOMECHANICAL
#define MATERIAL_MODEL_CHEMOMECHANICAL

#include "Definitions.h"

namespace MaterialModels {

class PhaseFieldBatteryModel2D {

/* Implementation of chemomechanical LiFePO4-FePO4 model
 * by N. Nadkarni
 */

    public:

      typedef Matrix<double, 8, 1> Strain;       // contains eps11,eps22,eps12, c, c,1,c,2, mu,1,mu,2
      typedef Matrix<double, 8, 1> Stress;       // contains sig11,sig22,sig12, xi, pi1, pi2, j1,j2
      typedef Matrix<double, 8, 8> TangentMatrix; // gradients of the stress matrix with respect to strain
			typedef EmptyInternalVariables InternalVariables;

      struct MaterialParameters {                 // this container is used to define all the material parameters
        double R, T, cmax, Omega;  // for entropic and mixing free energy
        double kappa; // Cahn-Hilliard phase field parameter
        double C11, C22, C33, C44, C55, C66, C12, C31, C23; // averaged elastic stiffness constants
        double eps0aa, eps0bb, eps0cc; // chemical volumetric strain components
        double Da, Db; // diffusivity in the two directions
      };
      
      
      // component names
      static
        vector<string>
        getStressComponentNames() {

          vector<string> stressComponentNames;
          stressComponentNames.push_back(string("MechanicalStress_XX"));
          stressComponentNames.push_back(string("MechanicalStress_YY"));
          stressComponentNames.push_back(string("MechanicalStress_XY"));
          stressComponentNames.push_back(string("ConcentrationStress"));
          stressComponentNames.push_back(string("ConcentrationGradientStress_1"));
          stressComponentNames.push_back(string("ConcentrationGradientStress_2"));
          stressComponentNames.push_back(string("ChemicalPotentialStress_1"));
          stressComponentNames.push_back(string("ChemicalPotentialStress_2"));

          return stressComponentNames;

        }

      static
        vector<string>
        getStrainComponentNames() {
          vector<string> strainComponentNames;
          strainComponentNames.push_back(string("Strain_XX"));
          strainComponentNames.push_back(string("Strain_YY"));
          strainComponentNames.push_back(string("2*Strain_XY"));
          strainComponentNames.push_back(string("Concentration"));
          strainComponentNames.push_back(string("ConcentrationGradient_1"));
          strainComponentNames.push_back(string("ConcentrationGradient_2"));
          strainComponentNames.push_back(string("ChemicalPotentialGradient_1"));
          strainComponentNames.push_back(string("ChemicalPotentialGradient_2"));

          return strainComponentNames;

        }

      PhaseFieldBatteryModel2D(const MaterialParameters materialParameters) :
        _materialParameters(materialParameters){
        };
	
      // converting long vector of variables into individual components
      void
        convertGradientIntoVectors(const Strain & gradients, Vector3d & strain, Vector1d & concentration,
            Vector2d & concentrationGradient, Vector2d & chemicalPotentialGradient) const {
          for (unsigned int i=0; i<3; i++) {
            strain(i) = gradients(i);
          }
            concentration(0) = gradients(3);
          for (unsigned int i=0; i<2; i++) {
            concentrationGradient(i) = gradients(4+i);
            chemicalPotentialGradient(i) = gradients(6+i);
          }
        }
      
      // converting individual components into a long vector of variabes
      void
        convertVectorsIntoGradient(Strain & gradients, const Vector3d & strain, const Vector1d & concentration,
            const Vector2d & concentrationGradient, const Vector2d & chemicalPotentialGradient) const {
          for (unsigned int i=0; i<3; i++) {
            gradients(i) = strain(i);
          }
            gradients(3) = concentration(0);
          for (unsigned int i=0; i<2; i++) {
            gradients(4+i) = concentrationGradient(i);
            gradients(6+i) = chemicalPotentialGradient(i);
          }
        }


      Stress
        computeStress(const Strain & gradients,
            const InternalVariables & oldIVS,
            const double time) const {   
          ignoreUnusedVariables(oldIVS,time);
          // split up the received vector into its individual components
          Vector3d strains;
          Vector1d concentration;
          Vector2d concentrationGradient, chemicalPotentialGradient;
          convertGradientIntoVectors(gradients, strains, concentration, concentrationGradient, chemicalPotentialGradient);
					
					// material parameters
					double R = _materialParameters.R;
					double T = _materialParameters.T;
					double cmax = _materialParameters.cmax;
					double Omega = _materialParameters.Omega;
					double kappa = _materialParameters.kappa;
					double C11 = _materialParameters.C11;
					double C22 = _materialParameters.C22;
					double C33 = _materialParameters.C33;
					double C44 = _materialParameters.C44;
					double C55 = _materialParameters.C55;
					double C66 = _materialParameters.C66;
					double C12 = _materialParameters.C12;
					double C23 = _materialParameters.C23;
					double C31 = _materialParameters.C31;
					double eps0aa = _materialParameters.eps0aa;
					double eps0bb = _materialParameters.eps0bb;
					double eps0cc = _materialParameters.eps0cc;
					double Da = _materialParameters.Da;
					double Db = _materialParameters.Db;
				
					
					ignoreUnusedVariables(C55);
					ignoreUnusedVariables(C66);
					// now compute the individual components to be returned

          Vector3d stresses;
					stresses.fill(0.);
          Vector1d concentrationStress;
					concentrationStress.fill(0.);
          Vector2d concentrationGradientStress;
					concentrationGradientStress.fill(0.);
          Vector2d chemicalPotentialStress;
					chemicalPotentialStress.fill(0.);
					
          // implementation for the plane strain case (a-b-c) -> (1-2-3) in simulations
            
          double elstrain11 = strains(0) - eps0aa*concentration(0)/cmax;
          double elstrain22 = strains(1) - eps0bb*concentration(0)/cmax;
					double elstrain33 = - eps0cc*concentration(0)/cmax;
          double elstrain12 = strains(2);
					
					//cout << eps0aa << endl;
					//cout << eps0bb << endl;
					//cout << cmax << endl;
					//cout << elstrain33 << endl;
					
					// stresses
					stresses(0) = C11*elstrain11 + C12*elstrain22 + C31*elstrain33;
					stresses(1) = C12*elstrain11 + C22*elstrain22 + C23*elstrain33;
					stresses(2) = 2*C44*elstrain12;
					double stress33 = C31*elstrain11 + C23*elstrain22 + C33*elstrain33;

          // concentration stress
					concentrationStress(0) = R*T*log(concentration(0)/(cmax-concentration(0))) + Omega*(1.-2.*concentration(0)/cmax) - (stresses(0)*eps0aa + stresses(1)*eps0bb + stress33*eps0cc)/cmax - 0.01*R*T*(cmax/concentration(0)) + 0.01*R*T*(cmax/(cmax-concentration(0)));
					
					
					// concentration gradient stresses
					concentrationGradientStress(0) = kappa*concentrationGradient(0);
					concentrationGradientStress(1) = 1.*kappa*concentrationGradient(1);

          // chemical potential gradient stresses or the flux
					chemicalPotentialStress(0) = -Da/(R*T)*concentration(0)*(1-concentration(0)/cmax)*chemicalPotentialGradient(0);
					chemicalPotentialStress(1) = -Db/(R*T)*concentration(0)*(1-concentration(0)/cmax)*chemicalPotentialGradient(1);
					
          Stress stressVector;
          convertVectorsIntoGradient(stressVector, stresses, concentrationStress, concentrationGradientStress, chemicalPotentialStress);

          return stressVector;
        };


      TangentMatrix
        computeTangentMatrix(const Strain & gradients,
            const InternalVariables & oldIVS,
            const double time) const {
					
					ignoreUnusedVariables(oldIVS,time);
					// split into individual vectors
					Vector3d strains;
					Vector1d concentration;
					Vector2d concentrationGradient, chemicalPotentialGradient;
					convertGradientIntoVectors(gradients, strains, concentration, concentrationGradient, chemicalPotentialGradient);
					TangentMatrix tangentMatrix;
          tangentMatrix.setZero();
          // compute the tangent matrix based on the above vectors
          // ...
					
					// implementation for the plane strain case (a-b-c) -> (1-2-3) in simulations
					
					//double elstrain11 = strains(0) - eps0aa*concentration(0)/cmax;
					//double elstrain22 = strains(1) - eps0bb*concentration(0)/cmax;
					//double elstrain33 = - eps0cc*concentration(0)/cmax;
					//double elstrain12 = strains(2);
					
					// material parameters
					double R = _materialParameters.R;
					double T = _materialParameters.T;
					double cmax = _materialParameters.cmax;
					double Omega = _materialParameters.Omega;
					double kappa = _materialParameters.kappa;
					double C11 = _materialParameters.C11;
					double C22 = _materialParameters.C22;
					double C33 = _materialParameters.C33;
					double C44 = _materialParameters.C44;
					double C55 = _materialParameters.C55;
					double C66 = _materialParameters.C66;
					double C12 = _materialParameters.C12;
					double C23 = _materialParameters.C23;
					double C31 = _materialParameters.C31;
					double eps0aa = _materialParameters.eps0aa;
					double eps0bb = _materialParameters.eps0bb;
					double eps0cc = _materialParameters.eps0cc;
					double Da = _materialParameters.Da;
					double Db = _materialParameters.Db;
					
					ignoreUnusedVariables(C55);
					ignoreUnusedVariables(C66);
					
					// stress11
					tangentMatrix(0,0) = C11;
					tangentMatrix(0,1) = C12;
					tangentMatrix(0,3) = -(1/cmax)*(C11*eps0aa+C12*eps0bb+C31*eps0cc);
					
					// stress22
					tangentMatrix(1,0) = C12;
					tangentMatrix(1,1) = C22;
					tangentMatrix(1,3) = -(1/cmax)*(C12*eps0aa+C22*eps0bb+C23*eps0cc);
					
					// stress12
					tangentMatrix(2,2) = 2*C44;
					
					// xi
					tangentMatrix(3,0) = -(1/cmax)*(C11*eps0aa+C12*eps0bb+C31*eps0cc);
					tangentMatrix(3,1) = -(1/cmax)*(C12*eps0aa+C22*eps0bb+C23*eps0cc);
					
					tangentMatrix(3,3) = R*T*cmax/(concentration(0)*(cmax-concentration(0))) - 2.*Omega/cmax + (C11*eps0aa*eps0aa + C22*eps0bb*eps0bb + C33*eps0cc*eps0cc + 2.*C12*eps0aa*eps0bb + 2.*C23*eps0bb*eps0cc + 2.*C31*eps0cc*eps0aa)/(cmax*cmax) + 0.01*R*T*(cmax/pow(concentration(0),2)) + 0.01*R*T*(cmax/pow(cmax-concentration(0),2));
					
					// pi1
					tangentMatrix(4,4) = kappa;
					
					// pi2
					tangentMatrix(5,5) = 1.*kappa;
					
					// j1
					tangentMatrix(6,6) = -Da/(R*T)*concentration(0)*(1-concentration(0)/cmax);
					tangentMatrix(6,3) = -Da/(R*T)*(1-2.*concentration(0)/cmax)*chemicalPotentialGradient(0);
					
					// j2
					tangentMatrix(7,7) = -Db/(R*T)*concentration(0)*(1-concentration(0)/cmax);
					tangentMatrix(7,3) = -Db/(R*T)*(1-2.*concentration(0)/cmax)*chemicalPotentialGradient(1);

          return tangentMatrix;
        };

      InternalVariables
        computeNewInternalVariables(const Strain & strain,
            const InternalVariables & oldIVS,
            const double time) const {
          ignoreUnusedVariables(strain,time);
          return oldIVS;
        };


    private:
      MaterialParameters _materialParameters;

  };
	
class PhaseFieldBatteryModelPlaneStress2D {
		
		/* Implementation of chemomechanical LiFePO4-FePO4 model
		 * by N. Nadkarni
		 */
		
	public:
		
		typedef Matrix<double, 8, 1> Strain;       // contains eps11,eps22,eps12, c, c,1,c,2, mu,1,mu,2
		typedef Matrix<double, 8, 1> Stress;       // contains sig11,sig22,sig12, xi, pi1, pi2, j1,j2
		typedef Matrix<double, 8, 8> TangentMatrix; // gradients of the stress matrix with respect to strain
		typedef EmptyInternalVariables InternalVariables;
		
		struct MaterialParameters {                 // this container is used to define all the material parameters
			double R, T, cmax, Omega;  // for entropic and mixing free energy
			double kappa; // Cahn-Hilliard phase field parameter
			double C11, C22, C33, C44, C55, C66, C12, C31, C23; // averaged elastic stiffness constants
			double eps0aa, eps0bb, eps0cc; // chemical volumetric strain components
			double Da, Db; // diffusivity in the two directions
		};
		
		
		// component names
		static
		vector<string>
		getStressComponentNames() {
			
			vector<string> stressComponentNames;
			stressComponentNames.push_back(string("MechanicalStress_XX"));
			stressComponentNames.push_back(string("MechanicalStress_YY"));
			stressComponentNames.push_back(string("MechanicalStress_XY"));
			stressComponentNames.push_back(string("ConcentrationStress"));
			stressComponentNames.push_back(string("ConcentrationGradientStress_1"));
			stressComponentNames.push_back(string("ConcentrationGradientStress_2"));
			stressComponentNames.push_back(string("ChemicalPotentialStress_1"));
			stressComponentNames.push_back(string("ChemicalPotentialStress_2"));
			
			return stressComponentNames;
			
		}
		
		static
		vector<string>
		getStrainComponentNames() {
			vector<string> strainComponentNames;
			strainComponentNames.push_back(string("Strain_XX"));
			strainComponentNames.push_back(string("Strain_YY"));
			strainComponentNames.push_back(string("2*Strain_XY"));
			strainComponentNames.push_back(string("Concentration"));
			strainComponentNames.push_back(string("ConcentrationGradient_1"));
			strainComponentNames.push_back(string("ConcentrationGradient_2"));
			strainComponentNames.push_back(string("ChemicalPotentialGradient_1"));
			strainComponentNames.push_back(string("ChemicalPotentialGradient_2"));
			
			return strainComponentNames;
			
		}
		
		PhaseFieldBatteryModelPlaneStress2D(const MaterialParameters materialParameters) :
		_materialParameters(materialParameters){
		};
		
		// converting long vector of variables into individual components
		void
		convertGradientIntoVectors(const Strain & gradients, Vector3d & strain, Vector1d & concentration,
															 Vector2d & concentrationGradient, Vector2d & chemicalPotentialGradient) const {
			for (unsigned int i=0; i<3; i++) {
				strain(i) = gradients(i);
			}
			concentration(0) = gradients(3);
			for (unsigned int i=0; i<2; i++) {
				concentrationGradient(i) = gradients(4+i);
				chemicalPotentialGradient(i) = gradients(6+i);
			}
		}
		
		// converting individual components into a long vector of variabes
		void
		convertVectorsIntoGradient(Strain & gradients, const Vector3d & strain, const Vector1d & concentration,
															 const Vector2d & concentrationGradient, const Vector2d & chemicalPotentialGradient) const {
			for (unsigned int i=0; i<3; i++) {
				gradients(i) = strain(i);
			}
			gradients(3) = concentration(0);
			for (unsigned int i=0; i<2; i++) {
				gradients(4+i) = concentrationGradient(i);
				gradients(6+i) = chemicalPotentialGradient(i);
			}
		}
		
		
		Stress
		computeStress(const Strain & gradients,
									const InternalVariables & oldIVS,
									const double time) const {
			ignoreUnusedVariables(oldIVS,time);
			// split up the received vector into its individual components
			Vector3d strains;
			Vector1d concentration;
			Vector2d concentrationGradient, chemicalPotentialGradient;
			convertGradientIntoVectors(gradients, strains, concentration, concentrationGradient, chemicalPotentialGradient);
			
			// material parameters
			double R = _materialParameters.R;
			double T = _materialParameters.T;
			double cmax = _materialParameters.cmax;
			double Omega = _materialParameters.Omega;
			double kappa = _materialParameters.kappa;
			double C11 = _materialParameters.C11;
			double C22 = _materialParameters.C22;
			double C33 = _materialParameters.C33;
			double C44 = _materialParameters.C44;
			double C55 = _materialParameters.C55;
			double C66 = _materialParameters.C66;
			double C12 = _materialParameters.C12;
			double C23 = _materialParameters.C23;
			double C31 = _materialParameters.C31;
			double eps0aa = _materialParameters.eps0aa;
			double eps0bb = _materialParameters.eps0bb;
			double eps0cc = _materialParameters.eps0cc;
			double Da = _materialParameters.Da;
			double Db = _materialParameters.Db;
			
			
			ignoreUnusedVariables(C55);
			ignoreUnusedVariables(C66);
			ignoreUnusedVariables(eps0cc);
			
			// now compute the individual components to be returned
			
			Vector3d stresses;
			stresses.fill(0.);
			Vector1d concentrationStress;
			concentrationStress.fill(0.);
			Vector2d concentrationGradientStress;
			concentrationGradientStress.fill(0.);
			Vector2d chemicalPotentialStress;
			chemicalPotentialStress.fill(0.);
			
			// implementation for the plane strain case (a-b-c) -> (1-2-3) in simulations
			
			double elstrain11 = strains(0) - eps0aa*concentration(0)/cmax;
			double elstrain22 = strains(1) - eps0bb*concentration(0)/cmax;
			double elstrain33 = -(C31*elstrain11 + C23*elstrain22)/C33;
			double elstrain12 = strains(2);
			
			//cout << eps0aa << endl;
			//cout << eps0bb << endl;
			//cout << cmax << endl;
			//cout << elstrain33 << endl;
			
			// stresses
			stresses(0) = C11*elstrain11 + C12*elstrain22 + C31*elstrain33;
			stresses(1) = C12*elstrain11 + C22*elstrain22 + C23*elstrain33;
			stresses(2) = 2*C44*elstrain12;
			
			// concentration stress
			concentrationStress(0) = R*T*log(concentration(0)/(cmax-concentration(0))) + Omega*(1-2*concentration(0)/cmax) - (stresses(0)*eps0aa + stresses(1)*eps0bb)/cmax- 0.01*R*T*(cmax/concentration(0)) + 0.01*R*T*(cmax/(cmax-concentration(0)));
			
			// concentration gradient stresses
			concentrationGradientStress(0) = kappa*concentrationGradient(0);
			concentrationGradientStress(1) = 1.*kappa*concentrationGradient(1);
			
			// chemical potential gradient stresses or the flux
			chemicalPotentialStress(0) = -Da/(R*T)*concentration(0)*(1-concentration(0)/cmax)*chemicalPotentialGradient(0);
			chemicalPotentialStress(1) = -Db/(R*T)*concentration(0)*(1-concentration(0)/cmax)*chemicalPotentialGradient(1);
			
			Stress stressVector;
			convertVectorsIntoGradient(stressVector, stresses, concentrationStress, concentrationGradientStress, chemicalPotentialStress);
			
			return stressVector;
		};
		
		
		TangentMatrix
		computeTangentMatrix(const Strain & gradients,
												 const InternalVariables & oldIVS,
												 const double time) const {
			
			ignoreUnusedVariables(oldIVS,time);
			// split into individual vectors
			Vector3d strains;
			Vector1d concentration;
			Vector2d concentrationGradient, chemicalPotentialGradient;
			convertGradientIntoVectors(gradients, strains, concentration, concentrationGradient, chemicalPotentialGradient);
			TangentMatrix tangentMatrix;
			tangentMatrix.setZero();
			// compute the tangent matrix based on the above vectors
			// ...
			
			// implementation for the plane strain case (a-b-c) -> (1-2-3) in simulations
			
			//double elstrain11 = strains(0) - eps0aa*concentration(0)/cmax;
			//double elstrain22 = strains(1) - eps0bb*concentration(0)/cmax;
			//double elstrain33 = - eps0cc*concentration(0)/cmax;
			//double elstrain12 = strains(2);
			
			// material parameters
			double R = _materialParameters.R;
			double T = _materialParameters.T;
			double cmax = _materialParameters.cmax;
			double Omega = _materialParameters.Omega;
			double kappa = _materialParameters.kappa;
			double C11 = _materialParameters.C11;
			double C22 = _materialParameters.C22;
			double C33 = _materialParameters.C33;
			double C44 = _materialParameters.C44;
			double C55 = _materialParameters.C55;
			double C66 = _materialParameters.C66;
			double C12 = _materialParameters.C12;
			double C23 = _materialParameters.C23;
			double C31 = _materialParameters.C31;
			double eps0aa = _materialParameters.eps0aa;
			double eps0bb = _materialParameters.eps0bb;
			double eps0cc = _materialParameters.eps0cc;
			double Da = _materialParameters.Da;
			double Db = _materialParameters.Db;
			
			ignoreUnusedVariables(C55);
			ignoreUnusedVariables(C66);
			ignoreUnusedVariables(eps0cc);
			
			// stress11
			tangentMatrix(0,0) = C11-C31*C31/C33;
			tangentMatrix(0,1) = C12-C31*C23/C33;
			tangentMatrix(0,3) = -(1/cmax)*(C11*eps0aa+C12*eps0bb-C31*(C31*eps0aa/C33+C23*eps0bb/C33));
			
			// stress22
			tangentMatrix(1,0) = C12-C23*C31/C33;
			tangentMatrix(1,1) = C22-C23*C23/C33;
			tangentMatrix(1,3) = -(1/cmax)*(C12*eps0aa+C22*eps0bb-C23*(C31*eps0aa/C33+C23*eps0bb/C33));
			
			// stress12
			tangentMatrix(2,2) = 2*C44;
			
			// xi
			tangentMatrix(3,0) = (-tangentMatrix(0,0)*eps0aa - tangentMatrix(1,0)*eps0bb)/cmax;
			tangentMatrix(3,1) = (-tangentMatrix(1,0)*eps0aa - tangentMatrix(1,1)*eps0bb)/cmax;
			tangentMatrix(3,3) = R*T*cmax/(concentration(0)*(cmax-concentration(0))) - 2.*Omega/cmax - (tangentMatrix(0,3)*eps0aa + tangentMatrix(1,3)*eps0bb)/cmax + 0.01*R*T*(cmax/pow(concentration(0),2)) + 0.01*R*T*(cmax/pow(cmax-concentration(0),2));
			
			// pi1
			tangentMatrix(4,4) = kappa;
			
			// pi2
			tangentMatrix(5,5) = 1.*kappa;
			
			// j1
			tangentMatrix(6,6) = -Da/(R*T)*concentration(0)*(1-concentration(0)/cmax);
			tangentMatrix(6,3) = -Da/(R*T)*(1-2.*concentration(0)/cmax)*chemicalPotentialGradient(0);
			
			// j2
			tangentMatrix(7,7) = -Db/(R*T)*concentration(0)*(1-concentration(0)/cmax);
			tangentMatrix(7,3) = -Db/(R*T)*(1-2.*concentration(0)/cmax)*chemicalPotentialGradient(1);
			
			return tangentMatrix;
		};
		
		InternalVariables
		computeNewInternalVariables(const Strain & strain,
																const InternalVariables & oldIVS,
																const double time) const {
			ignoreUnusedVariables(strain,time);
			return oldIVS;
		};
		
		
	private:
		MaterialParameters _materialParameters;
		
	};
	
}
#endif // MATERIAL_MODEL_CHEMOMECHANICAL

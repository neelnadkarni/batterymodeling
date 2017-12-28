// -*- C++ -*-
#include "Definitions.h"
#include "Quadrature.h"
#include "Utilities.h"
#include "ChemoElectroMechanicalTriangle.h"
#include "MaterialModelLithiumIronPhosphate.h"
#include "SurfaceEnergyElement.h"
#include "SurfaceFluxElement.h"
#include "Assemblers.h"
#include "MeshUtilities.h"
#include "SolverNewtonRaphson.h"
#include "../../core/PostProcessorVtk.h"
#include <Eigen/LU>
#include <Eigen/SuperLUSupport>

const unsigned int numberOfQuadraturePointsForTriangle                  = 3;
const unsigned int numberOfQuadraturePointsForSurfaceFlux								= 5;
typedef MaterialModels::PhaseFieldBatteryModelPlaneStress2D																																MaterialModel;
typedef MaterialModels::PhaseFieldBatteryModelPlaneStress2D::MaterialParameters																						MaterialParameters;
typedef MaterialModels::EmptyInternalVariables								 InternalVariables;
typedef MaterialModel::Strain                                  Strain;
typedef MaterialModel::Stress																	 Stress;
typedef MaterialModel::TangentMatrix													 TangentMatrix;
typedef Eigen::SparseMatrix<double>                            SparseMatrix;
typedef Eigen::SuperLU<SparseMatrix>                           SuperLUofSparseMatrix;

//static const size_t TotalDegreesOfFreedom = VolumeElement::TotalDegreesOfFreedom;
//static const size_t MechanicalDegreesOfFreedom = VolumeElement::MechanicalDegreesOfFreedom;


int main(int arc, char *argv[]) {

  ignoreUnusedVariables(arc,argv);
	double L_scale = 1.e-10; // m
	double mu_scale = 298*8.3144598; // J/mol
	double F_scale = 96485.33289; // A.s/mol
	double c_scale = 2.29e3; // mol/m^3
	double k_scale = 1.e-2; // A/m^2
	double T_scale = F_scale*c_scale*L_scale/k_scale; // s
	
	// input all material parameters
	MaterialParameters materialParameters;
	materialParameters.R = 8.3144598/mu_scale; // J/(mol.K)
	materialParameters.T = 298; // K
	materialParameters.cmax = 2.29e4/c_scale; // mol/m^3
	materialParameters.Omega = 11.2e3/mu_scale; // J/mol
	materialParameters.kappa = (0.022e-12/2.29e4)/(mu_scale*L_scale*L_scale/c_scale); // (Jm^2/mol)/(mol/m^3)
	materialParameters.C11 = 157.4e9/(mu_scale*c_scale); // Pa
	materialParameters.C22 = 175.8e9/(mu_scale*c_scale); // Pa
	materialParameters.C33 = 154.0e9/(mu_scale*c_scale); // Pa
	materialParameters.C44 = 37.8e9/(mu_scale*c_scale); // Pa
	materialParameters.C55 = 49.05e9/(mu_scale*c_scale); // Pa
	materialParameters.C66 = 51.6e9/(mu_scale*c_scale); // Pa
	materialParameters.C12 = 51.2e9/(mu_scale*c_scale); // Pa
	materialParameters.C23 = 53.25e9/(mu_scale*c_scale); //  Pa
	materialParameters.C31 = 32.7e9/(mu_scale*c_scale); // Pa
	materialParameters.eps0aa = 0.05;
	materialParameters.eps0bb = 0.028;
	materialParameters.eps0cc = -0.025;
	materialParameters.Da = 1.e-18/(L_scale*L_scale/T_scale); // m^2/s
	materialParameters.Db = 1.e-18/(L_scale*L_scale/T_scale); // m^2/s
	
	MaterialModel materialModel(materialParameters);

	Strain strains;
	InternalVariables oldIVs;
	
	for(unsigned int i = 0; i < 8; i++){
			strains(i) = 1.*(rand() % 10);
	}
	
	cout << "The strains are :" << endl << strains << endl;
	
	double perturbation = 1e-7;
	double tolerance = 1e-7;
	double timeStep = 0.;
	
	const Stress F0 = materialModel.computeStress(strains,oldIVs,timeStep);
	const TangentMatrix K0 = materialModel.computeTangentMatrix(strains,oldIVs,timeStep);
	TangentMatrix stiffnessMatrix;
	stiffnessMatrix.fill(0.);
	
	for (unsigned int i = 0; i < 8; i++) {
		for (unsigned int j = 0; j < 8; j++) {
			Strain perturbedStrains = strains;
			perturbedStrains(j) += perturbation;
			const Stress F1 = materialModel.computeStress(perturbedStrains,oldIVs,timeStep);
			
			//cout << "The perturbed stress component " << i << " for a perturbation in strain component " << j << " is:" << F1(i) << endl;
			//cout << "The unperturbed stress component " << i << " for a perturbation in strain component " << j << " is:" << F0(i) << endl << endl;
			
			stiffnessMatrix(i,j) = (F1(i)-F0(i))/perturbation;
			
			cout << "The numerical tangent matrix element (" << i << "," << j << ") is: " << stiffnessMatrix(i,j) << endl;
			cout << "The exact tangent matrix element (" << i << "," << j << ") is: " << K0(i,j) << endl << endl;
		}
	}
	
	
	cout << "computed stiffness matrix is: \n" << stiffnessMatrix << endl;
	
	cout << "stiffness matrix from perturbations is: \n" << K0 << endl;
	
	cout << "difference in the stiffness matrix is: \n" << K0 - stiffnessMatrix << endl;
	
	const double errorStiffness = (stiffnessMatrix - K0).norm()/K0.norm();
	
	cout << "error of method computeTangentMatrix = " << errorStiffness << endl;

	if (errorStiffness >= tolerance || std::isfinite(errorStiffness) == false) {
		cout << "Warning: element derivatives test failed." << endl << endl;
	}
	else{
		cout << "Element test derivatives passed for volume element." << endl;
	}
	
	return 0;

}
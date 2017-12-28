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
typedef MaterialModels::PhaseFieldBatteryModel2D																																MaterialModel;
typedef MaterialModels::PhaseFieldBatteryModel2D::MaterialParameters																						MaterialParameters;
typedef MaterialModels::EmptyInternalVariables												InternalVariables;
typedef Elements::TriangleForBatterySimulations::LinearChemoMechanical<MaterialModel,
                                                                      numberOfQuadraturePointsForTriangle>			VolumeElement;
typedef Elements::TriangleForBatterySimulations::Properties																											VolumeElementProperties;
typedef Elements::SurfaceGammaElement::LinearTwoNodeSurfaceEnergyElement																				SurfaceEnergyElement;
typedef Elements::SurfaceGammaElement::Properties																																SurfaceEnergyElementProperties;
typedef Elements::SurfaceFluxElement::LinearTwoNodeSurfaceFluxElement<numberOfQuadraturePointsForSurfaceFlux>		SurfaceFluxElement;
typedef Elements::SurfaceFluxElement::Properties																																SurfaceFluxElementProperties;
typedef SingleElementMesh<VolumeElement>                             Mesh;
typedef Assemblers::AssemblerChemoMechanicalProblem<VolumeElement,
																										SurfaceEnergyElement,
																										SurfaceFluxElement>																					Assembler;
typedef VolumeElement::Node                                    Node;
typedef VolumeElement::PotentialVector												 PotentialVector;
typedef VolumeElement::Forces																	 Forces;
typedef VolumeElement::StiffnessMatrix												 StiffnessMatrix;
typedef VolumeElement::Point                                   Point;
typedef VolumeElement::PrimitiveVariables											 PrimitiveVariablesForVolumeElement;
typedef SurfaceFluxElement::PrimitiveVariables								 PrimitiveVariablesForSurfaceElement;
typedef VolumeElement::TotalVariableVector										 TotalVariableVector;
typedef VolumeElement::NTNMatrix															 NTNMatrix;
typedef VolumeElement::Stress                                  Stress;
typedef VolumeElement::Strain                                  Strain;
typedef MaterialModel::TangentMatrix													 TangentMatrix;
typedef Eigen::SparseMatrix<double>                            SparseMatrix;
typedef Eigen::SuperLU<SparseMatrix>                           SuperLUofSparseMatrix;
typedef Solvers::NewtonRaphsonEigen<Assembler, SurfaceFluxElement> SolverNewtonRaphson;

static const size_t TotalDegreesOfFreedom = VolumeElement::TotalDegreesOfFreedom;
//static const size_t MechanicalDegreesOfFreedom = VolumeElement::MechanicalDegreesOfFreedom;


int main(int arc, char *argv[]) {

	ignoreUnusedVariables(arc);
	
	// input all material parameters
	MaterialParameters materialParameters;
	materialParameters.R = 8.3144598; // J/(mol.K)
	materialParameters.T = 298; // K
	materialParameters.cmax = 2.29e4; // mol/m^3
	materialParameters.Omega = 11.2e3; // J/mol
	materialParameters.kappa = 0.022; // Jm^2/mol
	materialParameters.C11 = 157.4; // Pa
	materialParameters.C22 = 175.8; // Pa
	materialParameters.C33 = 154.9; // Pa
	materialParameters.C44 = 37.; // Pa
	materialParameters.C55 = 49.05; // Pa
	materialParameters.C66 = 51.6; // Pa
	materialParameters.C12 = 51.2; // Pa
	materialParameters.C23 = 53.25; //  Pa
	materialParameters.C31 = 32.7; // Pa
	materialParameters.eps0aa = 0.05;
	materialParameters.eps0bb = 0.028;
	materialParameters.eps0cc = -0.025;
	materialParameters.Da = 10; // m^2/s
	materialParameters.Db = 10; // m^2/s
	
	MaterialModel materialModel(materialParameters);
	
	double thickness = 1.; // m
	double cmax = 2.29e4; // mol/m^3
	double gamma_ac_FePO4 = 0.;
	double gamma_ac_LiFePO4 = 0.26e6; // J/m^2
	double gamma_bc_FePO4 = 0.;
	double gamma_bc_LiFePO4 = -0.4e6; // J/m^2
	double R = 8.3144598; // J/(mol.K)
	double T = 298; // K
	double F = 96485.33289; // C/mol (Faraday constant)
	double k = 1.e-2; // A/m^2
	
	VolumeElementProperties volElementProperties(thickness,cmax);
	SurfaceEnergyElementProperties surfaceEnergyElementPropertiesAC(thickness,cmax,gamma_ac_FePO4,gamma_ac_LiFePO4);
	SurfaceEnergyElementProperties surfaceEnergyElementPropertiesBC(thickness,cmax,gamma_bc_FePO4,gamma_bc_LiFePO4);
	SurfaceFluxElementProperties surfaceFluxElementProperties(thickness,R,T,F,cmax,k);
	
	const QuadratureRule<2, numberOfQuadraturePointsForTriangle> quadratureRuleForVolume =
	Quadrature::buildSimplicialQuadrature<2, numberOfQuadraturePointsForTriangle>();
	const QuadratureRule<1, numberOfQuadraturePointsForSurfaceFlux> quadratureRuleForSurface =
	Quadrature::buildGaussianQuadrature<1, numberOfQuadraturePointsForSurfaceFlux>();
	
	Node node1;
	Node node2;
	Node node3;
	
	node1._position[0] = 0.;
	node1._position[1] = 0.;
	node2._position[0] = 0.;
	node2._position[1] = 100.;
	node3._position[0] = 100.;
	node3._position[1] = 0.;
	
	node1._id = 0;
	node2._id = 1;
	node3._id = 2;
	
	array<Node, 3> elementNodes;
	elementNodes[0] = node1;
	elementNodes[1] = node2;
	elementNodes[2] = node3;
	
	array<Node, 2> elementNodesForSurface1;
	elementNodesForSurface1[0] = node1;
	elementNodesForSurface1[1] = node2;
	
	array<Node, 2> elementNodesForSurface2;
	elementNodesForSurface2[0] = node2;
	elementNodesForSurface2[1] = node3;
	
	SurfaceEnergyElement element1(elementNodesForSurface1,surfaceEnergyElementPropertiesAC);
	SurfaceEnergyElement element2(elementNodesForSurface2,surfaceEnergyElementPropertiesAC);
	
	const unsigned int totalNumberOfDOFs = 12;
	double perturbation = 1.e-2;
	
	Assembler assembler(3);
	vector<EssentialBoundaryCondition> essentialBCs;
	
	assembler.addElement(element1);
	assembler.addElement(element2);
	
	Eigen::SparseMatrix<double> tangentMatrix(totalNumberOfDOFs+1,totalNumberOfDOFs+1);
	Eigen::SparseMatrix<double> perturbedTangentMatrix(totalNumberOfDOFs+1,totalNumberOfDOFs+1);
	
	VectorXd force(totalNumberOfDOFs+1);
	force.fill(0.);
	
	VectorXd allPrimitives(totalNumberOfDOFs+1);
	VectorXd perturbedPrimitives(totalNumberOfDOFs+1);
	allPrimitives.fill(0.);
	perturbedPrimitives.fill(0.);
	
	for (unsigned int i = 0; i<totalNumberOfDOFs+1; i++) {
		allPrimitives(i) = (rand() % 20000);
		perturbedPrimitives(i) = allPrimitives(i);
	}
	
	VectorXd prim(allPrimitives.size()-1);
	
	for (unsigned int i = 0; i < allPrimitives.size()-1; i++){
		prim(i) = allPrimitives(i);
	}
	
	vector<TotalVariableVector> nodalPrimitives = Utilities::distributeGlobalVectorToLocalVectors<Assembler>(prim);
	
	double boundaryPotential = allPrimitives(totalNumberOfDOFs);
	double time = 100.;
	//MatrixXd tangentMatrix(totalNumberOfDOFs+1,totalNumberOfDOFs+1);
	//VectorXd force(totalNumberOfDOFs+1);
	
	assembler.assembleForceVector(nodalPrimitives, boundaryPotential, time, &force);
	assembler.assembleStiffnessMatrix(nodalPrimitives, boundaryPotential, time,  &tangentMatrix);
	
	//cout << "The determinant of the tangent matrix is: " << tangentMatrix.determinant() << endl;
	
	VectorXd perturbedForce(totalNumberOfDOFs+1);
	perturbedForce.fill(0.);
	
	for(unsigned int i = 0; i < totalNumberOfDOFs+1; i++){
		for (unsigned int j = 0; j < totalNumberOfDOFs+1; j++) {
			
			perturbedPrimitives = allPrimitives;
			perturbedPrimitives(j) += perturbation;
			
			VectorXd perturbedPrim(allPrimitives.size()-1);
			
			for (unsigned int k = 0; k < allPrimitives.size()-1; k++){
				perturbedPrim(k) = perturbedPrimitives(k);
			}
			
			vector<TotalVariableVector> perturbedNodalPrimitives = Utilities::distributeGlobalVectorToLocalVectors<Assembler>(perturbedPrim);
			
			double perturbedBoundaryPotential = perturbedPrimitives(totalNumberOfDOFs);
			
			perturbedForce.fill(0.);
			
			try {
				assembler.assembleForceVector(perturbedNodalPrimitives,perturbedBoundaryPotential,time,&perturbedForce);
			} catch (std::exception & e) {
				errorStatement("exception caught trying to assemble the forces");
				throw;
			}
			
			perturbedTangentMatrix.coeffRef(i,j) = (perturbedForce(i)-force(i))/perturbation;
			
			if (j==totalNumberOfDOFs){
				cout << "The numerical tangent matrix component (" << i << "," << j << ") in this case is: " << perturbedTangentMatrix.coeffRef(i,j) << endl;
				cout << "The exact tangent matrix component (" << i << "," << j << ") in this case is: " << tangentMatrix.coeffRef(i,j) << endl;
			}
			
		}
	}
	
	cout << "The primitives are: " << endl << allPrimitives << endl;
	
	cout << "The tangent matrix is: " << tangentMatrix << endl;
	 
	cout << "The perturbed tangent matrix is: " << endl << perturbedTangentMatrix << endl;
	 
	 /*cout << "The force vector is: " << endl << force << endl;
	 
	 cout << "The primitives are:" << endl << allPrimitives << endl;*/
	
	const double errorStiffness = (tangentMatrix-perturbedTangentMatrix).norm()/tangentMatrix.norm();
	
	cout << "error of method computeStiffnessMatrix = " << errorStiffness << endl;
	double tolerance = 1.e-6;
	
	if (errorStiffness >= tolerance || std::isfinite(errorStiffness) == false) {
		cout << "Warning: element derivatives test failed." << endl << endl;
	}
	else{
		cout << "Element test derivatives passed for assembler." << endl;
	}
	
}
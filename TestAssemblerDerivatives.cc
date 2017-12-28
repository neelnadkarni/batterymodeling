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
SurfaceFluxElement>																						 Assembler;
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
	
	// scales
	double L_scale = 1.e-11; // m
	double mu_scale = 298*8.314; // J/mol
	double F_scale = 96485.33289; // A.s/mol
	double c_scale = 2.29e4; // mol/m^3
	double k_scale = 1.e-2; // A/m^2
	double T_scale = F_scale*c_scale*L_scale/k_scale; // s
	
	// mesh properties
	Mesh mesh;
	array<double, 2> sideLengths;
	sideLengths[0] = 100.e-9/L_scale; // m
	sideLengths[1] = 30.e-9/L_scale; // m
	
	array<size_t, 2> elementsAlongTheLength;
	elementsAlongTheLength[0] = atof(argv[1]);
	elementsAlongTheLength[1] = atof(argv[2]);
	
	MeshUtilities::buildRectangularTriangleMesh(sideLengths, elementsAlongTheLength, &mesh);
	const size_t numberOfNodes = mesh._nodes.size();
	const size_t numberOfElements = mesh._connectivity.size();
	
	Assembler assemble(numberOfNodes);
	vector<EssentialBoundaryCondition> essentialBCs;
	
	// input all material parameters
	MaterialParameters materialParameters;
	materialParameters.R = 8.3144598/mu_scale; // J/(mol.K)
	materialParameters.T = 298; // K
	materialParameters.cmax = 2.29e4/c_scale; // mol/m^3
	materialParameters.Omega = 11.2e3/mu_scale; // J/mol
	materialParameters.kappa = 0.022e-12/(mu_scale*L_scale*L_scale); // Jm^2/mol
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
	materialParameters.Da = 1.e-13/(L_scale*L_scale/T_scale); // m^2/s
	materialParameters.Db = 1.e-18/(L_scale*L_scale/T_scale); // m^2/s
	
	MaterialModel materialModel(materialParameters);
	
	double thickness = 1.; // m
	double cmax = 2.29e4/c_scale; // mol/m^3
	double gamma_ac_FePO4 = 0.;
	double gamma_ac_LiFePO4 = 0.26/(mu_scale*c_scale*L_scale); // J/m^2
	double gamma_bc_FePO4 = 0.;
	double gamma_bc_LiFePO4 = -0.4/(mu_scale*c_scale*L_scale); // J/m^2
	double R = 8.3144598/mu_scale; // J/(mol.K)
	double T = 298; // K
	double F = 96485.33289/F_scale; // C/mol (Faraday constant)
	double k = 1.e-2/k_scale; // A/m^2
	
	VolumeElementProperties volElementProperties(thickness,cmax);
	SurfaceEnergyElementProperties surfaceEnergyElementPropertiesAC(thickness,cmax,gamma_ac_FePO4,gamma_ac_LiFePO4);
	SurfaceEnergyElementProperties surfaceEnergyElementPropertiesBC(thickness,cmax,gamma_bc_FePO4,gamma_bc_LiFePO4);
	SurfaceFluxElementProperties surfaceFluxElementProperties(thickness,R,T,F,cmax,k);
	
	const QuadratureRule<2, numberOfQuadraturePointsForTriangle> quadratureRuleForVolume =
	Quadrature::buildSimplicialQuadrature<2, numberOfQuadraturePointsForTriangle>();
	const QuadratureRule<1, numberOfQuadraturePointsForSurfaceFlux> quadratureRuleForSurface =
	Quadrature::buildGaussianQuadrature<1, numberOfQuadraturePointsForSurfaceFlux>();
	
	// assemble all volume elements
	for (unsigned int index = 0; index < numberOfElements; index++){
		
		array<Node, VolumeElement::NumberOfNodes> elementNodes = MeshUtilities::getNodesFromMesh<VolumeElement>(mesh, index);
		
		assemble.addElement(VolumeElement(elementNodes,volElementProperties,&quadratureRuleForVolume,&materialModel));
	}
	
	// assemble all surface elements
	for (unsigned int index = 0; index < elementsAlongTheLength[0]; index++) {
		
		array<Node, 2> surfaceNodesDown;
		surfaceNodesDown[0] = mesh._nodes[index];
		surfaceNodesDown[1] = mesh._nodes[index+1];
		
		array<Node, 2> surfaceNodesUp;
		surfaceNodesUp[0] = mesh._nodes[elementsAlongTheLength[1]*(elementsAlongTheLength[0]+1)+index];
		surfaceNodesUp[1] = mesh._nodes[elementsAlongTheLength[1]*(elementsAlongTheLength[0]+1)+index+1];
		
		assemble.addElement(SurfaceEnergyElement(surfaceNodesUp,surfaceEnergyElementPropertiesAC));
		assemble.addElement(SurfaceEnergyElement(surfaceNodesDown,surfaceEnergyElementPropertiesAC));
		assemble.addElement(SurfaceFluxElement(surfaceNodesDown,surfaceFluxElementProperties,&quadratureRuleForSurface));
		assemble.addElement(SurfaceFluxElement(surfaceNodesUp,surfaceFluxElementProperties,&quadratureRuleForSurface));
		
	}
	
	for (unsigned int index = 0; index < elementsAlongTheLength[1]; index++) {
		
		array<Node, 2> surfaceNodesLeft;
		surfaceNodesLeft[0] = mesh._nodes[(elementsAlongTheLength[0]+1)*(index)];
		surfaceNodesLeft[1] = mesh._nodes[(elementsAlongTheLength[0]+1)*(index+1)];
		
		array<Node, 2> surfaceNodesRight;
		surfaceNodesRight[0] = mesh._nodes[(elementsAlongTheLength[0]+1)*(index)+elementsAlongTheLength[0]];
		surfaceNodesRight[1] = mesh._nodes[(elementsAlongTheLength[0]+1)*(index+1)+elementsAlongTheLength[0]];
		
		assemble.addElement(SurfaceEnergyElement(surfaceNodesLeft,surfaceEnergyElementPropertiesBC));
		assemble.addElement(SurfaceEnergyElement(surfaceNodesRight,surfaceEnergyElementPropertiesBC));
		
	}
	
	SolverNewtonRaphson solver(10000,1e-4,false);
	
	const unsigned int totalNumberOfDOFs = numberOfNodes*TotalDegreesOfFreedom;
	double boundaryCurrent = 0.1;
	double time = 100.;
	double penalty = 100.;
	double tolerance = 1.e-5;
	double perturbation = 1.e-9;
	
	/*VectorXd allPrimitives(totalNumberOfDOFs+1);
	VectorXd perturbedPrimitives(totalNumberOfDOFs+1);
	allPrimitives.fill(0.5);
	
	for (unsigned int i = 0; i<totalNumberOfDOFs+1; i++) {
		allPrimitives(i) = (rand() % 100000)/100000.;
		perturbedPrimitives(i) = allPrimitives(i);
	}*/
	
	Eigen::SparseMatrix<double> tangentMatrix(totalNumberOfDOFs+1,totalNumberOfDOFs+1);
	
	/*try {
		solver.computeAllStiffnesses(allPrimitives,boundaryCurrent,time,assemble,penalty,&tangentMatrix);
	} catch (std::exception & e) {
		errorStatement("exception caught trying to assemble the stiffness matrix");
		throw;
	}*/
	
	Eigen::SparseMatrix<double> perturbedTangentMatrix(totalNumberOfDOFs+1,totalNumberOfDOFs+1);
	//perturbedTangentMatrix.fill(0.);
	
	VectorXd force(totalNumberOfDOFs+1);
	force.fill(0.);
	
	VectorXd allPrimitives(totalNumberOfDOFs+1);
	VectorXd perturbedPrimitives(totalNumberOfDOFs+1);
	allPrimitives.fill(0.);
	perturbedPrimitives.fill(0.);
	
	for (unsigned int i = 0; i<numberOfNodes; i++) {
		allPrimitives(i*TotalDegreesOfFreedom) = (rand() % 20000)/200.;
		allPrimitives(i*TotalDegreesOfFreedom+1) = (rand() % 20000)/200.;
		allPrimitives(i*TotalDegreesOfFreedom+2) = (rand() % 20000)/20000.;
		allPrimitives(i*TotalDegreesOfFreedom+3) = (rand() % 20000)/2000.;
		perturbedPrimitives(i*TotalDegreesOfFreedom) = allPrimitives(i*TotalDegreesOfFreedom);
		perturbedPrimitives(i*TotalDegreesOfFreedom+1) = allPrimitives(i*TotalDegreesOfFreedom+1);
		perturbedPrimitives(i*TotalDegreesOfFreedom+2) = allPrimitives(i*TotalDegreesOfFreedom+2);
		perturbedPrimitives(i*TotalDegreesOfFreedom+3) = allPrimitives(i*TotalDegreesOfFreedom+3);
	}
	
	allPrimitives(totalNumberOfDOFs) = 10.;
	perturbedPrimitives(totalNumberOfDOFs) = 10.;
	
	VectorXd prim(allPrimitives.size()-1);
	
	for (unsigned int i = 0; i < allPrimitives.size()-1; i++){
		prim(i) = allPrimitives(i);
	}
	
	vector<TotalVariableVector> nodalPrimitives = Utilities::distributeGlobalVectorToLocalVectors<Assembler>(prim);
	
	double boundaryPotential = allPrimitives(totalNumberOfDOFs);
	
	//MatrixXd tangentMatrix(totalNumberOfDOFs+1,totalNumberOfDOFs+1);
	//VectorXd force(totalNumberOfDOFs+1);

	assemble.assembleForceVector(nodalPrimitives, boundaryPotential, time, &force);
	assemble.assembleStiffnessMatrix(nodalPrimitives, boundaryPotential, time,  &tangentMatrix);
	
	//cout << "The determinant of the tangent matrix is: " << tangentMatrix.determinant() << endl;

	VectorXd perturbedForce(totalNumberOfDOFs+1);
	perturbedForce.fill(0.);
	
	double numericalSum = 0.;
	double analyticalSum = 0.;
	
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
				assemble.assembleForceVector(perturbedNodalPrimitives,perturbedBoundaryPotential,time,&perturbedForce);
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
	
	/*cout << "The tangent matrix is: " << tangentMatrix << endl;
	
	cout << "The perturbed tangent matrix is: " << endl << perturbedTangentMatrix << endl;
	
	cout << "The force vector is: " << endl << force << endl;*/
	
	cout << "The primitives are:" << endl << allPrimitives << endl;
	
	const double errorStiffness = (tangentMatrix-perturbedTangentMatrix).norm()/tangentMatrix.norm();
	
	cout << "error of method computeStiffnessMatrix = " << errorStiffness << endl;
	
	
	if (errorStiffness >= tolerance || std::isfinite(errorStiffness) == false) {
		cout << "Warning: element derivatives test failed." << endl << endl;
	}
	else{
		cout << "Element test derivatives passed for assembler." << endl;
	}
	
	// Test for solver

	force.fill(0.);
	
	try {
		solver.computeAllForces(allPrimitives,boundaryCurrent,time,assemble,penalty,&force);
	 } catch (std::exception & e) {
		errorStatement("exception caught trying to assemble the forces");
		throw;
	 }
	
	Eigen::SparseMatrix<double> tangentMatrixSolver(totalNumberOfDOFs+1,totalNumberOfDOFs+1);
	
	Eigen::SparseMatrix<double> perturbedTangentMatrixSolver(totalNumberOfDOFs+1,totalNumberOfDOFs+1);
	
	try {
		solver.computeAllStiffnesses(allPrimitives,boundaryCurrent,time,assemble,penalty,&tangentMatrixSolver);
	} catch (std::exception & e) {
		errorStatement("exception caught trying to assemble the forces");
		throw;
	}
	
	for (unsigned int j = 0; j < totalNumberOfDOFs+1; j++) {
	
		perturbedPrimitives = allPrimitives;
		perturbedPrimitives(j) += perturbation;
	 
		perturbedForce.fill(0.);
		
		try {
			solver.computeAllForces(perturbedPrimitives,boundaryCurrent,time,assemble,penalty,&perturbedForce);
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble the forces");
			throw;
		}
		
		for(unsigned int i = 0; i < totalNumberOfDOFs+1; i++){
		
			perturbedTangentMatrixSolver.coeffRef(i,j) = (perturbedForce(i)-force(i))/perturbation;
			
			cout << "The numerical tangent matrix component (" << i << "," << j << ") in this case is: " << perturbedTangentMatrixSolver.coeffRef(i,j) << endl;
			cout << "The exact tangent matrix component (" << i << "," << j << ") in this case is: " << tangentMatrixSolver.coeffRef(i,j) << endl << endl;
			
			numericalSum += perturbedTangentMatrix.coeffRef(i,j);
			analyticalSum += tangentMatrix.coeffRef(i,j);
		}
		
	}
	
	cout << "The numerical sum is: " << numericalSum << endl;
	cout << "The analytical sum is: " << analyticalSum << endl;
	
	cout << "The force vector is: " << endl << force << endl;
	
	cout << "The stiffness matrix is: " << endl << tangentMatrixSolver << endl;
	
	cout << "The perturbed stiffness matrix is: " << endl << perturbedTangentMatrixSolver << endl;
	
	 /*cout << "The perturbed tangent matrix compiled by the solver is: " << endl << perturbedTangentMatrixSolver << endl;
	 
	 cout << "The force vector compiled by the solver is: " << endl << force << endl;
	 
	 cout << "The primitives are:" << endl << allPrimitives << endl;*/
	 
	 const double errorStiffnessSolver = (perturbedTangentMatrixSolver - tangentMatrixSolver).norm()/tangentMatrixSolver.norm();
	 
	 cout << "error of method computeStiffnessMatrix for the solver is = " << errorStiffnessSolver << endl;
	 
	if (errorStiffnessSolver >= tolerance || std::isfinite(errorStiffnessSolver) == false) {
		cout << "Warning: element derivatives test failed." << endl << endl;
	 }
	 else{
		cout << "Element test derivatives passed for solver assembly." << endl;
	 }
	
	return 0;
	
}
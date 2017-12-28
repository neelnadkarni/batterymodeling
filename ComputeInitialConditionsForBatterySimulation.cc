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
const double cMin = 1.e-9;
typedef MaterialModels::PhaseFieldBatteryModelPlaneStress2D																											MaterialModel;
typedef MaterialModels::PhaseFieldBatteryModelPlaneStress2D::MaterialParameters																	MaterialParameters;
typedef MaterialModels::EmptyInternalVariables												InternalVariables;
typedef Elements::TriangleForBatterySimulations::LinearChemoMechanical<MaterialModel,
                                                                      numberOfQuadraturePointsForTriangle>			VolumeElement;
typedef Elements::TriangleForBatterySimulations::Properties																											VolumeElementProperties;
typedef Elements::SurfaceGammaElement::LinearTwoNodeSurfaceEnergyElement<numberOfQuadraturePointsForSurfaceFlux>																				SurfaceEnergyElement;
typedef Elements::SurfaceGammaElement::Properties																																SurfaceEnergyElementProperties;
typedef Elements::SurfaceFluxElementDoyle::LinearTwoNodeSurfaceFluxElement<numberOfQuadraturePointsForSurfaceFlux>		SurfaceFluxElement;
typedef Elements::SurfaceFluxElement::Properties																																SurfaceFluxElementProperties;
typedef SingleElementMesh<VolumeElement>                             Mesh;
typedef Assemblers::AssemblerChemoMechanicalProblem<VolumeElement,
																										SurfaceEnergyElement,
																										SurfaceFluxElement>																					Assembler;
typedef VolumeElement::Node                                    Node;
typedef VolumeElement::PotentialVector												 PotentialVector;
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
	
	// For 60C, use cscale = 2.29e2 and timestep 0.0001. Seems to work.
	
	
	// scales
	double L_scale = 1.e-10; // m
	double mu_scale = 298*8.3144598; // J/mol
	double F_scale = 96485.33289; // A.s/mol
	double c_scale = 2.29e3; // mol/m^3
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
	
	// input all material parameters
	MaterialParameters materialParameters;
	materialParameters.R = 8.3144598/mu_scale; // J/(mol.K)
	materialParameters.T = 298; // K
	materialParameters.cmax = 2.29e4/c_scale; // mol/m^3
	materialParameters.Omega = 11.2e3/mu_scale; // J/mol
	materialParameters.kappa = (0.022e-12/2.29e4)/(mu_scale*L_scale*L_scale/c_scale); // Jm^2/mol
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
	materialParameters.Db = 1.e-13/(L_scale*L_scale/T_scale); // m^2/s
	
	MaterialModel materialModel(materialParameters);
	
	double cmax = 2.29e4/c_scale; // mol/m^3
	double thickness = 1.; // m
	VolumeElementProperties volElementProperties(thickness,cmax);
	
	const QuadratureRule<2, numberOfQuadraturePointsForTriangle> quadratureRuleForVolume =
	Quadrature::buildSimplicialQuadrature<2, numberOfQuadraturePointsForTriangle>();
	const QuadratureRule<1, numberOfQuadraturePointsForSurfaceFlux> quadratureRuleForSurface =
	Quadrature::buildGaussianQuadrature<1, numberOfQuadraturePointsForSurfaceFlux>();
	
	// assemble all volume elements
	for (unsigned int index = 0; index < numberOfElements; index++){
		
		array<Node, VolumeElement::NumberOfNodes> elementNodes = MeshUtilities::getNodesFromMesh<VolumeElement>(mesh, index);
		
		assemble.addElement(VolumeElement(elementNodes,volElementProperties,&quadratureRuleForVolume,&materialModel));
	}
	
	// solver
	double timeStep = atof(argv[3]);
	double boundaryPotential = 0.;
	double totalVolume = sideLengths[0]*sideLengths[1]*thickness;
	double boundaryCurrentOverF = 0.;//-cmax*totalVolume*10/(3600./T_scale);
	SolverNewtonRaphson solver(10000,1e-4,false,cmax);
	//vector<TotalVariableVector> allPrimitives;
	
	vector<TotalVariableVector> newPrimitives;

	for (unsigned int i = 0; i<numberOfNodes; i++) {
		TotalVariableVector vector;
		vector.fill(0.);
		vector(2) = 0.1*cmax;
		newPrimitives.push_back(vector);
	}
	
	// ******* sample preparation ******* //
	// assign essential boundary conditions for concentration and solve step-wise
	
	vector<EssentialBoundaryCondition> essentialBCs;
	
	essentialBCs.push_back(EssentialBoundaryCondition(0,0,0.));
	essentialBCs.push_back(EssentialBoundaryCondition(0,1,0.));
	essentialBCs.push_back(EssentialBoundaryCondition(elementsAlongTheLength[0],1,0.));
	
	for (unsigned int iteration = 2; iteration < 11; iteration++) {
		
		for (unsigned int st = 1; st < 3; st++){
			
			cout << "for iteration: " << iteration << endl;
			
			double gamma_ac_FePO4 = 0.;
			double gamma_ac_LiFePO4 = 0.26/(mu_scale*c_scale*L_scale); // J/m^2
			double gamma_bc_FePO4 = 0.;
			double gamma_bc_LiFePO4 = -iteration*0.04/(mu_scale*c_scale*L_scale); // J/m^2
			double R = 8.3144598/mu_scale; // J/(mol.K)
			double T = 298; // K
			double F = 96485.33289/F_scale; // C/mol (Faraday constant)
			double k = 1.e-2/k_scale; // A/m^2
			double totalTime = 0.;
			double barrierCoeff = 0.0001;
		
			SurfaceEnergyElementProperties surfaceEnergyElementPropertiesAC(thickness,cmax,gamma_ac_FePO4,gamma_ac_LiFePO4,barrierCoeff);
			SurfaceEnergyElementProperties surfaceEnergyElementPropertiesBC(thickness,cmax,gamma_bc_FePO4,gamma_bc_LiFePO4,barrierCoeff);
			SurfaceFluxElementProperties surfaceFluxElementProperties(thickness,R,T,F,cmax,k);
		
			// assemble all surface elements
			for (unsigned int index = 2; index < elementsAlongTheLength[0]-2; index++) {
			
				array<Node, 2> surfaceNodesDown;
				surfaceNodesDown[0] = mesh._nodes[index];
				surfaceNodesDown[1] = mesh._nodes[index+1];
			
				array<Node, 2> surfaceNodesUp;
				surfaceNodesUp[0] = mesh._nodes[elementsAlongTheLength[1]*(elementsAlongTheLength[0]+1)+index];
				surfaceNodesUp[1] = mesh._nodes[elementsAlongTheLength[1]*(elementsAlongTheLength[0]+1)+index+1];
			
				assemble.addElement(SurfaceEnergyElement(surfaceNodesUp,surfaceEnergyElementPropertiesAC,&quadratureRuleForSurface));
				assemble.addElement(SurfaceEnergyElement(surfaceNodesDown,surfaceEnergyElementPropertiesAC,&quadratureRuleForSurface));
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
			
				assemble.addElement(SurfaceEnergyElement(surfaceNodesLeft,surfaceEnergyElementPropertiesBC,&quadratureRuleForSurface));
				assemble.addElement(SurfaceEnergyElement(surfaceNodesRight,surfaceEnergyElementPropertiesBC,&quadratureRuleForSurface));
			
			}
		
			for (unsigned int index = 0; index < elementsAlongTheLength[1]+1; index++) {
				
				Node surfaceNodeLeft = mesh._nodes[(elementsAlongTheLength[0]+1)*(index)];
				Node surfaceNodeRight = mesh._nodes[(elementsAlongTheLength[0]+1)*(index)+elementsAlongTheLength[0]];
				
				essentialBCs.push_back(EssentialBoundaryCondition(surfaceNodeLeft._id,2,(iteration-1)*0.1*cmax));
				essentialBCs.push_back(EssentialBoundaryCondition(surfaceNodeRight._id,2,(iteration-1)*0.1*cmax));
				
			}
		
			newPrimitives = solver.computeSolution(essentialBCs,newPrimitives,boundaryCurrentOverF,boundaryPotential,timeStep,&assemble,totalTime);
			
			//cout << "The normalized concentration at the boundary now is: " << iteration << endl << endl;
		}
		
		string address = "../../../batterySimulations/initialConditionKappa10_";
		string length = argv[1];
		string space = "_";
		string breadth = argv[2];
		string text = ".txt";
		
		string folder;
		folder.append(address);
		folder.append(length);
		folder.append(space);
		folder.append(breadth);
		folder.append(text);
		
		std::ofstream primitivesFile;
		primitivesFile.open(folder);
		
		primitivesFile << newPrimitives << boundaryPotential;
	}
	
	//cout << "Applying Neumann boundary conditions now." << endl;
	
	/*
	 for (unsigned int index = 0; index < elementsAlongTheLength[1]+1; index++) {
	 
	 Node surfaceNodeLeft = mesh._nodes[(elementsAlongTheLength[0]+1)*(index)];
	 Node surfaceNodeRight = mesh._nodes[(elementsAlongTheLength[0]+1)*(index)+elementsAlongTheLength[0]];
	 
	 essentialBCs.push_back(EssentialBoundaryCondition(surfaceNodeLeft._id,2,0.98*cmax));
	 essentialBCs.push_back(EssentialBoundaryCondition(surfaceNodeRight._id,2,0.98*cmax));
	 
	}
	 for (unsigned int iteration = 1; iteration<30; iteration++) {
		
		vector<EssentialBoundaryCondition> essentialBCs;
		
		essentialBCs.push_back(EssentialBoundaryCondition(0,0,0.));
		essentialBCs.push_back(EssentialBoundaryCondition(0,1,0.));
		essentialBCs.push_back(EssentialBoundaryCondition(elementsAlongTheLength[0],1,0.));
		
		newPrimitives = solver.computeSolution(essentialBCs,newPrimitives,boundaryCurrentOverF,boundaryPotential,timeStep,&assemble,penalty);
		
		cout << "The normalized concentration at the boundary now is: " << iteration << endl << endl;
		
	 }
	 
	 for (unsigned int iteration = 1; iteration<10; iteration++) {
		
		vector<EssentialBoundaryCondition> essentialBCs;
		
		essentialBCs.push_back(EssentialBoundaryCondition(0,0,0.));
		essentialBCs.push_back(EssentialBoundaryCondition(0,1,0.));
		essentialBCs.push_back(EssentialBoundaryCondition(elementsAlongTheLength[0],1,0.));
		
		newPrimitives = solver.computeSolution(essentialBCs,newPrimitives,boundaryCurrent,boundaryPotential,timeStep,&assemble,penalty);
		
	}*/
	
	return 0;

}
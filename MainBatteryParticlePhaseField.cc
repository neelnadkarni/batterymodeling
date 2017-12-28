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
#include <fstream>

const unsigned int numberOfQuadraturePointsForTriangle                  = 1;
const unsigned int numberOfQuadraturePointsForSurfaceFlux								= 3;
typedef MaterialModels::PhaseFieldBatteryModelPlaneStress2D																											MaterialModel;
typedef MaterialModels::PhaseFieldBatteryModelPlaneStress2D::MaterialParameters																	MaterialParameters;
typedef MaterialModels::EmptyInternalVariables												InternalVariables;
typedef Elements::TriangleForBatterySimulations::LinearChemoMechanical<MaterialModel,
                                                                      numberOfQuadraturePointsForTriangle>			VolumeElement;
typedef Elements::TriangleForBatterySimulations::Properties																											VolumeElementProperties;
typedef Elements::SurfaceGammaElement::LinearTwoNodeSurfaceEnergyElement<numberOfQuadraturePointsForSurfaceFlux>																				SurfaceEnergyElement;
typedef Elements::SurfaceGammaElement::Properties																																SurfaceEnergyElementProperties;
typedef Elements::SurfaceFluxElement::LinearTwoNodeSurfaceFluxElement<numberOfQuadraturePointsForSurfaceFlux>		SurfaceFluxElement;
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


int main(int arc, char *argv[]) {

  ignoreUnusedVariables(arc);
	
	// For 60C, use c_scale = 2.29e2 and timestep = 0.0001. Seems to work.
	// For others use c_scale 2.29e3.
	
	// scales
	double kappa_scale = 0.022e-12/2.29e4; // (Jm^2/mol)/(mol/m^3)
	double mu_scale = 298*8.3144598; // J/mol
	double F_scale = 96485.33289; // A.s/mol
	double c_scale = 2.29e3; // mol/m^3
	double k_scale = 1.e-2; // A/m^2
	double L_scale = sqrt(kappa_scale*c_scale/mu_scale); // m
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
	materialParameters.Db = 1.e-13/(L_scale*L_scale/T_scale); // m^2/s
	
	MaterialModel materialModel(materialParameters);
	
	cout << materialParameters.kappa << endl;
	
	double thickness = 1.; // m
	double cmax = 2.29e4/c_scale; // mol/m^3
	double gamma_ac_FePO4 = 0.;
	double gamma_ac_LiFePO4 = 0.26/(mu_scale*c_scale*L_scale); // J/m^2
	double gamma_bc_FePO4 = 0.;
	double gamma_bc_LiFePO4 = -0.4/(mu_scale*c_scale*L_scale); // J/m^2
	double R = 8.3144598/mu_scale; // J/(mol.K)
	double T = 298.; // K
	double F = 96485.33289/F_scale; // C/mol (Faraday constant)
	double k = 1.e-2/k_scale; // A/m^2
	double totalTime = 0.;
	double barrierCoeff = 0.;//0.001;

	VolumeElementProperties volElementProperties(thickness,cmax);
	SurfaceEnergyElementProperties surfaceEnergyElementPropertiesAC(thickness,cmax,gamma_ac_FePO4,gamma_ac_LiFePO4,barrierCoeff);
	SurfaceEnergyElementProperties surfaceEnergyElementPropertiesBC(thickness,cmax,gamma_bc_FePO4,gamma_bc_LiFePO4,barrierCoeff);
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
	
	// intializing the solver with initial condition having surface wetting
	double timeStep = atof(argv[3]);
	double tolerance = 1.e-6;
	SolverNewtonRaphson solver(10000,tolerance,false,cmax);
	vector<TotalVariableVector> newPrimitives;
	double boundaryPotential = 0.;//atof(argv[7]);
	string text = ".txt";

	string address = "../../../batterySimulations/InitialConditions/initialConditionExtractionIsotropic_";
	string length = argv[1];
	string space = "_";
	string breadth = argv[2];

	string folder;
	folder.append(address);
	folder.append(length);
	folder.append(space);
	folder.append(breadth);
	folder.append(text);
	
	std::ifstream initialCondition(folder);
	TotalVariableVector nodalValues;
	nodalValues.fill(0.);
	
	const unsigned int totalNumberOfVariables =  (elementsAlongTheLength[0]+1)*(elementsAlongTheLength[1]+1)*TotalDegreesOfFreedom+2;
	VectorXd fileValues(totalNumberOfVariables);
	fileValues.fill(0.);
	
	unsigned int n = 0;
	while (n < fileValues.size() && initialCondition >> fileValues(n)) {
		++n;
	}
	
	//cout << fileValues << endl;
	
	for (unsigned int i = 0; i < (elementsAlongTheLength[0]+1)*(elementsAlongTheLength[1]+1); i++) {
		nodalValues.fill(0.);
		for (unsigned int j = 0; j < TotalDegreesOfFreedom; j++) {
			nodalValues(j) = fileValues(TotalDegreesOfFreedom*i+1+j);
		}
		newPrimitives.push_back(nodalValues);
	}
	
	boundaryPotential = fileValues(totalNumberOfVariables-1);
	
	/*// initializing the solver with initial condition without surface wetting
	for (unsigned int i = 0; i<numberOfNodes; i++) {
		TotalVariableVector vector;
		vector.fill(0.);
		vector(2) = 0.97*cmax;
		//vector(3) = 1.;
		newPrimitives.push_back(vector);
	}*/
 	
	//cout << " the initial primitives are: " << endl << newPrimitives << endl;
	
	// adding the boundary conditions: No translations and rotations
	
	double cRate = atof(argv[4]);
	double totalVolume = sideLengths[0]*sideLengths[1]*thickness;
	double boundaryCurrentOverF = -cmax*totalVolume*cRate/(3600./T_scale);
	
	// SOLVING THE PROBLEM
	// assign essential boundary conditions for concentration and solve step-wise
	
	cout << "the normalized current is: " << boundaryCurrentOverF << endl;
	
	// initial few timesteps are chosen to be small to stabilize the solution
	vector<EssentialBoundaryCondition> essentialBCs;
	essentialBCs.push_back(EssentialBoundaryCondition(0,0,0.));
	essentialBCs.push_back(EssentialBoundaryCondition(0,1,0.));
	//essentialBCs.push_back(EssentialBoundaryCondition(elementsAlongTheLength[0],0,0.));
	essentialBCs.push_back(EssentialBoundaryCondition(elementsAlongTheLength[0],1,0.));
	
	double smallTimeStep = 0.000001;
	cout << "Stabilizing the solution with a small normalized timestep of " << smallTimeStep << "." << endl;
	
	for (unsigned int iteration = 0; iteration < 3; ++iteration) {
		newPrimitives = solver.computeSolution(essentialBCs,newPrimitives,boundaryCurrentOverF,boundaryPotential,smallTimeStep,&assemble,totalTime);
	}
	
	cout << "Solution stablized. Using larger prescribed timestep to complete computation." << endl;
	
	int iteration = 0;
	
	double cAvg = 0.;
	assemble.getAverageConcentration(newPrimitives,timeStep,&cAvg);

	while (cAvg > (1.-0.97)*cmax) {
		
		iteration += 1;
		newPrimitives = solver.computeSolution(essentialBCs,newPrimitives,boundaryCurrentOverF,boundaryPotential,timeStep,&assemble,totalTime);
		
		cAvg = 0.;
		assemble.getAverageConcentration(newPrimitives,timeStep,&cAvg);
		
		cout << "iteration: " << iteration << endl;
		
		std::ofstream primitivesFile;
			
		string outputAddress = "../../../batterySimulations/ExtractionIsotropic/LimSimulationForCRate_";
		string rate = argv[4];
		string timeSt = "_iteration_";
		string tStep = to_string(iteration);
			
		string outputFolder;
		outputFolder.append(outputAddress);
		outputFolder.append(rate);
		outputFolder.append(timeSt);
		outputFolder.append(tStep);
			
		outputFolder.append(text);
		primitivesFile.open(outputFolder);
		primitivesFile << newPrimitives << boundaryPotential << endl << cAvg << endl << totalTime << endl;
			
		cout << "File saved." << endl;
		cout << "Simulation " << (cAvg/0.95)*100/cmax << "% done." << endl;
		
		if (cAvg < (1.-0.97)*cmax) {
			break;
		}

	}
	
	return 0;

}
// -*- C++ -*-
#ifndef SOLVER_NEWTON_RAPHSON_H
#define SOLVER_NEWTON_RAPHSON_H

#include "Definitions.h"
#include "Utilities.h"
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/SuperLUSupport>
#include <ads/timer.h>

namespace Solvers {

template <class Assembler, class SurfaceElementType2>
class NewtonRaphsonEigen {

  //typedef typename Assembler::MechanicalVector				MechanicalVector;
	typedef typename Assembler::TotalVariableVector			TotalVariableVector;
  //const unsigned int MechanicalDegreesOfFreedom = Assembler::MechanicalDegreesOfFreedom;
	const unsigned int TotalDegreesOfFreedom = Assembler::TotalDegreesOfFreedom;
	

public:

  NewtonRaphsonEigen(const unsigned int maxIterations = 10000,
                     const double tolerance = 1e-4,
                     const bool verbose = false,
										 const double cmax=1) :
    _maxIterations(maxIterations), _tolerance(tolerance), _verbose(verbose), _cmax(cmax){
  }

  vector<TotalVariableVector>
  computeSolution(const vector<EssentialBoundaryCondition> & essentialBCs,
									vector<TotalVariableVector> & initialPrimitives,
									const double boundaryCurrentOverF,
									double & boundaryPotential,
									const double time,
                  Assembler * assembler,
									double & totalTime) const {

    if (_verbose == true) {
      printf("\n\n\nNewton Raphson solver trying to achieve a tolerance of %e in %u "
             "maximum iterations\n", _tolerance, _maxIterations);
    }

    ads::Timer timer;
		
		double penalty = 1.;
    //const unsigned int mechNumberOfDOFs = assembler->getNumberOfNodes() * MechanicalDegreesOfFreedom;
		const unsigned int totalNumberOfDOFs = assembler->getNumberOfNodes() * TotalDegreesOfFreedom;
    VectorXd globalForces(totalNumberOfDOFs+1); ///******** Add 1 here for extra galvanostatic condition
		
    Eigen::SparseMatrix<double> tangentMatrix(totalNumberOfDOFs+1, totalNumberOfDOFs+1);
		
		VectorXd solution(totalNumberOfDOFs+1); // 1 extra for galvanostatic condition
		solution.fill(0.);
		
		double initialBoundaryPotential = boundaryPotential;
		
		VectorXd solutionNew(totalNumberOfDOFs+1); // 1 extra for galvanostatic condition
		solutionNew.fill(0.);
		
    // we start with an initial guess of zero displacements everywhere:
    VectorXd temp(totalNumberOfDOFs);
		temp.fill(0.);
    timer.tic();
    // fill in the initial guess
    for (size_t nodeIndex = 0; nodeIndex < initialPrimitives.size(); ++nodeIndex) {
      for (unsigned int dofIndex = 0; dofIndex < TotalDegreesOfFreedom; ++dofIndex) {
        temp(nodeIndex * TotalDegreesOfFreedom + dofIndex) =
          initialPrimitives[nodeIndex](dofIndex);
      }
    }
		
    if (_verbose == true){
      printf("Time to fill initial guess: %e\n",timer.toc());
    }

    timer.tic();
    // set essential boundary conditions
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      const size_t dofIndex = bc._nodeId * TotalDegreesOfFreedom + bc._coordinate;
      temp(dofIndex) = bc._constraint;
    }
		
    if (_verbose == true){
      printf("Time to set essential boundary conditions: %e\n",timer.toc());
    }

    // make nodal-wise displacements
    vector<TotalVariableVector> nodalPrimitives =
      Utilities::distributeGlobalVectorToLocalVectors<Assembler>(temp);
		
		vector<TotalVariableVector> initialNodalPrimitives =
		Utilities::distributeGlobalVectorToLocalVectors<Assembler>(temp);
		
		VectorXd globalForcesOfPreviousConcentration(totalNumberOfDOFs+1);
		globalForcesOfPreviousConcentration.fill(0.);
		
		try {
			assembler->assembleForceVectorOfPreviousConcentrations(nodalPrimitives,boundaryPotential,time,&globalForcesOfPreviousConcentration);
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble the concentration force vector");
			throw;
		}
		
    timer.tic();
		
		/*----------------------------------------------------------------------------------------*/
		// STEP 1: initialize global forces
		globalForces.fill(0.);

		// STEP 2: assemble forces
		globalForces = computeAllForcesVec(nodalPrimitives,boundaryPotential,boundaryCurrentOverF,time,assembler,penalty);
		
		// STEP 3: subtract earlier concentration forces
		globalForces -= globalForcesOfPreviousConcentration;
		
		// STEP 4: zero out entries of global forces of BCs
		for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
			const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
			const size_t dofIndex = bc._nodeId * TotalDegreesOfFreedom + bc._coordinate;
			globalForces(dofIndex) = 0.;
		}
		
		/*----------------------------------------------------------------------------------------*/
		
		double residue = globalForces.norm();
		
		if (_verbose == true) {
			printf("The initial residue is : %e\n",residue);
		}
		
		cout << "the initial residue is: " << globalForces.norm() << endl;
		/////////////////////////////////////////////
		
		//// ITERATIONS BEGIN HERE //////////////////
		
		/////////////////////////////////////////////
		
		double newTimeStep = time;
		
    // while residue > tolerance

    unsigned int iterationIndex = 0;
    while(residue > _tolerance && iterationIndex < _maxIterations) {
			
			cout << "Using a timestep of: " << newTimeStep << endl;
			
			if (_verbose == true) {
				cout << endl << endl << endl << "NEW ITERATION" << endl;
			}
			
			/*----------------------------------------------------------------------------------------*/
			// STEP 1: initialize global forces
			globalForces.fill(0.);
			
			// STEP 2: assemble forces
			globalForces = computeAllForcesVec(nodalPrimitives,boundaryPotential,boundaryCurrentOverF,newTimeStep,assembler,penalty);
			
			// STEP 3: subtract earlier concentration forces
			globalForces -= globalForcesOfPreviousConcentration;
			
			// STEP 4: zero out entries of global forces of BCs
			for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
				const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
				const size_t dofIndex = bc._nodeId * TotalDegreesOfFreedom + bc._coordinate;
				globalForces(dofIndex) = 0.;
			}
			
			/*----------------------------------------------------------------------------------------*/
			
			
			/*----------------------------------------------------------------------------------------*/
			// ASSEMBLING THE TANGENT MATRIX
			
			timer.tic();
			
			//cout << "computing the stiffness tensor" << endl;
			
			// STEP 1: assemble the global stiffness matrix
			tangentMatrix = computeAllStiffnessesVec(nodalPrimitives,boundaryPotential,boundaryCurrentOverF,newTimeStep,assembler,penalty);
			
			//cout << "stiffness matrix computed" << endl;
			
      // STEP 2: check for and report rows of zeros in the tangent matrix
      {
        timer.tic();
        vector<size_t> rowIndicesToZeroOut;
        //   apply boundary conditions to stiffness matrix
        for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
          const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
          const size_t dofIndex = bc._nodeId * TotalDegreesOfFreedom + bc._coordinate;
          rowIndicesToZeroOut.push_back(dofIndex);
        }
        tangentMatrix.prune(SparseEigenMatrixRowZeroer(rowIndicesToZeroOut));
        for (size_t rowIndex = 0; rowIndex < rowIndicesToZeroOut.size(); rowIndex++){
          // this is an okay coeffRef
          tangentMatrix.coeffRef(rowIndicesToZeroOut[rowIndex],
                                 rowIndicesToZeroOut[rowIndex]) = 1;
        }

      }
			
      {
        timer.tic();
        vector<double> collapseColumnsToOne;
        collapseColumnsToOne.reserve(tangentMatrix.rows());
        for (size_t rowIndex = 0; rowIndex < tangentMatrix.rows(); rowIndex++){
          collapseColumnsToOne.push_back(0);
        }
        for (int colIndex = 0; colIndex < tangentMatrix.outerSize(); ++colIndex){
          for (Eigen::SparseMatrix<double>::InnerIterator
                 iterator(tangentMatrix,colIndex);
               iterator; ++iterator){
            collapseColumnsToOne[iterator.row()] += std::abs(iterator.value());
          }
        }

        vector<size_t> zeroRows;
        for (size_t rowIndex = 0; rowIndex < tangentMatrix.rows(); rowIndex++){
          if (collapseColumnsToOne[rowIndex] == 0){
            zeroRows.push_back(rowIndex);
          }
        }
        if (zeroRows.size() > 0) {
          printf("Error: Newton Raphson solver found that on iteration %u, "
                 "the tangent matrix has %zu rows full of zeros.  The indices are: \n",
                 iterationIndex, zeroRows.size());
          for (size_t i = 0; i < zeroRows.size(); ++i) {
            printf("%zu, ", zeroRows[i]);
          }
          printf("\nSolver stops.\n");
          throwException("there were rows of zeros in the tangent matrix on "
                         "iteration %u", iterationIndex);
        }
        if (_verbose == true){
          printf("Time to find zeros in tangent matrix: %e\n",timer.toc());
        }
      }
      if (_verbose == true){
        printf("Time to apply boundary conditions to tangent matrix: %e\n",
               timer.toc());
      }
      // this seems a lot of work for just the determinant but it seems
      // to be the only way
      timer.tic();
      const Eigen::SuperLU<Eigen::SparseMatrix<double> >
        superLUofTangentMatrix(tangentMatrix);
      const double density =
        tangentMatrix.nonZeros() / float(tangentMatrix.rows() * tangentMatrix.cols());
      if (density > 0.9) {
        fprintf(stderr, "Error: on iteration %4u the density of the tangent "
                "matrix is very high: %%%5.1f\n", iterationIndex,
                100. * density);
      }
      if (_verbose == true) {
        if (density > 0.2) {
          fprintf(stderr, "Warning: on iteration %4u the density of the tangent "
                  "matrix is fairly high: %%%5.1f\n", iterationIndex,
                  100. * density);
        }
        printf("Time to compute LU decomposition: %e\n",timer.toc());
      }
      const Eigen::ComputationInfo infoAboutLUdecomposition = superLUofTangentMatrix.info();
      if (infoAboutLUdecomposition != Eigen::Success){
        errorStatement("LU decomposition failed\n\n");
        printf("The error code is: ");
        if (infoAboutLUdecomposition == Eigen::NumericalIssue){
          printf("Numerical Issue\n");
					exit(1);
        }
        else if (infoAboutLUdecomposition == Eigen::NoConvergence){
          printf("No Convergence\n");
				}
        else if (infoAboutLUdecomposition == Eigen::InvalidInput){
          printf("Invalid Input\n");
        }
      }
      else{
        if (_verbose == true){
          printf("LU decomposition successful!\n");
        }
      }
			
			// COMPUTING THE SOLUTION
			
			// STEP 1: assign nodal primitives to solution
			for (size_t nodeIndex = 0; nodeIndex < nodalPrimitives.size(); ++nodeIndex) {
				for (unsigned int dofIndex = 0; dofIndex < TotalDegreesOfFreedom; ++dofIndex) {
					solution(nodeIndex * TotalDegreesOfFreedom + dofIndex) =
					nodalPrimitives[nodeIndex](dofIndex);
				}
			}
			
			solution(totalNumberOfDOFs) = boundaryPotential;
			
      // STEP 2: update solution
      // solving via a sparse method (if assembler is large then this may help)
			timer.tic();
			
			solutionNew = solution - superLUofTangentMatrix.solve(globalForces);
			
			VectorXd allConcentrations(assembler->getNumberOfNodes());
			allConcentrations.fill(0.);
			
			for(unsigned int i = 0; i < assembler->getNumberOfNodes(); i++){
				allConcentrations(i) = solutionNew(TotalDegreesOfFreedom*i+2);
			}
			
			cout << "The normalized max concentration value is: " << allConcentrations.maxCoeff()/_cmax << endl;
			cout << "The normalized min concentration value is: " << allConcentrations.minCoeff()/_cmax << endl;
			
			while (allConcentrations.maxCoeff() > 1.*_cmax || allConcentrations.minCoeff() < 0.*_cmax) {
				
				newTimeStep = newTimeStep/2.;
				
				cout << "time-step updated to: " << newTimeStep << ". Redoing the simulation with new time-step." << endl;
				
				// STEP 1: assemble the global stiffness matrix
				tangentMatrix = computeAllStiffnessesVec(initialNodalPrimitives,initialBoundaryPotential,boundaryCurrentOverF,newTimeStep,assembler,penalty);
				
				// STEP 2: check for and report rows of zeros in the tangent matrix
				{
					timer.tic();
					vector<size_t> rowIndicesToZeroOut;
					//   apply boundary conditions to stiffness matrix
					for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
						const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
						const size_t dofIndex = bc._nodeId * TotalDegreesOfFreedom + bc._coordinate;
						rowIndicesToZeroOut.push_back(dofIndex);
					}
					tangentMatrix.prune(SparseEigenMatrixRowZeroer(rowIndicesToZeroOut));
					for (size_t rowIndex = 0; rowIndex < rowIndicesToZeroOut.size(); rowIndex++){
						// this is an okay coeffRef
						tangentMatrix.coeffRef(rowIndicesToZeroOut[rowIndex],
																	 rowIndicesToZeroOut[rowIndex]) = 1;
					}
					
				}
				
				{
					timer.tic();
					vector<double> collapseColumnsToOne;
					collapseColumnsToOne.reserve(tangentMatrix.rows());
					for (size_t rowIndex = 0; rowIndex < tangentMatrix.rows(); rowIndex++){
						collapseColumnsToOne.push_back(0);
					}
					for (int colIndex = 0; colIndex < tangentMatrix.outerSize(); ++colIndex){
						for (Eigen::SparseMatrix<double>::InnerIterator
								 iterator(tangentMatrix,colIndex);
								 iterator; ++iterator){
							collapseColumnsToOne[iterator.row()] += std::abs(iterator.value());
						}
					}
					
					vector<size_t> zeroRows;
					for (size_t rowIndex = 0; rowIndex < tangentMatrix.rows(); rowIndex++){
						if (collapseColumnsToOne[rowIndex] == 0){
							zeroRows.push_back(rowIndex);
						}
					}
					if (zeroRows.size() > 0) {
						printf("Error: Newton Raphson solver found that on iteration %u, "
									 "the tangent matrix has %zu rows full of zeros.  The indices are: \n",
									 iterationIndex, zeroRows.size());
						for (size_t i = 0; i < zeroRows.size(); ++i) {
							printf("%zu, ", zeroRows[i]);
						}
						printf("\nSolver stops.\n");
						throwException("there were rows of zeros in the tangent matrix on "
													 "iteration %u", iterationIndex);
					}
					if (_verbose == true){
						printf("Time to find zeros in tangent matrix: %e\n",timer.toc());
					}
				}
				if (_verbose == true){
					printf("Time to apply boundary conditions to tangent matrix: %e\n",
								 timer.toc());
				}
				// this seems a lot of work for just the determinant but it seems
				// to be the only way
				timer.tic();
				const Eigen::SuperLU<Eigen::SparseMatrix<double> >
				superLUofTangentMatrix(tangentMatrix);
				const double density =
				tangentMatrix.nonZeros() / float(tangentMatrix.rows() * tangentMatrix.cols());
				if (density > 0.9) {
					fprintf(stderr, "Error: on iteration %4u the density of the tangent "
									"matrix is very high: %%%5.1f\n", iterationIndex,
									100. * density);
				}
				if (_verbose == true) {
					if (density > 0.2) {
						fprintf(stderr, "Warning: on iteration %4u the density of the tangent "
										"matrix is fairly high: %%%5.1f\n", iterationIndex,
										100. * density);
					}
					printf("Time to compute LU decomposition: %e\n",timer.toc());
				}
				const Eigen::ComputationInfo infoAboutLUdecomposition = superLUofTangentMatrix.info();
				if (infoAboutLUdecomposition != Eigen::Success){
					errorStatement("LU decomposition failed\n\n");
					printf("The error code is: ");
					if (infoAboutLUdecomposition == Eigen::NumericalIssue){
						printf("Numerical Issue\n");
						exit(1);
					}
					else if (infoAboutLUdecomposition == Eigen::NoConvergence){
						printf("No Convergence\n");
					}
					else if (infoAboutLUdecomposition == Eigen::InvalidInput){
						printf("Invalid Input\n");
					}
				}
				else{
					if (_verbose == true){
						printf("LU decomposition successful!\n");
					}
				}
				
				// STEP 1: initialize global forces
				globalForces.fill(0.);
				
				// STEP 2: assemble forces
				globalForces = computeAllForcesVec(initialNodalPrimitives,initialBoundaryPotential,boundaryCurrentOverF,newTimeStep,assembler,penalty);
				
				// STEP 3: subtract earlier concentration forces
				globalForces -= globalForcesOfPreviousConcentration;
				
				// STEP 4: zero out entries of global forces of BCs
				for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
					const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
					const size_t dofIndex = bc._nodeId * TotalDegreesOfFreedom + bc._coordinate;
					globalForces(dofIndex) = 0.;
				}
				
				// COMPUTING THE SOLUTION
				
				// STEP 1: assign nodal primitives to solution
				for (size_t nodeIndex = 0; nodeIndex < nodalPrimitives.size(); ++nodeIndex) {
					for (unsigned int dofIndex = 0; dofIndex < TotalDegreesOfFreedom; ++dofIndex) {
						solution(nodeIndex * TotalDegreesOfFreedom + dofIndex) =
						initialNodalPrimitives[nodeIndex](dofIndex);
					}
				}
				
				solution(totalNumberOfDOFs) = initialBoundaryPotential;
				
				// STEP 2: update solution
				// solving via a sparse method (if assembler is large then this may help)
				timer.tic();
				int count = 0;
				solutionNew = solution - superLUofTangentMatrix.solve(globalForces);
				
				VectorXd allConcentrations(assembler->getNumberOfNodes());
				allConcentrations.fill(0.);
				
				for(unsigned int i = 0; i < assembler->getNumberOfNodes(); i++){
					allConcentrations(i) = solutionNew(TotalDegreesOfFreedom*i+2);
				}
				
				cout << "The normalized max concentration value is: " << allConcentrations.maxCoeff()/_cmax << endl;
				cout << "The normalized min concentration value is: " << allConcentrations.minCoeff()/_cmax << endl;

				if (allConcentrations.maxCoeff() < 1.*_cmax && allConcentrations.minCoeff() > 0.*_cmax) {
					break;
				}
				
			}
			
			solution = solutionNew;
			
      // reassemble primitives and boundary potential from solution
			for (size_t nodeIndex = 0; nodeIndex < nodalPrimitives.size(); ++nodeIndex) {
				for (unsigned int dofIndex = 0; dofIndex < TotalDegreesOfFreedom; ++dofIndex) {
					nodalPrimitives[nodeIndex](dofIndex) =
					solution(nodeIndex * TotalDegreesOfFreedom + dofIndex);
				}
			}
			
			boundaryPotential = solution(totalNumberOfDOFs);
			
			if (_verbose == true) {
        printf("Time to compute LU decomposition: %e\n",timer.toc());
      }

      // COMPUTING GLOBAL FORCES AGAIN
			// STEP 1: zero out the forces first
			globalForces.fill(0);
			
			// STEP 2: assemble forces
			globalForces = computeAllForcesVec(nodalPrimitives,boundaryPotential,boundaryCurrentOverF,newTimeStep,assembler,penalty);

			// STEP 3: subtract earlier concentration forces
			globalForces -= globalForcesOfPreviousConcentration;
			
			// STEP 4: zero out entries of global forces of BCs
			for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
				const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
				const size_t dofIndex = bc._nodeId * TotalDegreesOfFreedom + bc._coordinate;
				globalForces(dofIndex) = 0.;
			}
						
			cout << "global force residue is: " << globalForces.norm() << endl; //globalForces.lpNorm<Infinity>();
			residue = globalForces.norm();
			
      if (_verbose == true) {
				printf("Newton Raphson iteration %4u, residue = %8.3e\n",
               iterationIndex, residue);
      }

      iterationIndex++;
    }
		
		totalTime += newTimeStep;
		
		//////////////// ITERATION ENDS HERE //////////////
		
		if (iterationIndex == _maxIterations) {
      throwException("Error: Newton Raphson solver could not converge in %u "
                     "iterations.\n", _maxIterations);
    }
		
		else {
			if (_verbose == true) {
				printf("Converged to a stable solution! \n");
			}
		}

    assembler->updateInternalVariables(nodalPrimitives, newTimeStep);
		
		//cout << globalForcesOfPreviousConcentration << endl;
		
    return nodalPrimitives;
  }
	
	void
	computeAllForces(const VectorXd primitives,
									 const double boundaryCurrentOverF,
									 const double time,
									 Assembler assembler,
									 const double penalty,
									 VectorXd * forces) const{
		
		ignoreUnusedVariable(penalty);
		
		const unsigned int totalNumberOfDOFs = assembler.getNumberOfNodes() * TotalDegreesOfFreedom;
		
		VectorXd prim(primitives.size()-1);
		
		for (unsigned int i = 0; i < primitives.size()-1; i++){
			prim(i) = primitives(i);
		}
		
		vector<TotalVariableVector> nodalPrimitives = Utilities::distributeGlobalVectorToLocalVectors<Assembler>(prim);
		
		double boundaryPotential = primitives(totalNumberOfDOFs);
		
		VectorXd globalForces(totalNumberOfDOFs+1);
		globalForces.fill(0.);
		
		// STEP 2: compute global forces
		try {
			assembler.assembleForceVector(nodalPrimitives, boundaryPotential, time, &globalForces);
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble the force vector before "
										 "starting solver iterations\n");
			throw;
		}
		
		// STEP 3: compute forces of the surface flux elements
		double constraint = 0.;
		vector<SurfaceElementType2> allSurfaceElements = assembler.getTertiaryPhysicalElements();
		for(unsigned int gbcIndex = 0; gbcIndex < allSurfaceElements.size(); ++gbcIndex) {
			
			array<size_t, 2> nodeIds = allSurfaceElements[gbcIndex].getNodeIds();
			array<TotalVariableVector, 2> primitives;
			
			primitives[0] = nodalPrimitives[nodeIds[0]];
			primitives[1] = nodalPrimitives[nodeIds[1]];
			
			double force = allSurfaceElements[gbcIndex].computeForceForCurrentEquation(primitives,boundaryPotential,time);
			
			constraint += force;
			
		}
		
		constraint -= boundaryCurrentOverF;
		
		globalForces(totalNumberOfDOFs) = constraint;
		
		(*forces) = globalForces;
	}
	
	
	void
	computeAllStiffnesses(const VectorXd & primitives,
												const double boundaryCurrentOverF,
												const double time,
												Assembler assembler,
												const double penalty,
												Eigen::SparseMatrix<double> * tangent) const{
		
		const unsigned int totalNumberOfDOFs = assembler.getNumberOfNodes() * TotalDegreesOfFreedom;
		
		ignoreUnusedVariables(penalty,boundaryCurrentOverF);
		
		VectorXd prim(primitives.size()-1);
		
		for (unsigned int i = 0; i < primitives.size()-1; i++){
			prim(i) = primitives(i);
		}
		
		vector<TotalVariableVector> nodalPrimitives = Utilities::distributeGlobalVectorToLocalVectors<Assembler>(prim);
		
		double boundaryPotential = primitives(totalNumberOfDOFs);
		
		// STEP 1: initialize the tangent matrix
		Eigen::SparseMatrix<double> tangentMatrix(totalNumberOfDOFs+1, totalNumberOfDOFs+1);

		// STEP 2: assembler assembles the problem
		try {
			assembler.assembleStiffnessMatrix(nodalPrimitives, boundaryPotential, time,  &tangentMatrix);
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble the stiffness matrix on "
										 "iteration\n");
			throw;
		}
		
		// STEP 2: apply galvanostatic boundary conditions
		vector<SurfaceElementType2> allSurfaceElements = assembler.getTertiaryPhysicalElements();
		for(unsigned int gbcIndex = 0; gbcIndex < allSurfaceElements.size(); ++gbcIndex) {
			
			array<size_t, 2> nodeIds = allSurfaceElements[gbcIndex].getNodeIds();
			array<TotalVariableVector, 2> primitives;
			
			primitives[0] = nodalPrimitives[nodeIds[0]];
			primitives[1] = nodalPrimitives[nodeIds[1]];
			
			Matrix<double, 5, 1> matrixElements = allSurfaceElements[gbcIndex].computeStiffnessMatrixComponentsForCurrentEquation(primitives,boundaryPotential,time);
			
			Matrix<double, 2, 1> phiStiffness = allSurfaceElements[gbcIndex].computePhiStiffness(primitives,boundaryPotential,time);
			
			tangentMatrix.coeffRef(totalNumberOfDOFs,TotalDegreesOfFreedom*nodeIds[0]+3) += matrixElements(0);
			tangentMatrix.coeffRef(totalNumberOfDOFs,TotalDegreesOfFreedom*nodeIds[1]+3) += matrixElements(1);
			tangentMatrix.coeffRef(totalNumberOfDOFs,totalNumberOfDOFs) += matrixElements(2);
			tangentMatrix.coeffRef(totalNumberOfDOFs,TotalDegreesOfFreedom*nodeIds[0]+2) += matrixElements(3);
			tangentMatrix.coeffRef(totalNumberOfDOFs,TotalDegreesOfFreedom*nodeIds[1]+2) += matrixElements(4);
			
			tangentMatrix.coeffRef(TotalDegreesOfFreedom*nodeIds[0]+3,totalNumberOfDOFs) += phiStiffness(0);
			tangentMatrix.coeffRef(TotalDegreesOfFreedom*nodeIds[1]+3,totalNumberOfDOFs) += phiStiffness(1);
			
		}
		
		for (unsigned int i = 0; i < totalNumberOfDOFs+1; i++) {
			for (unsigned int j = 0; j < totalNumberOfDOFs+1; j++) {
				(*tangent).coeffRef(i,j) = tangentMatrix.coeffRef(i,j);
			}
		}
		
	}
	
	VectorXd
	computeAllForcesVec(const vector<TotalVariableVector> nodalPrimitives,
									 const double boundaryPotential,
									 const double boundaryCurrentOverF,
									 const double time,
									 Assembler * assembler,
									 const double penalty) const{
		
		const unsigned int totalNumberOfDOFs = assembler->getNumberOfNodes() * TotalDegreesOfFreedom;
		
		ignoreUnusedVariable(penalty);
		
		VectorXd globalForces(totalNumberOfDOFs+1);
		globalForces.fill(0.);
		
		// STEP 2: compute global forces
		try {
			assembler->assembleForceVector(nodalPrimitives, boundaryPotential, time, &globalForces);
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble the force vector before "
										 "starting solver iterations\n");
			throw;
		}
		
		// STEP 3: compute forces of the surface flux elements
		double constraint = 0.;
		vector<SurfaceElementType2> allSurfaceElements = assembler->getTertiaryPhysicalElements();
		for(unsigned int gbcIndex = 0; gbcIndex < allSurfaceElements.size(); ++gbcIndex) {
			
			array<size_t, 2> nodeIds = allSurfaceElements[gbcIndex].getNodeIds();
			array<TotalVariableVector, 2> primitives;
			
			primitives[0] = nodalPrimitives[nodeIds[0]];
			primitives[1] = nodalPrimitives[nodeIds[1]];
			
			double force = allSurfaceElements[gbcIndex].computeForceForCurrentEquation(primitives,boundaryPotential,time);
			
			constraint += force;
			
		}
		
		//cout << "the computed current from the boundary elements for this iteration is: " << constraint << endl;
		
		constraint -= boundaryCurrentOverF;
		
		globalForces(totalNumberOfDOFs) = constraint;
		
		return globalForces;
	}
	
	
	Eigen::SparseMatrix<double>
	computeAllStiffnessesVec(const vector<TotalVariableVector> nodalPrimitives,
												const double boundaryPotential,
												const double boundaryCurrentOverF,
												const double time,
												Assembler * assembler,
												const double penalty) const{
		
		const unsigned int totalNumberOfDOFs = assembler->getNumberOfNodes() * TotalDegreesOfFreedom;
		
		ignoreUnusedVariables(penalty,boundaryCurrentOverF);

		// STEP 1: initialize the tangent matrix
		Eigen::SparseMatrix<double> tangentMatrix(totalNumberOfDOFs+1, totalNumberOfDOFs+1);
		
		// STEP 2: assembler assembles the problem
		try {
			assembler->assembleStiffnessMatrix(nodalPrimitives, boundaryPotential, time,  &tangentMatrix);
		} catch (std::exception & e) {
			errorStatement("exception caught trying to assemble the stiffness matrix on "
										 "iteration\n");
			throw;
		}
		
		//cout << "assembled the stiffness matrix by calling function." << endl;
		
		// STEP 2: apply galvanostatic boundary conditions
		vector<SurfaceElementType2> allSurfaceElements = assembler->getTertiaryPhysicalElements();
		for(unsigned int gbcIndex = 0; gbcIndex < allSurfaceElements.size(); ++gbcIndex) {
			
			array<size_t, 2> nodeIds = allSurfaceElements[gbcIndex].getNodeIds();
			array<TotalVariableVector, 2> primitives;
			
			primitives[0] = nodalPrimitives[nodeIds[0]];
			primitives[1] = nodalPrimitives[nodeIds[1]];
			
			Matrix<double, 5, 1> matrixElements = allSurfaceElements[gbcIndex].computeStiffnessMatrixComponentsForCurrentEquation(primitives,boundaryPotential,time);
			
			Matrix<double, 2, 1> phiStiffness = allSurfaceElements[gbcIndex].computePhiStiffness(primitives,boundaryPotential,time);
			
			//cout << "assembling surface element stiffness components." << totalNumberOfDOFs << endl;
			
			//cout << tangentMatrix.coeffRef(totalNumberOfDOFs,0) << endl;
			
			/*tangentMatrix.coeffRef(totalNumberOfDOFs,TotalDegreesOfFreedom*nodeIds[0]+3) += matrixElements(0);
			tangentMatrix.coeffRef(totalNumberOfDOFs,TotalDegreesOfFreedom*nodeIds[1]+3) += matrixElements(1);
			tangentMatrix.coeffRef(totalNumberOfDOFs,totalNumberOfDOFs) += matrixElements(2);
			tangentMatrix.coeffRef(totalNumberOfDOFs,TotalDegreesOfFreedom*nodeIds[0]+2) += matrixElements(3);
			tangentMatrix.coeffRef(totalNumberOfDOFs,TotalDegreesOfFreedom*nodeIds[1]+2) += matrixElements(4);
			
			tangentMatrix.coeffRef(TotalDegreesOfFreedom*nodeIds[0]+3,totalNumberOfDOFs) += phiStiffness(0);
			tangentMatrix.coeffRef(TotalDegreesOfFreedom*nodeIds[1]+3,totalNumberOfDOFs) += phiStiffness(1);
			
		}*/
	
		return tangentMatrix;
	}

	/*
	
	int count = 1;
	
	if (maxC > _cmax || minC < 0.) {
		
		count = count + 1;
		//	cout << "The iteration is: " << count << endl;
		
		for(unsigned int i = 0; i < assembler->getNumberOfNodes(); i++){
			if(solutionNew(TotalDegreesOfFreedom*i+2) < 0.0*_cmax) {
				solutionNew(TotalDegreesOfFreedom*i+2) = abs(solutionNew(TotalDegreesOfFreedom*i+2));
			}
			else if(solutionNew(TotalDegreesOfFreedom*i+2) > 1.0*_cmax) {
				solutionNew(TotalDegreesOfFreedom*i+2) = _cmax - abs(_cmax-solutionNew(TotalDegreesOfFreedom*i+2));
			}
		}
		
		maxC = allConcentrations.maxCoeff();
		minC = allConcentrations.minCoeff();
		
		//	cout << "The max concentration value is: " << maxC << endl;
		//	cout << "The min concentration value is: " << minC << endl;
	 //count += 1;
	 //cout << "internal iteration " << count << endl;
	 //lambda = lambda/2.;
	 const double a = 120.; // smoothening parameter
	 cout << "applying postconditioning to the concentrations with smoothening parameter = " << a << endl;
	 for(unsigned int i = 0; i < assembler->getNumberOfNodes(); i++){
	 double x = 2.*solutionNew(TotalDegreesOfFreedom*i+2)/_cmax-1.;
	 x = (-0.5*log(cosh(a*(1.-x))) + 0.5*log(cosh(a*(1.+x))))/a;
	 solutionNew(TotalDegreesOfFreedom*i+2) = abs(1.e-6*_cmax + (1-2.e-6)*_cmax*(1.+x)/2.);
	 allConcentrations(i) = solutionNew(TotalDegreesOfFreedom*i+2);
	 }
	 
	 //solutionNew = solution - lambda*superLUofTangentMatrix.solve(globalForces);
	 
	 //VectorXd allConcentrations(assembler->getNumberOfNodes());
	 //allConcentrations.fill(0.);
	 
	 //for(unsigned int i = 0; i < assembler->getNumberOfNodes(); i++){
	 //allConcentrations(i) = solutionNew(TotalDegreesOfFreedom*i+2);
	 //}
	 
	 cout << "The new normalized max concentration value is: " << allConcentrations.maxCoeff()/_cmax << endl;
	 cout << "The new normalized min concentration value is: " << allConcentrations.minCoeff()/_cmax << endl;
	 
	 //if (allConcentrations.maxCoeff() < 1.*_cmax && allConcentrations.minCoeff() > 0.*_cmax){
	 //break;
	 //}
		
	}*/
	
private:
  const unsigned int _maxIterations;
  const double _tolerance;
  const bool _verbose;
	const double _cmax;
};
}

#endif  // SOLVER_NEWTON_RAPHSON_H

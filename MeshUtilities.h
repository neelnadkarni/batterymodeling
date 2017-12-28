// -*- C++ -*-
#ifndef MESH_UTILITIES_H
#define MESH_UTILITIES_H
#include "Definitions.h"
#include "Utilities.h"

#include "geom/mesh/iss/build.h"
#include <geom/mesh/iss/IndSimpSetIncAdj.h>
#include <geom/mesh/iss/IndSimpSet.h>

namespace MeshUtilities {

enum IndexingStyle {ZeroIndexed, OneIndexed};

Eigen::Matrix<double,3,3>
generateRotationMatrix(const double thetaRadians,
                       const Eigen::Vector3d  & axisOfRotation){

  const double sinTheta = sin(thetaRadians);
  const double cosTheta = cos(thetaRadians);

  Matrix<double,3,3> axisCrossProductMatrix;
  axisCrossProductMatrix.fill(0);
  axisCrossProductMatrix(0,1) = -axisOfRotation(2);
  axisCrossProductMatrix(0,2) =  axisOfRotation(1);
  axisCrossProductMatrix(1,0) =  axisOfRotation(2);
  axisCrossProductMatrix(1,2) = -axisOfRotation(0);
  axisCrossProductMatrix(2,0) = -axisOfRotation(1);
  axisCrossProductMatrix(2,1) =  axisOfRotation(0);

  Matrix<double,3,3> axisTensorProduct;
  axisTensorProduct.fill(0);
  axisTensorProduct(0,0) = axisOfRotation(0) * axisOfRotation(0);
  axisTensorProduct(0,1) = axisOfRotation(0) * axisOfRotation(1);
  axisTensorProduct(0,2) = axisOfRotation(0) * axisOfRotation(2);
  axisTensorProduct(1,0) = axisTensorProduct(0,1);
  axisTensorProduct(1,1) = axisOfRotation(1) * axisOfRotation(1);
  axisTensorProduct(1,2) = axisOfRotation(1) * axisOfRotation(2);
  axisTensorProduct(2,0) = axisTensorProduct(0,2);
  axisTensorProduct(2,1) = axisTensorProduct(1,2);
  axisTensorProduct(2,2) = axisOfRotation(2) * axisOfRotation(2);

  return cosTheta * Matrix<double,3,3>::Identity() + sinTheta * axisCrossProductMatrix +
    (1 - cosTheta) * axisTensorProduct;
}

template<class Mesh>
void
rotate2DMesh(const double thetaRadians,
             Mesh * mesh){

  const size_t numberOfNodes = mesh->_nodes.size();

  Eigen::Matrix<double,2,2> rotationMatrix;
  rotationMatrix << cos(thetaRadians), -sin(thetaRadians), sin(thetaRadians), cos(thetaRadians);

  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    (*mesh)._nodes[nodeIndex]._position = rotationMatrix * (*mesh)._nodes[nodeIndex]._position;
  }
}

template<class Mesh>
void
rotate3DMesh(const double thetaRadians,
             const Eigen::Vector3d & axisOfRotation,
             Mesh * mesh){

  const size_t numberOfNodes = mesh->_nodes.size();

  const Eigen::Matrix<double,3,3> rotationMatrix = generateRotationMatrix(thetaRadians,
                                                                          axisOfRotation);
  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    (*mesh)._nodes[nodeIndex]._position = rotationMatrix * (*mesh)._nodes[nodeIndex]._position;
  }
}


template<class Mesh>
void
checkForSpatiallyDuplicateNodes(const double spatialTolerance,
                                const Mesh & mesh){
  Eigen::Matrix<double,Mesh::SpatialDimension,1> node0Position;
  Eigen::Matrix<double,Mesh::SpatialDimension,1> node1Position;

  vector<array<Eigen::Matrix<double,Mesh::SpatialDimension,1>,2>> duplicateNodePositions;

  vector<array<size_t,2>> nodeIdsWithSameSpatialLocation;

  for (size_t node0Index = 0; node0Index < mesh._nodes.size() - 1; node0Index++){
    node0Position = mesh._nodes[node0Index]._position;
    for (size_t node1Index = node0Index + 1; node1Index < mesh._nodes.size(); node1Index++){
      node1Position = mesh._nodes[node1Index]._position;
      if ( (node0Position - node1Position).norm() < spatialTolerance){
        nodeIdsWithSameSpatialLocation.push_back(array<size_t,2>{{node0Index,node1Index}});
        duplicateNodePositions.push_back
          (array<Eigen::Matrix<double,Mesh::SpatialDimension,1>,2>{{node0Position,node1Position}});
      }
    }
  }
  if (nodeIdsWithSameSpatialLocation.size() > 0){
    errorStatement("found nodes which are spatially duplicate\n");
    for (size_t index = 0; index < nodeIdsWithSameSpatialLocation.size(); index++){
      for (size_t nodeIndex = 0; nodeIndex < 2; nodeIndex++){
        printf("%8lu: ",nodeIdsWithSameSpatialLocation[index][nodeIndex]);
        for (size_t spatialIndex = 0; spatialIndex < Mesh::SpatialDimension; spatialIndex++){
          printf("%1.5e, ", duplicateNodePositions[index][nodeIndex](spatialIndex));
        }
        printf("\n");

      }
      printf("Distance = %e\n",(duplicateNodePositions[index][0] -
                                duplicateNodePositions[index][1]).norm());
    }
    throwException("found nodes which are spatially duplicate\n");
  }
}

template<class Mesh, unsigned int SpatialDimension>
void
findBoundingBoxOfGeneralMesh(const Mesh & mesh,
                             Matrix<double,SpatialDimension,1> * minVector,
                             Matrix<double,SpatialDimension,1> * maxVector){
  if (mesh._nodes.size() == 0){
    errorStatement("Mesh contains zero nodes\n");
    exit(1);
  }

  const size_t numberOfNodes = mesh._nodes.size();
  const size_t spatialDimensions = mesh._nodes[0]._position.size();

  VectorXd localMinVector = mesh._nodes[0]._position;
  VectorXd localMaxVector = mesh._nodes[0]._position;
  for (size_t nodeIndex = 1; nodeIndex < numberOfNodes; nodeIndex++){
    const VectorXd nodePosition = mesh._nodes[nodeIndex]._position;
    for (size_t spatialIndex = 0; spatialIndex < spatialDimensions; spatialIndex++){
      if (nodePosition(spatialIndex) > localMaxVector(spatialIndex)){
        localMaxVector(spatialIndex) = nodePosition(spatialIndex);
      }
      if (nodePosition(spatialIndex) < localMinVector(spatialIndex)){
        localMinVector(spatialIndex) = nodePosition(spatialIndex);
      }
    }
  }
  *minVector = localMinVector;
  *maxVector = localMaxVector;
}

template<class Element>
size_t
findIndexOfNodeClosestToAPosition(const typename Element::Point desiredPosition,
                                  const vector<typename Element::Node> & allNodes){
  if (allNodes.size() == 0){
    throwException("there are no nodes compared to desired position\n");
  }
  double smallestDistanceToPosition = (allNodes[0]._position - desiredPosition).norm();
  size_t indexOfClosestNode = 0;
  for (size_t nodeIndex = 0; nodeIndex < allNodes.size(); nodeIndex++){
    const typename Element::Point nodePosition = allNodes[nodeIndex]._position;
    if ((nodePosition - desiredPosition).norm() < smallestDistanceToPosition){
      smallestDistanceToPosition = (nodePosition - desiredPosition).norm();
      indexOfClosestNode = nodeIndex;
    }
  }
  return indexOfClosestNode;
}

template < unsigned int SpatialDimension, class Node>
bool
computeIfNodesArePeriodicPair(Node node1, Node node2,
                              vector<Matrix<double,SpatialDimension,1>> extrema){
  typedef Matrix<double,SpatialDimension,1> Point;

  bool isPair = false;
  const double tol = ( pow(extrema[1](2) - extrema[0](2),2) +
                       pow(extrema[1](1) - extrema[0](1),2) +
                       pow(extrema[1](0) - extrema[0](0),2) ) / double(1e4);

  const Point point1 = node1._position;
  const Point point2 = node2._position;

  if ( ( ( std::abs(point1(0) - point2(0)) < tol) &&
         ( std::abs(point1(1) - point2(1)) < tol) &&
         ( std::abs( std::abs(point1(2) - point2(2)) -
                     std::abs(extrema[1](2) - extrema[0](2))) < tol) ) ) {
    isPair=true;
  }
  if ( ( (std::abs( point1(1)-point2(1))<tol) &&
         (std::abs( point1(2)-point2(2))<tol) &&
         (std::abs( std::abs(point1(0)-point2(0)) -
                    std::abs(extrema[1](0) - extrema[0](0))) < tol) ) ) {
    isPair=true;
  }
  if ( ( (std::abs(point1(0)-point2(0))<tol) &&
         (std::abs(point1(2)-point2(2))<tol) &&
         (std::abs( std::abs(point1(1)-point2(1)) -
                    std::abs(extrema[1](1) - extrema[0](1)))<tol) ) ) {
    isPair=true;
  }
  return isPair;
};

template <unsigned int SpatialDimension, class Node,class Mesh>
bool
computeIfNodesArePeriodicPair(Node node1, Node node2,
                              Mesh &mesh){
  Vector3d minVector;
  Vector3d maxVector;
  findBoundingBoxOfGeneralMesh<Mesh,SpatialDimension>(mesh,&minVector,&maxVector);
  vector<Matrix<double,SpatialDimension,1>> extrema;
  extrema.push_back(minVector);
  extrema.push_back(maxVector);
  return computeIfNodesArePeriodicPair<SpatialDimension>(node1,node2,extrema);

}
template <unsigned int SpatialDimension, class Mesh>
vector<array<size_t,2>> computeNodalPairsForPeriodicBCs(Mesh & mesh){

  const unsigned int numberOfNodes = mesh._nodes.size();
  vector<array<size_t,2>> idsOfNodePairs;
  for (unsigned int node1Index = 0; node1Index < numberOfNodes; node1Index++){
    for (unsigned int node2Index = node1Index; node2Index < numberOfNodes; node2Index++){
      const bool isPair = computeIfNodesArePeriodicPair<SpatialDimension>
        (mesh._nodes[node1Index],
         mesh._nodes[node2Index],
         mesh);
      if (isPair && node1Index != node2Index){
        idsOfNodePairs.push_back(array<size_t,2>
                                 {{mesh._nodes[node1Index]._id,
                                       mesh._nodes[node2Index]._id}});
      }
    }
  }
  return idsOfNodePairs;
}


template <class Element, class Mesh>
void
uniformRefineBarMesh(Mesh * mesh, size_t divisions){
  static const unsigned int D = Element::SpatialDimension;
  typedef Matrix<double, D ,1> Point;

  const size_t numberOfNodes = mesh->_nodes.size();
  const size_t numberOfElements = mesh->_connectivity.size();

  size_t node1Id = 0;
  size_t node2Id = 0;

  size_t availableNodeNumber = numberOfNodes;

  for (size_t elementIndex = 0; elementIndex < numberOfElements; elementIndex++){

    // grab the two nodes of each element
    node1Id = mesh->_connectivity[elementIndex][0];
    node2Id = mesh->_connectivity[elementIndex][1];

    Point node1Position = mesh->_nodes[node1Id]._position;
    Point node2Position = mesh->_nodes[node2Id]._position;

    // add new nodes
    Point positionDifference = node2Position - node1Position;
    Point positionIncrement = positionDifference / double(divisions);

    if (divisions == 1){
      errorStatement("WARNING: You specified a refinement of 1 which does nothing\n");
      break;
    }
    else if (divisions < 1){
      errorStatement("ERROR: You specified less than 1 divisions. Please enter an integer greater than 1\n");
      exit(1);
    }

    size_t newElementLeftNodeId = node1Id;
    size_t newElementRightNodeId = availableNodeNumber;

    array<size_t,2> connection;

    for (size_t addIndex = 1; addIndex < divisions ; addIndex++){
      // create new node
      Point point = node1Position + positionIncrement * addIndex;
      mesh->_nodes.push_back(typename Element::Node(availableNodeNumber, point));

      // add the connections between the new node(s) and to the left
      connection[0] = newElementLeftNodeId;
      connection[1] = newElementRightNodeId;
      mesh->_connectivity.push_back(connection);
      availableNodeNumber++;

      // if the last connection is to be made use the newElementRightNodeId
      if (addIndex == divisions - 1){
        newElementLeftNodeId = newElementRightNodeId;
        newElementRightNodeId = node2Id;
        connection[0] = newElementLeftNodeId;
        connection[1] = newElementRightNodeId;
        mesh->_connectivity.push_back(connection);
      }
      // if not the last element then properly update the left and right node Id
      else {
        newElementLeftNodeId = newElementRightNodeId;
        newElementRightNodeId = availableNodeNumber;
      }

    }
  }
  if (divisions > 1){
    // delete all the original elements / connectivities
    mesh->_connectivity.erase(mesh->_connectivity.begin(),mesh->_connectivity.begin() + numberOfElements);
  }
}

template <unsigned D, class Mesh>
void
uniformRefineBarMeshTwoElement(Mesh * mesh, size_t divisions){
  // this assumes the primary element is the beam
  typedef Matrix<double, D ,1> Point;

  const size_t numberOfNodes = mesh->_nodes.size();
  const size_t numberOfElements = mesh->_connectivity1.size();

  size_t node1Id = 0;
  size_t node2Id = 0;

  size_t availableNodeNumber = numberOfNodes;

  for (size_t elementIndex = 0; elementIndex < numberOfElements; elementIndex++){

    // grab the two nodes of each element
    node1Id = mesh->_connectivity1[elementIndex][0];
    node2Id = mesh->_connectivity1[elementIndex][1];

    Point node1Position = mesh->_nodes[node1Id]._position;
    Point node2Position = mesh->_nodes[node2Id]._position;

    // add new nodes
    Point positionDifference = node2Position - node1Position;
    Point positionIncrement = positionDifference / double(divisions);

    if (divisions == 1){
      errorStatement("WARNING: You specified a refinement of 1 which does nothing\n");
      break;
    }
    else if (divisions < 1){
      errorStatement("ERROR: You specified less than 1 divisions. Please enter an integer greater than 1\n");
      exit(1);
    }

    size_t newElementLeftNodeId = node1Id;
    size_t newElementRightNodeId = availableNodeNumber;

    array<size_t,2> connection;

    for (size_t addIndex = 1; addIndex < divisions ; addIndex++){
      // create new node
      Point point = node1Position + positionIncrement * addIndex;
      mesh->_nodes.push_back(NodeWithId<Point>(availableNodeNumber, point));

      // add the connections between the new node(s) and to the left
      connection[0] = newElementLeftNodeId;
      connection[1] = newElementRightNodeId;
      mesh->_connectivity1.push_back(connection);
      availableNodeNumber++;

      // if the last connection is to be made use the newElementRightNodeId
      if (addIndex == divisions - 1){
        newElementLeftNodeId = newElementRightNodeId;
        newElementRightNodeId = node2Id;
        connection[0] = newElementLeftNodeId;
        connection[1] = newElementRightNodeId;
        mesh->_connectivity1.push_back(connection);
      }
      // if not the last element then properly update the left and right node Id
      else {
        newElementLeftNodeId = newElementRightNodeId;
        newElementRightNodeId = availableNodeNumber;
      }

    }
  }
  if (divisions > 1){
    // delete all the original elements / connectivities
    mesh->_connectivity1.erase(mesh->_connectivity1.begin(),
                               mesh->_connectivity1.begin() + numberOfElements);
  }
}

template <unsigned D, class Mesh>
void
refineBarMeshToEqualizeElementLengths(const double maximumAcceptableLengthRatio,
                                      const unsigned int maximumNumberOfDivisions,
                                      Mesh * mesh,
                                      vector<double> * equalizedElementLengths){

  typedef Matrix<double, D, 1> Point;

  const size_t numberOfNodesInitially = mesh->_nodes.size();
  const size_t numberOfElementsInitially = mesh->_connectivity.size();

  if (numberOfElementsInitially == 0){
    errorStatement("Number of elements = 0\n");
  }

  vector<double> binsOfUniqueLength;
  vector<vector<size_t>> binsOfElementIndices;
  ////// sort the bars into bins //////
  for (size_t elementIndex = 0; elementIndex < numberOfElementsInitially; elementIndex++){
    const size_t node1Id = mesh->_connectivity[elementIndex][0];
    const size_t node2Id = mesh->_connectivity[elementIndex][1];
    const Point node1Position = mesh->_nodes[node1Id]._position;
    const Point node2Position = mesh->_nodes[node2Id]._position;
    const double lengthOfBar = (node1Position - node2Position).norm();

    // if the bins are empty
    if (binsOfUniqueLength.size() == 0){
      binsOfUniqueLength.push_back(lengthOfBar);
      binsOfElementIndices.push_back(vector<size_t>{elementIndex});
    }

    // need to search through bins and see if there are any other lengths
    //  within the acceptable range
    else{
      for (size_t binIndex = 0; binIndex < binsOfUniqueLength.size(); binIndex++){
        const double testBarLengthRatio =
          std::abs((lengthOfBar - binsOfUniqueLength[binIndex])/binsOfUniqueLength[binIndex]);
        if (testBarLengthRatio <= maximumAcceptableLengthRatio){
          binsOfElementIndices[binIndex].push_back(elementIndex);
          break;
        }

        // couldn't find a beam with the same length
        else if (testBarLengthRatio > maximumAcceptableLengthRatio &&
                 binIndex == binsOfUniqueLength.size() - 1){
          binsOfUniqueLength.push_back(lengthOfBar);
          binsOfElementIndices.push_back(vector<size_t>{elementIndex});
        }
      }
    }
  }

  /////// perform iterative scheme to find divisions ///////
  vector<size_t> divisionsForEachBin;
  divisionsForEachBin.reserve(binsOfUniqueLength.size());
  vector<double> dividedLengths;
  dividedLengths.reserve(binsOfUniqueLength.size());
  for (size_t binIndex = 0; binIndex < binsOfUniqueLength.size(); binIndex++){
    divisionsForEachBin.push_back(1);
    dividedLengths.push_back(binsOfUniqueLength[binIndex]);
  }

  while (binsOfUniqueLength.size() > 1 &&
         *max_element(dividedLengths.begin(), dividedLengths.end()) /
         *min_element(dividedLengths.begin(), dividedLengths.end()) >
         maximumAcceptableLengthRatio &&
         *max_element(divisionsForEachBin.begin(), divisionsForEachBin.end()) <
         maximumNumberOfDivisions){

    // find the bin with the largest length and increase the divisor
    const std::vector<double>::iterator binIteratorWithGreatestLength =
      max_element(dividedLengths.begin(), dividedLengths.end());
    const size_t binIndexWithGreatestLength =
      distance(dividedLengths.begin(),binIteratorWithGreatestLength);
    divisionsForEachBin[binIndexWithGreatestLength]++;

    for (size_t binIndex = 0; binIndex < dividedLengths.size(); binIndex++){
      dividedLengths[binIndex] = binsOfUniqueLength[binIndex] /
        double(divisionsForEachBin[binIndex]);
    }
  }

  *equalizedElementLengths = dividedLengths;
  /////////////// TODO FINISH THIS ///////////////
  for (size_t i = 0; i < 100; i++){
    errorStatement("function not completed!\n");
  }
}


template <class Element>
void
checkMeshForErrors(const SingleElementMesh<Element> & mesh){

  const size_t numberOfNodes = mesh._nodes.size();

  Eigen::Matrix<double,Element::SpatialDimension,1> nodePosition =
    mesh._nodes[0]._position;
  for (size_t nodeIndex = 0; nodeIndex < mesh._nodes.size(); ++nodeIndex){
    nodePosition = mesh._nodes[nodeIndex]._position;
    for (size_t spatialIndex = 0; spatialIndex < Element::SpatialDimension; spatialIndex++){
      if (std::isfinite(nodePosition(spatialIndex)) == false){
        throwException("node (%lu) has a non-finite position number"
                       " in spatial index %lu",nodeIndex,spatialIndex);
      }
    }
  }

  for (size_t elementIndex = 0; elementIndex < mesh._connectivity.size(); ++elementIndex){
    for (size_t nodeIndex = 0; nodeIndex < Element::NumberOfNodes; ++nodeIndex){
      const size_t requestedNodeIndex = mesh._connectivity[elementIndex][nodeIndex];
      if (requestedNodeIndex > numberOfNodes - 1){
        throwException("element (%lu) is asking for a node that doesn't exist (%lu)",
                       elementIndex,requestedNodeIndex);
      }
    }
  }
}


template <class Element>
void
readMeshFromFile(const string filename,
                 SingleElementMesh<Element> * mesh,
                 IndexingStyle indexingStyle = ZeroIndexed) {
  typedef Matrix<double, Element::SpatialDimension, 1> Point;
  ifstream file(filename.c_str());
  string line;
  bool doingConnectivity = false;
  size_t nextNodeId = 0;
  if (file.is_open()) {
    while (file.good()) {
      getline (file,line);
      if (strcmp(line.c_str(), "") == 0) {
        break;
      }
      //printf("read line \"%s\"\n", line.c_str());
      if (strcmp(line.c_str(), "connectivity") == 0) {
        doingConnectivity = true;
        //printf("switching to reading connectivity\n");
      } else {
        vector<string> tokens;
        Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
        if (doingConnectivity == true) {
          if (tokens.size() != Element::NumberOfNodes) {
            printf("Problem in readMeshFromFile, connectivity tokens.size "
                   "= %zu\n", tokens.size());
            exit(1);
          }
          array<size_t, Element::NumberOfNodes> connection;
          for (size_t i = 0; i < Element::NumberOfNodes; ++i) {
            size_t index = size_t(atoi(tokens[i].c_str()));
            if (indexingStyle == OneIndexed && index == 0) {
              printf("Error: you said that file %s was one-indexed but "
                     "line %s has a connectivity with a zero index\n",
                     filename.c_str(), line.c_str());
              exit(1);
            }
            if (indexingStyle == OneIndexed) {
              --index;
            }
            connection[i] = index;
            if (connection[i] < 0 || connection[i] >= mesh->_nodes.size()) {
              printf("Bad vertex number (%zu) in element %zu of mesh file\n",
                     connection[i], mesh->_connectivity.size());
            }
          }
          mesh->_connectivity.push_back(connection);
        } else {
          if (tokens.size() != Element::SpatialDimension) {
            printf("Problem in readMeshFromFile, point tokens.size = %zu, "
                   "read from line \"%s\"\n", tokens.size(), line.c_str());
            exit(1);
          }
          Point point;
          for (size_t i = 0; i < Element::SpatialDimension; ++i) {
            point[i] = atof(tokens[i].c_str());
          }
          mesh->_nodes.push_back(NodeWithId<Point>(nextNodeId, point));
          nextNodeId++;
        }
      }
    }
    file.close();
    if (doingConnectivity == false){
      throwException("failed to capture connectivity from mesh file\n");
      exit(1);
    }
  } else {
    printf("Unable to open file at %s\n", filename.c_str());
    exit(1);
  }
}

template <class Element1, class Element2>
void
readTwoElementMeshFromFile(const string filename,
                           TwoElementMesh<Element1,Element2> * mesh,
                           IndexingStyle indexingStyle = ZeroIndexed) {
  typedef Matrix<double, Element1::SpatialDimension, 1> Point;
  ifstream file(filename.c_str());
  string line;
  bool doingConnectivity = false;
  size_t nextNodeId = 0;
  if (file.is_open()) {
    while (file.good()) {
      getline (file,line);
      if (strcmp(line.c_str(), "") == 0) {
        break;
      }
      //printf("read line \"%s\"\n", line.c_str());
      if (strcmp(line.c_str(), "connectivity") == 0) {
        doingConnectivity = true;
        //printf("switching to reading connectivity\n");
      } else {
        vector<string> tokens;
        Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
        if (doingConnectivity == true) {
          if (tokens.size() != Element1::NumberOfNodes &&
              tokens.size() != Element2::NumberOfNodes) {
            printf("Problem in readMeshFromFile, connectivity tokens.size "
                   "= %zu\n", tokens.size());
            exit(1);
          }
          if (tokens.size() == Element1::NumberOfNodes){
            // start of first type of element
            array<size_t, Element1::NumberOfNodes> connection;
            for (size_t i = 0; i < Element1::NumberOfNodes; ++i) {
              size_t index = size_t(atoi(tokens[i].c_str()));
              if (indexingStyle == OneIndexed && index == 0) {
                printf("Error: you said that file %s was one-indexed but "
                       "line %s has a connectivity with a zero index\n",
                       filename.c_str(), line.c_str());
                exit(1);
              }
              if (indexingStyle == OneIndexed) {
                --index;
              }
              connection[i] = index;
              if (connection[i] < 0 || connection[i] >= mesh->_nodes.size()) {
                printf("Bad vertex number (%zu) in element %zu of mesh file\n",
                       connection[i], mesh->_connectivity1.size());
              }
            }
            mesh->_connectivity1.push_back(connection);
            // end of first type of element
          }
          if (tokens.size() == Element2::NumberOfNodes){
            // start of second type of element
            array<size_t, Element2::NumberOfNodes> connection;
            for (size_t i = 0; i < Element2::NumberOfNodes; ++i) {
              size_t index = size_t(atoi(tokens[i].c_str()));
              if (indexingStyle == OneIndexed && index == 0) {
                printf("Error: you said that file %s was one-indexed but "
                       "line %s has a connectivity with a zero index\n",
                       filename.c_str(), line.c_str());
                exit(1);
              }
              if (indexingStyle == OneIndexed) {
                --index;
              }
              connection[i] = index;
              if (connection[i] < 0 || connection[i] >= mesh->_nodes.size()) {
                printf("Bad vertex number (%zu) in element %zu of mesh file\n",
                       connection[i], mesh->_connectivity2.size());
              }
            }
            mesh->_connectivity2.push_back(connection);
            // end of second type of element
          }

        } else {
          if (tokens.size() != Element1::SpatialDimension) {
            printf("Problem in readMeshFromFile, point tokens.size = %zu, "
                   "read from line \"%s\"\n", tokens.size(), line.c_str());
            exit(1);
          }
          Point point;
          for (size_t i = 0; i < Element1::SpatialDimension; ++i) {
            point[i] = atof(tokens[i].c_str());
          }
          mesh->_nodes.push_back(NodeWithId<Point>(nextNodeId, point));
          nextNodeId++;
        }
      }
    }
    file.close();
    if (doingConnectivity == false){
      throwException("failed to capture connectivity from mesh file\n");
      exit(1);
    }
  } else {
    printf("Unable to open file at %s\n", filename.c_str());
    exit(1);
  }
}
template <class Element1, class Element2 >
void
buildBandGapAuxeticLattice3D(const double sparLength,
                             const double zHeight,
                             const double thetaDegrees,
                             TwoElementMesh<Element1,Element2> * mesh,
                             SingleElementMesh<Element1> * representativeMesh){
  typedef Vector3d Point;
  typedef array<size_t,2> Connection1;
  typedef array<size_t,1> Connection2;
  array<Point,6> points;
  const double thetaRadians = thetaDegrees / 180. * M_PI;
  const double xLength = sparLength * sin(thetaRadians);
  const double zLength = sparLength * cos(thetaRadians);

  points[0] <<  0       , 0        , 0;
  points[1] <<  xLength, 0        , zLength;
  points[2] << -xLength, 0        , zLength;
  points[3] <<  0       , xLength , zLength;
  points[4] <<  0       ,-xLength , zLength;
  points[5] <<  0       , 0       , zHeight;

  for (size_t nodeIndex = 0; nodeIndex < points.size(); nodeIndex++){
    mesh->_nodes.push_back(NodeWithId<Point>(nodeIndex,points[nodeIndex]));
    representativeMesh->_nodes.push_back(NodeWithId<Point>(nodeIndex,points[nodeIndex]));
  }

  // beams
  array<Connection1,5> connections;
  connections[0] = {{0,1}};
  connections[1] = {{0,2}};
  connections[2] = {{0,3}};
  connections[3] = {{0,4}};
  connections[4] = {{0,5}};
  for (size_t elementIndex = 0; elementIndex < connections.size(); elementIndex++){
    mesh->_connectivity1.push_back(connections[elementIndex]);
    representativeMesh->_connectivity.push_back(connections[elementIndex]);
  }
  // point masses
  array<Connection2,4> massConnections;
  massConnections[0] = {{1}};
  massConnections[1] = {{2}};
  massConnections[2] = {{3}};
  massConnections[3] = {{4}};
  for (size_t elementIndex = 0; elementIndex < massConnections.size(); elementIndex++){
    mesh->_connectivity2.push_back(massConnections[elementIndex]);
  }

}

template <class Element>
void
readMeshFromFemapFile(const string filename,
                      SingleElementMesh<Element> * mesh,
                      vector<size_t> * groupIdContainer,
                      vector<vector<size_t> > * elementContainer,
                      IndexingStyle indexingStyle = ZeroIndexed) {
  vector<vector<size_t> > nodeContainer;
  readMeshFromFemapFile<Element>(filename,
                                 mesh,
                                 groupIdContainer,
                                 elementContainer,
                                 &nodeContainer,
                                 indexingStyle);
}

template <class Element>
void
readMeshFromFemapFile(const string filename,
                      SingleElementMesh<Element> * mesh,
                      vector<size_t> * groupIdContainer,
                      vector<vector<size_t> > * elementContainer,
                      vector<vector<size_t> > * nodeContainer,
                      IndexingStyle indexingStyle = ZeroIndexed) {
  typedef Matrix<double, Element::SpatialDimension, 1> Point;
  //ignoring variable indexingStyle as it is never used in the function
  ignoreUnusedVariables(indexingStyle);
  
  if (Element::SpatialDimension == 3 && Element::NumberOfNodes >= 3 
                                     && Element::NumberOfNodes <= 4) {

    const std::string femapBlockEndMark = "   -1";
    std::map<size_t,size_t> nodeIdMap;
    std::map<size_t,size_t> elementIdMap;
    ifstream file(filename.c_str());
    string line, line2;
    int labelId = 0;
    if (file.is_open()) {
      while (file.good()) {
        getline(file,line);
        file >> labelId; // is this safe?
        getline(file,line);

        switch ( labelId ) {
          // read nodes
        case 403: {
          Point point;
          getline(file,line);
          line2.assign(line, 0, 5);
          size_t nextNodeId = 0;

          while (file.good() && line2 != femapBlockEndMark) {
            vector<string> tokens;
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t nodeId = size_t(atoi(tokens[0].c_str()));
            nodeIdMap.insert(std::make_pair(nodeId,nextNodeId));
            point[0] = atof(tokens[11].c_str());
            point[1] = atof(tokens[12].c_str());
            point[2] = atof(tokens[13].c_str());
            mesh->_nodes.push_back(typename Element::Node(nextNodeId, point));
            nextNodeId++;
            getline(file,line);
            line2.assign(line, 0, 5);
          }
        } break;

          // read elements
        case 404: {
          Point point;
          getline(file,line);
          line2.assign(line, 0, 5);
          size_t nextElementId = 0;

          while (file.good() && line2 != femapBlockEndMark) {
            vector<string> tokens;
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t elementId = size_t(atoi(tokens[0].c_str()));
            elementIdMap.insert(std::make_pair(elementId,nextElementId));
            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            array<size_t, Element::NumberOfNodes> localConnection;
            for (size_t i = 0; i < Element::NumberOfNodes; ++i) {
              localConnection[i] = size_t(atoi(tokens[i].c_str()));
              if (localConnection[i] == 0) {
                printf("Something went wrong reading element connectivities in line \"%s\"\n", line.c_str());
                exit(1);
              }
            }
            getline(file,line);
            getline(file,line);
            getline(file,line);
            getline(file,line);
            getline(file,line);

            array<size_t, Element::NumberOfNodes> connection;
            for (size_t i = 0; i < Element::NumberOfNodes; ++i) {
              connection[i] = nodeIdMap.find(localConnection[i])->second;
            }
            mesh->_connectivity.push_back(connection);

            nextElementId++;
            getline(file,line);
            line2.assign(line, 0, 5);
          }
        } break;

        case 408: {
          // read groups
          getline(file,line);
          vector<string> tokens;
          Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
          size_t groupId = size_t(atoi(tokens[0].c_str())); //starts counting at 1

          vector<size_t> groupElementList;
          vector<size_t> groupNodeList;

          while (groupId != -1) {
            groupIdContainer->push_back(groupId);
            getline(file,line);
            cout << "Group Id: " << groupId << " ---  Title: " << line << endl;

            // skip clipping info
            for (size_t i = 0; i < 21; i++) {
              getline(file,line);
            }

            // read group rules, don't care to store
            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t maxRule = size_t(atoi(tokens[0].c_str()));

            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t groupRuleType = size_t(atoi(tokens[0].c_str()));

            while (groupRuleType != -1) {
              if (groupRuleType < maxRule) {
                getline(file,line);
                tokens.resize(0);
                Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                size_t startId = size_t(atoi(tokens[0].c_str()));
                while (startId != -1) {
                  getline(file,line);
                  tokens.resize(0);
                  Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                  startId = size_t(atoi(tokens[0].c_str()));
                }
              } else {
                getline(file,line);
              }
              getline(file,line);
              tokens.resize(0);
              Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
              groupRuleType = size_t(atoi(tokens[0].c_str()));
            }

            // read group lists
            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t maxList = size_t(atoi(tokens[0].c_str()));

            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t groupListType = size_t(atoi(tokens[0].c_str()));

            while (groupListType != -1) {
              if (groupListType < maxList) {
                getline(file,line);
                tokens.resize(0);
                Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                size_t entityId = size_t(atoi(tokens[0].c_str()));
                // store element list
                while (entityId != -1 && groupListType == 8) {
                  groupElementList.push_back(elementIdMap.find(entityId)->second);
                  getline(file,line);
                  tokens.resize(0);
                  Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                  entityId = size_t(atoi(tokens[0].c_str()));
                }
                // store node list
                while (entityId != -1 && groupListType == 7) {
                  groupNodeList.push_back(nodeIdMap.find(entityId)->second);
                  getline(file,line);
                  tokens.resize(0);
                  Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                  entityId = size_t(atoi(tokens[0].c_str()));
                }
                while (entityId != -1 && groupListType != 8 && groupListType != 7) {
                  getline(file,line);
                  tokens.resize(0);
                  Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                  entityId = size_t(atoi(tokens[0].c_str()));
                }
              } else {
                getline(file,line);
              }
              getline(file,line);
              tokens.resize(0);
              Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
              groupListType = size_t(atoi(tokens[0].c_str()));
            }

            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            groupId = size_t(atoi(tokens[0].c_str()));

            elementContainer->push_back(groupElementList);
            groupElementList.resize(0);

            nodeContainer->push_back(groupNodeList);
            groupNodeList.resize(0);
          }
        } break;

        default:
          while (file.good() && line != femapBlockEndMark) {
            getline(file,line);
            line.assign(line, 0, 5);
          }
          break;
        }
      }
    }
    else {
      printf("Unable to open file at %s\n", filename.c_str());
      exit(1);
    }
  }
  else if (Element::SpatialDimension == 2 && Element::NumberOfNodes == 3) {

    const std::string femapBlockEndMark = "   -1";
    std::map<size_t,size_t> nodeIdMap;
    std::map<size_t,size_t> elementIdMap;
    ifstream file(filename.c_str());
    string line, line2;
    int labelId = 0;
    if (file.is_open()) {
      while (file.good()) {
        getline(file,line);
        file >> labelId; // is this safe?
        getline(file,line);

        switch ( labelId ) {
          // read nodes
        case 403: {
          Point point;
          getline(file,line);
          line2.assign(line, 0, 5);
          size_t nextNodeId = 0;

          while (file.good() && line2 != femapBlockEndMark) {
            vector<string> tokens;
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t nodeId = size_t(atoi(tokens[0].c_str()));
            nodeIdMap.insert(std::make_pair(nodeId,nextNodeId));
            point[0] = atof(tokens[11].c_str());
            point[1] = atof(tokens[12].c_str());
            mesh->_nodes.push_back(typename Element::Node(nextNodeId, point));
            nextNodeId++;
            getline(file,line);
            line2.assign(line, 0, 5);
          }
        } break;

          // read elements
        case 404: {
          Point point;
          getline(file,line);
          line2.assign(line, 0, 5);
          size_t nextElementId = 0;

          while (file.good() && line2 != femapBlockEndMark) {
            vector<string> tokens;
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t elementId = size_t(atoi(tokens[0].c_str()));
            elementIdMap.insert(std::make_pair(elementId,nextElementId));
            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            array<size_t, Element::NumberOfNodes> localConnection;
            for (size_t i = 0; i < Element::NumberOfNodes; ++i) {
              localConnection[i] = size_t(atoi(tokens[i].c_str()));
              if (localConnection[i] == 0) {
                printf("Something went wrong reading element connectivities in line \"%s\"\n", line.c_str());
                exit(1);
              }
            }
            getline(file,line);
            getline(file,line);
            getline(file,line);
            getline(file,line);
            getline(file,line);

            array<size_t, Element::NumberOfNodes> connection;
            for (size_t i = 0; i < Element::NumberOfNodes; ++i) {
              connection[i] = nodeIdMap.find(localConnection[i])->second;
            }
            mesh->_connectivity.push_back(connection);

            nextElementId++;
            getline(file,line);
            line2.assign(line, 0, 5);
          }
        } break;

        case 408: {
          // read groups
          getline(file,line);
          vector<string> tokens;
          Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
          size_t groupId = size_t(atoi(tokens[0].c_str())); //starts counting at 1

          vector<size_t> groupElementList;
          vector<size_t> groupNodeList;

          while (groupId != -1) {
            groupIdContainer->push_back(groupId);
            getline(file,line);
            cout << "Group Id: " << groupId << " ---  Title: " << line << endl;

            // skip clipping info
            for (size_t i = 0; i < 21; i++) {
              getline(file,line);
            }

            // read group rules, don't care to store
            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t maxRule = size_t(atoi(tokens[0].c_str()));

            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t groupRuleType = size_t(atoi(tokens[0].c_str()));

            while (groupRuleType != -1) {
              if (groupRuleType < maxRule) {
                getline(file,line);
                tokens.resize(0);
                Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                size_t startId = size_t(atoi(tokens[0].c_str()));
                while (startId != -1) {
                  getline(file,line);
                  tokens.resize(0);
                  Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                  startId = size_t(atoi(tokens[0].c_str()));
                }
              } else {
                getline(file,line);
              }
              getline(file,line);
              tokens.resize(0);
              Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
              groupRuleType = size_t(atoi(tokens[0].c_str()));
            }

            // read group lists
            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t maxList = size_t(atoi(tokens[0].c_str()));

            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            size_t groupListType = size_t(atoi(tokens[0].c_str()));

            while (groupListType != -1) {
              if (groupListType < maxList) {
                getline(file,line);
                tokens.resize(0);
                Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                size_t entityId = size_t(atoi(tokens[0].c_str()));
                // store element list
                while (entityId != -1 && groupListType == 8) {
                  groupElementList.push_back(elementIdMap.find(entityId)->second);
                  getline(file,line);
                  tokens.resize(0);
                  Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                  entityId = size_t(atoi(tokens[0].c_str()));
                }
                // store node list
                while (entityId != -1 && groupListType == 7) {
                  groupNodeList.push_back(nodeIdMap.find(entityId)->second);
                  getline(file,line);
                  tokens.resize(0);
                  Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                  entityId = size_t(atoi(tokens[0].c_str()));
                }
                while (entityId != -1 && groupListType != 8 && groupListType != 7) {
                  getline(file,line);
                  tokens.resize(0);
                  Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
                  entityId = size_t(atoi(tokens[0].c_str()));
                }
              } else {
                getline(file,line);
              }
              getline(file,line);
              tokens.resize(0);
              Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
              groupListType = size_t(atoi(tokens[0].c_str()));
            }

            getline(file,line);
            tokens.resize(0);
            Utilities::tokenize(line, ",", Utilities::Trim, &tokens);
            groupId = size_t(atoi(tokens[0].c_str()));

            elementContainer->push_back(groupElementList);
            groupElementList.resize(0);

            nodeContainer->push_back(groupNodeList);
            groupNodeList.resize(0);
          }
        } break;

        default:
          while (file.good() && line != femapBlockEndMark) {
            getline(file,line);
            line.assign(line, 0, 5);
          }
          break;
        }
      }
    }
    else {
      printf("Unable to open file at %s\n", filename.c_str());
      exit(1);
    }
  }
  else {
    printf("error: readMeshFromFemapFile can only handle 3D data of triangle and tetrahedral meshes\n");
  }
}

template <class Element>
void
buildType1AuxeticLattice2D(const double lowerEdgeLength,
                           const double height,
                           const double thetaDegrees,
                           SingleElementMesh<Element> * mesh) {

  typedef Vector3d Point;
  //typedef NodeWithId<Point> Node;
  typedef array<size_t,2> Connection;



  const double thetaRadians = thetaDegrees/180.*M_PI;

  array<Point,9> points;
  points[0] << 0                                                      ,0        ,0;
  points[1] << lowerEdgeLength                                        ,0        ,0;
  points[2] << lowerEdgeLength - (height/2. / tan(thetaRadians))      ,height/2.,0;
  points[3] << lowerEdgeLength                                        ,height   ,0;
  points[4] << 0                                                      ,height   ,0;
  points[5] << height/2. / tan(thetaRadians)                          ,height/2.,0;
  points[6] << 2 * (lowerEdgeLength - (height/2. / tan(thetaRadians))),0        ,0;
  points[7] << 2 * lowerEdgeLength - (height/2. / tan(thetaRadians))  ,height/2.,0;
  points[8] << 2 * (lowerEdgeLength - (height/2. / tan(thetaRadians))),height   ,0;

  if (points[5](0) > points[2](0)){
    errorStatement("The unit cell has overlapping elements based on the theta, height,"
                   "and lowerEdgeLength you have prescribed\n");
    exit(1);
  }

  for (size_t nodeIndex = 0; nodeIndex < points.size(); nodeIndex++){
    mesh->_nodes.push_back(NodeWithId<Point>(nodeIndex,points[nodeIndex]));
  }

  array<Connection,9> connections;
  connections[0] = {{0,1}};
  connections[1] = {{1,2}};
  connections[2] = {{2,3}};
  connections[3] = {{3,4}};
  connections[4] = {{4,5}};
  connections[5] = {{5,0}};
  connections[6] = {{2,7}};
  connections[7] = {{6,7}};
  connections[8] = {{7,8}};
  for (size_t elementIndex = 0; elementIndex < connections.size(); elementIndex++){
    mesh->_connectivity.push_back(connections[elementIndex]);
  }
}

template <class Element>
void
trimMeshWithBox(const array<Matrix<double,Element::SpatialDimension,1>,2> boundingBox,
                SingleElementMesh<Element> * mesh){
  size_t numberOfElements = mesh->_connectivity.size();
  size_t numberOfNodes = mesh->_nodes.size();

  // finding the bounding box values
  array<double,Element::SpatialDimension> maxValues;
  array<double,Element::SpatialDimension> minValues;
  for (size_t spatialIndex = 0; spatialIndex < Element::SpatialDimension; spatialIndex++){
    if (boundingBox[0](spatialIndex) > boundingBox[1](spatialIndex)){
      maxValues[spatialIndex] = boundingBox[0](spatialIndex);
      minValues[spatialIndex] = boundingBox[1](spatialIndex);
    }
    else{
      maxValues[spatialIndex] = boundingBox[1](spatialIndex);
      minValues[spatialIndex] = boundingBox[0](spatialIndex);
    }
  }

  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    const Matrix<double,Element::SpatialDimension,1> nodePosition = mesh->_nodes[nodeIndex]._position;
    size_t nodeFlag = 0;
    // trying to find a node outside the bounding box
    if (nodePosition(0) > maxValues[0] ||
        nodePosition(0) < minValues[0]){
      nodeFlag = 1;
    }
    if (Element::SpatialDimension >1){
      if (nodePosition(1) > maxValues[1] ||
          nodePosition(1) < minValues[1]){
        nodeFlag = 1;
      }
    }
    if (Element::SpatialDimension == 3){
      if (nodePosition(2) > maxValues[2] ||
          nodePosition(2) < minValues[2]){
        nodeFlag = 1;
      }
    }

    // if the node is outside the bounding box remove it and the element it belonged to
    if (nodeFlag == 1){
      mesh->_nodes.erase( mesh->_nodes.begin() + nodeIndex);
      numberOfNodes--;
      // decrement all the ids of nodes that followed the deleted node
      for (size_t index = nodeIndex; index < numberOfNodes; index++){
        mesh->_nodes[index]._id--;
      }

      // PRIMARY ELEMENTS
      // now go through the connectivity and see if any of the elements have that node
      for (size_t elementIndex = 0; elementIndex < numberOfElements; elementIndex++){
        for (size_t nodeIndex2 = 0; nodeIndex2 < Element::NumberOfNodes; nodeIndex2++){
          if (mesh->_connectivity[elementIndex][nodeIndex2] == nodeIndex){
            // erase the element having the deleted node
            mesh->_connectivity.erase( mesh->_connectivity.begin() + elementIndex);
            numberOfElements--;
            elementIndex--;
          }
        }
      }
      // reduce the IDS of nodes in the connectivity matrix
      for (size_t elementIndex2 = 0; elementIndex2 < numberOfElements; elementIndex2++){
        for (size_t elementNodeIndex2 = 0;
             elementNodeIndex2 < Element::NumberOfNodes;
             elementNodeIndex2++){
          if (mesh->_connectivity[elementIndex2][elementNodeIndex2] > nodeIndex){
            mesh->_connectivity[elementIndex2][elementNodeIndex2]--;
          }
        }
      }
      nodeIndex--;
    }
  }
}

template <class PrimaryElement, class SecondaryElement>
void
trimMeshWithBox(const array<Matrix<double,PrimaryElement::SpatialDimension,1>,2> boundingBox,
                TwoElementMesh<PrimaryElement,SecondaryElement> * mesh){
  size_t numberOfPrimaryElements = mesh->_connectivity1.size();
  size_t numberOfSecondaryElements = mesh->_connectivity2.size();
  size_t numberOfNodes = mesh->_nodes.size();

  // finding the bounding box values
  array<double,PrimaryElement::SpatialDimension> maxValues;
  array<double,PrimaryElement::SpatialDimension> minValues;
  for (size_t spatialIndex = 0; spatialIndex < PrimaryElement::SpatialDimension; spatialIndex++){
    if (boundingBox[0](spatialIndex) > boundingBox[1](spatialIndex)){
      maxValues[spatialIndex] = boundingBox[0](spatialIndex);
      minValues[spatialIndex] = boundingBox[1](spatialIndex);
    }
    else{
      maxValues[spatialIndex] = boundingBox[1](spatialIndex);
      minValues[spatialIndex] = boundingBox[0](spatialIndex);
    }
  }
  printf("maxValues: %e,%e,%e\n",maxValues[0],maxValues[1],maxValues[2]);
  printf("minValues: %e,%e,%e\n",minValues[0],minValues[1],minValues[2]);

  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    //cout << "nodeIndex" << nodeIndex << endl;
    const Matrix<double,PrimaryElement::SpatialDimension,1> nodePosition = mesh->_nodes[nodeIndex]._position;
    size_t nodeFlag = 0;
    // trying to find a node outside the bounding box
    if (nodePosition(0) > maxValues[0] ||
        nodePosition(0) < minValues[0]){
      nodeFlag = 1;
    }
    if (PrimaryElement::SpatialDimension >1){
      if (nodePosition(1) > maxValues[1] ||
          nodePosition(1) < minValues[1]){
        nodeFlag = 1;
      }
    }
    if (PrimaryElement::SpatialDimension == 3){
      if (nodePosition(2) > maxValues[2] ||
          nodePosition(2) < minValues[2]){
        nodeFlag = 1;
      }
    }


    // if the node is outside the bounding box remove it and the element it belonged to
    if (nodeFlag == 1){
      mesh->_nodes.erase( mesh->_nodes.begin() + nodeIndex);
      numberOfNodes--;
      // decrement all the ids of nodes that followed the deleted node
      for (size_t index = nodeIndex; index < numberOfNodes; index++){
        mesh->_nodes[index]._id--;
      }

      // PRIMARY ELEMENTS
      // now go through the connectivity and see if any of the elements have that node
      for (size_t elementIndex = 0; elementIndex < numberOfPrimaryElements; elementIndex++){
        for (size_t nodeIndex2 = 0; nodeIndex2 < PrimaryElement::NumberOfNodes; nodeIndex2++){
          if (mesh->_connectivity1[elementIndex][nodeIndex2] == nodeIndex){
            // erase the element having the deleted node
            mesh->_connectivity1.erase( mesh->_connectivity1.begin() + elementIndex);
            numberOfPrimaryElements--;
            elementIndex--;
          }
        }
      }
      // reduce the IDS of nodes in the connectivity matrix
      for (size_t elementIndex2 = 0; elementIndex2 < numberOfPrimaryElements; elementIndex2++){
        for (size_t elementNodeIndex2 = 0;
             elementNodeIndex2 < PrimaryElement::NumberOfNodes;
             elementNodeIndex2++){
          if (mesh->_connectivity1[elementIndex2][elementNodeIndex2] > nodeIndex){
            mesh->_connectivity1[elementIndex2][elementNodeIndex2]--;
          }
        }
      }

      // SECONDARY ELEMENTS
      // now go through the connectivity and see if any of the elements have that node
      for (size_t elementIndex = 0; elementIndex < numberOfSecondaryElements; elementIndex++){
        for (size_t nodeIndex2 = 0; nodeIndex2 < SecondaryElement::NumberOfNodes; nodeIndex2++){
          if (mesh->_connectivity2[elementIndex][nodeIndex2] == nodeIndex){
            // erase the element having the deleted node
            mesh->_connectivity2.erase( mesh->_connectivity2.begin() + elementIndex);
            numberOfSecondaryElements--;
            elementIndex--;
          }
        }
      }
      // reduce the IDS of nodes in the connectivity matrix
      for (size_t elementIndex2 = 0; elementIndex2 < numberOfSecondaryElements; elementIndex2++){
        for (size_t elementNodeIndex2 = 0;
             elementNodeIndex2 < SecondaryElement::NumberOfNodes;
             elementNodeIndex2++){
          if (mesh->_connectivity2[elementIndex2][elementNodeIndex2] > nodeIndex){
            mesh->_connectivity2[elementIndex2][elementNodeIndex2]--;
          }
        }
      }
      nodeIndex--;
    }
  }
}

template <class Element>
void
translateAndCopyMesh(const Matrix<double,Element::SpatialDimension,1> translation,
                     SingleElementMesh<Element> originalMesh,
                     SingleElementMesh<Element> * mesh){

  typedef typename Element::Point Point;

  const size_t numberOfNodes = originalMesh._nodes.size();
  const size_t numberOfElements = originalMesh._connectivity.size();

  const size_t numberOfNodesInMesh = mesh->_nodes.size();
  size_t freeNodeId = mesh->_nodes.size();
  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    const Point newPosition =
      originalMesh._nodes[nodeIndex]._position + translation;
    mesh->_nodes.push_back(typename Element::Node(freeNodeId,newPosition));
    freeNodeId++;
  }

  for (size_t elementIndex = 0; elementIndex < numberOfElements; elementIndex++){
    array<size_t,Element::NumberOfNodes>
      connection = originalMesh._connectivity[elementIndex];
    for (size_t nodeIndex = 0; nodeIndex < Element::NumberOfNodes; nodeIndex++){
      connection[nodeIndex] += numberOfNodesInMesh;
    }
    mesh->_connectivity.push_back(connection);
  }
}

template <class PrimaryElement, class SecondaryElement>
void
translateAndCopyMesh(const Matrix<double,PrimaryElement::SpatialDimension,1> translation,
                     TwoElementMesh<PrimaryElement,SecondaryElement> originalMesh,
                     TwoElementMesh<PrimaryElement,SecondaryElement> * mesh){

  typedef typename PrimaryElement::Point Point;

  const size_t numberOfNodes = originalMesh._nodes.size();
  const size_t numberOfPrimaryElements = originalMesh._connectivity1.size();
  const size_t numberOfSecondaryElements = originalMesh._connectivity2.size();

  const size_t numberOfNodesInMesh = mesh->_nodes.size();
  size_t freeNodeId = mesh->_nodes.size();
  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    const Point newPosition =
      originalMesh._nodes[nodeIndex]._position + translation;
    mesh->_nodes.push_back(NodeWithId<Point>(freeNodeId,newPosition));
    freeNodeId++;
  }

  for (size_t elementIndex = 0; elementIndex < numberOfPrimaryElements;
       elementIndex++){
    array<size_t,PrimaryElement::NumberOfNodes>
      connection = originalMesh._connectivity1[elementIndex];
    for (size_t nodeIndex = 0; nodeIndex < PrimaryElement::NumberOfNodes; nodeIndex++){
      connection[nodeIndex] += numberOfNodesInMesh;
    }
    mesh->_connectivity1.push_back(connection);
  }

  for (size_t elementIndex = 0; elementIndex < numberOfSecondaryElements; elementIndex++){
    array<size_t,SecondaryElement::NumberOfNodes>
      connection = originalMesh._connectivity2[elementIndex];
    for (size_t nodeIndex = 0; nodeIndex < SecondaryElement::NumberOfNodes; nodeIndex++){
      connection[nodeIndex] += numberOfNodesInMesh;
    }
    mesh->_connectivity2.push_back(connection);
  }
}

template <class Element>
void
removeOverlapingNodesAndElements(SingleElementMesh<Element> * mesh,
                                 const double spatialTolerance = 10){
  size_t numberOfNodes = mesh->_nodes.size();
  size_t numberOfElements = mesh->_connectivity.size();
  // find an acceptable tolerance to determine if one node is "on top of" another
  double maxValue = (mesh->_nodes[0]._position).norm();
  double minValue = (mesh->_nodes[0]._position).norm();
  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    if ((mesh->_nodes[nodeIndex]._position).norm() > maxValue){
      maxValue = (mesh->_nodes[nodeIndex]._position).norm();
    }
    if ((mesh->_nodes[nodeIndex]._position).norm() < minValue){
      minValue = (mesh->_nodes[nodeIndex]._position).norm();
    }
  }
  const double tolerance = spatialTolerance;

  // looking for nodes with the same position
  for (size_t nodeIndex1 = 0; nodeIndex1 < numberOfNodes; nodeIndex1++){
    const MatrixXd nodePosition1 = mesh->_nodes[nodeIndex1]._position;
    for (size_t nodeIndex2 = 0; nodeIndex2 < numberOfNodes; nodeIndex2++){
      const MatrixXd nodePosition2 = mesh->_nodes[nodeIndex2]._position;
      if (nodeIndex1 != nodeIndex2 &&
          (nodePosition1 - nodePosition2).norm() < tolerance){
        // erase the double node
        mesh->_nodes.erase( mesh->_nodes.begin() + nodeIndex2);
        numberOfNodes--;
        // change all the connectivities of double node to kept node // OK
        for (size_t elementIndex = 0; elementIndex < numberOfElements; elementIndex++){
          for (size_t elementNodeIndex = 0; elementNodeIndex < Element::NumberOfNodes; elementNodeIndex++){
            if (mesh->_connectivity[elementIndex][elementNodeIndex] == nodeIndex2){
              //printf("Replaced NodeIndex %u with NodeIndex %u in element %u  in the connectivity\n",nodeIndex2,nodeIndex1,elementIndex);
              mesh->_connectivity[elementIndex][elementNodeIndex] = nodeIndex1;
            }
          }
        }
        // reduce the IDs of all later nodes in the nodes themselves // OK
        for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
          if (mesh->_nodes[nodeIndex]._id > nodeIndex2){
            mesh->_nodes[nodeIndex]._id--;
          }
        }
        // reduce the IDS of all later nodes in the connectivity matrix
        for (size_t elementIndex = 0; elementIndex < numberOfElements; elementIndex++){
          for (size_t elementNodeIndex = 0; elementNodeIndex < Element::NumberOfNodes; elementNodeIndex++){
            if (mesh->_connectivity[elementIndex][elementNodeIndex] > nodeIndex2){
              mesh->_connectivity[elementIndex][elementNodeIndex]--;
            }
          }
        }
        if (nodeIndex2>0){
          nodeIndex2--;
        }
      }
    }
  }

  //printf("Mesh Connectivity After Node Editing\n");
  /*
    for (size_t elementIndex = 0; elementIndex < numberOfElements; elementIndex++){
    const array<size_t,2> connection = mesh->_connectivity[elementIndex];
    //printf("%u,%u\n",connection[0],connection[1]);
    }
  */

  vector<size_t> elementsMarkedForDeletion;
  for (size_t elementIndex = 0; elementIndex < numberOfElements; elementIndex++){
    elementsMarkedForDeletion.push_back(0);
  }

  vector<array<size_t,Element::NumberOfNodes>> possibleCombinations;
  // now let's look for elements with the same node numbers
  for (size_t elementIndex = 0; elementIndex < numberOfElements; elementIndex++){
    for (size_t comboIndex = 0; comboIndex < pow(Element::NumberOfNodes,Element::NumberOfNodes); comboIndex++){
      size_t number = comboIndex;
      vector<size_t> combo;
      while (number >= Element::NumberOfNodes){
        size_t remainder = number%(Element::NumberOfNodes);
        number = floor(number/(Element::NumberOfNodes));
        combo.push_back(remainder);
      }
      combo.push_back(number);
      while (combo.size() < Element::NumberOfNodes){
        combo.push_back(0);
      }
      // now let's make a fake connection out of the combo
      array<size_t,Element::NumberOfNodes> tempConnection;
      for (size_t comboSubIndex = 0; comboSubIndex < Element::NumberOfNodes; comboSubIndex++){
        if (comboSubIndex < combo.size()){
          tempConnection[comboSubIndex] = combo[comboSubIndex];
        }
        else{
          tempConnection[comboSubIndex] = 0;
        }
      }
      possibleCombinations.push_back(tempConnection);
    }
  }
  // now that we have all the possible node combinations for the same element
  // in space, let's find if any other elements in the matrix match this element
  for (size_t elementIndex1 = 0; elementIndex1 < numberOfElements; elementIndex1++){
    for (size_t elementIndex2 = elementIndex1; elementIndex2 < numberOfElements; elementIndex2++){
      for (size_t comboIndex = 0; comboIndex < pow(Element::NumberOfNodes,Element::NumberOfNodes); comboIndex++){
        size_t elementFlag = 0;
        // turn possibleCombination into node number of elementIndex2
        array<size_t,Element::NumberOfNodes> possibleConnection;
        for (size_t nodeIndex = 0; nodeIndex < Element::NumberOfNodes; nodeIndex++){
          const size_t nodeToUse = possibleCombinations[comboIndex][nodeIndex];
          possibleConnection[nodeIndex] = mesh->_connectivity[elementIndex1][nodeToUse];
        }
        // check if the possibleConnection is the same as the 1st element
        for (size_t nodeIndex = 0; nodeIndex < Element::NumberOfNodes; nodeIndex++){
          if (possibleConnection[nodeIndex] == mesh->_connectivity[elementIndex2][nodeIndex]){
            elementFlag++;
          }
        }
        if (elementFlag == Element::NumberOfNodes &&
            elementIndex1 != elementIndex2){
          elementsMarkedForDeletion[elementIndex2] = 1;
        }
      }
      for (size_t checkIndex = 0; checkIndex < numberOfElements; checkIndex++){
      }
    }
    for (size_t elementIndex = 0; elementIndex < numberOfElements; elementIndex++){
      if (elementsMarkedForDeletion[elementIndex] == 1){
        mesh->_connectivity.erase(mesh->_connectivity.begin() + elementIndex);
        elementsMarkedForDeletion.erase(elementsMarkedForDeletion.begin() + elementIndex);
        numberOfElements--;
        if (elementIndex > 0){
          elementIndex--;
        }
      }
    }
  }
}

template <class PrimaryElement, class SecondaryElement>
void
removeOverlapingNodesAndElements(TwoElementMesh<PrimaryElement,SecondaryElement> * mesh){
  size_t numberOfNodes = mesh->_nodes.size();
  size_t numberOfPrimaryElements = mesh->_connectivity1.size();
  size_t numberOfSecondaryElements = mesh->_connectivity2.size();
  // find an acceptable tolerance to determine if one node is "on top of" another
  double maxValue = (mesh->_nodes[0]._position).norm();
  double minValue = (mesh->_nodes[0]._position).norm();
  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    if ((mesh->_nodes[nodeIndex]._position).norm() > maxValue){
      maxValue = (mesh->_nodes[nodeIndex]._position).norm();
    }
    if ((mesh->_nodes[nodeIndex]._position).norm() < minValue){
      minValue = (mesh->_nodes[nodeIndex]._position).norm();
    }
  }
  const double tolerance = (maxValue - minValue)/1e3;

  // looking for nodes with the same position
  for (size_t nodeIndex1 = 0; nodeIndex1 < numberOfNodes; nodeIndex1++){
    const MatrixXd nodePosition1 = mesh->_nodes[nodeIndex1]._position;
    for (size_t nodeIndex2 = 0; nodeIndex2 < numberOfNodes; nodeIndex2++){
      const MatrixXd nodePosition2 = mesh->_nodes[nodeIndex2]._position;
      if (nodeIndex1 != nodeIndex2 &&
          (nodePosition1 - nodePosition2).norm() < tolerance){
        // erase the double node
        mesh->_nodes.erase( mesh->_nodes.begin() + nodeIndex2);
        numberOfNodes--;
        // change all the connectivities of double node to kept node
        for (size_t elementIndex = 0; elementIndex < numberOfPrimaryElements;
             elementIndex++){
          for (size_t elementNodeIndex = 0;
               elementNodeIndex < PrimaryElement::NumberOfNodes; elementNodeIndex++){
            if (mesh->_connectivity1[elementIndex][elementNodeIndex] == nodeIndex2){
              mesh->_connectivity1[elementIndex][elementNodeIndex] = nodeIndex1;
            }
          }
        }
        for (size_t elementIndex = 0; elementIndex < numberOfSecondaryElements;
             elementIndex++){
          for (size_t elementNodeIndex = 0;
               elementNodeIndex < SecondaryElement::NumberOfNodes; elementNodeIndex++){
            if (mesh->_connectivity2[elementIndex][elementNodeIndex] == nodeIndex2){
              mesh->_connectivity2[elementIndex][elementNodeIndex] = nodeIndex1;
            }
          }
        }
        // reduce the IDs of all later nodes in the nodes themselves // OK
        for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
          if (mesh->_nodes[nodeIndex]._id > nodeIndex2){
            mesh->_nodes[nodeIndex]._id--;
          }
        }
        // reduce the IDS of all later nodes in the connectivity matrix
        for (size_t elementIndex = 0; elementIndex < numberOfPrimaryElements;
             elementIndex++){
          for (size_t elementNodeIndex = 0;
               elementNodeIndex < PrimaryElement::NumberOfNodes; elementNodeIndex++){
            if (mesh->_connectivity1[elementIndex][elementNodeIndex] > nodeIndex2){
              mesh->_connectivity1[elementIndex][elementNodeIndex]--;
            }
          }
        }
        // reduce the IDS of all later nodes in the connectivity matrix
        for (size_t elementIndex = 0; elementIndex < numberOfSecondaryElements;
             elementIndex++){
          for (size_t elementNodeIndex = 0;
               elementNodeIndex < SecondaryElement::NumberOfNodes; elementNodeIndex++){
            if (mesh->_connectivity2[elementIndex][elementNodeIndex] > nodeIndex2){
              mesh->_connectivity2[elementIndex][elementNodeIndex]--;
            }
          }
        }
        if (nodeIndex1>0){
        }
        if (nodeIndex2>0){
          nodeIndex2--;
        }
      }
    }
  }
  // PRIMARY ELEMENTS //
  vector<size_t> primaryElementsMarkedForDeletion;
  // just fill the vector to keep track of who to delete later
  for (size_t elementIndex = 0; elementIndex < numberOfPrimaryElements; elementIndex++){
    primaryElementsMarkedForDeletion.push_back(0);
  }

  vector<array<size_t,PrimaryElement::NumberOfNodes>> possibleCombinationsPrimary;
  // now let's look for elements with the same node numbers
  for (size_t elementIndex = 0; elementIndex < numberOfPrimaryElements; elementIndex++){
    for (size_t comboIndex = 0; comboIndex < pow(PrimaryElement::NumberOfNodes,PrimaryElement::NumberOfNodes); comboIndex++){
      size_t number = comboIndex;
      vector<size_t> combo;
      while (number >= PrimaryElement::NumberOfNodes){
        size_t remainder = number%(PrimaryElement::NumberOfNodes);
        number = floor(number/(PrimaryElement::NumberOfNodes));
        combo.push_back(remainder);
      }
      combo.push_back(number);
      while (combo.size() < PrimaryElement::NumberOfNodes){
        combo.push_back(0);
      }
      // now let's make a fake connection out of the combo
      array<size_t,PrimaryElement::NumberOfNodes> tempConnection;
      for (size_t comboSubIndex = 0; comboSubIndex < PrimaryElement::NumberOfNodes; comboSubIndex++){
        if (comboSubIndex < combo.size()){
          tempConnection[comboSubIndex] = combo[comboSubIndex];
        }
        else{
          tempConnection[comboSubIndex] = 0;
        }
      }
      possibleCombinationsPrimary.push_back(tempConnection);
    }
  }
  // now that we have all the possible node combinations for the same element
  // in space, let's find if any other elements in the matrix match this element
  for (size_t elementIndex1 = 0; elementIndex1 < numberOfPrimaryElements;
       elementIndex1++){
    for (size_t elementIndex2 = elementIndex1; elementIndex2 < numberOfPrimaryElements;
         elementIndex2++){
      for (size_t comboIndex = 0; comboIndex <
             pow(PrimaryElement::NumberOfNodes,PrimaryElement::NumberOfNodes); comboIndex++){
        size_t elementFlag = 0;
        // turn possibleCombination into node number of elementIndex2
        array<size_t,PrimaryElement::NumberOfNodes> possibleConnection;
        for (size_t nodeIndex = 0; nodeIndex < PrimaryElement::NumberOfNodes; nodeIndex++){
          const size_t nodeToUse = possibleCombinationsPrimary[comboIndex][nodeIndex];
          possibleConnection[nodeIndex] = mesh->_connectivity1[elementIndex1][nodeToUse];
        }
        // check if the possibleConnection is the same as the 1st element
        for (size_t nodeIndex = 0; nodeIndex < PrimaryElement::NumberOfNodes; nodeIndex++){
          if (possibleConnection[nodeIndex] == mesh->_connectivity1[elementIndex2][nodeIndex]){
            elementFlag++;
          }
        }
        if (elementFlag == PrimaryElement::NumberOfNodes &&
            elementIndex1 != elementIndex2){
          primaryElementsMarkedForDeletion[elementIndex2] = 1;
        }
      }
      for (size_t checkIndex = 0; checkIndex < numberOfPrimaryElements; checkIndex++){
      }
    }
    for (size_t elementIndex = 0; elementIndex < numberOfPrimaryElements; elementIndex++){
      if (primaryElementsMarkedForDeletion[elementIndex] == 1){
        mesh->_connectivity1.erase(mesh->_connectivity1.begin() + elementIndex);
        primaryElementsMarkedForDeletion.erase(primaryElementsMarkedForDeletion.begin()
                                               + elementIndex);
        numberOfPrimaryElements--;
        if (elementIndex > 0){
          elementIndex--;
        }
      }
    }
  }


  // SECONDARY ELEMENT //
  vector<size_t> secondaryElementsMarkedForDeletion;
  // just fill the vector to keep track of who to delete later
  for (size_t elementIndex = 0; elementIndex < numberOfSecondaryElements; elementIndex++){
    secondaryElementsMarkedForDeletion.push_back(0);
  }

  vector<array<size_t,SecondaryElement::NumberOfNodes>> possibleCombinationsSecondary;
  // now let's look for elements with the same node numbers
  for (size_t elementIndex = 0; elementIndex < numberOfSecondaryElements; elementIndex++){
    for (size_t comboIndex = 0;
         comboIndex < pow(SecondaryElement::NumberOfNodes,
                          SecondaryElement::NumberOfNodes); comboIndex++){
      size_t number = comboIndex;
      vector<size_t> combo;
      while (number >= SecondaryElement::NumberOfNodes){
        size_t remainder = number%(SecondaryElement::NumberOfNodes);
        number = floor(number/(SecondaryElement::NumberOfNodes));
        combo.push_back(remainder);
      }
      combo.push_back(number);
      while (combo.size() < SecondaryElement::NumberOfNodes){
        combo.push_back(0);
      }
      // now let's make a fake connection out of the combo
      array<size_t,SecondaryElement::NumberOfNodes> tempConnection;
      for (size_t comboSubIndex = 0; comboSubIndex < SecondaryElement::NumberOfNodes; comboSubIndex++){
        if (comboSubIndex < combo.size()){
          tempConnection[comboSubIndex] = combo[comboSubIndex];
        }
        else{
          tempConnection[comboSubIndex] = 0;
        }
      }
      possibleCombinationsSecondary.push_back(tempConnection);
    }
  }
  // now that we have all the possible node combinations for the same element
  // in space, let's find if any other elements in the matrix match this element
  for (size_t elementIndex1 = 0; elementIndex1 < numberOfSecondaryElements;
       elementIndex1++){
    for (size_t elementIndex2 = elementIndex1; elementIndex2 < numberOfSecondaryElements;
         elementIndex2++){
      for (size_t comboIndex = 0; comboIndex <
             pow(SecondaryElement::NumberOfNodes,SecondaryElement::NumberOfNodes); comboIndex++){
        size_t elementFlag = 0;
        // turn possibleCombination into node number of elementIndex2
        array<size_t,SecondaryElement::NumberOfNodes> possibleConnection;
        for (size_t nodeIndex = 0; nodeIndex < SecondaryElement::NumberOfNodes; nodeIndex++){
          const size_t nodeToUse = possibleCombinationsSecondary[comboIndex][nodeIndex];
          possibleConnection[nodeIndex] = mesh->_connectivity2[elementIndex1][nodeToUse];
        }
        // check if the possibleConnection is the same as the 1st element
        for (size_t nodeIndex = 0; nodeIndex < SecondaryElement::NumberOfNodes; nodeIndex++){
          if (possibleConnection[nodeIndex] == mesh->_connectivity2[elementIndex2][nodeIndex]){
            elementFlag++;
          }
        }
        if (elementFlag == SecondaryElement::NumberOfNodes &&
            elementIndex1 != elementIndex2){
          secondaryElementsMarkedForDeletion[elementIndex2] = 1;
        }
      }
    }
    for (size_t elementIndex = 0; elementIndex < numberOfSecondaryElements; elementIndex++){
      if (secondaryElementsMarkedForDeletion[elementIndex] == 1){
        mesh->_connectivity2.erase(mesh->_connectivity2.begin() + elementIndex);
        secondaryElementsMarkedForDeletion.erase(secondaryElementsMarkedForDeletion.begin() + elementIndex);
        numberOfSecondaryElements--;
        if (elementIndex > 0){
          elementIndex--;
        }
      }
    }
  }
}

template <class Element>
void
buildPeriodicMeshFromUnitCell(const array<Matrix<double,Element::SpatialDimension,1>
                              ,Element::SpatialDimension> patternVectors,
                              const array<size_t,Element::SpatialDimension> instances,
                              SingleElementMesh<Element> * mesh,
                              const double spatialTolerance = 1e-3){

  SingleElementMesh<Element> originalMesh = (*mesh);
  if (Element::SpatialDimension == 3){
    for (size_t indexX = 0; indexX <= instances[0]; indexX++){
      for (size_t indexY = 0; indexY <= instances[1]; indexY++){
        for (size_t indexZ = 0; indexZ <= instances[2]; indexZ++){
          Matrix<double,Element::SpatialDimension,1> translation =
            patternVectors[0] * indexX;
            translation += patternVectors[1] * indexY;
            translation += patternVectors[2] * indexZ;
          // TODO need to fix for 3d
          translateAndCopyMesh(translation,originalMesh,mesh);
        }
      }
    }
  }
  else if (Element::SpatialDimension == 2){
    for (size_t indexX = 0; indexX <= instances[0]; indexX++){
      for (size_t indexY = 0; indexY <= instances[1]; indexY++){
        Matrix<double,Element::SpatialDimension,1> translation =
          patternVectors[0] * indexX;
        cout << indexX << "," << indexY << endl;
        translation += patternVectors[1] * indexY;
        translateAndCopyMesh(translation,originalMesh,mesh);
      }
    }
  }

  // search for overlapping nodes and remove overlapping elements
  removeOverlapingNodesAndElements(mesh,spatialTolerance);
}

template <class PrimaryElement, class SecondaryElement>
void
buildPeriodicMeshFromUnitCell(const array<Matrix<double,PrimaryElement::SpatialDimension,1>
                              ,PrimaryElement::SpatialDimension> patternVectors,
                              const array<size_t,PrimaryElement::SpatialDimension> instances,
                              TwoElementMesh<PrimaryElement,SecondaryElement> * mesh){

  TwoElementMesh<PrimaryElement,SecondaryElement> originalMesh = (*mesh);
  if (PrimaryElement::SpatialDimension == 2){
    for (size_t indexX = 0; indexX <= instances[0]; indexX++){
      for (size_t indexY = 0; indexY <= instances[1]; indexY++){
        Matrix<double,PrimaryElement::SpatialDimension,1> translation =
          patternVectors[0] * indexX;
        if (patternVectors.size() == 2){
          translation += patternVectors[1] * indexY;
        }
        // TODO need to fix for 3d
        cout << "translating and copying mesh" << endl;
        translateAndCopyMesh(translation,originalMesh,mesh);
      }
    }
  }
  if (PrimaryElement::SpatialDimension == 3){
    for (size_t indexX = 0; indexX <= instances[0]; indexX++){
      for (size_t indexY = 0; indexY <= instances[1]; indexY++){
        for (size_t indexZ = 0; indexZ <= instances[2]; indexZ++){
          Matrix<double,PrimaryElement::SpatialDimension,1> translation =
            patternVectors[0] * indexX;
          if(patternVectors.size() == 3){
            translation += patternVectors[1] * indexY;
            translation += patternVectors[2] * indexZ;
          }
          // TODO need to fix for 3d
          cout << "translating and copying mesh" << endl;
          translateAndCopyMesh(translation,originalMesh,mesh);
        }
      }
    }
  }

  // search for overlapping nodes and remove overlapping elements
  cout << "removing overalpping nodes" << endl;
  removeOverlapingNodesAndElements(mesh);
}


template <class Element>
void
buildBrickMeshedCuboid(const array<double, 3> & sideLengths,
                       const array<size_t, 3> & numberOfElements,
                       SingleElementMesh<Element> * mesh) {
  typedef typename Element::Node Node;

  const double a = sideLengths[0];
  const double b = sideLengths[1];
  const double c = sideLengths[2];
  const size_t n = numberOfElements[0];
  const size_t m = numberOfElements[1];
  const size_t l = numberOfElements[2];

  double dx = a/double(n);
  double dy = b/double(m);
  double dz = c/double(l);

  size_t nextNodeId = 0;
  Eigen::Vector3d point;
  for (size_t k=0; k<=l; k++) {
    for (size_t j=0; j<=m; j++) {
      for (size_t i=0; i<=n; i++) {
        point[0] = double(i * dx);
        point[1] = double(j * dy);
        point[2] = double(k * dz);
        mesh->_nodes.push_back(Node(nextNodeId, point));
        ++nextNodeId;
      }
    }
  }

  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<m; j++) {
      for (size_t k=0; k<l; k++) {
        array<size_t, 8> connection;
        size_t bottomLeftCorner = k*(n+1)*(m+1)+j*(n+1)+i;
        connection[0] = bottomLeftCorner;
        connection[1] = bottomLeftCorner+1;
        connection[2] = bottomLeftCorner+(n+1)+1;
        connection[3] = bottomLeftCorner+(n+1);
        connection[4] = bottomLeftCorner+(n+1)*(m+1);
        connection[5] = bottomLeftCorner+(n+1)*(m+1)+1;
        connection[6] = bottomLeftCorner+(n+1)*(m+1)+(n+1)+1;
        connection[7] = bottomLeftCorner+(n+1)*(m+1)+(n+1);
        mesh->_connectivity.push_back(connection);
      }
    }
  }
}

template <class Element>
void
buildTetMeshedCuboid(const array<double, 3> & sideLengths,
                     const array<size_t, 3> & numberOfElements,
                     SingleElementMesh<Element> * mesh) {

  if (Element::SpatialDimension != 3) {
    cout << "Error: Using wrong dimension element for building a tet mesh."<< endl;
    exit(1);
  }

  typedef typename Element::Point Point;
  typedef typename Element::Node Node;

  const double a = sideLengths[0];
  const double b = sideLengths[1];
  const double c = sideLengths[2];
  const size_t n = numberOfElements[0];
  const size_t m = numberOfElements[1];
  const size_t L = numberOfElements[2];

  double dx = a/double(n);
  double dy = b/double(m);
  double dz = c/double(L);

  size_t nextNodeId = 0;
  Point point;
  // add the corner nodes
  for (size_t k=0; k<=L; k++) {
    for (size_t j=0; j<=m; j++) {
      for (size_t i=0; i<=n; i++) {
        point[0] = double(i * dx);
        point[1] = double(j * dy);
        point[2] = double(k * dz);
        mesh->_nodes.push_back(Node(nextNodeId, point));
        ++nextNodeId;
      }
    }
  }

  // add the center nodes
  for (size_t k=0; k<L; k++) {
    for (size_t j=0; j<m; j++) {
      for (size_t i=0; i<n; i++) {
        point[0] = double(i * dx + 0.5 * dx);
        point[1] = double(j * dy + 0.5 * dy);
        point[2] = double(k * dz + 0.5 * dz);
        mesh->_nodes.push_back(Node(nextNodeId, point));
        ++nextNodeId;
      }
    }
  }

  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<m; j++) {
      for (size_t k=0; k<L; k++) {

        size_t bottomLeftCorner = k*(n+1)*(m+1)+j*(n+1)+i;
        size_t centerNode = (n+1)*(m+1)*(L+1) + k*n*m+j*n+i;
        size_t b1 = bottomLeftCorner;
        size_t b2 = bottomLeftCorner+1;
        size_t b3 = bottomLeftCorner+(n+1)+1;
        size_t b4 = bottomLeftCorner+(n+1);
        size_t t1 = bottomLeftCorner+(n+1)*(m+1);
        size_t t2 = bottomLeftCorner+(n+1)*(m+1)+1;
        size_t t3 = bottomLeftCorner+(n+1)*(m+1)+(n+1)+1;
        size_t t4 = bottomLeftCorner+(n+1)*(m+1)+(n+1);
        size_t c = centerNode;

        // bottom two tetrahedrons
        mesh->_connectivity.push_back((array<size_t, 4>) {{b1,b2,b3,c}});
        mesh->_connectivity.push_back((array<size_t, 4>) {{b4,b1,b3,c}});

        // top two tetrahedrons
        mesh->_connectivity.push_back((array<size_t, 4>) {{t3,t2,t1,c}});
        mesh->_connectivity.push_back((array<size_t, 4>) {{t4,t3,t1,c}});

        // front two tetrahedrons
        mesh->_connectivity.push_back((array<size_t, 4>) {{t2,b2,b1,c}});
        mesh->_connectivity.push_back((array<size_t, 4>) {{t1,t2,b1,c}});

        // back two tetrahedrons
        mesh->_connectivity.push_back((array<size_t, 4>) {{b3,t3,b4,c}});
        mesh->_connectivity.push_back((array<size_t, 4>) {{t3,t4,b4,c}});

        // left two tetrahedrons
        mesh->_connectivity.push_back((array<size_t, 4>) {{b1,b4,t4,c}});
        mesh->_connectivity.push_back((array<size_t, 4>) {{t1,b1,t4,c}});

        // right two tetrahedrons
        mesh->_connectivity.push_back((array<size_t, 4>) {{b3,b2,t3,c}});
        mesh->_connectivity.push_back((array<size_t, 4>) {{t3,b2,t2,c}});

        // if eigen complains, use these
        /*
          array<size_t, 4> temp;
          temp[0] = b1;
          temp[1] = b2;
          temp[2] = b4;
          temp[3] = c;
          mesh->_connectivity.push_back(temp);
          temp[0] = b2;
          temp[1] = b3;
          temp[2] = b4;
          temp[3] = c;
          mesh->_connectivity.push_back(temp);
          temp[0] = t2;
          temp[1] = t4;
          temp[2] = t3;
          temp[3] = c;
          mesh->_connectivity.push_back(temp);
          temp[0] = t2;
          temp[1] = t1;
          temp[2] = t4;
          temp[3] = c;
          mesh->_connectivity.push_back(temp);
        */

      }
    }
  }
}


/*
  template <class Element>
  void
  buildTetMeshedCuboid(const array<double, 3> & sideLengths,
  const array<size_t, 3> & numberOfElements,
  SingleElementMesh<Element> * mesh) {

  if (Element::SpatialDimension != 3) {
  cout << "Error: Using wrong dimension element for building a tet mesh."<< endl;
  exit(1);
  }

  typedef typename Element::Point Point;
  typedef NodeWithId<Point> Node;

  const double a = sideLengths[0];
  const double b = sideLengths[1];
  const double c = sideLengths[2];
  const size_t n = numberOfElements[0];
  const size_t m = numberOfElements[1];
  const size_t l = numberOfElements[2];

  const double dx = a/double(n);
  const double dy = b/double(m);
  const double dz = c/double(l);

  size_t nextNodeId = 0;
  Point point;
  // add corner nodes
  for (size_t i=0; i<=n; i++) {
  for (size_t j=0; j<=m; j++) {
  for (size_t k=0; k<=l; k++) {
  point[0] = double(i * dx);
  point[1] = double(j * dy);
  point[2] = double(k * dz);
  mesh->_nodes.push_back(Node(nextNodeId, point));
  ++nextNodeId;
  }
  }
  }

  // connectivities
  for (size_t i=0; i<n; i++) {
  for (size_t j=0; j<m; j++) {
  for (size_t k=0; k<l; k++) {
  const size_t bottomLeftCorner = i*(l+1)*(m+1)+j*(l+1)+k;
  const size_t b1 = bottomLeftCorner;
  const size_t t1 = bottomLeftCorner+1;
  const size_t b4 = bottomLeftCorner+(l+1);
  const size_t t4 = bottomLeftCorner+(l+1)+1;
  const size_t b2 = bottomLeftCorner+(l+1)*(m+1);
  const size_t t2 = bottomLeftCorner+(l+1)*(m+1)+1;
  const size_t b3 = bottomLeftCorner+(l+1)*(m+1)+(l+1);
  const size_t t3 = bottomLeftCorner+(l+1)*(m+1)+(l+1)+1;

  mesh->_connectivity.push_back((array<size_t, 4>) {{b1,t1,t2,t3}});
  mesh->_connectivity.push_back((array<size_t, 4>) {{b1,t2,b2,t3}});
  mesh->_connectivity.push_back((array<size_t, 4>) {{b1,b2,b3,t3}});
  mesh->_connectivity.push_back((array<size_t, 4>) {{b1,b3,b4,t3}});
  mesh->_connectivity.push_back((array<size_t, 4>) {{b1,b4,t4,t3}});
  mesh->_connectivity.push_back((array<size_t, 4>) {{b1,t1,t3,t4}});
  }
  }
  }
  }
*/
template <class Element>
void
buildRectangularQuadMesh(const array<double, 2> & sideLengths,
                         const array<size_t, 2> & numberOfElements,
                         SingleElementMesh<Element> * mesh) {

  if (Element::SpatialDimension != 2){
    cout << "Error: Using wrong dimension element for building a rectangular quad mesh." << endl;
    exit(1);
  }

  typedef typename Element::Node  Node;

  const double a = sideLengths[0];
  const double b = sideLengths[1];
  const size_t n = numberOfElements[0];
  const size_t m = numberOfElements[1];

  double dx = a/double(n);
  double dy = b/double(m);

  size_t nextNodeId = 0;
  Eigen::Vector2d point;
  for (size_t j=0; j<=m; j++) {
    for (size_t i=0; i<=n; i++) {
      point[0] = double(i * dx);
      point[1] = double(j * dy);
      mesh->_nodes.push_back(Node(nextNodeId, point));
      ++nextNodeId;
    }
  }
  for (size_t j=0; j<m; j++) {
    for (size_t i=0; i<n; i++) {
      array<size_t, 4> connection;
      size_t bottomLeftCorner = j*(n+1)+i;
      connection[0] = bottomLeftCorner;
      connection[1] = bottomLeftCorner+1;
      connection[2] = bottomLeftCorner+(n+1)+1;
      connection[3] = bottomLeftCorner+(n+1);
      /*  cout << mesh->_nodes[bottomLeftCorner]._position << endl;
          cout << mesh->_nodes[bottomLeftCorner+1]._position << endl;
          cout << mesh->_nodes[bottomLeftCorner+n+2]._position << endl;
          cout << mesh->_nodes[bottomLeftCorner+n+1]._position << endl << endl; */
      mesh->_connectivity.push_back(connection);
    }
  }
}

template <class Element>
void
buildRectangularTriangleMesh(const array<double, 2> & sideLengths,
                             const array<size_t, 2> & numberOfElements,
                             SingleElementMesh<Element> * mesh) {

  typedef Matrix<double, 2, 1> Point;
  typedef NodeWithId<Point> Node;

  const float a = sideLengths[0];
  const float b = sideLengths[1];
  const size_t n = numberOfElements[0];
  const size_t m = numberOfElements[1];

  float dx = a/float(n);
  float dy = b/float(m);

  size_t nextNodeId = 0;
  Eigen::Vector2d point;
  for (size_t j=0; j<=m; j++) {
    for (size_t i=0; i<=n; i++) {
      point[0] = float(i * dx);
      point[1] = float(j * dy);
      mesh->_nodes.push_back(Node(nextNodeId, point));
      ++nextNodeId;
    }
  }
  for (size_t j=0; j<m; j++) {
    for (size_t i=0; i<n; i++) {
      array<size_t, 3> connection;
      size_t bottomLeftCorner = j*(n+1)+i;
      connection[0] = bottomLeftCorner;
      connection[1] = bottomLeftCorner+1;
      connection[2] = bottomLeftCorner+(n+1)+1;
      mesh->_connectivity.push_back(connection);
      connection[0] = bottomLeftCorner+(n+1)+1;
      connection[1] = bottomLeftCorner+(n+1);
      connection[2] = bottomLeftCorner;
      mesh->_connectivity.push_back(connection);
    }
  }
}

template <class Element>
void
buildRectangularTriangleMeshIn3D(const array<double, 2> & sideLengths,
                             const array<size_t, 2> & numberOfElements,
                             SingleElementMesh<Element> * mesh) {

  typedef Matrix<double, 3, 1> Point;
  typedef NodeWithId<Point> Node;

  const float a = sideLengths[0];
  const float b = sideLengths[1];
  const size_t n = numberOfElements[0];
  const size_t m = numberOfElements[1];

  float dx = a/float(n);
  float dy = b/float(m);

  size_t nextNodeId = 0;
  Eigen::Vector3d point;
  for (size_t j=0; j<=m; j++) {
    for (size_t i=0; i<=n; i++) {
      point[0] = float(i * dx);
      point[1] = float(j * dy);
      point[2] = 0.0;
      mesh->_nodes.push_back(Node(nextNodeId, point));
      ++nextNodeId;
    }
  }
  for (size_t j=0; j<m; j++) {
    for (size_t i=0; i<n; i++) {
      array<size_t, 3> connection;
      size_t bottomLeftCorner = j*(n+1)+i;
      connection[0] = bottomLeftCorner;
      connection[1] = bottomLeftCorner+1;
      connection[2] = bottomLeftCorner+(n+1)+1;
      mesh->_connectivity.push_back(connection);
      connection[0] = bottomLeftCorner+(n+1)+1;
      connection[1] = bottomLeftCorner+(n+1);
      connection[2] = bottomLeftCorner;
      mesh->_connectivity.push_back(connection);
    }
  }
}

// NOT TESTED AT ALL!!!!

#if 0
void
buildBowTieMesh(const double L, const double Theta, const int ElementsPerLength, SingleElementMesh<NodeWithId<Vector3d>, 2> * mesh){

  typedef NodeWithId<Vector3d> Node;

  vector<Vector3d> Directions;
  vector<Vector3d> startingPoints;
  double a = sqrt(2.) * sin(Theta) * L;
  //double c = 2*cos(Theta)*L;
  double dz = L *(1-cos(Theta));
  Vector3d P1, P2, P3, P4,P5;
  P1<< 0,a,L-dz;
  startingPoints.push_back(P1);
  P2<< a,a,L-dz;
  startingPoints.push_back(P2);
  P3<< a,0,L-dz;
  startingPoints.push_back(P3);
  P4<< 0,0,L-dz;
  startingPoints.push_back(P4);
  P5<< a/2.,a/2.,L;
  startingPoints.push_back(P5);
  Vector3d disp;
  disp <<L,L,L;
  Vector3d center;
  center <<a/2.,a/2,0;



  Vector3d Direction1, Direction2, Direction3, Direction4,Direction5;
  Direction1=center-P1;
  Directions.push_back(Direction1);
  Direction2=center-P2;
  Directions.push_back(Direction2);
  Direction3=center-P3;
  Directions.push_back(Direction3);
  Direction4=center-P4;
  Directions.push_back(Direction4);
  Direction5=center-P5;
  Directions.push_back(Direction5);
  int nextNodeId = 0;

  for(int i = 0; i < 5 ; i++){
    Vector3d dir = Directions[i];
    Vector3d start= startingPoints[i];
    for(int j = 0 ; j < ElementsPerLength+1; j++){
      if(i==0){
        Vector3d Point = start + j * dir / ElementsPerLength;

        mesh->_nodes.push_back(Node(nextNodeId, Point));
        ++nextNodeId;
      }
      if(i>0){
        Vector3d Point = start + j * dir / ElementsPerLength;
        //if((Point-center).norm() < 1e-5){cout << "found another center node"<<endl;}
        if((Point-center).norm() > 1e-5){

          mesh->_nodes.push_back(Node(nextNodeId, Point));
          ++nextNodeId;
        }
        //else{cout << "This is a center Point " << endl<<Point<<endl<<endl;}
      }
    }

  }


  //Build Connectivity
  int ID = 0;
  for(unsigned int k = 0; k < 5; k++){

    ID = ElementsPerLength * k + 1;
    if(k==0){ID=0;}
    int count = 0;
    for(int j = 0; j < ElementsPerLength;j++){
      if(k==0){
        mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ID+1}});

        ID++;

      }
      else{

        if(count == (ElementsPerLength-1)){
          mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ElementsPerLength}});
          //mesh->_connectivity.push_back((array<size_t, 2>) {{ElementsPerLength, ID+1}});

          //cout << "Connecting Nodes " << ElementsPerLength << " " <<ID+1 <<endl;
          ID++;
          j++;
          count += 1;
        }
        else{

          mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ID+1}});

          ID++;
          count++;
        }
      }
    }


  }

}
#endif


// NOT TESTED AT ALL!!!!
#if 0

void
buildMicrolatticeMesh(const double dx, const double dy, const double dz, const unsigned int ElementsPerLength, SingleElementMesh<NodeWithId<Vector3d>, 2> * mesh){
  typedef NodeWithId<Vector3d> Node;

  vector<Vector3d> Directions;
  vector<Vector3d> startingPoints;

  Vector3d P1, P2, P3, P4;
  P1<< 0.,0.,0.;
  startingPoints.push_back(P1);
  P2<< dx,0.,0.;
  startingPoints.push_back(P2);
  P3<< dx,dy,0.;
  startingPoints.push_back(P3);
  P4<< 0.,dy,0.;
  startingPoints.push_back(P4);

  Vector3d center;
  center <<dx/2.,dy/2.,dz/2.;

  Vector3d Direction1, Direction2, Direction3, Direction4;
  Direction1<<dx,dy,dz;
  Directions.push_back(Direction1);
  Direction2<<-dx,dy,dz;
  Directions.push_back(Direction2);
  Direction3<<-dx,-dy,dz;
  Directions.push_back(Direction3);
  Direction4<<dx,-dy,dz;
  Directions.push_back(Direction4);
  int nextNodeId = 0;

  for(int i = 0; i < 4 ; i++){
    Vector3d dir = Directions[i];
    Vector3d start= startingPoints[i];
    for(unsigned int j = 0 ; j < ElementsPerLength+1; j++){
      if(i==0){
        Vector3d Point = start + j * dir / ElementsPerLength;
        //cout << "Point in i = 0 Loop " << endl<<Point<<endl<<endl;
        mesh->_nodes.push_back(Node(nextNodeId, Point));
        ++nextNodeId;
      }
      else{
        Vector3d Point = start + j * dir / ElementsPerLength;
        if(Point != center){
          cout << " This is not the center Point at j got to push this back" << j <<endl;
          mesh->_nodes.push_back(Node(nextNodeId, Point));
          ++nextNodeId;
        }
        //else{cout << "This is a center Point " << endl<<Point<<endl<<endl;}
      }
    }

  }
  cout << "We have overall: " << mesh->_nodes.size()<<" Nodes"<<endl;
  cout << "These are the Nodes: " << endl;
  for (unsigned int h = 0; h < mesh->_nodes.size();h++){
    Node node = mesh->_nodes[h];
    Vector3d pos=node._position;
    cout <<"ID "<<node._id<< " x " << pos(0) << " y " << pos(1) << " z " << pos(2) <<endl;
  }

  //Build Connectivity
  unsigned int ID = 0;
  for(unsigned int k = 0; k < 4; k++){
    cout << "beginning k = " << k <<endl;
    ID = ElementsPerLength * k + 1;
    if(k==0){ID=0;}
    int count = 0;
    for(unsigned int j = 0; j < ElementsPerLength;j++){
      if(k==0){
        mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ID+1}});
        cout << "Connecting Nodes in k=0 loop: " << ID << " and " << ID+1<<endl;
        ID++;

      }
      else{

        if(count == (ElementsPerLength/2-1)){
          mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ElementsPerLength/2}});
          mesh->_connectivity.push_back((array<size_t, 2>) {{ElementsPerLength/2, ID+1}});
          cout << "Connecting Nodes in k loop with center Nodes:" << ID << " and " << ElementsPerLength/2<<endl;
          cout << "Connecting Nodes in k loop with center Nodes:" << ElementsPerLength/2 << " and " <<ID+1 <<endl;
          ID++;
          j++;
          count += 2;
        }
        else{
          cout << "Connecting Nodes in k loop without center Nodes:" << ID << " and " << ID+1<<endl;
          mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ID+1}});
          ID++;
          count++;
        }
      }
    }


  }

}
// NOT TESTED AT ALL!!!!
#endif
// students, do NOT use this; it doesn't even work.
#if 0
void
buildOctahedronMesh(const double L, const double Theta, const int ElementsPerLength, SingleElementMesh<NodeWithId<Vector3d>, 2> * mesh){

  typedef NodeWithId<Vector3d> Node;

  vector<Vector3d> Directions;
  vector<Vector3d> startingPoints;
  double a = sqrt(2.) * cos(Theta) * L;
  double c = 2*sin(Theta)*L;
  Vector3d P1, P2, P3, P4,P5, P6,P7, P8;
  P1<< 0.,0.,0.;
  startingPoints.push_back(P1);
  P2<< a,0.,0.;
  startingPoints.push_back(P2);
  P3<< a,a,0.;
  startingPoints.push_back(P3);
  P4<< 0.,a,0.;
  startingPoints.push_back(P4);
  P5=P1;
  startingPoints.push_back(P5);
  P6=P2;
  startingPoints.push_back(P6);
  P7=P3;
  startingPoints.push_back(P7);
  P8=P4;
  startingPoints.push_back(P8);

  Vector3d center;
  center <<a/2.,a/2.,c/2.;

  Vector3d Direction1, Direction2, Direction3, Direction4,Direction5, Direction6,Direction7, Direction8;
  Direction1<<a/2,a/2,c/2.;
  Directions.push_back(Direction1);
  Direction2<<-a/2,a/2,c/2;
  Directions.push_back(Direction2);
  Direction3<<-a/2,-a/2,c/2;
  Directions.push_back(Direction3);
  Direction4<<a/2,-a/2,c/2;
  Directions.push_back(Direction4);
  Direction5<<a,0,0;
  Directions.push_back(Direction5);
  Direction6<<0,a,0;
  Directions.push_back(Direction6);
  Direction7<<-a,0,0;
  Directions.push_back(Direction7);
  Direction8<<0,-a,0;
  Directions.push_back(Direction8);


  int nextNodeId = 0;

  for(int i = 0; i < 8 ; i++){
    Vector3d dir = Directions[i];
    Vector3d start= startingPoints[i];
    for(int j = 0 ; j < ElementsPerLength+1; j++){
      if(i==0){
        Vector3d Point = start + j * dir / ElementsPerLength;

        mesh->_nodes.push_back(Node(nextNodeId, Point));
        cout << "This is a Point " << endl<<Point.transpose()<<endl<<" at ID " << nextNodeId<< endl;
        ++nextNodeId;
      }
      if((i>0) && (i < 4)){
        Vector3d Point = start + j * dir / ElementsPerLength;
        //if((Point-center).norm() < 1e-5){cout << "found another center node"<<endl;}
        if((Point-center).norm() > 1e-5){

          mesh->_nodes.push_back(Node(nextNodeId, Point));
          cout << "This is a Point " << endl<<Point.transpose()<<endl<<" at ID " << nextNodeId<< endl;
          ++nextNodeId;
        }
        //else{cout << "This is a center Point " << endl<<Point<<endl<<endl;}
      }
      if(i>=4){
        Vector3d Point = start + j * dir / ElementsPerLength;
        if(((Point-startingPoints[0]).norm() > 1e-5) && ((Point-startingPoints[1]).norm() > 1e-5) && ((Point-startingPoints[2]).norm() > 1e-5) && ((Point-startingPoints[3]).norm() > 1e-5)  ){

          mesh->_nodes.push_back(Node(nextNodeId, Point));
          cout << "This is a Point " << endl<<Point.transpose()<<endl<<" at ID " << nextNodeId<< endl;
          ++nextNodeId;
        }

      }
    }

  }


  //Build Connectivity
  int ID = 0;
  for(unsigned int k = 0; k < 4; k++){

    ID = ElementsPerLength * k + 1;
    if(k==0){ID=0;}
    int count = 0;
    for(int j = 0; j < ElementsPerLength;j++){
      if(k==0){
        mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ID+1}});
        cout << "Connecting Nodes " << ID << " " <<ID+1 <<endl;
        ID++;

      }
      else{

        if(count == (ElementsPerLength-1)  ){
          mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ElementsPerLength}});
          //mesh->_connectivity.push_back((array<size_t, 2>) {{ElementsPerLength, ID+1}});

          cout << "Connecting Nodes " << ID << " " <<ElementsPerLength <<endl;
          ID++;
          j++;
          count += 1;
        }
        else{

          mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ID+1}});
          cout << "Connecting Nodes " << ID << " " <<ID +1 <<endl;
          ID++;
          count++;
        }
      }
    }


  }
  for(unsigned int k = 4; k < 8; k++){
    //ID = ElementsPerLength * k + 1;
    int count = 0;

    for(int j = 0; j < ElementsPerLength;j++){
      if(count == 0){
        int IDprime = (ElementsPerLength)*(k-4)+1;
        if (IDprime == 1){IDprime=0;};
        mesh->_connectivity.push_back((array<size_t, 2>) {{IDprime,ID}});
        //mesh->_connectivity.push_back((array<size_t, 2>) {{ElementsPerLength, ID+1}});

        cout << "Connecting Nodes " << IDprime << " first loop " <<ID <<endl;
        //ID++;
        j++;
        count += 1;
      }
      if(count == (ElementsPerLength-1)){
        int IDprime = (ElementsPerLength)*(k-3)+1;
        if (IDprime==(ElementsPerLength)*(7-3)+1){IDprime = 0;};
        mesh->_connectivity.push_back((array<size_t, 2>) {{ID, IDprime}});
        //mesh->_connectivity.push_back((array<size_t, 2>) {{ElementsPerLength, ID+1}});

        cout << "Connecting Nodes " << ID << "  second loop " <<IDprime <<endl;
        //ID++;
        j++;
        count += 1;
      }
      else{

        mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ID+1}});
        cout << "Connecting Nodes " << ID << " third loop " <<ID+1 <<endl;
        ID++;
        count++;
      }
    }ID++;
  }
}
#endif

// NOT TESTED AT ALL!!!!
#if 0
void
buildMicrolatticeLocalResonanceMesh(const double L, const double theta,const int ElementsPerLength, SingleElementMesh<NodeWithId<Vector3d>, 2> * mesh){


  typedef NodeWithId<Vector3d> Node;

  vector<Vector3d> Directions;
  vector<Vector3d> startingPoints;
  double a = sqrt(2.) * cos(theta) * L;
  //double c = 2*sin(theta)*L;
  // double dz = L *(1-cos(theta));
  Vector3d P1, P2, P3, P4,P5;
  P1<< 0,a,sin(theta)*L;
  startingPoints.push_back(P1);
  P2<< a,a,sin(theta)*L;
  startingPoints.push_back(P2);
  P3<< a,0,sin(theta)*L;
  startingPoints.push_back(P3);
  P4<< 0,0,sin(theta)*L;
  startingPoints.push_back(P4);
  P5<< a/2.,a/2.,sin(theta)*L;
  startingPoints.push_back(P5);
  Vector3d disp;
  disp <<L,L,L;
  Vector3d center;
  center <<a/2.,a/2,0;



  Vector3d Direction1, Direction2, Direction3, Direction4,Direction5;
  Direction1=center-P1;
  Directions.push_back(Direction1);
  Direction2=center-P2;
  Directions.push_back(Direction2);
  Direction3=center-P3;
  Directions.push_back(Direction3);
  Direction4=center-P4;
  Directions.push_back(Direction4);
  Direction5=center-P5;
  Directions.push_back(Direction5);
  int nextNodeId = 0;

  for(int i = 0; i < 5 ; i++){
    Vector3d dir = Directions[i];
    Vector3d start= startingPoints[i];
    for(int j = 0 ; j < ElementsPerLength+1; j++){
      if(i==0){
        Vector3d Point = start + j * dir / ElementsPerLength;

        mesh->_nodes.push_back(Node(nextNodeId, Point));
        ++nextNodeId;
      }
      if(i>0){
        Vector3d Point = start + j * dir / ElementsPerLength;
        //if((Point-center).norm() < 1e-5){cout << "found another center node"<<endl;}
        if((Point-center).norm() > 1e-5){

          mesh->_nodes.push_back(Node(nextNodeId, Point));
          ++nextNodeId;
        }
        //else{cout << "This is a center Point " << endl<<Point<<endl<<endl;}
      }
    }

  }


  //Build Connectivity
  int ID = 0;
  for(unsigned int k = 0; k < 5; k++){

    ID = ElementsPerLength * k + 1;
    if(k==0){ID=0;}
    int count = 0;
    for(int j = 0; j < ElementsPerLength;j++){
      if(k==0){
        mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ID+1}});

        ID++;

      }
      else{

        if(count == (ElementsPerLength-1)){
          mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ElementsPerLength}});
          //mesh->_connectivity.push_back((array<size_t, 2>) {{ElementsPerLength, ID+1}});

          //cout << "Connecting Nodes " << ElementsPerLength << " " <<ID+1 <<endl;
          ID++;
          j++;
          count += 1;
        }
        else{

          mesh->_connectivity.push_back((array<size_t, 2>) {{ID, ID+1}});

          ID++;
          count++;
        }
      }
    }


  }

}
#endif

#if 0
void
buildRectangularQuadMesh(const array<double, 2> & sideLengths,
                         const array<size_t, 2> & numberOfElements,
                         SingleElementMesh<NodeWithId<Matrix<double, 2, 1> >, 4> * mesh,
                         vector<array<NodeWithId<Matrix<double, 2, 1> >, 2> > * boundary) {
  typedef Matrix<double, 2, 1> Point;
  typedef NodeWithId<Point> Node;

  const double a = sideLengths[0];
  const double b = sideLengths[1];
  const size_t n = numberOfElements[0];
  const size_t m = numberOfElements[1];

  double dx = a/double(n);
  double dy = b/double(m);

  size_t nextNodeId = 0;
  Eigen::Vector2d point;
  for (size_t j=0; j<=m; j++) {
    for (size_t i=0; i<=n; i++) {
      point[0] = double(i * dx);
      point[1] = double(j * dy);
      mesh->_nodes.push_back(Node(nextNodeId, point));
      ++nextNodeId;
    }
  }
  for (size_t j=0; j<m; j++) {
    for (size_t i=0; i<n; i++) {
      array<size_t, 4> connection;
      size_t bottomLeftCorner = j*(n+1)+i;
      connection[0] = bottomLeftCorner;
      connection[1] = bottomLeftCorner+1;
      connection[2] = bottomLeftCorner+(n+1)+1;
      connection[3] = bottomLeftCorner+(n+1);
      /*  cout << mesh->_nodes[bottomLeftCorner]._position << endl;
          cout << mesh->_nodes[bottomLeftCorner+1]._position << endl;
          cout << mesh->_nodes[bottomLeftCorner+n+2]._position << endl;
          cout << mesh->_nodes[bottomLeftCorner+n+1]._position << endl << endl; */
      mesh->_connectivity.push_back(connection);
    }
  }
  for (size_t i=0; i<n; i++) {
    array<Node, 2> boundaryConnection;
    boundaryConnection[0] = mesh->_nodes[i];
    boundaryConnection[1] = mesh->_nodes[i+1];
    boundary->push_back(boundaryConnection);
    boundaryConnection[0] = mesh->_nodes[(n+1)*(m+1)-i-1];
    boundaryConnection[1] = mesh->_nodes[(n+1)*(m+1)-(i+1)-1];
    boundary->push_back(boundaryConnection);
  }
}
#endif

#if 0
void
buildRectangularQuadMesh(const array<double, 2> & sideLengths,
                         const array<size_t, 2> & numberOfElements,
                         SingleElementMesh<NodeWithId<Matrix<double, 2, 1> >, 4> * mesh) {
  vector<array<NodeWithId<Matrix<double, 2, 1> >, 2> > boundary;
  buildRectangularQuadMesh(sideLengths, numberOfElements, mesh, &boundary);
}
#endif




#if 0
void
buildRectangularTriangleMesh(const array<double, 2> & sideLengths,
                             const array<size_t, 2> & numberOfElements,
                             SingleElementMesh<NodeWithId<Matrix<double, 2, 1> >, 3> * mesh,
                             vector<array<NodeWithId<Matrix<double, 2, 1> >, 2> > * boundary) {
  typedef Matrix<double, 2, 1> Point;
  typedef NodeWithId<Point> Node;

  const double a = sideLengths[0];
  const double b = sideLengths[1];
  const size_t n = numberOfElements[0];
  const size_t m = numberOfElements[1];

  double dx = a/double(n);
  double dy = b/double(m);

  size_t nextNodeId = 0;
  Eigen::Vector2d point;
  for (size_t j=0; j<=m; j++) {
    for (size_t i=0; i<=n; i++) {
      point[0] = double(i * dx);
      point[1] = double(j * dy);
      mesh->_nodes.push_back(Node(nextNodeId, point));
      ++nextNodeId;
    }
  }
  for (size_t j=0; j<m; j++) {
    for (size_t i=0; i<n; i++) {
      array<size_t, 3> connection;
      size_t bottomLeftCorner = j*(n+1)+i;
      // oriented one way
      /*
        connection[0] = bottomLeftCorner;
        connection[1] = bottomLeftCorner+1;
        connection[2] = bottomLeftCorner+(n+1)+1;
        mesh->_connectivity.push_back(connection);
        connection[0] = bottomLeftCorner+(n+1)+1;
        connection[1] = bottomLeftCorner+(n+1);
        connection[2] = bottomLeftCorner;
      */
      // oriented the other
      connection[0] = bottomLeftCorner;
      connection[1] = bottomLeftCorner+1;
      connection[2] = bottomLeftCorner+(n+1);
      mesh->_connectivity.push_back(connection);
      connection[0] = bottomLeftCorner+1;
      connection[1] = bottomLeftCorner+(n+1)+1;
      connection[2] = bottomLeftCorner+(n+1);
      mesh->_connectivity.push_back(connection);
    }
  }
  for (size_t i=0; i<n; i++) {
    array<Node, 2> boundaryConnection;
    boundaryConnection[0] = mesh->_nodes[i];
    boundaryConnection[1] = mesh->_nodes[i+1];
    boundary->push_back(boundaryConnection);
    boundaryConnection[0] = mesh->_nodes[(n+1)*(m+1)-i-1];
    boundaryConnection[1] = mesh->_nodes[(n+1)*(m+1)-(i+1)-1];
    boundary->push_back(boundaryConnection);
  }
}
#endif

#if 0
void
buildRectangularTriangleMesh(const array<double, 2> & sideLengths,
                             const array<size_t, 2> & numberOfElements,
                             SingleElementMesh<NodeWithId<Matrix<double, 2, 1> >, 3> * mesh) {
  vector<array<NodeWithId<Matrix<double, 2, 1> >, 2> > boundary;
  buildRectangularTriangleMesh(sideLengths, numberOfElements, mesh, &boundary);
}
#endif

template <class Element>
array<typename Element::Node, Element::NumberOfNodes>
getNodesFromMesh(const SingleElementMesh<Element> & mesh,
                 const size_t index) {
  array<typename Element::Node, Element::NumberOfNodes> elementNodes;
  for (size_t nodeIndex = 0; nodeIndex < Element::NumberOfNodes; ++nodeIndex) {
    elementNodes[nodeIndex] =
      mesh._nodes[mesh._connectivity[index][nodeIndex]];
  }
  return elementNodes;
}

bool
checkIfNodesAreEdgePeriodic(const Eigen::Matrix<double,3,1> & node1,
                            const Eigen::Matrix<double,3,1> & node2,
                            const size_t & coordinate,
                            const double & tolerance = 1.e-4) {

  return (std::abs(node1(coordinate) - node2(coordinate)) < tolerance);
}

bool
checkIfNodesAreFacePeriodic(const Eigen::Matrix<double,3,1> & node1,
                            const Eigen::Matrix<double,3,1> & node2,
                            const size_t & coordinate1,
                            const size_t & coordinate2,
                            const double & tolerance = 1.e-4) {

  return (std::abs(node1(coordinate1) - node2(coordinate1)) < tolerance &&
          std::abs(node1(coordinate2) - node2(coordinate2)) < tolerance);
}

template <class Mesh>
vector<size_t>
computeNodalAssignmentsForPeriodicBCs(Mesh & mesh,
                                      size_t * numberOfReducedNodes,
                                      const double & tolerance = 1.e-4,
                                      const bool & verbose = false){
  const unsigned int spatialDimension = 3;
  Matrix<double,spatialDimension,1> minVector, maxVector;
  findBoundingBoxOfGeneralMesh<Mesh,spatialDimension>(mesh,
                                                      &minVector,
                                                      &maxVector);


  const unsigned int numberOfNodes = mesh._nodes.size();
  vector<size_t> nodalAssignments(numberOfNodes);

  vector<size_t> xEdge, yEdge, zEdge;
  vector<size_t> x0Edge, y0Edge, z0Edge;
  vector<size_t> origin;
  vector<size_t> vertices;
  vector<size_t> x0Face, y0Face, z0Face;
  vector<size_t> xFace, yFace, zFace;
  vector<size_t> inner;

  for (unsigned int nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){

    Vector3d position = mesh._nodes[nodeIndex]._position;

    if (std::abs(position(0)-minVector(0)) < tolerance &&
        std::abs(position(1)-minVector(1)) < tolerance &&
        std::abs(position(2)-minVector(2)) < tolerance) {
      origin.push_back(nodeIndex);
    }
    else if (std::abs(position(0)-maxVector(0)) < tolerance &&
             ((std::abs(position(1)-minVector(1)) < tolerance ||
               std::abs(position(1)-maxVector(1)) < tolerance) &&
              (std::abs(position(2)-minVector(2)) < tolerance ||
               std::abs(position(2)-maxVector(2)) < tolerance)) &&
             std::abs(position(0)-minVector(0)) > tolerance) {
      vertices.push_back(nodeIndex);
    }
    else if (std::abs(position(1)-maxVector(1)) < tolerance &&
             ((std::abs(position(0)-minVector(0)) < tolerance ||
               std::abs(position(0)-maxVector(0)) < tolerance) &&
              (std::abs(position(2)-minVector(2)) < tolerance ||
               std::abs(position(2)-maxVector(2)) < tolerance)) &&
             std::abs(position(1)-minVector(1)) > tolerance) {
      vertices.push_back(nodeIndex);
    }
    else if (std::abs(position(2)-maxVector(2)) < tolerance &&
             ((std::abs(position(1)-minVector(1)) < tolerance ||
               std::abs(position(1)-maxVector(1)) < tolerance) &&
              (std::abs(position(0)-minVector(0)) < tolerance ||
               std::abs(position(0)-maxVector(0)) < tolerance)) &&
             std::abs(position(2)-minVector(2)) > tolerance) {
      vertices.push_back(nodeIndex);
    }
    else if (std::abs(position(0)-minVector(0)) < tolerance &&
             std::abs(position(1)-minVector(1)) < tolerance &&
             position(2) < maxVector(2) &&
             position(2) > minVector(2)) {
      z0Edge.push_back(nodeIndex);
    }
    else if (std::abs(position(0)-minVector(0)) < tolerance &&
             std::abs(position(2)-minVector(2)) < tolerance &&
             position(1) < maxVector(1) &&
             position(1) > minVector(1)) {
      y0Edge.push_back(nodeIndex);
    }
    else if (std::abs(position(1)-minVector(1)) < tolerance &&
             std::abs(position(2)-minVector(2)) < tolerance &&
             position(0) < maxVector(0) &&
             position(0) > minVector(0)) {
      x0Edge.push_back(nodeIndex);
    }
    else if ((std::abs(position(0)-minVector(0)) < tolerance &&      // edges in the z-direction
              std::abs(position(1)-maxVector(1)) < tolerance) ||
             (std::abs(position(0)-maxVector(0)) < tolerance &&
              std::abs(position(1)-maxVector(1)) < tolerance) ||
             (std::abs(position(0)-maxVector(0)) < tolerance &&
              std::abs(position(1)-minVector(1)) < tolerance)) {
      zEdge.push_back(nodeIndex);
    }
    else if ((std::abs(position(0)-minVector(0)) < tolerance &&      // edges in the y-direction
              std::abs(position(2)-maxVector(2)) < tolerance) ||
             (std::abs(position(0)-maxVector(0)) < tolerance &&
              std::abs(position(2)-maxVector(2)) < tolerance) ||
             (std::abs(position(0)-maxVector(0)) < tolerance &&
              std::abs(position(2)-minVector(2)) < tolerance)) {
      yEdge.push_back(nodeIndex);
    }
    else if ((std::abs(position(2)-minVector(2)) < tolerance &&      // edges in the x-direction
              std::abs(position(1)-maxVector(1)) < tolerance) ||
             (std::abs(position(2)-maxVector(2)) < tolerance &&
              std::abs(position(1)-maxVector(1)) < tolerance) ||
             (std::abs(position(2)-maxVector(2)) < tolerance &&
              std::abs(position(1)-minVector(1)) < tolerance)) {
      xEdge.push_back(nodeIndex);
    }
    else if (std::abs(position(0)-minVector(0)) < tolerance) {
      x0Face.push_back(nodeIndex);
    }
    else if (std::abs(position(0)-maxVector(0)) < tolerance) {
      xFace.push_back(nodeIndex);
    }
    else if (std::abs(position(1)-minVector(1)) < tolerance) {
      y0Face.push_back(nodeIndex);
    }
    else if (std::abs(position(1)-maxVector(1)) < tolerance) {
      yFace.push_back(nodeIndex);
    }
    else if (std::abs(position(2)-minVector(2)) < tolerance) {
      z0Face.push_back(nodeIndex);
    }
    else if (std::abs(position(2)-maxVector(2)) < tolerance) {
      zFace.push_back(nodeIndex);
    }
    else {
      inner.push_back(nodeIndex);
    }
  }

  if (verbose == true) {
    cout << "nodes on x-edges: " << xEdge.size() << endl;
    cout << "nodes on y-edges: " << yEdge.size() << endl;
    cout << "nodes on z-edges: " << zEdge.size() << endl;
    cout << "nodes on x0-edge: " << x0Edge.size() << endl;
    cout << "nodes on y0-edge: " << y0Edge.size() << endl;
    cout << "nodes on z0-edge: " << z0Edge.size() << endl;
    cout << "vertices: " << vertices.size() << endl;
    cout << "nodes on x0-face: " << x0Face.size() << endl;
    cout << "nodes on y0-face: " << y0Face.size() << endl;
    cout << "nodes on z0-face: " << z0Face.size() << endl;
    cout << "nodes on x-face: " << xFace.size() << endl;
    cout << "nodes on y-face: " << yFace.size() << endl;
    cout << "nodes on z-face: " << zFace.size() << endl;
    cout << "inner nodes: " << inner.size() << endl;
  }

  if (origin.size() != 1) {
    cout << "error: more than one origin detected" << endl;
    exit(1);
  }
  if (x0Face.size() != xFace.size()) {
    cout << "error: nodes on x-sides: " << x0Face.size() << " != " << xFace.size() << endl;
    exit(1);
  }
  if (y0Face.size() != yFace.size()) {
    cout << "error: nodes on y-sides: " << y0Face.size() << " != " << yFace.size() << endl;
    exit(1);
  }
  if (z0Face.size() != zFace.size()) {
    cout << "error: nodes on z-sides: " << z0Face.size() << " != " << zFace.size() << endl;
    exit(1);
  }
  if (vertices.size() != 7) {
    cout << "error: " << vertices.size() << " vertices deteceted" << endl;
    exit(1);
  }
  if (xEdge.size() != 3 * x0Edge.size()) {
    cout << "error: nodes on x-edges: " << 3 * x0Edge.size() << " != " << xEdge.size() << endl;
    exit(1);
  }
  if (yEdge.size() != 3 * y0Edge.size()) {
    cout << "error: nodes on y-edges: " << 3 * y0Edge.size() << " != " << yEdge.size() << endl;
    exit(1);
  }
  if (zEdge.size() != 3 * z0Edge.size()) {
    cout << "error: nodes on z-edges: " << 3 * z0Edge.size() << " != " << zEdge.size() << endl;
    exit(1);
  }
  if (inner.size() == 0) {
    cout << "error: no inner nodes detected" << endl;
    exit(1);
  }

  // set nodal assignments (i.e. for each node, assign the unique reduced node)
  vector<size_t> reducedNumbering(numberOfNodes);
  size_t nodeCount=0;
  reducedNumbering[origin[0]] = nodeCount;
  nodalAssignments[origin[0]] = nodeCount;
  for (size_t nodeID = 0; nodeID < vertices.size(); nodeID++) {
    reducedNumbering[vertices[nodeID]] = nodeCount;
    nodalAssignments[vertices[nodeID]] = nodeCount;
  }
  nodeCount++;
  for (size_t nodeID = 0; nodeID < x0Edge.size(); nodeID++) {
    reducedNumbering[x0Edge[nodeID]] = nodeCount;
    nodalAssignments[x0Edge[nodeID]] = nodeCount;
    nodeCount++;
  }
  for (size_t nodeID = 0; nodeID < y0Edge.size(); nodeID++) {
    reducedNumbering[y0Edge[nodeID]] = nodeCount;
    nodalAssignments[y0Edge[nodeID]] = nodeCount;
    nodeCount++;
  }
  for (size_t nodeID = 0; nodeID < z0Edge.size(); nodeID++) {
    reducedNumbering[z0Edge[nodeID]] = nodeCount;
    nodalAssignments[z0Edge[nodeID]] = nodeCount;
    nodeCount++;
  }
  for (size_t nodeID = 0; nodeID < x0Face.size(); nodeID++) {
    reducedNumbering[x0Face[nodeID]] = nodeCount;
    nodalAssignments[x0Face[nodeID]] = nodeCount;
    nodeCount++;
  }
  for (size_t nodeID = 0; nodeID < y0Face.size(); nodeID++) {
    reducedNumbering[y0Face[nodeID]] = nodeCount;
    nodalAssignments[y0Face[nodeID]] = nodeCount;
    nodeCount++;
  }
  for (size_t nodeID = 0; nodeID < z0Face.size(); nodeID++) {
    reducedNumbering[z0Face[nodeID]] = nodeCount;
    nodalAssignments[z0Face[nodeID]] = nodeCount;
    nodeCount++;
  }
  for (size_t nodeID = 0; nodeID < inner.size(); nodeID++) {
    reducedNumbering[inner[nodeID]] = nodeCount;
    nodalAssignments[inner[nodeID]] = nodeCount;
    nodeCount++;
  }
  *numberOfReducedNodes = nodeCount;

  for (size_t nodeID = 0; nodeID < xEdge.size(); nodeID++) {
    for (size_t nodeID2 = 0; nodeID2 < x0Edge.size(); nodeID2++) {
      if (checkIfNodesAreEdgePeriodic(mesh._nodes[xEdge[nodeID]]._position,mesh._nodes[x0Edge[nodeID2]]._position,0,tolerance)) {
        nodalAssignments[xEdge[nodeID]] = reducedNumbering[x0Edge[nodeID2]];
        break;
      }
    }
  }
  for (size_t nodeID = 0; nodeID < yEdge.size(); nodeID++) {
    for (size_t nodeID2 = 0; nodeID2 < y0Edge.size(); nodeID2++) {
      if (checkIfNodesAreEdgePeriodic(mesh._nodes[yEdge[nodeID]]._position,mesh._nodes[y0Edge[nodeID2]]._position,1,tolerance)) {
        nodalAssignments[yEdge[nodeID]] = reducedNumbering[y0Edge[nodeID2]];
        break;
      }
    }
  }
  for (size_t nodeID = 0; nodeID < zEdge.size(); nodeID++) {
    for (size_t nodeID2 = 0; nodeID2 < z0Edge.size(); nodeID2++) {
      if (checkIfNodesAreEdgePeriodic(mesh._nodes[zEdge[nodeID]]._position,mesh._nodes[z0Edge[nodeID2]]._position,2,tolerance)) {
        nodalAssignments[zEdge[nodeID]] = reducedNumbering[z0Edge[nodeID2]];
        break;
      }
    }
  }
  for (size_t nodeID = 0; nodeID < xFace.size(); nodeID++) {
    for (size_t nodeID2 = 0; nodeID2 < x0Face.size(); nodeID2++) {
      if (checkIfNodesAreFacePeriodic(mesh._nodes[xFace[nodeID]]._position,mesh._nodes[x0Face[nodeID2]]._position,1,2,tolerance)) {
        nodalAssignments[xFace[nodeID]] = reducedNumbering[x0Face[nodeID2]];
        break;
      }
    }
  }
  for (size_t nodeID = 0; nodeID < yFace.size(); nodeID++) {
    for (size_t nodeID2 = 0; nodeID2 < y0Face.size(); nodeID2++) {
      if (checkIfNodesAreFacePeriodic(mesh._nodes[yFace[nodeID]]._position,mesh._nodes[y0Face[nodeID2]]._position,0,2,tolerance)) {
        nodalAssignments[yFace[nodeID]] = reducedNumbering[y0Face[nodeID2]];
        break;
      }
    }
  }
  for (size_t nodeID = 0; nodeID < zFace.size(); nodeID++) {
    for (size_t nodeID2 = 0; nodeID2 < z0Face.size(); nodeID2++) {
      if (checkIfNodesAreFacePeriodic(mesh._nodes[zFace[nodeID]]._position,mesh._nodes[z0Face[nodeID2]]._position,0,1,tolerance)) {
        nodalAssignments[zFace[nodeID]] = reducedNumbering[z0Face[nodeID2]];
        break;
      }
    }
  }

  return nodalAssignments;
}

template <class Element>
vector<array<size_t, Element::SpatialDimension> >
computeBoundaryElementConnectivity(const SingleElementMesh<Element> & mesh) {
  static const unsigned int SpatialDimension = Element::SpatialDimension;
  static const unsigned int NumberOfNodes    = Element::NumberOfNodes;
  if (NumberOfNodes != SpatialDimension + 1) {
    throwException("Cannot do computeBoundaryElementConnectivity on anything "
                   "but simplicial elements");
  }

  vector<array<double, SpatialDimension> > vertices(mesh._nodes.size());
  for (unsigned int nodeIndex = 0; nodeIndex < mesh._nodes.size();
       ++nodeIndex) {
    for (unsigned int coordinate = 0; coordinate < SpatialDimension;
         ++coordinate) {
      vertices[nodeIndex][coordinate] =
        mesh._nodes[nodeIndex]._position(coordinate);
    }
  }

  // create the stlib mesh data structures
  geom::IndSimpSetIncAdj<SpatialDimension, SpatialDimension,
                         double> simplicialMesh;
  geom::build<SpatialDimension, SpatialDimension,
              double>(&simplicialMesh, vertices, mesh._connectivity);

  // run the stlib function to extract the boundary
  geom::IndSimpSetIncAdj<SpatialDimension, SpatialDimension - 1,
                         double> boundaryMesh;
  geom::buildBoundaryWithoutPacking<SpatialDimension, SpatialDimension,
                                    double>(simplicialMesh, &boundaryMesh);

  return boundaryMesh.indexedSimplices;
}
template<class Node, int Dimension>
struct BoundaryNodesWithNormal{
  vector<Node> _nodes;
  Eigen::Matrix<double,Dimension,1> _boundaryNormal;
};

template<class Mesh, class Node, class BoundaryNodesWithNormal>
void
fill2DHexagonMeshWithTriangularNodePlacement(Mesh * mesh,
                                             const int numberOfBars,
                                             const double barLength,
                                             vector<BoundaryNodesWithNormal> * vectorOfBoundaryNodesWithNormal){
  BoundaryNodesWithNormal tempBoundaryNodesWithNormal;
  size_t freeNodeIndex = 0;
  for (size_t sideIndex = 0; sideIndex < 12; sideIndex++){
    vectorOfBoundaryNodesWithNormal->push_back(tempBoundaryNodesWithNormal);
  }
  for (int rowIndex = -numberOfBars; rowIndex <= numberOfBars; rowIndex++){
    const int maxCols = (2 * numberOfBars) - std::abs(rowIndex);
    for (int colIndex = 0; colIndex <= maxCols; colIndex++){
      Node tempNode;
      tempNode._position << double(std::abs(rowIndex)) * 1/2. * barLength + double(colIndex) * barLength, double(rowIndex) * sqrt(3.)/2. * barLength;
      tempNode._id = freeNodeIndex;
      freeNodeIndex++;
      mesh->_nodes.push_back(tempNode);
      if (rowIndex == -numberOfBars && colIndex !=0 && colIndex != maxCols){
        (*vectorOfBoundaryNodesWithNormal)[0]._nodes.push_back(tempNode);
      }
      else if (rowIndex == -numberOfBars && colIndex == 0){
        (*vectorOfBoundaryNodesWithNormal)[1]._nodes.push_back(tempNode);
      }
      else if (rowIndex < 0 && rowIndex != -numberOfBars && colIndex == 0){
        (*vectorOfBoundaryNodesWithNormal)[2]._nodes.push_back(tempNode);
      }
      else if (rowIndex == 0 && colIndex == 0){
        (*vectorOfBoundaryNodesWithNormal)[3]._nodes.push_back(tempNode);
      }
      else if (rowIndex > 0 && rowIndex != numberOfBars && colIndex == 0){
        (*vectorOfBoundaryNodesWithNormal)[4]._nodes.push_back(tempNode);
      }
      else if (rowIndex == numberOfBars && colIndex == 0){
        (*vectorOfBoundaryNodesWithNormal)[5]._nodes.push_back(tempNode);
      }
      else if (rowIndex == numberOfBars && colIndex != 0 && colIndex != maxCols){
        (*vectorOfBoundaryNodesWithNormal)[6]._nodes.push_back(tempNode);
      }
      else if (rowIndex == numberOfBars && colIndex == maxCols){
        (*vectorOfBoundaryNodesWithNormal)[7]._nodes.push_back(tempNode);
      }
      else if (rowIndex > 0 && rowIndex != numberOfBars && colIndex == maxCols){
        (*vectorOfBoundaryNodesWithNormal)[8]._nodes.push_back(tempNode);
      }
      else if (rowIndex == 0 && colIndex == maxCols){
        (*vectorOfBoundaryNodesWithNormal)[9]._nodes.push_back(tempNode);
      }
      else if (rowIndex < 0 && rowIndex != -numberOfBars && colIndex == maxCols){
        (*vectorOfBoundaryNodesWithNormal)[10]._nodes.push_back(tempNode);
      }
      else if (rowIndex == -numberOfBars && colIndex == maxCols){
        (*vectorOfBoundaryNodesWithNormal)[11]._nodes.push_back(tempNode);
      }
    }
  }
  double theta = -M_PI/2.;
  for (size_t sideIndex = 0; sideIndex < 12; sideIndex++){
    (*vectorOfBoundaryNodesWithNormal)[sideIndex]._boundaryNormal << cos(theta) ,sin(theta);
    theta -= M_PI/6.;
  }
}

void
initializeBeamsOfHexagonMesh(const size_t & numberOfBars,
                             array<vector<array<size_t,2>>,3> * vectorsOfNodeSetsForEachBeamType){
  (*vectorsOfNodeSetsForEachBeamType)[0].push_back(array<size_t,2>{{0,1}});
  (*vectorsOfNodeSetsForEachBeamType)[1].push_back(array<size_t,2>{{1,numberOfBars + 2}});
  (*vectorsOfNodeSetsForEachBeamType)[2].push_back(array<size_t,2>{{numberOfBars + 2,0}});
}

template<class Node, class Mesh>
void
fractileTypeFillEqualateralTriangularMeshWithThreeTypesOfBeams(const double barLength,
                                                               const int numberOfBars,
                                                               Mesh * mesh,
                                                               array<vector<array<size_t,2>>,3> * vectorsOfNodeSetsForEachBeamType){
  vector<Node> masterVectorsTipNode;
  vector<Node> masterVectorsTailNode;
  for (size_t beamTypeIndex = 0; beamTypeIndex < 3; beamTypeIndex++){
    // provide initial vector for fractile setup
    const size_t vectorTipIndex = (*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex][0][0];
    const size_t vectorTailIndex = (*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex][0][1];
    masterVectorsTipNode.push_back(mesh->_nodes[vectorTipIndex]);
    masterVectorsTailNode.push_back(mesh->_nodes[vectorTailIndex]);
    mesh->_connectivity.push_back(array<size_t,2>{{vectorTailIndex,vectorTipIndex}});
    size_t iterations = 0;
    const size_t maximumIterations = 6 * numberOfBars;
    while (masterVectorsTipNode.size() > 0 && iterations < maximumIterations){
      iterations++;
      size_t numberOfMasterVectorsAddedThisWhileIteration = 0;
      for (size_t masterVectorsIndex = 0; masterVectorsIndex < masterVectorsTipNode.size() ; masterVectorsIndex++){
        for (size_t pointIndex = 0; pointIndex < mesh->_nodes.size(); pointIndex++){
          Vector2d comparisonVector = mesh->_nodes[pointIndex]._position - masterVectorsTipNode[masterVectorsIndex]._position;
          Vector2d masterVectorToCompare = masterVectorsTipNode[masterVectorsIndex]._position -
            masterVectorsTailNode[masterVectorsIndex]._position;
          const double dotProductResult = masterVectorToCompare.dot(comparisonVector)/(barLength * barLength);
          if (dotProductResult < (1+1e-2) * cos(M_PI/3.) &&
              dotProductResult > (1-1e-2) * cos(M_PI/3.) &&
              comparisonVector.norm() < (1+1e-2) * barLength){
            string unique = "true";
            for (size_t uniqueIndex = 0; uniqueIndex < (*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex].size(); uniqueIndex++){
              if (std::min(masterVectorsTipNode[masterVectorsIndex]._id, pointIndex) ==
                  std::min((*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex][uniqueIndex][0],
                           (*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex][uniqueIndex][1])){
                if (std::max(masterVectorsTipNode[masterVectorsIndex]._id,pointIndex) ==
                    std::max((*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex][uniqueIndex][0],
                             (*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex][uniqueIndex][1])){
                  unique = "false";

                }
              }
            }
            if (unique.compare("true") == 0){
              mesh->_connectivity.push_back(array<size_t,2>{{masterVectorsTipNode[masterVectorsIndex]._id,pointIndex}});
              (*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex].push_back({{(masterVectorsTipNode[masterVectorsIndex]._id),pointIndex}});
              numberOfMasterVectorsAddedThisWhileIteration++;
            }
          }
        }
      }
      masterVectorsTipNode.clear();
      masterVectorsTailNode.clear();
      for(size_t addedVectorIndex = 1; addedVectorIndex <= numberOfMasterVectorsAddedThisWhileIteration; addedVectorIndex++){
        const size_t lengthOfBeamNodeSet = (*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex].size();
        Node tempTipNode;
        tempTipNode._id = (*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex][lengthOfBeamNodeSet - addedVectorIndex][1];
        tempTipNode._position = mesh->_nodes[tempTipNode._id]._position;
        masterVectorsTipNode.push_back(tempTipNode);
        Node tempTailNode;
        tempTailNode._id = (*vectorsOfNodeSetsForEachBeamType)[beamTypeIndex][lengthOfBeamNodeSet - addedVectorIndex][0];
        tempTailNode._position = mesh->_nodes[tempTailNode._id]._position;
        masterVectorsTailNode.push_back(tempTailNode);
      }
      if (iterations == maximumIterations - 1){
        errorStatement("reached iteration limit on beamTypeIndex %lu\n",beamTypeIndex);
      }
    }
  }

}

template<class Mesh, class Node, class BoundaryNodesWithNormal>
void
buildHexagonMeshWithThreeTypesOfBarsInHoneycombPattern(Mesh * mesh,
                                                       const int numberOfBars,
                                                       const double barLength,
                                                       vector<BoundaryNodesWithNormal> * vectorOfBoundaryNodesWithNormal,
                                                       array<vector<array<size_t,2>>,3> * vectorsOfNodeSetsForEachBeamType){

  fill2DHexagonMeshWithTriangularNodePlacement<Mesh,Node,BoundaryNodesWithNormal>(mesh,
                                                                                  numberOfBars,
                                                                                  barLength,
                                                                                  vectorOfBoundaryNodesWithNormal);

  initializeBeamsOfHexagonMesh(numberOfBars,
                               vectorsOfNodeSetsForEachBeamType);


  fractileTypeFillEqualateralTriangularMeshWithThreeTypesOfBeams<Node,Mesh>(barLength,
                                                                            numberOfBars,
                                                                            mesh,
                                                                            vectorsOfNodeSetsForEachBeamType);

}

template <class Element>
vector<array<size_t, 2> >
computeAdjacentPairsOfSimplicialMesh(const SingleElementMesh<Element> & mesh) {
  static const unsigned int SpatialDimension = Element::SpatialDimension;
  static const unsigned int NumberOfNodes    = Element::NumberOfNodes;
  if (NumberOfNodes != SpatialDimension + 1) {
    throwException("Cannot do computeBoundaryElementConnectivity on anything "
                   "but simplicial elements");
  }

  vector<array<double, SpatialDimension> > vertices(mesh._nodes.size());
  for (unsigned int nodeIndex = 0; nodeIndex < mesh._nodes.size();
       ++nodeIndex) {
    for (unsigned int coordinate = 0; coordinate < SpatialDimension;
         ++coordinate) {
      vertices[nodeIndex][coordinate] =
        mesh._nodes[nodeIndex]._position(coordinate);
    }
  }

  // create the stlib mesh data structures
  geom::IndSimpSetIncAdj<SpatialDimension, SpatialDimension,
                         double> simplicialMesh;
  geom::build<SpatialDimension, SpatialDimension,
              double>(&simplicialMesh, vertices, mesh._connectivity);
  const unsigned int numberOfElements = mesh._connectivity.size();

  Eigen::SparseMatrix<int> adjacentPairs(numberOfElements, numberOfElements);

  for (unsigned int elementIndex = 0; elementIndex < numberOfElements;
       ++elementIndex) {
    for (unsigned int nodeNumber = 0; nodeNumber < NumberOfNodes;
         ++nodeNumber) {
      const size_t adjacentElementIndex =
        simplicialMesh.adjacent[elementIndex][nodeNumber];
      if (adjacentElementIndex <= numberOfElements) {
        array<size_t, 2> adjacentPair;
        adjacentPair[0] = elementIndex;
        adjacentPair[1] = adjacentElementIndex;
        if (adjacentPair[1] < adjacentPair[0]) {
          const size_t temp = adjacentPair[0];
          adjacentPair[0] = adjacentPair[1];
          adjacentPair[1] = temp;
        }
        adjacentPairs.coeffRef(adjacentPair[0], adjacentPair[1]) += 1;
      }
    }
  }

  vector<array<size_t, 2> > returnValue;
  returnValue.reserve(1e4);
  // for each nonzero matrix entry
  for (unsigned int outerIndex=0;
       outerIndex < adjacentPairs.outerSize(); ++outerIndex) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(adjacentPairs,
                                                    outerIndex);
         it; ++it) {

      // what's this entry's row, column, and value?
      const unsigned int row = it.row();
      const unsigned int col = it.col();
      const unsigned int matrixEntry = it.value();
      if (matrixEntry != 2) {
        throwException("something went wrong in computeAdjacentPairsOfSimplicialMesh "
                       "because the matrix entry for pair %u %u is not 2: %u",
                       row, col, matrixEntry);
      }
      returnValue.push_back((array<size_t, 2>) {{row, col}});
    }
  }
  return returnValue;
}

template <class Mesh, class Element1, class Element2>
void
printTwoElementMeshToFile (const Mesh & mesh, const string & filename) {
  static const size_t  SpatialDimensions = Element1::SpatialDimension;
  static const size_t     NumberOfNodes1 = Element1::NumberOfNodes;
  static const size_t     NumberOfNodes2 = Element2::NumberOfNodes;
  
  std::ofstream file;
  file.open(filename.c_str());
  for (unsigned int nodeIndex = 0; nodeIndex < mesh._nodes.size() ; nodeIndex++) {
    file << mesh._nodes[nodeIndex]._position(0);
    for (unsigned int dimensionIndex = 1; dimensionIndex < SpatialDimensions; dimensionIndex++) {
      file << ", " << mesh._nodes[nodeIndex]._position(dimensionIndex);
    }
    file << endl;
  }
  file << "connectivity1" << endl;
  for (unsigned int elementIndex = 0; elementIndex < mesh._connectivity1.size() ; elementIndex++) {
    array<size_t, NumberOfNodes1> nodeIds = mesh._connectivity1[elementIndex];
    file << nodeIds[0];
    for (unsigned int nodeIndex = 1; nodeIndex < NumberOfNodes1; nodeIndex++) {
      file << ", " << nodeIds[nodeIndex];
    }
    file << endl;
  }
  file << "connectivity2" << endl;
  for (unsigned int elementIndex = 0; elementIndex < mesh._connectivity2.size() ; elementIndex++) {
    array<size_t, NumberOfNodes2> nodeIds = mesh._connectivity2[elementIndex];
    file << nodeIds[0];
    for (unsigned int nodeIndex = 1; nodeIndex < NumberOfNodes2; nodeIndex++) {
      file << ", " << nodeIds[nodeIndex];
    }
    file << endl;
  }
  file.close();
}

}

#endif//MESH_UTILITIES_H

#include <algorithm>
#include <limits>
#include "../lib/VectorMath.h"
#include "TriangTree.h"

TriangTree::TriangTree(double *triangs, double *aabb, int numTriangs,
					   double *rootBox, int maxPrimsInLeaf)
{
	// allocate
	P0x = new double[numTriangs];
	P0y = new double[numTriangs];
	P0z = new double[numTriangs];
	P1x = new double[numTriangs];
	P1y = new double[numTriangs];
	P1z = new double[numTriangs];
	P2x = new double[numTriangs];
	P2y = new double[numTriangs];
	P2z = new double[numTriangs];
	xMin = new double[numTriangs];
	yMin = new double[numTriangs];
	zMin = new double[numTriangs];
	xMax = new double[numTriangs];
	yMax = new double[numTriangs];
	zMax = new double[numTriangs];
	E1x = new double[numTriangs];
	E1y = new double[numTriangs];
	E1z = new double[numTriangs];
	E2x = new double[numTriangs];
	E2y = new double[numTriangs];
	E2z = new double[numTriangs];

	// copy triangs and aabb
	std::copy(triangs + 0 * numTriangs, triangs + 1 * numTriangs, P0x);
	std::copy(triangs + 1 * numTriangs, triangs + 2 * numTriangs, P0y);
	std::copy(triangs + 2 * numTriangs, triangs + 3 * numTriangs, P0z);
	std::copy(triangs + 3 * numTriangs, triangs + 4 * numTriangs, P1x);
	std::copy(triangs + 4 * numTriangs, triangs + 5 * numTriangs, P1y);
	std::copy(triangs + 5 * numTriangs, triangs + 6 * numTriangs, P1z);
	std::copy(triangs + 6 * numTriangs, triangs + 7 * numTriangs, P2x);
	std::copy(triangs + 7 * numTriangs, triangs + 8 * numTriangs, P2y);
	std::copy(triangs + 8 * numTriangs, triangs + 9 * numTriangs, P2z);
	std::copy(aabb + 0 * numTriangs, aabb + 1 * numTriangs, xMin);
	std::copy(aabb + 1 * numTriangs, aabb + 2 * numTriangs, yMin);
	std::copy(aabb + 2 * numTriangs, aabb + 3 * numTriangs, zMin);
	std::copy(aabb + 3 * numTriangs, aabb + 4 * numTriangs, xMax);
	std::copy(aabb + 4 * numTriangs, aabb + 5 * numTriangs, yMax);
	std::copy(aabb + 5 * numTriangs, aabb + 6 * numTriangs, zMax);

	// compute E1, E2
	VectorMath::Sub(P0x, P1x, E1x, numTriangs);
	VectorMath::Sub(P0y, P1y, E1y, numTriangs);
	VectorMath::Sub(P0z, P1z, E1z, numTriangs);
	VectorMath::Sub(P0x, P2x, E2x, numTriangs);
	VectorMath::Sub(P0y, P2y, E2y, numTriangs);
	VectorMath::Sub(P0z, P2z, E2z, numTriangs);

	// copy root box
	rootNode = new Node3D;
	for (int i = 0; i < 6; i++)
		rootNode->box[i] = rootBox[i];

	// fill triangs numbers
	int *primInds = new int[numTriangs];
	for (int i = 0; i < numTriangs; i++)
		primInds[i] = i;
	rootNode->primInds = primInds;
	rootNode->numPrims = numTriangs;

	// run tree building
	int curDepth = 0;
	buildRecur(rootNode, maxPrimsInLeaf, curDepth);
}

// recursively builds nodes
void TriangTree::buildRecur(Node3D *parNode, int maxNumPrimsInLeaf, int curDepth)
{
	int *primInds = parNode->primInds;
	int numPrims = parNode->numPrims;
	double *box = parNode->box;

	// preset Node3D fields
	parNode->leftNode = NULL;
	parNode->rightNode = NULL;

	// max depth, finalize
	if (curDepth == MAX_DEPTH)
		return;

	// small number of primitives in leaf, finalize
	if (numPrims <= maxNumPrimsInLeaf)
		return;

	// choose the longest axis as the axis for splitting
	int curAxis = 0;
	double maxLen = box[3] - box[0];
	for (int i = 1; i < 3; i++)
		if (box[3+i] - box[i] > maxLen)
		{
			maxLen = box[3+i] - box[i];
			curAxis = i;
		}

		// choose axis
		double *cMin, *cMax;
		switch (curAxis)
		{
		case 0:
			cMin = xMin;
			cMax = xMax;
			break;
		case 1:
			cMin = yMin;
			cMax = yMax;
			break;
		case 2:
			cMin = zMin;
			cMax = zMax;
			break;
		default:
			cMin = xMin;
			cMax = xMax;
			break;
		}

		// find split and make child primInds
		int lNum, rNum;
		double splitCoord = getOptimalSplit(parNode, cMin, cMax, &lNum, &rNum);
		int* leftPrimInds = new int[lNum];
		int* rightPrimInds = new int[rNum];
		int lInd = 0;
		int rInd = 0;
		for (int i = 0; i < numPrims; i++)
		{
			if (cMin[primInds[i]] <= splitCoord)
			{
				leftPrimInds[lInd] = primInds[i];
				lInd++;
			}
			if (cMax[primInds[i]] >= splitCoord)
			{
				rightPrimInds[rInd] = primInds[i];
				rInd++;
			}
		}

		// make left node
		Node3D *leftNode = new Node3D;
		for (int i = 0; i < 6; i++)
			leftNode->box[i] = parNode->box[i];
		leftNode->box[3 + curAxis] = splitCoord;
		leftNode->primInds = leftPrimInds;
		leftNode->numPrims = lNum;

		// make right node
		Node3D *rightNode = new Node3D;
		for (int i = 0; i < 6; i++)
			rightNode->box[i] = parNode->box[i];
		rightNode->box[0 + curAxis] = splitCoord;
		rightNode->primInds = rightPrimInds;
		rightNode->numPrims = rNum;

		// update parent node
		parNode->leftNode = leftNode;
		parNode->rightNode = rightNode;
		parNode->splitAxis = curAxis;
		parNode->splitCoord = splitCoord;

		// run recursive build
		buildRecur(leftNode, maxNumPrimsInLeaf, curDepth + 1);
		buildRecur(rightNode, maxNumPrimsInLeaf, curDepth + 1);
}

// computes good splitCoordinate
double TriangTree::getOptimalSplit(Node3D *node, double *cMin, double *cMax,
								   int *lNum, int *rNum)
{
	const int MAX_BIN_SEARCH = 10;
	int* primInds = node->primInds;
	int numPrims = node->numPrims;
	double cMinMin, cMaxMax, curSplit, bestSplit;
	double curMeas, bestMeas;

	// find minimum and maximum
	cMinMin = std::numeric_limits<double>::infinity();
	cMaxMax = -std::numeric_limits<double>::infinity();
	for (int i = 0; i < numPrims; i++)
	{
		if (cMin[primInds[i]] < cMinMin)
			cMinMin = cMin[primInds[i]];
		if (cMax[primInds[i]] > cMaxMax)
			cMaxMax = cMax[primInds[i]];
	}

	// 0 iteration
	curSplit = (cMinMin + cMaxMax) / 2.0;
	bestSplit = curSplit;
	countPrimsInChilds(node, cMin, cMax, curSplit, lNum, rNum);
	bestMeas = abs((*rNum) * (curSplit - cMinMin) - (*lNum) * (cMaxMax - curSplit));

	// get split coordinate - binary search
	int curBinDepth = 0;
	while ((*rNum) - (*lNum) != 0 && curBinDepth != MAX_BIN_SEARCH)
	{
		// find new split and compute new numbers of primitives
		if ((*rNum) > (*lNum))
			cMinMin = curSplit;
		else
			cMaxMax = curSplit;
		curSplit = (cMinMin + cMaxMax) / 2.0;
		countPrimsInChilds(node, cMin, cMax, curSplit, lNum, rNum);
		curMeas = abs((*rNum) * (curSplit - cMinMin) - (*lNum) * (cMaxMax - curSplit));

		// update bestSplit if necessary
		if (curMeas < bestMeas)
			bestSplit = curSplit;

		curBinDepth++;
	}

	// update lNum, rNum and exit
	countPrimsInChilds(node, cMin, cMax, bestSplit, lNum, rNum);
	return bestSplit;
}

// computes number of primitives in child nodes provided by splitCoord
void TriangTree::countPrimsInChilds(Node3D *node, double *cMin, double *cMax,
									double splitCoord, int *lNum, int *rNum)
{
	int* primInds = node->primInds;
	int numPrims = node->numPrims;

	(*lNum) = 0;
	(*rNum) = 0;
	for (int i = 0; i < numPrims; i++)
	{
		if (cMin[*primInds] <= splitCoord)
			(*lNum)++;
		if (cMax[*primInds] >= splitCoord)
			(*rNum)++;
		primInds++;
	}
}

// computes intersections
void TriangTree::getInter(double *ray, int numRays, double *t, double *ind)
{
	// inf const
	double const Inf = std::numeric_limits<double>::infinity();

	// allocate and declare
	double curTMin;
	double curTMax;
	double *curRay = new double[7];

	// process ray-by-ray
	for (int i = 0; i < numRays; i++)
	{
		curTMin = -Inf;
		curTMax = Inf;
		curRay[0] = ray[i + 0 * numRays];
		curRay[1] = ray[i + 1 * numRays];
		curRay[2] = ray[i + 2 * numRays];
		curRay[3] = ray[i + 3 * numRays];
		curRay[4] = ray[i + 4 * numRays];
		curRay[5] = ray[i + 5 * numRays];
		curRay[6] = ray[i + 6 * numRays];
		t[i] = Inf;
		// if tree can be intersected
		if (isRootBoxIntersected(curRay, &curTMin, &curTMax))
		{
			t[i] = getInterRecur(curRay, rootNode, curTMin, curTMax, ind + i);
		}
	}
	delete [] curRay;
}

// checks intersection of ray and box
bool TriangTree::isRootBoxIntersected(double *ray, double *tMin, double *tMax)
{
	// An Efficient and Robust Ray-Box Intersection Algorithm
	// Amy Williams, Steve Barrus, ...    
	double tyMin, tyMax, tzMin, tzMax;
	double *rootBox = rootNode->box;

	if (ray[3] >= 0)
	{
		(*tMin) = (rootBox[0] - ray[0]) / ray[3];
		(*tMax) = (rootBox[3] - ray[0]) / ray[3];
	}
	else
	{
		(*tMin) = (rootBox[3] - ray[0]) / ray[3];
		(*tMax) = (rootBox[0] - ray[0]) / ray[3];
	}
	if (ray[4] >= 0)
	{
		tyMin = (rootBox[1] - ray[1]) / ray[4];
		tyMax = (rootBox[4] - ray[1]) / ray[4];
	}
	else
	{
		tyMin = (rootBox[4] - ray[1]) / ray[4];
		tyMax = (rootBox[1] - ray[1]) / ray[4];
	}
	if (((*tMin) > tyMax) || (tyMin > (*tMax)))
		return false;
	if (tyMin > (*tMin))
		(*tMin) = tyMin;
	if (tyMax < (*tMax))
		(*tMax) = tyMax;
	if (ray[5] >= 0)
	{
		tzMin = (rootBox[2] - ray[2]) / ray[5];
		tzMax = (rootBox[5] - ray[2]) / ray[5];
	}
	else
	{
		tzMin = (rootBox[5] - ray[2]) / ray[5];
		tzMax = (rootBox[2] - ray[2]) / ray[5];
	}
	if (((*tMin) > tzMax) || (tzMin > (*tMax)))
		return false;
	if (tzMin > (*tMin))
		(*tMin) = tzMin;
	if (tzMax < (*tMax))
		(*tMax) = tzMax;
	return true;
}

// recursive search of intersection for one ray
double TriangTree::getInterRecur(double *ray, Node3D *node, double tMin, double tMax, double *ind)
{
	if (node->leftNode == NULL || node->rightNode == NULL)
	{
		// leaf
		int numPrims = node->numPrims;
		int primInd;
		double Px, Py, Pz;
		double deter, invDeter;
		double Tx, Ty, Tz;
		double Qx, Qy, Qz;
		double u, v;
		double tTemp;
		tMin = std::numeric_limits<double>::infinity();

		// cycle on primitives
		for (int i = 0; i < numPrims; i++)
		{
			primInd = node->primInds[i];

			// P = cross(D, E2)
			Px = ray[4] * E2z[primInd] - ray[5] * E2y[primInd];
			Py = ray[5] * E2x[primInd] - ray[3] * E2z[primInd];
			Pz = ray[3] * E2y[primInd] - ray[4] * E2x[primInd];
			// deter = dot(E1, P)
			deter = E1x[primInd] * Px + E1y[primInd] * Py + E1z[primInd] * Pz;
			// check 1
			if (deter > -1e-15 && deter < 1e-15)
				continue;

			invDeter = 1.0 / deter;
			// T = O - P0
			Tx = ray[0] - P0x[primInd];
			Ty = ray[1] - P0y[primInd];
			Tz = ray[2] - P0z[primInd];
			// u = dot(T, P) * invDet
			u = (Tx * Px + Ty * Py + Tz * Pz) * invDeter;
			// check 2
			if (u < 0.0 || u > 1.0)
				continue;

			// Q = cross(T, E1)
			Qx = Ty * E1z[primInd] - Tz * E1y[primInd];
			Qy = Tz * E1x[primInd] - Tx * E1z[primInd];
			Qz = Tx * E1y[primInd] - Ty * E1x[primInd];
			// v = dot(D, Q) * invDet
			v = (ray[3] * Qx + ray[4] * Qy + ray[5] * Qz) * invDeter;
			// check 3
			if (v < 0.0 || u + v > 1.0)
				continue;

			// fin t = dot(E2, Q) * invDet
			tTemp = (E2x[primInd] * Qx + E2y[primInd] * Qy + E2z[primInd] * Qz) * invDeter;
			if (tTemp > 1e-6 && tTemp < tMin)
			{
				tMin = tTemp;
				(*ind) = primInd + 1;
			}
		}

		return tMin;
	}
	else
	{
		// find split t
		int splitAxis = node->splitAxis;
		double splitCoord = node->splitCoord;
		double tSplit = (splitCoord - ray[splitAxis]) / ray[splitAxis + 3];

		// get near and far nodes
		Node3D *nearNode;
		Node3D *farNode;
		if (ray[splitAxis + 3] >= 0)
		{
			nearNode = node->leftNode;
			farNode = node->rightNode;
		}
		else
		{
			nearNode = node->rightNode;
			farNode = node->leftNode;
		}

		// run recursive intersection search
		if (tSplit <= tMin)
			return getInterRecur(ray, farNode, tMin, tMax, ind);
		else
			if (tSplit >= tMax)
				return getInterRecur(ray, nearNode, tMin, tMax, ind);
			else
			{
				double ind1, ind2;
				double t1, t2;
				t1 = getInterRecur(ray, nearNode, tMin, tSplit, &ind1);
				t2 = getInterRecur(ray, farNode, tSplit, tMax, &ind2);
				if (t1 < t2)
				{
					(*ind) = ind1;
					return t1;
				}
				else
				{
					(*ind) = ind2;
					return t2;
				}
			}
	}
}

TriangTree::~TriangTree()
{
	delete [] P0x;
	delete [] P0y;
	delete [] P0z;
	delete [] P1x;
	delete [] P1y;
	delete [] P1z;
	delete [] P2x;
	delete [] P2y;
	delete [] P2z;
	delete [] xMin;
	delete [] yMin;
	delete [] zMin;
	delete [] xMax;
	delete [] yMax;
	delete [] zMax;
	delete [] E1x;
	delete [] E1y;
	delete [] E1z;
	delete [] E2x;
	delete [] E2y;
	delete [] E2z;
	deleteNode(rootNode);
}

// deletes node
void TriangTree::deleteNode(Node3D *node)
{
	if (node->leftNode != NULL)
		deleteNode(node->leftNode);
	if (node->rightNode != NULL)
		deleteNode(node->rightNode);
	delete [] node->primInds;
	delete node;
}

int TriangTree::numLeafs()
{
	return countLeafs(rootNode);
}

int TriangTree::countLeafs(Node3D *node)
{
	if (node->leftNode == NULL || node->rightNode == NULL)
		return 1;
	else
	{
		return countLeafs(node->leftNode) + countLeafs(node->rightNode);
	}
}

double TriangTree::averagePrimNumInLeaf()
{
	return computeAveragePrimNumInLeaf(rootNode);
}

double TriangTree::computeAveragePrimNumInLeaf(Node3D *node)
{
	if (node->leftNode == NULL || node->rightNode == NULL)
		return node->numPrims;
	else
	{
		return (computeAveragePrimNumInLeaf(node->leftNode) +
			computeAveragePrimNumInLeaf(node->rightNode)) / 2.0;
	}
}

int TriangTree::minPrimNumInLeaf()
{
	return computeMinPrimNumInLeaf(rootNode);
}

int TriangTree::computeMinPrimNumInLeaf(Node3D *node)
{
	if (node->leftNode == NULL || node->rightNode == NULL)
		return node->numPrims;
	else
	{
		return std::min(computeMinPrimNumInLeaf(node->leftNode),
			computeMinPrimNumInLeaf(node->rightNode));
	}
}

int TriangTree::maxPrimNumInLeaf()
{
	return computeMaxPrimNumInLeaf(rootNode);
}

int TriangTree::computeMaxPrimNumInLeaf(Node3D *node)
{
	if (node->leftNode == NULL || node->rightNode == NULL)
		return node->numPrims;
	else
	{
		return std::max(computeMaxPrimNumInLeaf(node->leftNode),
			computeMaxPrimNumInLeaf(node->rightNode));
	}
}

int TriangTree::minDepth()
{
	return computeMinDepth(rootNode);
}

int TriangTree::computeMinDepth(Node3D *node)
{
	if (node->leftNode == NULL || node->rightNode == NULL)
		return 0;
	else
	{
		return std::min(computeMinDepth(node->leftNode),
			computeMinDepth(node->rightNode)) + 1;
	}
}

int TriangTree::maxDepth()
{
	return computeMaxDepth(rootNode);
}

int TriangTree::computeMaxDepth(Node3D *node)
{
	if (node->leftNode == NULL || node->rightNode == NULL)
		return 0;
	else
	{
		return std::max(computeMaxDepth(node->leftNode),
			computeMaxDepth(node->rightNode)) + 1;
	}
}
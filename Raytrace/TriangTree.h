#ifndef _TRIANGTREE_H
#define _TRIANGTREE_H

#include "mex.h"

struct Node3D
{
	int *primInds;
	int numPrims;
	double box[6];
	int splitAxis;
	double splitCoord;
	Node3D *leftNode, *rightNode;
};

class TriangTree
{
private:

	double *P0x, *P0y, *P0z;
	double *P1x, *P1y, *P1z;
	double *P2x, *P2y, *P2z;
	double *E1x, *E1y, *E1z;
	double *E2x, *E2y, *E2z;
	double *xMin, *xMax;
	double *yMin, *yMax;
	double *zMin, *zMax;
	Node3D *rootNode;
	static const int MAX_DEPTH = 15;

public:

	// constructor
	TriangTree(double *triangs, double *aabb, int numTriangs,
		double *rootBox, int maxPrimsInLeaf);

	// computes intersections
	void getInter(double *ray, int numRays, double *t, double *ind);

	// information funcs
	int numLeafs();
	double averagePrimNumInLeaf();
	int minPrimNumInLeaf();
	int maxPrimNumInLeaf();
	int minDepth();
	int maxDepth();

	// destructor
	~TriangTree();

private:

	// recursively builds nodes
	void buildRecur(Node3D *parNode, int maxNumPrimsInLeaf, int curDepth);

	// computes good splitCoordinate
	double getOptimalSplit(Node3D *node, double *cMin, double *cMax,
		int *lNum, int *rNum);

	// computes number of primitives in child nodes provided by splitCoord
	void countPrimsInChilds(Node3D *node, double *cMin, double *cMax,
		double splitCoord, int *lNum, int *rNum);

	// checks intersection of ray and box
	bool isRootBoxIntersected(double *ray, double *tMin, double *tMax);

	// recursive search of intersection for one ray
	double getInterRecur(double *ray, Node3D *node, double tMin, double tMax, double *ind);

	// utilites for information funcs
	int countLeafs(Node3D *node);
	double computeAveragePrimNumInLeaf(Node3D *node);
	int computeMinPrimNumInLeaf(Node3D *node);
	int computeMaxPrimNumInLeaf(Node3D *node);
	int computeMinDepth(Node3D *node);
	int computeMaxDepth(Node3D *node);

	// deletes node
	void deleteNode(Node3D *node);
};

#endif /* _TRIANGTREE_H */

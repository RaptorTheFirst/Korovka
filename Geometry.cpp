#include <fstream>
#include <vector>
#include <cstdlib>
#include "Task.cpp"
using namespace std;

struct POINT {
	double x, y, z;
};

struct TRIANGLE {
	int p1, p2, p3, n1, n2, n3;
};

POINT crossLine(POINT & A, POINT & B, POINT & vec, POINT & pntOfView) {
	POINT crossP;
	double bXa = B.x - A.x, bYa = B.y - A.y, bZa = B.z - A.z;
	double t;

	if (vec.x / bXa != vec.y / bYa) {
		double topPart = vec.x * (A.y - pntOfView.y) - vec.y * (A.x - pntOfView.x);
		double bottomPart = vec.y * bXa - vec.x * bYa;
		t = topPart / bottomPart;
	}
	else if (vec.y / bYa != vec.z / bZa) {
		double topPart = vec.y * (A.z - pntOfView.z) - vec.z * (A.y - pntOfView.y);
		double bottomPart = vec.z * bYa - vec.y * bZa;
		t = topPart / bottomPart;
	}
	/*
	double btmPart = vec.x - bXa;
	double topPart = A.x - pntOfView.x;
	if (!btmPart) {
	if (!topPart) crossP.x = crossP.y = crossP.z = 1000;
	else crossP.x = crossP.y = crossP.z = -1000;
	}
	double t = topPart / btmPart;
	*/
	crossP.x = A.x + bXa * t;
	crossP.y = A.y + bYa * t;
	crossP.z = A.z + bZa * t;
	if (bZa) {
		if (A.z < B.z) {
			if (crossP.z < A.z || crossP.z > B.z) crossP.x = crossP.y = crossP.z = -1000;
		}
		else if (crossP.z > A.z || crossP.z < B.z) crossP.x = crossP.y = crossP.z = -1000;
	}
	else if (bXa) {
		if (A.x < B.x) {
			if (crossP.x < A.x || crossP.x > B.x) crossP.x = crossP.y = crossP.z = -1000;
		}
		else if (crossP.x > A.x || crossP.x < B.x) crossP.x = crossP.y = crossP.z = -1000;
	}
	else if (bYa) {
		if (A.y < B.y) {
			if (crossP.y < A.y || crossP.y > B.y) crossP.x = crossP.y = crossP.z = -1000;
		}
		else if (crossP.y > A.y || crossP.y < B.y) crossP.x = crossP.y = crossP.z = -1000;
	}
	else crossP.x = crossP.y = crossP.z = -1000;

	return crossP;
}

POINT crossSheet(POINT & A, POINT & B, POINT & C, POINT & vec, POINT & pntOfView) {
	double bXa = B.x - A.x, bYa = B.y - A.y, bZa = B.z - A.z;
	double cXa = C.x - A.x, cYa = C.y - A.y, cZa = C.z - A.z;
	double bcYZ = bYa * cZa - bZa * cYa, bcXZ = bZa * cXa - bXa * cZa, bcXY = bXa * cYa - bYa * cXa;
	double scalarMult = (pntOfView.x - A.x) * bcYZ + (pntOfView.y - A.y) * bcXZ + (pntOfView.z - A.z) * bcXY;
	double scalarConst = vec.x * bcYZ + vec.y * bcXZ + vec.z * bcXY;

	POINT D;

	if (!scalarConst) {
		if (!scalarMult) {
			POINT cross1, cross2, cross3;

			cross1 = crossLine(A, B, vec, pntOfView);
			cross2 = crossLine(B, C, vec, pntOfView);
			cross3 = crossLine(C, A, vec, pntOfView);

			if (cross1.x > 999) return cross1;
			if (cross2.x > 999) return cross2;
			if (cross3.x > 999) return cross3;

			POINT *vec1 = nullptr;
			POINT *vec2 = nullptr;
			POINT *vec3 = nullptr;
			if (cross1.x > -999) (*vec1).x = cross1.x - pntOfView.x, (*vec1).y = cross1.y - pntOfView.y, (*vec1).z = cross1.z - pntOfView.z;
			if (cross2.x > -999) (*vec2).x = cross2.x - pntOfView.x, (*vec2).y = cross2.y - pntOfView.y, (*vec2).z = cross2.z - pntOfView.z;
			if (cross3.x > -999) (*vec3).x = cross3.x - pntOfView.x, (*vec3).y = cross3.y - pntOfView.y, (*vec3).z = cross3.z - pntOfView.z;

			double len1 = 0, len2 = 0, len3 = 0;
			if (vec1) len1 = sqrt(pow((*vec1).x, 2) + pow((*vec1).y, 2) + pow((*vec1).z, 2));
			if (vec2) len2 = sqrt(pow((*vec2).x, 2) + pow((*vec2).y, 2) + pow((*vec2).z, 2));
			if (vec3) len3 = sqrt(pow((*vec3).x, 2) + pow((*vec3).y, 2) + pow((*vec3).z, 2));

			if (len1 < len2 && len1 < len3) return len2 < len3 ? cross3 : cross2;
			else if (len2 < len1 && len2 < len3) return len1 < len3 ? cross3 : cross1;
			else if (len3 < len1 && len3 < len2) return len1 < len2 ? cross2 : cross1;
			else return cross1;
		}
		else D.x = D.y = D.z = -1000;
		return D;
	}

	double lambda = -scalarMult / scalarConst;

	D.x = pntOfView.x + vec.x * lambda;
	D.y = pntOfView.y + vec.y * lambda;
	D.z = pntOfView.z + vec.z * lambda;

	return D;
}

bool ñrossTrngl(vector <POINT> & vertex, TRIANGLE & trngl, POINT & vec, POINT & pntOfView) {
	POINT a, b, c;
	a.x = vertex[trngl.p1].x;
	b.x = vertex[trngl.p2].x;
	c.x = vertex[trngl.p3].x;
	a.y = vertex[trngl.p1].y;
	b.y = vertex[trngl.p2].y;
	c.y = vertex[trngl.p3].y;
	a.z = vertex[trngl.p1].z;
	b.z = vertex[trngl.p2].z;
	c.z = vertex[trngl.p3].z;

	POINT d;
	d = crossSheet(a, b, c, vec, pntOfView);

	double res1, res2, res3;

	if (a.z == b.z && b.z == c.z) {
		res1 = (a.x - d.x)*(b.y - a.y) - (b.x - a.x)*(a.y - d.y);
		res2 = (b.x - d.x)*(c.y - b.y) - (c.x - b.x)*(b.y - d.y);
		res3 = (c.x - d.x)*(a.y - c.y) - (a.x - c.x)*(c.y - d.y);
	}
	else {
		res1 = (a.z - d.z)*(b.y - a.y) - (b.z - a.z)*(a.y - d.y);
		res2 = (b.z - d.z)*(c.y - b.y) - (c.z - b.z)*(b.y - d.y);
		res3 = (c.z - d.z)*(a.y - c.y) - (a.z - c.z)*(c.y - d.y);
	}

	if ((res1 > 0 && res2 > 0 && res3 > 0) || (res1 < 0 && res2 < 0 && res3 < 0)) return true;
	else return false;
}

POINT crossParall(double x1, double x2, double y1, double y2, double z1, double z2, POINT & vec, POINT & pntOfView) {
	if ((x1 == x2 && vec.x == 0) || (y1 == y2 && vec.y == 0) || (z1 == z2 && vec.z == 0)) {
		bool xInCube = (x1 == x2) && x1 == pntOfView.x,
			yInCube = (y1 == y2) && y1 == pntOfView.y,
			zInCube = (z1 == z2) && z1 == pntOfView.z;
		if (xInCube || yInCube || zInCube) {
			return { -1000, -1000, -1000 };
		}
		return { 1000, 1000, 1000 };
	}
	POINT res;
	if (x1 == x2) {
		res.x = x1;
		res.y = (x1 - pntOfView.x) * vec.y / vec.x + pntOfView.y;
		res.z = (x1 - pntOfView.x) * vec.z / vec.x + pntOfView.z;
	} else 
		if (y1 == y2) {
			res.x = (y1 - pntOfView.y) * vec.x / vec.y + pntOfView.x;
			res.y = y1;
			res.z = (y1 - pntOfView.y) * vec.z / vec.y + pntOfView.z;
		}
		else {
			res.x = (z1 - pntOfView.z) * vec.x / vec.z + pntOfView.x;
			res.y = (z1 - pntOfView.z) * vec.y / vec.z + pntOfView.y;
			res.z = z1;
		}
	return res;
}

bool crossFace(POINT & A, POINT & B, POINT & C, double x1, double x2, double y1, double y2, double z1, double z2, POINT & vec, POINT & pntOfView) {
	POINT crossPnt = crossParall(x1, x2, y1, y2, z1, z2, vec, pntOfView);
	if (crossPnt.x == 1000) {
		return false;
	}
	if (crossPnt.x == -1000) {
		return true;
	}
	bool xInCube = (x1 == x2) ||((crossPnt.x >= x1 - 0.01) && (crossPnt.x <= x2 + 0.01)),
		yInCube = (y1 == y2) || ((crossPnt.y >= y1 - 0.01) && (crossPnt.y <= y2 + 0.01)),
		zInCube = (z1 == z2) || ((crossPnt.z >= z1 - 0.01) && (crossPnt.z <= z2 + 0.01));
	return xInCube && yInCube && zInCube;
}

bool crossCube(Node* curNode, POINT & vec, POINT & pntOfView) {
	POINT blt = { curNode->back, curNode->left, curNode->top },
		blb = { curNode->back, curNode->left, curNode->bottom },
		brt = { curNode->back, curNode->right, curNode->top },
		brb = { curNode->back, curNode->right, curNode->bottom },
		flt = { curNode->front, curNode->left, curNode->top },
		flb = { curNode->front, curNode->left, curNode->bottom },
		frt = { curNode->front, curNode->right, curNode->top },
		frb = { curNode->front, curNode->right, curNode->bottom };
	bool backCrossed = crossFace(blt, blb, brt, curNode->back, curNode->back, curNode->left, curNode->right, curNode->bottom, curNode->top, vec, pntOfView),
		frontCrossed = crossFace(flt, flb, frt, curNode->front, curNode->front, curNode->left, curNode->right, curNode->bottom, curNode->top, vec, pntOfView),
		leftCrossed = crossFace(blt, blb, flt, curNode->back, curNode->front, curNode->left, curNode->left, curNode->bottom, curNode->top, vec, pntOfView),
		rightCrossed = crossFace(brt, brb, frt, curNode->back, curNode->front, curNode->right, curNode->right, curNode->bottom, curNode->top, vec, pntOfView),
		bottomCrossed = crossFace(blb, brb, flb, curNode->back, curNode->front, curNode->left, curNode->right, curNode->bottom, curNode->bottom, vec, pntOfView),
		topCrossed = crossFace(blt, brt, flt, curNode->back, curNode->front, curNode->left, curNode->right, curNode->top, curNode->top, vec, pntOfView);
	if (backCrossed || frontCrossed || leftCrossed || rightCrossed || bottomCrossed || topCrossed) {
		return true;
	}
	return false;
}

int clsst(int trngl1, int trngl2, vector<TRIANGLE> &triangle, vector<POINT> &vertex, POINT & vec, POINT & pntOfView) {
	POINT a1, b1, c1;
	POINT a2, b2, c2;

	a1.x = vertex[triangle[trngl1].p1].x;
	b1.x = vertex[triangle[trngl1].p2].x;
	c1.x = vertex[triangle[trngl1].p3].x;
	a1.y = vertex[triangle[trngl1].p1].y;
	b1.y = vertex[triangle[trngl1].p2].y;
	c1.y = vertex[triangle[trngl1].p3].y;
	a1.z = vertex[triangle[trngl1].p1].z;
	b1.z = vertex[triangle[trngl1].p2].z;
	c1.z = vertex[triangle[trngl1].p3].z;

	a2.x = vertex[triangle[trngl2].p1].x;
	b2.x = vertex[triangle[trngl2].p2].x;
	c2.x = vertex[triangle[trngl2].p3].x;
	a2.y = vertex[triangle[trngl2].p1].y;
	b2.y = vertex[triangle[trngl2].p2].y;
	c2.y = vertex[triangle[trngl2].p3].y;
	a2.z = vertex[triangle[trngl2].p1].z;
	b2.z = vertex[triangle[trngl2].p2].z;
	c2.z = vertex[triangle[trngl2].p3].z;

	POINT d1, d2;
	d1 = crossSheet(a1, b1, c1, vec, pntOfView);
	d2 = crossSheet(a2, b2, c2, vec, pntOfView);

	POINT vec1, vec2;

	vec1.x = d1.x - pntOfView.x, vec1.y = d1.y - pntOfView.y, vec1.z = d1.z - pntOfView.z;
	vec2.x = d2.x - pntOfView.x, vec2.y = d2.y - pntOfView.y, vec2.z = d2.z - pntOfView.z;

	double len1 = sqrt(pow(vec1.x, 2) + pow(vec1.y, 2) + pow(vec1.z, 2));
	double len2 = sqrt(pow(vec2.x, 2) + pow(vec2.y, 2) + pow(vec2.z, 2));

	return len1 < len2 ? trngl1 : trngl2;
}

POINT findVec(POINT & vecOfView, POINT & pntOfView, int xp, int yp) {
	POINT vecToPix, centerOfScreen;

	double Rvp = sqrt(pow(vecOfView.x, 2) + pow(vecOfView.y, 2) + pow(vecOfView.z, 2));
	double k = 1 / Rvp;
	centerOfScreen.x = pntOfView.x + vecOfView.x * k;
	centerOfScreen.y = pntOfView.y + vecOfView.y * k;
	centerOfScreen.z = pntOfView.z + vecOfView.z * k;

	POINT basicX, basicY, basicZ;

	basicZ.x = -vecOfView.x * k, basicZ.y = -vecOfView.y * k, basicZ.z = -vecOfView.z * k;
	POINT crossZ;
	crossZ.z = (basicZ.x * centerOfScreen.x + basicZ.y * centerOfScreen.y + basicZ.z * centerOfScreen.z) / basicZ.z;
	crossZ.x = crossZ.y = 0;

	basicY.x = centerOfScreen.x - crossZ.x, basicY.y = centerOfScreen.y - crossZ.y, basicY.z = centerOfScreen.z - crossZ.z;
	k = 1 / sqrt(pow(basicY.x, 2) + pow(basicY.y, 2) + pow(basicY.z, 2));
	basicY.x *= k, basicY.y *= k, basicY.z *= k;

	basicX.x = basicZ.y * basicY.z - basicZ.z * basicY.y;
	basicX.y = basicZ.z * basicY.x - basicZ.x * basicY.z;
	basicX.z = basicZ.x * basicY.y - basicZ.y * basicY.x;
	k = 1 / sqrt(pow(basicX.x, 2) + pow(basicX.y, 2) + pow(basicX.z, 2));
	basicX.x *= k, basicX.y *= k, basicX.z *= k;

	POINT cnvrtP;
	cnvrtP.x = ((double)xp - 960) / 1000;
	cnvrtP.y = (540 - (double)yp) / 1000;
	cnvrtP.z = 0;

	double m00 = basicX.x,
		m01 = basicY.x,
		m02 = basicZ.x,
		m10 = basicX.y,
		m11 = basicY.y,
		m12 = basicZ.y,
		m20 = basicX.z,
		m21 = basicY.z,
		m22 = basicZ.z;

	POINT pixel;
	pixel.x = m00 * cnvrtP.x + m01 * cnvrtP.y + m02 * cnvrtP.z + centerOfScreen.x;
	pixel.y = m10 * cnvrtP.x + m11 * cnvrtP.y + m12 * cnvrtP.z + centerOfScreen.y;
	pixel.z = m20 * cnvrtP.x + m21 * cnvrtP.y + m22 * cnvrtP.z + centerOfScreen.z;

	vecToPix.x = pixel.x - pntOfView.x;
	vecToPix.y = pixel.y - pntOfView.y;
	vecToPix.z = pixel.z - pntOfView.z;
	return vecToPix;
}

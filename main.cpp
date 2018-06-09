#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <cstring>
#include "string.h"
#include <vector>
#include <cstdlib>
#include "Task.cpp"
using namespace std;

#pragma pack(push, 2)
typedef struct {
	int8_t id[2];            // Завжди дві літери 'B' і 'M'
	int32_t filesize;        // Розмір файла в байтах
	int16_t reserved[2];     // 0, 0
	int32_t headersize;      // 54L для 24-бітних зображень
	int32_t infoSize;        // 40L для 24-бітних зображень
	int32_t width;           // ширина зображення в пікселях
	int32_t depth;           // висота зображення в пікселях
	int16_t biPlanes;        // 1 (для 24-бітних зображень)
	int16_t bits;            // 24 (для 24-бітних зображень)
	int32_t biCompression;   // 0L
	int32_t biSizeImage;     // Можна поставити в 0L для зображень без компрессії (наш варіант)
	int32_t biXPelsPerMeter; // Рекомендована кількість пікселів на метр, можна 0L
	int32_t biYPelsPerMeter; // Те саме, по висоті
	int32_t biClrUsed;       // Для індексованих зображень, можна поставити 0L
	int32_t biClrImportant;  // Те саме
} BMPHEAD;

typedef struct {
	int8_t redComponent;
	int8_t greenComponent;
	int8_t blueComponent;
} PIXELDATA;
#pragma pack(pop)

struct POINT {
	double x, y, z;
};

struct TRIANGLE {
	int p1, p2, p3, n1, n2, n3;
};

void readObj(vector<POINT> &vertex, vector<POINT> &normal, vector<TRIANGLE> &triangle) {
	ifstream in("../cow.obj");
char ch[256], trash[256];
in >> ch;
while (strcmp(ch, "v") && strcmp(ch, "vn") && strcmp(ch, "f")) {
	in.getline(trash, 256, '\n');
	in >> ch;
}
while (!strcmp(ch, "v")) {
	double x, y, z;
	in >> x >> y >> z;
	vertex.push_back({ x, y ,z });
	in >> ch;
}
while (strcmp(ch, "vn") && strcmp(ch, "f")) {
	in.getline(&trash[0], 256, '\n');
	in >> ch;
}
while (!strcmp(ch, "vn")) {
	double x, y, z;
	in >> x >> y >> z;
	normal.push_back({ x, y ,z });
	in >> ch;
}
while (strcmp(ch, "f")) {
	in.getline(&trash[0], 256, '\n');
	in >> ch;
}
while (!strcmp(ch, "f")) {
	int p1, p2, p3, n1, n2, n3, trashInt;
	char trashCh;
	in >> p1 >> trashCh >> trashCh >> n1;
	in >> p2 >> trashCh >> trashCh >> n2;
	in >> p3 >> trashCh >> trashCh >> n3;
	triangle.push_back({ --p1, --p2, --p3, --n1, --n2, --n3 });
	in >> ch;
}
in.close();
}

void readBMP(FILE *inputFile, BMPHEAD &myBMPHead, vector<vector<PIXELDATA>> &pixels) {
	int8_t nullByte = 0;
	fread(&myBMPHead, sizeof(myBMPHead), 1, inputFile);
	pixels.resize(myBMPHead.depth);
	int plusToWidth = 4 - (myBMPHead.width * 3) % 4;
	if (myBMPHead.width % 4 == 0) {
		plusToWidth = 0;
	}
	for (int i = 0; i < myBMPHead.depth; i++) {
		pixels[i].resize(myBMPHead.width);
		for (int j = 0; j < myBMPHead.width; j++) {
			fread(&pixels[i][j], sizeof(pixels[i][j]), 1, inputFile);
		}
		fread(&nullByte, 1, plusToWidth, inputFile);
	}
}

void getSpace(vector<POINT> &vertex, Node* root) {
	if (vertex.empty()) {
		return;
	}
	root->back = vertex[0].x;
	root->front = vertex[0].x;
	root->left = vertex[0].y;
	root->right = vertex[0].y;
	root->bottom = vertex[0].z;
	root->top = vertex[0].z;
	for (int i = 1; i < vertex.size(); i++) {
		double x = vertex[i].x, y = vertex[i].y, z = vertex[i].z;
		if (x < root->back) {
			root->back = x;
		}
		else
			if (x > root->front) {
				root->front = x;
			}
		if (y < root->left) {
			root->left = y;
		}
		else
			if (y > root->right) {
				root->right = y;
			}
		if (z < root->bottom) {
			root->bottom = z;
		}
		else
			if (z > root->top) {
				root->top = z;
			}
	}
}

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

bool сrossTrngl(vector <POINT> & vertex, TRIANGLE & trngl, POINT & vec, POINT & pntOfView) {
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

int clsstInLeaf(Node* curNode, POINT vec, POINT pntOfView, vector<TRIANGLE> &triangle, vector<POINT> &vertex) {
	int res = -1;
	for (int i = 0; i < curNode->numOfTrngls; i++) {
		if (сrossTrngl(vertex, triangle[curNode->triangleCount[i]], vec, pntOfView)) {
			if (res == -1) {
				res = curNode->triangleCount[i];
			}
			else {
				res = clsst(res, curNode->triangleCount[i], triangle, vertex, vec, pntOfView);
			}
		}
	}
	return res;
}

int clsstTrngl(Node* curNode, POINT vec, POINT pntOfView, vector<TRIANGLE> &triangle, vector<POINT> &vertex) {
	if (!crossCube(curNode, vec, pntOfView)) {
		return -1;
	}
	if (curNode->numOfSons == 0) {
		return clsstInLeaf(curNode, vec, pntOfView, triangle, vertex);
	}
	int res = -1, i = 0;
	while (res == -1 && i < 8) {
		if (curNode->Sons[i]) {
			res = clsstTrngl(curNode->Sons[i], vec, pntOfView, triangle, vertex);
		}
		i++;
	}
	i++;
	for (; i < 8; i++) {
		if (curNode->Sons[i]) {
			int nextTrngl = clsstTrngl(curNode->Sons[i], vec, pntOfView, triangle, vertex);
			if (nextTrngl != -1) {
				res = clsst(res, nextTrngl, triangle, vertex, vec, pntOfView);
			}
		}
	}
	return res;
}

bool pInCube(Node* curNode, POINT &p) {
	bool xInCube = (p.x >= curNode->back) && (p.x <= curNode->front);
	bool yInCube = (p.y >= curNode->left) && (p.y <= curNode->right);
	bool zInCube = (p.z >= curNode->bottom) && (p.z <= curNode->top);
	return xInCube && yInCube && zInCube;
}

int difSideTrngl(POINT p1, POINT p2, POINT p3, POINT pn1, POINT pn2, POINT pn3) {
	double x21 = p2.x - p1.x,
		y21 = p2.y - p1.y,
		z21 = p2.z - p1.z,
		x31 = p3.x - p1.x,
		y31 = p3.y - p1.y,
		z31 = p3.z - p1.z,
		A = (y21 * z31 - z21 * y31),
		B = (z21 * x31 - x21 * z31),
		C = (x21 * y31 - y21 * x31),
		D = -p1.x * A - p1.y * B - p1.z * C,
		sign1 = (A * pn1.x + B * pn1.y + C * pn1.z + D),
		sign2 = (A * pn2.x + B * pn2.y + C * pn2.z + D),
		sign3 = (A * pn3.x + B * pn3.y + C * pn3.z + D);
	return (sign1 * sign2 <= 0) || (sign1 * sign3 <= 0);
}

int difSideRctngl(POINT p1, POINT p2, POINT p3, POINT pn1, POINT pn2, POINT pn3, POINT pn4) {
	double x21 = p2.x - p1.x,
		y21 = p2.y - p1.y,
		z21 = p2.z - p1.z,
		x31 = p3.x - p1.x,
		y31 = p3.y - p1.y,
		z31 = p3.z - p1.z,
		A = (y21 * z31 - z21 * y31),
		B = (z21 * x31 - x21 * z31),
		C = (x21 * y31 - y21 * x31),
		D = -p1.x * A - p1.y * B - p1.z * C,
		sign1 = (A * pn1.x + B * pn1.y + C * pn1.z + D),
		sign2 = (A * pn2.x + B * pn2.y + C * pn2.z + D),
		sign3 = (A * pn3.x + B * pn3.y + C * pn3.z + D),
		sign4 = (A * pn4.x + B * pn4.y + C * pn4.z + D);
	return (sign1 * sign2 <= 0) || (sign2 * sign3 <= 0) || (sign3 * sign4 <= 0) || (sign4 * sign1 <= 0);
}

bool trnglInCube(Node* curNode, TRIANGLE &trngl, vector<POINT> &vertex) {
	POINT p1 = vertex[trngl.p1],
		p2 = vertex[trngl.p2],
		p3 = vertex[trngl.p3];
	if (pInCube(curNode, p1) || pInCube(curNode, p2) || pInCube(curNode, p3)) {
		return true;
	}
	POINT blt = { curNode->back, curNode->left, curNode->top },
		blb = { curNode->back, curNode->left, curNode->bottom },
		brt = { curNode->back, curNode->right, curNode->top },
		brb = { curNode->back, curNode->right, curNode->bottom },
		flt = { curNode->front, curNode->left, curNode->top },
		flb = { curNode->front, curNode->left, curNode->bottom },
		frt = { curNode->front, curNode->right, curNode->top },
		frb = { curNode->front, curNode->right, curNode->bottom };
	bool backCrossed = difSideTrngl(blt, blb, brt, p1, p2, p3) && difSideRctngl(p1, p2, p3, blt, blb, brt, brb),
		frontCrossed = difSideTrngl(flt, flb, frt, p1, p2, p3) && difSideRctngl(p1, p2, p3, flt, flb, frt, frb),
		leftCrossed = difSideTrngl(blt, blb, flt, p1, p2, p3) && difSideRctngl(p1, p2, p3, blt, blb, flt, flb),
		rightCrossed = difSideTrngl(brt, brb, frt, p1, p2, p3) && difSideRctngl(p1, p2, p3, brt, brb, frt, frb),
		bottomCrossed = difSideTrngl(blb, brb, flb, p1, p2, p3) && difSideRctngl(p1, p2, p3, blb, brb, flb, frb),
		topCrossed = difSideTrngl(blt, brt, flt, p1, p2, p3) && difSideRctngl(p1, p2, p3, blt, brt, flt, frt);
	if (backCrossed || frontCrossed || leftCrossed || rightCrossed || bottomCrossed || topCrossed) {
		return true;
	}
	return false;
}

void disassembleTrngls(Node* &curNode, vector<TRIANGLE> &triangle, vector<POINT> &vertex, int depth, int maxDepth) {
	curNode->numOfTrngls = curNode->triangleCount.size();
	for (int i = 0; i < 8; i++) {
		Node *temp = createNode(curNode, i);
		for (int j = 0; j < curNode->numOfTrngls; j++) {
			if (trnglInCube(curNode->Sons[i], triangle[curNode->triangleCount[j]], vertex)) {
				curNode->Sons[i]->triangleCount.push_back(curNode->triangleCount[j]);
			}
		}
		if (curNode->Sons[i]->triangleCount.size() == 0) {
			remNode(curNode, i);
		}
	}
	if (curNode->numOfSons > 0) {
		curNode->triangleCount.clear();
	}
	for (int i = 0; i < 8; i++) {
		if (curNode->Sons[i]) {
			if (depth > maxDepth) {
				curNode->Sons[i]->numOfTrngls = curNode->Sons[i]->triangleCount.size();
			}
			else {
				disassembleTrngls(curNode->Sons[i], triangle, vertex, depth + 1, maxDepth);
			}
		}
	}
}

PIXELDATA getColor(TRIANGLE &trngl, POINT &vecLight, vector<TRIANGLE> &triangle, vector<POINT> &vertex, Node* root) {
	POINT p1, p2, p3, pm, po, v1, v2, n;
	p1 = vertex[trngl.p1];
	p2 = vertex[trngl.p2];
	p3 = vertex[trngl.p3];
	pm.x = (p1.x + p2.x) / 2;
	pm.y = (p1.y + p2.y) / 2;
	pm.z = (p1.z + p2.z) / 2;
	po.x = (pm.x + 2 * p3.x) / 3;
	po.y = (pm.y + 2 * p3.y) / 3;
	po.z = (pm.z + 2 * p3.z) / 3;
	v1.x = p1.x - p2.x;
	v1.y = p1.y - p2.y;
	v1.z = p1.z - p2.z;
	v2.x = p2.x - p3.x;
	v2.y = p2.y - p3.y;
	v2.z = p2.z - p3.z;
	n.x = v1.y * v2.z - v1.z * v2.y;
	n.y = v1.z * v2.x - v1.x * v2.z;
	n.z = v1.x * v2.y - v1.y * v2.x;
	double lngth = sqrt(pow(n.x, 2) + pow(n.y, 2) + pow(n.z, 2));
	n.x /= lngth;
	n.y /= lngth;
	n.z /= lngth;
	if (n.y > 0) {
		n.x *= -1;
		n.y *= -1;
		n.z *= -1;
	}
	PIXELDATA res = { abs(n.y) * 255, abs(n.z) * 255, abs(n.x) * 255 };
	/*PIXELDATA res;
	if (clsstTrngl(root, n, po, triangle, vertex) == -1) {
		double cos = (vecLight.x * n.x + vecLight.y * n.y + vecLight.z * n.z) / sqrt(pow(vecLight.x, 2) + pow(vecLight.y, 2) + pow(vecLight.z, 2));
		double shadow = 255 * abs(cos);
		res.redComponent = 255 - shadow;
		res.greenComponent = 255 - shadow;
		res.blueComponent = 255 - shadow;
	}
	else {
		res.redComponent  = 0;
		res.greenComponent  = 0;
		res.blueComponent  = 0;
	}*/
	return res;
}

void calcTrnglsColor(vector<TRIANGLE> &triangle, vector<POINT> &vertex, vector<PIXELDATA> &trnglColor, POINT &pntOfView, POINT &vecOfView, POINT &vecLight, Node* root) {
	for (int i = 0; i < triangle.size(); i++) {
		trnglColor.push_back(getColor(triangle[i], vecLight, triangle, vertex, root));
	}
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

PIXELDATA getPixelColor(int xp, int yp, POINT &pntOfView, POINT &vecOfView, vector< vector<PIXELDATA> > &pixels, vector<PIXELDATA> &trnglColor, Node* root, vector<TRIANGLE> &triangle, vector<POINT> &vertex) {
	PIXELDATA res;
	POINT vec = findVec(vecOfView, pntOfView, yp, xp);
	int numOfTrngl = clsstTrngl(root, findVec(vecOfView, pntOfView, yp, xp), pntOfView, triangle, vertex);
	if (numOfTrngl == -1) {
		res = pixels[xp][yp];
	}
	else {
		res = trnglColor[numOfTrngl];
	}
	return res;
}

void getBMP(FILE *outputFile, BMPHEAD &myBMPHead, vector<vector<PIXELDATA>> &pixels, POINT &vecLight, POINT &pntOfView, POINT &vecOfView, Node* root, vector<TRIANGLE> &triangle, vector<POINT> &vertex, vector<PIXELDATA> &trnglColor) {
	int8_t nullByte = 0;
	int plusToWidth = 4 - (int(abs(myBMPHead.width)) * 3) % 4;
	if (int(abs(myBMPHead.width)) % 4 == 0) {
		plusToWidth = 0;
	}
	fwrite(&myBMPHead, sizeof(myBMPHead), 1, outputFile);
	for (int i = 0; i < myBMPHead.depth; i++) {
		for (int j = 0; j < myBMPHead.width; j++) {
			PIXELDATA newPixel = getPixelColor(i, j, pntOfView, vecOfView, pixels, trnglColor, root, triangle, vertex);
			fwrite(&newPixel, sizeof(newPixel), 1, outputFile);
		}
		if (plusToWidth > 0) {
			fwrite(&nullByte, 1, plusToWidth, outputFile);
		}
	}
}

void createBMP(POINT &vecLight, POINT &pntOfView, POINT &vecOfView, Node* root, vector<TRIANGLE> &triangle, vector<POINT> &vertex, vector<PIXELDATA> &trnglColor) {
	BMPHEAD myBMPHead;
	vector<vector<PIXELDATA>> pixels;
	FILE *inputFile = fopen("images/background.bmp", "rb"), *outputFile = fopen("images/out.bmp", "wb");
	readBMP(inputFile, myBMPHead, pixels);
	getBMP(outputFile, myBMPHead, pixels, vecLight, pntOfView, vecOfView, root, triangle, vertex, trnglColor);
	fclose(inputFile);
	fclose(outputFile);
}

int main() {
	vector<POINT> vertex, normal;
	vector<TRIANGLE> triangle;
	vector<PIXELDATA> trnglColor;
	POINT pntOfView = { 0.1, -1.1, 0.1 }, vecOfView = { -0.1, 1, -0.1 }, vecLight = { 1, 1, 1 };
	int maxDepth = 0;
	readObj(vertex, normal, triangle);
	Node* root = createNode(NULL, 0);
	for (int i = 0; i < triangle.size(); i++) {
		root->triangleCount.push_back(i);
	}
	getSpace(vertex, root);
	calcTrnglsColor(triangle, vertex, trnglColor, pntOfView, vecOfView, vecLight, root);
	disassembleTrngls(root, triangle, vertex, 0, maxDepth);
	cout << 1;
	createBMP(vecLight, pntOfView, vecOfView, root, triangle, vertex, trnglColor);
	system("pause");
}
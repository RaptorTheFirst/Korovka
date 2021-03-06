#include <iostream>
#include <vector>
using namespace std;

struct Node
{
	int numOfSons = 0, numOfTrngls = 0;
	double left, right, front, back, top, bottom;
	vector<int> triangleCount;
	vector<Node*> Sons;
};

Node* createNode(Node *Parent, int y, char axis)
{
	Node *Son = new Node;
	Son->Sons.resize(8, nullptr);
	if (Parent == NULL) {
		Parent = Son;
		return Son;
	}
	Son->back = Parent->back;
	Son->front = Parent->front;
	Son->bottom = Parent->bottom;
	Son->top = Parent->top;
	Son->right = Parent->right;
	Son->left = Parent->left;
	if (axis == 'x') {
		if (y == 0) {
			Son->front = (Parent->back + Parent->front) / 2;
		}
		else {
			Son->back = (Parent->back + Parent->front) / 2;
		}
	}
	else if (axis == 'y') {
		if (y == 0) {
			Son->right = (Parent->left + Parent->right) / 2;
		}
		else {
			Son->left = (Parent->left + Parent->right) / 2;
		}
	}
	else if (axis == 'z') {
		if (y == 0) {
			Son->top = (Parent->bottom + Parent->top) / 2;
		}
		else {
			Son->bottom = (Parent->bottom + Parent->top) / 2;
		}
	}
	Parent->numOfSons++;
	Parent->Sons[y] = Son;
	return Son;
}

void delAll(Node *Parent)
{
	if (Parent->numOfSons > 0) {
		for (int i = 0; i < 8; i++)
		{
			if (Parent->Sons[i] != NULL) delAll(Parent->Sons[i]);
		}
	}
	delete Parent;
}

void remNode(Node *Parent, int y)
{
	if (Parent->Sons[y] != NULL)
	{
		delAll(Parent->Sons[y]);
	}
	Parent->Sons[y] = NULL;
	Parent->numOfSons--;
}
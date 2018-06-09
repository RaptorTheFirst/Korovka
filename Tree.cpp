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

Node* createNode(Node *Parent, int y)
{
	Node *Son = new Node;
	Son->Sons.resize(8, nullptr);
	if (Parent == NULL) {
		Parent = Son;
		return Son;
	}
	if (y == 5 || y == 6 || y == 7 || y == 8) {
		Son->back = Parent->back;
		Son->front = (Parent->back + Parent->front) / 2;
	}
	else {
		Son->back = (Parent->back + Parent->front) / 2;
		Son->front = Parent->front;
	}
	if (y == 1 || y == 3 || y == 5 || y == 7) {
		Son->left = Parent->left;
		Son->right = (Parent->left + Parent->right) / 2;
	}
	else {
		Son->left = (Parent->left + Parent->right) / 2;
		Son->right = Parent->right;
	}
	if (y == 3 || y == 4 || y == 7 || y == 8) {
		Son->bottom = Parent->bottom;
		Son->top = (Parent->bottom + Parent->top) / 2;
	}
	else {
		Son->bottom = (Parent->bottom + Parent->top) / 2;
		Son->top = Parent->top;
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
/*
 * parser.h
 *
 *  Created on: Jul 15, 2012
 *      Author: winnen
 */

#ifndef PARSER_H_
#define PARSER_H_
using namespace std;
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include "Atom.h"
//#include "set.h" //TODO: Use this from source
/*
template <class T>
class Node{
protected:
	Node *head;  //Previous node, if null this is the head
	Node *tail;  //Next node
	T *data;
public:
	Node();
	Node(Node* head, Node* tail, T* obj);
	T* getItem(int counter, int desired) const;
	bool append(T*);
	//Precondition: T is an object of some sort to be stored.
	//Postcondition: T is appended to the
};
template <class T>
class List{
private:
	Node<T> *first;
	Node<T> *last;
	int count;
public:
	List();//Creates an empty list
	int add(T*);
	bool insert(T*, int pos);
		//Inserts the object at the indicated position
	T* operator [](int index);

};*/

//void readPDB(const string);
//string clean(const string&);


#endif /* PARSER_H_ */

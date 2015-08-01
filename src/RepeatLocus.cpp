/*
 * RepeatLocus.cpp
 *
 *  Created on: Oct 30, 2012
 *      Author: spinelli
 */

#include "RepeatLocus.h"

using namespace std;

RepeatLocus::RepeatLocus( string _chromosom, int _start, int _end, string _classe, string _family) {

	chromosom = _chromosom;
	start = _start;
	end = _end;
	classe = _classe;
	family = _family;

}

RepeatLocus::~RepeatLocus() {
	// TODO Auto-generated destructor stub
}

std::string RepeatLocus::getChromosom() const {
	return chromosom;
}

std::string RepeatLocus::getClasse() const {
	return classe;
}

int RepeatLocus::getEnd() const {
	return end;
}

std::string RepeatLocus::getFamily() const {
	return family;
}

int RepeatLocus::getStart() const {
	return start;
}



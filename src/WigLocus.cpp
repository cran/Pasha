/*
 * WigLocus.cpp
 *
 *  Created on: Oct 30, 2012
 *      Author: spinelli
 */

#include "WigLocus.h"

using namespace std;

WigLocus::WigLocus( string chromosom, int start, int step, double value) {

	this->chromosom = chromosom;
	this->start = start;
	this->step = step;
	this->value = value;
}

WigLocus::~WigLocus() {
	// TODO Auto-generated destructor stub
}

std::string WigLocus::getChromosom() const {
	return chromosom;
}

int WigLocus::getStart() const {
	return start;
}

int WigLocus::getStep() const {
	return step;
}

int WigLocus::getValue() const {
	return value;
}

void WigLocus::setValue( double value) {
	this->value = value;
}

int WigLocus::getEnd(){

	return start+step;
}



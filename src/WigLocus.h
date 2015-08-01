/*
 * WigLocus.h
 *
 *  Created on: Oct 30, 2012
 *      Author: spinelli
 */

#ifndef WIGLOCUS_H_
#define WIGLOCUS_H_

#include <string>

class WigLocus {
public:
	WigLocus( std::string chromosom, int start, int step, double value);
	virtual ~WigLocus();
	std::string getChromosom() const;
	int getStart() const;
	int getStep() const;
	int getValue() const;
	void setValue( double value);
	int getEnd();

private:
	std::string chromosom;
	int start;
	int step;
	int value;
};

#endif /* WIGLOCUS_H_ */

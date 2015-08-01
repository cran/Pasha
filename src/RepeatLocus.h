/*
 * RepeatLocus.h
 *
 *  Created on: Oct 30, 2012
 *      Author: spinelli
 */

#ifndef REPEATLOCUS_H_
#define REPEATLOCUS_H_

#include <string>

class RepeatLocus {
public:
	RepeatLocus( std::string chrom, int start, int end, std::string classe, std::string family);
	virtual ~RepeatLocus();
	std::string getChromosom() const;
	std::string getClasse() const;
	int getEnd() const;
	std::string getFamily() const;
	int getStart() const;

private:
	std::string chromosom;
	int start;
	int end;
	std::string classe;
	std::string family;
};

#endif /* REPEATLOCUS_H_ */

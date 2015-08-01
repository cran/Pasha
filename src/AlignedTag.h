/*
 * AlignedTag.h
 *
 *  Created on: Sep 5, 2012
 *      Author: spinelli
 */

#include <sstream>
#include <string>

#ifndef ALIGNEDTAG_H_
#define ALIGNEDTAG_H_

class AlignedTag {
public:
	AlignedTag( int, char, int, int);
	virtual ~AlignedTag();
	int getChromID() const;
	int getPosition() const;
	char getStrand() const;
	int getTagID() const;
	const char* output() const;
	std::string toString();
	int getRepeat() const;
	void setRepeat( int multiread);

private:
	int tagID;
	char strand;
	int chromID;
	int position;
	int repeat;
};

#endif /* ALIGNEDTAG_H_ */

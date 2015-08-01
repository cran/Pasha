/*
 * AlignedTag.cpp
 *
 *  Created on: Sep 5, 2012
 *      Author: spinelli
 */

#include "AlignedTag.h"

using namespace std;

AlignedTag::AlignedTag( int _tag_id, char _strand, int _chrom_id, int _position) {

	tagID = _tag_id;
	strand = _strand;
	chromID = _chrom_id;
	position = _position;
	repeat = 0;
}

AlignedTag::~AlignedTag() {
	// TODO Auto-generated destructor stub
}

int AlignedTag::getChromID() const {
	return chromID;
}

int AlignedTag::getPosition() const {
	return position;
}

char AlignedTag::getStrand() const {
	return strand;
}

int AlignedTag::getTagID() const {
	return tagID;
}

int AlignedTag::getRepeat() const {
	return repeat;
}

void AlignedTag::setRepeat( int multiread) {
	this->repeat = multiread;
}


string AlignedTag::toString() {

	stringstream result;
	result << "Tokens = " << tagID << ";" << strand << ";" << chromID << ";" << position << ";" << position;

	return result.str();
}




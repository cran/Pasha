#ifndef POSPAIR
#define POSPAIR

struct PosPair {
  char cid; // chromosone id
  int pos;  // position

  PosPair() { cid = pos = 0; }

  PosPair(char cid, int pos) {
    this->cid = cid;
    this->pos = pos;
  }

  bool operator< (const PosPair& o) const {
    return cid < o.cid || cid == o.cid && pos < o.pos;
  }

  bool operator== (const PosPair& o) const {
    return cid == o.cid && pos == o.pos;
  }
};

#endif

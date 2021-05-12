#include <iostream>
#include <fstream>
#include <cassert>
#include <bit>
using namespace std;

#ifndef NUM_VERTS
#define NUM_VERTS 9
#endif

#define bit(i) (((Set)1) << (i))

typedef unsigned int Set;
typedef int Num;

void writed6(ostream &o, Set *grf, int n) {
  int i, j, b, val;
  o.put('&');
  if (n <= 62)
    o.put(63 + n);
  else if (n <= 63 << 12) {
    o.put('~');
    o.put(63 + (n >> 12));
    o.put(63 + ((n >> 6) & 63));
    o.put(63 + (n & 63));
  } else {
    o.put('~');
    o.put('~');
    o.put(63 + (n >> 30));
    o.put(63 + ((n >> 24) & 63));
    o.put(63 + ((n >> 18) & 63));
    o.put(63 + ((n >> 12) & 63));
    o.put(63 + ((n >> 6) & 63));
    o.put(63 + (n & 63));
  }
  b = 6;
  val = 63;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      b--;
      if (b < 0) {
	o.put(val);
	val = 63;
	b = 5;
      }
      if (grf[i] & bit(j))
	val += 1 << b;
    }
  }
  if (n > 1)
    o.put(val);
}

constexpr Num calcImbalanceMod(Num n) {
  Num val = 1;
  while ((val * val) % n != 0)
    val++;
  return val;
}

static const Num IMBAL_MOD = calcImbalanceMod(NUM_VERTS);
static const Num NUM_I_GROUPS = 2 * NUM_VERTS / IMBAL_MOD;
#define SIGN(v) ((v) >= 0 ? 1 : -1)
#define SQR(v) ((v) * (v))

static_assert(NUM_VERTS % IMBAL_MOD == 0);

class ImbalGen {
public:
  Num imbal[NUM_VERTS], imbalGrp[NUM_VERTS], grpCnt[NUM_I_GROUPS], grpVal[NUM_I_GROUPS], imbalClass;
  Set grpMsk[NUM_I_GROUPS];
  bool imbalValid, imbalDone;
  ImbalGen() : imbalValid(false), imbalDone(false) {}
  void setClass(Num imbalClass_){
    imbalValid = imbalDone = false;
    imbalClass=imbalClass_;
    signICls = SIGN(imbalClass);
    rdcImbal[0] = -NUM_VERTS/IMBAL_MOD*signICls;
    for (Num i = 1; i<NUM_VERTS/IMBAL_MOD; i++){
      rdcImbal[2*i-1] = (NUM_VERTS/IMBAL_MOD-i)*signICls;
      rdcImbal[2*i] = -(NUM_VERTS/IMBAL_MOD-i)*signICls;
    }
    rdcImbal[NUM_I_GROUPS-1] = 0;
    for (Num i=0; i < NUM_I_GROUPS;i++)
      grpVal[i]=rdcImbal[i]*IMBAL_MOD+imbalClass;
    curDelta[0]=imbalClass*NUM_VERTS/IMBAL_MOD;
    level=0;
    grpCnt[0]=-1;
    numRem[0]=NUM_VERTS;
    grpMsk[0]=0;
    assert(rdcImbal[NEAR_END] != 0);
    assert(rdcImbal[NEAR_END+1] == 0);
  }
  bool next(){
    if (NEAR_END == 0 && !imbalValid)
      return finish();
    if(level>=NEAR_END)
      return up();
    grpCnt[level]++;
    if(grpCnt[level]>0) {
      if (numRem[level] == 0)
        return up();
      grpMsk[level]|=bit(numRem[level]-1);
      imbal[numRem[level]-1]=grpVal[level];
      imbalGrp[numRem[level]-1]=level;
      curDelta[level]+=rdcImbal[level];
      numRem[level]--;
    }
    if (((numRem[level]*rdcImbal[level+1]+curDelta[level])*SIGN(rdcImbal[level+1])<0) || (grpCnt[level]+min(abs(grpVal[level]),NUM_VERTS)>NUM_VERTS))
      return up();
    return down();
  }
  protected:
  Num rdcImbal[NUM_I_GROUPS], signICls, curDelta[NUM_I_GROUPS], level, numRem[NUM_I_GROUPS];
  static const Num NEAR_END = NUM_I_GROUPS-2;
  bool up() {
    if (level == 0) {
      imbalValid = false;
      imbalDone = true;
      return false;
    }
    level--;
    return next();
  }
  bool down() {
    level++;
    grpCnt[level]=-1;
    grpMsk[level]=0;
    curDelta[level]=curDelta[level-1];
    numRem[level]=numRem[level-1];
    if (level == NEAR_END)
      return finish();
    return next();
  }
  bool finish() {
    grpCnt[NEAR_END]=-curDelta[NEAR_END]/rdcImbal[NEAR_END];
    if (grpCnt[NEAR_END]<0) {
      return up();
    }
    assert(grpCnt[NEAR_END]<=numRem[NEAR_END]);
    grpMsk[NEAR_END+1]=0;
    numRem[NEAR_END]-=grpCnt[NEAR_END];
    curDelta[NEAR_END]+=rdcImbal[NEAR_END]*grpCnt[NEAR_END];
    assert(curDelta[NEAR_END]==0);
    for(Num i=numRem[NEAR_END];i<(NEAR_END == 0 ? NUM_VERTS : numRem[NEAR_END-1]);i++) {
      imbal[i]=grpVal[NEAR_END];
      imbalGrp[i]=NEAR_END;
      grpMsk[NEAR_END]|=bit(i);
    }
    grpCnt[NEAR_END+1]=numRem[NEAR_END];
    for(Num i=0; i<numRem[NEAR_END]; i++) {
      imbal[i]=grpVal[NEAR_END+1];
      imbalGrp[i]=NEAR_END+1;
      grpMsk[NEAR_END+1]|=bit(i);
    }
    imbalValid = true;
    return true;
  }
};
class ImbalGen2 : public ImbalGen{
public:
  ImbalGen2() : ImbalGen(), feasValid(false), feasDone(false) {}
  bool next() { feasValid = feasDone = false; return ImbalGen::next(); }
  bool calcFeas() {
    feasDone = true;
    for (Num grp=0;grp<NUM_I_GROUPS;grp++){
      feasCnt[grp]=0;
      feasLoc[grp][1 << NUM_VERTS] = 0;
      if (grpCnt[grp]==0) continue;
      Num gDelta[NUM_I_GROUPS],baseHigh=SQR(NUM_VERTS-grpVal[grp]),baseLow=SQR(-grpVal[grp]),iDelta=-(baseHigh-baseLow)*grpVal[grp];
      for (Num i=0;i<NUM_I_GROUPS;i++){
        Num lDiff=SQR(-grpVal[i])-baseLow,hDiff=SQR(NUM_VERTS-grpVal[i])-baseHigh;
        gDelta[i]=hDiff-lDiff;
        iDelta+=lDiff*grpCnt[i];
      }
      for (Set s=0;s<(1 << NUM_VERTS);s++){
        Num sDelta=iDelta;
        for (Num i = 0; i < NUM_I_GROUPS; i++)
          sDelta += gDelta[i] * popcount(s & grpMsk[i]);
        feasLoc[grp][s] = feasCnt[grp];
        if (sDelta == 0)
          feasVal[grp][feasCnt[grp]++] = s;
      }
      feasLoc[grp][1 << NUM_VERTS] = feasCnt[grp];
      if (feasCnt[grp] == 0)
        return false;
    }
    feasValid = true;
    return true;
  }
  Set feasVal[NUM_I_GROUPS][1 << NUM_VERTS];
  Num feasLoc[NUM_I_GROUPS][(1 << NUM_VERTS) + 1], feasCnt[NUM_I_GROUPS];
  bool feasValid, feasDone;
};

class GrfGen {
public:
  ImbalGen2 ig;
  Set grf[NUM_VERTS];
  bool grfInited, grfValid, grfDone;
  GrfGen() : ig(), grfInited(false), grfValid(false), grfDone(false) {
    initStatics();
  }
  void setClass(Num imbalClass) { invalidate(); ig.setClass(imbalClass); }
  bool nextImbalSet() { invalidate(); return ig.next(); }
  bool nextGrf() {
    assert(grfInited);
    return next();
  }

  void initGrf() {
    assert(ig.feasValid);
    vert = NUM_VERTS - 1;
    for (Num i = 0; i < NUM_VERTS; i++) {
      colBase[i] = rowBase[i] = 0;
      rowIter[i] = colIter[i] = 0;
    }
    rowIter[vert] = -1;
    colIter[vert] = 0;
    stopRowIter[vert] = ig.feasLoc[ig.imbalGrp[vert]][1 << (NUM_VERTS - 1)];
    for (Set s = 0; s < (1 << NUM_VERTS); s++) {
      Num cSum = 0;
      for (Num i = 0; i < NUM_I_GROUPS; i++)
        cSum += NUM_VERTS * ig.grpVal[i] * popcount(s & ig.grpMsk[i]);
      rowImbalProd[s] = cSum;
    }
    imbalBase0 = 0;
    for (Num i = 0; i < NUM_I_GROUPS; i++)
      imbalBase0 += ig.grpCnt[i] * SQR(ig.grpVal[i]);
    for (Num i = 0; i < NUM_I_GROUPS; i++) {
      for (Num j = 0; j < NUM_I_GROUPS; j++) {
        imbalBase[i][j] = imbalBase0 + NUM_VERTS * ig.grpVal[i] * ig.grpVal[j];
      }
    }
    grfInited = true;
  }

  void invalidate() { grfInited = grfValid = grfDone = false; }

protected:
  static Num cBound[NUM_VERTS];
  static Set colWithBits[(1 << (NUM_VERTS - 1)) + NUM_VERTS];
  static bool staticsInited;
  static void initStatics() {
    if (staticsInited)
      return;
    Num p = 0;
    colWithBits[0] = 5;
    for (Num nBits = 0; nBits <= NUM_VERTS - 1; nBits++) {
      cBound[nBits] = p;
      for (Set s = 0; s < (1 << (NUM_VERTS - 1)); s++) {
        if (popcount(s) == nBits)
          colWithBits[p++] = s;
      }
      colWithBits[p++] = 1 << NUM_VERTS;
    }
    staticsInited = true;
  }

  Num vert, rowIter[NUM_VERTS], stopRowIter[NUM_VERTS], colIter[NUM_VERTS];
  Set colBase[NUM_VERTS], rowBase[NUM_VERTS], cGrf[NUM_VERTS];
  Num rowImbalProd[1 << NUM_VERTS], imbalBase0, imbalBase[NUM_I_GROUPS][NUM_I_GROUPS], sRowBits[NUM_VERTS];
  Num reqShrBits[NUM_VERTS][NUM_VERTS];
  bool lessBits[NUM_VERTS];

  // I know, gotos are considered bad style (or something).  They used to be tail calls, but tail-call elimination isn't guaranteed, and the code was triggering stack overflows. :(
  bool next() {
  next:
    grfValid = true;
    if (vert == 0)
      goto up;
    colIter[vert]++;
    Set newCol;
    while (colWithBits[colIter[vert]] >= (1 << vert)) {
      rowIter[vert]++;
      if (rowIter[vert] >= stopRowIter[vert])
        goto up;
      Set newRow = ig.feasVal[ig.imbalGrp[vert]][rowIter[vert]];
      Num rowBits = popcount(newRow);
      Num moreColBits = rowBits - popcount(colBase[vert]) + ig.imbal[vert];
      if (moreColBits > vert || moreColBits < 0 || (vert < NUM_VERTS - 1 && ((lessBits[vert + 1] && rowBits <= sRowBits[vert + 1]) || (rowBits < sRowBits[vert + 1] && ig.imbalGrp[vert] == ig.imbalGrp[vert+1] && (rowBase[vert] & ~bit(vert+1)) == rowBase[vert+1]))))
        continue;
      bool badRow = false;
      for (Num i = vert + 1; i < NUM_VERTS; i++) {
        Num testVal = imbalBase[ig.imbalGrp[vert]][ig.imbalGrp[i]]-rowImbalProd[newRow]-rowImbalProd[grf[i]]+NUM_VERTS*(sRowBits[i]*ig.imbal[i]+rowBits*ig.imbal[vert]);
        Num rTestVal = testVal / SQR(NUM_VERTS);
        if (rTestVal * SQR(NUM_VERTS) != testVal) {
          badRow = true;
          break;
        }
        rTestVal += popcount(newRow & grf[i]) + (((newRow >> i) & 1) - ((grf[i] >> vert) & 1)) * (rowBits - sRowBits[i]);
        if (rTestVal < 0 || rTestVal > ig.imbal[i] + sRowBits[i] || rTestVal > ig.imbal[vert] + rowBits) {
          badRow = true;
          break;
        }
        reqShrBits[vert][i] = rTestVal;
      }
      if (badRow)
        continue;
      grf[vert] = newRow;
      sRowBits[vert] = rowBits;
      colIter[vert] = cBound[moreColBits];
      for (Num i = 0; i < vert; i++) {
        Set cBit = (newRow >> i) & 1;
        colBase[i] = (colBase[i] & ~bit(vert)) | (cBit << vert);
      }
    }
    newCol = colBase[vert] | colWithBits[colIter[vert]];
    for (Num i = 0; i < vert; i++) {
      Set cBit = (newCol >> i) & 1;
      rowBase[i] = (rowBase[i] & ~bit(vert)) | (cBit << vert);
      if (ig.imbalGrp[i] == ig.imbalGrp[i+1] && (rowBase[i] & ~bit(vert)) == (rowBase[i+1] & ~bit(vert)) && cBit == 0 && ((newCol >> (i + 1)) & 1) != 0)
        goto next;
      if (ig.feasLoc[ig.imbalGrp[i]][rowBase[i]] == ig.feasLoc[ig.imbalGrp[i]][rowBase[i]+bit(vert)])
        goto next;
    }
    for (Num i = vert + 1; i < NUM_VERTS; i++) {
      Num val = popcount(cGrf[i] & newCol);
      if (val != reqShrBits[vert][i])
        goto next;
    }
    lessBits[vert] = (ig.imbalGrp[vert-1] == ig.imbalGrp[vert] && (rowBase[vert-1] & ~bit(vert)) == rowBase[vert] && ((newCol >> (vert - 1)) & 1) != 0 && ((grf[vert] >> (vert - 1)) & 1) == 0);
    cGrf[vert] = newCol;

  down:
    vert--;
    if (vert == 0)
      goto finish;
    colIter[vert] = 0;
    rowIter[vert] = ig.feasLoc[ig.imbalGrp[vert]][rowBase[vert]] - 1;
    stopRowIter[vert] = ig.feasLoc[ig.imbalGrp[vert]][rowBase[vert]+bit(vert)];
    goto next;

  up:
    vert++;
    if (vert == NUM_VERTS) {
      grfValid = false;
      grfDone = true;
      return false;
    }
    for (Num i = 0; i < vert - 1; i++) {
      rowBase[i] &= ~bit(vert - 1);
      colBase[i] &= ~bit(vert - 1);
    }
    goto next;

  finish:
    grf[0] = rowBase[0];
    return true;
  }
};
bool GrfGen::staticsInited = false;
Num GrfGen::cBound[NUM_VERTS];
Set GrfGen::colWithBits[(1 << (NUM_VERTS - 1)) + NUM_VERTS];

int main() {
  GrfGen gg;
  for (Num imbalClass = -IMBAL_MOD/2; imbalClass<=(IMBAL_MOD-1)/2; imbalClass++) {
    gg.setClass(imbalClass);
    while (gg.ig.next()) {
#ifdef SKIP_NORMAL
      if (imbalClass == 0 && gg.ig.grpCnt[NUM_I_GROUPS - 1] == NUM_VERTS)
        continue;
#endif
#ifdef ONLY_NORMAL
      if (imbalClass != 0 || gg.ig.grpCnt[NUM_I_GROUPS - 1] != NUM_VERTS)
        continue;
#endif
#ifdef SKIP_JOINS
      Num sumA = 0, sumB = 0;
      bool skip = false;
      for (Num i = 0; i < NUM_I_GROUPS; i += 2) {
        sumA += gg.ig.grpCnt[i];
        sumB += gg.ig.grpCnt[i+1];
        if ((sumA + min(abs(gg.ig.grpVal[i]), NUM_VERTS - 1) > NUM_VERTS - 1) || (sumB + min(abs(gg.ig.grpVal[i+1]), NUM_VERTS - 1) > NUM_VERTS - 1)) {
          skip = true;
          break;
        }
      }
      if (skip)
        continue;
#endif
      if (!gg.ig.calcFeas())
        continue;
      gg.initGrf();
      while (gg.nextGrf()) {
        writed6(cout, gg.grf, NUM_VERTS);
        cout << "\n";
      }
    }
  }
  return 0;
}

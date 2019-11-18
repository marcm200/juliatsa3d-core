#include "math.h"
#include "stdio.h"
#include "string.h"
#include "quadmath.h"
#include "stdint.h"


// /////////////////////////////////////////

// comment in or out depending on the desired floating
// point datatype

#define _DOUBLE
//#define _LONGDOUBLE
//#define _QUADPREC

// /////////////////////////////////////////

typedef uint16_t *PDBYTE;

#ifdef _QUADPREC
typedef __float128 NTYP;
const char NNTYPSTR[]="f128_";
const int32_t _BITSNTYP=113;
#endif
#ifdef _LONGDOUBLE
typedef long double NTYP;
const char NNTYPSTR[]="ld_";
const int32_t _BITSNTYP=63;
#endif
#ifdef _DOUBLE
typedef double NTYP;
const char NNTYPSTR[]="";
const int32_t _BITSNTYP=53;
#endif

enum { CMD_CALC=1,CMD_PERIODICITY=2 };

enum {
	FUNC_BAIRD=0,FUNC_BRISTOR=1,FUNC_MAKIN=2,FUNC_SMITH=3,
	FUNC_BAIRD5=4,FUNC_TRICZ3B=5,FUNC_TRICZ4B=6,
	FUNC_MAKINEXP4=7,FUNC_MAKINEXP4B=8,FUNC_TRICZ5D=9,
	FUNC_MAKINEXP5=10
};

const int32_t funcanz=FUNC_MAKINEXP5+1;
const int32_t MAXPTR=2048;

const char funcname[][32] = {
	"BAIRD","BRISTOR","MAKIN","SMITH",
	"BAIRD5","TRICZ3B","TRICZ4B",
	"MAKINEXP4","MAKINEXP4B","TRICZ5D","MAKINEXP5"
};

int32_t SCREENWIDTH=256;
NTYP RANGE0=-4;
NTYP RANGE1= 4;

const int32_t DBYTEMAX=(1 << 16);
const int32_t MAXFATOUCOMPONENTS=65500;
const int32_t MAXCYCLES=110;
const int32_t FATOUCOMPONENTCOLOROFFSET=24;
// for 64bit operating systems
const int64_t CHUNKSIZE=(int64_t)1 << 30; 

// for 32bit operating systems
//const int64_t CHUNKSIZE=(int64_t)1 << 29; 

// color constants for cells
const uint16_t CC_GRAY=0;
const uint16_t CC_WHITE=0b01;
const uint16_t CC_BLACK=0b10;
const uint16_t CC_GRAYPOTW=0b11;


#define LOGMSG(TT) {\
	printf(TT);\
	fprintf(flog,TT); fflush(flog); \
}

#define LOGMSG2(TT,AA) {\
	printf(TT,AA);\
	fprintf(flog,TT,AA); fflush(flog); \
}

#define TRIM(WW)\
	if ( WW < 0) WW=0;\
	else if (WW >= SCREENWIDTH) WW=SCREENWIDTH-1;

#define LOGMSG3(TT,AA,BB) {\
	printf(TT,AA,BB);\
	fprintf(flog,TT,AA,BB); fflush(flog); \
}

#define LOGMSG4(TT,AA,BB,CC) {\
	printf(TT,AA,BB,CC);\
	fprintf(flog,TT,AA,BB,CC); fflush(flog); \
}

#define LOGMSG5(TT,AA,BB,CC,DD) {\
	printf(TT,AA,BB,CC,DD);\
	fprintf(flog,TT,AA,BB,CC,DD); fflush(flog); \
}


// a voxel representing a tiny cube in R^3-space
// interval arithmetics boundaries
// x-coordinate [x0..x1], y in [y0..y1] etc.
// all values finite
// zero is NOT in any interval, at most it is one/more of the border values
struct SpaceCube {
	NTYP x0,x1,y0,y1,z0,z1;
};

// memory manager for allocating DBYTE values
// low memory fragmentation, in chunks
struct ArrayMgrDByte {
	uint16_t *current;
	int32_t freeFromIdx;
	int32_t allocateHowMany;
	PDBYTE ptr[MAXPTR];
	int32_t countptr;
	
	ArrayMgrDByte();
	virtual ~ArrayMgrDByte();
	uint16_t *getMem(const int32_t);
};

// coordinates in cube
struct ScreenCube {
	int32_t x0,x1,y0,y1,z0,z1;
};

// an attracting periodic cycle
struct Cycle {
	int32_t len;
	uint16_t immediateBasinColorIdx;
	uint16_t attractionBasinColorIdx;
	// currently not used
	uint16_t fatouidx0,fatouidx1; // Cycle between those indices in array allcycles[]
};

// a Fatou component
// either an immediate basin or some other part of the attraction basin
struct FatouComponent {
	ScreenCube scrc;
	uint16_t currentOrbitColorIdxTemp; // < 256 => set, >= 256: temporary
	int32_t inCycleNbr; // < 0 => none
	// currently not used
	char isimmediate; // flag
};

// for internal viewer: one 2D screen pixel and its
// pointed-to voxel coordinate x,y,z in floating type
struct  Coord {
	float x,y,z; // low memory consumption: double not necessary
	float dimmingFactor; // dimming of RGB values dependent on distance and normal vector
	uint8_t R,G,B;
	int32_t intersectsWithVisibleObjectAtStep; // already intersected with object
		
	Coord(const Coord&);
	Coord(const double,const double,const double);
	Coord();
	double norm(void);
	double normQ(void);
	Coord& operator+=(const Coord&);
	friend Coord operator+(Coord lhs,const Coord& rhs) { lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z; return lhs; }
	friend Coord operator-(Coord lhs,const Coord& rhs) { lhs.x -= rhs.x; lhs.y -= rhs.y; lhs.z -= rhs.z; return lhs; }

	Coord& operator=(const Coord&);
	void normiere(void);
	void mult(const double);
};

typedef Coord* PCoord;

// memory manager for allocating Coords in screen
struct ArrayMgrCoord {
	Coord *current;
	int32_t freeFromIdx;
	int32_t allocateHowMany;
	PCoord ptr[MAXPTR];
	int32_t countptr;
	
	ArrayMgrCoord();
	virtual ~ArrayMgrCoord();
	Coord *getMem(const int32_t);
};

// internal viewer. based on my cube-viewer project
// fixed observer position
struct Screen {
	int32_t lenx,leny;
	double abstandDelta;
	double normVektor;
	PCoord* coordsY; 
	Coord vektor;
	Coord scrxv,scryv;
	Coord observer;
		
	Screen();
	virtual ~Screen();
	void setScreenSize(const int32_t,const int32_t);
	void setObserver(Coord&,Coord&,Coord&,Coord&);
	void slide(void);
	void dimm(const double dimmingFactor,const int32_t s0,const int32_t s1);
	void dimmByNormal(void);
	double DistanceToScreen(Coord&);
};

struct RGB {
	uint8_t R,G,B;
};

// helper object for speeding up finding Fatou components
// storing voxels to follow in flood fill
struct Int3 {
	int32_t x,y,z;
};

struct Zeile {
	int32_t memx0,memx1;
	uint16_t *werte; // allocated by memory manager
	int32_t gb0,gb1; // index region where gray cells are
	
	void setToNonExistent(void);
};

// main struct that computes the 3D object
struct Data3Ddyn {
	// every gray or black cells is contained in this cube
	int32_t encx0,encx1,ency0,ency1,encz0,encz1;
	Zeile *planesZY; // indexed by Z*SCREENWIDTH+Y

	Data3Ddyn();
	virtual ~Data3Ddyn();
	
	uint16_t getVoxel(const int32_t,const int32_t,const int32_t);
	void setVoxel(const int32_t,const int32_t,const int32_t,const uint16_t);
	void setGrayregionYZ(const int32_t,const int32_t,const int32_t,const int32_t);
	void convert_to_BitCube(const char*);
	void initComplete(void);
	int readRaw(void);
	void saveRaw(const char*);

	// main routines to compute a cell's color
	void search_specext(void);
	void search_directWhite(void);
	int coloring_black(void);
};


// global variables
double scaleRangePerPixel,scalePixelPerRange;
int64_t DENOM225=1 << 25;
// coefficients for functions
NTYP tricC0x,tricC1x,tricC0y,tricC1y,tricC0z,tricC1z;
NTYP tricAx,tricAy,tricAz;
NTYP tricBx,tricBy,tricBz;
NTYP varE=0.0,varF=0.0,varG=0.0,varH=0.0;
NTYP varJ=0.0,varK=0.0,varL=0.0,varM=0.0;
NTYP varN=0.0,varO=0.0,varP=0.0,varQ=0.0;
// main object
Data3Ddyn* data;
// pointer to current used function
int _FUNC=FUNC_MAKIN;
void (*getBoundingBoxfA)(SpaceCube&,SpaceCube&) = NULL;
FILE* flog=NULL;
int64_t ctrweiss=0,ctrgrau=0,ctrpotwgrau=0,ctrschwarz=0;
// memory managers
ArrayMgrDByte mgrVoxel;
ArrayMgrCoord mgrCoord;
int bitpower=1;
int GRAY_IS_TRANSPARENT=1;
int _TWDBITCUBE=-1;
int CMD=CMD_PERIODICITY;
int REFINEMENTLEVEL=8;
RGB pal256[256];
// Cycles if found
Cycle *cycles=NULL;
int32_t anzcycles=0;


// functions

inline int32_t minimumI(const int32_t a,const int32_t b) {
	if (a < b) return a;
	return b;
}

inline int32_t maximumI(const int32_t a,const int32_t b) {
	if (a > b) return a;
	return b;
}

char* upper(char* s) {
	if (!s) return 0;
	for(int32_t i=(strlen(s)-1);i>=0;i--) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
	}

	return s;
}

// converter function: triplex number in virtual screen coordinates
inline int32_t scrcoord0(const NTYP a) {
	return (int32_t)floor( (a - RANGE0) * scalePixelPerRange );
}

inline NTYP maximumD(const NTYP a,const NTYP b,const NTYP c,const NTYP d) {
	NTYP m=a;
	if (b > m) m=b;
	if (c > m) m=c;
	if (d > m) m=d;
	return m;
}

inline NTYP minimumD(const NTYP a,const NTYP b,const NTYP c,const NTYP d) {
	NTYP m=a;
	if (b < m) m=b;
	if (c < m) m=c;
	if (d < m) m=d;
	return m;
}

inline NTYP minimumD(const NTYP a,const NTYP b) {
	if (a < b) return a;
	return b;
}

inline NTYP maximumD(const NTYP a,const NTYP b) {
	if (a > b) return a;
	return b;
}

// implemented 3d functions
// bounding cubes in interval arithmetics
// formulas generated by a symbolic parser
inline void getBoundingBoxfA_tricz5d(SpaceCube& A,SpaceCube& fA) {
	fA.x0=tricC0x+minimumD(varF*A.x0,varF*A.x1)+minimumD(varH*minimumD(A.z0*A.z0,A.z1*A.z1),varH*maximumD(A.z0*A.z0,A.z1*A.z1))+A.x0*A.x0*A.x0*A.x0*A.x0-(2*(5*maximumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*minimumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1));
	fA.x1=tricC1x+maximumD(varF*A.x0,varF*A.x1)+maximumD(varH*minimumD(A.z0*A.z0,A.z1*A.z1),varH*maximumD(A.z0*A.z0,A.z1*A.z1))+A.x1*A.x1*A.x1*A.x1*A.x1-(2*(5*minimumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*maximumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1));
	fA.y0=tricC0y+minimumD(varE*A.y0,varE*A.y1)+5*minimumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y0*A.y0*A.y0*A.y0*A.y0;
	fA.y1=tricC1y+maximumD(varE*A.y0,varE*A.y1)+5*maximumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y1*A.y1*A.y1*A.y1*A.y1;
	fA.z0=tricC0z+minimumD(A.z0*A.z0,A.z1*A.z1)-maximumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1);
	fA.z1=tricC1z+maximumD(A.z0*A.z0,A.z1*A.z1)-minimumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1);
}

void getBoundingBoxfA_bristor(SpaceCube& A,SpaceCube& fA) {
	fA.x0=minimumD(A.x0*A.x0,A.x0*A.x1,A.x1*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y0*A.y1,A.y1*A.y0,A.y1*A.y1)-maximumD(A.z0*A.z0,A.z0*A.z1,A.z1*A.z0,A.z1*A.z1)+tricC0x;
	fA.x1=maximumD(A.x0*A.x0,A.x0*A.x1,A.x1*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y0*A.y1,A.y1*A.y0,A.y1*A.y1)-minimumD(A.z0*A.z0,A.z0*A.z1,A.z1*A.z0,A.z1*A.z1)+tricC1x;
	fA.y0=2*minimumD(A.y0*A.x0,A.y0*A.x1,A.y1*A.x0,A.y1*A.x1)-maximumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1)+tricC0y;
	fA.y1=2*maximumD(A.y0*A.x0,A.y0*A.x1,A.y1*A.x0,A.y1*A.x1)-minimumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1)+tricC1y;
	fA.z0=tricC0z+2*minimumD(A.z0*A.x0,A.z0*A.x1,A.z1*A.x0,A.z1*A.x1)+minimumD(A.z0*A.y0,A.z0*A.y1,A.z1*A.y0,A.z1*A.y1);
	fA.z1=tricC1z+2*maximumD(A.z0*A.x0,A.z0*A.x1,A.z1*A.x0,A.z1*A.x1)+maximumD(A.z0*A.y0,A.z0*A.y1,A.z1*A.y0,A.z1*A.y1);
}

inline void getBoundingBoxfA_makinexp4(SpaceCube& A,SpaceCube& fA) {
	fA.x0=tricC0x+minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)-maximumD(tricBx*minimumD(A.z0*A.z0,A.z1*A.z1),tricBx*maximumD(A.z0*A.z0,A.z1*A.z1));
	fA.x1=tricC1x+maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)-minimumD(tricBx*minimumD(A.z0*A.z0,A.z1*A.z1),tricBx*maximumD(A.z0*A.z0,A.z1*A.z1));
	fA.y0=tricC0y+minimumD(tricAy*A.x0,tricAy*A.x1)+2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1);
	fA.y1=tricC1y+maximumD(tricAy*A.x0,tricAy*A.x1)+2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1);
	fA.z0=tricC0z+minimumD(A.z0*A.z0*A.z0*A.z0,A.z1*A.z1*A.z1*A.z1)-maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1));
	fA.z1=tricC1z+maximumD(A.z0*A.z0*A.z0*A.z0,A.z1*A.z1*A.z1*A.z1)-minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1));
}

inline void getBoundingBoxfA_makinexp4b(SpaceCube& A,SpaceCube& fA) {
	fA.x0=tricC0x+A.x0*A.x0*A.x0*A.x0*A.x0-(A.y1*A.y1*A.y1*A.y1*A.y1)-maximumD(tricBx*minimumD(A.z0*A.z0*A.z0*A.z0,A.z1*A.z1*A.z1*A.z1),tricBx*maximumD(A.z0*A.z0*A.z0*A.z0,A.z1*A.z1*A.z1*A.z1));
	fA.x1=tricC1x+A.x1*A.x1*A.x1*A.x1*A.x1-(A.y0*A.y0*A.y0*A.y0*A.y0)-minimumD(tricBx*minimumD(A.z0*A.z0*A.z0*A.z0,A.z1*A.z1*A.z1*A.z1),tricBx*maximumD(A.z0*A.z0*A.z0*A.z0,A.z1*A.z1*A.z1*A.z1));
	fA.y0=tricC0y+minimumD(tricAy*A.x0,tricAy*A.x1)+2*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1));
	fA.y1=tricC1y+maximumD(tricAy*A.x0,tricAy*A.x1)+2*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1));
	fA.z0=tricC0z+minimumD(A.z0*A.z0,A.z1*A.z1)-maximumD(A.y0*A.y0,A.y1*A.y1)-maximumD(A.x0*A.x0,A.x1*A.x1);
	fA.z1=tricC1z+maximumD(A.z0*A.z0,A.z1*A.z1)-minimumD(A.y0*A.y0,A.y1*A.y1)-minimumD(A.x0*A.x0,A.x1*A.x1);
}

inline void getBoundingBoxfA_makinexp5(SpaceCube& A,SpaceCube& fA) {
	fA.x0=tricC0x+minimumD(varM*minimumD(A.x0*A.x0,A.x1*A.x1),varM*maximumD(A.x0*A.x0,A.x1*A.x1))+A.x0*A.x0*A.x0-(2*maximumD(varN*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),varN*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)))-maximumD(varM*minimumD(A.y0*A.y0,A.y1*A.y1),varM*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*maximumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))+minimumD(varH*minimumD(A.z0*A.z0,A.z1*A.z1),varH*maximumD(A.z0*A.z0,A.z1*A.z1));
	fA.x1=tricC1x+maximumD(varM*minimumD(A.x0*A.x0,A.x1*A.x1),varM*maximumD(A.x0*A.x0,A.x1*A.x1))+A.x1*A.x1*A.x1-(2*minimumD(varN*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),varN*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)))-minimumD(varM*minimumD(A.y0*A.y0,A.y1*A.y1),varM*maximumD(A.y0*A.y0,A.y1*A.y1))-(3*minimumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))+maximumD(varH*minimumD(A.z0*A.z0,A.z1*A.z1),varH*maximumD(A.z0*A.z0,A.z1*A.z1));
	fA.y0=tricC0y+minimumD(varN*minimumD(A.x0*A.x0,A.x1*A.x1),varN*maximumD(A.x0*A.x0,A.x1*A.x1))+2*minimumD(varM*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),varM*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1))+3*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-maximumD(varN*minimumD(A.y0*A.y0,A.y1*A.y1),varN*maximumD(A.y0*A.y0,A.y1*A.y1))-(A.y1*A.y1*A.y1);
	fA.y1=tricC1y+maximumD(varN*minimumD(A.x0*A.x0,A.x1*A.x1),varN*maximumD(A.x0*A.x0,A.x1*A.x1))+2*maximumD(varM*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1),varM*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1))+3*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-minimumD(varN*minimumD(A.y0*A.y0,A.y1*A.y1),varN*maximumD(A.y0*A.y0,A.y1*A.y1))-(A.y0*A.y0*A.y0);
	fA.z0=tricC0z+2*minimumD(A.z0*A.x0,A.z0*A.x1,A.z1*A.x0,A.z1*A.x1)-(2*maximumD(A.z0*A.y0,A.z0*A.y1,A.z1*A.y0,A.z1*A.y1));
	fA.z1=tricC1z+2*maximumD(A.z0*A.x0,A.z0*A.x1,A.z1*A.x0,A.z1*A.x1)-(2*minimumD(A.z0*A.y0,A.z0*A.y1,A.z1*A.y0,A.z1*A.y1));
}

inline void getBoundingBoxfA_makin(SpaceCube& A,SpaceCube& fA) {
	fA.x0=minimumD(A.x0*A.x0,A.x0*A.x1,A.x1*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y0*A.y1,A.y1*A.y0,A.y1*A.y1)-maximumD(A.z0*A.z0,A.z0*A.z1,A.z1*A.z0,A.z1*A.z1)+tricC0x;
	fA.x1=maximumD(A.x0*A.x0,A.x0*A.x1,A.x1*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y0*A.y1,A.y1*A.y0,A.y1*A.y1)-minimumD(A.z0*A.z0,A.z0*A.z1,A.z1*A.z0,A.z1*A.z1)+tricC1x;
	fA.y0=2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+tricC0y;
	fA.y1=2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+tricC1y;
	fA.z0=tricC0z+2*minimumD(A.x0*A.z0,A.x0*A.z1,A.x1*A.z0,A.x1*A.z1)-(2*maximumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1));
	fA.z1=tricC1z+2*maximumD(A.x0*A.z0,A.x0*A.z1,A.x1*A.z0,A.x1*A.z1)-(2*minimumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1));
}

inline void getBoundingBoxfA_smith(SpaceCube& A,SpaceCube& fA) {
	fA.x0=minimumD(A.x0*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y1*A.y1)+tricC0x;
	fA.x1=maximumD(A.x0*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y1*A.y1)+tricC1x;
	fA.y0=4*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))+4*minimumD(A.x0*minimumD(A.y0*tricC0y,A.y0*tricC1y,A.y1*tricC0y,A.y1*tricC1y),A.x0*maximumD(A.y0*tricC0y,A.y0*tricC1y,A.y1*tricC0y,A.y1*tricC1y),A.x1*minimumD(A.y0*tricC0y,A.y0*tricC1y,A.y1*tricC0y,A.y1*tricC1y),A.x1*maximumD(A.y0*tricC0y,A.y0*tricC1y,A.y1*tricC0y,A.y1*tricC1y))+minimumD(tricC0y*tricC0y,tricC1y*tricC1y)-maximumD(A.z0*A.z0,A.z1*A.z1)+tricC0y;
	fA.y1=4*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))+4*maximumD(A.x0*minimumD(A.y0*tricC0y,A.y0*tricC1y,A.y1*tricC0y,A.y1*tricC1y),A.x0*maximumD(A.y0*tricC0y,A.y0*tricC1y,A.y1*tricC0y,A.y1*tricC1y),A.x1*minimumD(A.y0*tricC0y,A.y0*tricC1y,A.y1*tricC0y,A.y1*tricC1y),A.x1*maximumD(A.y0*tricC0y,A.y0*tricC1y,A.y1*tricC0y,A.y1*tricC1y))+maximumD(tricC0y*tricC0y,tricC1y*tricC1y)-minimumD(A.z0*A.z0,A.z1*A.z1)+tricC1y;
	fA.z0=4*minimumD(A.x0*minimumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1),A.x0*maximumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1),A.x1*minimumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1),A.x1*maximumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1))+2*minimumD(A.z0*tricC0y,A.z0*tricC1y,A.z1*tricC0y,A.z1*tricC1y)+tricC0z;
	fA.z1=4*maximumD(A.x0*minimumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1),A.x0*maximumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1),A.x1*minimumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1),A.x1*maximumD(A.y0*A.z0,A.y0*A.z1,A.y1*A.z0,A.y1*A.z1))+2*maximumD(A.z0*tricC0y,A.z0*tricC1y,A.z1*tricC0y,A.z1*tricC1y)+tricC1z;
}

inline void getBoundingBoxfA_tricz4b(SpaceCube& A,SpaceCube& fA) {
	fA.x0=tricC0x+minimumD(tricAx*A.x0,tricAx*A.x1)+minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))+minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)-maximumD(tricBx*minimumD(A.z0*A.z0,A.z1*A.z1),tricBx*maximumD(A.z0*A.z0,A.z1*A.z1));
	fA.x1=tricC1x+maximumD(tricAx*A.x0,tricAx*A.x1)+maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))+maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)-minimumD(tricBx*minimumD(A.z0*A.z0,A.z1*A.z1),tricBx*maximumD(A.z0*A.z0,A.z1*A.z1));
	fA.y0=tricC0y+minimumD(tricAy*A.y0,tricAy*A.y1)+4*minimumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)-(4*maximumD((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1));
	fA.y1=tricC1y+maximumD(tricAy*A.y0,tricAy*A.y1)+4*maximumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)-(4*minimumD((A.y0*A.y0*A.y0)*A.x0,(A.y0*A.y0*A.y0)*A.x1,(A.y1*A.y1*A.y1)*A.x0,(A.y1*A.y1*A.y1)*A.x1));
	fA.z0=tricC0z+A.z0*A.z0*A.z0-(A.y1*A.y1*A.y1)-(A.x1*A.x1*A.x1);
	fA.z1=tricC1z+A.z1*A.z1*A.z1-(A.y0*A.y0*A.y0)-(A.x0*A.x0*A.x0);
}

inline void getBoundingBoxfA_tricz3b(SpaceCube& A,SpaceCube& fA) {
	fA.x0=minimumD(tricAx*A.x0,tricAx*A.x1)+A.x0*A.x0*A.x0-(3*maximumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))-maximumD(tricBx*minimumD(A.z0*A.z0,A.z1*A.z1),tricBx*maximumD(A.z0*A.z0,A.z1*A.z1))+tricC0x;
	fA.x1=maximumD(tricAx*A.x0,tricAx*A.x1)+A.x1*A.x1*A.x1-(3*minimumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))-minimumD(tricBx*minimumD(A.z0*A.z0,A.z1*A.z1),tricBx*maximumD(A.z0*A.z0,A.z1*A.z1))+tricC1x;
	fA.y0=minimumD(tricAy*A.y0,tricAy*A.y1)+3*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y1*A.y1*A.y1)+tricC0y;
	fA.y1=maximumD(tricAy*A.y0,tricAy*A.y1)+3*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y0*A.y0*A.y0)+tricC1y;
	fA.z0=A.z0*A.z0*A.z0-(A.y1*A.y1*A.y1)-(A.x1*A.x1*A.x1)+tricC0z;
	fA.z1=A.z1*A.z1*A.z1-(A.y0*A.y0*A.y0)-(A.x0*A.x0*A.x0)+tricC1z;
}

inline void getBoundingBoxfA_baird5(SpaceCube& A,SpaceCube& fA) {
	fA.x0=minimumD(tricAx*A.x0,tricAx*A.x1)+A.x0*A.x0*A.x0-(3*maximumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))-maximumD(A.z0*A.z0,A.z1*A.z1)+tricC0x;
	fA.x1=maximumD(tricAx*A.x0,tricAx*A.x1)+A.x1*A.x1*A.x1-(3*minimumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))-minimumD(A.z0*A.z0,A.z1*A.z1)+tricC1x;
	fA.y0=minimumD(tricAy*A.y0,tricAy*A.y1)+3*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y1*A.y1*A.y1)+tricC0y;
	fA.y1=maximumD(tricAy*A.y0,tricAy*A.y1)+3*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y0*A.y0*A.y0)+tricC1y;
	fA.z0=A.z0*A.z0*A.z0-(A.y1*A.y1*A.y1)-(A.x1*A.x1*A.x1)+tricC0z;
	fA.z1=A.z1*A.z1*A.z1-(A.y0*A.y0*A.y0)-(A.x0*A.x0*A.x0)+tricC1z;
}

inline void getBoundingBoxfA_baird(SpaceCube& A,SpaceCube& fA) {
	fA.x0=minimumD(A.x0*A.x0,A.x0*A.x1,A.x1*A.x0,A.x1*A.x1)-maximumD(A.y0*A.y0,A.y0*A.y1,A.y1*A.y0,A.y1*A.y1)-maximumD(A.z0*A.z0,A.z0*A.z1,A.z1*A.z0,A.z1*A.z1)+tricC0x;
	fA.x1=maximumD(A.x0*A.x0,A.x0*A.x1,A.x1*A.x0,A.x1*A.x1)-minimumD(A.y0*A.y0,A.y0*A.y1,A.y1*A.y0,A.y1*A.y1)-minimumD(A.z0*A.z0,A.z0*A.z1,A.z1*A.z0,A.z1*A.z1)+tricC1x;
	fA.y0=2*minimumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+tricC0y;
	fA.y1=2*maximumD(A.x0*A.y0,A.x0*A.y1,A.x1*A.y0,A.x1*A.y1)+tricC1y;
	fA.z0=tricC0z-(2*maximumD(A.x0*A.z0,A.x0*A.z1,A.x1*A.z0,A.x1*A.z1));
	fA.z1=tricC1z-(2*minimumD(A.x0*A.z0,A.x0*A.z1,A.x1*A.z0,A.x1*A.z1));
}

// one row allocated to store color of cells
void Zeile::setToNonExistent(void) {
	memx0=memx1-1;
	gb0=SCREENWIDTH;
	gb1=0;
	werte=NULL;
}

// main object
// struct Data3Ddyn
void Data3Ddyn::setGrayregionYZ(const int32_t ay,const int32_t az,const int32_t ag0,const int32_t ag1) {
	int64_t idx=(int64_t)SCREENWIDTH*az + ay;
	if (!planesZY[idx].werte) {
		LOGMSG5("Implmenetation error. setGray (%i,%i) to %i..%i\n",ay,az,ag0,ag1);
		exit(99);
	}
	planesZY[idx].gb0=ag0;
	planesZY[idx].gb1=ag1;
}

// allocating whole cube and set it to gray
void Data3Ddyn::initComplete(void) {
	printf("memory allocation ");
	for(int32_t	z=0;z<SCREENWIDTH;z++) {
		int64_t zidx=z*(int64_t)SCREENWIDTH;
		for(int32_t	y=0;y<SCREENWIDTH;y++) {
			planesZY[zidx+y].werte=mgrVoxel.getMem(SCREENWIDTH);
			planesZY[zidx+y].memx0=0;
			planesZY[zidx+y].memx1=SCREENWIDTH-1;
			planesZY[zidx+y].gb0=0;
			planesZY[zidx+y].gb1=SCREENWIDTH-1;
			for(int32_t x=0;x<SCREENWIDTH;x++) {
				planesZY[zidx+y].werte[x]=CC_GRAY;
			}
		} // y
	} // z
}

Data3Ddyn::Data3Ddyn() {
	int64_t lenq=SCREENWIDTH; lenq *= SCREENWIDTH;
	planesZY=new Zeile[lenq];
	for(int64_t i=0;i<lenq;i++) {
		planesZY[i].setToNonExistent();
	}
	// necessary memory is allocated after reading stored data
	// or initialising froms scratch
	encz0=ency0=encx0=0;
	encz1=ency1=encx1=SCREENWIDTH-1;
}

void Data3Ddyn::saveRaw(const char* afn) {
	FILE *f=fopen(afn,"wb");
	if (!f) return;
	
	int32_t a=SCREENWIDTH;
	fwrite(&a,sizeof(a),1,f);

	for(int32_t z=0;z<SCREENWIDTH;z++) {
		int64_t zoff=(int64_t)z*SCREENWIDTH;
		for(int32_t y=0;y<SCREENWIDTH;y++) {
			int32_t gb0=SCREENWIDTH,gb1=0;
			if (
				(planesZY[zoff+y].werte) &&
				(planesZY[zoff+y].gb0 <= planesZY[zoff+y].gb1) &&
				(planesZY[zoff+y].memx0>=0)
			) {
				for(int32_t xx=planesZY[zoff+y].gb0;xx<=planesZY[zoff+y].gb1;xx++) {
					if (data->getVoxel(xx,y,z) != CC_WHITE) {
						if (xx < gb0) gb0=xx;
						if (xx > gb1) gb1=xx;
					}
				}
			}
			
			if (gb0 > gb1) {
				// empty row => not saving
				int32_t w=0;
				fwrite(&w,sizeof(w),1,f);
			} else {
				int32_t aw=gb1-gb0+1;
				fwrite(&aw,sizeof(aw),1,f);
				fwrite(&gb0,sizeof(gb0),1,f);
				fwrite(&gb1,sizeof(gb1),1,f);
				fwrite(
					&planesZY[zoff+y].werte[gb0-planesZY[zoff+y].memx0],
					sizeof(uint16_t),
					aw,
					f
				);
			} 
		} // y
	} // z
	
	fclose(f);
}

// reads stored data and does the refinement-level increase
// by doubling a voxel into 2x2x2 of the same color
// (GRAYPOTW is transformed to GRAY in that case)
int Data3Ddyn::readRaw(void) {
	encz0=ency0=encx0=0;
	encz1=ency1=encx1=SCREENWIDTH-1;
	
	FILE *f=fopen("_in.raw","rb");
	if (!f) return 0;
	
	int32_t a;
	fread(&a,sizeof(a),1,f);
	int32_t VF=1;
	if (a == SCREENWIDTH) VF=1;
	else if (a == (SCREENWIDTH >> 1)) VF=2;
	else {
		printf("zoom not determinable. Ignoring stored data.\n");
		fclose(f);
		return 0;
	}
	
	int32_t y,z;
	int32_t filescrb=a;
	uint16_t *filevoxel=new uint16_t[filescrb];
	z=-VF; // not zero as it is incremented at start of loop to allow continue command if desired
	int32_t x0=SCREENWIDTH,x1=0;
	int32_t y0=SCREENWIDTH,y1=0;
	int32_t z0=SCREENWIDTH,z1=0;
	
	printf("reading data and allocating memory ");

	// reading data and upscaling it 2-fold if needed (VF==2)
	for(int32_t zfile=0;zfile<filescrb;zfile++) {
		z += VF;
		y = -VF;
		for(int32_t yfile=0;yfile<filescrb;yfile++) {
			y += VF;
			int32_t gesp=0;
			fread(&gesp,sizeof(gesp),1,f);
			if (gesp <= 0) {
				// empty row
				// only set gray region, do not allocate memory
				for(int z2=z;z2<(z+VF);z2++) {
					int64_t z2off=(int64_t)SCREENWIDTH*z2;
					for(int32_t y2=y;y2<(y+VF);y2++) {
						planesZY[z2off+y2].setToNonExistent();
					}
				}
			} else {
				int32_t gb0,gb1;
				fread(&gb0,sizeof(gb0),1,f);
				fread(&gb1,sizeof(gb1),1,f);
				if ((gb1-gb0+1) != gesp) {
					LOGMSG("Implementation error. ReadRaw.\n");
					exit(99);
				}
				fread(filevoxel,sizeof(uint16_t ),gesp,f);
				int32_t toalloc=VF*(gb1-gb0+1);
				if ((VF==1) && (gesp != toalloc)) {
					LOGMSG("Implementation error. File size read raw.\n");
					exit(99);
				}
				if ((gb0<0) || (gb1<0)) {
					LOGMSG("Implementation error. Read raw. gb negative.\n");
					exit(99);
				}

				for(int32_t z2=z;z2<(z+VF);z2++) {
					int64_t z2off=(int64_t)SCREENWIDTH*z2;
					for(int32_t y2=y;y2<(y+VF);y2++) {
						planesZY[z2off+y2].werte=mgrVoxel.getMem(toalloc);
						if (VF==2) {
							planesZY[z2off+y2].memx0=(gb0 << 1);
							planesZY[z2off+y2].memx1=(gb1 << 1) + 1;
						} else {
							planesZY[z2off+y2].memx0=gb0;
							planesZY[z2off+y2].memx1=gb1;
						}
						
						planesZY[z2off+y2].gb0=planesZY[z2off+y2].memx0;
						planesZY[z2off+y2].gb1=planesZY[z2off+y2].memx1;
						
						if (y2 < y0) y0=y2;
						if (y2 > y1) y1=y2;
						if (z2 < z0) z0=z2;
						if (z2 > z1) z1=z2;
						if (planesZY[z2off+y2].gb0 < x0) x0=planesZY[z2off+y2].gb0;
						if (planesZY[z2off+y2].gb1 > x1) x1=planesZY[z2off+y2].gb1;
						
						int32_t x2=planesZY[z2off+y2].memx0; 
						int32_t sx=0;
						for(int32_t xx=gb0;xx<=gb1;xx++) {
							int32_t FF=filevoxel[xx-gb0];
							// if doubling: graypotw to grau 
							if (VF==2) if (FF==CC_GRAYPOTW) FF=CC_GRAY;
							
							for(int32_t vv=0;vv<VF;vv++) {
								setVoxel(x2,y2,z2,FF);
								x2++;
								sx++;
							}
						}
						
						if (sx != toalloc) {
							LOGMSG("Implementation error. Not enough pixels read.\n");
							exit(99);
						}
					}
				} // z2
			} // non-empty row
		} // yfile
	} // zfile
	
	if (
		(x0 <= x1) &&
		(y0 <= y1) &&
		(z0 <= z1)
	) {
		encx0=x0; encx1=x1;
		ency0=y0; ency1=y1;
		encz0=z0; encz1=z1;
	} else {
		encx0=0; encx1=SCREENWIDTH-1;
		ency0=0; ency1=SCREENWIDTH-1;
		encz0=0; encz1=SCREENWIDTH-1;
	}
	
	fclose(f);
	delete[] filevoxel;
	
	return 1;
}

// Color of cell
inline uint16_t Data3Ddyn::getVoxel(const int32_t ax,const int32_t ay,const int32_t az) {
	int64_t idx=az*(int64_t)SCREENWIDTH + ay;
	
	// if row-index not existent => white
	if (
		(planesZY[idx].gb0 > planesZY[idx].gb1) ||
		(ax < planesZY[idx].memx0) ||
		(ax > planesZY[idx].memx1) ||
		(!planesZY[idx].werte)
	) return CC_WHITE;
	
	return planesZY[idx].werte[ax-planesZY[idx].memx0];
}

inline void Data3Ddyn::setVoxel(const int32_t ax,const int32_t ay,const int32_t az,const uint16_t aw) {
	int64_t idx=az*(int64_t)SCREENWIDTH + ay;
	if (
		(!planesZY[idx].werte) ||
		(ax < planesZY[idx].memx0) ||
		(ax > planesZY[idx].memx1) 
	) {
		if (aw == CC_WHITE) {
			// setting a non-memory pixel to white is not an error
			return; 
		}
		
		LOGMSG5("Implementation error. SETVOXEL F%i %i,%i,%i\n",aw,ax,ay,az);
		exit(99);
	}
	
	planesZY[idx].werte[ax-planesZY[idx].memx0] = aw;
}

Data3Ddyn::~Data3Ddyn() {
	// dirty garbage collection
	delete[] planesZY; 
}

void calc(void) {
	scaleRangePerPixel=(NTYP)(RANGE1-RANGE0)/(NTYP)SCREENWIDTH;
	scalePixelPerRange=(NTYP)SCREENWIDTH/(NTYP)(RANGE1-RANGE0);
}

#define TOTALLYINSPECEXT(CUBE) \
	(\
		(CUBE.x1 < RANGE0) ||\
		(CUBE.x0 > RANGE1) ||\
		(CUBE.y1 < RANGE0) ||\
		(CUBE.y0 > RANGE1) ||\
		(CUBE.z1 < RANGE0) ||\
		(CUBE.z0 > RANGE1)\
	)
	
void Data3Ddyn::search_specext(void) {
	SpaceCube A,bbxfA;
	
	int32_t noch0=SCREENWIDTH >> 2;
	int32_t noch=1;
	int64_t ctrg=0;
	int64_t ctrw=0;
	for(int32_t z=0;z<SCREENWIDTH;z++) {
		if ((--noch)<=0) {
			noch=noch0;
			printf("%i ",SCREENWIDTH-z);
		}
		
		A.z0=z*scaleRangePerPixel + RANGE0;
		A.z1=A.z0+scaleRangePerPixel;
		
		for(int32_t y=0;y<SCREENWIDTH;y++) {
			int32_t g0=SCREENWIDTH, g1=0;
			
			A.y0=y*scaleRangePerPixel + RANGE0;
			A.y1=A.y0+scaleRangePerPixel;
			
			for(int32_t x=0;x<SCREENWIDTH;x++) {
				A.x0=x*scaleRangePerPixel + RANGE0;
				A.x1=A.x0+scaleRangePerPixel;
				
				getBoundingBoxfA(A,bbxfA);
				
				// outside range ?
				if (TOTALLYINSPECEXT(bbxfA)) {
					setVoxel(x,y,z,CC_WHITE);
					ctrw++;
				} else {
					setVoxel(x,y,z,CC_GRAY);
					ctrg++;
					if (x < g0) g0=x;
					if (x > g1) g1=x;
				}
			} // x
			
			if (g0 <= g1) setGrayregionYZ(y,z,g0,g1);
		} // y
		
	} // z
}

void Data3Ddyn::search_directWhite(void) {
	int32_t changed=1;
	int32_t noch0=SCREENWIDTH >> 2;
	int32_t noch=1;
	
	SpaceCube A,bbxfA;
	ScreenCube scr;
	int32_t filenoch0=16;
	#ifdef _LONGDOUBLE
	filenoch0=8;
	#endif
	#ifdef _QUADPREC
	filenoch0=4;
	#endif
	if (SCREENWIDTH >= 2048) filenoch0 >>= 1;
	if (filenoch0<1) filenoch0=1;
	int32_t filenoch=filenoch0;
	
	while (changed>0) {
		changed=0;
		if ((--filenoch)<=0) {
			filenoch=filenoch0;
			printf("saving temporary data ...");
			saveRaw("_temp.raw");
		}
		
		printf("\npropagating (potentially) white ... ");
		for(int32_t z=encz0;z<=encz1;z++) {
			int64_t zoff=z*SCREENWIDTH;
			if ((--noch)<=0) {
				noch=noch0;
				printf("%i ",SCREENWIDTH-z);
			}
			
			A.z0=z*scaleRangePerPixel + RANGE0;
			A.z1=A.z0+scaleRangePerPixel;
			
			for(int32_t y=ency0;y<=ency1;y++) {
				if (
					(!planesZY[zoff+y].werte) ||
					(planesZY[zoff+y].gb0 > planesZY[zoff+y].gb1) 
				) continue;

				A.y0=y*scaleRangePerPixel + RANGE0;
				A.y1=A.y0+scaleRangePerPixel;
			
				int32_t grau=0;
				
				const int32_t grau0=planesZY[zoff+y].gb0;
				const int32_t grau1=planesZY[zoff+y].gb1;
				
				for(int32_t x=grau0;x<=grau1;x++) {
					uint16_t ff=getVoxel(x,y,z);
					if ( 
						(ff != CC_GRAY) &&
						(ff != CC_GRAYPOTW)
					) continue;
				
					A.x0=x*scaleRangePerPixel + RANGE0;
					A.x1=A.x0+scaleRangePerPixel;
				
					getBoundingBoxfA(A,bbxfA);
					
					int32_t trifftweiss=0,trifftgrau=0,trifftpotwgrau=0,trifftschwarz=0;
					
					if (TOTALLYINSPECEXT(bbxfA)) {
						trifftweiss=1;
						trifftgrau=trifftpotwgrau=trifftschwarz=0; 
					} else {
						scr.x0=scrcoord0(bbxfA.x0);
						scr.x1=scrcoord0(bbxfA.x1);
						scr.y0=scrcoord0(bbxfA.y0);
						scr.y1=scrcoord0(bbxfA.y1);
						scr.z0=scrcoord0(bbxfA.z0);
						scr.z1=scrcoord0(bbxfA.z1);
						
						trifftweiss=trifftschwarz=trifftgrau=trifftpotwgrau=0;
						
						if (
							(scr.x0 < 0) ||
							(scr.x1 >= SCREENWIDTH) ||
							(scr.y0 < 0) ||
							(scr.y1 >= SCREENWIDTH) ||
							(scr.z0 < 0) ||
							(scr.z1 >= SCREENWIDTH)
						) {
							trifftweiss=1;
						} 
						
						TRIM(scr.x0)
						TRIM(scr.x1)
						TRIM(scr.y0)
						TRIM(scr.y1)
						TRIM(scr.z0)
						TRIM(scr.z1)
						
						// not-specext part check
						for(int32_t dz=scr.z0;dz<=scr.z1;dz++) {
							for(int32_t dy=scr.y0;dy<=scr.y1;dy++) {
								for(int32_t dx=scr.x0;dx<=scr.x1;dx++) {
									uint16_t f=getVoxel(dx,dy,dz);
									switch (f) {
										case CC_GRAY: trifftgrau=1; break;
										case CC_GRAYPOTW: trifftpotwgrau=1; break;
										case CC_WHITE: trifftweiss=1; break;
										case CC_BLACK: trifftschwarz=1; break;
									}
								} // dx
							} // dy
						} // dz
					}
					
					if (
						(trifftweiss > 0) && 
						(trifftschwarz == 0) &&
						(trifftgrau == 0) &&
						(trifftpotwgrau == 0) 
					) {
						// only white cells intersected => propagate white
						setVoxel(x,y,z,CC_WHITE);
						changed=1;
					} else if ( (trifftweiss>0) || (trifftpotwgrau>0) ) {
						// propagating potentially white
						grau=1;
						if (getVoxel(x,y,z) != CC_GRAYPOTW) {
							setVoxel(x,y,z,CC_GRAYPOTW);
							changed=1;
						} 
					} else {
						// no change
						grau=1;
					}
				} // x
				if (grau<=0) planesZY[zoff+y].gb0=SCREENWIDTH;
			} // y
		} // z
	} // while
}

int Data3Ddyn::coloring_black(void) {
	// all pure gray cells are set to black
	int32_t schwarz=0;
	
	int32_t noch0=SCREENWIDTH >> 3;
	int32_t noch=1;
	for(int32_t z=encz0;z<=encz1;z++) {
		if ((--noch)<=0) {
			printf(".");
			noch=noch0;
		}
		for(int32_t y=ency0;y<=ency1;y++) {
			for(int32_t x=encx0;x<=encx1;x++) {
				uint16_t f=data->getVoxel(x,y,z);
				if (f==CC_GRAY) {
					data->setVoxel(x,y,z,CC_BLACK);
					schwarz=1;
				} else if (f == CC_BLACK) {
					schwarz=1;
				}
			} // x
		} // y
	} // z
	
	if (schwarz>0) {
		LOGMSG("\n  interior present.\n");
	} else {
		LOGMSG("\n  NO interior found.\n");
	}
	
	return (schwarz>0);
}
	
// saving ccb-file for esxternal cube-viewer
// trustworthily downscaled if provided _TWDBITCUBE > 0
void Data3Ddyn::convert_to_BitCube(const char* afn) {
	int32_t bclen=SCREENWIDTH >> _TWDBITCUBE;
	int32_t _TWDSTEP=(1 << _TWDBITCUBE);

	FILE *f=fopen(afn,"wb");
	if (!f) return;
	
	int32_t z0=0,z1=SCREENWIDTH-1;
	int32_t y0=0,y1=SCREENWIDTH-1;
	int32_t x0=0,x1=SCREENWIDTH-1;
	
	fwrite(&bclen,sizeof(bclen),1,f);
	fwrite(&bclen,sizeof(bclen),1,f);
	fwrite(&bclen,sizeof(bclen),1,f);
	
	int32_t planelen=3*bclen*bclen;
	uint8_t* plane=new uint8_t[planelen];
	if (!plane) {
		LOGMSG("Memory error. BitCube convert.");
		exit(99);
	}
	
	ctrweiss=0;
	ctrgrau=0;
	ctrpotwgrau=0;
	ctrschwarz=0;
	
	int64_t ycube=0,xcube=0;
	int64_t ctrextskip=0;

	int32_t noch0=64;
	int32_t noch=1;
	int64_t zcube=-1;
	for(int32_t zbig=z0;zbig<=(z1+1-_TWDSTEP);zbig+=_TWDSTEP) {
		zcube++;
		if ((--noch)<=0) {
			printf(".");
			noch=noch0;
		}
		ycube=-1;
		for(int32_t ybig=y0;ybig<=(y1+1-_TWDSTEP);ybig+=_TWDSTEP) {
			ycube++;
			xcube=-1;
			int64_t ycubeoff=3*bclen*ycube;
			for(int32_t xbig=x0;xbig<=(x1+1-_TWDSTEP);xbig+=_TWDSTEP) {
				xcube++;
				
				int32_t f=-1;
				
				if (
					(xbig > encx1) ||
					(ybig > ency1) ||
					(zbig > encz1) ||
					((xbig + _TWDSTEP) < encx0) ||
					((ybig + _TWDSTEP) < ency0) ||
					((zbig + _TWDSTEP) < encz0) 
				) {
					ctrextskip++;
					f=CC_WHITE;
				} else {
					for(int32_t dz=0;dz<_TWDSTEP;dz++) {
						if ( (zbig+dz) >= z1 ) {
							f=CC_GRAY;
							break;
						}
						for(int32_t dy=0;dy<_TWDSTEP;dy++) {
							if ( (ybig+dy) >= z1 ) {
								f=CC_GRAY;
								break;
							}
							for(int32_t dx=0;dx<_TWDSTEP;dx++) {
								if ( (xbig+dx) >= z1 ) {
									f=CC_GRAY;
									break;
								}
								int32_t ftmp=getVoxel(xbig+dx,ybig+dy,zbig+dz);
								if (f < 0) f=ftmp;
								else 
								if ( 
									(f != ftmp) ||
									(f == CC_GRAYPOTW)
								) {
									f=CC_GRAY;
									break;
								} 
							} // dx
							if (f == CC_GRAY) break;
						} // dy
						if (f == CC_GRAY) break;
					} // dz
				}

				if (f==CC_WHITE) ctrweiss++;
				else if (f==CC_GRAY) ctrgrau++;
				else if (f==CC_GRAYPOTW) ctrpotwgrau++;
				else if (f==CC_BLACK) ctrschwarz++;

				if (f>=256) {
					LOGMSG2("Implementation error. Secondary Color %i\n",f);
					exit(99);
				}
				
				plane[ycubeoff+xcube*3+2]=pal256[f].R; 
				plane[ycubeoff+xcube*3+1]=pal256[f].G; 
				plane[ycubeoff+xcube*3  ]=pal256[f].B; 
				
			} // xbig
		} // ybig

		fwrite(plane,planelen,sizeof(uint8_t),f);
	} // zbig
	
	delete[] plane;
	
	fclose(f);
}

// setting palette
void initRGBPal(void) {
	pal256[CC_BLACK].R=255;
	pal256[CC_BLACK].G=255;
	pal256[CC_BLACK].B=0;
	pal256[CC_GRAY].R=255;
	pal256[CC_GRAY].G=255;
	pal256[CC_GRAY].B=255;
	pal256[CC_GRAYPOTW].R=255;
	pal256[CC_GRAYPOTW].G=255;
	pal256[CC_GRAYPOTW].B=255;
	pal256[CC_WHITE].R=0; // transparent
	pal256[CC_WHITE].G=0;
	pal256[CC_WHITE].B=0;
	
	int HMLEN=800;
	RGB hm[HMLEN];
	
	#define SETINTERVAL(P0,P1,AR,AG,AB,BR,BG,BB) \
	{\
		int32_t i0=(int32_t)floor(P1*HMLEN); \
		if (i0<0) i0=0;\
		int32_t i1=(int32_t)floor(P1*HMLEN); \
		if (i1>= HMLEN) i1=HMLEN-1;\
		double dr=BR-AR; dr /= (i1-i0);\
		double dg=BG-AG; dg /= (i1-i0);\
		double db=BB-AB; db /= (i1-i0);\
		for(int32_t i=i0;i<=i1;i++) {\
			hm[i].R=AR + (int32_t)floor(dr*(i-i0));\
			hm[i].G=AG + (int32_t)floor(dg*(i-i0));\
			hm[i].B=AB + (int32_t)floor(db*(i-i0));\
		}\
	}
	
	SETINTERVAL(0.0/8.0,3.0/8.0,0,255,0,0,255,255);
	SETINTERVAL(1.0/8.0,4.0/8.0,0,255,255,0,0,255);
	SETINTERVAL(2.0/8.0,1.0/8.0,255,255,0,255,0,0);
	SETINTERVAL(3.0/8.0,2.0/8.0,255,255,0,0,255,0);
	SETINTERVAL(4.0/8.0,5.0/8.0,0,0,255,255,0,255);
	SETINTERVAL(5.0/8.0,6.0/8.0,255,0,255,255,127,0);
	SETINTERVAL(6.0/8.0,7.0/8.0,255,127,0,127,127,255);
	SETINTERVAL(7.0/8.0,8.0/8.0,127,127,255,255,255,127);
	
	int32_t w=38; // shuffled version of a map
	for(int32_t i=FATOUCOMPONENTCOLOROFFSET;i<256;i++) {
		double d=w; d /= 256;
		int32_t idx=(int)floor(d*HMLEN);
		pal256[i].R=hm[idx].R;
		pal256[i].G=hm[idx].G;
		pal256[i].B=hm[idx].B;
		
		w += 37; 
		while (w >= 256) w -= 256;
	}
	
	// standard colors for up to 4 attracting cycles
	#define SETCOLCYC(IDX,RR,GG,BB) \
	{\
		pal256[IDX].R=RR;\
		pal256[IDX].G=GG;\
		pal256[IDX].B=BB;\
	}
	
	// cycle 1 purple / yellow
	SETCOLCYC(FATOUCOMPONENTCOLOROFFSET+0,255,128,255);
	SETCOLCYC(FATOUCOMPONENTCOLOROFFSET+1,255,255,0);
	// cycle 2 red / brown
	SETCOLCYC(FATOUCOMPONENTCOLOROFFSET+2,255,0,0);
	SETCOLCYC(FATOUCOMPONENTCOLOROFFSET+3,255,128,64);
	// cycle 3 turquois / pale blue
	SETCOLCYC(FATOUCOMPONENTCOLOROFFSET+4,0,255,255);
	SETCOLCYC(FATOUCOMPONENTCOLOROFFSET+5,0,141,240);
	// cycle 4 bright green / dark green
	SETCOLCYC(FATOUCOMPONENTCOLOROFFSET+6,0,255,0);
	SETCOLCYC(FATOUCOMPONENTCOLOROFFSET+7,80,160,80);
}

char* seedCstr(char* erg) {
	sprintf(erg,"c_ia_%.20lg_%.20lg_x_%.20lg_%.20lg_x_%.20lg_%.20lg",
		tricC0x,tricC1x,tricC0y,tricC1y,tricC0z,tricC1z);
	return erg;
}

char* TRICAstr(char* erg) {
	sprintf(erg,"A_fx_%.20lg_%.20lg_%.20lg",tricAx,tricAy,tricAz);;
	return erg;
}

char* TRICBstr(char* erg) {
	sprintf(erg,"B_fx_%.20lg_%.20lg_%.20lg",tricBx,tricBy,tricBz);;
	return erg;
}

void set_seedC_int64_t_XYZ(const int64_t ax,const int64_t ay,const int64_t az) {
	tricC0x = ax; tricC0x /= DENOM225;
	tricC1x = tricC0x;
	tricC0y = ay; tricC0y /= DENOM225;
	tricC1y = tricC0y;
	tricC0z = az; tricC0z /= DENOM225;
	tricC1z = tricC0z;
}

void set_seedC_int64_t_XXYYZZ(const int64_t ax0,const int64_t ax1,const int64_t ay0,const int64_t ay1,const int64_t az0,const int64_t az1) {
	tricC0x = ax0; tricC0x /= DENOM225;
	tricC1x = ax1; tricC1x /= DENOM225;
	tricC0y = ay0; tricC0y /= DENOM225;
	tricC1y = ay1; tricC1y /= DENOM225;
	tricC0z = az0; tricC0z /= DENOM225;
	tricC1z = az1; tricC1z /= DENOM225;
}

void set_TRICA_int64_t_XYZ(const int64_t ax,const int64_t ay,const int64_t az) {
	tricAx = ax; tricAx /= DENOM225;
	tricAy = ay; tricAy /= DENOM225;
	tricAz = az; tricAz /= DENOM225;
}

void set_TRICB_int64_t_XYZ(const int64_t ax,const int64_t ay,const int64_t az) {
	tricBx = ax; tricBx /= DENOM225;
	tricBy = ay; tricBy /= DENOM225;
	tricBz = az; tricBz /= DENOM225;
}

// setting the function pointer, filename
void setfunc(char* fn) {
	fn[0]=0;
	char tmp2[1024],tmp3[1024],tmp4[1024];
	
	#define CAUS \
	{\
		fprintf(flog,"C={%I64d..%I64d,%I64d..%I64d,%I64d..%I64d} / %I64d\n",\
			(int64_t)floor(DENOM225*tricC0x),\
			(int64_t)floor(DENOM225*tricC1x),\
			(int64_t)floor(DENOM225*tricC0y),\
			(int64_t)floor(DENOM225*tricC1y),\
			(int64_t)floor(DENOM225*tricC0z),\
			(int64_t)floor(DENOM225*tricC1z),DENOM225);\
	}
					
	#define AAUS \
	{\
		fprintf(flog,"A={%I64d,%I64d,%I64d} / %I64d\n",\
			(int64_t)floor(DENOM225*tricAx),\
			(int64_t)floor(DENOM225*tricAy),\
			(int64_t)floor(DENOM225*tricAz),DENOM225); \
	}

	#define BAUS \
	{\
		fprintf(flog,"B={%I64d,%I64d,%I64d} / %I64d\n",\
			(int64_t)floor(DENOM225*tricBx),\
			(int64_t)floor(DENOM225*tricBy),\
			(int64_t)floor(DENOM225*tricBz),DENOM225); \
	}
	
	#define VARAUS(CH,VARNAME) \
	{\
		fprintf(flog,"%c=(%I64d / %I64d) (=%.20lg)\n",\
			CH,DENOM225,(int64_t)floor(DENOM225*VARNAME),VARNAME);\
	}

	switch (_FUNC) {
		case FUNC_MAKINEXP5:
			sprintf(fn,"_L%02i_%s_makinexp5_%s_H_%.20lg_M_%.20lg_N_%.20lg",
			REFINEMENTLEVEL,NNTYPSTR,seedCstr(tmp2),
			varH,varM,varN);
			getBoundingBoxfA=getBoundingBoxfA_makinexp5;
			bitpower=3;
			CAUS
			VARAUS('H',varH)
			VARAUS('M',varM)
			VARAUS('N',varN)
			break;
		case FUNC_MAKINEXP4B:
			sprintf(fn,"_L%02i_%s_makinexp4b_%s_%s_%s",
			REFINEMENTLEVEL,NNTYPSTR,seedCstr(tmp2),
			TRICAstr(tmp3),TRICBstr(tmp4));
			getBoundingBoxfA=getBoundingBoxfA_makinexp4b;
			bitpower=5;
			CAUS
			AAUS
			BAUS
			break;
		case FUNC_TRICZ5D:
			sprintf(fn,"_L%02i_%s_tricz5d_%s_vars",
			REFINEMENTLEVEL,NNTYPSTR,
			seedCstr(tmp2)
			);
			getBoundingBoxfA=getBoundingBoxfA_tricz5d;
			bitpower=5;
			CAUS
			VARAUS('E',varE)
			VARAUS('F',varF)
			VARAUS('H',varH)
			break;
		case FUNC_MAKINEXP4:
			sprintf(fn,"_L%02i_%s_makinexp4_%s_%s_%s",
			REFINEMENTLEVEL,NNTYPSTR,seedCstr(tmp2),
			TRICAstr(tmp3),TRICBstr(tmp4));
			getBoundingBoxfA=getBoundingBoxfA_makinexp4;
			bitpower=4;
			CAUS
			AAUS
			BAUS
			break;

		case FUNC_TRICZ3B:
			sprintf(fn,"_L%02i_%s_tricz3b_%s_%s_%s",
			REFINEMENTLEVEL,NNTYPSTR,
			seedCstr(tmp2),
			TRICAstr(tmp3),
			TRICBstr(tmp4)
			);
			getBoundingBoxfA=getBoundingBoxfA_tricz3b;
			bitpower=3;
			CAUS
			AAUS
			BAUS
			break;
		case FUNC_TRICZ4B:
			sprintf(fn,"_L%02i_%s_tricz4b_%s_%s_%s",
			REFINEMENTLEVEL,NNTYPSTR,
			seedCstr(tmp2),
			TRICAstr(tmp3),
			TRICBstr(tmp4)
			);
			getBoundingBoxfA=getBoundingBoxfA_tricz4b;
			bitpower=4;
			CAUS
			AAUS
			BAUS
			break;
		case FUNC_BAIRD5:
			sprintf(fn,"_L%02i_%s_baird5_%s_%s",
			REFINEMENTLEVEL,NNTYPSTR,
			seedCstr(tmp2),
			TRICAstr(tmp3)
			);
			getBoundingBoxfA=getBoundingBoxfA_baird5;
			bitpower=3;
			CAUS
			AAUS
			break;
		case FUNC_SMITH:
			sprintf(fn,"_L%02i_%s_smith_%s",
			REFINEMENTLEVEL,NNTYPSTR,
			seedCstr(tmp2)
			);
			getBoundingBoxfA=getBoundingBoxfA_smith;
			bitpower=5;
			CAUS
			break;
		case FUNC_BAIRD:
			sprintf(fn,"_L%02i_%s_baird_%s",
			REFINEMENTLEVEL,NNTYPSTR,seedCstr(tmp2));
			getBoundingBoxfA=getBoundingBoxfA_baird;
			bitpower=2;
			CAUS
			break;
		case FUNC_BRISTOR:
			sprintf(fn,"_L%02i_%s_bristor_%s",
			REFINEMENTLEVEL,NNTYPSTR,seedCstr(tmp2));
			getBoundingBoxfA=getBoundingBoxfA_bristor;
			bitpower=2;
			break;
			CAUS
		default:
			_FUNC=FUNC_MAKIN;
			sprintf(fn,"_L%02i_%s_makin_%s",
			REFINEMENTLEVEL,NNTYPSTR,seedCstr(tmp2));
			getBoundingBoxfA=getBoundingBoxfA_makin;
			bitpower=2;
			CAUS
			break;
	}
}

// detect periodic basins
void periodicity(void) {
	// Cycles
	// the temporary components in one orbit until cyclic component met
	FatouComponent *oneorbit=new FatouComponent[MAXFATOUCOMPONENTS];
	int32_t anzfatouinorbit=0;
	// routine-global Components of immediate basins
	FatouComponent* allcycles=new FatouComponent[MAXFATOUCOMPONENTS];
	int32_t anzfatouinCycleNbr=0;
	cycles=new Cycle[MAXCYCLES];
	anzcycles=0;
	
	const int32_t MINTEMPCOLOR=256;
	int32_t orbitfcnbr=MINTEMPCOLOR; // >= 256
	int32_t cyclesetnbrimmediate=FATOUCOMPONENTCOLOROFFSET;
	int32_t cyclesetnbrattraction=FATOUCOMPONENTCOLOROFFSET+1;
	uint16_t blobaktiv=FATOUCOMPONENTCOLOROFFSET-1;
	
	int32_t noch0=SCREENWIDTH >> 4;
	int32_t noch=1;
	const int64_t MAXINLISTE=( (int64_t)1 << 26); // a 1 GB Speicher
	Int3 *liste=new Int3[MAXINLISTE];
	int64_t anzliste=0;
	printf("searching for cycles ");
	if (!liste) {
		printf("\n  periodicity search using cube only (slow routine)\n");
	} 

	#define ADDLISTE(XX,YY,ZZ) \
	{\
		if ( (liste) && (anzliste < (MAXINLISTE-8)) ) {\
			liste[anzliste].x=XX;\
			liste[anzliste].y=YY;\
			liste[anzliste].z=ZZ;\
			anzliste++;\
		}\
	}
	
	#define STRAHLEN(BX,BY,BZ) \
	{\
		for(int32_t xx=(BX-1);xx>=0;xx--) {\
			if (data->getVoxel(xx,BY,BZ) == CC_BLACK) {\
				ADDLISTE(xx,BY,BZ)\
				if (xx < bx0) bx0=xx;\
				if (xx > bx1) bx1=xx;\
				data->setVoxel(xx,BY,BZ,blobaktiv);\
			} else break;\
		}\
		\
		for(int32_t xx=(BX+1);xx<SCREENWIDTH;xx++) {\
			if (data->getVoxel(xx,BY,BZ) == CC_BLACK) {\
				ADDLISTE(xx,BY,BZ)\
				if (xx < bx0) bx0=xx;\
				if (xx > bx1) bx1=xx;\
				data->setVoxel(xx,BY,BZ,blobaktiv);\
			} else break;\
		}\
		\
		for(int32_t yy=(BY+1);yy<SCREENWIDTH;yy++) {\
			if (data->getVoxel(BX,yy,BZ) == CC_BLACK) {\
				ADDLISTE(BX,yy,BZ)\
				if (yy < by0) by0=yy;\
				if (yy > by1) by1=yy;\
				data->setVoxel(BX,yy,BZ,blobaktiv);\
			} else break;\
		}\
		\
		for(int32_t yy=(BY-1);yy>=0;yy--) {\
			if (data->getVoxel(BX,yy,BZ) == CC_BLACK) {\
				ADDLISTE(BX,yy,BZ)\
				if (yy < by0) by0=yy;\
				if (yy > by1) by1=yy;\
				data->setVoxel(BX,yy,BZ,blobaktiv);\
			} else break;\
		}\
		\
		for(int32_t zz=(BZ+1);zz<SCREENWIDTH;zz++) {\
			if (data->getVoxel(BX,BY,zz) == CC_BLACK) {\
				ADDLISTE(BX,BY,zz)\
				if (zz < bz0) bz0=zz;\
				if (zz > bz1) bz1=zz;\
				data->setVoxel(BX,BY,zz,blobaktiv);\
			} else break;\
		}\
		\
		for(int32_t zz=(BZ-1);zz>=0;zz--) {\
			if (data->getVoxel(BX,BY,zz) == CC_BLACK) {\
				ADDLISTE(BX,BY,zz)\
				if (zz < bz0) bz0=zz;\
				if (zz > bz1) bz1=zz;\
				data->setVoxel(BX,BY,zz,blobaktiv);\
			} else break;\
		}\
	}
	
	// exchange colors from temporary to final
	#define COLORFC(FC,QF,ZF) \
	{\
		for(int32_t zz=FC.scrc.z0;zz<=FC.scrc.z1;zz++) {\
			for(int32_t yy=FC.scrc.y0;yy<=FC.scrc.y1;yy++) {\
				for(int32_t xx=FC.scrc.x0;xx<=FC.scrc.x1;xx++) {\
					if (data->getVoxel(xx,yy,zz) == QF) {\
						data->setVoxel(xx,yy,zz,ZF);\
					}\
				}\
			}\
		}\
	}
	
	int32_t maxorbitlen=0;
	
	// search for BLACK cube to detect a new Fatou component
	for(int32_t zb=0;zb<SCREENWIDTH;zb++) {
		if ((--noch)<=0) {
			printf("%i ",SCREENWIDTH-zb);
			noch=noch0;
		}
		for(int32_t yb=0;yb<SCREENWIDTH;yb++) {
			for(int32_t xb=0;xb<SCREENWIDTH;xb++) {
				if (data->getVoxel(xb,yb,zb) != CC_BLACK) continue;
				
				int32_t x=xb,y=yb,z=zb;
				orbitfcnbr=MINTEMPCOLOR;
				anzfatouinorbit=0; // starten
				
				// no Voxel: color >= MINTEMPCOLOR
				// no voxel: color blobaktiv
				
				while (1) {
					// current orbit ongoing
					// starts at x,y,z
				
					int32_t bx0=x,bx1=x,by0=y,by1=y,bz0=z,bz1=z;
					int32_t currorbitidx=anzfatouinorbit;
					oneorbit[currorbitidx].currentOrbitColorIdxTemp=orbitfcnbr;
					oneorbit[currorbitidx].inCycleNbr=-1; 
					oneorbit[currorbitidx].isimmediate=0;
					anzfatouinorbit++;
			
					if (anzliste != 0) {
						LOGMSG3("Implementation error. new blob#%i, but list with %I64d elements\n",anzfatouinorbit,anzliste);
						exit(99);
					
					}
					ADDLISTE(x,y,z)
					bx0=bx1=x;
					by0=by1=y;
					bz0=bz1=z;
					data->setVoxel(x,y,z,blobaktiv);
				
					// Floodfill 
					int32_t changed=1;
					while (changed>0) {
						changed=0;
					
						int32_t ey0=by0,ey1=by1,ex0=bx0,ex1=bx1;
						int32_t ez0=bz0,ez1=bz1;
					
						if (anzliste>0) {
							int32_t lx=liste[anzliste-1].x;
							int32_t ly=liste[anzliste-1].y;
							int32_t lz=liste[anzliste-1].z;
							if (lx < bx0) bx0=lx;
							if (lx > bx1) bx1=lx;
							if (ly < by0) by0=ly;
							if (ly > by1) by1=ly;
							if (lz < bz0) bz0=lz;
							if (lz > bz1) bz1=lz;
							anzliste--;
							if (data->getVoxel(lx,ly,lz) == blobaktiv) {
								data->setVoxel(lx,ly,lz,orbitfcnbr);
								STRAHLEN(lx,ly,lz)
							}
						
							// go on till list is empty
							changed=1;
							continue;
						};
					
						// if list is empty or has changed the cube
						// check for voxels to follow
						if ((anzliste>0) || (changed>0)) {
							LOGMSG2("Implementation error. list mistake %I64d elements, 0 expected\n",anzliste);
							exit(99);
						} else {
							// per cube
							for(int32_t bz=ez0;bz<=ez1;bz++) {
								for(int32_t by=ey0;by<=ey1;by++) {
									for(int32_t bx=ex0;bx<=ex1;bx++) {
										if (data->getVoxel(bx,by,bz) != blobaktiv) continue;
								
										data->setVoxel(bx,by,bz,orbitfcnbr);
										changed=1;
										if (bx < bx0) bx0=bx;
										if (bx > bx1) bx1=bx;
										if (by < by0) by0=by;
										if (by > by1) by1=by;
										if (bz < bz0) bz0=bz;
										if (bz > bz1) bz1=bz;
								
										STRAHLEN(bx,by,bz)
									} // bx
								} // by
							} // bz
						} // per
					} // while
			
					// enclosement cube of current Fatou component
					// SETCOL
					oneorbit[currorbitidx].scrc.x0=bx0;
					oneorbit[currorbitidx].scrc.x1=bx1;
					oneorbit[currorbitidx].scrc.y0=by0;
					oneorbit[currorbitidx].scrc.y1=by1;
					oneorbit[currorbitidx].scrc.z0=bz0;
					oneorbit[currorbitidx].scrc.z1=bz1;
					
					// follow ONE arbitrary voxel to
					// find the target Fatou component
					SpaceCube A,bbxfA;
					A.x0=x*scaleRangePerPixel+RANGE0;
					A.x1=A.x0+scaleRangePerPixel;
					A.y0=y*scaleRangePerPixel+RANGE0;
					A.y1=A.y0+scaleRangePerPixel;
					A.z0=z*scaleRangePerPixel+RANGE0;
					A.z1=A.z0+scaleRangePerPixel;
					
					getBoundingBoxfA(A,bbxfA);
					
					if (TOTALLYINSPECEXT(bbxfA) > 0) {
						LOGMSG("Implementation error. No target Fatou component.\n");
						exit(99);
					}
				
					ScreenCube scr;
					scr.x0=scrcoord0(bbxfA.x0);
					scr.x1=scrcoord0(bbxfA.x1);
					scr.y0=scrcoord0(bbxfA.y0);
					scr.y1=scrcoord0(bbxfA.y1);
					scr.z0=scrcoord0(bbxfA.z0);
					scr.z1=scrcoord0(bbxfA.z1);
					
					if (
						(scr.x0<0) || (scr.x0 >= SCREENWIDTH) ||
						(scr.y0<0) || (scr.y0 >= SCREENWIDTH) ||
						(scr.z0<0) || (scr.z0 >= SCREENWIDTH)
					) {
						LOGMSG("Implementation error. BbxfA is inconsistent.\n");
						exit(99);
					}
					
					int32_t vf=data->getVoxel(scr.x0,scr.y0,scr.z1);
					
					if (vf == CC_BLACK) {
						// new Fatou component found
						x=scr.x0;
						y=scr.y0;
						z=scr.z1;
						orbitfcnbr++;
						
						continue;
					} else if (vf >= MINTEMPCOLOR) {
						// new cycle found
						LOGMSG("\n  Cycle ");
						
						int32_t o0=-1;
						for(int32_t oi=0;oi<=currorbitidx;oi++) {
							if (oneorbit[oi].currentOrbitColorIdxTemp == vf) {
								o0=oi;
								break;
							}
						} // oi
						if (o0 < 0) {
							LOGMSG("Implementation error. No cycle in orbit\n");
							LOGMSG4("anzcycles %i currorbitidx %i currentOrbitColorIdxTemp %i\n",
								anzcycles,currorbitidx,oneorbit[currorbitidx].currentOrbitColorIdxTemp);
							exit(99);
						}
						
						cycles[anzcycles].len=currorbitidx-o0+1;

						LOGMSG2("len %i found\n",cycles[anzcycles].len);
						// up to first cfyclic point => attraction basin
						for(int32_t oi=0;oi<o0;oi++) {
							COLORFC(oneorbit[oi],oneorbit[oi].currentOrbitColorIdxTemp,cyclesetnbrattraction);
						}

						cycles[anzcycles].fatouidx0=anzfatouinCycleNbr;

						for(int32_t oi=o0;oi<=currorbitidx;oi++) {
							// immediate basins
							allcycles[anzfatouinCycleNbr].scrc.x0=oneorbit[oi].scrc.x0;
							allcycles[anzfatouinCycleNbr].scrc.x1=oneorbit[oi].scrc.x1;
							allcycles[anzfatouinCycleNbr].scrc.y0=oneorbit[oi].scrc.y0;
							allcycles[anzfatouinCycleNbr].scrc.y1=oneorbit[oi].scrc.y1;
							allcycles[anzfatouinCycleNbr].scrc.z0=oneorbit[oi].scrc.z0;
							allcycles[anzfatouinCycleNbr].scrc.z1=oneorbit[oi].scrc.z1;
							allcycles[anzfatouinCycleNbr].currentOrbitColorIdxTemp=cyclesetnbrimmediate;
							allcycles[anzfatouinCycleNbr].inCycleNbr=anzcycles;
							allcycles[anzfatouinCycleNbr].isimmediate=1;
							anzfatouinCycleNbr++;
							COLORFC(oneorbit[oi],oneorbit[oi].currentOrbitColorIdxTemp,cyclesetnbrimmediate);
						} // copy
						
						cycles[anzcycles].fatouidx1=-1 + anzfatouinCycleNbr;
						cycles[anzcycles].immediateBasinColorIdx=cyclesetnbrimmediate;
						cycles[anzcycles].attractionBasinColorIdx=cyclesetnbrattraction;
						anzcycles++;
						if (anzcycles >= (MAXCYCLES-2)) {
							LOGMSG("Not possible: Too many cycles (yet to be implemented).\n");
							exit(99);
						}
						
						cyclesetnbrimmediate+=2;;
						cyclesetnbrattraction+=2;
						if (cyclesetnbrattraction >= MINTEMPCOLOR) {
							LOGMSG("Implementation error. Too many cycles.\n");
							exit(99);
						}
						
						if (currorbitidx > maxorbitlen) maxorbitlen=currorbitidx;
						break; // beendet
					} else {
						// orbit lands in already found cycle
						int32_t zyklus=-1;
						for(int32_t cyc=0;cyc<anzcycles;cyc++) {
							if (
								(cycles[cyc].immediateBasinColorIdx == vf) ||
								(cycles[cyc].attractionBasinColorIdx == vf)
							) {
								zyklus=cyc;
								break;
							}
						} // cyc
						
						if (zyklus<0) {
							LOGMSG("Implementation error. Found cycle not detected in orbit\n");
							exit(99);
						}
						
						// whole orbit is attraction basin
						for(int32_t oi=0;oi<=currorbitidx;oi++) {
							COLORFC(oneorbit[oi],oneorbit[oi].currentOrbitColorIdxTemp,cycles[zyklus].attractionBasinColorIdx);
						}
						
						break;
					}
				} // while in Orbit
				
			} // xb
		} // yb
	} // zb
	
	LOGMSG3("\n%i cycles (max. orbit length %i)\n",anzcycles,maxorbitlen);
	for(int32_t i=0;i<anzcycles;i++) {
		LOGMSG2("  Cycle #%i: ",i);
		LOGMSG5("len=%i immediate RGB(%i,%i,%i) ",cycles[i].len,
			pal256[cycles[i].immediateBasinColorIdx].R,
			pal256[cycles[i].immediateBasinColorIdx].G,
			pal256[cycles[i].immediateBasinColorIdx].B);
		LOGMSG4("attraction RGB(%i,%i,%i)\n",
			pal256[cycles[i].attractionBasinColorIdx].R,
			pal256[cycles[i].attractionBasinColorIdx].G,
			pal256[cycles[i].attractionBasinColorIdx].B);
	}
	
	delete[] oneorbit;
	delete[] allcycles;
}

// datatype sufficient mantissa precision
void bitTest(void) {
	if (bitpower < 0) {
		LOGMSG("Implementation error: Exponent missing.\n");
		exit(99);
	}

	SpaceCube A,bbxfA;

	double d=0.0;
	
	#define GROESSER(DD) \
	{\
		double dwert=fabs(DD);\
		if (dwert > d2) d2=dwert;\
	}
	
	#define MAXVALUE(AA,BB,CC,DD,EE,FF) \
	{\
		A.x0=AA;\
		A.x1=BB;\
		A.y0=CC;\
		A.y1=DD;\
		A.z0=EE;\
		A.z1=FF;\
		getBoundingBoxfA(A,bbxfA);\
		double d2=0.0;\
		GROESSER(bbxfA.x0)\
		GROESSER(bbxfA.x1)\
		GROESSER(bbxfA.y0)\
		GROESSER(bbxfA.y1)\
		GROESSER(bbxfA.z0)\
		GROESSER(bbxfA.z1)\
		if (d2 > d) d=d2;\
	}

	// Octant 1-8
	int32_t z0,z1,y0,y1,x0,x1;
	for(int z=0;z<2;z++) {
		if (z==0) { z0=RANGE0; z1=0; }
		else if (z==0) { z0=0; z1=RANGE1-1; }
		for(int y=0;y<2;y++) {
			if (y==0) { y0=RANGE0; y1=0; }
			else if (y==0) { y0=0; y1=RANGE1-1; }
			for(int x=0;x<2;x++) {
				if (x==0) { x0=RANGE0; x1=0; }
				else if (x==0) { x0=0; x1=RANGE1-1; }
				
				MAXVALUE(x0,x1,y0,y1,z0,z1);
			}
		}
	}
	
	d=ceil(d);
	
	int bitsd=(int)ceil(log(d)/log(2.0));
	double sk=(RANGE1-RANGE0); sk /= SCREENWIDTH;
	int bitsnachkomma=(int)ceil(fabs(log(sk)/log(2.0)));
	int bitsused=bitsd + bitpower*bitsnachkomma;
	LOGMSG("mantissa precision ");
	if (bitsused < (_BITSNTYP-1)) {
		LOGMSG3("sufficient (%i needed, %i provided) ",
			bitsused,_BITSNTYP);
	} else {
		LOGMSG3("NOT SUFFICIENT (buffer 2 Bit) (%i needed, %i provided) ",
			bitsused,_BITSNTYP);
	}
	
	LOGMSG2("maximum's decimal part %.0lf\n",d);
}

// memory manager for DBYTE-Array
ArrayMgrDByte::~ArrayMgrDByte() {
	for(int32_t i=0;i<countptr;i++) {
		delete[] ptr[i];
	}
	countptr=0;
}

uint16_t* ArrayMgrDByte::getMem(const int32_t aanz) {
	if (
		(!current) ||
		( (allocateHowMany - freeFromIdx-8) < aanz)
	) {
		// allocated new
		if (countptr >= (MAXPTR-8)) {
			LOGMSG("Memory error. ArrayMgrDByte\n");
			exit(99);
		}
		printf("x");
		ptr[countptr]=current=new uint16_t[allocateHowMany];
		countptr++;
		freeFromIdx=0;
		if (!current) {
			LOGMSG("Memory error. Manager DBYTE\n");
			exit(99);
		}
	}
	
	uint16_t *p=&current[freeFromIdx];
	freeFromIdx += aanz;
	
	return p;
}

ArrayMgrDByte::ArrayMgrDByte() {
	current=NULL;
	int64_t memblock=CHUNKSIZE; // 1 GB-Chunks
	allocateHowMany=memblock / sizeof(uint16_t);
	countptr=0;
}

// memory manager for Coord array
Coord* ArrayMgrCoord::getMem(const int32_t aanz) {
	if (
		(!current) ||
		( (allocateHowMany - freeFromIdx-8) < aanz)
	) {
		if (countptr >= (MAXPTR-8)) {
			LOGMSG("Memory error. ArrayMgrCoord\n");
			exit(99);
		}
		ptr[countptr]=current=new Coord[allocateHowMany];
		countptr++;
		freeFromIdx=0;
		if (!current) {
			LOGMSG("Not enough memory for coords.n");
			exit(99);
		}
	}
	
	Coord *p=&current[freeFromIdx];
	freeFromIdx += aanz;
	
	return p;
}

ArrayMgrCoord::ArrayMgrCoord() {
	current=NULL;
	int64_t memblock=CHUNKSIZE; // 1 GB-Chunks
	allocateHowMany=memblock / sizeof(Coord);
	countptr=0;
}

ArrayMgrCoord::~ArrayMgrCoord() {
	for(int32_t i=0;i<countptr;i++) {
		delete[] ptr[i];
	}
	countptr=0;
}

// writing bitmap data to file
unsigned char dez(const char c) {
	if ((c>='A')&&(c<='F')) { return c-'A'+10; }
	if ((c>='a')&&(c<='f')) { return c-'a'+10; }
	if ((c>='0')&&(c<='9')) { return c-'0'; }
	
	return 0;
}

void writehex(FILE *f,const char* s) {
	for(uint32_t i=0;i<strlen(s);i+=2) {
		unsigned char c = 16*dez(s[i]) + dez(s[i+1]);
		fwrite(&c,sizeof(c),1,f);
	}
}

// internal viewer
// base version of my cube-viewer project with a fixed view point
void internalViewer(const char* afn) {
	// standard-viewer: above one corner looking into the whole cube's center
	
	Screen scr;
	scr.setScreenSize(SCREENWIDTH << 1,SCREENWIDTH << 1); // should be enough for the whole object

	const int32_t hb=SCREENWIDTH >> 1;
	
	Coord cubeM=Coord(hb,hb,hb);
	scr.observer=Coord(-hb,SCREENWIDTH+hb,SCREENWIDTH+hb); 
	Coord obsUL;
	
	Coord s=cubeM - scr.observer;
	
	Coord rx=Coord(-s.y,s.x,0); 
	Coord ry=Coord(
		s.z/(-s.y*s.y/s.x-s.x),
		s.y*s.z/(-s.y*s.y-s.x*s.x),
		1
	);
	
	if ( scr.observer.z < -0.01 ) { 
		rx.mult(-1); // von HINTEN sieht man in die ANDERE x-Richtung
	}

	rx.normiere(); rx.mult(scr.lenx);
	ry.normiere(); ry.mult(scr.leny);
	
	obsUL.x = scr.observer.x - 0.5*rx.x - 0.5*ry.x;	
	obsUL.y = scr.observer.y - 0.5*rx.y - 0.5*ry.y;	
	obsUL.z = scr.observer.z - 0.5*rx.z - 0.5*ry.z;	
	
	scr.setObserver(obsUL,rx,ry,cubeM);
	
	printf("sending virtual screen through cube ...\n");
	scr.slide();
	
	// introducing small colored rectangles in the lower
	// left to indicate the colors of the found cycles
	
	#define COLORCYC(XX0,YY0,FF) \
	{\
		if ( ( (FF) < 0) || ( (FF) >= 256) ) {\
			LOGMSG2("Implementation error. COLORCYC %i\n",(FF));\
			exit(99);\
		}\
		for(int y=(YY0);y<((YY0)+cycwidth);y++) {\
			for(int x=(XX0);x<((XX0)+cycwidth);x++) {\
				scr.coordsY[y][x].R=pal256[FF].R;\
				scr.coordsY[y][x].G=pal256[FF].G;\
				scr.coordsY[y][x].B=pal256[FF].B;\
			}\
		}\
	}

	int32_t cycx=0,cycy=0;
	int32_t cycwidth=(int32_t)( (double)SCREENWIDTH / 640.0 );
	if (cycwidth < 1) cycwidth=1;
	if (SCREENWIDTH <= 512) cycwidth=8; else cycwidth <<= 4;
	for(int cyc=0;cyc<anzcycles;cyc++) {
		COLORCYC(cycx,cycy,cycles[cyc].immediateBasinColorIdx);
		COLORCYC(cycx+cycwidth,cycy,cycles[cyc].attractionBasinColorIdx);
		cycx += (3*cycwidth);
		if (cycx +(3*cycwidth) > scr.lenx) {
			cycx=0;
			cycy += (2*cycwidth);
		}
	}
	
	// saving as 24-RGB-Bitmap
	printf("saving bitmap ...\n");
	char tt[1024];
	sprintf(tt,"%s.bmp",afn);
	FILE *fbmp=fopen(tt,"wb");
	
	int32_t entries_per_row=scr.lenx;
	uint8_t* rgbz=new uint8_t[3*entries_per_row];

	writehex(fbmp,"424DF6C62D00000000003600000028000000");
	fwrite(&scr.lenx,sizeof(scr.lenx),1,fbmp);
	fwrite(&scr.leny,sizeof(scr.leny),1,fbmp);
	writehex(fbmp,"0100180000000000C0C62D00C40E0000C40E00000000000000000000");
	
	// direct values
	
	for(int y=0;y<scr.leny;y++) {
		for(int x=0;x<scr.lenx;x++) {
			rgbz[3*x+0]=scr.coordsY[y][x].B;
			rgbz[3*x+1]=scr.coordsY[y][x].G;
			rgbz[3*x+2]=scr.coordsY[y][x].R;
		}
		fwrite(rgbz,sizeof(uint8_t),entries_per_row*3,fbmp);
	}
		
	fclose(fbmp);
	
	delete[] rgbz;
}

// struct Screen
double Screen::DistanceToScreen(Coord& p) {
	return ( (vektor.x*p.x + vektor.y*p.y + vektor.z*p.z - abstandDelta) / normVektor);
}

Screen::Screen() {
	coordsY=NULL;
	lenx=leny;
}

Screen::~Screen() {
	if ((lenx>0) && (coordsY) ) { 
		delete[] coordsY; // dirty memory free up
		coordsY=NULL;
		lenx=leny=0;
	}
}

void Screen::setScreenSize(const int32_t axl,const int32_t ayl) {
	coordsY=new PCoord[ayl];
	if (!coordsY) {
		LOGMSG("Memory not sufficient: Screen.setsize\n");
		exit(99);
	}
	lenx=axl;
	leny=ayl;
	for(int32_t y=0;y<ayl;y++) {
		coordsY[y]=mgrCoord.getMem(axl);
		if (!coordsY[y]) {
			LOGMSG("Memory not sufficient/2: Screen.setsize\n");
			exit(99);
		}
	}
}

void Screen::setObserver(Coord& obsUL,Coord& rvx,Coord& rvy,Coord& Mcube) {
	rvx.normiere();
	rvy.normiere();
	scrxv=rvx;
	scryv=rvy;
	
	const double x0=0.5*lenx;
	const double y0=0.5*leny;
	
	vektor.x=(Mcube.x-(obsUL.x+x0*rvx.x+y0*rvy.x));
	vektor.y=(Mcube.y-(obsUL.y+x0*rvx.y+y0*rvy.y));
	vektor.z=(Mcube.z-(obsUL.z+x0*rvx.z+y0*rvy.z));
	vektor.normiere();
	
	for(int y=0;y<leny;y++) {
		for(int x=0;x<lenx;x++) {
			coordsY[y][x].x=obsUL.x+x*rvx.x+y*rvy.x;
			coordsY[y][x].y=obsUL.y+x*rvx.y+y*rvy.y;
			coordsY[y][x].z=obsUL.z+x*rvx.z+y*rvy.z;
			coordsY[y][x].intersectsWithVisibleObjectAtStep=-1; // noch nicht getroffen
		}
	}
	
	abstandDelta=vektor.x*obsUL.x + vektor.y*obsUL.y + vektor.z*obsUL.z;
	normVektor=vektor.norm();
}

void Screen::dimm(const double dimmingFactor,const int32_t s0,const int32_t s1) {
	int32_t breite=s0-s1; // erste ja GRö?ER im Index
	double df=1.0; 
	if (breite<=0) breite=1;
	
	df /= breite;
	df *= dimmingFactor;

	for(int32_t y=0;y<leny;y++) {
		for(int32_t x=0;x<lenx;x++) {
			if (coordsY[y][x].intersectsWithVisibleObjectAtStep>0) { 
				double dd=1-df*(s0 - coordsY[y][x].intersectsWithVisibleObjectAtStep);
				coordsY[y][x].R=(int)floor(dd*coordsY[y][x].R);
				coordsY[y][x].G=(int)floor(dd*coordsY[y][x].G);
				coordsY[y][x].B=(int)floor(dd*coordsY[y][x].B);
				coordsY[y][x].dimmingFactor=dd;
			}
		}
	}
}

void getRGBColor(const int32_t cubecol,RGB& ff) {
	if (cubecol <= 255) {
		ff.R=pal256[cubecol].R;
		ff.G=pal256[cubecol].G;
		ff.B=pal256[cubecol].B;
	} else {
		ff.R=ff.G=ff.B=127;
		LOGMSG("cubecol possible error. Color set to gray.\n");
	}
}

void Screen::slide(void) {
	double len=sqrt(3)*SCREENWIDTH;
	int32_t schritte=2*(int)len;
	int32_t MAXSCR=(int)ceil(1.2*len);
	
	int32_t maxanzschritte=2*(int)len;
	schritte=maxanzschritte;
	int32_t intersectidx0=-1,intersectidx1=-1; // invers, d.h. je kleiner Index, desto näher an Cube
	int32_t erster=1;
	const int32_t XLEN=lenx;
	const int32_t YLEN=leny;
	
	const int32_t NOCH0=64;
	int32_t noch=NOCH0;
	
	for(int32_t y=0;y<YLEN;y++) {
		for(int32_t x=0;x<XLEN;x++) {
			coordsY[y][x].intersectsWithVisibleObjectAtStep=-1;
		}
	}
	
	double mina=len*len;
	// distance of the eigth corners to screen
	for(int32_t z=0;z<2;z++) {
		for(int32_t y=0;y<2;y++) {
			for(int32_t x=0;x<2;x++) {
				Coord temp=Coord(
					x*SCREENWIDTH,
					y*SCREENWIDTH,
					z*SCREENWIDTH
				);

				double d=DistanceToScreen(temp);
				if (d < mina) mina=d;
			}
		}
	}
		
	// jumping forth and back to find the world cube in larger steps
	int32_t sprung=(int32_t)floor(0.9*mina);
	
	int32_t x0=0,x1=XLEN-1,y0=0,y1=YLEN-1;
	int32_t ersterPunkt=1;
	int32_t ersterSchritt=schritte;
	int32_t schrittErsterTreffer=-1;
	
	const int32_t SPRUNGSCHRITT=8;
	const int32_t ZURUECKSCHRITTE=64;
	
	int32_t inkubus=1;
	Coord add;
	
	#define WITHIN(XX,YY,ZZ) \
	(\
		(XX >= 0) && (XX < SCREENWIDTH) &&\
		(YY >= 0) && (YY < SCREENWIDTH) &&\
		(ZZ >= 0) && (ZZ < SCREENWIDTH)\
	)
	
	#define WITHINRGB(XX,YY,ZZ,FF,DRIN) \
	{\
		DRIN=0;\
		if (WITHIN(XX,YY,ZZ)) {\
			DRIN=1;\
			FF=data->getVoxel(XX,YY,ZZ);\
		}\
	}

	while (inkubus>0) {
		inkubus=0;
		add=vektor;
		add.mult(sprung);
		for(int32_t y=0;(y<YLEN);y++) {
			for(int32_t x=0;(x<XLEN);x++) {
				if (
					(coordsY[y][x].intersectsWithVisibleObjectAtStep<0)
				) { 
					Coord temp=coordsY[y][x];
					temp = temp + add;
					if (WITHIN(temp.x,temp.y,temp.z)) {
						sprung -= ZURUECKSCHRITTE;
						inkubus=1;
						break;
					}
				}
			} // for x
			if (inkubus>0) break;
		} // for y
	} // while
	
	int32_t zaehle=schritte;
	while (zaehle>0) {
		int32_t treffer=0;
		
		add=vektor; add.mult(sprung);
		for(int32_t y=0;y<YLEN;y++) {
			for(int32_t x=0;x<XLEN;x++) {
				if (
					(coordsY[y][x].intersectsWithVisibleObjectAtStep<0)
				) { 
					Coord temp=coordsY[y][x];
					temp = temp + add;
					if (WITHIN(temp.x,temp.y,temp.z)) {
						treffer=1;
						break;
					}
				}
			}
			if (treffer>0) {
				break;
			}
		} // fors

		if (treffer>0) {
			sprung -= SPRUNGSCHRITT;
			break;
		}
		
		sprung += SPRUNGSCHRITT;
		zaehle -= SPRUNGSCHRITT;
	}
	
	if ( (sprung > 0) && (sprung < 32000) ) {
		Coord vsprung=vektor;
		vsprung.mult(sprung);
		for(int32_t y=0;y<YLEN;y++) {
			for(int32_t x=0;x<XLEN;x++) {
				coordsY[y][x] += vsprung;
			}
		}
	}

	schritte -= sprung; // eigentlich egal, wenn nix mehr kommt, hoert er eh auf)
	
	while (schritte>0) {
		
		if ( (noch--)== 0) {
			printf("%i ",schritte);
			noch=NOCH0;
		}
		schritte--;
		int32_t scrschneidet=0;
		
		for(int32_t y=y0;y<=y1;y++) {
			for(int32_t x=x0;x<=x1;x++) {
				int32_t cubecol;
				
				if (
					(coordsY[y][x].intersectsWithVisibleObjectAtStep<0)
				) { 
					int32_t drin=0;
					WITHINRGB(
						coordsY[y][x].x,coordsY[y][x].y,coordsY[y][x].z,
						cubecol,
						drin
					);
					
					if (drin>0) {
						scrschneidet=1;
						x0=maximumI(x0,x-MAXSCR);
						x1=minimumI(x1,x+MAXSCR);
						y0=maximumI(y0,y-MAXSCR);
						y1=minimumI(y1,y+MAXSCR);
						if (ersterPunkt>0) {
							schrittErsterTreffer=schritte;
							ersterPunkt=0;						
						}
						
						if (
							(cubecol==CC_WHITE) ||
							(
								(
									(cubecol==CC_GRAY) ||
									(cubecol==CC_GRAYPOTW)
								) &&
								(GRAY_IS_TRANSPARENT>0)
							)
						) {
							coordsY[y][x] += vektor;
						} else {
							if (intersectidx0<0) intersectidx0=schritte;
							intersectidx1=schritte;
							
							coordsY[y][x].intersectsWithVisibleObjectAtStep=schritte;
							
							RGB ff;
							getRGBColor(cubecol,ff);
							coordsY[y][x].R=ff.R;
							coordsY[y][x].G=ff.G;
							coordsY[y][x].B=ff.B;
						}
					} else {
						// move screen point
						coordsY[y][x] += vektor;
					}
				} 
			} // for x
		} // for y
		
		if (scrschneidet<=0) {
			if (erster==0) {
				// exit loop: no part of the cube (not the object) intersected with
				break;
			}
		} else if (erster>0) erster=0;
		
	} // while
	
	// dimming by distance and normal vector
	printf("\ndimming screen pixels ...\n");
	dimm(0.75,intersectidx0,intersectidx1);
	dimmByNormal();
	
	if (schrittErsterTreffer>0) if (schrittErsterTreffer==ersterSchritt) {
		LOGMSG("Cube error. Not able to find object.\n");
		exit(99);
	}
}

void swap(Coord& a,Coord& b) {
	Coord c=a;
	a=b;
	b=c;
}

inline void crossproduct(Coord& a,Coord& b,Coord& e) {
	e.x=a.y*b.z-a.z*b.y;
	e.y=a.z*b.x-a.x*b.z;
	e.z=a.x*b.y-a.y*b.x;
}

inline double scalarproduct(Coord& a,Coord& b) {
	return (
		a.x*b.x + a.y*b.y + a.z*b.z
	);
}

// EXPERIMENTAL
// routine to determine the normal vector
// and its angle to the viewing/light axis
void Screen::dimmByNormal(void) {
	const double VNORM=vektor.norm();
	
	for(int32_t y=1;y<leny;y++) {
		for(int32_t x=1;x<lenx;x++) {
			int32_t punkte=0;
			int32_t nidx=-1;
			double mindQ=-1;
			Coord ps[4];
			for(int32_t dy=-1;dy<=0;dy++) {
				for(int32_t dx=-1;dx<=0;dx++) {
					if (coordsY[y+dy][x+dx].intersectsWithVisibleObjectAtStep>=0) {
						ps[punkte]=coordsY[y+dy][x+dx];
						punkte++;
						Coord dist=observer - coordsY[y+dy][x+dx];
						double d=dist.normQ();
						if ((d < mindQ) || (mindQ<-0.5)) {
							mindQ=d;
							nidx=(punkte-1);
						}
					}
				}
			}
			
			if (punkte <= 2) continue;
			
			Coord nvektor;
			if (punkte==3) {
				Coord v1,v2;
				if (nidx != 0) swap(ps[0],ps[nidx]);
				nidx=0;
				v1=ps[1]-ps[0];
				v2=ps[2]-ps[0];
				// crossproduct
				crossproduct(v1,v2,nvektor);
			} else if (punkte==4) {
				Coord v1,v2,e;
				nvektor=Coord(0,0,0);
				if (nidx != 0) swap(ps[0],ps[nidx]);
				nidx=0;
				
				#define EBENE(AA,BB) \
				v1=ps[AA]-ps[0];\
				v2=ps[BB]-ps[0];\
				crossproduct(v1,v2,e);\
				nvektor += e;\
				
				EBENE(1,2)
				EBENE(1,3)
				EBENE(2,3)
			} // 4 points
			
			double nvn=nvektor.norm();
			if (fabs(nvn)<1E-5) {
				// nullvector
				continue;
			}

			double alpha=
				fabs(
				scalarproduct(vektor,nvektor)
				/ (VNORM*nvn)
			);
			// not dimming to full black
			alpha=1-0.2*alpha;
			for(int32_t dy=-1;dy<=0;dy++) {
				for(int32_t dx=-1;dx<=0;dx++) {
					if (coordsY[y+dy][x+dx].intersectsWithVisibleObjectAtStep<0) continue;

					coordsY[y+dy][x+dx].R=(int)floor(alpha*coordsY[y+dy][x+dx].R);
					coordsY[y+dy][x+dx].G=(int)floor(alpha*coordsY[y+dy][x+dx].G);
					coordsY[y+dy][x+dx].B=(int)floor(alpha*coordsY[y+dy][x+dx].B);
					coordsY[y+dy][x+dx].dimmingFactor*=alpha;
				}
			}
		}
	} // y
}

// struct Coord

double Coord::norm(void) { 
	return sqrt(x*x+y*y+z*z); 
}

double Coord::normQ(void) { 
	return (x*x+y*y+z*z); 
}

Coord& Coord::operator=(const Coord& a) {
	if (this != &a) {
		x=a.x;
		y=a.y;
		z=a.z;
	}
	return *this;
}
 
void Coord::normiere(void) {
	double d=1.0; d /= norm();
	mult(d);
}

void Coord::mult(const double d) {
	x *= d;
	y *= d;
	z *= d;
}

Coord::Coord() {
	x=y=z=0.0;
}

Coord::Coord(const double ax,const double ay,const double az) {
	x=ax;
	y=ay;
	z=az;
}

Coord::Coord(const Coord& a) {
	x=a.x;
	y=a.y;
	z=a.z;
}
bool operator==(Coord const& lhs,Coord const& rhs) {
	if (fabs(lhs.x-rhs.x) > 1E-300) return 0;
	if (fabs(lhs.y-rhs.y) > 1E-300) return 0;
	if (fabs(lhs.z-rhs.z) > 1E-300) return 0;

	return 1;
}

Coord& Coord::operator+=(const Coord& a) {
	x += a.x;
	y += a.y;
	z += a.z;

	return *this; 
}

int main(int argc,char** argv) {
	getBoundingBoxfA=getBoundingBoxfA_makin;
	_FUNC=FUNC_MAKIN;

	printf("juliatza3dcore\n  CMD=CALC or CMD=PERIOD\n  FUNC=string / GRAY / TWDB=n / LEN=n / RANGE=n\n");
	printf("  C=x,y,z or x,x,y,y,z,z or CI=nx,ny,nz or CI=nx0,nx1,ny0,ny1,nz0,nz1\n");
	printf("  A,B=x,y,z or AI,BI=nx,ny,nz / E,F,...Q=r\n");
	flog=fopen("juliatsa3dcore.log","at");
	fprintf(flog,"\n\n----------------------\n\n");
	
	initRGBPal();
	CMD=CMD_PERIODICITY;
	char fn[1024];

	if (argc>1) {
		for(int i=1;i<argc;i++) {
			upper(argv[i]);
			if (strstr(argv[i],"CMD=")==argv[i]) {
				if (
					(strstr(argv[i],"CMD=PERIOD")==argv[i]) ||
					(strstr(argv[i],"CMD=BASIN")==argv[i])
				) {
					CMD=CMD_PERIODICITY;
				} else
				if (strstr(argv[i],"CMD=CALC")==argv[i]) CMD=CMD_CALC;
			} else
			if (strstr(argv[i],"GRAY")==argv[i]) {
				GRAY_IS_TRANSPARENT=0;
			} else
			if (strstr(argv[i],"FUNC=")==argv[i]) {
				_FUNC=FUNC_MAKIN;
				for(int k=0;k<funcanz;k++) {
					if (!strcmp(&argv[i][5],funcname[k])) {
						_FUNC=k;
						break;
					}
				}
			} 
			else if (strstr(argv[i],"C=")==argv[i]) {
				double vx0,vx1,vy0,vy1,vz0,vz1;
				if (sscanf(&argv[i][2],"%lf,%lf,%lf,%lf,%lf,%lf",&vx0,&vx1,&vy0,&vy1,&vz0,&vz1) == 6) {
					set_seedC_int64_t_XXYYZZ(
						(int64_t)floor(vx0*DENOM225),(int64_t)floor(vx1*DENOM225),
						(int64_t)floor(vy0*DENOM225),(int64_t)floor(vy1*DENOM225),
						(int64_t)floor(vz0*DENOM225),(int64_t)floor(vz1*DENOM225)
					);
				} else
				if (sscanf(&argv[i][2],"%lf,%lf,%lf",&vx0,&vy0,&vz0) == 3) {
					set_seedC_int64_t_XYZ(
						(int64_t)floor(vx0*DENOM225),
						(int64_t)floor(vy0*DENOM225),
						(int64_t)floor(vz0*DENOM225)
					);
				}
			} 
			else if (strstr(argv[i],"CI=")==argv[i]) {
				int64_t vx0,vx1,vy0,vy1,vz0,vz1;
				if (sscanf(&argv[i][3],"%I64d,%I64d,%I64d,%I64d,%I64d,%I64d",&vx0,&vx1,&vy0,&vy1,&vz0,&vz1) == 6) {
					set_seedC_int64_t_XXYYZZ(vx0,vx1,vy0,vy1,vz0,vz1);
				} else
				if (sscanf(&argv[i][3],"%I64d,%I64d,%I64d",&vx0,&vy0,&vz0) == 3) {
					set_seedC_int64_t_XYZ(vx0,vy0,vz0);
				}
			} 
			else if (strstr(argv[i],"TWDB")==argv[i]) {
				int a;
				if (sscanf(&argv[i][5],"%i",&a) == 1) _TWDBITCUBE=a;
				else _TWDBITCUBE=0; // original size
			} else if (strstr(argv[i],"A=")==argv[i]) {
				double vx,vy,vz;
				if (sscanf(&argv[i][2],"%lf,%lf,%lf",&vx,&vy,&vz) == 3) {
					set_TRICA_int64_t_XYZ(
						(int64_t)floor(vx*DENOM225),
						(int64_t)floor(vy*DENOM225),
						(int64_t)floor(vz*DENOM225)
					);
				}
			} else if (strstr(argv[i],"AI=")==argv[i]) {
				int64_t vx,vy,vz;
				if (sscanf(&argv[i][3],"%I64d,%I64d,%I64d",&vx,&vy,&vz) == 3) {
					set_TRICA_int64_t_XYZ(vx,vy,vz);
				}
			} 
			else if (strstr(argv[i],"B=")==argv[i]) {
				double vx,vy,vz;
				if (sscanf(&argv[i][2],"%lf,%lf,%lf",&vx,&vy,&vz) == 3) {
					set_TRICB_int64_t_XYZ(
						(int64_t)floor(vx*DENOM225),
						(int64_t)floor(vy*DENOM225),
						(int64_t)floor(vz*DENOM225)
					);
				}
			} 
			else if (strstr(argv[i],"BI=")==argv[i]) {
				int64_t vx,vy,vz;
				if (sscanf(&argv[i][3],"%I64d,%I64d,%I64d",&vx,&vy,&vz) == 3) {
					set_TRICB_int64_t_XYZ(vx,vy,vz);
				}
			} else
			if (strstr(argv[i],"LEN=")==argv[i]) {
				int a;
				if (sscanf(&argv[i][4],"%i",&a) == 1) {
					if (a >= 256) SCREENWIDTH=a;
					else if (a < 32) SCREENWIDTH=(1 << a);
				}
			} else
			if (strstr(argv[i],"RANGE=")==argv[i]) {
				int64_t a;
				if (sscanf(&argv[i][6],"%I64d",&a) == 1) {
					if (a<0) a=-a;
					RANGE0=-a;
					RANGE1= a;
				}
			} else
			if (argv[i][1]=='=') {
				// eine andere Variable
				double a;
				if (sscanf(&argv[i][2],"%lf",&a) == 1) {
					switch (argv[i][0]) {
						case 'E': varE=a;break;
						case 'F': varF=a;break;
						case 'G': varG=a;break;
						case 'H': varH=a;break;
						case 'J': varJ=a;break;
						case 'K': varK=a;break;
						case 'L': varL=a;break;
						case 'M': varM=a;break;
						case 'N': varN=a;break;
						case 'O': varO=a;break;
						case 'P': varP=a;break;
						case 'Q': varQ=a;break;
					}
				}
			}
		} // i
	}
	
	REFINEMENTLEVEL=(int)ceil(log(SCREENWIDTH)/log(2.0));
	// calculate the width of a voxel
	calc();
	
	// sets function pointers and filename's principal part
	setfunc(fn);
	LOGMSG2("\n%s\n",fn);
	
	// check for sufficient datatype precision
	bitTest();
	
	printf("initialising main object ...\n");
	data=new Data3Ddyn;
	// reading and increasing data
	int32_t geladen=data->readRaw();
	if (geladen <= 0) {
		// if file not there => start with a totally gray world
		data->initComplete();
		printf("\nsearching for special exterior ... ");
		data->search_specext();
	}
	
	//////////////////////////////////////////////
	
	// ain routine to detect new white and potnetially white cells
	data->search_directWhite();

	//////////////////////////////////////////////

	printf("\ncoloring interior ");
	// everything gray but not potentially white has no path
	// to white and is therefore bounded
	int32_t schwarz=data->coloring_black();
	
	printf("\nsaving final data ... ");
	data->saveRaw("_3d.raw");
	fflush(flog);
	printf("done. Safe to close window.\n");
	
	if (
		(CMD == CMD_CALC) || 
		(schwarz<=0)
	) {
		if (schwarz <= 0) GRAY_IS_TRANSPARENT=0; // overall shape
		
		if (_TWDBITCUBE >= 0) {
			printf("\nconvert to bitcube ...\n");
			data->convert_to_BitCube("bc1.ccb");
			LOGMSG4("%I64d interior, %I64d exterior, %I64d gray cells\n",
				ctrschwarz,ctrweiss,ctrgrau+ctrpotwgrau);
		} else {
			// if no interior present (or GRAY was specified in
			// the command line) => show gray shape
			// otherwise show interior cells in yellow
			internalViewer(fn);
		}
	}

	if (
		(CMD==CMD_PERIODICITY) &&
		(schwarz>0)
	) {
		// detects and colors periodic cycles and their basins
		
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// NOTE: the data-structure struct data3d is no longer
		// containing only CC_xxx colors
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		periodicity();
		if (_TWDBITCUBE >= 0) {
			printf("\nconverting to bitcube ...\n");
			data->convert_to_BitCube("bc1per.ccb");
		} else {
			internalViewer(fn);
		}
	}
	
	fclose(flog);
	delete data;
	if (cycles) delete[] cycles;
	
    return 0;
}

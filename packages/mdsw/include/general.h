//Last Modified : Sun Dec  2 19:08:06 2007
#ifndef _GENERAL_H
#define _GENERAL_H


//=============================================================
//General headers
//#ifdef _ENDIAN_H //Clearn <endian.h> definition
//#undef __LITTLE_ENDIAN
//#undef __BIG_ENDIAN
//#else
//#define _ENDIAN_H
//#endif //_ENDIAN_H

#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <fcntl.h>

#include <unistd.h>
#include <errno.h>
#include <signal.h>

#include <time.h>
#include <stdarg.h>
#include <sys/stat.h>

#include <math.h>

//=============================================================================
//Scaler stuff
//All math constants starts with M_ (only non-composit numbers included)
//Extended precision to 60 digits (using Maple)
//This list will increase upon usage
#ifndef M_SQRT2
#define M_SQRT2 (1.41421356237309504880168872420969807856967187537694807317668)
#endif 
#ifndef M_SQRT3 //sqrt(3)
#define M_SQRT3 (1.73205080756887729352744634150587236694280525381038062805581)
#endif
#ifndef M_SQRT5
#define M_SQRT5 (2.23606797749978969640917366873127623544061835961152572427090)
#endif
#ifndef M_SQRT7
#define M_SQRT7 (2.64575131106459059050161575363926042571025918308245018036833)
#endif

#ifndef M_CBRT2 //2^(1/3)
#define M_CBRT2 (1.25992104989487316476721060727822835057025146470150798008198)
#endif
#ifndef M_CBRT3
#define M_CBRT3 (1.44224957030740838232163831078010958839186925349935057754642)
#endif
#ifndef M_CBRT5
#define M_CBRT5 (1.70997594667669698935310887254386010986805511054305492438286)
#endif
#ifndef M_CBRT7
#define M_CBRT7 (1.91293118277238910119911683954876028286243905034587576621065)
#endif

#ifndef M_E
#define M_E (2.71828182845904523536028747135266249775724709369995957496697)
#endif
#ifndef M_PI
#define M_PI    (3.14159265358979323846264338327950288419716939937510582097494)
#endif
#ifndef M_SQRTPI //sqrt(pi)
#define M_SQRTPI 1.77245385090551602729816748334114518279754945612238712821381
#endif
#ifndef M_CBRTPI
#define M_CBRTPI 1.46459188756152326302014252726379039173859685562793717435726
#endif

//Numerical cutoff
#ifndef MIN_DOUBLE
#define MIN_DOUBLE (2.220446e-16)
#endif
// IEEE 754 32-bit floating-point standard: 1.175494e-38 to 3.402823e38 
// IEEE 754 64-bit floating-point standard: 2.225074e-308 to 1.797693e308 
#define DOUBLE_PRECISION_INFINITY (1e308)

//=============================================================================
//Physics Stuff
//All physical constants starts with P_
#define P_C 2.99792458e8  // (m/s) the speed of light (EXACT)
#define P_HBAR 1.0545919e-34 // (kg*m^2/s) hbar plank constant
#define P_KB 1.380662e-23  // (kg*m^2/s^2/K) Boltzmann constant
#define P_E 1.6021892e-19  // (C) electron charge
#define P_ME 9.109558e-31  // (kg) electron mass
#define P_G 6.6732e-11      // (m^3/s^2/kg) Gravitational constant
#define P_NA 6.022169e23   // (1/mol) Avogadro constant
#define P_U 1.660531e-27    // (kg) atomic mass unit
#define P_MU0 (4*M_PI*1e-7) // (C^2*s^4/kg/m^5)
#define P_EPSILON0 (1/P_MU0/P_C/P_C) // (kg*m^3/C^2/s^2)
//Equality sqrt(P_MU0*P_EPSILON0)*P_C=1

//=============================================================
//Site definitions
//=============================================================
//Checking for necessary conditions (only for new CC compiler)


//#if !defined(__LITTLE_ENDIAN) && !defined(__BIG_ENDIAN)
//
//#warning Define __LITTLE_ENDIAN or __BIG_ENDIAN(def) !
//#define __LITTLE_ENDIAN
//
//#endif


//============================================================
//General declearations and definitions


//-----------------------------------
//Color codes $FIXME$ NameSpace pollution

#define ESC	"\x1b"

//  Foreground Colors

#define BLK ESC"[30m"          // Black
#define RED ESC"[31m"          // Red
#define GRN ESC"[32m"          // Green
#define YEL ESC"[33m"          // Yellow
#define BLU ESC"[34m"          // Blue
#define MAG ESC"[35m"          // Magenta
#define CYN ESC"[36m"          // Cyan
#define WHT ESC"[37m"          // White

//   Hi Intensity Foreground Colors

#define HIK ESC"[1;30m"	       // Black
#define HIR ESC"[1;31m"        // Red
#define HIG ESC"[1;32m"        // Green
#define HIY ESC"[1;33m"        // Yellow
#define HIB ESC"[1;34m"        // Blue
#define HIM ESC"[1;35m"        // Magenta
#define HIC ESC"[1;36m"        // Cyan
#define HIW ESC"[1;37m"        // White

// High Intensity Background Colors

#define HBRED ESC"[41;1m"       // Red
#define HBGRN ESC"[42;1m"       // Green
#define HBYEL ESC"[43;1m"       // Yellow
#define HBBLU ESC"[44;1m"       // Blue
#define HBMAG ESC"[45;1m"       // Magenta
#define HBCYN ESC"[46;1m"       // Cyan
#define HBWHT ESC"[47;1m"       // White

//  Background Colors

#define BBLK ESC"[40m"          // Black
#define BRED ESC"[41m"          // Red
#define BGRN ESC"[42m"          // Green
#define BYEL ESC"[43m"          // Yellow
#define BBLU ESC"[44m"          // Blue
#define BMAG ESC"[45m"          // Magenta
#define BCYN ESC"[46m"          // Cyan
#define BWHT ESC"[47m"          // White

#define NOR ESC"[2;37;0m"      // Puts everything back to normal

//   Additional controls

// #define UDL ESC"[2m"
// #define BOLD ESC"[1m"          // Turn on bold mode
// #define CLR ESC"[2J"           // Clear the screen
// #define HOME ESC"[H"           // Send cursor to home position
// #define REF CLR HOME            // Clear screen and home cursor
// #define BIGTOP ESC"#3"         /* Dbl height characters, top half */
// #define BIGBOT ESC"#4"         /* Dbl height characters, bottem half */
// #define SAVEC ESC"[s"           // Save cursor position
// #define REST ESC"[u"            // Restore cursor to saved position
// #define REVINDEX ESC"M"        // Scroll screen in opposite direction
// #define SINGW ESC"#5"          /* Normal, single-width characters */
// #define DBL ESC"#6"            /* Creates double-width characters */
// #define FRTOP ESC"[2;25r"      // Freeze top line
// #define FRBOT ESC"[1;24r"      // Freeze bottom line
// #define UNFR ESC"[r"           // Unfreeze top and bottom lines
// #define BLINK ESC"[5m"         // Initialize blink mode
// #define U ESC"[4m"             // Initialize underscore mode
// #define REV ESC"[7m"           // Turns reverse video mode on
// #define HIREV ESC"[1,7m"       // Hi intensity reverse video


#ifdef __cplusplus

//These ENDIAN dependent subroutines should go to a separate file
//#ifdef __LITTLE_ENDIAN
//inline char LoByte(unsigned short s)
//{ return *(char *)&s; }
//inline char HiByte(unsigned short s)
//{ return *(1+(char *)&s); }
//inline char Byte0(unsigned u)
//{ return *(char *)&u; }
//inline char Byte1(unsigned u)
//{ return *(1+(char *)&u); }
//inline char Byte2(unsigned u)
//{ return *(2+(char *)&u); }
//inline char Byte3(unsigned u)
//{ return *(3+(char *)&u); }

//#else //defined(__BIG_ENDIAN)
//inline char LoByte(unsigned short s)
//{ return *(1+(char *)&s); }
//inline char HiByte(unsigned short s)
//{ return *(char *)&s; }
//inline char Byte3(unsigned u)
//{ return *(char *)&u; }
//inline char Byte2(unsigned u)
//{ return *(1+(char *)&u); }
//inline char Byte1(unsigned u)
//{ return *(2+(char *)&u); }
//inline char Byte0(unsigned u)
//{ return *(3+(char *)&u); }
//#endif

//===========================================================
// Debuging

#ifdef _DEBUG
#define bad_alloc //What is the REAL standard??

//Overloading the global new and delete operators
void *operator new[](size_t) throw(bad_alloc);
void *operator new[](size_t, const char *, unsigned) throw(bad_alloc);
void operator delete[](void *) throw();

void *operator new(size_t)throw(bad_alloc);
void *operator new(size_t, const char *, unsigned)throw(bad_alloc);
void operator delete(void *) throw();

#define new new(__FILE__, __LINE__) //How to solve placement new?
#endif


//==================================================================
//Inline helper functions
#include <sys/time.h>
#include <sys/resource.h>

inline void udelay(long usec)
{
    //Define to use timeval without including <sys/time.h>
//    struct {time_t tv_sec, tv_usec;} tv;
    struct timeval tv;
    tv.tv_sec=usec/1000000;
    tv.tv_usec=usec%1000000;
//    select(0,NULL,NULL,NULL,(timeval *)&tv);
    select(0,NULL,NULL,NULL,&tv);    
}
inline void delay(double sec){ udelay((long)(sec*1e6)); }
inline double CPUTime()
{
    struct rusage ru;
    getrusage(RUSAGE_SELF,&ru);
    return ru.ru_utime.tv_sec+1e-6*ru.ru_utime.tv_usec;
}


unsigned UpdateCRC(unsigned crc, unsigned char *buf, int len);

inline unsigned CRC(unsigned char *buf, int len)
{ return UpdateCRC(0, buf, len); }

template <class T>inline T Abs(T a)
{ return a>=0?a:-a; }
template <class T>inline T Max(T a, T b)
{ return a>b? a:b; }
template <class T>inline T Max(T a, T b, T c)
{ return Max(a, Max(b, c)); }
template <class T>inline T Max(T a, T b, T c, T d)
{ return Max(Max(a,b),Max(c,d)); }
template <class T>inline T Min(T a, T b)
{ return a>b? b:a; }
template <class T>inline T Min(T a, T b, T c)
{ return Min(Min(a,b),c); }
template <class T>inline void Swap(T &a, T &b)
{ T c;c=a;a=b;b=c; }
template <class T>inline void Swap(T &a, T &b, T &c)
{ T d;d=a;a=b;b=c;c=d; }
template <class T>inline T Square(T v)
{ return v*v; }
template <class T>inline T Cube(T v)
{ return v*v*v; }
template <class T>inline T iMod(T a, T b)
{ return a>=0 ? a%b : (b-1)-(-a-1)%b; }
//Mainly for int, unsigned, long, float, double, and maybe char
template <class T>inline T Factorial(T n)
{
    T nf;
    for(nf=1;n>1;nf*=n--);
    return nf;
}

inline double drand() { return (rand()+0.5)/(RAND_MAX+1.0); }
//Returns normal distribution of (0, 1)
inline double randnorm() 
{
    static double w=0;
    double t;
    if(w==0) {
        double r1, r2;
        r1=drand();
        r2=drand();
        t=sqrt(-log(r1)*2);
        w=t*cos(M_PI*2*r2);
        return t*sin(M_PI*2*r2);
    } else {
        t=w; w=0;
        return t;
    }
}
//Returns normal distribution of (xbar, sigma)
inline double RandNorm(double xbar, double sigma)
{ return randnorm()*sigma+xbar; }

inline double randnorm48()
{
    static double w=0;
    double t;
    if (w==0) {
        double r1,r2;
        r1=drand48();
        r2=drand48();
        t=sqrt(-log(r1)*2);
        w=t*cos(M_PI*2*r2);
        return t*sin(M_PI*2*r2);
    } else {
        t=w; w=0;
        return t;
    }
}
//Returns normal distribution of (xbar, sigma)
inline double RandNorm48(double xbar, double sigma)
{ return randnorm48()*sigma+xbar; }

/////////////////////////////////////////////////////////////
// The QuickSort function
// Function QuickSort(T *data, int ndata)
//    Sort the arry by swapping, store in the same place
// Function IndexQuickSort(int *idx, T const *data, int ndata)
//    Sort the data by swapping the index, the index must be prepared with
//    ndata indices, which refers to data array (not necessarily distinct,
//    consecutive or ordered). The sort is done by swapping the index.

#define PUSH(low, high)	((void) ((top->lo = (low)), (top->hi = (high)), ++top))
#define	POP(low, high)	((void) (--top, (low = top->lo), (high = top->hi)))
template <class T> void QuickSort(T *data, int ndata)
{
    const unsigned MaxThresh=4;
    static struct TStack { T *lo, *hi; } stack[8*sizeof(int)];

    T *base_ptr = data;
    T pivot;

    if(ndata == 0) return;

    if((unsigned)ndata > MaxThresh)
    {
        T *lo = data;
        T *hi = &lo[ndata-1];
        
        TStack *top = stack + 1;
        
        while(stack < top)
        {
            T *left_ptr;
            T *right_ptr;
            
            /* Select median value from among LO, MID, and HI. Rearrange
               LO and HI so the three values are sorted. This lowers the
               probability of picking a pathological pivot value and
               skips a comparison for both the LEFT_PTR and RIGHT_PTR. */
            
            T *mid = lo +  ((hi - lo) >> 1);
            
            if(*mid<*lo) Swap(*mid, *lo);
            if(*hi<*mid) Swap(*mid, *hi);
            else goto jump_over;
            if(*mid<*lo) Swap(*mid, *lo);
        jump_over:
            pivot=*mid;
            
            left_ptr  = lo+1;
            right_ptr = hi-1;
            
            /* Here's the famous ``collapse the walls'' section of quicksort.
               Gotta like those tight inner loops!  They are the main reason
               that this algorithm runs much faster than others. */
            do
            {
                while (*left_ptr<pivot) left_ptr++;
                
                while (pivot<*right_ptr) right_ptr--;
                
                if (left_ptr < right_ptr)
                {
                    Swap(*left_ptr, *right_ptr);
                    left_ptr++;
                    right_ptr--;
                }
                else if (left_ptr == right_ptr)
                {
                    left_ptr ++;
                    right_ptr --;
                    break;
                }
            }
            while (left_ptr <= right_ptr);
            
            /* Set up pointers for next iteration.  First determine whether
               left and right partitions are below the threshold size.  If so,
               ignore one or both.  Otherwise, push the larger partition's
               bounds on the stack and continue sorting the smaller one. */
            
            if ((size_t) (right_ptr - lo) <= MaxThresh)
            {
                if ((size_t) (hi - left_ptr) <= MaxThresh)
                    /* Ignore both small partitions. */
                    POP(lo, hi);
                else
                    /* Ignore small left partition. */
                    lo = left_ptr;
            }
            else if ((size_t) (hi - left_ptr) <= MaxThresh)
                /* Ignore small right partition. */
                hi = right_ptr;
            else if ((right_ptr - lo) > (hi - left_ptr))
            {
                /* Push larger left partition indices. */
                PUSH(lo, right_ptr);
                lo = left_ptr;
            }
            else
            {
                /* Push larger right partition indices. */
                PUSH(left_ptr, hi);
                hi = right_ptr;
            }
        }
    }
    
    /* Once the BASE_PTR array is partially sorted by quicksort the rest
       is completely sorted using insertion sort, since this is efficient
       for partitions below MAX_THRESH size. BASE_PTR points to the beginning
       of the array to sort, and END_PTR points at the very last element in
       the array (*not* one beyond it!). */

    {
        T * const end_ptr = &data[ndata-1];
        T *tmp_ptr = base_ptr;
        T *thresh = Min(end_ptr, base_ptr + MaxThresh);
        register T *run_ptr; //key word register could be a helpful hint
        
        /* Find smallest element in first threshold and place it at the
           array's beginning.  This is the smallest array element,
           and the operation speeds up insertion sort's inner loop. */

        for (run_ptr = tmp_ptr+1; run_ptr <= thresh; run_ptr ++)
            if (*run_ptr < *tmp_ptr)
                tmp_ptr = run_ptr;

        if (tmp_ptr != base_ptr) Swap(*tmp_ptr, *base_ptr);

        /* Insertion sort, running from left-hand-side up to right-hand-side.*/

        run_ptr = base_ptr+1;
        while (++run_ptr <= end_ptr)
        {
            tmp_ptr = run_ptr-1;
            while (*run_ptr<*tmp_ptr)
                tmp_ptr --;

            tmp_ptr ++;
            if (tmp_ptr != run_ptr)
            {
                T *trav;

                trav = run_ptr+1;
                while (--trav >= run_ptr)
                {
                    T c = *trav;
                    T *hi, *lo;

                    for (hi = lo = trav; --lo >= tmp_ptr; hi = lo)
                        *hi = *lo;
                    *hi = c;
                }
            }
        }
    }
}
#undef PUSH
#undef POP
template <class T> void QuickSortTo(T *dest, const T *src, int ndata)
{
    memcpy(dest, src, ndata*sizeof(T)); //Should have better way to do this
    QuickSort(dest, ndata);
}
/*
static void * TIndex_data;
template <class T> class TIndex
{
public:
    int id;
    //static T const * data;
    bool operator < (class TIndex<T> b)
    {
        return ((T const *)TIndex_data)[id] < ((T const *)TIndex_data)[b.id];
    }
};

template <class T> void IndexQuickSort(T * const data, int *index, int ndata)
{
    ASSERT(sizeof(TIndex<T>)==sizeof(int));
    TIndex_data=data;
    QuickSort((TIndex<T> *)index, ndata);
}
*/
//=======================================================================
// Base of Stream

// Note: The operators(<, >) all try to call Read/Write only once!
//       the operators(<=, >=) Endian-independent binary read/write
//       the operators(<<,>>) maybe write multiple times.

template <class T>
inline void ReverseTo(T *dest, const T *src, int n)
{ for(int i=0;i<n;i++) dest[i]=src[n-i-1]; }

template <class T>
inline void Reverse(T *data, int n)
{ for(int i=0;i<n/2;i++) Swap(data[i], data[n-i-1]); }

class LIStream
{
public:
    enum { BUFLEN=500 };
private:
//    static char *inbuf;
public:
    // All input stream must have a read function, pure virtual in base class
    virtual int Read(char *buffer, int size)=0;//Should throw $ not return
    // The lock function, do nothing now
    virtual void Lock() {}
    // Returns strlen, -1 means error
    virtual ~LIStream() {}

    int ReadLine(char *buffer, int n)
    {
        //Let's ignore \r completely
        int i, nr;
        i=0;
        do nr=Read(&buffer[i++],1); while(nr==1 && i<n && buffer[i-1]!='\n');
        buffer[i-1]=0;
        if(nr==0 && i==1) return -1;
        return i;
    }

    // ReadWord, read until the first occurrance of `delim' or length n
    // pd    : pointer to return the punctuation. (it get dumped)
    // RETURN: save as read()
    int ReadWord(char *buffer, int n, const char *punct=0,
                char *pd=0, const char *delim=" \n\t")
    {
        int i, nr;
        char ch;

        //Skipping leading deliminators
        do nr=Read(&ch,1);
        while(nr==1 && delim!=0 && strchr(delim, ch));
        if(nr!=1) return -1;

        buffer[0]=ch;
        i=1;
        if(punct!=0 && strchr(punct, ch))
        {
            if(pd!=0) *pd=ch;
            i--;
        }
        else for(;i<n;i++)
        {
            nr=Read(buffer+i,1);
            if(nr==-1) return -1;
            if( (nr==0) ||
               ((delim!=0) && strchr(delim, buffer[i])) ||
               ((punct!=0) && strchr(punct, buffer[i])) )
            {
                if(pd!=0) *pd=(nr==0)?0:buffer[i];
                break;
            }
        }
        buffer[i]='\0'; //terminate the string
        return i;
    }

    LIStream & operator >(unsigned &a)
    { Read((char *)&a, sizeof(a)); return *this; }
//    LIStream & operator >=(unsigned &a)
//    {
//#ifdef __LITTLE_ENDIAN
//        char c[sizeof(unsigned)];
//        Read(c, sizeof(unsigned));
//        Reverse((char *)&a, c, sizeof(unsigned));
//#else
//        Read((char *)&a, sizeof(a));
//#endif
//        return *this;
//    }
    LIStream & operator >(int &a)
    { Read((char *)&a, sizeof(a)); return *this; }
//    LIStream & operator >=(int &a)
//    {
//#ifdef __LITTLE_ENDIAN
//        char c[sizeof(int)];
//        Read(c, sizeof(int));
//        Reverse((char *)&a, c, sizeof(int));
//#else
//        Read((char *)&a, sizeof(a));
//#endif
//        return *this;
//    }
    LIStream & operator >(double &a)
    { Read((char *)&a, sizeof(a)); return *this; }
//    LIStream & operator >=(double &a)
//    {
//#ifdef __LITTLE_ENDIAN
//        char c[sizeof(double)];
//        Read(c, sizeof(double));
//        Reverse((char *)&a, c, sizeof(double));
//#else
//        Read((char *)&a, sizeof(a));
//#endif
//        return *this;
//    }
    LIStream & operator >(char &a)
    { Read((char *)&a, sizeof(a)); return *this; }
//new formatted input function
    LIStream &operator >>(int &a);
    LIStream &operator >>(double &a);
    LIStream &operator >>(char *a);
};


class LOStream
{
public:
//    enum { BUFLEN=LIStream::BUFLEN };
private:
//    static char *outbuf;
public:
    //Functions
    virtual int Write(const char *data, int size)=0; //$FIXME$ throw
    // The lock function
    virtual void Lock() {}
    virtual ~LOStream() {}

    //Binary outputs, using operator <
    LOStream & operator <(unsigned a)
    { Write((const char *)&a,sizeof(a)); return *this; }
    LOStream & operator <(int a)
    { Write((const char *)&a,sizeof(a)); return *this; }
    LOStream & operator <(double a)
    { Write((const char *)&a,sizeof(a)); return *this; }
    LOStream & operator <(char a)
    { Write((const char *)&a,sizeof(a)); return *this; }
    LOStream & operator <(const char *s)
    {
        int l=strlen(s);
        Write(s,l);
        return *this;
    }

    //Write with error check and retry, nretry=-1 means retry infinitely
    //udelay is delay in usec, default to 100 seconds
    int SafeWrite(const char *data, int size,
                  int nretry=-1, long udelay=100000000L);
    
    void FmtOutput(long l, unsigned format=0, int radix=10);
    void FmtOutput(double d/*,unsigned format=0, int radix=10*/);
public:
    enum Formats
    {
        NoPrefix=1,
        ZeroFill=2,
        Truncate=4,
        Hex=8,
        Oct=0x10,
        Bin=0x20,
        Unsigned=0x40,
        FillSign=0x80
        //Persist=0x8000,
    };
    int Printf(const char *fmt...)
    {
        int l;
        va_list ap;
        va_start(ap,fmt);
        l=vPrintf(fmt,ap);
        va_end(ap);
        return l;
    }

    int vPrintf(const char *fmt, va_list ap);
    //Text outputs, using operator <<
    //----bool----
    LOStream & operator <<(bool b)
    {
        const char strue[]="(true)";
        const char sfalse[]="(false)";
        if(b) Write(strue, sizeof(strue)-1);
        else Write(sfalse, sizeof(sfalse)-1);
        return *this;
    }
    //----char----
    LOStream & operator <<(char c)
    {
        Write(&c, sizeof(c));
        return *this;
    }
    //LOStream & operator <<(signed char c){ return operator<<((char)c);}
    //LOStream & operator <<(unsigned char c){ return operator<<((char)c);}
    //----short----
    //----int----
    LOStream & operator <<(unsigned u)
    {
        FmtOutput((long)u,Unsigned);
        return *this;
    }
    LOStream & operator <<(unsigned long u)
    {
        FmtOutput((long)u,Unsigned);
        return *this;
    }
    //----float----(is there any one who still use this?)
    //----double----
    LOStream & operator <<(double d)
    {
        FmtOutput(d);
        return *this;
    }
    LOStream & operator <<(int i)
    {
        FmtOutput((long)i);
        return *this;
    }
    LOStream & operator <<(long i)
    {
        FmtOutput(i);
        return *this;
    }
    LOStream & operator <<(const char *s)
    {
        if(s)
        {
            int l=strlen(s);
            Write(s, l);
        }
        return *this;
    }
    LOStream & operator <<(void *p);
};

//-----------------------------------------------------------------------
//Begin formated output (no input now)
//-----------------------------------------------------------------------

//Maybe I should put it inside LOStream,
// but to save typing of user, he just need _IO << LFmt(ldkfjsd)
/*
class LFmt
{
};


class LOStream : public LOStr
{
private:
    int radix; //A number from 2 to 36(all the number + all the digit)
    //showpos:  0--put sign only necessary(-), 1--always put sign
    //showbase: 0--don't put prefix to radix, 1--put prefix,
    //          format of prefix: 2:(0b), 8:(0), 10:(null) 16:(0x)
    //             others: n:(0nr)(e.g. 013r3AC==3*13^2+10*13+12)
    //             Notice here for octals, zerofill may get confused
    //fixed:
    //left/right/center:
    //scientific:
    enum
    {
        ShowBase=0x0010;
        Align=0x0003;
    };
    enum
    {
        AlignLeft=0x0001;
        AlignRight=0x0002;
        AlignCenter=0x0003;
    };
    unsigned flag;

    char fill;//default ' '
    int precision;//default 6
    int width;//default 0 (as many as you want)

    char toChar(unsigned d) const
    {
        {
          '0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F',
            'G','H','I','J','K','L',
        };
        return tc[(d%radix)];
    }
    char *writeUnsigned(unsigned val)
    {
        char *ptr=&outbuf[BUFLEN];
        do *--ptr=toChar(u); while((u/=radix)!=0 && ptr>outbuf);
        if(u) *outbuf='$';//means incorrect.
        return ptr;
    }
    void doWrite(char *ptr)
    {
        ASSERT(ptr>=outbuf && ptr<outbuf+BUFLEN);
        Write(ptr,outbuf+BUFLEN-ptr);
    }
};
*/

// This is the only place a multi-inherence is induced, we can redefine this
// class to eliminate multi-inherence.
class LStream : public LOStream, public LIStream {};

class LFDStream : public LStream
{
    int ifd, ofd;
    int oldifd, oldofd;
public:
    LFDStream():ifd(-1),ofd(-1){}
    LFDStream(int If, int Of):ifd(If),ofd(Of){}
    virtual int Write(const char *data, int size)
    { return ofd < 0?0:write(ofd, data, size); }
    virtual int Read(char *buff, int size)
    { return ifd < 0?0:read(ifd, buff, size); }
    virtual void Lock()
    {
//        lock(ofd);
//        if(ifd!=ofd)
//            lock(ifd);
    }
    void Flush(){ fsync(ofd); }
    void Redirect(int i, int o){ oldifd=ifd; oldofd=ofd; ifd=i, ofd=o; }
    void RedirectOut(int o){ oldofd=ofd; ofd=o; }
    void Restore() { ifd=oldifd; ofd=oldofd; }
};

extern LFDStream _IO;
extern LFDStream _Error;
extern LFDStream _Null;


class LSys
{
#ifdef _DEBUG
private:
    static int m_Init;
public:
    static const char *sDebug, *sDebugEnd;
    static void sig_handler(int, siginfo_t *sig, void *);
    LSys(){ }; //Debugger initialization
    ~LSys();
#endif

public:
    static const char *sInfo, *sWarning, *sError, *sFatal,
        *sInfoEnd, *sWarningEnd, *sErrorEnd, *sFatalEnd;

    static LOStream *_pSys;

    static LOStream *Redirect(LOStream *os)
    {
        LOStream *p=_pSys;
        _pSys=os;
        return p;
    }
    static void Restore()
    {
        _pSys=&_Error;
        sInfo=GRN "[I]" NOR" "; sInfoEnd=GRN "[i]" NOR"\n";
        sWarning=YEL "[W]" NOR" "; sWarningEnd=YEL "[w]" NOR"\n";
        sError=RED "[E]" NOR" "; sErrorEnd=HIR "[e]" NOR"\n";
        sError=RED "[E]" NOR" "; sErrorEnd=HIR "[e]" NOR"\n";
        sFatal=HIR "[F]" NOR" "; sFatalEnd=HIR "[f]" NOR"\n";

#ifdef _DEBUG
        sDebug=HIC"[D]"NOR" "; sDebugEnd=HIC"[d]"NOR"\n";
#endif
			}
			static void Abort()
			{
#ifdef _DEBUG
			m_Init=0; //Disabling memory leak checking
#endif
			exit(1);
			}
			};

			//extern LSys _Sys;
#define INFO(d) (void)(*LSys::_pSys << LSys::sInfo << d << LSys::sInfoEnd)
#define INFO_(d) (void)(*LSys::_pSys << LSys::sInfo << d)
#define _INFO_(d) (void)(*LSys::_pSys << d)
#define _INFO(d) (void)(*LSys::_pSys << d << LSys::sInfoEnd)
#define WARNING(d) (void)(*LSys::_pSys<<LSys::sWarning<< d <<LSys::sWarningEnd)
#define WARNING_(d) (void)(*LSys::_pSys << LSys::sWarning << d)
#define _WARNING(d) (void)(*LSys::_pSys << d << LSys::sWarningEnd)
#define ERROR(d) (void)(*LSys::_pSys << LSys::sError << d << LSys::sErrorEnd)
#define ERROR_(d) (void)(*LSys::_pSys << LSys::sError << d <<)
#define _ERROR(d) (void)(*LSys::_pSys << d << LSys::sErrorEnd)
#define FATAL(d) (void)((*LSys::_pSys << LSys::sFatal << d <<                 \
			LSys::sFatalEnd), exit(1))
#define SYSERROR(d) (void)(*LSys::_pSys << LSys::sError << d <<               \
		" ("<<strerror(errno)<<")"<<LSys::sErrorEnd)
#define SYSFATAL(d) (void)((*LSys::_pSys << LSys::sError << d <<              \
			" ("<<strerror(errno)<<")"<<LSys::sErrorEnd), exit(1))

inline void INFO_Dump(const char *fmt...)
{
	va_list ap;
	va_start(ap,fmt);
	*LSys::_pSys << LSys::sInfo;
	LSys::_pSys->vPrintf(fmt,ap);
	*LSys::_pSys << LSys::sInfoEnd;
}
inline void WARNING_Dump(const char *fmt...)
{
	va_list ap;
	va_start(ap,fmt);
	*LSys::_pSys << LSys::sWarning;
	LSys::_pSys->vPrintf(fmt,ap);
	*LSys::_pSys << LSys::sWarningEnd;
}
inline void ERROR_Dump(const char *fmt...)
{
	va_list ap;
	va_start(ap,fmt);
	*LSys::_pSys << LSys::sError;
	LSys::_pSys->vPrintf(fmt,ap);
	*LSys::_pSys << LSys::sErrorEnd;
}
inline void FATAL_Dump(const char *fmt...)
{
	va_list ap;
	va_start(ap,fmt);
	*LSys::_pSys << LSys::sFatal;
	LSys::_pSys->vPrintf(fmt,ap);
	*LSys::_pSys << LSys::sFatalEnd;
	LSys::Abort();
}

#define INFO_Printf LSys::_pSys->Printf
#define INFO_Flush LSys::_pSys->Flush

#ifdef _DEBUG
//This named object is not used, just to invoke the initialization process
static class LSys _Internal_Debugging_Object; 

#define DUMP(d) (void)(*LSys::_pSys << LSys::sDebug << d << LSys::sDebugEnd)
#define DUMP_(d) (void)(*LSys::_pSys << LSys::sDebug << d)
#define _DUMP(d) (void)(*LSys::_pSys << d << LSys::sDebugEnd)
inline void DEBUG_Dump(const char *fmt...)
{
	va_list ap;
	va_start(ap,fmt);
	*LSys::_pSys << LSys::sDebug;
	LSys::_pSys->vPrintf(fmt,ap);
	*LSys::_pSys << LSys::sDebugEnd;
}

#define ASSERT(p) if(p) ; else DUMP( "Assertion Failed (file "                \
	<< __FILE__ << ", line " << __LINE__ << ')'),LSys::Abort()
#define VERIFY(p) ASSERT(p)
#define ASSERTMSG(p, strm) if(p); else DUMP( "Assertion Failed (file "        \
	<< __FILE__ << ", line "<< __LINE__ << "):" << strm), LSys::Abort()
#define NOTREACHED() DUMP("Illegal branch (file "                             \
	<< __FILE__ << ", line "<<__LINE__ << ')'), LSys::Abort()
#define NOTREACHEDMSG(d) DUMP("Illegal branch (file "                         \
	<< __FILE__ << ", line "<<__LINE__ << "):" << d), LSys::Abort()

#else //no _DEBUG

#define DUMP(d) ((void)0)
#define DUMP_(d) ((void)0)
#define _DUMP(d) ((void)0)
#define ASSERT(p) ((void)0)
#define ASSERTMSG(p, strm) ((void)0)
#define NOTREACHED() FATAL("Internal error: Illegal branch")
#define NOTREACHEDMSG(d) FATAL("Internal error: " << d)
#define VERIFY(p) (p)

#endif //_DEBUG

#define NOTIMPLEMENTED() NOTREACHEDMSG("Feature not impleneted")

//============================================================
// String functions
const unsigned _cU=0x0001;  //Upper case
const unsigned _cL=0x0002;  //Lower case
const unsigned _cN=0x0004;  //Numeral (digit)
const unsigned _cS=0x0008;  //Space
const unsigned _cP=0x0010;  //Printable
const unsigned _cX=0x0080;  //Hexdecimal digit

extern const unsigned _CharTable[];

inline bool IsDigit(char ch){  return _CharTable[(unsigned char)ch]&_cN; }
inline bool IsHexDigit(char ch){  return _CharTable[(unsigned char)ch]&_cX; }
inline bool IsSpace(char ch){  return _CharTable[(unsigned char)ch]&_cS; }
inline bool IsUpper(char ch){  return _CharTable[(unsigned char)ch]&_cU; }
inline bool IsLower(char ch){  return _CharTable[(unsigned char)ch]&_cL; }
inline bool IsLetter(char ch){ return _CharTable[(unsigned char)ch]&(_cL|_cU);}
inline bool IsPrintable(char ch){  return _CharTable[(unsigned char)ch]&_cP; }

inline char _ToUpper(char c){ ASSERT(IsLower(c)); return c&0xdf;}
inline char _ToLower(char c){ ASSERT(IsUpper(c)); return c|0x20;}
inline char ToUpper(char c){  return IsLower(c) ? _ToUpper(c):c; }
inline char ToLower(char c){  return IsUpper(c) ? _ToLower(c):c; }

inline char HexToChar(unsigned h)
{
	static const char htc[]=
	{
		'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'
	};
	return htc[h&0x000F];
}

inline unsigned CharToHex(char ch)
{
	ASSERT(IsHexDigit(ch));
	return (ch>'9')? _ToUpper(ch)-'A'+10 : ch-'0';
}

// Aliases of string.h functions
// The functions below (1) has ASSERT, (2) works fine with StrStream etc.
#ifdef _DEBUG
inline void *StrMemCpy(void *dest, const void *src, unsigned len)
{ ASSERT(dest!=0 && src!=0); return memcpy(dest, src, len);}
inline void *StrMemMove(void *dest, const void *src, unsigned len)
{ ASSERT(dest!=0 && src!=0); return memmove(dest, src, len);}
inline char *StrCpy(char *dest, const char *src)
{ ASSERT(dest!=0 && src!=0); return strcpy(dest, src);}
inline char *StrNCpy(char *dest, const char *src, unsigned len)
{ ASSERT(dest!=0 && src!=0); return strncpy(dest, src, len);}
inline char *StrCat(char *dest, const char *src)
{ ASSERT(dest!=0 && src!=0); return strcat(dest, src);}
inline char *StrNCat(char *dest, const char *src, unsigned len)
{ ASSERT(dest!=0 && src!=0); return strncat(dest, src, len);}
inline int StrMemCmp(const void *m1, const void *m2, unsigned len)
{ ASSERT(m1!=0 && m2!=0);  return memcmp(m1, m2, len);}
inline int StrCmp(const char *m1, const char *m2)
{ ASSERT(m1!=0 && m2!=0);  return strcmp(m1, m2);}
inline int StrCaseCmp(const char *m1, const char *m2)
{ ASSERT(m1!=0 && m2!=0);  return strcasecmp(m1, m2);}
inline int StrNCmp(const char *m1, const char *m2, unsigned len)
{ ASSERT(m1!=0 && m2!=0); return strncmp(m1, m2, len);}
inline int StrNCaseCmp(const char *m1, const char *m2, unsigned len)
{ ASSERT(m1!=0 && m2!=0); return strncasecmp(m1, m2, len);}
inline void *StrMemChr(const void *s, char c, unsigned l)
{ ASSERT(s!=0); return memchr(s, c, l); }
inline char *StrChr(const char *s, char c)
{ ASSERT(s!=0); return strchr(s, c); }
inline const char *StrStr(const char *s1, const char *s2)
{ ASSERT(s1!=0 && s2!=0); return strstr(s1, s2); }
inline char *StrStr(char *s1, const char *s2)
{ ASSERT(s1!=0 && s2!=0); return strstr(s1, s2); }

inline void *StrMemSet(void *s, char b, unsigned l)
{ ASSERT(s!=0);  return memset(s, b, l); }
inline void StrBZero(void *s, int n)
{ ASSERT(s!=0); memset(s,0,n); }
inline unsigned StrLen(const char *s)
{  ASSERT(s!=0); return strlen(s); }
#else //Release
#define StrMemCpy(d,s,l) memcpy(d,s,l)
#define StrMemMove(d,s,l) memmove(d,s,l)
#define StrCpy(d,s) strcpy(d,s)
#define StrNCpy(d,s,l) strncpy(d,s,l)
#define StrCat(d,s) strcat(d,s)
#define StrNCat(d,s,l) strncat(d,s,l)
#define StrMemCmp(d,s,l) memcpy(d,s,l)
#define StrCmp(d,s) strcmp(d,s)
#define StrCaseCmp(d,s) strcasecmp(d,s)
#define StrNCmp(d,s,l) strncmp(d,s,l)
#define StrNCaseCmp(d,s,l) strncasecmp(d,s,l)
#define StrMemChr(s,c,l) memchr(s,c,l)
#define StrChr(s,c) strchr(s,c)
#define StrStr(a,b) strstr(a,b)
#define StrMemSet(s,b,l) memset(s,b,l)
#define StrBZero(s,l) ((void)memset(s,0,l))
#define StrLen(s) strlen(s)
#endif //_DEBUG
//StrCaseStr(): strstr while case insensitive
const char *StrCaseStr(const char *s1, const char *s2);
//StrSkipSpaces(): Return the first non-space charactor from str
const char *StrSkipSpaces(const char *str);
inline char *StrSkipSpaces(char *str)
{ return (char *)StrSkipSpaces((const char *)str); }
//StrSnipSpaces(): Put \0 to truncate string, make sure no trailing spaces
void StrSnipSpaces(char *str);


// string formatting functions.
bool StrScanInt(const char *s, int &i);
bool StrScanLong(const char *s, long &i);
bool StrScanDouble(const char *s, double &d);
int StrPrintf(char *str, const char *fmt...);
int StrNPrintf(char *str, int s, const char *fmt...);
//#ifdef __GNUC__
//int StrScanf(const char *str, const char *ftm...);
//#else
//extern "C" int sscanf(const char *str, const char *format...);
#define StrScanf sscanf
//#endif

#else //not __cplusplus

//#define Abs(a) ((a)>=0?(a):-(a))
//#define Max(a, b) ((a)>(b)?(a):(b))
//#define Min(a, b) ((a)<(b)?(a):(b))
//#define Square(a) ((a)*(a))
//#define Cube(a) ((a)*(a)*(a))

//===========================================================
//Debugging for C (without C++ compiler)
static void Warning(const char* fmt, ...)
{
	va_list ap;
	va_start(ap,fmt);
	printf(BRED HIW"*Warning*"NOR": ");
	vprintf(fmt, ap);
	printf("\n");
}

static void Error(const char* fmt, ...)
{
	va_list ap;
	va_start(ap,fmt);
	printf(BRED HIW"*Error*"NOR": ");
	vprintf(fmt, ap);
	printf("\n");
	exit(1);
}

#ifdef _DEBUG

#include <stdarg.h>
static const char* DB_File;
static int DB_Line;
static void DB_RegLocus(const char* file, int line)
{
	DB_File=file;
	DB_Line=line;
}
static void DB_Assert(int c)
{
	if(!c) printf("Assertion Failed at %s:%d\n",DB_File,DB_Line);
}
static void DB_AssertMsg(int c, const char *fmt, ...)
{
	if(!c)
	{
		va_list ap;
		va_start(ap,fmt);
		printf("Assertion Failed at %s:%d--",DB_File,DB_Line);
		vprintf(fmt,ap);
		printf("\n");
		va_end(ap);
	}
}

#define ASSERT(p) DB_RegLocus(__FILE__,__LINE__),DB_Assert(p);
#define VERIFY(p) ASSERT(p)
#define ASSERTMSG DB_RegLocus(__FILE__,__LINE__),DB_AssertMsg
#define NOTREACHED() ASSERTMSG(0,"Not Reachable")

#else //no _DEBUG
#define DUMP(d) ((void)0)
#define ASSERT(p) ((void)0)
#define ASSERTMSG 1?(void)0:(void)
#define NOTREACHED() ((void)0)
#define NOTREACHEDMSG(d) ((void)0)
#define VERIFY(p) ((void)(p))
#endif //_DEBUG


#endif //__cplusplus

//---------
// Debug for memory leaks in Pool management
#ifndef _DEBUG
inline void *MAlloc(size_t size)
{
	//    void *ptr;
	//    ptr=malloc(size);
	//    INFO("Alloc("<<size<<")="<<ptr);
	//    return ptr;
	return malloc(size);
}
inline void *Realloc(void *ptr, size_t size)
{
	//    void *ptr2;
	//    ptr2=realloc(ptr, size);
	//    INFO("Realloc("<<ptr<<','<<size<<")="<<ptr2);
	//    return ptr2;
	return realloc(ptr, size);
}

inline void Free(void *ptr)
{
	//    INFO("Free()="<<ptr);
	free(ptr);
}
#else
enum MemType
{
	M_Array,
	M_Obj,
	M_Alloc
};
void *MemAlloc(size_t size, const char *file, unsigned line, enum MemType mt);
void *MemRealloc(void *p,size_t size,
		const char *file, unsigned line, enum MemType mt);
void MemFree(void *p, const char *f, unsigned l, enum MemType mt);
void MemDump(const char *f, unsigned l);


#define MAlloc(s) MemAlloc(s,__FILE__,__LINE__, M_Alloc)
#define Realloc(p,s) MemRealloc(p,s,__FILE__,__LINE__, M_Alloc)
#define Free(p) MemFree(p,__FILE__,__LINE__, M_Alloc)
#define DUMP_MEM() MemDump(__FILE__, __LINE__)

#endif //_DEBUG

// Array allocation helpers
inline char **ArrAlloc2(unsigned rsize, int dim1, int dim2)
{
	int i;
	char **p=(char **)MAlloc((sizeof(void *)+rsize*dim2)*dim1);
	p[0]=(char *)p+dim1*sizeof(void *);

	for(i=1;i<dim1;i++)
		p[i]=p[i-1]+rsize;
	return p;
}
inline char **ReArrAlloc2(void *ptr, unsigned rsize, int dim1, int dim2)
{
	int i;
	char **p=(char **)Realloc(ptr, (sizeof(void *)+rsize*dim2)*dim1);
	p[0]=(char *)p+dim1*sizeof(void *);

	for(i=1;i<dim1;i++)
		p[i]=p[i-1]+rsize*dim2;
	return p;
}
inline char *StrDup(const char *s)
{
	//    char *p=new char[StrLen(s)];
	char *p=(char *)MAlloc(StrLen(s)+1);
	return StrCpy(p,s);
}

inline bool Equal(double a, double b, double torlerate=1e-13)
{
	return Abs(a-b) < torlerate;
}

#endif // _GENERAL_H

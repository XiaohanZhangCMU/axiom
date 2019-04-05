// Last Modified : Sat Sep 22 13:50:29 2001

#include "general.h"

/////////////////////////////////////////////////////////////
//StreamBase variables and functions
LOStream *LSys::_pSys=&_Error;

const char *LSys::sInfo= GRN "[I]" NOR " ",*LSys::sInfoEnd="\n",//" "GRN"[i]"NOR"\n",
    *LSys::sWarning= YEL "[W]" NOR " ",
    *LSys::sWarningEnd="\n",//" "YEL"[w]"NOR"\n",
    *LSys::sError= RED "[E]" NOR " ",*LSys::sErrorEnd="\n",//" "RED"[e]"NOR"\n",
    *LSys::sFatal= HIR "[F]" NOR " ",*LSys::sFatalEnd="\n";//" "HIR"[f]"NOR"\n";
#ifdef _DEBUG
char *LSys::sDebug= HIC "[D]" NOR " ",*LSys::sDebugEnd="\n";//" "HIC"[d]"NOR"\n";
#endif

LFDStream _IO(0,1);
LFDStream _Error(-1,2);
LFDStream _Null;

int LOStream::SafeWrite(const char *data, int size, int nretry, long udel)
{
    int retval, retry, nwritten;
    
    retry=nwritten=0;
 l_retry:
    retval=Write(data+nwritten, size-nwritten);
    if(retval<0 && (nretry==-1 || ++retry < nretry)) {//Error, need retry
        SYSERROR("SafeWrite: Failed to write "<<size-nwritten<<
                 " bytes to file, retrying...");
        udelay(udel);
        goto l_retry;
    }
    else if(retval < size-nwritten){ //Not all written
        nwritten+=retval;
        if(nretry==-1 || ++retry < nretry) {
            ERROR("SafeWrite: Written "<<nwritten<<" out of "<<size<<
                  " bytes, trying to write remaining...");
            goto l_retry;
        }
    }
    return retval;
}
/////////////////////////////////////////////////////////////
//The CRCTable
unsigned UpdateCRC(unsigned crc, unsigned char *buf, int len)
{
    static unsigned _CRCTable[256]=
    {
        0x00000000, 0x77073096, 0xEE0E612C, 0x990951BA,
        0x076DC419, 0x706AF48F, 0xE963A535, 0x9E6495A3,
        0x0EDB8832, 0x79DCB8A4, 0xE0D5E91E, 0x97D2D988,
        0x09B64C2B, 0x7EB17CBD, 0xE7B82D07, 0x90BF1D91,
        0x1DB71064, 0x6AB020F2, 0xF3B97148, 0x84BE41DE,
        0x1ADAD47D, 0x6DDDE4EB, 0xF4D4B551, 0x83D385C7,
        0x136C9856, 0x646BA8C0, 0xFD62F97A, 0x8A65C9EC,
        0x14015C4F, 0x63066CD9, 0xFA0F3D63, 0x8D080DF5,
        0x3B6E20C8, 0x4C69105E, 0xD56041E4, 0xA2677172,
        0x3C03E4D1, 0x4B04D447, 0xD20D85FD, 0xA50AB56B,
        0x35B5A8FA, 0x42B2986C, 0xDBBBC9D6, 0xACBCF940,
        0x32D86CE3, 0x45DF5C75, 0xDCD60DCF, 0xABD13D59,
        0x26D930AC, 0x51DE003A, 0xC8D75180, 0xBFD06116,
        0x21B4F4B5, 0x56B3C423, 0xCFBA9599, 0xB8BDA50F,
        0x2802B89E, 0x5F058808, 0xC60CD9B2, 0xB10BE924,
        0x2F6F7C87, 0x58684C11, 0xC1611DAB, 0xB6662D3D,
        0x76DC4190, 0x01DB7106, 0x98D220BC, 0xEFD5102A,
        0x71B18589, 0x06B6B51F, 0x9FBFE4A5, 0xE8B8D433,
        0x7807C9A2, 0x0F00F934, 0x9609A88E, 0xE10E9818,
        0x7F6A0DBB, 0x086D3D2D, 0x91646C97, 0xE6635C01,
        0x6B6B51F4, 0x1C6C6162, 0x856530D8, 0xF262004E,
        0x6C0695ED, 0x1B01A57B, 0x8208F4C1, 0xF50FC457,
        0x65B0D9C6, 0x12B7E950, 0x8BBEB8EA, 0xFCB9887C,
        0x62DD1DDF, 0x15DA2D49, 0x8CD37CF3, 0xFBD44C65,
        0x4DB26158, 0x3AB551CE, 0xA3BC0074, 0xD4BB30E2,
        0x4ADFA541, 0x3DD895D7, 0xA4D1C46D, 0xD3D6F4FB,
        0x4369E96A, 0x346ED9FC, 0xAD678846, 0xDA60B8D0,
        0x44042D73, 0x33031DE5, 0xAA0A4C5F, 0xDD0D7CC9,
        0x5005713C, 0x270241AA, 0xBE0B1010, 0xC90C2086,
        0x5768B525, 0x206F85B3, 0xB966D409, 0xCE61E49F,
        0x5EDEF90E, 0x29D9C998, 0xB0D09822, 0xC7D7A8B4,
        0x59B33D17, 0x2EB40D81, 0xB7BD5C3B, 0xC0BA6CAD,
        0xEDB88320, 0x9ABFB3B6, 0x03B6E20C, 0x74B1D29A,
        0xEAD54739, 0x9DD277AF, 0x04DB2615, 0x73DC1683,
        0xE3630B12, 0x94643B84, 0x0D6D6A3E, 0x7A6A5AA8,
        0xE40ECF0B, 0x9309FF9D, 0x0A00AE27, 0x7D079EB1,
        0xF00F9344, 0x8708A3D2, 0x1E01F268, 0x6906C2FE,
        0xF762575D, 0x806567CB, 0x196C3671, 0x6E6B06E7,
        0xFED41B76, 0x89D32BE0, 0x10DA7A5A, 0x67DD4ACC,
        0xF9B9DF6F, 0x8EBEEFF9, 0x17B7BE43, 0x60B08ED5,
        0xD6D6A3E8, 0xA1D1937E, 0x38D8C2C4, 0x4FDFF252,
        0xD1BB67F1, 0xA6BC5767, 0x3FB506DD, 0x48B2364B,
        0xD80D2BDA, 0xAF0A1B4C, 0x36034AF6, 0x41047A60,
        0xDF60EFC3, 0xA867DF55, 0x316E8EEF, 0x4669BE79,
        0xCB61B38C, 0xBC66831A, 0x256FD2A0, 0x5268E236,
        0xCC0C7795, 0xBB0B4703, 0x220216B9, 0x5505262F,
        0xC5BA3BBE, 0xB2BD0B28, 0x2BB45A92, 0x5CB36A04,
        0xC2D7FFA7, 0xB5D0CF31, 0x2CD99E8B, 0x5BDEAE1D,
        0x9B64C2B0, 0xEC63F226, 0x756AA39C, 0x026D930A,
        0x9C0906A9, 0xEB0E363F, 0x72076785, 0x05005713,
        0x95BF4A82, 0xE2B87A14, 0x7BB12BAE, 0x0CB61B38,
        0x92D28E9B, 0xE5D5BE0D, 0x7CDCEFB7, 0x0BDBDF21,
        0x86D3D2D4, 0xF1D4E242, 0x68DDB3F8, 0x1FDA836E,
        0x81BE16CD, 0xF6B9265B, 0x6FB077E1, 0x18B74777,
        0x88085AE6, 0xFF0F6A70, 0x66063BCA, 0x11010B5C,
        0x8F659EFF, 0xF862AE69, 0x616BFFD3, 0x166CCF45,
        0xA00AE278, 0xD70DD2EE, 0x4E048354, 0x3903B3C2,
        0xA7672661, 0xD06016F7, 0x4969474D, 0x3E6E77DB,
        0xAED16A4A, 0xD9D65ADC, 0x40DF0B66, 0x37D83BF0,
        0xA9BCAE53, 0xDEBB9EC5, 0x47B2CF7F, 0x30B5FFE9,
        0xBDBDF21C, 0xCABAC28A, 0x53B39330, 0x24B4A3A6,
        0xBAD03605, 0xCDD70693, 0x54DE5729, 0x23D967BF,
        0xB3667A2E, 0xC4614AB8, 0x5D681B02, 0x2A6F2B94,
        0xB40BBE37, 0xC30C8EA1, 0x5A05DF1B, 0x2D02EF8D,
    };
    unsigned c=crc^0xFFFFFFFF;
    int n;
    
    for (n = 0; n < len; n++)
        c = _CRCTable[(c ^ buf[n]) & 0xff] ^ (c >> 8);
    return c ^ 0xFFFFFFFF;
}
/* The above program is made by the following:
static void MakeCRCTable(void)
{
    unsigned c;
    unsigned n, k;
    for (n = 0; n < 256; n++)
    {
        c = n;
        for (k = 0; k < 8; k++)
            if(c & 1) c = 0xEDB88320 ^ (c >> 1);
            else c >>= 1;
        _CRCTable[n] = c;
    }
}
void PrintCRCTable()
{
    MakeCRCTable();
    _IO << "\nstatic unsigned _CRCTable[256]=\n{";
    for(int i=0;i<256;i++)
    {
        if(i%4==0) _IO << "\n    ";
        _IO << LOStream::Hex(_CRCTable[i]) << ", ";
    }
    _IO << "\n};\n";
}
*/

/////////////////////////////////////////////////////////////
// The Stream Base
/*
    class Hex
    {
        unsigned h;
        unsigned fmt;
        int param;
        friend LOStream;
    public:
        Hex(unsigned d,unsigned format=0,int parameter=-1)
            :h(d),fmt(format),param(parameter)
        {
            if(fmt&No0x)
            {
                if(param==-1) param=8;
            }
        }
        friend LOStream & operator <<(LOStream &os, Hex h)
        {
            //We know that the largest unsigned int is 0xFFFFFFFF
            char buf[10];
            int i=sizeof(buf);
            do
            {
                buf[--i]=HexToChar(h.h);
                h.h >>=4;
            } while(h.h);
            buf[--i]='x', buf[--i]='0';
            os.Write(buf+i, sizeof(buf)-i);
            return os;
        }
    };
*/

LOStream & LOStream::operator <<(void *p)
{
    (*this << '^').FmtOutput((long)p,Hex|NoPrefix|ZeroFill);
    return *this;
}

void LOStream::FmtOutput(long l, unsigned format, int /*radix*/)
{
    if((format & Hex) && (format & ZeroFill))
        Printf("%08X", l);
    else Printf("%d",l);
//    if(!(format&Unsigned))
//        if(l<0)
//        {
//            Write("-",1);
//            l=-l;
//        }
//        else if(format&FillSign)
//            Write("+",1);
//    //Output the radix prefix
//    if((format&Hex)&&(radix==10))
//        radix=16;
//    if(radix < 2 || radix > 36)
//        WARNING("(Stream):radix out of range\n");
//    if((!(format&NoPrefix))&&radix!=10)
//    {
//        outbuf[0]='0';
//        if(radix==2)
//            outbuf[1]='b';
//        else if(radix==16)
//            outbuf[1]='x';
//        else outbuf[1]='@';
//        Write(outbuf,2);
//    }
//    //Output the unsigned
//    {
//        char *ptr=&outbuf[BUFLEN];
//        do *--ptr=toChar(l%radix); while((l/=radix)!=0 && ptr>outbuf);
//        if(l) *outbuf='$';//means incorrect.
//        Write(ptr,outbuf+BUFLEN-ptr);
//    }
}

#include <stdio.h>
void LOStream::FmtOutput(double d/*,unsigned format=0, int radix*/)
{
    Printf("%.14g",d);
    /*
    double e=1e-15; //Try to display at most 16 digits
    unsigned n;
    if(d<0)
    {
        *this << '-';
        d=-d;
    }
    n=(unsigned)d;
    *this << n << '.';
    d-=n;
    do
    {
        d*=10; e*=10;
        n=(unsigned)(d+e);
        *this << (unsigned)n;
        d-=n;
    }
    while(d > e);
    */
}
int LOStream::vPrintf(const char *fmt, va_list ap)
{
    char outbuf[LIStream::BUFLEN];
    int l;
    l=vsprintf(outbuf, fmt, ap);
    Write(outbuf, l);
    return l;
}
//int LIStream::Scanf(const char *fmt...)
//{
//    
//}
LIStream &LIStream::operator >>(int &a)
{
    char inbuf[100];
    ReadWord(inbuf,BUFLEN);
    StrScanInt(inbuf,a);
    return *this;
}

LIStream &LIStream::operator >>(double &a)
{
    char inbuf[100];
    ReadWord(inbuf,BUFLEN);
    StrScanDouble(inbuf,a);
    return *this;
}

LIStream &LIStream::operator >>(char *a)
{ ReadWord(a,BUFLEN); return *this;}

/////////////////////////////////////////////////////////////
// The Strings, class
const char *StrCaseStr(const char *s1, const char *s2) 
{
    ASSERT(s1!=NULL && s2!=NULL);
    int sl1, sl2, i;
    sl1=StrLen(s1); sl2=StrLen(s2);
    
    for(i=0;i<=sl1-sl2;i++) if(StrNCaseCmp(s1+i, s2, sl2)==0) return s1+i;
    return NULL;
}
const char *StrSkipSpaces(const char *str)
{
    ASSERT(str!=NULL);
    const char *ret=str;
    while(IsSpace(*ret)) ret++;
    return ret;
}
void StrSnipSpaces(char *str)
{
    char *ret;
    for(ret=str+StrLen(str)-1; ret>=str && IsSpace(*ret); ret--);
    ret[1]='\0';
}

bool StrScanInt(const char *s, int &i)
{
    return sscanf(s,"%d",&i)==1;
}    
bool StrScanLong(const char *s, long &i)
{
    return sscanf(s,"%ld",&i)==1;
}    
bool StrScanDouble(const char *s, double &d)
{
    return sscanf(s,"%lf",&d)==1;    
}
#ifndef NO_VSNPRINTF
int StrNPrintf(char *str, int size, const char *fmt...)
{
    int ret;
    va_list ap;
    va_start(ap, fmt);
    ret=vsnprintf(str,size, fmt, ap);
    va_end(ap);
    return ret;
}
#endif
int StrPrintf(char *str, const char *fmt...)
{
    int ret;
    va_list ap;
    va_start(ap, fmt);
    ret=vsprintf(str, fmt, ap);
    va_end(ap);
    return ret;
}

//Following line extracted from <stdio.h>
//#ifdef __GNUC__
//extern "C" int vsscanf(const char *str, const char *format, va_list ap);
//
//int StrScanf(const char *str, const char *fmt...)
//{    
//    int ret;
//    va_list ap;
//    va_start(ap, fmt);
//    ret=vsscanf(str, fmt, ap);
//    va_end(ap);
//    return ret;
//}
//#endif

//CharTable contents:
//  Space (_cS) ^L,\t,\n,\r,(sp)
const unsigned _CharTable[256] =
{
    0/*`^@(nul)'*/, 0/*`^A(soh)'*/, 0/*`^B(stx)'*/, 0/*`^C(etx)'*/, 
    0/*`^D(eot)'*/, 0/*`^E(enq)'*/, 0/*`^F(ack)'*/, 0/*`^G(bel)'*/, 
    0/*`^H(bs)'*/, _cS|_cP/*`\t(ht)'*/, _cS|_cP/*`\n(nl)'*/, 0/*`^K(vt)'*/, 
    _cS/*`^L(np)'*/, _cS/*`\r(cr)'*/, 0/*`^N(so)'*/, 0/*`^O(si)'*/, 
    0/*`^P(dle)'*/, 0/*`^Q(dc1)'*/, 0/*`^R(dc2)'*/, 0/*`^S(dc3)'*/, 
    0/*`^T(dc4)'*/, 0/*`^U(nak)'*/, 0/*`^V(syn)'*/, 0/*`^W(etb)'*/, 
    0/*`^X(can)'*/, 0/*`^Y(em)'*/, 0/*`^Z(sum)'*/, 0/*`^[(esc)'*/, 
    0/*`^\(fs)'*/, 0/*`^](gs)'*/, 0/*`^^(rs)'*/, 0/*`^_(us)'*/, 
    _cS|_cP/*` (sp)'*/, _cP/*`!'*/, _cP/*`"'*/, _cP/*`#'*/, 
    _cP/*`$'*/, _cP/*`%'*/, _cP/*`&'*/, _cP/*`''*/, 
    _cP/*`('*/, _cP/*`)'*/, _cP/*`*'*/, _cP/*`+'*/, 
    _cP/*`,'*/, _cP/*`-'*/, _cP/*`.'*/, _cP/*`/'*/, 
    _cN|_cX|_cP/*`0'*/, _cN|_cX|_cP/*`1'*/,
    _cN|_cX|_cP/*`2'*/, _cN|_cX|_cP/*`3'*/, 
    _cN|_cX|_cP/*`4'*/, _cN|_cX|_cP/*`5'*/,
    _cN|_cX|_cP/*`6'*/, _cN|_cX|_cP/*`7'*/, 
    _cN|_cX|_cP/*`8'*/, _cN|_cX|_cP/*`9'*/, _cP/*`:'*/, _cP/*`;'*/, 
    _cP/*`<'*/, _cP/*`='*/, _cP/*`>'*/, _cP/*`?'*/, 
    _cP/*`@'*/, _cU|_cX|_cP/*`A'*/, _cU|_cX|_cP/*`B'*/, _cU|_cX|_cP/*`C'*/, 
    _cU|_cX|_cP/*`D'*/, _cU|_cX|_cP/*`E'*/,
    _cU|_cX|_cP/*`F'*/, _cU|_cP/*`G'*/, 
    _cU|_cP/*`H'*/, _cU|_cP/*`I'*/, _cU|_cP/*`J'*/, _cU|_cP/*`K'*/, 
    _cU|_cP/*`L'*/, _cU|_cP/*`M'*/, _cU|_cP/*`N'*/, _cU|_cP/*`O'*/, 
    _cU|_cP/*`P'*/, _cU|_cP/*`Q'*/, _cU|_cP/*`R'*/, _cU|_cP/*`S'*/, 
    _cU|_cP/*`T'*/, _cU|_cP/*`U'*/, _cU|_cP/*`V'*/, _cU|_cP/*`W'*/, 
    _cU|_cP/*`X'*/, _cU|_cP/*`Y'*/, _cU|_cP/*`Z'*/, _cP/*`['*/, 
    _cP/*`\'*/, _cP/*`]'*/, _cP/*`^'*/, _cP/*`_'*/, 
    _cP/*``'*/, _cL|_cX|_cP/*`a'*/, _cL|_cX|_cP/*`b'*/, _cL|_cX|_cP/*`c'*/, 
    _cL|_cX|_cP/*`d'*/, _cL|_cX|_cP/*`e'*/,
    _cL|_cX|_cP/*`f'*/, _cL|_cP/*`g'*/, 
    _cL|_cP/*`h'*/, _cL|_cP/*`i'*/, _cL|_cP/*`j'*/, _cL|_cP/*`k'*/, 
    _cL|_cP/*`l'*/, _cL|_cP/*`m'*/, _cL|_cP/*`n'*/, _cL|_cP/*`o'*/, 
    _cL|_cP/*`p'*/, _cL|_cP/*`q'*/, _cL|_cP/*`r'*/, _cL|_cP/*`s'*/, 
    _cL|_cP/*`t'*/, _cL|_cP/*`u'*/, _cL|_cP/*`v'*/, _cL|_cP/*`w'*/, 
    _cL|_cP/*`x'*/, _cL|_cP/*`y'*/, _cL|_cP/*`z'*/, _cP/*`{'*/, 
    _cP/*`|'*/, _cP/*`}'*/, _cP/*`~'*/, 0/*`'*/, 
    0/*`\x80'*/, 0/*`\x81'*/, 0/*`\x82'*/, 0/*`\x83'*/, 
    0/*`\x84'*/, 0/*`\x85'*/, 0/*`\x86'*/, 0/*`\x87'*/, 
    0/*`\x88'*/, 0/*`\x89'*/, 0/*`\x8A'*/, 0/*`\x8B'*/, 
    0/*`\x8C'*/, 0/*`\x8D'*/, 0/*`\x8E'*/, 0/*`\x8F'*/, 
    0/*`\x90'*/, 0/*`\x91'*/, 0/*`\x92'*/, 0/*`\x93'*/, 
    0/*`\x94'*/, 0/*`\x95'*/, 0/*`\x96'*/, 0/*`\x97'*/, 
    0/*`\x98'*/, 0/*`\x99'*/, 0/*`\x9A'*/, 0/*`\x9B'*/, 
    0/*`\x9C'*/, 0/*`\x9D'*/, 0/*`\x9E'*/, 0/*`\x9F'*/, 
    0/*`\xA0'*/, 0/*`\xA1'*/, 0/*`\xA2'*/, 0/*`\xA3'*/, 
    0/*`\xA4'*/, 0/*`\xA5'*/, 0/*`\xA6'*/, 0/*`\xA7'*/, 
    0/*`\xA8'*/, 0/*`\xA9'*/, 0/*`\xAA'*/, 0/*`\xAB'*/, 
    0/*`\xAC'*/, 0/*`\xAD'*/, 0/*`\xAE'*/, 0/*`\xAF'*/, 
    0/*`\xB0'*/, 0/*`\xB1'*/, 0/*`\xB2'*/, 0/*`\xB3'*/, 
    0/*`\xB4'*/, 0/*`\xB5'*/, 0/*`\xB6'*/, 0/*`\xB7'*/, 
    0/*`\xB8'*/, 0/*`\xB9'*/, 0/*`\xBA'*/, 0/*`\xBB'*/, 
    0/*`\xBC'*/, 0/*`\xBD'*/, 0/*`\xBE'*/, 0/*`\xBF'*/, 
    0/*`\xC0'*/, 0/*`\xC1'*/, 0/*`\xC2'*/, 0/*`\xC3'*/, 
    0/*`\xC4'*/, 0/*`\xC5'*/, 0/*`\xC6'*/, 0/*`\xC7'*/, 
    0/*`\xC8'*/, 0/*`\xC9'*/, 0/*`\xCA'*/, 0/*`\xCB'*/, 
    0/*`\xCC'*/, 0/*`\xCD'*/, 0/*`\xCE'*/, 0/*`\xCF'*/, 
    0/*`\xD0'*/, 0/*`\xD1'*/, 0/*`\xD2'*/, 0/*`\xD3'*/, 
    0/*`\xD4'*/, 0/*`\xD5'*/, 0/*`\xD6'*/, 0/*`\xD7'*/, 
    0/*`\xD8'*/, 0/*`\xD9'*/, 0/*`\xDA'*/, 0/*`\xDB'*/, 
    0/*`\xDC'*/, 0/*`\xDD'*/, 0/*`\xDE'*/, 0/*`\xDF'*/, 
    0/*`\xE0'*/, 0/*`\xE1'*/, 0/*`\xE2'*/, 0/*`\xE3'*/, 
    0/*`\xE4'*/, 0/*`\xE5'*/, 0/*`\xE6'*/, 0/*`\xE7'*/, 
    0/*`\xE8'*/, 0/*`\xE9'*/, 0/*`\xEA'*/, 0/*`\xEB'*/, 
    0/*`\xEC'*/, 0/*`\xED'*/, 0/*`\xEE'*/, 0/*`\xEF'*/, 
    0/*`\xF0'*/, 0/*`\xF1'*/, 0/*`\xF2'*/, 0/*`\xF3'*/, 
    0/*`\xF4'*/, 0/*`\xF5'*/, 0/*`\xF6'*/, 0/*`\xF7'*/, 
    0/*`\xF8'*/, 0/*`\xF9'*/, 0/*`\xFA'*/, 0/*`\xFB'*/, 
    0/*`\xFC'*/, 0/*`\xFD'*/, 0/*`\xFE'*/, 0/*`\xFF'*/, 
};

            
/////////////////////////////////////////////////////////////
// The memory leak control block
#ifdef _DEBUG

#undef new
#include <malloc.h>

static const unsigned MEM_HEAD =0xA1A3A5A7;
static const unsigned MEM_START=0xC1C3C5C7;
static const unsigned MEM_END  =0xF1F3F5F7;
static const char MEM_NEW='?';
static const char MEM_DEL='!';

struct MemBlk
{
    unsigned m_bHead; //Must = MEM_HEAD
    unsigned m_nSize;
    unsigned m_nLine;
    enum MemType m_bType;
    char m_cFile[0x20];
    MemBlk *m_pNext;
    unsigned m_bStart; //Must = MEM_START
    unsigned char m_Data[sizeof(unsigned)]; //Memory Data, MEM_NEW
};

static LOStream &operator <<(LOStream &os, MemBlk &bl)
{
    os << "[" << (bl.m_bType==M_Array?"Array": bl.m_bType==M_Obj? "Object":
                  bl.m_bType==M_Alloc?"Alloc":"Unknown")
       << " size:" << bl.m_nSize << ", file: " << bl.m_cFile
       << ", line: " << bl.m_nLine << "]";
    return os;
}

inline unsigned WB(unsigned u)
{ return (u+sizeof(unsigned)-1)&~(sizeof(unsigned)-1); }

static MemBlk *m_pFirstBlk=0;
struct MemBlk *MemBlkPtr(void *pmem)
{
    return (struct MemBlk *)(((char *)pmem)-sizeof(MemBlk)+sizeof(unsigned));
}

static void MemValid(struct MemBlk *p, const char *file, int line)
{
    if(p->m_bHead!=MEM_HEAD)
    {
        if((char)p->m_bHead==MEM_NEW)
            FATAL("Not a block (inside another new block?) ("
                  << file<<','<<line<<")");
        if((char)p->m_bHead==MEM_DEL)
            FATAL("Not a block (already deleted?) ("
                  << file<<','<<line<<")");
        FATAL("Not a block (variables on stack?)("
              << file<<','<<line<<")");
    }
    if(p->m_bStart!=MEM_START)
        WARNING("Subscript underflow" << *p<<"("<<
                  file<<','<<line<<")");
    if(*(unsigned *)(p->m_Data+WB(p->m_nSize))!=MEM_END)
        WARNING("Subscript overflow" << *p<<"("<<
                  file<<','<<line<<")");
}

void MemChk(const char *file, unsigned line) //Checking for all memory blocks
{
    MemBlk **q;
    for(q=&m_pFirstBlk; *q!=0; q=&((*q)->m_pNext)) MemValid(*q, file, line);
}
void MemDump(const char *file, unsigned line)
{
    MemBlk **q;
    DUMP("Begin memory block dumping");
    for(q=&m_pFirstBlk; *q!=0; q=&((*q)->m_pNext)) {
        MemValid(*q, file, line);
        *LSys::_pSys << **q << '\n';
    }
    DUMP("End memory block dumping");
}


void *MemAlloc(size_t size, const char *file, unsigned line, enum MemType mt)
{
    MemChk(file, line);
    MemBlk *p=(MemBlk *)malloc(sizeof(struct MemBlk)+WB(size));
    //DUMP("MEMAlloc("<<size<<','<<file<<','<<line<<")="<<p);
    ASSERTMSG(p!=0, "Out of memory: File " << file << ", Line " << line);
    p->m_bType=mt;
    p->m_nSize=size;
    p->m_nLine=line;
    StrNCpy(p->m_cFile, file, sizeof(p->m_cFile)-1);
    p->m_cFile[sizeof(p->m_cFile)-1]='\0';
    p->m_bHead=MEM_HEAD;
    p->m_bStart=MEM_START;
    //p->m_bEnd()=MEM_END;
    *(unsigned *)(p->m_Data+WB(p->m_nSize))=MEM_END;
    StrMemSet(p->m_Data,MEM_NEW,size);
    p->m_pNext=m_pFirstBlk;
    m_pFirstBlk=p;
    return p->m_Data;
}

void MemFree(void *pp, const char *file, unsigned line, enum MemType mt)
{
    MemChk(file, line);
    if(pp==0) return;
    MemBlk *p=MemBlkPtr(pp);
    //DUMP("MEMFree("<<*p<<','<<file<<','<<line<<')');

    MemBlk **q;
    MemValid(p, file, line);
    if(mt!=p->m_bType)
        DUMP("Memory type inconsistant during MemFree("<<
             file<<','<<line<<")"<<*p);
    for(q=&m_pFirstBlk; *q!=0; q=&((*q)->m_pNext))
    {
        if(*q==p)
        {
            *q=p->m_pNext;
            StrMemSet(p, MEM_DEL, sizeof(MemBlk)+WB(p->m_nSize));
            free(p);
            return;
        }
        else MemValid(*q, file, line);
    }
    NOTREACHEDMSG("Block not found during MemFree("<<
                  file<<','<<line<<")" << *p);
}

void *MemRealloc(void *p, size_t s, const char *file, unsigned line,
                 enum MemType mt)
{
    MemChk(file, line);
    if(p==0) return MemAlloc(s,file,line,mt);

    void *pmem;
    MemBlk *pb=MemBlkPtr(p);

    MemValid(pb, file, line);
    
    //DUMP("MEMRealloc("<<*pb<<','<<s<<','<<file<<','<<line<<")");
    if(s > pb->m_nSize)
    {
        pmem=MemAlloc(s, file, line, mt);
        StrMemCpy(pmem, p, pb->m_nSize);
        MemFree(p,file,line,M_Alloc);
        return pmem;
    }
    else if(s < pb->m_nSize) {//else, shrinking        
        pb=(MemBlk *)realloc(pb,sizeof(struct MemBlk)+WB(s));
        ASSERTMSG(pb!=0, "Out of memory: File " << file << ", Line " << line);
        //Need to refill the headers/footers
        pb->m_bType=mt;
        pb->m_nSize=s;
        pb->m_nLine=line;
        StrNCpy(pb->m_cFile, file, sizeof(pb->m_cFile)-1);
        pb->m_cFile[sizeof(pb->m_cFile)-1]='\0';
        //pb->m_bHead=MEM_HEAD;
        //pb->m_bStart=MEM_START; //No need to reset HEAD and START
        //p->m_bEnd()=MEM_END;
        *(unsigned *)(pb->m_Data+WB(pb->m_nSize))=MEM_END;
    }
    return pb->m_Data;
}


static void MemInit()
{
}

static void MemCleanup()
{
    struct MemBlk *p;
    while((p=m_pFirstBlk)!=0)
    {
        DUMP("Block not free" << *p);
        MemFree(p->m_Data,"<cleanup>",0,p->m_bType);
    }
}

int LSys::m_Init=0;

static void iDUMP(const char *s) //Can't use DUMP yet
{
    write(2,LSys::sDebug, strlen(LSys::sDebug));
    write(2,s,strlen(s));
    write(2,"\n",1);
    fsync(2);
}

static const char *print_siginfo(pid_t pid, siginfo_t *siginfo)
{
    static char printbuf[200];
    const char *signame;
    const char *si_code_name;
    static char info[100];

    si_code_name=0;
    info[0]=0;
    
    switch (siginfo->si_code) {
    case SI_USER:
        si_code_name="(SI_USER) kill, sigsend or raise";
        sprintf(info, "uid=%d, pid=%d",
                (int)siginfo->si_uid,(int)siginfo->si_pid);
        break;
//    case SI_KERNEL:  si_code_name="(SI_KERNEL) kernel";
    case SI_QUEUE:
        si_code_name="(SI_QUEUE) sigqueue";
        sprintf(info, "uid=%d, pid=%d",
                (int)siginfo->si_uid, (int)siginfo->si_pid);
        break;
    case SI_TIMER:
        si_code_name="(SI_TIMER) timer expired";
        break;
    case SI_MESGQ:
        si_code_name="(SI_MESGQ) mesq state changed";
        break;
    case SI_ASYNCIO:
        si_code_name="(SI_ASYNCIO) AIO completed";
        break;
//    case SI_SIGIO:   si_code_name="(SI_SIGIO) queued SIGIO";
    }
    switch (siginfo->si_signo) {
    case SIGHUP:   signame="SIGHUP"; break;
    case SIGINT:   signame="SIGINT"; break;
    case SIGQUIT:  signame="SIGQUIT"; break;
    case SIGILL:   signame="SIGILL";
        sprintf(info, "addr=^%X", (unsigned)(siginfo->si_addr));
        switch (siginfo->si_code) {
        case ILL_ILLOPC: si_code_name="(ILL_ILLOPC) illegal opcode";
            break;
        case ILL_ILLOPN: si_code_name="(ILL_ILLOPN) illegal operand";
            break;
        case ILL_ILLADR: si_code_name="(ILL_ILLADR) illegal addressing mode";
            break;
        case ILL_ILLTRP: si_code_name="(ILL_ILLTRP) illegal trap";
            break;
        case ILL_PRVOPC: si_code_name="(ILL_PRVOPC) privileged opcode";
            break;
        case ILL_PRVREG: si_code_name="(ILL_PRVREG) privileged register";
            break;
        case ILL_COPROC: si_code_name="(ILL_COPROC) coprocessor error";
            break;
        case ILL_BADSTK: si_code_name="(ILL_BADSTK) internal stack error";
            break;
        }
        break;
    case SIGTRAP:  signame="SIGTRAP"; break;
    case SIGABRT:  signame="SIGABRT"; break;
    case SIGFPE:   signame="SIGFPE"; 
        sprintf(info, "addr=^%X", (unsigned)(siginfo->si_addr));
        switch (siginfo->si_code) {
        case FPE_INTDIV:
            si_code_name="(FPE_INTDIV) integer divide by zero";
            break;
        case FPE_INTOVF:
            si_code_name="(FPE_INTOVF) integer overflow";
            break;
        case FPE_FLTDIV:
            si_code_name="(FPE_FLTDIV) floating point divide by zero";
            break;
        case FPE_FLTOVF:
            si_code_name="(FPE_FLTOVF) floating point overflow";
            break;
        case FPE_FLTUND:
            si_code_name="(FPE_FLTUND) floating point underflow";
            break;
        case FPE_FLTRES:
            si_code_name="(FPE_FLTRES) floating point inexact result";
            break;
        case FPE_FLTINV:
            si_code_name="(FPE_FLTINV) floating point invalid operation";
            break;
        case FPE_FLTSUB:
            si_code_name="(FPE_FLTSUB) subscript out of range";
            break;
        }
        break;
    case SIGKILL:  signame="SIGKILL"; break;
    case SIGBUS:   signame="SIGBUS"; 
        sprintf(info, "addr=^%X", (unsigned)(siginfo->si_addr));
        switch (siginfo->si_code) {
        case BUS_ADRALN:
            si_code_name="(BUS_ADRALN) invalid address alignment";
            break;
        case BUS_ADRERR:
            si_code_name="(BUS_ADRERR) non-existant physical address";
            break;
        case BUS_OBJERR:
            si_code_name="(BUS_OBJERR) object specific hardware error";
            break;
        }
        break;
    case SIGSEGV:  signame="SIGSEGV"; 
        sprintf(info, "addr=^%X", (unsigned)(siginfo->si_addr));
        switch (siginfo->si_code) {
        case SEGV_MAPERR:
            si_code_name="(SEGV_MAPERR) address not mapped to object";
            break;
        case SEGV_ACCERR:
            si_code_name="(SEGV_ACCERR) invalid permissions for mapped object";
            break;
        }
        break;
    case SIGSYS:   signame="SIGSYS"; break;
    case SIGPIPE:  signame="SIGPIPE"; break;
    case SIGALRM:  signame="SIGALRM"; break;
    case SIGTERM:  signame="SIGTERM"; break;
    case SIGCHLD:  signame="SIGCHLD"; 
        sprintf(info, "pid=%d, status=%x", (int)siginfo->si_pid,
                (int)siginfo->si_status);
        switch (siginfo->si_code) {
        case CLD_EXITED:
            //Normal, get back without doing anything
            return 0;
            //si_code_name="(CLD_EXITED) child has exited";
            //break;
        case CLD_KILLED:
            si_code_name="(CLD_KILLED) child was killed";
            break;
        case CLD_DUMPED:
            si_code_name="(CLD_DUMPED) child terminated abnormally";
            break;
        case CLD_TRAPPED:
            si_code_name="(CLD_TRAPPED) traced child has trapped";
            break;
        case CLD_STOPPED:
            si_code_name="(CLD_STOPPED) child has stopped";
            break;
        case CLD_CONTINUED:
            si_code_name="(CLD_CONTINUED) stopped child has continued";
            break;
        }
        break;
    case SIGUSR1:  signame="SIGUSR1"; break;
    case SIGUSR2:  signame="SIGUSR2"; break;
    default:
        signame="SIG???";
    }
    if(si_code_name==0)
    {
        static char si_code_buf[20];
        sprintf(si_code_buf, "(%d) Unknown si_code", siginfo->si_code);
        si_code_name=si_code_buf;
    }
    if(info[0])
        sprintf(printbuf, HIR"[%d]%s%s[%s](%s)"NOR"\n", (int)pid,signame,
                si_code_name, info, strerror(siginfo->si_errno));
    else sprintf(printbuf, HIR"[%d]%s%s(%s)"NOR"\n", (int)pid,signame,
                 si_code_name, strerror(siginfo->si_errno));
    return printbuf;
}

void LSys::sig_handler(int, siginfo_t *sig, void *)
{
    const char *str=print_siginfo(getpid(), sig);
    if(str) {
        write(2, str, strlen(str));
        fsync(2);
    }
    
    switch(sig->si_signo) {
        //Terminate
    case SIGSEGV:
    case SIGBUS:
    case SIGILL:
        _exit(1);
        //Ignore
    case SIGCHLD:
    default:;
    }
}

LSys::LSys()
{
    if(!m_Init)
    {
        struct sigaction sa;
        sa.sa_flags=SA_SIGINFO;
        sa.sa_sigaction=sig_handler;
        //Setup Signal handler
        //No SIGCHLD because it is a NORMAL signal, while we only interested in
        //the disastrous signals
        //if(sigaction(SIGCHLD, &sa, 0)!=0) iDUMP("Error in sigaction(SIGCHLD)");
        if(sigaction(SIGFPE, &sa, 0)!=0) iDUMP("Error in sigaction(SIGFPE)");
        if(sigaction(SIGBUS, &sa, 0)!=0) iDUMP("Error in sigaction(SIGBUS)");
        if(sigaction(SIGSEGV, &sa, 0)!=0) iDUMP("Error in sigaction(SIGSEGV)");
//        unsigned u=1;
//        if(sizeof(unsigned)!=4)
//            iDUMP("Warning:Word size should be 32bit.");
//#ifdef __LITTLE_ENDIAN
//        if((*(char *)&u)!=1)
//            iDUMP("Warning:Byte-order should be big-endian(file:"__FILE__")");
//#else //__BIG_ENDIAN
//        if((*(char *)&u)!=0)
//            iDUMP("Warning:Byte-order should be little-endian(file:"__FILE__
//                  ")");
//#endif
        MemInit();
        iDUMP("Program starting...");
    }
    m_Init++;
}
    
#include <sys/time.h>
#include <sys/resource.h>
LSys::~LSys()
{
    if(!--m_Init)
    {
        struct rusage ru;
        MemCleanup();
        
        getrusage(RUSAGE_SELF,&ru);
        DUMP("CPU time spent(in seconds):"<<
             ru.ru_utime.tv_sec+1e-6*ru.ru_utime.tv_usec);
    }
}

void *operator new[](size_t size) throw(bad_alloc)
{
    return MemAlloc(size, "<new[]>", 0, M_Array);
}
void *operator new[](size_t size, const char *file, unsigned line)
    throw(bad_alloc)
{
    return MemAlloc(size, file, line, M_Array);
}
//void *operator new[](size_t size, void * ptr)throw(bad_alloc)
//{
//    return (new(ptr, size, "<Unknown>", 0, true) LMemoryBlock())->Data();
//}
//void *operator new[](size_t size, void * ptr,
//                     const char *file, unsigned line)
//    throw(bad_alloc)
//{
//    return (new(ptr, size, file, line, true) LMemoryBlock())->Data();
//}

void operator delete[](void *pmem) throw()
{
    if(pmem!=0) MemFree(pmem,"<delete[]>",0, M_Array);
}

void *operator new(size_t size)throw(bad_alloc)
{
    return MemAlloc(size, "<new>", 0, M_Obj);
}
void *operator new(size_t size, const char *file, unsigned line)
    throw(bad_alloc)
{
    return MemAlloc(size, file, line, M_Obj);
}
//void *operator new(size_t size, void * ptr)throw(bad_alloc)
//{
//    return (new(ptr, size, "<Unknown>", 0) LMemoryBlock())->Data();
//}
//void *operator new(size_t size, void * ptr, const char *file, unsigned line)
//    throw(bad_alloc)
//{
//    return (new(ptr, size, file, line) LMemoryBlock())->Data();
//}
void operator delete(void *pmem)throw()
{
    if(pmem!=0) MemFree(pmem,"<delete>",0,M_Obj);
}

#endif //_DEBUG

#ifdef GENERAL_TEST
#define new new(__FILE__,__LINE__)

int main()
{
    {
        _IO << "=*= Message system\n";
        INFO("This is a test message");
        WARNING("This is a test warning");
        ERROR("This is a test error");
        //FATAL("This is a test fatal");
        DUMP("This is a test debug dump");
    }

    {
        char str[]="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        _IO << "=*= CRC checking\n";
        _IO << "CRC of " << str << " is "
            << CRC((unsigned char *)str, sizeof(str)) << '\n';
    }
    {
        double data[10]={12.5,2,13,1,331,4,31,21,35,12};
        int index[10];
        _IO << "=*=QuickSort test\n";
        _IO << "[ ";
        for(int i=0;i<10;i++)
            _IO << data[i] << ' ';
        _IO << "];\n";
        IndexQuickSort(data,index, 10);
        _IO << "[ ";
        for(int i=0;i<10;i++)
            _IO << index[i] << ' ';
        _IO << "];\n";
    }
    {
        char strbuf[200];
        double sec= 231342;
        _IO << "=*=LTime demo\n";
        _IO << "PrettyPrint of "<<sec<<" time seconds:\n";
        LTime::PrettyPrint(strbuf, sizeof(strbuf), sec);
        _IO << strbuf << '\n';
        sec=3892740.233;
        _IO << "PrettyPrint of "<<sec<<" time seconds:\n";
        LTime::PrettyPrint(strbuf, sizeof(strbuf), sec);
        _IO << strbuf << '\n';
    }
        
//    int *a=new int[322];
//    delete a;
//    delete []a;
    return 0;
}

#endif //_TEST

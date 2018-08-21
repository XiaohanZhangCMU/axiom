//Last Modified : Wed Aug  1 18:38:54 2007

#ifndef _FILECLS_H
#define _FILECLS_H

#include <stdio.h>
#include <stdlib.h>
#include "general.h"

#ifndef MPLIB_SINGLE
class LFile : public LFDStream
#else
class LFile : public LStream
#endif
{
public:
    //Type members
    enum {//Seek parameters
        SeekBegin=0,
        SeekRel=1,
        SeekEnd=2
    };
    enum { //Open modes
        O_Read=O_RDONLY,
        O_Write=O_WRONLY,
        O_AccMode=O_ACCMODE,
        O_ReadWrite=O_RDWR,
        O_Create=O_CREAT,
        O_Exclusive=O_EXCL,
        O_Truncate=O_TRUNC,
        O_Append=O_APPEND,
        
        O_WriteNew=O_WRONLY|O_TRUNC|O_CREAT
    };
    enum { //Supported file formats
        FILE_UNKNOWN,
        FILE_GZ,
        FILE_BZ2
    };
protected:
    //Data members
    struct stat m_stat;
    //Function members
    LFile(int h) :handle(h)
    {
        if(fstat(handle, &m_stat) < 0) {
            SYSERROR("Can't stat file !");
            LSys::Abort();
        }
    }
public:
    int handle; 
    LFile(const char *name, unsigned attr, mode_t mode=00644)
    {
        handle=open(name, attr, mode); // c open function
        if(handle < 0) {
            SYSERROR("File open error : "<< name );
            LSys::Abort();
        }
        Stat();
    }
    virtual ~LFile()
    {
        close(handle);
    }
    virtual int Read(char *buffer, int size)
    {
        ASSERT(buffer!=NULL || size==0);
        ASSERT(handle > 0);
        return read(handle, buffer, size);
    }
    virtual int Write(const char *data, int size)
    {
        ASSERT(data!=NULL || size==0);
        ASSERT(handle > 0);
        return write(handle, data, size);
    }
    
    virtual int Seek(int gap, int mode)
    { return lseek(handle, gap, mode); }

    int SeekFromBegin(int gap)
    { return Seek(gap, SeekBegin); }

    int SeekFromEnd(int gap)
    { return Seek(gap, SeekEnd); }

    int SeekRelatively(int gap)
    { return Seek(gap, SeekRel); }
    void Stat()
    {
        if(fstat(handle, &m_stat) < 0) {
            SYSERROR("LFile::Stat(): Can't stat file");
            LSys::Abort();
        }
    }
    //Ignore types of (long long), assume file size less than unsigned long
    unsigned long Length() { return (unsigned long)m_stat.st_size; }

    static LFile *MakeTemp(char *tpl)
    {
        int handle;
        handle=mkstemp(tpl);
        if(handle > 0) return new LFile(handle);
        return NULL;
    }
    static LFile *Open(const char *name, unsigned attr, mode_t mode=00644)
    {
        int handle;
        handle=open(name, attr, mode);
        if(handle > 0) return new LFile(handle);
        return NULL;
    }
    static void Close(LFile *&f)
    {
        if(f!=NULL) {
            ASSERT(f->handle>0);
            delete f;
            f=NULL;
        }
    }
    static int Unlink(const char * name)
    {
        int ret=unlink(name);
        return ret==-1?-errno:ret;
    }
    static int Rename(const char *oldname, const char *newname)
    {
        int ret=rename(oldname, newname);
        return ret==-1?-errno:ret;
    }
    static int Symlink(const char* oldpath, const char *newpath)
    {
        int ret=symlink(oldpath, newpath);
        return ret==-1?-errno:ret;
    }
    static int Stat(const char *name, struct stat *buf=0)
    {
        static struct stat buff;
        if(buf==0) buf=&buff; //Means user does not want the stat_buf
        int ret=stat(name, buf);
        return ret==-1?-errno:ret;
    }
    //LoadFileToString(fname, buffer, size): Load whole file into a string
    //  up to `size' byes. If size==0, use malloc to get the buffer, return
    //  number of bytes loaded.
    static int LoadToString(const char *fname, char *&buffer, int size);
    static void BZipFile(const char *fname, bool bg=false);
    static void GZipFile(const char *fname, bool bg=false);

    static void SubHomeDir(const char *fname, char *extname);
    
    int FileType();
};

#endif //_FILECLS_H

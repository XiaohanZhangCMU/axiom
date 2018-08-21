// Last Modified : Wed Apr 14 11:42:51 2010

#include "filecls.h"

int LFile::FileType()
{
    char buf[3];
    if(Read(buf, 3)!=3) return false;
    
    if(buf[0]=='B' && buf[1]=='Z' && buf[2]=='h') return FILE_BZ2;
    if(buf[0]=='\037' && buf[1]=='\213') return FILE_GZ;
    return FILE_UNKNOWN;
}
    
static LFile *Unzip(const char *fname, const char *exec)
{
    char cmd[1000];
    char tmp[]="/tmp/.TMPXXXXXX";
    LFile *fp;
    
    fp=LFile::MakeTemp(tmp);
    if(fp==NULL) ERROR("Can't create tempfile");
    StrPrintf(cmd, "%s -cd \"%s\" > \"%s\"", exec, fname, tmp);
    system(cmd);
    fp->Stat();
    LFile::Unlink(tmp); 
    return fp;
}

void LFile::BZipFile(const char *fname, bool bg)
{
    char buffer[1000]="bzip2 -9 -f \"";

    StrCat(buffer, fname);
    StrCat(buffer, "\"");
    if(bg) StrCat(buffer, " &");
    //INFO(buffer); 
    system(buffer);
}

void LFile::GZipFile(const char *fname, bool bg)
{
    char buffer[1000]="gzip -9 -f  \"";

    StrCat(buffer, fname);
    StrCat(buffer, "\"");
    if(bg) StrCat(buffer, " &");
    /* This INFO causes segmentation fault by writing buffer to the logfile,
       which was already closed in Organizer::closefiles(int, bool) in organizer.cpp. */
    //INFO(buffer);
    system(buffer);
}

int LFile::LoadToString(const char *fname, char *&buffer, int size)
{
    int s, ret;
    LFile *file, *f2;
    file=Open(fname, LFile::O_Read); // c open function
    if(file==NULL) {
        SYSERROR("LoadConfig: Can not open "YEL<<fname<<NOR);
        return 0;
    }
    switch(file->FileType())
    {
    case FILE_BZ2:
        f2=Unzip(fname, "bzip2");
        goto l_unzipped;
    case FILE_GZ:
        f2=Unzip(fname, "gzip");
    l_unzipped:
        if(f2==NULL) return 0;
        Close(file);
        file=f2;
        //fall through (for a quick break)
    case FILE_UNKNOWN:
        break;
    default: NOTREACHED();
    }
    //Now file is the pointer to the true file.
    if(size==0) {
        s=(int)file->Length();
        buffer=(char *)MAlloc(s+1);
        if(buffer==NULL) {
            ERROR("LoadFileToString: Out of memory.");
            return 0;
        }
    }
    else s=size;
    file->SeekFromBegin(0);
    ret=file->Read(buffer, s);
    Close(file);
    if(ret!=s) {
        ERROR("LoadFileToString: File size changed!");
        s=ret;
    }
    buffer[s]=0;//Terminate
    return s;
}

void LFile::SubHomeDir(const char *fname, char *extname)
{
    char Home[400], *ph, *pt;
    char fname_local[400];
//    INFO("getenv");
    ph=getenv("HOME"); 
    strcpy(fname_local, fname);
    pt=strchr(fname_local, '~');
//    INFO("ph = "<<ph);
    if((pt==NULL)||(ph==NULL))
    {
       if(fname_local!=extname) strcpy(extname,fname_local);
//       INFO_Printf("~ not found in %s\n",fname);
       return;
    }
    *pt=0;
    strcpy(Home,fname_local);
    strcpy(Home+strlen(Home),ph);       
    strcpy(Home+strlen(Home),pt+1);
    strcpy(extname,Home);
    if(fname_local!=extname)  *pt='~';
    INFO_Printf("replace %s by %s\n",fname_local,extname);
}

#ifdef _ZIPCAT_TEST
int main(int argc, char *argv[]) //ZipCat
{
    char *buffer;
    int nb;
    if(argc==1) _Error<< "Usage: zipcat [file] ...\n"
                    "dump contents of bzip2 or gzip compressed file.\n";    
    for(int i=1;i<argc;i++) {
        nb=LFile::LoadToString(argv[i], buffer, 0);
        _IO << buffer;
        Free(buffer);
    }
}

#endif //_ZIPCAT_TEST

//const int LEN=30000;
//int main(int argc, char *argv[])
//{
//    char buf[LEN+1];
//    char buf2[LEN+1];
//    int n;
//    if(argc != 2)
//    {
//        _Error <<
//            "Usage : dup <filename>\n"
//            "  Converts '\\r' to '\\n'\n";
//    }
//    else
//    {
//        LFile file(argv[1]);
//        while((n=file.GetTextString(buf, LEN)))
//        {
//            if((buf[1]=='V')&&(buf[2]=='2')&&(buf[3]=='3')&&(buf[4]=='4'))
//            {
//                file.GetTextString(buf2,LEN);
//                while((buf2[1]=='V')&&(buf2[2]=='2')&&
//                      (buf2[3]=='3')&&(buf2[4]=='4'))
//                {
//                    StrCpy(buf,buf2);
//                    file.GetTextString(buf2,LEN);
//                }
//                _IO <<buf+9<<' ';
//                _IO << buf2 +9<<' ';
//                file.GetTextString(buf,LEN);
//                _IO << buf +23<<' ';
//                file.GetTextString(buf,LEN);
//                _IO << buf+20<<' ';
//                file.GetTextString(buf,LEN);
//                buf[27]=0;
//                _IO<<buf +8<<buf +35<<'\n';
//            }
//            buf[n]='\0';
//            _IO << buf;
//        }
//        ASSERT(n==0);
//    }
//    return 0;
//}

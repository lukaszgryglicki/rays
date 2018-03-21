#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>

#define MAXSIZ 512

typedef struct _BMPTag
{
 char ident[2];
 int fsize;
 int dummy;
 int offset;
 int dummy2;
 int bm_x;
 int bm_y;
 short planes;
 short bpp;
 int compress;
 int nbytes;
 int no_matter[4];
} BMPTag;

void error(char* fmt, ...)
{
 va_list lst;
 va_start(lst,fmt);
 printf("Error: \t");
 vprintf(fmt,lst);
 printf("\n\tClient halted\n");
 fflush(stdout);
 va_end(lst);
 exit(1);
}

void get_host(char* to, char* name)
{
 struct hostent* hose;
 struct in_addr addr;
 printf("DNS resolving IP for: %s\n", name);
 hose = gethostbyname(name);
 if (!hose)
   {
    perror("gethostbyname");
    error("Cannot get IP for hostname: %s\n", name);
   }
 memcpy((void*)&addr,(void*)(hose->h_addr_list[0]),sizeof(struct in_addr));
 strcpy(to, inet_ntoa(addr));
 printf("Got IP %s for hostname %s\n", to, name);
 fflush(stdout);
}

void init_bmp(BMPTag* b)
{
 int i;
 b->ident[0]='B';
 b->ident[1]='M';
 b->fsize=0;
 b->dummy=0;
 b->offset=sizeof(BMPTag);
 b->bm_x=b->bm_y=0x20;
 b->dummy2=40;
 b->bpp=0x18;
 b->planes=1;
 b->compress=0;
 b->nbytes=3*32*32;
 for (i=0;i<4;i++) b->no_matter[i]=0;
}

void write_bmp(int s, char* fn)
{
 FILE* plik;
 BMPTag bm_handle;
 int x,y,n,i,j;
 char buff[MAXSIZ+1];
 char* lbuff;
 if (!(plik = fopen(fn, "w"))) error("cannot write to file: %s", fn);
 if ((n=read(s,buff,2*sizeof(int)))==-1) error("read from socket");
 if (n != 2*sizeof(int)) error("bad server answer, read %d bytes", n);
 memcpy(&x, buff, sizeof(int));
 memcpy(&y, (void*)(buff+sizeof(int)), sizeof(int));
 if (x < 0 || x > 30000) error("bad x value: %d", x);
 if (y < 0 || y > 20000) error("bad y value: %d", y);
 init_bmp(&bm_handle);
 fprintf(plik,"%c%c",'B', 'M');
 bm_handle.bm_y = y;
 bm_handle.bm_x = x;
 printf("Dimnesions: (%dx%d)\n", bm_handle.bm_x, bm_handle.bm_y);
 bm_handle.fsize = sizeof(BMPTag)+(bm_handle.bm_y*bm_handle.bm_x*3);
 fwrite(&bm_handle.fsize,4,1,plik);
 fwrite(&bm_handle.dummy,4,1,plik);
 bm_handle.offset=sizeof(BMPTag);
 bm_handle.planes=1;
 bm_handle.bpp=24;
/* bm_handle.nbytes = s->x * s->y * 3;*/
 fwrite(&bm_handle.offset,4,1,plik);
 fwrite(&bm_handle.dummy2,4,1,plik);
 fwrite(&bm_handle.bm_x,4,1,plik);
 fwrite(&bm_handle.bm_y,4,1,plik);
 fwrite(&bm_handle.planes,2,1,plik);
 fwrite(&bm_handle.bpp,2,1,plik);
 fwrite(&bm_handle.compress,4,1,plik);
 fwrite(&bm_handle.nbytes,4,1,plik);
 for (i=0;i<4;i++)  fwrite(&bm_handle.no_matter[i],4,1,plik);
 fseek(plik,bm_handle.offset,SEEK_SET);
 lbuff = (char*)malloc(bm_handle.bm_x*3+1);
  for (i=0;i<bm_handle.bm_y;i++)  
    {
     j = 0;
     while (j < bm_handle.bm_x*3)
      {
       n = read(s, lbuff+j, bm_handle.bm_x*3-j);
       j += n;
/*       printf("partial: %d\n", j);*/
      }
/*     printf("line %d finihed\n", i);*/
     for (j=0;j<bm_handle.bm_x*3;j++) fprintf(plik,"%c", lbuff[j]);
    }
 fclose(plik);
 printf("Writen bitmap: %s\n", fn);
}

void client(char* ipstr, int port, char* out)
{
 int sock;
 struct sockaddr_in server;
 char buff[MAXSIZ+1];
 printf("Starting up client out=%s\n", out);
 if (port<1 || port>0xffff) error("Port: %d is invalid", port);
 server.sin_family = AF_INET;
 server.sin_port = htons(port);
 server.sin_addr.s_addr = inet_addr(ipstr);
 if (server.sin_addr.s_addr==(in_addr_t)(-1)) error("IP %s is invalid.\n", ipstr);
 printf("Connecting to %s, at port %d...\n", ipstr, port);
 sock = socket(AF_INET, SOCK_STREAM, 0);
 if (sock==-1)
      { perror("socket"); error("create socket failed."); }
 if (connect(sock,(struct sockaddr*)&server, sizeof(struct sockaddr_in))==-1)
      { perror("connect"); error("connect to %s:%d failed",ipstr,port); }
 printf("connected.\n");
 fflush(stdout);
 sprintf(buff, "get\n");
 if (write(sock,buff,strlen(buff))==-1) error("write to socket");
 write_bmp(sock, out);
 close(sock);
}

void help()
{
 printf("options are: -p PORT, -i IPNUM, -s HOSTNAME -o BMPNAME\n");
}


int main(int lb, char** par)
{
 char u;
 char ip[16];
 char port[12];
 char out[256];
 strcpy(ip,"127.0.0.1");
 strcpy(port,"2500");
 strcpy(out,"output.bmp");
 printf("Starting client.., use -h to see available options...\n");
 while ((u = getopt(lb,par,"i:p:s:o:h"))!=-1)
 {
  switch (u)
   {
    case 'i': if (strlen(optarg)<16)  strcpy(ip, optarg);      break;
    case 'o': if (strlen(optarg)<256) strcpy(out, optarg);     break;
    case 'p': if (strlen(optarg)<12)  strcpy(port, optarg);    break;
    case 's': get_host(ip, optarg);                            break;
    case 'h': help(); return 0;
    default: printf("Unrecognized option\n"); return 1;
   }
 }
 client(ip,atoi(port), out);
 return 0;
}


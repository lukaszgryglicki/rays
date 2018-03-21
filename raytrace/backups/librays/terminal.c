#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <sys/socket.h>

#define MSGSIZ 2048


int client(char* ipstr, int port, char* cmd)
{
 int sock;
 char buff[MSGSIZ+1];
 struct sockaddr_in server;
 if (port<1 || port>0xffff) return -1;
 server.sin_family = AF_INET;
 server.sin_port = (in_port_t)(htons(port));
 server.sin_addr.s_addr = inet_addr(ipstr);
 if (server.sin_addr.s_addr==(in_addr_t)(-1)) return -1;
 printf("Connecting to %s, at port %d\n", ipstr, port);
 sock = socket(AF_INET, SOCK_STREAM, 0);
 if (sock==-1) return -1;
 if (connect(sock,(struct sockaddr*)&server, sizeof(struct sockaddr_in))==-1) return -1;
 printf("Connect OK\n");
 strcpy(buff,"");
 if (write(sock, cmd, strlen(cmd)+1)==-1) 
    { 
     printf("write error!\n"); 
     return -1; 
    }
 strcpy(buff,"");
 if (read(sock, buff, MSGSIZ)==-1) 
    { 
     printf("read error!\n"); 
     return -1; 
    }
 printf("SERVER: %s\n", buff);
 close(sock);
 printf("Socket closed.\n");
 return 0;
}


int main(int lb, char** par)
{
 char u;
 char ip[16];
 char cmd[MSGSIZ+1];
 char port[12];
 printf("Options are: -p PORTNUM, -i IPNUM -c 'CMD'\n");
 strcpy(ip,"127.0.0.1");
 strcpy(port,"2500");
 strcpy(cmd, "rem");
 while ((u = getopt(lb,par,"i:p:c:"))!=-1)
 {
  switch (u)
   {
    case 'i': if (strlen(optarg)<16) strcpy(ip, optarg);      break;
    case 'c': if (strlen(optarg)<MSGSIZ) strcpy(cmd, optarg); break;
    case 'p': if (strlen(optarg)<12) strcpy(port, optarg);    break;
    default: printf("Unrecognized option\n"); return 1;
   }
 }
 if (client(ip,atoi(port), cmd)==-1) printf("Error.\n");
 return 0;
}

/* end, written by MorgothDBMA: morgothdbma@o2.pl, tel +48693582014 */
/* license: BSD */

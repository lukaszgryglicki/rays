#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void get_char(char *s)
{
 int i;
 char str[32];
 i = fscanf(stdin, "%s", str);
 if (i != 1) { printf("Cannot scanf CHAR\n"); exit(1); }
 *s = str[0];
 fflush(stdin);
 /*printf("read CHAR \"%c(%d)\"\n", *s, *s);*/
}

void get_int(int *s, int def)
{
 int i;
 char str[32];
 *s = 0;
 i = fscanf(stdin, "%s", str);
 if (i != 1) { printf("Cannot scanf CHAR\n"); exit(1); }
 if (str[0] < '0' || str[0] > '9') *s = def;
 else *s = atoi(str);
 fflush(stdin);
 /*printf("read CHAR \"%c(%d)\"\n", *s, *s);*/
}

void get_double(double *s, double def)
{
 int i;
 char str[32];
 *s = 0.0;
 i = fscanf(stdin, "%s", str);
 if (i != 1) { printf("Cannot scanf CHAR\n"); exit(1); }
 if ((str[0] < '0' || str[0] > '9') && str[0] != '.') *s = def;
 else *s = atof(str);
 fflush(stdin);
 /*printf("read CHAR \"%c(%d)\"\n", *s, *s);*/
}

void get_string(char* s)
{
 int i;
 i = fscanf(stdin, "%s", s);
 if (i != 1) { printf("Cannot scanf STRING\n"); exit(1); }
 fflush(stdin);
 /*printf("read STRING \"%s\"\n", s);*/
}

void wrapper()
{
 char s,ss;
 int nr;
 int os_unix;
 char dircmd[32];
 char diropt[16];
 char rayscmd[2048];
 char temp[2048];
 char wdir[1024];
 double dd;
 time_t tm;
 s = nr = 0;
 getcwd(wdir, 512);
 if (wdir[strlen(wdir)-1] != '/') strcat(wdir, "/");
 printf("WorkDir = %s\n", wdir);
 printf("RAYS wrapper\n");
 printf("Usage: (Option1-shortcut1,..,OptionN-shortcut-n) [default]\n");
 printf("use shortcut1-n or 'd' for default\n");
 printf("Operating system (UNIX-u,WINDOWS-w) [UNIX]: ");
 get_char(&s);
 if (s == 'w') { os_unix = 0; printf("OS set to: WINDOWS\n"); }
 else { os_unix = 1; printf("OS set to UNIX\n"); }
 strcpy(rayscmd, wdir);
 if (os_unix)
   {
    printf("DirCommand (\"ls\"-l,other-o) [ls]: ");
    get_char(&s);
    if (s == 'o')
      {
       printf("Enter new dir command for UNIX OS: ");
       get_string(dircmd);
       strcat(dircmd, " ");
       printf("Command options? (\"none\"-n,set-s) [none]: ");
       get_char(&s);
       if (s == 's')
         { 
          printf("Dir command (%s) option: ", dircmd);
	  get_string(diropt);
          strcat(dircmd, diropt);
	  strcat(dircmd, " ");
         }
      }
    else strcpy(dircmd,"ls ");
    printf("Dir command is: %s\n", dircmd);
    printf("RAYS command: (\"rays.debug\"-g,\"rays.fast\"-f,other-o) [rays.fast]: ");
    get_char(&s);
    if (s == 'g') strcpy(temp, "rays.debug ");
    else if (s == 'o')
      {
       printf("Enter new RAYS cmd for UNIX OS: ");
       get_string(temp);
       strcat(temp, " ");
      }
    else strcpy(temp, "rays.fast ");
    strcat(rayscmd, temp);
    printf("RAYScmd: %s\n", rayscmd);
   }
 else
   {
    printf("DirCommand (\"dir\"-d,other-o) [dir]: ");
    get_char(&s);
    if (s == 'o')
      {
       printf("Enter new dir command for WINDOWS OS: ");
       get_string(dircmd);
       strcat(dircmd, " ");
       printf("Command options? (\"none\"-n,set-s) [none]: ");
       get_char(&s);
       if (s == 's')
         { 
          printf("Dir command (%s) option: ", dircmd);
	  get_string(diropt);
          strcat(dircmd, diropt);
	  strcat(dircmd, " ");
         }
      }
    else strcpy(dircmd,"dir ");
    printf("Dir command is: %s\n", dircmd);
    printf("RAYS command: (\"cyg_rays.exe\"-c,\"mingw_rays.exe\"-m,other-o) [cyg_rays.exe]: ");
    get_char(&s);
    if (s == 'm') strcpy(temp, "mingw_rays.exe ");
    else if (s == 'o')
      {
       printf("Enter new RAYS cmd for WINDOWS OS: ");
       get_string(temp);
       strcat(temp, " ");
      }
    else strcpy(temp, "cyg_rays.exe ");
    strcat(rayscmd, temp);
    printf("RAYScmd: %s\n", rayscmd);
   }
 /* we have dir command in dircmd and rays command in rayscmd */
 printf("Want help? (\"no\"-d,\"yes\"-y) [\"no\"]: ");
 get_char(&s);
 if (s == 'y')
   {
    sprintf(temp, "%s -H", rayscmd);
    system(temp);	/* help */
   }
 printf("Input file selection:\n");
 printf("Scenefile to use (\"dat/table.dat\"-d,other-o) [dat/table.dat]: ");
 get_char(&s);
 if (s == 'd')
   {
    strcat(rayscmd, "-i ");
    strcat(rayscmd, wdir);
    strcat(rayscmd, "dat/table.dat ");
   }
 else
   {
    sprintf(temp, "%s %sdat/", dircmd, wdir);
    printf("Executing system command: %s\n", temp);
    system(temp);
    printf("Enter secenefile name: without path and dat/ directory: ");
    get_string(temp);
    strcat(rayscmd, "-i ");
    strcat(rayscmd, wdir);
    strcat(rayscmd, "dat/");
    strcat(rayscmd, temp);
    strcat(rayscmd, " ");
   }
 printf("Texture directory prefix (\"use current dir\"-d,other-o) [use current]: ");
 get_char(&s);
 if (s == 'o')
   {
    printf("Enter texture prefix dir (current is \"%s\"): ", wdir);
    get_string(temp);
    if (temp[strlen(temp)-1] != '/') strcat(temp, "/");
    strcat(rayscmd, "-T ");
    strcat(rayscmd, temp);
    strcat(rayscmd, " ");
   }
 else
   {
    strcat(rayscmd, "-T ");
    strcat(rayscmd, wdir);
    strcat(rayscmd, " ");
   }
 printf("Output file selection: (\"bmp/screen.bmp\"-d,other-o) [bmp/screen.bmp]: ");
 get_char(&s);
 if (s == 'd')
   {
    strcat(rayscmd, "-o ");
    strcat(rayscmd, wdir);
    strcat(rayscmd, "bmp/screen.bmp ");
   }
 else
   {
    printf("Enter output name: without path and bmp/ directory\n");
    get_string(temp);
    strcat(rayscmd, "-o ");
    strcat(rayscmd, wdir);
    strcat(rayscmd, "bmp/");
    strcat(rayscmd, temp);
    strcat(rayscmd, " ");
   }
 printf("Disable signal handlers? (\"no\"-d,\"yes\"-y) [\"no\"]: ");
 get_char(&ss);
 if (ss == 'y') 
   {
    strcat(rayscmd, "-L ");
   }
 else
   {
    printf("Demand file selection: (\"bmp/demand.bmp\"-d,other-o) [bmp/demand.bmp]: ");
    get_char(&s);
    if (s == 'd')
      {
       strcat(rayscmd, "-D ");
       strcat(rayscmd, wdir);
       strcat(rayscmd, "bmp/demand.bmp ");
      }
    else
      {
       printf("Enter demand name: without path and bmp/ directory\n");
       get_string(temp);
       strcat(rayscmd, "-D ");
       strcat(rayscmd, wdir);
       strcat(rayscmd, "bmp/");
       strcat(rayscmd, temp);
       strcat(rayscmd, " ");
      }
    printf("Signalled file selection: (\"bmp/signal.bmp\"-d,other-o) [bmp/signal.bmp]: ");
    get_char(&s);
    if (s == 'd')
      {
       strcat(rayscmd, "-S ");
       strcat(rayscmd, wdir);
       strcat(rayscmd, "bmp/signal.bmp ");
      }
    else
      {
       printf("Enter signal name: without path and bmp/ directory\n");
       get_string(temp);
       strcat(rayscmd, "-S ");
       strcat(rayscmd, wdir);
       strcat(rayscmd, "bmp/");
       strcat(rayscmd, temp);
       strcat(rayscmd, " ");
      }
   }
 printf("Panic file selection: (\"bmp/panic.bmp\"-d,other-o) [bmp/panic.bmp]: ");
 get_char(&s);
 if (s == 'd')
   {
    strcat(rayscmd, "-P ");
    strcat(rayscmd, wdir);
    strcat(rayscmd, "bmp/panic.bmp ");
   }
 else
   {
    printf("Enter panic name: without path and bmp/ directory\n");
    get_string(temp);
    strcat(rayscmd, "-P ");
    strcat(rayscmd, wdir);
    strcat(rayscmd, "bmp/");
    strcat(rayscmd, temp);
    strcat(rayscmd, " ");
   }
 printf("Partial file selection: (\"bmp/demand.bmp\"-d,other-o) [bmp/part.bmp]: ");
 get_char(&s);
 if (s == 'd')
   {
    strcat(rayscmd, "-p ");
    strcat(rayscmd, wdir);
    strcat(rayscmd, "bmp/part.bmp ");
   }
 else
   {
    printf("Enter partial name: without path and bmp/ directory\n");
    get_string(temp);
    strcat(rayscmd, "-p ");
    strcat(rayscmd, wdir);
    strcat(rayscmd, "bmp/");
    strcat(rayscmd, temp);
    strcat(rayscmd, " ");
   }
 printf("Want to recover from a file (\"no\"-d,\"yes\"-y) [\"no\"]: ");
 get_char(&s);
 if (s == 'y')
   {
    sprintf(temp, "%s %sbmp/*.bmp", dircmd, wdir);
    printf("Files *_ab.bmp are antialiased backups\n");
    printf("Executing system command: %s\n", temp);
    system(temp);
    printf("Enter file name to recover from: without path and bmp/ directory: ");
    get_string(temp);
    strcat(rayscmd, "-R ");
    strcat(rayscmd, wdir);
    strcat(rayscmd, "bmp/");
    strcat(rayscmd, temp);
    strcat(rayscmd, " ");
   }		 
 printf("Recursion level: (6-d,N) [6]: ");
 get_int(&nr, 6);
 sprintf(temp,"-r %d ", nr);
 strcat(rayscmd, temp);
 printf("Backup per lines: (64-d,N) [64]: ");
 get_int(&nr, 64);
 sprintf(temp,"-b %d ", nr);
 strcat(rayscmd, temp);
 printf("Override screen size? (\"no\"-d,\"yes\"-y) [\"no\"]: ");
 get_char(&s);
 if (s == 'y')
   {
    printf("Screen X: (400-d,N) [400]: ");
    get_int(&nr, 400);
    sprintf(temp,"-x %d ", nr);
    strcat(rayscmd, temp);
    printf("Screen Y: (300-d,N) [300]: ");
    get_int(&nr, 300);
    sprintf(temp,"-y %d ", nr);
    strcat(rayscmd, temp);
   }
 printf("Override distorbers (ADVANCED OPTION)? (\"no\"-d,\"yes\"-y) [\"no\"]: ");
 get_char(&s);
 if (s == 'y')
   {
    printf("DistorberX: (1/17-d,N) [1/17]: ");
    get_double(&dd, 1./17.);
    sprintf(temp,"-X %f ", dd);
    strcat(rayscmd, temp);
    printf("DistorberY: (1/19-d,N) [1/19]: ");
    get_double(&dd, 1./19.);
    sprintf(temp,"-y %f ", dd);
    strcat(rayscmd, temp);
   }
 printf("Use normal distorbers? (\"yes\"-d,\"no\"-n) [\"yes\"]: ");
 get_char(&s);
 if (s != 'n')
   {
    strcat(rayscmd, "-N ");
    printf("Use global normal distorber? (\"no\"-d,\"yes\"-n) [\"no\"]: ");
    get_char(&s);
    if (s == 'y')
      {
       printf("Global normal distorber (in %%): (0-d,percent) [0]: ");
       get_double(&dd, 0.);
       sprintf(temp,"-n %f ", dd);
       strcat(rayscmd, temp);
      }
    printf("Random Seed: (current_time-d,N) [current_time]: ");
    time(&tm);
    get_int(&nr, (int)tm);
    sprintf(temp,"-s %d ", nr);
    strcat(rayscmd, temp);
   }
 printf("Use OpenGL GUI? (\"yes\"-d,\"no\"-n) [\"yes\"]: ");
 get_char(&s);
 if (s != 'n')
   {
    strcat(rayscmd, "-G ");
    printf("GL flush period (in ms): (2000-d,N) [2000]: ");
    get_int(&nr, 2000);
    sprintf(temp,"-t %d ", nr);
    strcat(rayscmd, temp);
   }
 printf("Ambient light: (0.25-d,N) [0.25]: ");
 get_double(&dd, 0.25);
 sprintf(temp,"-a %f ", dd);
 strcat(rayscmd, temp);
 printf("Minimal shadow cast: (0.2-d,N) [0.2]: ");
 get_double(&dd, 0.2);
 sprintf(temp,"-m %f ", dd);
 strcat(rayscmd, temp);
 printf("Use JPEG library? (\"yes\"-d,\"no\"-n) [\"yes\"]: ");
 get_char(&s);
 if (s != 'n')
   {
    strcat(rayscmd, "-J ");
    printf("JPEG quality: (90-d,N) [90]: ");
    get_int(&nr, 90);
    sprintf(temp,"-q %d ", nr);
    strcat(rayscmd, temp);
    printf("Save grayscale JPEG? (\"no\"-d,\"yes\"-y) [\"no\"]: ");
    get_char(&s);
    if (s == 'y') strcat(rayscmd, "-g ");
   }
 printf("Save binary scene snaphot? (\"no\"-d,\"yes\"-y) [\"no\"]: ");
 get_char(&s);
 if (s == 'y') strcat(rayscmd, "-B ");
 printf("Enable antialiased output? (\"no\"-d,\"yes\"-y) [\"no\"]: ");
 get_char(&s);
 if (s == 'y')
   {
    strcat(rayscmd, "-A ");
    printf("Double resolution for antialiasing algorithm? (\"yes\"-d,\"no\"-n) [\"yes\"]: ");
    get_char(&s);
    if (s != 'n') strcat(rayscmd, "-2 ");
    printf("Enable antialiased full backup? (\"yes\"-d,\"no\"-n) [\"yes\"]: ");
    get_char(&s);
    if (s != 'n') strcat(rayscmd, "-K ");
   }
 printf("SYSTEMCMD: \"%s\"\nPress d & ENTER to RUN", rayscmd);
 get_char(&s);
 system(rayscmd);
 fflush(stdout);
}

int main()
{
 wrapper();
 return 0;
}

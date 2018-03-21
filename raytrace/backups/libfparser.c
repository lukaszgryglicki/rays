#include "libfparser.h"
#define ID_LEN 32

char* buffer;
int position = 0;
char ch = 0;
int maxpos;
double term();
double factor();
double expression(); 
double exponential();
double x, y, z;
static int err=0;
int nvar=1;

void exception(char* reason)
{
 printf("Parser Exception: %s\n\nBuffer: '%s'", reason,buffer);
 err=1;
}

void fpar_info()
{
 printf("supported functions: sin,cos,tan,tg,ctan,ctg,exp,ln,log,sqrt,abs\n"
        "supported functions: asin,acos,asinh,acosh,sinh,cosh,cbrt,ceil,sgn\n"
        "supported functions: tanh,atan,atanh,tgh,atg,atgh\n"
        "supported operators: +,-,*,/,^, unary -\n"
        "blank characters are skipped\n");
}

int fpar_ok()
{
 nvar=1;
 fpar_f(0.);
 return !err;
}

int fpar2_ok()
{
 nvar=2;
 fpar2_f(0., 0.);
 return !err;
}

int fpar3_ok()
{
 nvar=3;
 fpar3_f(0., 0., 0.);
 return !err;
}


int fpar_function(char* tab)
{
 int i;
 if (!tab) return 1;
 if (buffer) { free(buffer); buffer = 0; }
 buffer = (char*)malloc(strlen(tab)+2);
 if (!buffer) return 1;
 for (i=0;i<strlen(tab);i++) if (tab[i] >='A' && tab[i] <= 'Z') tab[i] += 0x20;
 strcpy(buffer,tab);
 buffer[strlen(buffer)]   = 10;
 buffer[strlen(buffer)] = 0;
 maxpos = strlen(buffer);
 return 0;
}

void fpar_free()
{
 if (buffer) { free(buffer); buffer = 0; }
}


void read_next_char()
{
 if (position < maxpos && ch != 10) ch = buffer[position++];
}

void skipblanks()
{
 if (ch!=10)
  {
   while (isspace(ch) && ch!=10) read_next_char();
  }
}

void read_id(char *ident)
{
 int cnt=0;
 skipblanks();
 if (isalpha(ch))
   {
    while (isalpha(ch) || isdigit(ch))
      {
       if (cnt < ID_LEN-1) ident[cnt++] = ch;
       read_next_char();
      }
    ident[cnt] = 0;
   }
 else exception("Expected: function name or variable.\n");
}

double factor()
{
 double f, minus,tmp;
 char ident[ID_LEN],c;
 minus = 1.0;
 read_next_char();
 skipblanks();
 while (ch=='+' || ch=='-')
   {
     if (ch == '-') minus *= -1;
     read_next_char();
   }
 if (isdigit(ch) || ch=='.')
   {
    buffer[--position] = ch;
    sscanf(buffer+position, "%lf%c", &f, &ch);
    c=ch;
    do read_next_char();
    while (ch != c);
   }
 else if (ch == '(')
   {
    f = expression();
    skipblanks();
    if (ch == ')') read_next_char();
    else exception("Expected: '('.\n");
   }
 else
   {
    read_id(ident);
    if (!strcmp(ident, "pi")) f = 3.1415926353;
    if (!strcmp(ident, "x") || !strcmp(ident, "a")) f = x;
    else if (!strcmp(ident, "y") || !strcmp(ident, "b")) 
      {
       if (nvar < 2) exception("'y' found in one variable function.\n");
       else f = y;
      }
    else if (!strcmp(ident, "z") || !strcmp(ident, "c")) 
      {
       if (nvar < 3) exception("'z' found in one or two variable function.\n");
       else f = z;
      }
    else if (!strcmp(ident, "sin"))
      {
	skipblanks();
	if (ch=='(')
	  {
	   buffer[--position] = ch;
	   f = sin(exponential());
	  }
	else exception("Expected: '(' after sin.\n");
      }
     else if (!strcmp(ident, "cos"))
      {
       skipblanks();
       if (ch=='(')
	  {
	   buffer[--position] = ch;
	   f = cos(exponential());
	  }
       else exception("Expected: '(' after cos.\n");
      }
    else if (!strcmp(ident, "sinh"))
      {
	skipblanks();
	if (ch=='(')
	  {
	   buffer[--position] = ch;
	   f = sinh(exponential());
	  }
	else exception("Expected: '(' after sinh.\n");
      }
     else if (!strcmp(ident, "cosh"))
      {
       skipblanks();
       if (ch=='(')
	  {
	   buffer[--position] = ch;
	   f = cosh(exponential());
	  }
       else exception("Expected: '(' after cosh.\n");
      }
     else if (!strcmp(ident, "acos"))
      {
       skipblanks();
       if (ch=='(')
	  {
	   buffer[--position] = ch;
	   f = acos(exponential());
	  }
       else exception("Expected: '(' after acos.\n");
      }
     else if (!strcmp(ident, "asin"))
      {
       skipblanks();
       if (ch=='(')
	  {
	   buffer[--position] = ch;
	   f = asin(exponential());
	  }
       else exception("Expected: '(' after asin.\n");
      }
     else if (!strcmp(ident, "asinh"))
      {
       skipblanks();
       if (ch=='(')
	  {
	   buffer[--position] = ch;
	   f = asinh(exponential());
	  }
       else exception("Expected: '(' after asinh.\n");
      }
     else if (!strcmp(ident, "acosh"))
      {
       skipblanks();
       if (ch=='(')
	  {
	   buffer[--position] = ch;
	   f = acosh(exponential());
	  }
       else exception("Expected: '(' after acosh.\n");
      }
     else if (!strcmp(ident, "tg") || !strcmp(ident, "tan"))
      {
       skipblanks();
       if (ch=='(')
	 {
	  buffer[--position] = ch;
	  f = tan(exponential());
	 }
        else exception("Expected: '(' after tg or tan.\n");
       }
     else if (!strcmp(ident, "atg") || !strcmp(ident, "atan"))
      {
       skipblanks();
       if (ch=='(')
	 {
	  buffer[--position] = ch;
	  f = atan(exponential());
	 }
        else exception("Expected: '(' after atg or atan.\n");
       }
     else if (!strcmp(ident, "tgh") || !strcmp(ident, "tanh"))
      {
       skipblanks();
       if (ch=='(')
	 {
	  buffer[--position] = ch;
	  f = tanh(exponential());
	 }
        else exception("Expected: '(' after tgh or tanh.\n");
       }
     else if (!strcmp(ident, "atgh") || !strcmp(ident, "atanh"))
      {
       skipblanks();
       if (ch=='(')
	 {
	  buffer[--position] = ch;
	  f = atanh(exponential());
	 }
        else exception("Expected: '(' after atgh or atanh.\n");
       }
     else if (!strcmp(ident, "ctg") || !strcmp(ident,"ctan"))
       {
	skipblanks();
	if (ch=='(')
	  {
	   buffer[--position] = ch;
	   f = 1/tan(exponential());
	  }
	else exception("Expected: '(' after ctg or ctan.\n");
       }
     else if (!strcmp(ident, "exp") || !strcmp(ident,"e"))
       {
	skipblanks();
	if (ch=='(')
	  {
	   buffer[--position] = ch;
	   f = exp(exponential());
	  }
	else exception("Expected: '(' after exp or e.\n");
	}
	else if (!strcmp(ident, "ln"))
	  {
	   skipblanks();
	   if (ch=='(')
	     {
	      buffer[--position] = ch;
	      f = log(exponential());
	     }
	   else exception("Expected: '(' after ln.\n");
	  }
	else if (!strcmp(ident, "log"))
	  {
	   skipblanks();
	   if (ch=='(')
	     {
	      buffer[--position] = ch;
	      f = log10(exponential());
	     }
	   else exception("Expected: '(' after log.\n");
	  }
	else if (!strcmp(ident, "cbrt"))
	  {
	   skipblanks();
	   if (ch=='(')
	     {
	      buffer[--position] = ch;
	      f = cbrt(exponential());
	     }
	   else exception("Expected: '(' after cbrt.\n");
	  }
	else if (!strcmp(ident, "sqrt"))
	  {
	   skipblanks();
	   if (ch=='(')
	     {
	      buffer[--position] = ch;
	      f = cbrt(exponential());
	     }
	   else exception("Expected: '(' after cbrt.\n");
	  }
	else if (!strcmp(ident, "ceil"))
	  {
	   skipblanks();
	   if (ch=='(')
	     {
	      buffer[--position] = ch;
	      f = ceil(exponential());
	     }
	   else exception("Expected: '(' after ceil.\n");
	  }
	else if (!strcmp(ident, "sgn"))
	  {
	   skipblanks();
	   if (ch=='(')
	     {
	      buffer[--position] = ch;
	      tmp = exponential();
	      if (tmp > 0.0) f = 1.;
	      else if (tmp < 0.0) f = -1.;
	      else f = 0.;
	     }
	   else exception("Expected: '(' after sgn.\n");
	  }
	else if (!strcmp(ident, "abs"))
	  {
	   skipblanks();
	   if (ch=='(')
	     {
	      buffer[--position] = ch;
	      f = fabs(exponential());
	     }
	   else exception("Expected: '(' after abs.\n");
	  }
	else exception("Unknown identifier.\n");
	}
	skipblanks();
	return minus*f;		
}

double term()
{
 double f1;
 f1 = exponential();
 while(1)
   {
    switch(ch)
      {
       case '*': f1 *= exponential(); break;
       case '/': f1 /= exponential(); break;
       default: return f1;
      }
   }
}

double expression()
{
 double t1;
 t1 = term();
 while(1)
   {
    switch(ch)
      {
       case '+': t1 += term(); break;
       case '-': t1 -= term(); break;
       default: return t1;
      }
   }
}

double exponential()
{
 double f = factor();
 while (ch == '^') f = pow(f, exponential());
 return f;
}


double fpar_f(double X)
{
 double e;
 x=X;
 y=0.;
 z=0.;
 position=0;
 ch=0;
 e=expression();
 if (ch != 10 && ch != ';') exception("Garbage in function expression.\n");
 if (err) 
   {
    exception("Value returned MAY BE invalid, exceptions occured.\n");
    return e;
   }
 return e;
}

double fpar2_f(double X, double Y)
{
 double e;
 x=X;
 y=Y;
 z=0.;
 position=0;
 ch=0;
 e=expression();
 if (ch != 10 && ch != ';') exception("Garbage in function expression.\n");
 if (err) 
   {
    exception("Value returned MAY BE invalid, exceptions occured.\n");
    return e;
   }
 return e;
}

double fpar3_f(double X, double Y, double Z)
{
 double e;
 x=X;
 y=Y;
 z=Z;
 position=0;
 ch=0;
 e=expression();
 if (ch != 10 && ch != ';') exception("Garbage in function expression.\n");
 if (err) 
   {
    exception("Value returned MAY BE invalid, exceptions occured.\n");
    return e;
   }
 return e;
}


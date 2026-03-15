#include "define.h"
#include "global.h"

void initparams(PARAMETERS *params)
{
  params->seed = 0;

  params->cooling = LOGARITHMIC;
  params->rate = 48.5165;

  params->outdim = 3;
  params->indim = 2;

  params->T0 = -1;
  params->Tmin = -1;

  params->iteration = 0;

  params->retry = 10000000;

  params->maxfail = -1;
  params->maxtick = -1;
  params->maxtock = -1;

  params->acceptance = SYMMETRIC;
  params->distribution = NORMAL;
  params->movement = LOCAL;

  params->sdminimum = 0.0;
  params->sdrate = 1.0;
}

static void param_accept(char *value, PARAMETERS *params);
static void param_cool(char *value, PARAMETERS *params);
static void param_count(char *value, PARAMETERS *params);
static void param_distr(char *value, PARAMETERS *params);
static void param_fail(char *value, PARAMETERS *params);
static void param_indim(char *value, PARAMETERS *params);
static void param_iter(char *value, PARAMETERS *params);
static void param_mintemp(char *value, PARAMETERS *params);
static void param_move(char *value, PARAMETERS *params);
static void param_outdim(char *value, PARAMETERS *params);
static void param_rate(char *value, PARAMETERS *params);
static void param_report(char *value, PARAMETERS *params);
static void param_save(char *value, PARAMETERS *params);
static void param_sdmin(char *value, PARAMETERS *params);
static void param_sdrate(char *value, PARAMETERS *params);
static void param_seed(char *value, PARAMETERS *params);
static void param_temp(char *value, PARAMETERS *params);
static void param_tries(char *value, PARAMETERS *params);

typedef struct parameter PARAMETER;

struct parameter {
  char *name;
  int length;
  void (*action)();
  char *description;
};

static PARAMETER table[] = {
  { "acceptance",  1, param_accept, "[string]  Set acceptance formula:\n\t\t\t\t      b(aised):     min(1,1/exp(dE/T))\n\t\t\t\t      s(ymmetric):  1/(1+exp(dE/T))\n" },
  { "cooling",     2, param_cool,   "[string]  Select cooling rate:\n\t\t\t\t      g(eometric):    Same as exponential\n\t\t\t\t      e(xponential):  T = exp(-k*t)\n\t\t\t\t      li(near):       T = 1/(k*t)\n\t\t\t\t      lo(garithmic):  T = 1/log(1+k*t)\n" },
  { "distribution",1, param_distr,  "[string]  Set random distribution type:\n\t\t\t\t      u(niform):  Uniform distribution\n\t\t\t\t      n(ormal):   Normal distribution\n" },
  { "fail",        1, param_fail,   "[integer] Set maximum # of failures\n" },
  { "indim",       2, param_indim,  "[integer] Dimension of inside cube\n" },
  { "iteration",   2, param_iter,   "[integer] Set current iteration number\n" },
  { "mintemp",     2, param_mintemp,"[real]    Set minimum temperature\n" },
  { "movement",    2, param_move,   "[string]  Set type of movement:\n\t\t\t\t      l(ocal):   Move one point\n" },
  { "outdim",      1, param_outdim, "[integer] Dimension of outside cube\n" },
  { "rate",        2, param_rate,   "[real]    Set cooling rate, k\n" },
  { "report",      2, param_report, "[integer] Set interval between reports\n" },
  { "save",        2, param_save,   "[integer] Set interval between saves\n" },
  { "sdminimum",   3, param_sdmin,  "[real]    Set minimum standard deviation\n" },
  { "sdrate",      3, param_sdrate, "[real]    Set rate standard deviation decreased\n" },
  { "seed",        2, param_seed,   "[integer] Set seed of random number generator\n" },
  { "temperature", 2, param_temp,   "[real]    Set initial temperature\n" },
  { "tries",       2, param_tries,  "[integer] Set tries before restarting cooling\n" },
};

static int table_size = sizeof(table) / sizeof(table[0]);

static int first_char(register char *string)
{
  register int cur;

  cur = *string;

  while (cur != '\0' && (cur == ' ' || cur == '\t')) {
    string++;
    cur = *string;
  }

  return(cur);
}

void readparams(char *paramname, PARAMETERS *params)
{
  FILE *fd;
  char key[8192];
  char value[8192];
  char line[8192];
  int i;

  fd = fopen(paramname,"r");

  if (fd == NULL) {
    fprintf(stderr,"Unable to open '%s' to read parameters\n",paramname);
    exit(1);
  }

  while (fgets(line,8192,fd) != NULL) {
    line[8191] = '\0';
    line[strlen(line)-1] = '\0';

    if (first_char(line) != '#') {
      if (sscanf(line,"%s%s",key,value) != 2) {
        fprintf(stderr,"Bad line: '%s'\n",line);
        exit(1);
      }

      for (i = 0; i < table_size; i++) {
        if (strncasecmp(key,table[i].name,table[i].length) == 0) {
          (table[i].action)(value,params);
          break;
        }
      }

      if (i == table_size) {
        fprintf(stderr,"Unknown parameter: '%s'\n",line);
        exit(1);
      }
    }
  }
}

void printparams(FILE *file, PARAMETERS *params)
{
  switch (params->acceptance) {
    case BIASED:
      fprintf(file,"acceptance  biased\n");
      break;

    case SYMMETRIC:
      fprintf(file,"acceptance  symmetric\n");
      break;
  }
  switch (params->cooling) {
    case GEOMETRIC:
      fprintf(file,"cooling   geometric\n");
      break;

    case LINEAR:
      fprintf(file,"cooling   linear\n");
      break;

    case LOGARITHMIC:
      fprintf(file,"cooling   logarithmic\n");
      break;
  }
  switch (params->distribution) {
    case UNIFORM:
      fprintf(file,"distribution  uniform\n");
      break;

    case NORMAL:
      fprintf(file,"distribution  normal\n");
      break;
  }
  fprintf(file,"fail    %.0lf\n",params->maxfail);
  fprintf(file,"indim   %d\n",params->indim);
  fprintf(file,"iteration %.0lf\n",params->iteration);
  fprintf(file,"mintemp   %.27Le\n",params->Tmin);
  switch (params->movement) {
    case LOCAL:
      fprintf(file,"movement  local\n");
      break;
  }
  fprintf(file,"outdim    %d\n",params->outdim);
  fprintf(file,"rate    %.27Le\n",params->rate);
  fprintf(file,"report    %.0lf\n",params->maxtick);
  fprintf(file,"save    %.0lf\n",params->maxtock);
  fprintf(file,"sdminimum %.27Le\n",params->sdminimum);
  fprintf(file,"sdrate    %.27Le\n",params->sdrate);
  fprintf(file,"seed    %d\n",params->seed);
  fprintf(file,"temperature %.27Le\n",params->T0);
  fprintf(file,"tries   %.0lf\n",params->retry);
  fprintf(file,"\n");
}

void helpparams()
{
  int i,j;
  char *name;
  int length;
  int maxlength;

  maxlength = 0;

  for (i = 0; i < table_size; i++) {
    name = table[i].name;
    length = table[i].length;

    if (length < (int)strlen(name)) {
      if (maxlength < (int)strlen(name) + 2) {
        maxlength = strlen(name) + 2;
      }
    } else {
      if (maxlength < (int)strlen(name)) {
        maxlength = strlen(name);
      }
    }
  }

  maxlength = (maxlength/8 + 1) * 8;

  printf("\n");
  printf("Parameters:\n");

  for (i = 0; i < table_size; i++) {
    name = table[i].name;
    length = table[i].length;

    printf("\t");

    if (length < (int)strlen(name)) {
      for (j = 0; j < length; j++) {
        putchar(name[j]);
      }

      putchar('(');

      for ( ; j < (int)strlen(name); j++) {
        putchar(name[j]);
      }

      putchar(')');

      length = strlen(name) + 2;
      printf("\t");
      length = (length/8 + 1) * 8;

      while (length < maxlength) {
        printf("\t");
        length += 8;
      }
    } else {
      for (j = 0; j < (int)strlen(name); j++) {
        putchar(name[j]);
      }

      length = strlen(name);
      printf("\t");
      length = (length/8 + 1) * 8;

      while (length < maxlength) {
        printf("\t");
        length += 8;
      }
    }

    printf("%s",table[i].description);
  }

  printf("\n");
}

static void param_accept(char *value, PARAMETERS *params)
{
  if (strncasecmp(value,"b",1) == 0) {
    params->acceptance = BIASED;
  } else
  if (strncasecmp(value,"s",1) == 0) {
    params->acceptance = SYMMETRIC;
  } else {
    fprintf(stderr,"Unknown type of acceptance: %s\n",value);
    exit(1);
  }
}

static void param_cool(char *value, PARAMETERS *params)
{
  if (strncasecmp(value,"g",1) == 0) {
    params->cooling = GEOMETRIC;
  } else
  if (strncasecmp(value,"e",1) == 0) {
    params->cooling = EXPONENTIAL;
  } else
  if (strncasecmp(value,"li",2) == 0) {
    params->cooling = LINEAR;
  } else
  if (strncasecmp(value,"lo",2) == 0) {
    params->cooling = LOGARITHMIC;
  } else {
    fprintf(stderr,"Unknown type of cooling: %s\n",value);
    exit(1);
  }
}

static void param_distr(char *value, PARAMETERS *params)
{
  if (strncasecmp(value,"u",1) == 0) {
    params->distribution = UNIFORM;
  } else
  if (strncasecmp(value,"n",1) == 0) {
    params->distribution = NORMAL;
  } else {
    fprintf(stderr,"Unknown type of distribution: %s\n",value);
    exit(1);
  }
}

static void param_fail(char *value, PARAMETERS *params)
{
  params->maxfail = atof(value);
  if (params->maxfail <= 0) {
    fprintf(stderr,"Maximum failures, %.0lf, must be > 0\n",params->maxfail);
    exit(1);
  }
}

static void param_indim(char *value, PARAMETERS *params)
{
  params->indim = atoi(value);
  if (params->indim < 0) {
    fprintf(stderr,"Inner cube dimension, %d, must be >= 0\n",params->indim);
    exit(1);
  }
}

static void param_iter(char *value, PARAMETERS *params)
{
  params->iteration = atof(value);
  if (params->iteration < 0) {
    fprintf(stderr,"Current iteration, %.0lf, must be >= 0\n",params->iteration);
    exit(1);
  }
}

static void param_mintemp(char *value, PARAMETERS *params)
{
  params->Tmin = atof(value);
  if (params->Tmin <= 0) {
    fprintf(stderr,"Minimum temperature, %Lf, must be > 0\n",params->Tmin);
    exit(1);
  }
}

static void param_move(char *value, PARAMETERS *params)
{
  if (strncasecmp(value,"l",1) == 0) {
    params->movement = LOCAL;
  } else {
    fprintf(stderr,"Unknown type of movement: %s\n",value);
    exit(1);
  }
}

static void param_outdim(char *value, PARAMETERS *params)
{
  params->outdim = atoi(value);
  if (params->outdim < 0) {
    fprintf(stderr,"Outer cube dimension, %d, must be >= 0\n",params->outdim);
    exit(1);
  }
}

static void param_rate(char *value, PARAMETERS *params)
{
  params->rate = atof(value);
  if (params->rate <= 0) {
    fprintf(stderr,"Rate, %Lf, must be > 0\n",params->rate);
    exit(1);
  }
}

static void param_report(char *value, PARAMETERS *params)
{
  params->maxtick = atof(value);
  if (params->maxtick <= 0) {
    fprintf(stderr,"Report interval, %.0lf, must be > 0\n",params->maxtick);
    exit(1);
  }
}

static void param_save(char *value, PARAMETERS *params)
{
  params->maxtock = atof(value);
  if (params->maxtock <= 0) {
    fprintf(stderr,"Save interval, %.0lf, must be > 0\n",params->maxtock);
    exit(1);
  }
}

static void param_sdmin(char *value, PARAMETERS *params)
{
  params->sdminimum = atof(value);
  if (params->sdminimum < 0) {
    fprintf(stderr,"SD minimum, %Lf, must be >= 0\n",params->sdminimum);
    exit(1);
  }
}

static void param_sdrate(char *value, PARAMETERS *params)
{
  params->sdrate = atof(value);
  if (params->sdrate < 0) {
    fprintf(stderr,"SD rate, %Lf, must be >= 0\n",params->sdrate);
    exit(1);
  }
}

static void param_seed(char *value, PARAMETERS *params)
{
  params->seed = atoi(value);
}

static void param_temp(char *value, PARAMETERS *params)
{
  params->T0 = atof(value);
  if (params->T0 <= 0) {
    fprintf(stderr,"Initial temperature, %Lf, must be > 0\n",params->T0);
    exit(1);
  }
}

static void param_tries(char *value, PARAMETERS *params)
{
  params->retry = atof(value);
  if (params->retry <= 0) {
    fprintf(stderr,"Tries, %.0lf, must be > 0\n",params->retry);
    exit(1);
  }
}

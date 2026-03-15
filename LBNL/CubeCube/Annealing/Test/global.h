extern void initparams(PARAMETERS *);
extern void helpparams(void);
extern void readparams(char *, PARAMETERS *);
extern void printparams(FILE *, PARAMETERS *);
extern void newchain(CHAIN *, int, int);
extern void initchain(CHAIN *, int *, PARAMETERS *);
extern int next(CHAIN *, CHAIN *, int, REAL, PARAMETERS *, REAL *);
extern int approve(REAL, REAL, REAL, PARAMETERS *);
extern void copychain(CHAIN *, CHAIN *, int, int);
extern REAL energy(CHAIN *, int, PARAMETERS *);

extern void printchain(FILE *, CHAIN *, int, PARAMETERS *);
extern void debugchain(FILE *, CHAIN *, int, int);

extern void initpanic(PARAMETERS *, CHAIN *, CHAIN *, int *, REAL *,
		      REAL *);
extern void panic(int);
extern void endpanic(void);

extern void check(FILE *, CHAIN *, int, PARAMETERS *);

extern int sleep(unsigned);

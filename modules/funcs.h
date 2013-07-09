
void init_tconvert(void);
void assign_type(void);
int read_header(void);
//int errorcheck(unsigned int, unsigned int, char);
void skip_blocks(int);

void skiprho(void);
void skipgas(void);

void gadget_readtmax(void);
void gadget_readpotential(void);
void gadget_posvel(void);
void readpid(void);
void gadget_readnspawn(void);
void gadget_readZ(void);
void gadget_readmetals(void);
void gadget_mass(void);
void gas_props(void);
void gadget_readage(void);
void tipsy_posvel(void);
void tipsy_bin(void);
void tipsy_auxmetals(void);
void tipsy_aux(void);
void read_tipsy(void);
void read_tipsy_envira(void);
void read_tipsy_future(int, int);
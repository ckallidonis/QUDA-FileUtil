/* Christos Kallidonis                                    */
/* July 2015                                              */
/* This code takes as input a .h5 file with a three-point */
/* function in position-space, it performs the Fourier    */
/* transform and writes the result in an ASCII file       */


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>

#define MAX_MOM 5000
#define Ntype 3

void usage(char exe[]){
  printf("%s:\n<.h5-file>\n<output_dir>\n <thrp-type (ultra_local, noether, oneD)>\n<proj>\n<conf_traj>\n<src_x>\n<src_y>\n<src_z>\n<src_t>\n<tsink>\n<L>\n<T>\n<Qsq>\n",exe);
  exit(-1);
}
//=================================================================================================

int createMomenta(int **mom, int L, int Qsq){

  int totMom = 0;
  for(int iQ = 0 ; iQ <= Qsq ; iQ++){
    for(int nx = iQ ; nx >= -iQ ; nx--)
      for(int ny = iQ ; ny >= -iQ ; ny--)
        for(int nz = iQ ; nz >= -iQ ; nz--){
          if( nx*nx + ny*ny + nz*nz == iQ ){
	    if(totMom==MAX_MOM){
	      fprintf(stderr,"Total number of momenta exceeded! Exiting.\n");
	      exit(-1);
	    }
            mom[totMom][0] = nx;
            mom[totMom][1] = ny;
            mom[totMom][2] = nz;
	    printf("Mom %d: %+d %+d %+d\n",totMom+1,mom[totMom][0],mom[totMom][1],mom[totMom][2]);
            totMom++;
          }
        }
  }

  return totMom;
}
//=================================================================================================

int main(int argc, char *argv[]){

  if(argc!=14) usage(argv[0]);

  char thrp_type[3][256];
  int M,DIR;

  strcpy(thrp_type[0],"ultra_local");
  strcpy(thrp_type[1],"noether");
  strcpy(thrp_type[2],"oneD");

  enum THRP_TYPE{THRP_LOCAL,THRP_NOETHER,THRP_ONED};

  char *h5_file, *outdir, *thrp, *conf,*proj;
  int src[4];
  asprintf(&h5_file,"%s",argv[1]);
  asprintf(&outdir ,"%s",argv[2]);
  asprintf(&thrp   ,"%s",argv[3]);
  asprintf(&proj   ,"%s",argv[4]);
  asprintf(&conf   ,"%s",argv[5]);
  src[0] = atoi(argv[6]);
  src[1] = atoi(argv[7]);
  src[2] = atoi(argv[8]);
  src[3] = atoi(argv[9]);
  int tsink = atoi(argv[10]);
  int L     = atoi(argv[11]);
  int T     = atoi(argv[12]);
  int Qsq   = atoi(argv[13]);

  bool thrpOK = false;
  int dt;
  for(int i=0;(i<Ntype && !thrpOK);i++)
    if(strcmp(thrp,thrp_type[i])==0){
      thrpOK = true;
      dt = i;
    }
  THRP_TYPE type = (THRP_TYPE) dt;

  if(!thrpOK){
    printf("Error: Thrp must be one of:\n");
    for(int i=0;i<Ntype;i++) printf(" %s\n",thrp_type[i]);
    exit(-1);
  }
  
  if( type==THRP_LOCAL || type==THRP_ONED ) M = 16;
  else if (type==THRP_NOETHER) M = 4;
  else{
    printf("Undefined thrp type!\n");
    exit(-1);
  }
  
  if( type==THRP_LOCAL || type==THRP_NOETHER ) DIR = 1;
  else if (type==THRP_ONED) DIR = 4;

  printf("Got the following input:\n");
  printf("h5_file: %s\n",h5_file);
  printf("outdir: %s\n",outdir);
  printf("thrp: %d - %s\n",dt,thrp);
  printf("Projector: %s\n",proj);
  printf("conf traj: %s\n",conf);
  printf("src [x,y,z,t] = [%02d,%02d,%02d,%02d]\n",src[0],src[1],src[2],src[3]);
  printf("Tsink = %02d\n",tsink);
  printf("L     = %02d\n",L);
  printf("T     = %02d\n",T);
  printf("Qsq   = %02d\n",Qsq);
  //-----------------------------------------

  //-Open the h5 file
  hid_t file_id = H5Fopen(h5_file, H5F_ACC_RDONLY, H5P_DEFAULT);

  //-Read the dataset
  hid_t group_id;
  hid_t dset_id;
  char *group_dir;
  char *dset_name;

  int sV = L*L*L;
  int V = sV*T;

  float *thrpBuf = (float*) malloc(DIR*M*2*V*2*sizeof(float));
  if(thrpBuf == NULL){
    fprintf(stderr,"Cannot allocate thrpBuf. Exiting\n");
    exit(-1);
  }

  asprintf(&dset_name,"threep");

  asprintf(&group_dir,"/conf_%s/sx%02dsy%02dsz%02dst%02d/tsink_%d/proj_%s/%s",conf,src[0],src[1],src[2],src[3],tsink,proj,thrp);
  group_id = H5Gopen(file_id, group_dir, H5P_DEFAULT);
  dset_id = H5Dopen(group_id, dset_name, H5P_DEFAULT);
  
  herr_t status = H5Dread( dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, thrpBuf);
  
  if (status<0){
    fprintf (stderr, "Dataset read failed!\n");
    status = H5Gclose(group_id);
    status = H5Dclose(dset_id);
    status = H5Fclose(file_id);
    free(thrpBuf);
    exit(-1);
  }
  H5Dclose(dset_id);
  H5Gclose(group_id);

  printf("Dataset read successfully\n");
  //----------------------------------------------

  //-Perform the FT
  float two_pi = 4.0*asin(1.0);
  int **mom;
  
  mom = (int**) malloc(MAX_MOM*sizeof(int*));
  if(mom==NULL){
    fprintf(stderr,"Cannot allocate mom, top-level\n");
    exit(-1);
  }
  for(int i=0;i<MAX_MOM;i++){
    mom[i] = (int*) malloc(3*sizeof(int));
    if(mom[i] == NULL){
      fprintf(stderr,"Cannot allocate mom[%d]\n",i);
      exit(-1);
    }
    for(int j=0;j<3;j++) mom[i][j] = 0;
  }

  int Nmoms = createMomenta(mom,L,Qsq);
  printf("Created %d Momenta\n",Nmoms);

  float *thrpMom = (float*) malloc(DIR*M*2*T*Nmoms*2*sizeof(float));
  if(thrpMom == NULL){
    fprintf(stderr,"Cannot allocate thrpMom. Exiting\n");
    exit(-1);
  }
  memset(thrpMom,0,DIR*M*2*T*Nmoms*2*sizeof(float));
  
  for(int uORd = 0;uORd<2;uORd++){
    for(int id=0;id<DIR;id++){
      for(int ip=0;ip<Nmoms;ip++){
	int px = mom[ip][0];
	int py = mom[ip][1];
	int pz = mom[ip][2];
	
	int v = 0;
	for(int z=0;z<L;z++){
	  for(int y=0;y<L;y++){
	    for(int x=0;x<L;x++){
	      float expn = two_pi*(px*(x-src[0]) + py*(y-src[1]) + pz*(z-src[2])) / (float)L;
	      float phase[2];
	      phase[0] = cos(expn);
	      phase[1] = sin(expn);
	      for(int t=0;t<T;t++){
		for(int gm=0;gm<M;gm++){
		  thrpMom[ 0 + 2*gm + 2*M*ip + 2*M*Nmoms*t + 2*M*Nmoms*T*id + 2*M*Nmoms*T*DIR*uORd] += 
		    thrpBuf[0 + 2*v + 2*sV*t + 2*sV*T*uORd + 2*sV*T*2*gm + 2*sV*T*2*M*id]*phase[0] - thrpBuf[1 + 2*v + 2*sV*t + 2*sV*T*uORd + 2*sV*T*2*gm + 2*sV*T*2*M*id]*phase[1];
		  
		  thrpMom[ 1 + 2*gm + 2*M*ip + 2*M*Nmoms*t + 2*M*Nmoms*T*id + 2*M*Nmoms*T*DIR*uORd] += 
		    thrpBuf[0 + 2*v + 2*sV*t + 2*sV*T*uORd + 2*sV*T*2*gm + 2*sV*T*2*M*id]*phase[1] + thrpBuf[1 + 2*v + 2*sV*t + 2*sV*T*uORd + 2*sV*T*2*gm + 2*sV*T*2*M*id]*phase[0];			   
		}//-gm
	      }//-t
	      v++;
	    }//-x
	  }//-y
	}//-z
	printf("imom = %d done\n",ip);
      }//-ip
    }//-id
  }//-uORd
  printf("Fourier Transform completed successfully\n");
  //----------------------------------------------
   
  //-Write the output files
  char *file[2];
  char *Parts[2];

  asprintf(&(Parts[0]),"up");
  asprintf(&(Parts[1]),"down");

  for(int i=0;i<2;i++) asprintf(&(file[i]),"%s/threep.%s_tsink%d_proj%s.neutron.%s.%s.SS.%02d.%02d.%02d.%02d.dat",outdir, conf, tsink, proj, Parts[i], thrp, src[0], src[1], src[2], src[3]);

  FILE *outp[2];
  for(int i=0;i<2;i++){
    if( (outp[i]=fopen(file[i],"w")) == NULL ){
      fprintf(stderr,"Cannot open %s for writing\n",file[i]);
      exit(-1);
    }
  }

  if(type==THRP_LOCAL || type==THRP_NOETHER){
    int id = 0;
    for(int uORd=0;uORd<2;uORd++){    
      for(int ip=0;ip<Nmoms;ip++){
	for(int t=0;t<=tsink;t++){
	  int ts = (t + src[3])%T;
	  for(int gm=0;gm<M;gm++){
	    fprintf(outp[uORd],"%02d\t%02d\t%+d %+d %+d\t%+e %+e\n",t,gm,mom[ip][0],mom[ip][1],mom[ip][2],
		    thrpMom[ 0 + 2*gm + 2*M*ip + 2*M*Nmoms*ts + 2*M*Nmoms*T*id + 2*M*Nmoms*T*DIR*uORd], thrpMom[ 1 + 2*gm + 2*M*ip + 2*M*Nmoms*ts + 2*M*Nmoms*T*id + 2*M*Nmoms*T*DIR*uORd]);
	  }}}
    }
  }
  else if(type==THRP_ONED){
    for(int uORd=0;uORd<2;uORd++){    
      for(int ip=0;ip<Nmoms;ip++){
	for(int id=0;id<DIR;id++){
	  for(int t=0;t<=tsink;t++){
	    int ts = (t + src[3])%T;
	    for(int gm=0;gm<M;gm++){
	      fprintf(outp[uORd],"%02d\t%02d\t%02d\t%+d %+d %+d\t%+e %+e\n",id,t,gm,mom[ip][0],mom[ip][1],mom[ip][2],
		      thrpMom[ 0 + 2*gm + 2*M*ip + 2*M*Nmoms*ts + 2*M*Nmoms*T*id + 2*M*Nmoms*T*DIR*uORd], thrpMom[ 1 + 2*gm + 2*M*ip + 2*M*Nmoms*ts + 2*M*Nmoms*T*id + 2*M*Nmoms*T*DIR*uORd]);
	    }}}}
    }
  }
  

  free(thrpBuf);
  free(thrpMom);
  for(int i=0;i<Nmoms;i++) free(mom[i]);
  free(mom);
  for(int i=0;i<2;i++) fclose(outp[i]);

  printf("Extracting Three-point function and FT completed successfully.\n");

  return 0;
}

/* Christos Kallidonis                                  */
/* July 2015                                            */
/* This code takes as input a .h5 file with a baryon    */
/* two-point function in position-space, it performs    */
/* the Fourier transform and writes the result in an    */
/* ASCII file                                           */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>

#define Ntype 10
#define M 16
#define MAX_MOM 5000

void usage(char exe[]){
  printf("%s:\n<.h5-file>\n<output_dir>\n<baryon-type>\n  nucl_nucl\n  nucl_roper\n  roper_nucl\n  roper_roper\n  deltapp_deltamm_11\n  deltapp_deltamm_22\n  deltapp_deltamm_33\n  deltap_deltaz_11\n  deltap_deltaz_22\n  deltap_deltaz_33\n  all\n",exe);
  printf("<conf_traj>\n<src_x>\n<src_y>\n<src_z>\n<src_t>\n<L>\n<T>\n<Qsq>\n");
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

// int momentaCounter(int Qsq){
//   int counter = 0;

//   for(int iQ = 0 ; iQ <= Qsq ; iQ++){
//     for(int nx = iQ ; nx >= -iQ ; nx--)
//       for(int ny = iQ ; ny >= -iQ ; ny--)
//         for(int nz = iQ ; nz >= -iQ ; nz--){
//           if( nx*nx + ny*ny + nz*nz == iQ ){
//             counter++;
//           }
//         }
//   }

//   return counter;
// }
//=================================================================================================

// void createMomenta(int **mom, int L, int Qsq){

//   int totMom = 0;
//   int Gmom[3];

//   for(int pz = 0; pz < L; pz++){
//     for(int py = 0; py < L; py++){
//       for(int px = 0; px < L; px++){
//         if(px < L/2) Gmom[0] = px;
//         else         Gmom[0] = px - L;

//         if(py < L/2) Gmom[1] = py;
//         else         Gmom[1] = py - L;

//         if(pz < L/2) Gmom[2] = pz;
//         else         Gmom[2] = pz - L;

//         if( (Gmom[0]*Gmom[0] + Gmom[1]*Gmom[1] + Gmom[2]*Gmom[2]) <= Qsq ){
// 	  for(int i=0;i<3;i++) mom[totMom][i] = Gmom[i];
//           printf("Mom %d: %+d %+d %+d\n",totMom+1,mom[totMom][0],mom[totMom][1],mom[totMom][2]);
//           totMom++;
// 	}
//       }}
//   }

// }
//=================================================================================================

int main(int argc, char *argv[]){

  if(argc!=12) usage(argv[0]);

  char twop_type[Ntype+1][256];

  strcpy(twop_type[0],"nucl_nucl");
  strcpy(twop_type[1],"nucl_roper");
  strcpy(twop_type[2],"roper_nucl");
  strcpy(twop_type[3],"roper_roper");
  strcpy(twop_type[4],"deltapp_deltamm_11");
  strcpy(twop_type[5],"deltapp_deltamm_22");
  strcpy(twop_type[6],"deltapp_deltamm_33");
  strcpy(twop_type[7],"deltap_deltaz_11");
  strcpy(twop_type[8],"deltap_deltaz_22");
  strcpy(twop_type[9],"deltap_deltaz_33");
  strcpy(twop_type[10],"all");

  char *h5_file, *outdir, *twop, *conf;
  int src[4];
  asprintf(&h5_file,"%s",argv[1]);
  asprintf(&outdir ,"%s",argv[2]);
  asprintf(&twop   ,"%s",argv[3]);
  asprintf(&conf   ,"%s",argv[4]);
  src[0] = atoi(argv[5]);
  src[1] = atoi(argv[6]);
  src[2] = atoi(argv[7]);
  src[3] = atoi(argv[8]);
  int L   = atoi(argv[9]);
  int T   = atoi(argv[10]);
  int Qsq = atoi(argv[11]);

  bool twopOK = false;
  int dt;
  for(int i=0;(i<=Ntype && !twopOK);i++)
    if(strcmp(twop,twop_type[i])==0){
      twopOK = true;
      dt = i;
    }

  if(!twopOK){
    printf("Error: Twop must be one of:\n");
    for(int i=0;i<Ntype;i++) printf(" %s\n",twop_type[i]);
    exit(-1);
  }

  printf("Got the following input:\n");
  printf("h5_file: %s\n",h5_file);
  printf("outdir: %s\n",outdir);
  printf("twop: %d - %s\n",dt,twop);
  printf("conf traj: %s\n",conf);
  printf("src [x,y,z,t] = [%02d,%02d,%02d,%02d]\n",src[0],src[1],src[2],src[3]);
  printf("L   = %02d\n",L);
  printf("T   = %02d\n",T);
  printf("Qsq = %02d\n",Qsq);
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
  int Np, Ns;
  bool READ_ALL;
  if(strcmp(twop,"all")==0){
    READ_ALL = true;
    Np = Ntype;
    Ns = 0;
  }
  else{
    READ_ALL= false;
    Np = 1;
    Ns = dt;
  }

  float *twopBuf = (float*) malloc(Np*2*V*M*2*sizeof(float));
  if(twopBuf == NULL){
    fprintf(stderr,"Cannot allocate twopBuf. Exiting\n");
    exit(-1);
  }

  asprintf(&dset_name,"twop-baryon");

  for(int bar=Ns;bar<(Np+Ns);bar++){
    asprintf(&group_dir,"/conf_%s/sx%02dsy%02dsz%02dst%02d/%s",conf,src[0],src[1],src[2],src[3],twop_type[bar]);
    group_id = H5Gopen(file_id, group_dir, H5P_DEFAULT);
    dset_id = H5Dopen(group_id, dset_name, H5P_DEFAULT);
    
    herr_t status = H5Dread( dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(twopBuf[(bar%Np)*2*V*M*2]) );

    if (status<0){
      fprintf (stderr, "Dataset read failed!\n");
      status = H5Gclose(group_id);
      status = H5Dclose(dset_id);
      status = H5Fclose(file_id);
      free(twopBuf);
      exit(-1);
    }
    H5Dclose(dset_id);
    H5Gclose(group_id);
  }
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

  float *twopMom = (float*) malloc(2*Np*Nmoms*T*M*2*sizeof(float));
  if(twopMom == NULL){
    fprintf(stderr,"Cannot allocate twopMom. Exiting\n");
    exit(-1);
  }
  memset(twopMom,0,2*Np*Nmoms*T*M*2*sizeof(float));
  
  for(int bar=Ns;bar<(Np+Ns);bar++){
    int bidx = (bar%Np);
    printf("%d - Transforming baryon %s...\n",bidx,twop_type[bar]);
    for(int pr=0;pr<2;pr++){
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
	      phase[0] =  cos(expn);
	      phase[1] = -sin(expn);
	      for(int t=0;t<T;t++){
		for(int gm=0;gm<M;gm++){
		  twopMom[ 0 + 2*gm + 2*M*ip + 2*M*Nmoms*t + 2*M*Nmoms*T*pr + 2*M*Nmoms*T*2*bidx] += 
		    twopBuf[0 + 2*gm + 2*M*v + 2*M*sV*t + 2*M*sV*T*pr + 2*M*sV*T*2*bidx]*phase[0] - twopBuf[1 + 2*gm + 2*M*v + 2*M*sV*t + 2*M*sV*T*pr + 2*M*sV*T*2*bidx]*phase[1];
		  
		  twopMom[ 1 + 2*gm + 2*M*ip + 2*M*Nmoms*t + 2*M*Nmoms*T*pr + 2*M*Nmoms*T*2*bidx] += 
		    twopBuf[0 + 2*gm + 2*M*v + 2*M*sV*t + 2*M*sV*T*pr + 2*M*sV*T*2*bidx]*phase[1] + twopBuf[1 + 2*gm + 2*M*v + 2*M*sV*t + 2*M*sV*T*pr + 2*M*sV*T*2*bidx]*phase[0];			   
		}//-gm
	      }//-t
	      v++;
	    }//-x
	  }//-y
	}//-z
      }//-ip
      printf(" flav %d\n",pr+1);
    }//-pr
    printf("Baryon %s done.\n",twop_type[bar]);
  }//-bar
  printf("Fourier Transform completed successfully\n");
  //----------------------------------------------
  
  
  //-Write the output file
  char *file;
  asprintf(&file,"%s/twop.%s.%s.SS.%02d.%02d.%02d.%02d.dat",outdir, conf, READ_ALL ? "baryons" : twop_type[dt], src[0], src[1], src[2], src[3]);

  FILE *outp;
  if( (outp=fopen(file,"w")) == NULL ){
    fprintf(stderr,"Cannot open %s for writing\n",file);
    exit(-1);
  }


  for(int bar=Ns;bar<(Np+Ns);bar++){
    int bidx = (bar%Np);
    for(int ip=0;ip<Nmoms;ip++){
      for(int t=0;t<T;t++){
	int ts = (t + src[3])%T;
	for(int ga=0;ga<4;ga++){
	  for(int gap=0;gap<4;gap++){
	    int gm = gap+4*ga;
	    fprintf(outp,"%d\t%02d\t%+d %+d %+d\t%d %d\t%+e %+e\t%+e %+e\n",bar,t,mom[ip][0],mom[ip][1],mom[ip][2],ga,gap,
		    twopMom[ 0 + 2*gm + 2*M*ip + 2*M*Nmoms*ts + 2*M*Nmoms*T*0 + 2*M*Nmoms*T*2*bidx], twopMom[ 1 + 2*gm + 2*M*ip + 2*M*Nmoms*ts + 2*M*Nmoms*T*0 + 2*M*Nmoms*T*2*bidx],
		    twopMom[ 0 + 2*gm + 2*M*ip + 2*M*Nmoms*ts + 2*M*Nmoms*T*1 + 2*M*Nmoms*T*2*bidx], twopMom[ 1 + 2*gm + 2*M*ip + 2*M*Nmoms*ts + 2*M*Nmoms*T*1 + 2*M*Nmoms*T*2*bidx]);
	  }}}}
  }


  free(twopBuf);
  free(twopMom);
  for(int i=0;i<Nmoms;i++) free(mom[i]);
  free(mom);
  fclose(outp);

  printf("Extracting Two-point function and FT completed successfully.\n");

  return 0;
}

/* Christos Kallidonis                              */
/* February 2015                                    */
/* This code takes as input a .h5 file and extracts */
/* the desired loop data into an ASCII file         */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>

#define M 16
#define DIR 4

void usage(char exe[]){
  printf("%s: <.h5-file> <output_file> <loop (Scalar, dOp, Loops, LoopsCv, LpsDw, LpsDwCv)> <momx> <momy> <momz> <conf_traj> <Nvec> <HighPrec:0, LowPrec:1> <T>\n",exe);
  exit(-1);
}

int main(int argc, char *argv[]){

  if(argc!=11) usage(argv[0]);

  char loop_type[6][256];
  bool loop_oneD[6];
  strcpy(loop_type[0],"Scalar");   loop_oneD[0] = false;   // std-ultra_local
  strcpy(loop_type[1],"dOp");      loop_oneD[1] = false;   // gen-ultra_local
  strcpy(loop_type[2],"Loops");    loop_oneD[2] = true;    // std-one_derivative
  strcpy(loop_type[3],"LoopsCv");  loop_oneD[3] = true;    // gen-one_derivative
  strcpy(loop_type[4],"LpsDw");    loop_oneD[4] = true;    // std-conserved current
  strcpy(loop_type[5],"LpsDwCv");  loop_oneD[5] = true;    // gen-conserved current

  char *h5_file, *out_file, *loop,*mom[3],*conf;
  asprintf(&h5_file ,"%s",argv[1]);
  asprintf(&out_file,"%s",argv[2]);
  asprintf(&loop    ,"%s",argv[3]);
  asprintf(&mom[0]  ,"%s",argv[4]);
  asprintf(&mom[1]  ,"%s",argv[5]);
  asprintf(&mom[2]  ,"%s",argv[6]);
  asprintf(&conf    ,"%s",argv[7]);  

  int Ns = atoi(argv[8]);
  bool LowPrecSum = (atoi(argv[9]) == 1 ? true : false);
  int T = atoi(argv[10]);

  bool loopOK = false;
  int dt;
  for(int i=0;(i<6 && !loopOK);i++)
    if(strcmp(loop,loop_type[i])==0){
      loopOK = true;
      dt = i;
    }

  if(!loopOK){
    printf("Error: Loop must be one of:\n");
    for(int i=0;i<6;i++) printf(" %s\n",loop_type[i]);
    exit(-1);
  }
  
  printf("Got the following input:\n");
  printf("h5_file: %s\n",h5_file);
  printf("out_file: %s\n",out_file);
  printf("loop: %d - %s\n",dt,loop);
  printf("momentum: %s %s %s\n",mom[0],mom[1],mom[2]);
  printf("conf traj: %s\n",conf);
  printf("Nvec: %04d, %s\n",Ns, LowPrecSum ? "LowPrecSum" : "HighPrecSum");
  printf("T = %d\n",T);
  //-----------------------------------------

  //-Open the h5 file
  hid_t file_id = H5Fopen(h5_file, H5F_ACC_RDONLY, H5P_DEFAULT);

  //-Open the desired dataset and write it in the ASCII file
  hid_t group_id;
  hid_t dset_id;
  char *group_dir;

  double *loopBuf = (double*) malloc(T*M*2*sizeof(double));

  FILE *out = fopen(out_file,"w");
  if(out==NULL){
    fprintf(stderr,"Cannot open %s for writing. Exiting.\n",out_file);
    exit(-1);
  }

  char *vec_str;
  if(LowPrecSum) asprintf(&vec_str,"NLP");
  else asprintf(&vec_str,"NHP");

  
  if(loop_oneD[dt]){
    for(int mu=0;mu<DIR;mu++){
      asprintf(&group_dir,"/conf_%s/%s_%04d/%s/mom_xyz_%s_%s_%s/dir_%02d",conf,vec_str,Ns,loop,mom[0],mom[1],mom[2],mu);
      group_id = H5Gopen(file_id, group_dir, H5P_DEFAULT);
      dset_id = H5Dopen(group_id, "loop", H5P_DEFAULT);
    
      herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, loopBuf);

      if (status<0){
	fprintf (stderr, "Dataset read failed!\n");
	status = H5Gclose(group_id);
	status = H5Dclose(dset_id);
	status = H5Fclose(file_id);
	fclose(out);
	free(loopBuf);
	exit(-1);
      }

      for(int t=0;t<T;t++)
	for(int m=0;m<M;m++) fprintf(out,"%02d %02d %02d %s %s %s %+16.15e %+16.15e\n",t,m,mu,mom[0],mom[1],mom[2],loopBuf[0+2*m+2*M*t],loopBuf[1+2*m+2*M*t]); 

      H5Dclose(dset_id);
      H5Gclose(group_id);
    }//-mu
  }
  else{
    asprintf(&group_dir,"/conf_%s/%s_%04d/%s/mom_xyz_%s_%s_%s",conf,vec_str,Ns,loop,mom[0],mom[1],mom[2]);
    group_id = H5Gopen(file_id, group_dir, H5P_DEFAULT);    
    dset_id = H5Dopen(group_id, "loop", H5P_DEFAULT);
    
    herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, loopBuf);

    if (status<0){
      fprintf (stderr, "Dataset read failed!\n");
      status = H5Gclose(group_id);
      status = H5Dclose(dset_id);
      status = H5Fclose(file_id);
      fclose(out);
      free(loopBuf);
      exit(-1);
    }

    for(int t=0;t<T;t++)
      for(int m=0;m<M;m++) fprintf(out,"%02d %02d %s %s %s %+16.15e %+16.15e\n",t,m,mom[0],mom[1],mom[2],loopBuf[0+2*m+2*M*t],loopBuf[1+2*m+2*M*t]); 

    H5Dclose(dset_id);
    H5Gclose(group_id);
  }

  H5Fclose(file_id);

  fclose(out);
  free(loopBuf);

  printf("Extracting Loop completed successfully.\n");

  return 0;
}

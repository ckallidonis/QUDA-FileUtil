/* Christos Kallidonis                                */
/* April 2017                                         */
/* This program reads disconnected loops written      */
/* in ASCII format and converts them into HDF5 format */
/* This is for the exact part of the loop!            */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#define M 16
#define DIR 4
#define LTYPE 6

typedef struct{
  int Nmoms,Qsq,T,conf;
  char loop_type[LTYPE][256];
  bool loop_oneD[LTYPE];
} Info;

void usage(char exe[]){
  printf("%s: <.h5-file prefix> <ASCII loop prefix> <mom_list> <NGPU> <Nmoms> <Qsq> <T> <conf-trajectory>\n",exe);
  exit(-1);
}
//=============================================================

void errorMsg(const char msg[]){
  fprintf(stderr,"%s",msg);
  exit(-1);
}
//=============================================================

void getWriteBuf(double *writeBuf, double *loopBuf, int T, int Nmoms, int imom){

  for(int it=0;it<T;it++){
    for(int im=0;im<M;im++){
      writeBuf[0+2*im+2*M*it] = loopBuf[0+2*imom+2*Nmoms*it+2*Nmoms*T*im];
      writeBuf[1+2*im+2*M*it] = loopBuf[1+2*imom+2*Nmoms*it+2*Nmoms*T*im];
    }
  }
}
//=============================================================

void writeLoops_HDF5(double *buf_std_uloc, double *buf_gen_uloc, double **buf_std_oneD, double **buf_std_csvC, double **buf_gen_oneD, double **buf_gen_csvC, char *file_pref, int **momQsq, Info loopInfo){

  char fname[512];

  sprintf(fname,"%s.%04d_exact_Qsq%d.h5",file_pref,loopInfo.conf,loopInfo.Qsq);

  double *loopBuf = NULL;
  double *writeBuf = (double*) malloc(loopInfo.T*M*2*sizeof(double));

  //  hsize_t start[3]  = {0,0,0};

  // Dimensions of the dataspace
  hsize_t dims[3]  = {(hsize_t)loopInfo.T, M, 2}; // Global
  //  hsize_t ldims[3] = {loopInfo.T, M, 2}; // Local

  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  char *group1_tag;
  asprintf(&group1_tag,"conf_%04d",loopInfo.conf);
  hid_t group1_id = H5Gcreate(file_id, group1_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t group2_id;
  hid_t group3_id;
  hid_t group4_id;
    
  for(int it=0;it<LTYPE;it++){
    char *group2_tag;
    asprintf(&group2_tag,"%s",loopInfo.loop_type[it]);
    group2_id = H5Gcreate(group1_id, group2_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for(int imom=0;imom<loopInfo.Nmoms;imom++){
      char *group3_tag;
      asprintf(&group3_tag,"mom_xyz_%+d_%+d_%+d",momQsq[imom][0],momQsq[imom][1],momQsq[imom][2]);
      group3_id = H5Gcreate(group2_id, group3_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      if(loopInfo.loop_oneD[it]){
	for(int mu=0;mu<4;mu++){
	  if(strcmp(loopInfo.loop_type[it],"Loops")==0)   loopBuf = buf_std_oneD[mu];
	  if(strcmp(loopInfo.loop_type[it],"LoopsCv")==0) loopBuf = buf_std_csvC[mu];
	  if(strcmp(loopInfo.loop_type[it],"LpsDw")==0)   loopBuf = buf_gen_oneD[mu];
	  if(strcmp(loopInfo.loop_type[it],"LpsDwCv")==0) loopBuf = buf_gen_csvC[mu];

	  char *group4_tag;
	  asprintf(&group4_tag,"dir_%02d",mu);
	  group4_id = H5Gcreate(group3_id, group4_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	  hid_t filespace  = H5Screate_simple(3, dims, NULL);
	  hid_t dataset_id = H5Dcreate(group4_id, "loop", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  //            filespace = H5Dget_space(dataset_id);

	  getWriteBuf(writeBuf, loopBuf, loopInfo.T, loopInfo.Nmoms, imom);

	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, filespace, H5P_DEFAULT, writeBuf);
	  if(status<0) errorMsg("File writing NOT SUCCESSFUL.\n");

	  H5Dclose(dataset_id);
	  H5Sclose(filespace);
	  H5Gclose(group4_id);
	}//-mu
      }//-if
      else{
	if(strcmp(loopInfo.loop_type[it],"Scalar")==0) loopBuf = buf_std_uloc;
	if(strcmp(loopInfo.loop_type[it],"dOp")==0)    loopBuf = buf_gen_uloc;

	hid_t filespace  = H5Screate_simple(3, dims, NULL);
	hid_t dataset_id = H5Dcreate(group3_id, "loop", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//          filespace = H5Dget_space(dataset_id);

	getWriteBuf(writeBuf, loopBuf, loopInfo.T, loopInfo.Nmoms, imom);

	herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, filespace, H5P_DEFAULT, writeBuf);
	if(status<0) errorMsg("File writing NOT SUCCESSFUL.\n");

	H5Dclose(dataset_id);
	H5Sclose(filespace);
      }
      H5Gclose(group3_id);
    }//-imom
    H5Gclose(group2_id);
  }//-it

  H5Gclose(group1_id);
  H5Fclose(file_id);

  free(writeBuf);
}
//=============================================================


int main(int argc, char *argv[]){

  if(argc!=9) usage(argv[0]);

  //-Read and define input
  char loop_type[LTYPE][256];
  bool loop_oneD[LTYPE];
  strcpy(loop_type[0],"Scalar");   loop_oneD[0] = false;   // std-ultra_local
  strcpy(loop_type[1],"dOp");      loop_oneD[1] = false;   // gen-ultra_local
  strcpy(loop_type[2],"Loops");    loop_oneD[2] = true;    // std-one_derivative
  strcpy(loop_type[3],"LoopsCv");  loop_oneD[3] = true;    // std-conserved current
  strcpy(loop_type[4],"LpsDw");    loop_oneD[4] = true;    // gen-one_derivative
  strcpy(loop_type[5],"LpsDwCv");  loop_oneD[5] = true;    // gen-conserved current

  char *h5_file,*in_loop,*mom_list;
  asprintf(&h5_file,"%s" ,argv[1]);
  asprintf(&in_loop,"%s"  ,argv[2]);
  asprintf(&mom_list,"%s",argv[3]);
  int NGPU   = atoi(argv[4]);
  int Nmoms  = atoi(argv[5]);
  int Qsq    = atoi(argv[6]);
  int T      = atoi(argv[7]);
  int conf   = atoi(argv[8]);

  if(T%NGPU != 0) errorMsg("Error: NGPU MUST divide T exactly. Exiting.\n");
  int locT = T/NGPU;

  printf("Got the following input:\n");
  printf("h5_file prefix: %s\n",h5_file);
  printf("ASCII in_loop: %s\n",in_loop);
  printf("momenta list: %s\n",mom_list);
  printf("NGPU   = %d\n",NGPU);
  printf("Nmoms  = %d\n",Nmoms);
  printf("Qsq    = %d\n",Qsq);
  printf("T      = %d\n",T);
  printf("conf   = %d\n",conf);
  printf("locT   = %d\n",locT);
  //-------------------------------------

  //-Allocate the write buffers
  double *buf_std_uloc,*buf_gen_uloc,*tmp_buf_uloc;
  double **buf_std_oneD,**buf_gen_oneD,**buf_std_csvC,**buf_gen_csvC,**tmp_buf_oneD;

  size_t buf_bytes = 2*Nmoms*T*M*sizeof(double);

  if( (tmp_buf_uloc = (double*) malloc(buf_bytes))==NULL ) errorMsg("Allocation of buffer tmp_buf_uloc failed.\n");
  if( (buf_std_uloc = (double*) malloc(buf_bytes))==NULL ) errorMsg("Allocation of buffer buf_std_uloc failed.\n");
  if( (buf_gen_uloc = (double*) malloc(buf_bytes))==NULL ) errorMsg("Allocation of buffer buf_gen_uloc failed.\n");

  if( (tmp_buf_oneD = (double**) malloc(DIR*sizeof(double*)))==NULL ) errorMsg("Allocation of buffer tmp_buf_oneD failed.\n");
  if( (buf_std_oneD = (double**) malloc(DIR*sizeof(double*)))==NULL ) errorMsg("Allocation of buffer buf_std_oneD failed.\n");
  if( (buf_gen_oneD = (double**) malloc(DIR*sizeof(double*)))==NULL ) errorMsg("Allocation of buffer buf_gen_oneD failed.\n");
  if( (buf_std_csvC = (double**) malloc(DIR*sizeof(double*)))==NULL ) errorMsg("Allocation of buffer buf_std_csvC failed.\n");
  if( (buf_gen_csvC = (double**) malloc(DIR*sizeof(double*)))==NULL ) errorMsg("Allocation of buffer buf_gen_csvC failed.\n");
  for(int mu = 0; mu < DIR ; mu++){
    if( (tmp_buf_oneD[mu] = (double*) malloc(buf_bytes))==NULL ) errorMsg("Allocation of buffer tmp_buf_oneD inside mu-loop failed.\n");
    if( (buf_std_oneD[mu] = (double*) malloc(buf_bytes))==NULL ) errorMsg("Allocation of buffer buf_std_oneD inside mu-loop failed.\n");
    if( (buf_gen_oneD[mu] = (double*) malloc(buf_bytes))==NULL ) errorMsg("Allocation of buffer buf_gen_oneD inside mu-loop failed.\n");
    if( (buf_std_csvC[mu] = (double*) malloc(buf_bytes))==NULL ) errorMsg("Allocation of buffer buf_std_csvC inside mu-loop failed.\n");
    if( (buf_gen_csvC[mu] = (double*) malloc(buf_bytes))==NULL ) errorMsg("Allocation of buffer buf_gen_csvC inside mu-loop failed.\n");
  }
  //-------------------------------------

  //-Allocate and read the momenta
  int **mom;
  if( (mom = (int**) malloc(Nmoms*sizeof(int*)))==NULL ) errorMsg("Error in allocating mom.\n");
  for(int ip=0; ip<Nmoms; ip++)
    if( (mom[ip] = (int*) malloc(3*sizeof(int)))==NULL ) errorMsg("Error in allocating momQsq inside ip-loop.\n");

  FILE *p_mom;
  if( (p_mom=fopen(mom_list,"r"))==NULL ) errorMsg("Cannot open momenta list for reading. Exiting.\n");

  for(int imom=0;imom<Nmoms;imom++) fscanf(p_mom,"%d %d %d\n",&mom[imom][0],&mom[imom][1],&mom[imom][2]);
  fclose(p_mom);
  //-------------------------------------

  //-Read the loops from the ASCII files
  int idum;
  char *ASCII_file;
  FILE *p_loop;

  printf("Reading ASCII Loops...\n");
  for(int lt=0;lt<LTYPE;lt++){

    for(int iGPU=0;iGPU<NGPU;iGPU++){
      asprintf(&ASCII_file,"%s_%s.loop.%d_%d",in_loop,loop_type[lt],NGPU,iGPU);
      printf("Reading loop file %s\n",ASCII_file);
      if( (p_loop = fopen(ASCII_file,"r"))==NULL ) errorMsg("Cannot open loop file for reading. Exiting.\n");
      
      if(loop_oneD[lt]){
	for(int mu=0;mu<DIR;mu++){
	  for(int imom=0;imom<Nmoms;imom++){  
	    for(int it=0;it<locT;it++){
	      int gT = it + locT*iGPU;
	      for(int im=0;im<M;im++){
		fscanf(p_loop,"%d %d %d %d %d %d %lf %lf\n",&idum,&idum,&idum,&idum,&idum,&idum,
		       &tmp_buf_oneD[mu][0+2*imom+2*Nmoms*gT+2*Nmoms*T*im], &tmp_buf_oneD[mu][1+2*imom+2*Nmoms*gT+2*Nmoms*T*im]);
	      }//-im
	    }//-it
	  }//-imom
	}//-mu
      }
      else{
	for(int imom=0;imom<Nmoms;imom++){  
	  for(int it=0;it<locT;it++){
	    int gT = it + locT*iGPU;
	    for(int im=0;im<M;im++){
	      fscanf(p_loop,"%d %d %d %d %d %lf %lf\n",&idum,&idum,&idum,&idum,&idum,
		     &tmp_buf_uloc[0+2*imom+2*Nmoms*gT+2*Nmoms*T*im], &tmp_buf_uloc[1+2*imom+2*Nmoms*gT+2*Nmoms*T*im]);
	    }//-im
	  }//-it
	}//-imom
      }
      
    }//-iGPU
    if(strcmp(loop_type[lt],"Scalar")==0) memcpy(buf_std_uloc, tmp_buf_uloc, buf_bytes);
    if(strcmp(loop_type[lt],"dOp")==0)    memcpy(buf_gen_uloc, tmp_buf_uloc, buf_bytes);

    if(strcmp(loop_type[lt],"Loops")==0)   for(int mu=0;mu<DIR;mu++) memcpy(buf_std_oneD[mu], tmp_buf_oneD[mu], buf_bytes);
    if(strcmp(loop_type[lt],"LoopsCv")==0) for(int mu=0;mu<DIR;mu++) memcpy(buf_std_csvC[mu], tmp_buf_oneD[mu], buf_bytes);
    if(strcmp(loop_type[lt],"LpsDw")==0)   for(int mu=0;mu<DIR;mu++) memcpy(buf_gen_oneD[mu], tmp_buf_oneD[mu], buf_bytes);
    if(strcmp(loop_type[lt],"LpsDwCv")==0) for(int mu=0;mu<DIR;mu++) memcpy(buf_gen_csvC[mu], tmp_buf_oneD[mu], buf_bytes);
    printf("%s done.\n",loop_type[lt]);
  }//-lt      
  printf("Finished reading ASCII Loops\n");
  //-------------------------------------

  //-Write the loops in HDF5 format
  Info loopInfo;
  loopInfo.conf = conf; loopInfo.T = T;
  loopInfo.Nmoms = Nmoms; loopInfo.Qsq  = Qsq;

  strcpy(loopInfo.loop_type[0],loop_type[0]);  loopInfo.loop_oneD[0] = loop_oneD[0];   // std-ultra_local
  strcpy(loopInfo.loop_type[1],loop_type[1]);  loopInfo.loop_oneD[1] = loop_oneD[1];   // gen-ultra_local
  strcpy(loopInfo.loop_type[2],loop_type[2]);  loopInfo.loop_oneD[2] = loop_oneD[2];   // std-one_derivative
  strcpy(loopInfo.loop_type[3],loop_type[3]);  loopInfo.loop_oneD[3] = loop_oneD[3];   // std-conserved current
  strcpy(loopInfo.loop_type[4],loop_type[4]);  loopInfo.loop_oneD[4] = loop_oneD[4];   // gen-one_derivative
  strcpy(loopInfo.loop_type[5],loop_type[5]);  loopInfo.loop_oneD[5] = loop_oneD[5];   // gen-conserved current

  writeLoops_HDF5(buf_std_uloc, buf_gen_uloc, buf_std_oneD, buf_std_csvC, buf_gen_oneD, buf_gen_csvC, h5_file, mom, loopInfo);
  printf("Writing in HDF5 Completed Successfully\n");
  //-------------------------------------

  //-Free loop write buffers and momentum matrix
  free(tmp_buf_uloc);
  free(buf_std_uloc);
  free(buf_gen_uloc);
  for(int mu = 0 ; mu < 4 ; mu++){
    free(tmp_buf_oneD[mu]);
    free(buf_std_oneD[mu]);
    free(buf_std_csvC[mu]);
    free(buf_gen_oneD[mu]);
    free(buf_gen_csvC[mu]);
  }
  free(tmp_buf_oneD);
  free(buf_std_oneD);
  free(buf_std_csvC);
  free(buf_gen_oneD);
  free(buf_gen_csvC);

  for(int ip=0;ip<Nmoms;ip++) free(mom[ip]);
  free(mom);
  //-------------------------------------

  return 0;
}

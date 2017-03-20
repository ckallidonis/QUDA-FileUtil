/* Christos Kallidonis                              */
/* March 2016                                       */
/* This code takes as input a .h5 file and extracts */
/* the desired loop data into an ASCII file         */
/* High-Momenta Form!                               */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>

#define M 16
#define DIR 4
#define QSQ_STR "Qsq"
#define NMOM_STR "Nmoms"
#define MOM_DSET "Momenta_list_xyz"

void usage(char exe[]){
  printf("%s: <.h5-file> <output_file> <loop (Scalar, dOp, Loops, LoopsCv, LpsDw, LpsDwCv)> <conf_traj> <T>\n",exe);
  exit(-1);
}

int main(int argc, char *argv[]){

  if(argc!=6) usage(argv[0]);

  char loop_type[6][256];
  bool loop_oneD[6];
  strcpy(loop_type[0],"Scalar");   loop_oneD[0] = false;   // std-ultra_local
  strcpy(loop_type[1],"dOp");      loop_oneD[1] = false;   // gen-ultra_local
  strcpy(loop_type[2],"Loops");    loop_oneD[2] = true;    // std-one_derivative
  strcpy(loop_type[3],"LoopsCv");  loop_oneD[3] = true;    // gen-one_derivative
  strcpy(loop_type[4],"LpsDw");    loop_oneD[4] = true;    // std-conserved current
  strcpy(loop_type[5],"LpsDwCv");  loop_oneD[5] = true;    // gen-conserved current

  char *h5_file, *out_file, *loop,*conf;
  asprintf(&h5_file ,"%s",argv[1]);
  asprintf(&out_file,"%s",argv[2]);
  asprintf(&loop    ,"%s",argv[3]);
  asprintf(&conf    ,"%s",argv[4]);

  int T = atoi(argv[5]);

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
  printf("conf traj: %s\n",conf);
  printf("T = %d\n",T);
  //-----------------------------------------

  //-Open the h5 file
  hid_t file_id = H5Fopen(h5_file, H5F_ACC_RDONLY, H5P_DEFAULT);


  //-Get the momenta-related attributes
  hid_t  Qattr = H5Aopen (file_id, QSQ_STR , H5P_DEFAULT);
  hid_t  Mattr = H5Aopen (file_id, NMOM_STR, H5P_DEFAULT);
  hid_t  Qattr_type = H5Aget_type(Qattr);
  hid_t  Mattr_type = H5Aget_type(Mattr);
  size_t Qattr_dim  = H5Tget_size(Qattr_type);
  size_t Mattr_dim  = H5Tget_size(Mattr_type);

  hid_t type_id  = H5Tcopy(H5T_C_S1);
  herr_t Qstatus = H5Tset_size (type_id, Qattr_dim);
  herr_t Mstatus = H5Tset_size (type_id, Mattr_dim);

  char *cQsq = (char*) malloc(Qattr_dim*sizeof(char));
  char *cNmoms = (char*) malloc(Mattr_dim*sizeof(char));

  Qstatus = H5Aread (Qattr, type_id, &(cQsq[0]));
  Mstatus = H5Aread (Mattr, type_id, &(cNmoms[0]));

  if (Mstatus<0 || Qstatus<0){
    fprintf (stderr, "Momenta attributes read failed!\n");
    Qstatus = H5Aclose(Qattr);
    Mstatus = H5Aclose(Mattr);
    Qstatus = H5Tclose(Qattr_type);
    Mstatus = H5Tclose(Mattr_type);
    Mstatus = H5Fclose(file_id);
    exit(-1);
  }

  int Qsq   = atoi(cQsq);
  int Nmoms = atoi(cNmoms);
  printf("Momenta attributes: Qsq = %d , Nmoms = %d\n",Qsq,Nmoms);

  Qstatus = H5Aclose(Qattr);
  Mstatus = H5Aclose(Mattr);
  Qstatus = H5Tclose(Qattr_type);
  Mstatus = H5Tclose(Mattr_type);
  //------------------------------------------------------------------

  //-Open the momenta dataset from the file and read the momenta
  int *moms = (int*) malloc(Nmoms*3*sizeof(int));

  hid_t Mdset_id = H5Dopen(file_id, MOM_DSET, H5P_DEFAULT);

  Mstatus = H5Dread(Mdset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, moms);

  if (Mstatus<0){
    fprintf (stderr, "Momenta Dataset read failed!\n");
    Mstatus = H5Dclose(Mdset_id);
    Mstatus = H5Fclose(file_id);
    free(moms);
    exit(-1);
  }
  //------------------------------------------------------------------


  //-Open the desired dataset and write it in the ASCII file
  hid_t group_id;
  hid_t dset_id;
  char *group_dir;

  double *loopBuf,*DloopBuf[4];

  FILE *out = fopen(out_file,"w");
  if(out==NULL){
    fprintf(stderr,"Cannot open %s for writing. Exiting.\n",out_file);
    exit(-1);
  }
  
  if(loop_oneD[dt]){
    for(int mu=0;mu<DIR;mu++){
      DloopBuf[mu] = (double*) malloc(T*Nmoms*M*2*sizeof(double));
      if(DloopBuf[mu]==NULL){
	fprintf(stderr,"Cannot allocate DloopBuf[%d]. Exiting.\n",mu);
	exit(-1);
      }

      asprintf(&group_dir,"/conf_%s/%s/dir_%02d",conf,loop,mu);
      group_id = H5Gopen(file_id, group_dir, H5P_DEFAULT);
      dset_id = H5Dopen(group_id, "loop", H5P_DEFAULT);
    
      herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, DloopBuf[mu]);

      if (status<0){
	fprintf (stderr, "Dataset read failed!\n");
	status = H5Gclose(group_id);
	status = H5Dclose(dset_id);
	status = H5Fclose(file_id);
	fclose(out);
	free(DloopBuf[mu]);
	free(moms);
	exit(-1);
      }

      H5Dclose(dset_id);
      H5Gclose(group_id);
    }//-mu

    for(int imom=0;imom<Nmoms;imom++)
      for(int mu=0;mu<DIR;mu++)
	for(int t=0;t<T;t++)
	  for(int m=0;m<M;m++)
	    fprintf(out,"%02d %02d %02d %+d %+d %+d %+16.15e %+16.15e\n",t,m,mu,moms[0+3*imom],moms[1+3*imom],moms[2+3*imom],DloopBuf[mu][0+2*m+2*M*imom+2*M*Nmoms*t],DloopBuf[mu][1+2*m+2*M*imom+2*M*Nmoms*t]); 

    for(int mu=0;mu<DIR;mu++) free(DloopBuf[mu]);

  }
  else{
    loopBuf = (double*) malloc(T*Nmoms*M*2*sizeof(double));
    if(loopBuf==NULL){
      fprintf(stderr,"Cannot allocate loopBuf. Exiting.\n");
      exit(-1);
    }

    asprintf(&group_dir,"/conf_%s/%s",conf,loop);
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
      free(moms);
      exit(-1);
    }

    for(int imom=0;imom<Nmoms;imom++)
      for(int t=0;t<T;t++)
	for(int m=0;m<M;m++)
	  fprintf(out,"%02d %02d %+d %+d %+d %+16.15e %+16.15e\n",t,m,moms[0+3*imom],moms[1+3*imom],moms[2+3*imom],loopBuf[0+2*m+2*M*imom+2*M*Nmoms*t],loopBuf[1+2*m+2*M*imom+2*M*Nmoms*t]); 

    H5Dclose(dset_id);
    H5Gclose(group_id);
    
    free(loopBuf);
  }

  H5Fclose(file_id);

  fclose(out);
  free(moms);

  printf("Extracting Loop completed successfully.\n");

  return 0;
}

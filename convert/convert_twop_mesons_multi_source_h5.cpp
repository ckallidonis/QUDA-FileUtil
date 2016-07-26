/* Christos Kallidonis                                    */
/* July 2016                                              */
/* This program reads meson two-point functions written   */
/* in ASCII format and converts them into HDF5 format     */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <typeinfo>

typedef float Float;
//typedef double Float; // Uncomment this line and comment the previous one for double precision !!

typedef struct{
  int Nmoms,Qsq,T,conf,src_pos[4],N_PART;
  char twop_type[10][256];
} Info;

void errorMsg(const char msg[]){
  fprintf(stderr,"%s",msg);
  exit(-1);
}
//=============================================================

void usage(char exe[]){
  printf("%s: <.h5-file prefix> <twop ASCII in_dir> <mom_list> <src_list> <Nmoms> <Qsq> <Nsrc> <T> <conf-trajectory>\n",exe);
  exit(-1);
}
//=============================================================

void writeTwop_HDF5(Float *twopMesons, char *h5_file, int **mom, Info info){

  hid_t DATATYPE_H5;
  if( typeid(Float) == typeid(float)  )  DATATYPE_H5 = H5T_NATIVE_FLOAT;
  if( typeid(Float) == typeid(double) )  DATATYPE_H5 = H5T_NATIVE_DOUBLE;

  int T  = info.T;
  int Nmoms = info.Nmoms;

  char fname[512];
  sprintf(fname,"%s.%04d_mesons_Qsq%d_SS.%02d.%02d.%02d.%02d.h5",h5_file,info.conf,info.Qsq,info.src_pos[0],info.src_pos[1],info.src_pos[2],info.src_pos[3]);
  printf("The two-point function HDF5 file is: %s\n",fname);

  Float *writeTwopBuf;

  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  char *group1_tag;
  asprintf(&group1_tag,"conf_%04d",info.conf);
  hid_t group1_id = H5Gcreate(file_id, group1_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  char *group2_tag;
  asprintf(&group2_tag,"sx%02dsy%02dsz%02dst%02d",info.src_pos[0],info.src_pos[1],info.src_pos[2],info.src_pos[3]);
  hid_t group2_id = H5Gcreate(group1_id, group2_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t group3_id;
  hid_t group4_id;

  hsize_t dims[2] = {T,2}; // Size of the dataspace

  for(int mes=0;mes<info.N_PART;mes++){
    char *group3_tag;
    asprintf(&group3_tag,"%s",info.twop_type[mes]);
    group3_id = H5Gcreate(group2_id, group3_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for(int imom=0;imom<Nmoms;imom++){
      char *group4_tag;
      asprintf(&group4_tag,"mom_xyz_%+d_%+d_%+d",mom[imom][0],mom[imom][1],mom[imom][2]);
      group4_id = H5Gcreate(group3_id, group4_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      hid_t filespace  = H5Screate_simple(2, dims, NULL);

      for(int ip=0;ip<2;ip++){
	char *dset_tag;
	asprintf(&dset_tag,"twop_meson_%d",ip+1);

	hid_t dataset_id = H5Dcreate(group4_id, dset_tag, DATATYPE_H5, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	writeTwopBuf = &(twopMesons[2*T*imom + 2*T*Nmoms*mes + 2*T*Nmoms*info.N_PART*ip]);

	herr_t status = H5Dwrite(dataset_id, DATATYPE_H5, H5S_ALL, filespace, H5P_DEFAULT, writeTwopBuf);
        if(status<0) errorMsg("File writing NOT SUCCESSFUL.\n");

	H5Dclose(dataset_id);
      }//-ip                                                                                                                                                                                                                                 
      H5Sclose(filespace);
      H5Gclose(group4_id);
    }//-imom                                                                                                                                                                                                                                 
    H5Gclose(group3_id);
  }//-mes                                                                                                                                                                                                                                    

  H5Gclose(group2_id);
  H5Gclose(group1_id);
  H5Fclose(file_id);

}
//=============================================================

int main(int argc, char *argv[]){

  if(argc!=10) usage(argv[0]);

  char twop_type[10][256];
  strcpy(twop_type[0],"pseudoscalar");
  strcpy(twop_type[1],"scalar");
  strcpy(twop_type[2],"g5g1");
  strcpy(twop_type[3],"g5g2");
  strcpy(twop_type[4],"g5g3");
  strcpy(twop_type[5],"g5g4");
  strcpy(twop_type[6],"g1");
  strcpy(twop_type[7],"g2");
  strcpy(twop_type[8],"g3");
  strcpy(twop_type[9],"g4");

  char *h5_file,*in_dir,*mom_list,*src_list;
  asprintf(&h5_file   ,"%s",argv[1]);
  asprintf(&in_dir    ,"%s",argv[2]);
  asprintf(&mom_list  ,"%s",argv[3]);
  asprintf(&src_list  ,"%s",argv[4]);
  int Nmoms  = atoi(argv[5]);
  int Qsq    = atoi(argv[6]);
  int Nsrc   = atoi(argv[7]);
  int T      = atoi(argv[8]);
  int conf   = atoi(argv[9]);

  int N_MESONS = 10;

  printf("Got the following input:\n");
  printf("h5_file prefix: %s\n",h5_file);
  printf("ASCII in_dir: %s\n",in_dir);
  printf("momenta list: %s\n",mom_list);
  printf("Nmoms  = %d\n",Nmoms);
  printf("Qsq    = %d\n",Qsq);
  printf("Nsrc  = %d\n",Nsrc);
  printf("src list: %s\n",src_list);
  printf("T      = %d\n",T);
  printf("conf   = %d\n",conf);
  //-------------------------------------

  //-Allocate the write buffer
  Float *Twop_HDF5 = NULL;

  if( (Twop_HDF5 = (Float*) malloc(2*T*Nmoms*2*N_MESONS*sizeof(Float)))==NULL ) errorMsg("Cannot allocate memory for Twop_HDF5.\n");
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

  //-Read the source positions
  int **src_pos;
  if( (src_pos = (int**) malloc(Nsrc*sizeof(int*)))==NULL ) errorMsg("Error in allocating src_pos.\n");
  for(int is=0; is<Nsrc; is++)
    if( (src_pos[is] = (int*) malloc(4*sizeof(int)))==NULL ) errorMsg("Error in allocating src_pos inside is-loop.\n");

  FILE *p_src;
  if( (p_src=fopen(src_list,"r"))==NULL ) errorMsg("Cannot open source list for reading. Exiting.\n");

  for(int is=0;is<Nsrc;is++) fscanf(p_src,"%d %d %d %d\n",&src_pos[is][0],&src_pos[is][1],&src_pos[is][2],&src_pos[is][3]);
  fclose(p_src);
  //------------------------------------- 

  //-Read the two-point functions from the ASCII files
  int idum;
  char *ASCII_file;
  FILE *p_twop;

  for(int is=0;is<Nsrc;is++){
    asprintf(&ASCII_file,"%s/twop.%04d.mesons.SS.%02d.%02d.%02d.%02d.dat",in_dir,conf,src_pos[is][0],src_pos[is][1],src_pos[is][2],src_pos[is][3]);
    printf("Reading twop file %s\n",ASCII_file);
    if( (p_twop = fopen(ASCII_file,"r"))==NULL ) errorMsg("Cannot open twop file for reading. Exiting.\n");

    for(int mes=0;mes<N_MESONS;mes++){
      for(int it=0;it<T;it++){
	for(int imom=0;imom<Nmoms;imom++){
	  fscanf(p_twop,"%d %d %d %d %d %f %f %f %f\n",&idum,&idum,&idum,&idum,&idum,
		 &(Twop_HDF5[0 + 2*it + 2*T*imom + 2*T*Nmoms*mes + 2*T*Nmoms*N_MESONS*0]), &(Twop_HDF5[1 + 2*it + 2*T*imom + 2*T*Nmoms*mes + 2*T*Nmoms*N_MESONS*0]),
		 &(Twop_HDF5[0 + 2*it + 2*T*imom + 2*T*Nmoms*mes + 2*T*Nmoms*N_MESONS*1]), &(Twop_HDF5[1 + 2*it + 2*T*imom + 2*T*Nmoms*mes + 2*T*Nmoms*N_MESONS*1]));
	}}
    }
    fclose(p_twop);

    printf("Reading the two-point function in ASCII format for is = %d completed.\n",is);
    //------------------------------------- 

    //-Write the two-point function in HDF5 format
    Info twopInfo;

    twopInfo.Nmoms = Nmoms; twopInfo.Qsq = Qsq; twopInfo.T = T; twopInfo.conf = conf; twopInfo.N_PART = N_MESONS;
    for(int i=0;i<4;i++) twopInfo.src_pos[i] = src_pos[is][i];
    for(int i=0;i<N_MESONS;i++) strcpy(twopInfo.twop_type[i],twop_type[i]);

    writeTwop_HDF5(Twop_HDF5, h5_file, mom, twopInfo);
    printf("Two-point function for is = %d written in HDF5 format successfully\n",is);

  }//-is
  //------------------------------------- 


  for(int ip=0;ip<Nmoms;ip++) free(mom[ip]);
  free(mom);

  for(int is=0;is<Nsrc;is++) free(src_pos[is]);
  free(src_pos);

  free(Twop_HDF5);

  return 0;
}

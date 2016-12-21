/* Christos Kallidonis                                */
/* April 2016                                         */
/* This program reads three-point functions written   */
/* in ASCII format and converts them into HDF5 format */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <typeinfo>

typedef float Float;
//typedef double Float; // Uncomment this line and comment the previous one for double precision !!


typedef struct{
  int Nmoms,Qsq,T,conf,src_pos[4],Ntsink,tsink[10];
  char thrp_type[3][256];
  char thrp_proj_type[5][256];
  char thrp_part[2][256];
} Info;

enum THRP_TYPE{THRP_LOCAL,THRP_NOETHER,THRP_ONED};

void errorMsg(const char msg[]){
  fprintf(stderr,"%s",msg);
  exit(-1);
}
//=============================================================

void usage(char exe[]){
  printf("%s: <.h5-file prefix> <threep ASCII in_dir> <mom_list> <tsink_list> <Nmoms> <Qsq> <sx> <sy> <sz> <st> <T> <conf-trajectory>\n",exe);
  exit(-1);
}
//=============================================================

//template<typename Float>
void writeThrp_HDF5(Float *Thrp_local_HDF5, Float *Thrp_noether_HDF5, Float **Thrp_oneD_HDF5, char *h5_file, int **mom, Info info){

  hid_t DATATYPE_H5;
  if( typeid(Float) == typeid(float) ) DATATYPE_H5 = H5T_NATIVE_FLOAT;
  if( typeid(Float) == typeid(double)) DATATYPE_H5 = H5T_NATIVE_DOUBLE;

  char fname[512];
  sprintf(fname,"%s.%04d_neutron_Qsq%d_SS.%02d.%02d.%02d.%02d.h5",h5_file,info.conf,info.Qsq,info.src_pos[0],info.src_pos[1],info.src_pos[2],info.src_pos[3]);
  printf("The three-point function HDF5 file is: %s\n",fname);

  Float *writeThrpBuf;

  int Nsink = info.Ntsink;
  int Nmoms = info.Nmoms;
  int T  = info.T;
  int Mel;

  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  char *group1_tag;
  asprintf(&group1_tag,"conf_%04d",info.conf);
  hid_t group1_id = H5Gcreate(file_id, group1_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  char *group2_tag;
  asprintf(&group2_tag,"sx%02dsy%02dsz%02dst%02d",info.src_pos[0],info.src_pos[1],info.src_pos[2],info.src_pos[3]);
  hid_t group2_id = H5Gcreate(group1_id, group2_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t group3_id;
  hid_t group4_id;
  hid_t group5_id;
  hid_t group6_id;
  hid_t group7_id;
  hid_t group8_id;

  hsize_t dims[3];

  for(int its=0;its<Nsink;its++){
    int tsink = info.tsink[its];
    char *group3_tag;
    asprintf(&group3_tag,"tsink_%02d",tsink);
    group3_id = H5Gcreate(group2_id, group3_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    char *group4_tag;
    asprintf(&group4_tag,"proj_G4");
    group4_id = H5Gcreate(group3_id, group4_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for(int part=0;part<2;part++){
      char *group5_tag;
      asprintf(&group5_tag,"%s", (part==0) ? "up" : "down");
      group5_id = H5Gcreate(group4_id, group5_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      for(int thrp_int=0;thrp_int<3;thrp_int++){
	THRP_TYPE type = (THRP_TYPE) thrp_int;

	char *group6_tag;
	asprintf(&group6_tag,"%s", info.thrp_type[thrp_int]);
	group6_id = H5Gcreate(group5_id, group6_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//-Determine the global dimensions                                                                                                                                                                                                 
	if(type==THRP_LOCAL || type==THRP_ONED) Mel = 16;
	else if (type==THRP_NOETHER) Mel = 4;
	else{
	  printf("writeThrp_HDF5: Undefined three-point function type.\n");
	  exit(-1);
	}
	dims[0] = tsink+1;
	dims[1] = Mel;
	dims[2] = 2;

	for(int imom=0;imom<Nmoms;imom++){
	  char *group7_tag;
	  asprintf(&group7_tag,"mom_xyz_%+d_%+d_%+d",mom[imom][0],mom[imom][1],mom[imom][2]);
	  group7_id = H5Gcreate(group6_id, group7_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	  if(type==THRP_ONED){
	    for(int mu=0;mu<4;mu++){
	      char *group8_tag;
	      asprintf(&group8_tag,"dir_%02d",mu);
	      group8_id = H5Gcreate(group7_id, group8_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	      hid_t filespace  = H5Screate_simple(3, dims, NULL);
	      hid_t dataset_id = H5Dcreate(group8_id, "threep", DATATYPE_H5, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	      writeThrpBuf = &(Thrp_oneD_HDF5[mu][2*Mel*T*imom + 2*Mel*T*Nmoms*part + 2*Mel*T*Nmoms*2*its]);

	      herr_t status = H5Dwrite(dataset_id, DATATYPE_H5, H5S_ALL, filespace, H5P_DEFAULT, writeThrpBuf);
	      if(status<0) errorMsg("File writing NOT SUCCESSFUL.\n");

	      H5Dclose(dataset_id);
	      H5Sclose(filespace);
	      H5Gclose(group8_id);
	    }//-mu
	  }//-if
	  else{
	    Float *thrpBuf;
	    if(type==THRP_LOCAL) thrpBuf = Thrp_local_HDF5;
	    else if(type==THRP_NOETHER) thrpBuf = Thrp_noether_HDF5;

	    hid_t filespace  = H5Screate_simple(3, dims, NULL);
	    hid_t dataset_id = H5Dcreate(group7_id, "threep", DATATYPE_H5, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	    writeThrpBuf = &(thrpBuf[2*Mel*T*imom + 2*Mel*T*Nmoms*part + 2*Mel*T*Nmoms*2*its]);

	    herr_t status = H5Dwrite(dataset_id, DATATYPE_H5, H5S_ALL, filespace, H5P_DEFAULT, writeThrpBuf);
	    if(status<0) errorMsg("File writing NOT SUCCESSFUL.\n");

	    H5Dclose(dataset_id);
	    H5Sclose(filespace);
	  }//-else                                                                                                                                                                                                                         
	  H5Gclose(group7_id);
	}//-imom                                                                                                                                                                                                                           
	H5Gclose(group6_id);
      }//-thrp_int                                                                                                                                                                                                                         
      H5Gclose(group5_id);
    }//-part                                                                                                                                                                                                                               
    H5Gclose(group4_id);
    H5Gclose(group3_id);
  }//-its                                                                                                                                                                                                                                    

  H5Gclose(group2_id);
  H5Gclose(group1_id);
  H5Fclose(file_id);
  
}

//=============================================================

int main(int argc, char *argv[]){

  if(argc!=13) usage(argv[0]);

  char thrp_type[3][256];
  strcpy(thrp_type[0],"ultra_local");
  strcpy(thrp_type[1],"noether");
  strcpy(thrp_type[2],"oneD");

  char thrp_proj_type[5][256];
  strcpy(thrp_proj_type[0],"G4");
  strcpy(thrp_proj_type[1],"G5G123");
  strcpy(thrp_proj_type[2],"G5G1");
  strcpy(thrp_proj_type[3],"G5G2");
  strcpy(thrp_proj_type[4],"G5G3");

  char thrp_part[2][256];
  strcpy(thrp_part[0],"up");
  strcpy(thrp_part[1],"down");

  char *h5_file,*in_dir,*mom_list,*tsink_list;
  int src_pos[4];
  asprintf(&h5_file   ,"%s",argv[1]);
  asprintf(&in_dir    ,"%s",argv[2]);
  asprintf(&mom_list  ,"%s",argv[3]);
  asprintf(&tsink_list,"%s",argv[4]);
  int Nmoms  = atoi(argv[5]);
  int Qsq    = atoi(argv[6]);
  src_pos[0] = atoi(argv[7]);
  src_pos[1] = atoi(argv[8]);
  src_pos[2] = atoi(argv[9]);
  src_pos[3] = atoi(argv[10]);
  int T      = atoi(argv[11]);
  int conf   = atoi(argv[12]);

  FILE *t_ptr = fopen(tsink_list,"r");
  if(t_ptr==NULL) errorMsg("Cannot open tsink_list for reading\n");
  int Ntsink;
  fscanf(t_ptr,"%d\n",&Ntsink);
  int *tsink = (int*) malloc(Ntsink*sizeof(int));
  for(int i=0;i<Ntsink;i++) fscanf(t_ptr,"%d\n",&(tsink[i]));
  fclose(t_ptr);

  printf("Got the following input:\n");
  printf("h5_file prefix: %s\n",h5_file);
  printf("ASCII in_dir: %s\n",in_dir);
  printf("momenta list: %s\n",mom_list);
  printf("Ntsink = %d\n",Ntsink);
  for(int i=0;i<Ntsink;i++) printf("  tsink[%d] = %d\n",i,tsink[i]);
  printf("Nmoms  = %d\n",Nmoms);
  printf("Qsq    = %d\n",Qsq);
  printf("source position (x,y,z,t) = (%02d,%02d,%02d,%02d)\n",src_pos[0],src_pos[1],src_pos[2],src_pos[3]);
  printf("T      = %d\n",T);
  printf("conf   = %d\n",conf);
  //-------------------------------------

  //-Allocate the write buffers
  Float *Thrp_local_HDF5 = NULL;
  Float *Thrp_noether_HDF5 = NULL;
  Float **Thrp_oneD_HDF5 = NULL;

  if( (Thrp_local_HDF5   = (Float*) malloc(2*16*T*Nmoms*2*Ntsink*sizeof(Float)))==NULL ) errorMsg("Cannot allocate memory for Thrp_local_HDF5.\n");
  if( (Thrp_noether_HDF5 = (Float*) malloc(2* 4*T*Nmoms*2*Ntsink*sizeof(Float)))==NULL ) errorMsg("Cannot allocate memory for Thrp_noether_HDF5.\n");

  if( (Thrp_oneD_HDF5 = (Float**) malloc(4*sizeof(Float*))) == NULL ) errorMsg("Cannot allocate memory for Thrp_oneD_HDF5.\n");
  for(int mu=0;mu<4;mu++) if( (Thrp_oneD_HDF5[mu] = (Float*) malloc(2*16*T*Nmoms*2*Ntsink*sizeof(Float)))==NULL ) errorMsg("Cannot allocate memory for Thrp_oned_HDF5 inside mu-loop.\n");
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


  //-Read the three-point functions from the ASCII files
  int idum,Mel;
  char *ASCII_file;
  FILE *p_thrp;

  for(int its=0;its<Ntsink;its++){
    for(int thrp_int=0;thrp_int<3;thrp_int++){
      THRP_TYPE type = (THRP_TYPE) thrp_int;    
      if(type==THRP_LOCAL || type==THRP_ONED) Mel = 16;
      else if (type==THRP_NOETHER) Mel = 4;

      for(int uORd=0;uORd<2;uORd++){ // up and down
	asprintf(&ASCII_file,"%s/threep.%04d_tsink%d.neutron.%s.%s.SS.%02d.%02d.%02d.%02d.dat",in_dir,conf,tsink[its],thrp_part[uORd],thrp_type[thrp_int],src_pos[0],src_pos[1],src_pos[2],src_pos[3]);
        printf("Reading threep file %s\n",ASCII_file);
        if( (p_thrp = fopen(ASCII_file,"r"))==NULL ) errorMsg("Cannot open thrp file for reading. Exiting.\n");

	if(type==THRP_LOCAL){
	  for(int im=0;im<Mel;im++){
	    for(int it=0;it<T;it++){
	      for(int imom=0;imom<Nmoms;imom++){
		fscanf(p_thrp,"%d %d %d %d %d %f %f\n",&idum,&idum,&idum,&idum,&idum,
		       &(Thrp_local_HDF5[0 + 2*im + 2*Mel*it + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its]),
		       &(Thrp_local_HDF5[1 + 2*im + 2*Mel*it + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its]));
	      }}
	  }
	}//-if THRP_LOCAL
	else if(type==THRP_NOETHER){
	  for(int im=0;im<Mel;im++){
	    for(int it=0;it<T;it++){
	      for(int imom=0;imom<Nmoms;imom++){
		fscanf(p_thrp,"%d %d %d %d %d %f %f\n",&idum,&idum,&idum,&idum,&idum,
		       &(Thrp_noether_HDF5[0 + 2*im + 2*Mel*it + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its]),
		       &(Thrp_noether_HDF5[1 + 2*im + 2*Mel*it + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its]));
	      }}
	  }
	}//-if THRP_NOETHER
	else if(type==THRP_ONED){
	  for(int im=0;im<Mel;im++){
	    for(int mu=0;mu<4;mu++){
	      for(int it=0;it<T;it++){
		for(int imom=0;imom<Nmoms;imom++){
		  fscanf(p_thrp,"%d %d %d %d %d %d %f %f\n",&idum,&idum,&idum,&idum,&idum,&idum,
			 &(Thrp_oneD_HDF5[mu][0 + 2*im + 2*Mel*it + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its]),
			 &(Thrp_oneD_HDF5[mu][1 + 2*im + 2*Mel*it + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its]));
		}}}
	  }
	}//-if THRP_ONED
	
	fclose(p_thrp);
      }//-uORd
    }//-thrp_int
  }//-its
  printf("Reading the three-point function in ASCII format completed.\n");
  //------------------------------------- 


  //-Write the three-point function in HDF5 format
  Info thrpInfo;

  thrpInfo.Nmoms = Nmoms; thrpInfo.Qsq = Qsq; thrpInfo.T = T; thrpInfo.conf = conf; thrpInfo.Ntsink = Ntsink;
  for(int i=0;i<4;i++) thrpInfo.src_pos[i] = src_pos[i];
  for(int i=0;i<Ntsink;i++) thrpInfo.tsink[i] = tsink[i];
  for(int i=0;i<3;i++) strcpy(thrpInfo.thrp_type[i],thrp_type[i]);
  for(int i=0;i<5;i++) strcpy(thrpInfo.thrp_proj_type[i],thrp_proj_type[i]);
  for(int i=0;i<2;i++) strcpy(thrpInfo.thrp_part[i],thrp_part[i]);

  writeThrp_HDF5(Thrp_local_HDF5, Thrp_noether_HDF5, Thrp_oneD_HDF5, h5_file, mom, thrpInfo);
  printf("Three-point function writte in HDF5 format successfully\n");
  //------------------------------------- 


  for(int ip=0;ip<Nmoms;ip++) free(mom[ip]);
  free(mom);

  free(tsink);

  free(Thrp_local_HDF5);
  free(Thrp_noether_HDF5);
  for(int mu=0;mu<4;mu++) free(Thrp_oneD_HDF5[mu]);
  free(Thrp_oneD_HDF5);

  return 0;
}

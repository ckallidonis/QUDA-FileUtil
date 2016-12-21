/* Christos Kallidonis                                       */
/* July 2016                                                 */
/* This code takes as input a .h5 file with a three-point    */
/* function in position-space, it performs the Fourier       */
/* transform and writes the result in an ASCII or HDF5 file  */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <typeinfo>
#include <fftw3.h>

#define MAX_MOM 5000
#define MAX_PROJ 5
#define MAX_SINK 10
#define Ntype 3

typedef float Float;
//typedef double Float; // Uncomment this line and comment the previous one for double precision !!

typedef struct{
  int Nmoms,Qsq,T,L,conf,src_pos[4],Ntsink,tsink[MAX_SINK];
  int Nproj[MAX_SINK], proj_list[MAX_SINK][MAX_PROJ];
  char thrp_type[3][256];
  char thrp_proj_type[5][256];
  char thrp_part[2][256];
} Info;


enum THRP_TYPE{THRP_LOCAL,THRP_NOETHER,THRP_ONED};

void usage(char exe[]){
  printf("%s:\n<.h5-file>\n<output_prefix>\n<proj_list>\n<tsink_list>\n<Ntsink>\n<conf_traj>\n<src_x>\n<src_y>\n<src_z>\n<src_t>\n<L>\n<T>\n<Qsq>\n<output format: ASCII/HDF5>\n",exe);
  exit(-1);
}
//=================================================================================================

int createMomenta(int **mom, int L, int Qsq){

  bool nxOK,nyOK,nzOK;

  int totMom = 0;
  for(int iQ = 0 ; iQ <= Qsq ; iQ++){
    for(int nx = iQ ; nx >= -iQ ; nx--)
      for(int ny = iQ ; ny >= -iQ ; ny--)
        for(int nz = iQ ; nz >= -iQ ; nz--){
          nxOK = ( (nx>=(-L/2)) && (nx<L/2) ) ? true : false;
          nyOK = ( (ny>=(-L/2)) && (ny<L/2) ) ? true : false;
          nzOK = ( (nz>=(-L/2)) && (nz<L/2) ) ? true : false;

          if( ((nx*nx + ny*ny + nz*nz) == iQ) && nxOK && nyOK && nzOK ){
            if(totMom==MAX_MOM){
              fprintf(stderr,"Total number of momenta exceeded! Exiting.\n");
              exit(-1);
            }
            mom[totMom][0] = nx;
            mom[totMom][1] = ny;
            mom[totMom][2] = nz;
            totMom++;
          }
        }
  }

  return totMom;
}
//=================================================================================================

long int createALLMomenta(int **momALL, int L, int Qsq){

  long int momIdx = 0;
  for(int pz = 0; pz < L; pz++)
    for(int py = 0; py < L; py++)
      for(int px = 0; px < L; px++){
        if(px < L/2)
          momALL[momIdx][0]   = px;
        else
          momALL[momIdx][0]   = px - L;

        if(py < L/2)
          momALL[momIdx][1]   = py;
	else
          momALL[momIdx][1]   = py - L;

        if(pz < L/2)
          momALL[momIdx][2]   = pz;
        else
          momALL[momIdx][2]   = pz - L;

        momIdx++;
      }

  return momIdx;
}
//=================================================================================================

void writeThrp_HDF5(Float *Thrp_local_HDF5, Float *Thrp_noether_HDF5, Float **Thrp_oneD_HDF5, char *fname, int **mom, Info info){

  hid_t DATATYPE_H5;
  if( typeid(Float) == typeid(float) ) DATATYPE_H5 = H5T_NATIVE_FLOAT;
  if( typeid(Float) == typeid(double)) DATATYPE_H5 = H5T_NATIVE_DOUBLE;

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

    for(int ipr=0;ipr<info.Nproj[its];ipr++){
      char *group4_tag;
      asprintf(&group4_tag,"proj_%s",info.thrp_proj_type[info.proj_list[its][ipr]]);
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
		
		writeThrpBuf = &(Thrp_oneD_HDF5[mu][2*Mel*T*imom + 2*Mel*T*Nmoms*part + 2*Mel*T*Nmoms*2*its + 2*Mel*T*Nmoms*2*Nsink*ipr]);
		
		herr_t status = H5Dwrite(dataset_id, DATATYPE_H5, H5S_ALL, filespace, H5P_DEFAULT, writeThrpBuf);
		if(status<0){
		  printf("File writing NOT SUCCESSFUL.\n");
		  exit(-1);
		}

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
	      
	      writeThrpBuf = &(thrpBuf[2*Mel*T*imom + 2*Mel*T*Nmoms*part + 2*Mel*T*Nmoms*2*its + 2*Mel*T*Nmoms*2*Nsink*ipr]);
	      
	      herr_t status = H5Dwrite(dataset_id, DATATYPE_H5, H5S_ALL, filespace, H5P_DEFAULT, writeThrpBuf);
	      if(status<0){
		printf("File writing NOT SUCCESSFUL.\n");
		exit(-1);
	      }

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
    }//-proj
    H5Gclose(group3_id);
  }//-its

  H5Gclose(group2_id);
  H5Gclose(group1_id);
  H5Fclose(file_id);

}
//=================================================================================================



void copyThrpToWriteBuf(Float *thrpMom, Float *thrpFFT, int **momALL, int **mom, int its, int ipr, int mu, Info info, THRP_TYPE type){

  int T = info.T;
  int L = info.L;
  long int SpV = (long int)L*L*L;
  int x_src = info.src_pos[0];
  int y_src = info.src_pos[1];
  int z_src = info.src_pos[2];
  int t_src = info.src_pos[3];
  int Qsq = info.Qsq;
  int Nmoms = info.Nmoms;
  int Ntsink = info.Ntsink;
  int Mel;
  if(type==THRP_LOCAL || type==THRP_ONED) Mel = 16;
  else if(type==THRP_NOETHER) Mel = 4;
  else{
    printf("Undefined THRP_TYPE passed to copyThrpToWriteBuf.\n");
    exit(-1);
  }

  Float two_pi = 4.0*asin(1.0);

  for(int ip=0;ip<SpV;ip++){
    int px = momALL[ip][0];
    int py = momALL[ip][1];
    int pz = momALL[ip][2];
    
    if( (px*px + py*py + pz*pz) <= Qsq ){
      int imom = -1;
      for(int im=0;im<Nmoms;im++)                                          // With this we ensure that the momenta
	if(px==mom[im][0] && py==mom[im][1] && pz==mom[im][2]) imom = im;  // have the standard ordering
      
      if(imom==-1){
        printf("Check your momenta arrays. Exiting.\n");
        exit(-1);
      }

      Float expn = two_pi*(px*x_src + py*y_src + pz*z_src ) / (Float)L;
      Float phase[2];
      phase[0] =  cos(expn);
      phase[1] = -sin(expn);

      for(int uORd=0;uORd<2;uORd++){
	for(int t=0;t<T;t++){
	  int ts = (t + t_src)%T;
	  for(int im=0;im<Mel;im++){
	    unsigned long long Mompos_Re = (unsigned long long)(0ULL + 2*im + 2*Mel*t  + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its + 2*Mel*T*Nmoms*2*Ntsink*ipr);
	    unsigned long long Mompos_Im = (unsigned long long)(1ULL + 2*im + 2*Mel*t  + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its + 2*Mel*T*Nmoms*2*Ntsink*ipr);
	    unsigned long long FFTpos_Re = (unsigned long long)(0ULL + 2*ip + 2*SpV*ts + 2*SpV*T*uORd + 2*SpV*T*2*im       + 2*SpV*T*2*Mel*mu);
	    unsigned long long FFTpos_Im = (unsigned long long)(1ULL + 2*ip + 2*SpV*ts + 2*SpV*T*uORd + 2*SpV*T*2*im       + 2*SpV*T*2*Mel*mu);

	    thrpMom[ Mompos_Re ] = thrpFFT[ FFTpos_Re ]*phase[0] - thrpFFT[ FFTpos_Im ]*phase[1];
	    thrpMom[ Mompos_Im ] = thrpFFT[ FFTpos_Re ]*phase[1] + thrpFFT[ FFTpos_Im ]*phase[0];
	  }}}
      imom++;
    }//-if
  }

}
//=================================================================================================

int main(int argc, char *argv[]){

  if(argc!=15) usage(argv[0]);

  char thrp_type[3][256];

  strcpy(thrp_type[0],"ultra_local");
  strcpy(thrp_type[1],"noether");
  strcpy(thrp_type[2],"oneD");

  char thrp_part[2][256];
  strcpy(thrp_part[0],"up");
  strcpy(thrp_part[1],"down");

  char thrp_proj_type[5][256];
  strcpy(thrp_proj_type[0],"G4");
  strcpy(thrp_proj_type[1],"G5G123");
  strcpy(thrp_proj_type[2],"G5G1");
  strcpy(thrp_proj_type[3],"G5G2");
  strcpy(thrp_proj_type[4],"G5G3");

  char *h5_file, *outpre, *tsink_list,*proj_list_f, *outform;
  int src[4];
  asprintf(&h5_file    ,"%s",argv[1]);
  asprintf(&outpre     ,"%s",argv[2]);
  asprintf(&proj_list_f,"%s",argv[3]);
  asprintf(&tsink_list ,"%s",argv[4]);
  int Ntsink = atoi(argv[5]);
  int conf   = atoi(argv[6]);
  src[0] = atoi(argv[7]);
  src[1] = atoi(argv[8]);
  src[2] = atoi(argv[9]);
  src[3] = atoi(argv[10]);
  int L     = atoi(argv[11]);
  int T     = atoi(argv[12]);
  int Qsq   = atoi(argv[13]);

  asprintf(&outform,"%s",argv[14]);
  bool h5out;
  if( strcmp(outform,"HDF5")==0 ) h5out=true;
  else if( strcmp(outform,"ASCII")==0 ) h5out=false;
  else{
    printf("output format must be either ASCII or HDF5. Exiting\n");
    exit(-1);
  }

  //-Read the tsinks
  FILE *ptr_tsink;
  int tsink[Ntsink];
  if( (ptr_tsink = fopen(tsink_list,"r")) == NULL ){
    fprintf(stderr,"Error opening file %s for sink-source separations\n",tsink_list);
    exit(-1);
  }
  for(int it=0;it<Ntsink;it++) fscanf(ptr_tsink,"%d\n",&(tsink[it]));
  fclose(ptr_tsink);

  //-Read the projectors
  FILE *proj_ptr;
  char *proj_file;
  int Nproj[Ntsink], proj_list[Ntsink][MAX_PROJ];
  for(int it=0;it<Ntsink;it++){
    asprintf(&proj_file,"%s_tsink%d.txt",proj_list_f,tsink[it]);
    if( (proj_ptr = fopen(proj_file,"r")) == NULL ){
      fprintf(stderr,"Cannot open projector file %s for reading.\n Hint: Make sure that 1: it ends as _tsink%d.txt and 2: the input passed is truncated up to this string.\n",proj_file,tsink[it]);
      exit(-1);
    }
    fscanf(proj_ptr,"%d",&(Nproj[it]));
    for(int p=0;p<Nproj[it];p++) fscanf(proj_ptr,"%d\n",&(proj_list[it][p]));
    fclose(proj_ptr);
  }

  printf("Got the following input:\n");
  printf("h5_file: %s\n",h5_file);
  printf("output prefix: %s\n",outpre);
  printf("Got the following %d source-sink separations and projectors:\n",Ntsink);
  int NprojMax = 0;
  for(int its=0;its<Ntsink;its++){
    if(Nproj[its] >= NprojMax) NprojMax = Nproj[its];
    printf(" sink-source = %d:\n",tsink[its]);
    for(int p=0;p<Nproj[its];p++) printf("  %s\n",thrp_proj_type[proj_list[its][p]]);
  }
  printf("conf traj: %04d\n",conf);
  printf("src [x,y,z,t] = [%02d,%02d,%02d,%02d]\n",src[0],src[1],src[2],src[3]);
  printf("L     = %02d\n",L);
  printf("T     = %02d\n",T);
  printf("Qsq   = %02d\n",Qsq);
  printf("Output format is %s\n", h5out ? "HDF5" : "ASCII");
  //-----------------------------------------

  //-Define useful stuff
  long int SpV = (long int)L*L*L;
  long int V   = (long int)SpV*T;

  Info thrpInfo;
  thrpInfo.Qsq = Qsq; thrpInfo.T = T; thrpInfo.L = L; thrpInfo.Ntsink = Ntsink;
  thrpInfo.conf = conf;
  for(int i=0;i<4;i++) thrpInfo.src_pos[i] = src[i];
  for(int i=0;i<3;i++) strcpy(thrpInfo.thrp_type[i],thrp_type[i]);
  for(int i=0;i<5;i++) strcpy(thrpInfo.thrp_proj_type[i],thrp_proj_type[i]);
  for(int i=0;i<2;i++) strcpy(thrpInfo.thrp_part[i],thrp_part[i]);
  for(int i=0;i<Ntsink;i++){
    thrpInfo.tsink[i] = tsink[i];
    thrpInfo.Nproj[i] = Nproj[i];
    for(int p=0;p<Nproj[i];p++) thrpInfo.proj_list[i][p] = proj_list[i][p];
  }
  //-----------------------------------------

  //-Create the momenta arrays
  int **mom, **momALL;

  mom    = (int**) malloc(MAX_MOM*sizeof(int*));
  momALL = (int**) malloc(SpV*sizeof(int*));
  if(mom==NULL){
    fprintf(stderr,"Cannot allocate mom, top-level\n");
    exit(-1);
  }
  if(momALL==NULL){
    fprintf(stderr,"Cannot allocate momALL, top-level\n");
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
  for(int i=0;i<SpV;i++){
    momALL[i] = (int*) malloc(3*sizeof(int));
    if(momALL[i] == NULL){
      fprintf(stderr,"Cannot allocate momALL[%d]\n",i);
      exit(-1);
    }
    for(int j=0;j<3;j++) momALL[i][j] = 0;
  }

  int Nmoms = createMomenta(mom, L, Qsq);
  printf("Required Momenta array created, Nmoms = %d\n",Nmoms);
  thrpInfo.Nmoms = Nmoms;

  int Qsq_ALL = 3*L*L/4;
  long int Nmoms_ALL = createALLMomenta(momALL, L, Qsq_ALL);
  printf("All Momenta array created, Nmoms_ALL = %ld, Sp.Vol = %ld\n",Nmoms_ALL,SpV);
  //----------------------------------------------

  //-Open the h5 file
  hid_t file_id = H5Fopen(h5_file, H5F_ACC_RDONLY, H5P_DEFAULT);

  hid_t DATATYPE_H5;
  if( typeid(Float) == typeid(float)  )  DATATYPE_H5 = H5T_NATIVE_FLOAT;
  if( typeid(Float) == typeid(double) )  DATATYPE_H5 = H5T_NATIVE_DOUBLE;

  //-Read the dataset
  hid_t group_id;
  hid_t dset_id;
  char *group_dir;
  char *dset_name;
  asprintf(&dset_name,"threep");

  unsigned long long buflen_uloc = (unsigned long long)(2)*V*2*16  *Ntsink*NprojMax;
  unsigned long long buflen_noet = (unsigned long long)(2)*V*2* 4  *Ntsink*NprojMax;
  unsigned long long buflen_oneD = (unsigned long long)(2)*V*2*16*4*Ntsink*NprojMax;
  printf("Buffer lengths:\n");
  printf(" u-local = %lld\n",buflen_uloc);
  printf(" noether = %lld\n",buflen_noet);
  printf(" oneD    = %lld\n",buflen_oneD);

  Float *thrpBuf_uloc = (Float*) calloc(buflen_uloc,sizeof(Float));
  Float *thrpBuf_noet = (Float*) calloc(buflen_noet,sizeof(Float));
  Float *thrpBuf_oneD = (Float*) calloc(buflen_oneD,sizeof(Float));
  if( thrpBuf_uloc == NULL ){
    fprintf(stderr,"Cannot allocate thrpBuf_uloc. Exiting\n");
    exit(-1);
  }
  if( thrpBuf_noet == NULL ){
    fprintf(stderr,"Cannot allocate thrpBuf_noet. Exiting\n");
    exit(-1);
  }
  if( thrpBuf_oneD == NULL ){
    fprintf(stderr,"Cannot allocate thrpBuf_oneD. Exiting\n");
    exit(-1);
  }

  printf("Reading Dataset...\n");
  for(int its=0;its<Ntsink;its++){
    for(int ipr=0;ipr<Nproj[its];ipr++){
      unsigned long long posu = (unsigned long long)(2)*V*2*16  *its + (unsigned long long)(2)*V*2*16  *Ntsink*ipr;
      unsigned long long posn = (unsigned long long)(2)*V*2* 4  *its + (unsigned long long)(2)*V*2* 4  *Ntsink*ipr;
      unsigned long long poso = (unsigned long long)(2)*V*2*16*4*its + (unsigned long long)(2)*V*2*16*4*Ntsink*ipr;

      for(int itype=0;itype<Ntype;itype++){
	THRP_TYPE type = (THRP_TYPE) itype;

	asprintf(&group_dir,"/conf_%04d/sx%02dsy%02dsz%02dst%02d/tsink_%d/proj_%s/%s",conf,src[0],src[1],src[2],src[3],tsink[its],thrp_proj_type[proj_list[its][ipr]],thrp_type[itype]);

	group_id = H5Gopen(file_id, group_dir, H5P_DEFAULT);
	dset_id = H5Dopen(group_id, dset_name, H5P_DEFAULT);
  
	herr_t status;
	if(type==THRP_LOCAL)   status = H5Dread(dset_id, DATATYPE_H5, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(thrpBuf_uloc[posu]));
	if(type==THRP_NOETHER) status = H5Dread(dset_id, DATATYPE_H5, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(thrpBuf_noet[posn]));
	if(type==THRP_ONED)    status = H5Dread(dset_id, DATATYPE_H5, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(thrpBuf_oneD[poso]));
  
	if (status<0){
	  fprintf (stderr, "Dataset read failed!\n");
	  status = H5Gclose(group_id);
	  status = H5Dclose(dset_id);
	  status = H5Fclose(file_id);
	  free(thrpBuf_uloc);
	  free(thrpBuf_noet);
	  free(thrpBuf_oneD);
	  exit(-1);
	}
	
	H5Dclose(dset_id);
	H5Gclose(group_id);
      }
    }
  }
  printf("Dataset reading successfull\n");
  //----------------------------------------------

  //-Perform the FT
  int rank = 3;
  int nRank[3] = {L,L,L};
  int istride = 1; // array is contigious
  int ostride = 1; // in memory

  int *inembed = nRank, *onembed = nRank;
  int idist = 1;
  for(int i=0;i<rank;i++) idist *= nRank[i];
  int odist = idist;

  int howmany_u = T*2*16  ;
  int howmany_n = T*2* 4  ;
  int howmany_o = T*2*16*4;

  unsigned long long FFTlen_uloc = (unsigned long long)(2)*V*2*16  ;
  unsigned long long FFTlen_noet = (unsigned long long)(2)*V*2* 4  ;
  unsigned long long FFTlen_oneD = (unsigned long long)(2)*V*2*16*4;
  printf("FFT Buffer lengths:\n");
  printf(" u-local = %lld\n",FFTlen_uloc);
  printf(" noether = %lld\n",FFTlen_noet);
  printf(" oneD    = %lld\n",FFTlen_oneD);

  //-FFT Buffers
  Float *thrpFFT_uloc = (Float*) calloc(FFTlen_uloc,sizeof(Float));
  Float *thrpFFT_noet = (Float*) calloc(FFTlen_noet,sizeof(Float));
  Float *thrpFFT_oneD = (Float*) calloc(FFTlen_oneD,sizeof(Float));
  if( thrpFFT_uloc == NULL ){
    fprintf(stderr,"Cannot allocate thrpFFT_uloc. Exiting\n");
    exit(-1);
  }
  if( thrpFFT_noet == NULL ){
    fprintf(stderr,"Cannot allocate thrpFFT_noet. Exiting\n");
    exit(-1);
  }
  if( thrpFFT_oneD == NULL ){
    fprintf(stderr,"Cannot allocate thrpFFT_oneD. Exiting\n");
    exit(-1);
  }

  unsigned long long momlen_uloc = (unsigned long long)(2)*T*Nmoms*2*16*Ntsink*NprojMax;
  unsigned long long momlen_noet = (unsigned long long)(2)*T*Nmoms*2* 4*Ntsink*NprojMax;
  unsigned long long momlen_oneD = (unsigned long long)(2)*T*Nmoms*2*16*Ntsink*NprojMax;
  printf("Momfer lengths:\n");
  printf(" u-local = %lld\n",momlen_uloc);
  printf(" noether = %lld\n",momlen_noet);
  printf(" oneD    = %lld\n",momlen_oneD);

  Float *thrpMom_uloc  = (Float*)  calloc(momlen_uloc,sizeof(Float));
  Float *thrpMom_noet  = (Float*)  calloc(momlen_noet,sizeof(Float));
  Float **thrpMom_oneD = (Float**) calloc(4,sizeof(Float*));
  if( thrpMom_uloc == NULL ){
    fprintf(stderr,"Cannot allocate thrpMom_uloc. Exiting\n");
    exit(-1);
  }
  if( thrpMom_noet == NULL ){
    fprintf(stderr,"Cannot allocate thrpMom_noet. Exiting\n");
    exit(-1);
  }
  if( thrpMom_oneD == NULL ){
    fprintf(stderr,"Cannot allocate thrpMom_oneD. Exiting\n");
    exit(-1);
  }
  for(int mu=0;mu<4;mu++){
    thrpMom_oneD[mu] = (Float*) calloc(momlen_oneD,sizeof(Float));
    if( thrpMom_oneD[mu] == NULL ){
      fprintf(stderr,"Cannot allocate thrpMom_oneD[%d]. Exiting\n",mu);
      exit(-1);
    }
  }

  for(int its=0;its<Ntsink;its++){
    for(int ipr=0;ipr<Nproj[its];ipr++){
      unsigned long long posu = (unsigned long long)(2)*V*2*16  *its + (unsigned long long)(2)*V*2*16  *Ntsink*ipr;
      unsigned long long posn = (unsigned long long)(2)*V*2* 4  *its + (unsigned long long)(2)*V*2* 4  *Ntsink*ipr;
      unsigned long long poso = (unsigned long long)(2)*V*2*16*4*its + (unsigned long long)(2)*V*2*16*4*Ntsink*ipr;
      
      if( typeid(Float) == typeid(float) ){
	//-Create the FFT plans
	fftwf_plan FFTplanMany_u = fftwf_plan_many_dft(rank, nRank, howmany_u, (fftwf_complex*) thrpFFT_uloc, inembed, istride, idist, (fftwf_complex*) thrpFFT_uloc, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);
	fftwf_plan FFTplanMany_n = fftwf_plan_many_dft(rank, nRank, howmany_n, (fftwf_complex*) thrpFFT_noet, inembed, istride, idist, (fftwf_complex*) thrpFFT_noet, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);
	fftwf_plan FFTplanMany_o = fftwf_plan_many_dft(rank, nRank, howmany_o, (fftwf_complex*) thrpFFT_oneD, inembed, istride, idist, (fftwf_complex*) thrpFFT_oneD, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);

	memcpy(thrpFFT_uloc,&(thrpBuf_uloc[posu]),sizeof(Float)*FFTlen_uloc);
	memcpy(thrpFFT_noet,&(thrpBuf_noet[posn]),sizeof(Float)*FFTlen_noet);
	memcpy(thrpFFT_oneD,&(thrpBuf_oneD[poso]),sizeof(Float)*FFTlen_oneD);

	fftwf_execute(FFTplanMany_u);   //-Perform FFT
	fftwf_execute(FFTplanMany_n);   //-Perform FFT
	fftwf_execute(FFTplanMany_o);   //-Perform FFT
	fftwf_destroy_plan(FFTplanMany_u);   //-Destroy the plan
	fftwf_destroy_plan(FFTplanMany_n);   //-Destroy the plan
	fftwf_destroy_plan(FFTplanMany_o);   //-Destroy the plan
      }
      if( typeid(Float) == typeid(double) ){
	//-Create the FFT plans
	fftw_plan FFTplanMany_u = fftw_plan_many_dft(rank, nRank, howmany_u, (fftw_complex*) thrpFFT_uloc, inembed, istride, idist, (fftw_complex*) thrpFFT_uloc, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);
	fftw_plan FFTplanMany_n = fftw_plan_many_dft(rank, nRank, howmany_n, (fftw_complex*) thrpFFT_noet, inembed, istride, idist, (fftw_complex*) thrpFFT_noet, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);
	fftw_plan FFTplanMany_o = fftw_plan_many_dft(rank, nRank, howmany_o, (fftw_complex*) thrpFFT_oneD, inembed, istride, idist, (fftw_complex*) thrpFFT_oneD, onembed, ostride, odist, FFTW_BACKWARD, FFTW_MEASURE);

	memcpy(thrpFFT_uloc,&(thrpBuf_uloc[posu]),sizeof(Float)*FFTlen_uloc);
	memcpy(thrpFFT_noet,&(thrpBuf_noet[posn]),sizeof(Float)*FFTlen_noet);
	memcpy(thrpFFT_oneD,&(thrpBuf_oneD[poso]),sizeof(Float)*FFTlen_oneD);

	fftw_execute(FFTplanMany_u);   //-Perform FFT
	fftw_execute(FFTplanMany_n);   //-Perform FFT
	fftw_execute(FFTplanMany_o);   //-Perform FFT
	fftw_destroy_plan(FFTplanMany_u);   //-Destroy the plan
	fftw_destroy_plan(FFTplanMany_n);   //-Destroy the plan
	fftw_destroy_plan(FFTplanMany_o);   //-Destroy the plan
      }

      //-Copy the thrpFFT buffer to the thrpMom buffers
      printf("Copying to write buffer for tsink = %d, proj %s...\n",tsink[its],thrp_proj_type[proj_list[its][ipr]]);
      copyThrpToWriteBuf(thrpMom_uloc, thrpFFT_uloc, momALL, mom, its, ipr, 0, thrpInfo, THRP_LOCAL);
      copyThrpToWriteBuf(thrpMom_noet, thrpFFT_noet, momALL, mom, its, ipr, 0, thrpInfo, THRP_NOETHER);
      for(int mu=0;mu<4;mu++)
	copyThrpToWriteBuf(thrpMom_oneD[mu], thrpFFT_oneD, momALL, mom, its, ipr, mu, thrpInfo, THRP_ONED);

    }//-ipr
  }//-its

  printf("Fourier Transform completed successfully\n");
  //---------------------------------------------- 
   
  //-Write the output files
  char *file;

  if(h5out){
    asprintf(&file,"%s.%04d_neutron_Qsq%d_SS.%02d.%02d.%02d.%02d.h5",outpre, conf, Qsq, src[0], src[1], src[2], src[3]);
    writeThrp_HDF5(thrpMom_uloc, thrpMom_noet, thrpMom_oneD, file, mom, thrpInfo);
  }
  else{
    for(int itype=0;itype<3;itype++){
      THRP_TYPE type = (THRP_TYPE) itype;

      for(int its=0;its<Ntsink;its++){
	for(int proj=0;proj<Nproj[its];proj++){
	  for(int uORd=0;uORd<2;uORd++){
	    asprintf(&file,"%s.%04d_tsink%d_proj%s.neutron.%s.%s.SS.%02d.%02d.%02d.%02d.dat",outpre, conf, tsink[its], thrp_proj_type[proj_list[its][proj]], thrp_part[uORd], thrp_type[itype], src[0], src[1], src[2], src[3]);
	  
	    FILE *outp;
	    if( (outp=fopen(file,"w")) == NULL ){
	      fprintf(stderr,"Cannot open %s for writing\n",file);
	      exit(-1);
	    }

	    if(type==THRP_LOCAL){
	      int Mel = 16;
	      for(int imom=0;imom<Nmoms;imom++){
		for(int t=0;t<=tsink[its];t++){
		  for(int im=0;im<Mel;im++){
		    unsigned long long Mompos_Re = (unsigned long long)(0ULL + 2*im + 2*Mel*t  + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its + 2*Mel*T*Nmoms*2*Ntsink*proj);
		    unsigned long long Mompos_Im = (unsigned long long)(1ULL + 2*im + 2*Mel*t  + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its + 2*Mel*T*Nmoms*2*Ntsink*proj);

		    fprintf(outp,"%02d\t%02d\t%+d %+d %+d\t%+e %+e\n",t,im, mom[imom][0],mom[imom][1],mom[imom][2], thrpMom_uloc[Mompos_Re], thrpMom_uloc[Mompos_Im]);
		  }}}
	    }//-if THRP_LOCAL
	    else if(type==THRP_NOETHER){
	      int Mel = 4;
	      for(int imom=0;imom<Nmoms;imom++){
		for(int t=0;t<=tsink[its];t++){
		  for(int im=0;im<Mel;im++){
		    unsigned long long Mompos_Re = (unsigned long long)(0ULL + 2*im + 2*Mel*t  + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its + 2*Mel*T*Nmoms*2*Ntsink*proj);
		    unsigned long long Mompos_Im = (unsigned long long)(1ULL + 2*im + 2*Mel*t  + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its + 2*Mel*T*Nmoms*2*Ntsink*proj);

		    fprintf(outp,"%02d\t%02d\t%+d %+d %+d\t%+e %+e\n",t,im, mom[imom][0],mom[imom][1],mom[imom][2], thrpMom_noet[Mompos_Re], thrpMom_noet[Mompos_Im]);
		  }}}
	    }//-if THRP_NOETHER
	    else if(type==THRP_ONED){
	      int Mel = 16;
	      for(int imom=0;imom<Nmoms;imom++){
		for(int mu=0;mu<4;mu++){
		  for(int t=0;t<=tsink[its];t++){
		    for(int im=0;im<Mel;im++){
		      unsigned long long Mompos_Re = (unsigned long long)(0ULL + 2*im + 2*Mel*t  + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its + 2*Mel*T*Nmoms*2*Ntsink*proj);
		      unsigned long long Mompos_Im = (unsigned long long)(1ULL + 2*im + 2*Mel*t  + 2*Mel*T*imom + 2*Mel*T*Nmoms*uORd + 2*Mel*T*Nmoms*2*its + 2*Mel*T*Nmoms*2*Ntsink*proj);

		      fprintf(outp,"%02d\t%02d\t%02d\t%+d %+d %+d\t%+e %+e\n",mu, t,im, mom[imom][0],mom[imom][1],mom[imom][2], thrpMom_oneD[mu][Mompos_Re], thrpMom_oneD[mu][Mompos_Im]);
		    }}}}
	    }//-if THRP_ONED

	    fclose(outp);
	  }//-uORd
	}//-proj
      }//-its
    }//-itype

  }//-else h5out


  free(thrpBuf_uloc);
  free(thrpBuf_noet);
  free(thrpBuf_oneD);
  free(thrpFFT_uloc);
  free(thrpFFT_noet);
  free(thrpFFT_oneD);
  free(thrpMom_uloc);
  free(thrpMom_noet);
  for(int mu=0;mu<4;mu++) free(thrpMom_oneD[mu]);
  free(thrpMom_oneD);

  for(int i=0;i<MAX_MOM;i++) free(mom[i]);
  free(mom);
  for(int i=0;i<SpV;i++) free(momALL[i]);
  free(momALL);

  printf("Extracting Three-point function and FT completed successfully.\n");

  return 0;
}

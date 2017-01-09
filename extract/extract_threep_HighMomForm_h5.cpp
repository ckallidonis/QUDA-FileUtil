/* Christos Kallidonis                                  */
/* March 2015                                           */
/* This code takes as input a .h5 file and extracts the */
/* desired three-point function data into an ASCII file */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>

#define DIR 4
#define Ntype 3
#define QSQ_STR "Qsq"
#define NMOM_STR "Nmoms"
#define MOM_DSET "Momenta_list_xyz"


void usage(char exe[]){
  printf("%s: <.h5-file> <output_file> <thrp-type (ultra_local, noether, oneD)> <proj> <part (up,down)> <conf_traj> <src_pos (sx00sy00xz00st00)> <Tsink>\n",exe);
  exit(-1);
}

int main(int argc, char *argv[]){

  if(argc!=9) usage(argv[0]);

  char thrp_type[3][256];
  int M;

  strcpy(thrp_type[0],"ultra_local");
  strcpy(thrp_type[1],"noether");
  strcpy(thrp_type[2],"oneD");

  enum THRP_TYPE{THRP_LOCAL,THRP_NOETHER,THRP_ONED};

  char *h5_file, *out_file, *thrp, *proj, *part, *conf, *src_pos;
  asprintf(&h5_file ,"%s",argv[1]);
  asprintf(&out_file,"%s",argv[2]);
  asprintf(&thrp    ,"%s",argv[3]);
  asprintf(&proj    ,"%s",argv[4]);
  asprintf(&part    ,"%s",argv[5]);
  asprintf(&conf    ,"%s",argv[6]);
  asprintf(&src_pos ,"%s",argv[7]);

  int tsink = atoi(argv[8]);
  int T = tsink+1;

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

  printf("Got the following input:\n");
  printf("h5_file: %s\n",h5_file);
  printf("out_file: %s\n",out_file);
  printf("thrp: %d - %s\n",dt,thrp);
  printf("proj: %s\n",proj);
  printf("part: %s\n",part);
  printf("conf traj: %s\n",conf);
  printf("src_pos = %s\n",src_pos);
  printf("tsink = %02d\n",tsink);
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

  FILE *out = fopen(out_file,"w");
  if(out==NULL){
    fprintf(stderr,"Cannot open %s for writing. Exiting.\n",out_file);
    exit(-1);
  }
  
  if(type==THRP_ONED){
    float *thrpBuf[DIR];

    for(int mu=0;mu<DIR;mu++){
      thrpBuf[mu] = (float*) malloc(T*Nmoms*M*2*sizeof(float));

      asprintf(&group_dir,"/conf_%s/%s/tsink_%02d/proj_%s/%s/%s/dir_%02d",conf,src_pos,tsink,proj,part,thrp,mu);
      printf("Opening group: %s\n",group_dir);
      group_id = H5Gopen(file_id, group_dir, H5P_DEFAULT);
      dset_id = H5Dopen(group_id, "threep", H5P_DEFAULT);
    
      herr_t status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, thrpBuf[mu]);

      if (status<0){
	fprintf (stderr, "Dataset read failed!\n");
	status = H5Gclose(group_id);
	status = H5Dclose(dset_id);
	status = H5Fclose(file_id);
	fclose(out);
	free(thrpBuf[mu]);
	exit(-1);
      }

      H5Dclose(dset_id);
      H5Gclose(group_id);
    }//-mu

    for(int imom=0;imom<Nmoms;imom++){
      for(int mu=0;mu<DIR;mu++){
	for(int t=0;t<T;t++){
	  for(int m=0;m<M;m++){
	    fprintf(out,"%02d\t%02d\t%02d\t%+d %+d %+d\t%+e %+e\n",mu,t,m,moms[0+3*imom],moms[1+3*imom],moms[2+3*imom],
		    thrpBuf[mu][0 + 2*m + 2*M*imom + 2*M*Nmoms*t],thrpBuf[mu][1 + 2*m + 2*M*imom + 2*M*Nmoms*t]); 
	  }
	}
      }
    }//-imom

    for(int mu=0;mu<DIR;mu++) free(thrpBuf[mu]);
  }
  else{
    float *thrpBuf = (float*) malloc(T*Nmoms*M*2*sizeof(float));

    asprintf(&group_dir,"/conf_%s/%s/tsink_%02d/proj_%s/%s/%s",conf,src_pos,tsink,proj,part,thrp);
    printf("Opening group: %s\n",group_dir);
    group_id = H5Gopen(file_id, group_dir, H5P_DEFAULT);
    dset_id = H5Dopen(group_id, "threep", H5P_DEFAULT);
    
    herr_t status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, thrpBuf);

    if (status<0){
      fprintf (stderr, "Dataset read failed!\n");
      status = H5Gclose(group_id);
      status = H5Dclose(dset_id);
      status = H5Fclose(file_id);
      fclose(out);
      free(thrpBuf);
      exit(-1);
    }

    for(int imom=0;imom<Nmoms;imom++){
      for(int t=0;t<T;t++){
	for(int m=0;m<M;m++){
	  fprintf(out,"%02d\t%02d\t%+d %+d %+d\t%+e %+e\n",t,m,moms[0+3*imom],moms[1+3*imom],moms[2+3*imom],
		  thrpBuf[0 + 2*m + 2*M*imom + 2*M*Nmoms*t],thrpBuf[1 + 2*m + 2*M*imom + 2*M*Nmoms*t]); 
	}
      }
    }//-imom

    H5Dclose(dset_id);
    H5Gclose(group_id);

    free(thrpBuf);
  }

  H5Fclose(file_id);

  fclose(out);
  free(moms);

  printf("Extracting Three-point function completed successfully.\n");

  return 0;
}

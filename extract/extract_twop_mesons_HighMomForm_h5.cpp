/* Christos Kallidonis                                  */
/* March 2015                                           */
/* This code takes as input a .h5 file and extracts the */
/* desired two-point function data into an ASCII file   */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>

#define Ntype 10
#define QSQ_STR "Qsq"
#define NMOM_STR "Nmoms"
#define MOM_DSET "Momenta_list_xyz"

void usage(char exe[]){
  printf("%s:\n<.h5-file>\n<output_file>\n<meson-type>\n  pseudoscalar\n  scalar\n  g5g1\n  g5g2\n  g5g3\n  g5g4\n  g1\n  g2\n  g3\n  g4\n",exe);
  printf("<conf_traj>\n<src_pos (sx00sy00xz00st00)>\n<T>\n<extract_type, 0:both, 1:particle, 2:anti-particle>\n");
  exit(-1);
}

int main(int argc, char *argv[]){

  if(argc!=8) usage(argv[0]);

  char twop_type[Ntype][256];

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

  char *h5_file, *out_file, *twop, *conf, *src_pos;
  asprintf(&h5_file ,"%s",argv[1]);
  asprintf(&out_file,"%s",argv[2]);
  asprintf(&twop    ,"%s",argv[3]);
  asprintf(&conf    ,"%s",argv[4]);
  asprintf(&src_pos ,"%s",argv[5]);

  int T = atoi(argv[6]);
  int xtype = atoi(argv[7]);

  bool twopOK = false;
  int dt;
  for(int i=0;(i<Ntype && !twopOK);i++)
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
  printf("out_file: %s\n",out_file);
  printf("twop: %d - %s\n",dt,twop);
  printf("conf traj: %s\n",conf);
  printf("src pos: %s\n",src_pos);
  printf("T = %02d\n",T);
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
  char *dset_name;

  int Np,Ns;
  if(xtype==0) Np = 2;
  else Np = 1;

  Ns = xtype/2;

  float *twopBuf = (float*) malloc(Np*T*Nmoms*2*sizeof(float));

  FILE *out = fopen(out_file,"w");
  if(out==NULL){
    fprintf(stderr,"Cannot open %s for writing. Exiting.\n",out_file);
    exit(-1);
  }

  for(int ip=Ns;ip<(Np+Ns);ip++){
    asprintf(&dset_name,"twop_meson_%d",ip+1);

    asprintf(&group_dir,"/conf_%s/%s/%s",conf,src_pos,twop);
    group_id = H5Gopen(file_id, group_dir, H5P_DEFAULT);
    dset_id = H5Dopen(group_id, dset_name, H5P_DEFAULT);
    
    herr_t status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(twopBuf[(ip%Np)*T*Nmoms*2]));

    if (status<0){
      fprintf (stderr, "Dataset read failed!\n");
      status = H5Gclose(group_id);
      status = H5Dclose(dset_id);
      status = H5Fclose(file_id);
      fclose(out);
      free(twopBuf);
      exit(-1);
    }

    H5Dclose(dset_id);
    H5Gclose(group_id);
  }

  H5Fclose(file_id);
  //------------------------------------------------------------------

  if(xtype==0){
    for(int imom=0;imom<Nmoms;imom++){
      for(int t=0;t<T;t++){
	fprintf(out,"%d\t%02d\t%+d %+d %+d\t%+e %+e\t%+e %+e\n",dt,t,moms[0+3*imom],moms[1+3*imom],moms[2+3*imom],
		twopBuf[0 + 2*imom + 2*Nmoms*t + 2*Nmoms*T*0],twopBuf[1 + 2*imom + 2*Nmoms*t + 2*Nmoms*T*0],
		twopBuf[0 + 2*imom + 2*Nmoms*t + 2*Nmoms*T*1],twopBuf[1 + 2*imom + 2*Nmoms*t + 2*Nmoms*T*1]);
      }//-t
    }//-imom
  }
  else{
    for(int imom=0;imom<Nmoms;imom++){
      for(int t=0;t<T;t++){
	fprintf(out,"%d\t%02d\t%+d %+d %+d\t%+e %+e\n",dt,t,moms[0+3*imom],moms[1+3*imom],moms[2+3*imom],twopBuf[0 + 2*imom + 2*Nmoms*t],twopBuf[1 + 2*imom + 2*Nmoms*t]);
      }
    }
  }


  fclose(out);
  free(twopBuf);
  free(moms);

  printf("Extracting Two-point function completed successfully.\n");

  return 0;
}

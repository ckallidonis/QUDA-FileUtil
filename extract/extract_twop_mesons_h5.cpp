/* Christos Kallidonis                                  */
/* March 2015                                           */
/* This code takes as input a .h5 file and extracts the */
/* desired two-point function data into an ASCII file   */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>

#define Ntype 10

void usage(char exe[]){
  printf("%s:\n<.h5-file>\n<output_file>\n<meson-type>\n  pseudoscalar\n  scalar\n  g5g1\n  g5g2\n  g5g3\n  g5g4\n  g1\n  g2\n  g3\n  g4\n",exe);
  printf("<momx> <momy> <momz>\n<conf_traj>\n<src_pos (sx00sy00xz00st00)>\n<T>\n<extract_type, 0:both, 1:particle, 2:anti-particle>\n");
  exit(-1);
}

int main(int argc, char *argv[]){

  if(argc!=11) usage(argv[0]);

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

  char *h5_file, *out_file, *twop, *mom[3], *conf, *src_pos;
  asprintf(&h5_file ,"%s",argv[1]);
  asprintf(&out_file,"%s",argv[2]);
  asprintf(&twop    ,"%s",argv[3]);
  asprintf(&mom[0]  ,"%s",argv[4]);
  asprintf(&mom[1]  ,"%s",argv[5]);
  asprintf(&mom[2]  ,"%s",argv[6]);
  asprintf(&conf    ,"%s",argv[7]);
  asprintf(&src_pos ,"%s",argv[8]);

  int T = atoi(argv[9]);
  int xtype = atoi(argv[10]);

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
  printf("momentum: %s %s %s\n",mom[0],mom[1],mom[2]);
  printf("conf traj: %s\n",conf);
  printf("src pos: %s\n",src_pos);
  printf("T = %02d\n",T);
  //-----------------------------------------

  //-Open the h5 file
  hid_t file_id = H5Fopen(h5_file, H5F_ACC_RDONLY, H5P_DEFAULT);

  //-Open the desired dataset and write it in the ASCII file
  hid_t group_id;
  hid_t dset_id;
  char *group_dir;
  char *dset_name;

  int Np,Ns;
  if(xtype==0) Np = 2;
  else Np = 1;

  Ns = xtype/2;

  float *twopBuf = (float*) malloc(Np*T*2*sizeof(float));

  FILE *out = fopen(out_file,"w");
  if(out==NULL){
    fprintf(stderr,"Cannot open %s for writing. Exiting.\n",out_file);
    exit(-1);
  }

  for(int ip=Ns;ip<(Np+Ns);ip++){
    asprintf(&dset_name,"twop_meson_%d",ip+1);

    asprintf(&group_dir,"/conf_%s/%s/%s/mom_xyz_%s_%s_%s",conf,src_pos,twop,mom[0],mom[1],mom[2]);
    group_id = H5Gopen(file_id, group_dir, H5P_DEFAULT);
    dset_id = H5Dopen(group_id, dset_name, H5P_DEFAULT);
    
    herr_t status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(twopBuf[(ip%Np)*T*2]));

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

  if(xtype==0){
    for(int t=0;t<T;t++)
      fprintf(out,"%d\t%02d\t%s %s %s\t%+e %+e\t%+e %+e\n",dt,t,mom[0],mom[1],mom[2],
	      twopBuf[0 + 2*t + 2*T*0],twopBuf[1 + 2*t + 2*T*0],
	      twopBuf[0 + 2*t + 2*T*1],twopBuf[1 + 2*t + 2*T*1]);
  }
  else{
    for(int t=0;t<T;t++)
      fprintf(out,"%d\t%02d\t%s %s %s\t%+e %+e\n",dt,t,mom[0],mom[1],mom[2],twopBuf[0 + 2*t],twopBuf[1 + 2*t]);
  }
  H5Fclose(file_id);

  fclose(out);
  free(twopBuf);

  printf("Extracting Two-point function completed successfully.\n");

  return 0;
}

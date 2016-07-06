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

void usage(char exe[]){
  printf("%s: <.h5-file> <output_file> <thrp-type (ultra_local, noether, oneD)> <proj> <part (up,down)> <momx> <momy> <momz> <conf_traj> <src_pos (sx00sy00xz00st00)> <Tsink>\n",exe);
  exit(-1);
}

int main(int argc, char *argv[]){

  if(argc!=12) usage(argv[0]);

  char thrp_type[3][256];
  int M;

  strcpy(thrp_type[0],"ultra_local");
  strcpy(thrp_type[1],"noether");
  strcpy(thrp_type[2],"oneD");

  enum THRP_TYPE{THRP_LOCAL,THRP_NOETHER,THRP_ONED};

  char *h5_file, *out_file, *thrp, *proj, *part, *mom[3], *conf, *src_pos;
  asprintf(&h5_file ,"%s",argv[1]);
  asprintf(&out_file,"%s",argv[2]);
  asprintf(&thrp    ,"%s",argv[3]);
  asprintf(&proj    ,"%s",argv[4]);
  asprintf(&part    ,"%s",argv[5]);
  asprintf(&mom[0]  ,"%s",argv[6]);
  asprintf(&mom[1]  ,"%s",argv[7]);
  asprintf(&mom[2]  ,"%s",argv[8]);
  asprintf(&conf    ,"%s",argv[9]);
  asprintf(&src_pos ,"%s",argv[10]);

  int tsink = atoi(argv[11]);
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
  printf("momentum: %s %s %s\n",mom[0],mom[1],mom[2]);
  printf("conf traj: %s\n",conf);
  printf("src_pos = %s\n",src_pos);
  printf("tsink = %02d\n",tsink);
  //-----------------------------------------

  //-Open the h5 file
  hid_t file_id = H5Fopen(h5_file, H5F_ACC_RDONLY, H5P_DEFAULT);

  //-Open the desired dataset and write it in the ASCII file
  hid_t group_id;
  hid_t dset_id;
  char *group_dir;

  float *thrpBuf = (float*) malloc(T*M*2*sizeof(float));

  FILE *out = fopen(out_file,"w");
  if(out==NULL){
    fprintf(stderr,"Cannot open %s for writing. Exiting.\n",out_file);
    exit(-1);
  }
  
  if(type==THRP_ONED){
    for(int mu=0;mu<DIR;mu++){
      asprintf(&group_dir,"/conf_%s/%s/tsink_%02d/proj_%s/%s/%s/mom_xyz_%s_%s_%s/dir_%02d",conf,src_pos,tsink,proj,part,thrp,mom[0],mom[1],mom[2],mu);
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

      for(int t=0;t<T;t++)
	for(int m=0;m<M;m++) fprintf(out,"%02d\t%02d\t%02d\t%s %s %s\t%+e %+e\n",mu,t,m,mom[0],mom[1],mom[2],thrpBuf[0+2*m+2*M*t],thrpBuf[1+2*m+2*M*t]); 

      H5Dclose(dset_id);
      H5Gclose(group_id);
    }//-mu
  }
  else{
    asprintf(&group_dir,"/conf_%s/%s/tsink_%02d/proj_%s/%s/%s/mom_xyz_%s_%s_%s",conf,src_pos,tsink,proj,part,thrp,mom[0],mom[1],mom[2]);
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

    for(int t=0;t<T;t++)
      for(int m=0;m<M;m++) fprintf(out,"%02d\t%02d\t%s %s %s\t%+e %+e\n",t,m,mom[0],mom[1],mom[2],thrpBuf[0+2*m+2*M*t],thrpBuf[1+2*m+2*M*t]); 

    H5Dclose(dset_id);
    H5Gclose(group_id);
  }

  H5Fclose(file_id);

  fclose(out);
  free(thrpBuf);

  printf("Extracting Three-point function completed successfully.\n");

  return 0;
}

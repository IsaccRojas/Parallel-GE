/* Original author: Jeevan Joseph, ece/cs 566 student, S'13 */
/* modified by Shafigh Mehraeen, Chemical Engineering, UIC, 5/1/2018 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define ROWS 4096
#define COLS 4096
//#define ROWS 2048
//#define COLS 2048
//#define ROWS 1024
//#define COLS 1024
//#define ROWS 512
//#define COLS 512
//#define ROWS 32
//#define COLS 32
//#define ROWS 16
//#define COLS 16
//#define ROWS 8
//#define COLS 8
//#define ROWS 4
//#define COLS 4

float generate_rand(float file_r);

int main(int argc, char* argv[]) {
    int rank, size, num,i,j,k;
    float final_number, mant_r, exp_r, b_r, file_r ;
    int limit_per_no = pow(2,16), sign_exp, sign_mant;
    char op_file_name[100];
    FILE *fp;
    srand(time(NULL));
    
    int id=atoi(argv[1]);
    
    //allocate 2d ROWS x COLS array
    float *mat[ROWS];
    for (i=0; i<ROWS; i++) {mat[i] = (float *)malloc(COLS * sizeof(float));}
    
    //allocate 1d array of size ROWS (the vector)
    float *bMatrix;
    bMatrix=malloc(ROWS * sizeof(float));

    //generate a random float from 0.0f to 1.0f
    file_r = (rand()%100)/(float)100;
        
    for(i=0; i<ROWS;i++){
      //get valid random value and store in 2D matrix array
      for(j=0;j<COLS;j++){
        final_number = generate_rand(file_r);
        mat[i][j] = final_number;
      }

      //get valid random value and store in 1D vector array
      bMatrix[i] = generate_rand(file_r);
      if (abs(mat[i][i])<1e-6) {printf("diagonal element is small! Rerun the code!\n"); return 1;}
      else {

        //scale rest of row with diagonal index
        for (j=0; j<COLS; j++) {
          if (i!=j) {mat[i][j]=mat[i][j]/mat[i][i];}
        }

        //scale vector at i by diagonal
        bMatrix[i]=bMatrix[i]/mat[i][i];
        mat[i][i]=1;
      }
    }
    
    //open file for writing
    sprintf(op_file_name,"matrix_%d_%d.txt",ROWS,id);
    fp = fopen(op_file_name, "w");

    for(i=0; i<ROWS; i++){
      //write all values in 2D matrix array to file
      for(j=0; j<COLS; j++){
        fprintf(fp, "%011.4e", mat[i][j]);
        if((j+1)%COLS == 0 )
          fprintf(fp, "\n");
        else
          fprintf(fp, " ");
      }
      
      //write vector at end of file
      fprintf(fp, "%011.4e", bMatrix[i]);
      if(j!=COLS-1){
        fprintf(fp, "\n");
      }
    }
    
    fclose(fp);
    return 0;
}

float generate_rand(float file_r){
  float final_number, mant_r, exp_r;
  int sign_exp;
  do{
    mant_r = rand() % 65536;
    while(mant_r >10){
    mant_r = mant_r/(float)10;
        if (mant_r<10)
            break;
    }
    sign_exp = (rand() % 2 ==0)? -1:1 ;
    mant_r = mant_r * sign_exp;
    
    if(file_r <=.33){
        exp_r = (rand() % 7) -3;
    }
    else if(file_r <=.66){
        exp_r = (rand() % 7) -1;
    }
    else{
        exp_r = (rand() % 7) +2;
    }
    
    final_number = mant_r * pow(2, exp_r);
  }while(final_number == 0);
  return final_number;
}

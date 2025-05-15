
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("motif_region_file open error.\n");
    exit(1); //ファイルが開けなかった場合プログラムを終了
  }

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
    seq_num++;
  }
  return seq_num;
}

int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1);
  }

  while(fscanf(fp, "%s", buffer) != EOF){
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0';
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); 
    }else{
      strcpy(g_pro[gene_num].seq,buffer);
      gene_num++;
    }    
  }
  return gene_num;
}

void tablemaker(int seq_num, char A[MAX_SEQ_NUM][BUFSIZE])
{
    int i,j;
    float qx[4]={7519429.0/24374210,4637676.0/24374210,4637676.0/24374210,7519429.0/24374210};
    int a=strlen(A[0]);
    float oddsscore[4][BUFSIZE]={0};
    int number[4][BUFSIZE]={0};
    
    for(i=0;i<strlen(A[0]);i++)
    {
        for(j=0;j<seq_num;j++)
        {
          if(g_motif[j][i]=='T')
          {
            number[3][i]++;
          }
          if(g_motif[j][i]=='A')
          {
            number[0][i]++;
          }
          if(g_motif[j][i]=='C')
          {
            number[1][i]++;
          }
          if(g_motif[j][i]=='G')
          {
            number[2][i]++;
          }
        }
    }
for(i=0;i<4;i++)
    {
      for(j=0;j<strlen(A[0]);j++)
      {
        number[i][j]++;
      }
    }

    for(i=0;i<4;i++)
    {
      for(j=0;j<strlen(A[0]);j++)
      {
        printf("%d\n",number[i][j]);
      }

    }
    for(i=0;i<4;i++)
    {
      for(j=0;j<strlen(A[0]);j++)
      {
        oddsscore[i][j]=log(number[i][j]/(qx[i]*(seq_num+4)));
      }

    }
      
    
    for(i=0;i<4;i++)
    {
      for(j=0;j<strlen(A[0]);j++)
      {
        printf("%f\n",oddsscore[i][j]);
      }
      
    }
}

int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む

  printf("motif region:\n");
  for(int i = 0; i < seq_num; i++){
    printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
  }
  printf("\n");

  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  
  printf("promoter_sequence:\n");
  for(int i = 0; i < gene_num; i++){
    printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
    printf("%s\n", g_pro[i].seq);
  }

  tablemaker(seq_num, g_motif);


  return 0;
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define Base 4

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

void tablemaker(int seq_num, char A[MAX_SEQ_NUM][BUFSIZE], float oddsscore[Base][BUFSIZE])
{
    int i,j;
    float qx[Base]={7519429.0/24374210,4637676.0/24374210,4637676.0/24374210,7519429.0/24374210};
    int number[Base][BUFSIZE]={0};
    
    
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
for(i=0;i<Base;i++)
    {
      for(j=0;j<strlen(A[0]);j++)
      {
        number[i][j]++;
      }
    }

    for(i=0;i<Base;i++)
    {
      for(j=0;j<strlen(A[0]);j++)
      {
       // printf("%d\n",number[i][j]);
      }

    }
    for(i=0;i<Base;i++)
    {
      for(j=0;j<strlen(A[0]);j++)
      {
        oddsscore[i][j]=log(number[i][j]/(qx[i]*(seq_num+4)));
      }

    }
      
    
    for(i=0;i<Base;i++)
    {
      for(j=0;j<strlen(A[0]);j++)
      {
        printf("%f ",oddsscore[i][j]);
      }
      printf("\n");
    }
    printf("\n");
  }

void scansequence(float oddsscore[Base][BUFSIZE])
{
  int i,j,k,l;
  float sum[8][BUFSIZE]={0};
    for(j=0;j<8;j++)
    {
      float hit[10]={0};
      int order=0;
      for(k=0;k<strlen(g_pro[j].seq)-strlen(g_motif[0]);k++)
      {
        for(i=0;i<strlen(g_motif[0]);i++)
        {
          
          if(g_pro[j].seq[k+i]=='T')
          {
            sum[j][k]+=oddsscore[3][i];
          }
          if(g_pro[j].seq[k+i]=='A')
          {
            sum[j][k]+=oddsscore[0][i];
          }
          if(g_pro[j].seq[k+i]=='C')
          {
            sum[j][k]+=oddsscore[1][i];
          }
          if(g_pro[j].seq[k+i]=='G')
          {
            sum[j][k]+=oddsscore[2][i];
          }
          
        }
        if(sum[j][k]>=6)
          {

            printf("pro:%s\n",g_pro[j].name);
            printf("pos:%d\n",k);
            printf("hit(");
            for(l=k;l<k+strlen(g_motif[0]);l++)
            {
            printf("%c",g_pro[j].seq[l]);
            }
            printf(")= %f\n",sum[j][k]);
            printf("\n");
          }
      }

    }
    for(j=0;j<8;j++)
    {
      for(k=0;k<strlen(g_pro[j].seq)-strlen(g_motif[0]);k++)
      {
       // printf("%f ",sum[j][k]);
      }
      //printf("\n");
    }
}

int GetRandom(void)
{
  int a,i,min=0,max=24314210;
  char result[MAX_SEQ_NUM];
  for(i=0;i<strlen(g_motif[0]);i++)
  {
    a=min + (int)(rand() * (max - min + 1.0) / (1.0 + RAND_MAX));
    if(0<=a<7519429)
    {
      result[i]='A';
    }
    if(7519429<=a<15038858)
    {
      result[i]='T';
    }if(15038858<=a<19676534)
    {
      result[i]='G';
    }if(19676534<=a<24314210)
    {
      result[i]='C';
    }
  }
  printf("%s\n",result);
  
}

int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む
  float oddsscore[Base][BUFSIZE];
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

  tablemaker(seq_num, g_motif, oddsscore);
  scansequence(oddsscore);
  GetRandom();


  return 0;
}

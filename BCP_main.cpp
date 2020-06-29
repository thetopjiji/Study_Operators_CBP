//============================================================================
// Name        : BCP_main.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <cmath>
#include <math.h>
#include <complex>
#include <assert.h>
#include <iomanip>
#include <sys/times.h>
#include <unistd.h>
#include "BCP_mem.h"
using namespace std;
#define min(a,b) ((a)<(b)?(a):(b))

//data struct of solution
//typedef struct struct_individual {
//	int *permutation;
//	int *permutationNew;
//	double cbmp;
//} Individual;
//Individual *CurrentS;

//Variable on input
int alreadybest;
char *rep;
char *name_instance;		//-i
int seed;					//-seed
int pop;					//population
char resultsFile[100];
char name_final_result[256];
char file_G[256];
char file_E[256];
char file_D[256];
char file_S[256];
char benchmark[100];
/*******************************/
//variable for traiting the graph
int num_vex;
int num_unk;
int num_edge;
int max_deg;
int **start_end;
int **edge;
int *degree_node;
int **simple_edge;
//variable for generating the se and simple_edge
int *seq_node;
int **linked_list;
int num_read=0;
/*******************************/
//function declaration
void freedataStrucGraph(void);
void setdataStrucGraph(void);
void generate_linked_list(int node1, int node2);
int *get_vector(int size);
int **get_matrix(int num_row,int num_col);
int parameters(int count, char *arguments[]);
void read_fiche(void);
int get_CycD(int label1, int label2);
/*******************************/
int get_CycD(int label1, int label2){
	int absD;
	int cycD;
	absD=abs(label1-label2);
	cycD=min(absD, num_vex-absD);
	return cycD;
}


void freedataStrucGraph(void){
	//Free the **start_end, **edge, **simple_edge, *degree_node, CurrentS
	int i;
	for(i=0;i<2*num_edge;i++){
		free(edge[i]);
		edge[i]=NULL;
	}
	for(i=0;i<num_edge;i++){
		free(simple_edge[i]);
		simple_edge[i]=NULL;
	}
	free(start_end); 					start_end=NULL;
	free(edge); 						edge=NULL;
	free(simple_edge); 					simple_edge=NULL;
	free(degree_node);					degree_node=NULL;

}

void setdataStrucGraph(void){
	//Fill the **start_end, **edge, **simple_edge, *degree_node
	int i=0,j=0;
	int line_read=0;
	int p_now=0;
	int temp=0;
	/*****************************************/
	start_end=(int**)get_matrix(num_vex,2);
	edge=(int**)get_matrix(2*num_edge,2);
	simple_edge=(int**)get_matrix(num_edge,2);
	degree_node=(int*)get_vector(num_vex);
	/******************************************/
	for (i=0;i<num_vex;i++) degree_node[i]=0;
	for(i=0;i<num_vex;i++){
		p_now=seq_node[i];
		if (p_now!=-1){
		edge[line_read][0]=i;
		edge[line_read][1]=linked_list[p_now][0];
		start_end[i][0]=line_read;
		line_read++;
		p_now=linked_list[p_now][1];
		while(p_now!=-1){
			edge[line_read][0]=i;
			edge[line_read][1]=linked_list[p_now][0];
			line_read++;
			p_now=linked_list[p_now][1];
		}
		if (p_now==-1){
			start_end[i][1]=line_read-1;
		}
		}
		else {
			cout<<"there is a independent vertex"<<endl;
			exit(-1);
		}
	}
	/******************************************/
	j=0;
	for (i=0;i<2*num_edge;i++){
		if(edge[i][0]<edge[i][1]){
			simple_edge[j][0]=edge[i][0];
			simple_edge[j][1]=edge[i][1];
			degree_node[edge[i][0]]++;
			degree_node[edge[i][1]]++;
//			cout<<simple_edge[j][0]<<" "<<simple_edge[j][1]<<endl;
			j++;
		}
	}
	/******************************************/
	for (i=0;i<num_vex;i++){
		temp=start_end[i][1]-start_end[i][0]+1;
		if (temp>max_deg){
			max_deg=temp;
		}
	}
	/******************************************/
	//free  *seq_node, **link_list
	for (i=0;i<num_edge;i++){
		free(linked_list[i]);
		linked_list[i]=NULL;
	}
	free(seq_node);seq_node=NULL;

	for(i=0;i<2*num_edge;i++){
				free(linked_list[i]);
				linked_list[i]=NULL;
			}

	free(linked_list);linked_list=NULL;

}


void generate_linked_list(int node1, int node2){
	int j=0,k=0,temp=0;
	int *p_now=NULL;
	j=node1;
	k=node2;
	temp=num_read;
	if (seq_node[j-1]==-1){
		seq_node[j-1]=num_read;
		linked_list[num_read][0]=k-1;
		linked_list[num_read][1]=-1;
		num_read++;
	}
	else{
		p_now=&seq_node[j-1];
		while (*p_now!=-1){
			if (linked_list[*p_now][0]>(k-1)){
				linked_list[num_read][0]=k-1;
				linked_list[num_read][1]=*p_now;
				*p_now=temp;
				num_read++;
				break;
				}
			else{
				p_now=&linked_list[*p_now][1];
			}
		}

		if (*p_now==-1){
			linked_list[num_read][0]=k-1;
			linked_list[num_read][1]=-1;
			*p_now=temp;
			num_read++;
		}
	}
}


int *get_vector(int size){
	int *Poi;
	Poi= (int*)malloc(size*sizeof(int));
	if (Poi==NULL){
		cout<<"Memory error in get_vector"<<endl;
		exit(-1);
	}
	return Poi;
}

int **get_matrix(int num_row,int num_col){
	int **Poi,i;
	Poi=(int**)malloc(num_row*sizeof(int*));
	if (!Poi){
			cout<<"Memory error in get_matrix"<<endl;
			exit(-1);
		}
	for (i=0;i<num_row;i++){
		Poi[i]=(int*)malloc(num_col*sizeof(int));
		if (!Poi[i]){
					cout<<"Memory error in get_matrix"<<endl;
					exit(-1);
		}
	}
	return Poi;
}

int parameters(int count, char *arguments[])  {
	char *temp, filename[80] = "no";
	char *nf=filename;
//	char *wewant=filename;
//	char *token;
	pop=-1;
	alreadybest=-1;

	strcpy(resultsFile, filename);
	strcpy(benchmark, filename);
	strcpy(nf, filename);

	while (count != 1) {
		temp = arguments[count - 2];
		if (strcmp(temp,"-i") == 0) strcpy(benchmark, arguments[count - 1]);
		else
			if (strcmp(temp,"--seed") == 0) seed= atoi(arguments[count - 1]);
		else
			if (strcmp(temp,"-rep") == 0) rep=arguments[count - 1];
		else
			if (strcmp(temp,"-alb") == 0) alreadybest= atoi(arguments[count - 1]);
		else
			if (strcmp(temp,"-pop") == 0) pop=atoi(arguments[count - 1]);
		else {  // unknow parameter
			return 0;
		}
		count = count - 2;
	}
	if (strcmp(benchmark, "no") == 0 || (alreadybest==-1) || (pop==-1)) {
		printf("enter error\n");
		exit(-1);
	}
	strcpy(nf, benchmark);

//	token=strtok(nf,"/");
//	while(token!=NULL){
//		wewant=token;
//		token=strtok(NULL,"/");
//	}
//	if(strcmp(wewant,"path100.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"cycle650.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"cycle1000.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"T-dwt__592.mtx.rnd") == 0) alreadybest=7;
//	else if(strcmp(wewant,"cycle200.rnd") == 0) alreadybest=1;
//
//	else if(strcmp(wewant,"caterpillar29.rnd") == 0) alreadybest=24;
//	else if(strcmp(wewant,"hypercube11.rnd") == 0) alreadybest=526;
//	else if(strcmp(wewant,"cycle300.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"path825.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"path200.rnd") == 0) alreadybest=1;
//
//	else if(strcmp(wewant,"X-can__715.mtx.rnd") == 0) alreadybest=52;
//	else if(strcmp(wewant,"mesh2D20x50.rnd") == 0) alreadybest=20;
//	else if(strcmp(wewant,"Q-494_bus.mtx.rnd") == 0) alreadybest=5;
//	else if(strcmp(wewant,"W-685_bus.mtx.rnd") == 0) alreadybest=6;
//	else if(strcmp(wewant,"U-662_bus.mtx.rnd") == 0) alreadybest=5;
//
//	else if(strcmp(wewant,"path650.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"mesh2D5x25.rnd") == 0) alreadybest=5;
//	else if(strcmp(wewant,"tree2x9.rnd") == 0) alreadybest=57;
//	else if(strcmp(wewant,"P-can__445.mtx.rnd") == 0) alreadybest=6;
//	else if(strcmp(wewant,"mesh3D12x12x12.rnd") == 0) alreadybest=114;
//	else alreadybest=1;
	return 1;
}

void read_fiche(void){

	ifstream FIC;
	FIC.open(benchmark);
	 if ( FIC.fail() ){
		 cout << "### No way,check your fiche " << benchmark << endl;
		 exit(-1);
	 }
	 char line1[100];
     FIC.getline(line1,100,'\n'); /*ignore the first line*/
	 FIC>>num_vex>>num_unk>>num_edge; /*get the number of vertex and edge*/
	 if ( FIC.eof() ){
		 cout << "### Your fiche is empty " << benchmark << endl;
		 exit(-1);
	 }
//	     cout<<max_vertex1<<max_vertex2<<num_edge<<endl;
	 seq_node=(int*)get_vector(num_vex);
	 linked_list=(int**)get_matrix(2*num_edge,2);

	 int i;
	 for (i=0;i<num_vex;i++) seq_node[i]=-1;
	 for (i=0;i<2*num_edge;i++) {
		 linked_list[i][0]=-1;
		 linked_list[i][1]=-1;
	 }

	 int x1,x2,t;
	 int cnt = 0;
	 while(!FIC.eof() && cnt < num_edge){
		 FIC>>x1>>x2;
		 /*wonder if the node is out of range*/
		 if ((x2<1)||(x1<1)||(x1>num_vex)||(x2>num_vex)) {
			 cout<<"the number of node is out of range"<<endl;
			 exit(-1);
		 }
	    /*ensure the x1<x2*/
		if(x1!=x2){
			if (x1>x2){
				 t=x1;
				 x1=x2;
				 x2=t;
			 }
			generate_linked_list(x1,x2);
			generate_linked_list(x2,x1);
	    }
		 cnt++;
	 }
	 assert(cnt == num_edge);
	 if(num_read/2!=num_edge){
		 cout<<"the number of edge is not correct"<<endl;
		 exit(-1);
	 }
	 FIC.close();
}

int main(int argc, char **argv) {
	char final[50]="F";
	char Ginal[50]="G";
	char Einal[50]="E";
	char Dinal[50]="D";
	char Sinal[50]="S";
	int cb_final;
	if(parameters(argc, argv)==0) exit(-1);
	printf("graph_name=%s seed=%d best_found=%d population=%d \n",benchmark, seed, alreadybest ,pop);
	srand(seed);
//  srand((unsigned int)time(NULL));
	read_fiche();						//read graph and fill the link_list and seq_node
    strcpy(name_final_result,benchmark);
    strcat(name_final_result,final);
    strcat(name_final_result,rep);
    printf("%s\n",name_final_result);
/*******************************************************/
    strcpy(file_G,benchmark);
    strcat(file_G,Ginal);
    strcat(file_G,rep);
    printf("%s\n",file_G);
/*******************************************************/
	strcpy(file_E,benchmark);
	strcat(file_E,Einal);
	strcat(file_E,rep);
	printf("%s\n",file_E);
/*******************************************************/
/*******************************************************/
	strcpy(file_D,benchmark);
	strcat(file_D,Dinal);
	strcat(file_D,rep);
	printf("%s\n",file_D);
/*******************************************************/
/*******************************************************/
	strcpy(file_S,benchmark);
	strcat(file_S,Sinal);
	strcat(file_S,rep);
	printf("%s\n",file_S);
/*******************************************************/

	setdataStrucGraph();
	cb_final=memSearch(alreadybest, pop);
	freedataStrucGraph();


	cout << "program finished" << endl; // prints !!!Hello World!!!
	cout << "cb_final="<<cb_final << endl; // prints !!!Hello World!!!
	return 0;
}

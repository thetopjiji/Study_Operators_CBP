//============================================================================
// Name        : BCP_mem.cpp
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
#include "BCP_X.h"

using namespace std;
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

//extern function
extern int *get_vector(int size);
extern int **get_matrix(int num_row,int num_col);
extern int get_CycD(int label1, int label2);
//extern variable
extern int num_vex;
extern int num_unk;
extern int num_edge;
extern int max_deg;
extern int **start_end;
extern int **edge;
extern int *degree_node;
extern int **simple_edge;
extern char name_final_result[256];
extern char file_G[256];
extern char file_E[256];
extern char file_D[256];
extern char file_S[256];
//Struct of the Solution

typedef struct struct_individual {
	int *permutation;
	int *permutationNew;
	int cbmp;
	int *wc;
	int *cbnodes;
} Struc_Sol;
Struc_Sol *CurrentS,*BestS, *PopS, *ChS1, *ChS2;

// Variables time
struct tms glo_start;
struct tms glo_end;
struct tms midt;
double clockTicksPerSecond;
double startTimeSeconds;
double cpuTimeTheBest;

//Local variable
int *CV;
int **Candidatemove;
int **mat_sum;
//declaration of function
int memSearch(int value_alb, int num_pop);
void setdataStructSol(int pop_sol);
void freedataStructSol(int pop_sol);
void EvaSol(Struc_Sol *Sol);
void IniSol(Struc_Sol *Sol);
void CopySolution(Struc_Sol *SouSol, Struc_Sol *DesSol);
int LS(Struc_Sol *C, long long ite);
double get_time(void);
void get_newwc(Struc_Sol *C, int *wc, int u, int v);
void get_swap(Struc_Sol *C, int node1, int node2);
void update_swap(Struc_Sol *C, int node1, int node2, int *newwc);
void update_cbnodes(Struc_Sol *C, int node1, int node2);
void inten_Search(Struc_Sol *C);
void generate_child(int num_pop,Struc_Sol *C1,Struc_Sol *C2, int length);
void quality_update_pop(Struc_Sol *C, int length);
int get_bestCB_pop(int value_alb, int num_pop);
double get_entropy_pop(int value_alb, int num_pop);
double get_averagedistance_pop(int value_alb, int num_pop);
// Function

int get_bestCB_pop(int value_alb, int num_pop){
	int i;
	int bestreturn=num_vex/2+1;

	for(i=0;i<num_pop;i++) if(PopS[i].cbmp<bestreturn) bestreturn=PopS[i].cbmp;

	return (bestreturn);
}
double get_entropy_pop(int value_alb, int num_pop){

	int i,j,k;

	double entropy=0.0;
	for(i=0;i<num_vex;i++) for(j=0;j<num_vex;j++) mat_sum[i][j]=0;
	for(i=0;i<num_vex;i++)
		for(k=0;k<num_pop;k++)
			mat_sum[i][PopS[k].permutation[i]]++;
	for(i=0;i<num_vex;i++) for(j=0;j<num_vex;j++) {
		if(mat_sum[i][j]==0) continue;
		entropy=entropy+1.0*mat_sum[i][j]/num_pop*log(mat_sum[i][j]*1.0/num_pop);
	}
	entropy=-entropy/(num_vex*log(num_vex));
//	cout<<entropy<<endl;
	return entropy;
}
double get_averagedistance_pop(int value_alb, int num_pop){
	int i,j;
	double sum_dis=0.0;
	double ave_dis;
//	cout<<num_vex<<endl;
	for(i=0;i<num_pop;i++) for(j=i+1;j<num_pop;j++){
		sum_dis=sum_dis+get_dis_tsp(PopS[i].permutationNew,PopS[j].permutationNew,num_vex, mat_sum);
	}

	ave_dis=2*sum_dis/((num_pop-1)*num_pop);

	return ave_dis;
}

void setdataStructSol(int pop_sol){
	//initial *CV,*CurrentS, *BestS, **Candidatemove,*PopS, *ChS1, *ChS2;
	int ub_cb=num_vex/2+1;

	int i;
	CV=(int*)get_vector(num_vex);
	for(i=0;i<num_vex;i++) CV[i]=-1;
	/*************************************/
	CurrentS=(Struc_Sol*)malloc(sizeof(Struc_Sol));
	if(CurrentS==NULL){
		cout<<"Memory error in CurrentS"<<endl;
		exit(-1);
	}
	CurrentS->permutation=(int*)get_vector(num_vex);
	CurrentS->permutationNew=(int*)get_vector(num_vex);
	CurrentS->wc=(int*)get_vector(ub_cb);
	CurrentS->cbnodes=(int*)get_vector(num_vex);
	/*************************************/
	BestS=(Struc_Sol*)malloc(sizeof(Struc_Sol));
	if(BestS==NULL){
		cout<<"Memory error in BestS"<<endl;
		exit(-1);
	}
	BestS->permutation=(int*)get_vector(num_vex);
	BestS->permutationNew=(int*)get_vector(num_vex);
	BestS->wc=(int*)get_vector(ub_cb);
	BestS->cbnodes=(int*)get_vector(num_vex);
	/*************************************/
	Candidatemove=(int**)get_matrix(num_vex*num_vex,2);
	/*************************************/
	PopS=(Struc_Sol*)malloc(pop_sol*sizeof(Struc_Sol));
	if(PopS==NULL){
		cout<<"Memory error in PopS"<<endl;
		exit(-1);
	}
	for(i=0;i<pop_sol;i++){
		PopS[i].permutation=(int*)get_vector(num_vex);
		PopS[i].permutationNew=(int*)get_vector(num_vex);
		PopS[i].wc=(int*)get_vector(ub_cb);
		PopS[i].cbnodes=(int*)get_vector(num_vex);
	}
	/*************************************/
	ChS1=(Struc_Sol*)malloc(sizeof(Struc_Sol));
	if(ChS1==NULL){
		cout<<"Memory error in ChS1"<<endl;
		exit(-1);
	}
	ChS1->permutation=(int*)get_vector(num_vex);
	ChS1->permutationNew=(int*)get_vector(num_vex);
	ChS1->wc=(int*)get_vector(ub_cb);
	ChS1->cbnodes=(int*)get_vector(num_vex);
	/*************************************/
	/*************************************/
	ChS2=(Struc_Sol*)malloc(sizeof(Struc_Sol));
	if(ChS2==NULL){
		cout<<"Memory error in CurrentS"<<endl;
		exit(-1);
	}
	ChS2->permutation=(int*)get_vector(num_vex);
	ChS2->permutationNew=(int*)get_vector(num_vex);
	ChS2->wc=(int*)get_vector(ub_cb);
	ChS2->cbnodes=(int*)get_vector(num_vex);
	/**intitial the mat_s********************/
	mat_sum=(int**)get_matrix(num_vex,num_vex);

}

void freedataStructSol(int pop_sol){
	int i;
	//free the space of *CV, *CurrentS, *BestS, **Candidatemove
	free(CV);		CV=NULL;

	free(CurrentS->permutation);		CurrentS->permutation=NULL;
	free(CurrentS->permutationNew);		CurrentS->permutationNew=NULL;
	free(CurrentS->wc);					CurrentS->wc=NULL;
	free(CurrentS->cbnodes);			CurrentS->cbnodes=NULL;
	free(CurrentS);						CurrentS=NULL;

	free(BestS->permutation);			BestS->permutation=NULL;
	free(BestS->permutationNew);		BestS->permutationNew=NULL;
	free(BestS->wc);					BestS->wc=NULL;
	free(BestS->cbnodes);				BestS->cbnodes=NULL;
	free(BestS);						BestS=NULL;

	for(i=0;i<num_vex;i++){
		free(Candidatemove[i]);
		Candidatemove[i]=NULL;
	}
	free(Candidatemove);				Candidatemove=NULL;

	for(i=0;i<pop_sol;i++){
		free(PopS[i].permutation);			PopS[i].permutation=NULL;
		free(PopS[i].permutationNew);		PopS[i].permutationNew=NULL;
		free(PopS[i].wc);					PopS[i].wc=NULL;
		free(PopS[i].cbnodes);				PopS[i].cbnodes=NULL;
	}
	free(PopS);							PopS=NULL;

	free(ChS1->permutation);			ChS1->permutation=NULL;
	free(ChS1->permutationNew);			ChS1->permutationNew=NULL;
	free(ChS1->wc);						ChS1->wc=NULL;
	free(ChS1->cbnodes);				ChS1->cbnodes=NULL;
	free(ChS1);							ChS1=NULL;

	free(ChS2->permutation);			ChS2->permutation=NULL;
	free(ChS2->permutationNew);			ChS2->permutationNew=NULL;
	free(ChS2->wc);						ChS2->wc=NULL;
	free(ChS2->cbnodes);				ChS2->cbnodes=NULL;
	free(ChS2);							ChS2=NULL;

	for(i=0;i<num_vex;i++){
		free(mat_sum[i]);
		mat_sum[i]=NULL;
	}
	free(mat_sum);						mat_sum=NULL;

}

void IniSol(Struc_Sol *Sol){
	int i,temp;
	int expermu[num_vex];

	for (i=0;i<num_vex;i++) Sol->permutation[i]=i;
	for (i=0;i<num_vex;i++) expermu[i]=rand()%num_vex;
	for (i=0;i<num_vex;i++){
		temp=Sol->permutation[i];
		Sol->permutation[i]=Sol->permutation[expermu[i]];
		Sol->permutation[expermu[i]]=temp;
	}
	for (i=0;i<num_vex;i++) Sol->permutationNew[Sol->permutation[i]]=i;
	EvaSol(Sol);
}

void EvaSol(Struc_Sol *Sol){
	int ub_cb=num_vex/2+1;
	int i;
	int cd;

	Sol->cbmp=0;
	memset(Sol->wc,0,ub_cb*sizeof(int));
	memset(Sol->cbnodes,0,num_vex*sizeof(int));

	for (i=0;i<num_vex;i++) Sol->permutationNew[Sol->permutation[i]]=i;
	for (i=0;i<num_edge;i++){
		cd=get_CycD(Sol->permutation[simple_edge[i][0]],Sol->permutation[simple_edge[i][1]]);
		Sol->wc[cd]++;
		if(cd>Sol->cbnodes[simple_edge[i][0]]) Sol->cbnodes[simple_edge[i][0]]=cd;
		if(cd>Sol->cbnodes[simple_edge[i][1]]) Sol->cbnodes[simple_edge[i][1]]=cd;
		if(cd>Sol->cbmp) Sol->cbmp=cd;
	}
}



void CopySolution(Struc_Sol *SouSol, Struc_Sol *DesSol){
	memcpy(DesSol->permutation, SouSol->permutation, num_vex*sizeof(int));
	memcpy(DesSol->permutationNew, SouSol->permutationNew, num_vex*sizeof(int));
	memcpy(DesSol->wc, SouSol->wc, (num_vex/2+1)*sizeof(int));
	memcpy(DesSol->cbnodes, SouSol->cbnodes, num_vex*sizeof(int));
	DesSol->cbmp = SouSol->cbmp;
}

double get_time(void){
	double endTimeSeconds=0.0;
	times(&glo_end);
	endTimeSeconds = glo_end.tms_utime/clockTicksPerSecond;
	return (endTimeSeconds - startTimeSeconds);
}



void get_newwc(Struc_Sol *C, int *wc, int node1, int node2){
	int lab_oldnode1,lab_oldnode2;
	int lab_newnode1,lab_newnode2;
	int i;
	int ls,le;
	int cycd;

	for(i=0;i<num_vex/2+1;i++) wc[i]=C->wc[i];
	lab_newnode2=lab_oldnode1=C->permutation[node1];
	lab_newnode1=lab_oldnode2=C->permutation[node2];
	/**************************************************************/
	ls=start_end[node1][0];
	le=start_end[node1][1];

	for (i=ls;i<le+1;i++){
		if (lab_newnode1!=C->permutation[edge[i][1]]){
			/*delete the old cycd*/
			cycd=get_CycD(lab_oldnode1, C->permutation[edge[i][1]]);
			wc[cycd]--;
			/*get the new cycd*/
			cycd=get_CycD(lab_newnode1, C->permutation[edge[i][1]]);
			wc[cycd]++;
		}
	}
	/**************************************************************/
	ls=start_end[node2][0];
	le=start_end[node2][1];

	for (i=ls;i<le+1;i++){
		if (lab_newnode2!=C->permutation[edge[i][1]]){
			/*delete the old cycd*/
			cycd=get_CycD(lab_oldnode2, C->permutation[edge[i][1]]);
			wc[cycd]--;
			/*get the new cycd*/
			cycd=get_CycD(lab_newnode2, C->permutation[edge[i][1]]);
			wc[cycd]++;
		}
	}
	/**************************************************************/
}

void get_swap(Struc_Sol *C, int node1, int node2){
	int temp_node;
	int i;
	int newwc[num_vex/2+1];
	//update the newwc and cbnodes
	update_swap(C, node1, node2, newwc);

	//update the wc and cbmp
	for (i=0;i<num_vex/2+1;i++){
		C->wc[i]=newwc[i];
	}
	for (i=num_vex/2;newwc[i]==0;i--);
	C->cbmp=i;
	//update the permutation and permutationNew
	temp_node=C->permutation[node1];
	C->permutation[node1]=C->permutation[node2];
	C->permutation[node2]=temp_node;

	C->permutationNew[C->permutation[node1]]=node1;
	C->permutationNew[C->permutation[node2]]=node2;

	update_cbnodes(C,node1,node2);
}

void update_swap(Struc_Sol *C, int node1, int node2,int *wc){
	int lab_oldnode1,lab_oldnode2;
	int lab_newnode1,lab_newnode2;
	int i;
	int ls,le;
	int cycd;

	for(i=0;i<num_vex/2+1;i++) wc[i]=C->wc[i];
	lab_newnode2=lab_oldnode1=C->permutation[node1];
	lab_newnode1=lab_oldnode2=C->permutation[node2];
	/**************************************************************/
	ls=start_end[node1][0];
	le=start_end[node1][1];

	for (i=ls;i<le+1;i++){
		if (lab_newnode1!=C->permutation[edge[i][1]]){
			/*delete the old cycd*/
			cycd=get_CycD(lab_oldnode1, C->permutation[edge[i][1]]);
			wc[cycd]--;
			/*get the new cycd*/
			cycd=get_CycD(lab_newnode1, C->permutation[edge[i][1]]);
			wc[cycd]++;
		}
	}
	/**************************************************************/
	ls=start_end[node2][0];
	le=start_end[node2][1];

	for (i=ls;i<le+1;i++){
		if (lab_newnode2!=C->permutation[edge[i][1]]){
			/*delete the old cycd*/
			cycd=get_CycD(lab_oldnode2, C->permutation[edge[i][1]]);
			wc[cycd]--;
			/*get the new cycd*/
			cycd=get_CycD(lab_newnode2, C->permutation[edge[i][1]]);
			wc[cycd]++;
		}
	}
	/**************************************************************/
}

void update_cbnodes(Struc_Sol *C, int node1, int node2){
	int i,j,cbNi,cbNj,vl,ul,temp_node;

	cbNi=0;

	for(i=start_end[node1][0];i<=start_end[node1][1];i++){
		cbNj=0;
		for(j=start_end[edge[i][1]][0];j<=start_end[edge[i][1]][1];j++){
			vl=get_CycD(C->permutation[edge[i][1]], C->permutation[edge[j][1]]);
			if(vl>cbNj) cbNj=vl;
		}
		C->cbnodes[edge[i][1]]=cbNj;
		vl=get_CycD(C->permutation[node1], C->permutation[edge[i][1]]);
		if(vl>cbNi) cbNi=vl;
	}
	C->cbnodes[node1]=cbNi;

	cbNi=0;
	for(i=start_end[node2][0];i<=start_end[node2][1];i++){
		cbNj=0;
		for(j=start_end[edge[i][1]][0];j<=start_end[edge[i][1]][1];j++){
			ul=get_CycD(C->permutation[edge[i][1]], C->permutation[edge[j][1]]);
			if(ul>cbNj) cbNj=ul;
		}
		C->cbnodes[edge[i][1]]=cbNj;
		ul=get_CycD(C->permutation[node2], C->permutation[edge[i][1]]);
		if(ul>cbNi) cbNi=ul;
	}
	C->cbnodes[node2]=cbNi;

}

int LS(Struc_Sol *C, long long ite){
	int i,j,k;
	int cbmax=C->cbmp;
	int cbth=C->cbmp;
	int u,v;
	int num_wc=num_vex/2+1;
	int oldwc[num_wc];
	int newwc[num_wc];
	int bestwc[num_wc];
	int judge_f=0;
	int num_cand=0;
	int	index_best=0;


	memset(CV,0,num_vex*sizeof(int));
	for (i=0;i<num_wc;i++) oldwc[i]=C->wc[i];
	for (i=0;i<num_vex;i++) if(C->cbnodes[i]>=cbth) CV[i]=1;

	for (i=0;i<num_vex;i++) {
		if(CV[i]==0) continue;
		u=i;
		for(j=0;j<num_vex;j++){
			if(j==u) continue;
			v=j;
			get_newwc(C,newwc,u,v);
//			judge_f=judge_one(oldwc,newwc,num_wc);
			judge_f=judge_or(oldwc,newwc,num_wc);
			if(judge_f<0){
				Candidatemove[0][0]=u;
				Candidatemove[0][1]=v;
				num_cand=1;
				for (k=0;k<num_wc;k++) oldwc[k]=newwc[k];
			}
		}
	}
	if(num_cand) get_swap(C,Candidatemove[0][0],Candidatemove[0][1]);
	return num_cand;
}

void inten_Search(Struc_Sol *C){
	int flag_up=1;
	static long long ite=0;
	double time_total=0;

	while(flag_up){
		if(LS(C,ite)) ite++;
		else flag_up=0;

		time_total=get_time();
		if(C->cbmp<BestS->cbmp){
			CopySolution(C,BestS);
			cout<<"ite="<<ite<<" "<<BestS->cbmp<<endl;
			ofstream caout(name_final_result,ios::out|ios::app);
			if (caout.is_open()){
				caout<<ite<<" ";
				caout<<BestS->cbmp<<" ";
				caout<<time_total<<" ";
				caout<<endl;
				caout.close();
			}
		}
//		if(time_total>600) break;
	}

}

void generate_child(int num_pop, Struc_Sol *C1, Struc_Sol *C2, int length){
	int idx_f,idx_m;

	do{
		idx_f=rand()%num_pop;
		idx_m=rand()%num_pop;
	}while(idx_f==idx_m);
//	PMX(PopS[idx_f].permutation, PopS[idx_m].permutation, length, C1->permutation , C2->permutation);
//	CX(PopS[idx_f].permutationNew, PopS[idx_m].permutationNew, length, C1->permutationNew , C2->permutationNew);
//	OX(PopS[idx_f].permutation, PopS[idx_m].permutation, length, C1->permutation , C2->permutation);
//	OX2(PopS[idx_f].permutation, PopS[idx_m].permutation, length, C1->permutation , C2->permutation);
	while(DPX(PopS[idx_f].permutationNew, PopS[idx_m].permutationNew, length, C1->permutationNew , C2->permutationNew, mat_sum)==0){};
//	cout<<"enter the generate children"<<endl;
//	int veri[length];
//	for(int i=0;i<length;i++) veri[i]=0;
//	for(int i=0;i<length;i++) {
//		if(veri[C1->permutation[i]]==1) {
//			cout<<"error in X"<<endl;
//			exit(-1);
//		}
//		veri[C1->permutation[i]]=1;
//	}
//	for(int i=0;i<length;i++) cout<<C1->permutation[i]<<" ";
//	cout<<endl;
//	cout<<"length="<<length<<endl;
	for (int i=0;i<num_vex;i++) C1->permutation[C1->permutationNew[i]]=i;
	EvaSol(C1);
//	cout<<"finish the generate children"<<endl;
//	EvaSol(C2);
}

void quality_update_pop(Struc_Sol *C,int length){
	int i;
	int index_pire=-1, max_dif=0;
	Struc_Sol *temp;

	for(i=0;i<length;i++) if(PopS[i].cbmp-C->cbmp>max_dif) {
		index_pire=i;
		max_dif=PopS[i].cbmp-C->cbmp;
	}
	if(index_pire!=-1){
		if(max_dif>=0){
			temp=&PopS[index_pire];
			CopySolution(C,temp);
		}
	}
}

int memSearch(int value_alb, int num_pop){
	int i;
	int bc_final;
	double time_total=0;
	Struc_Sol *temp;
	int count=0;

    clockTicksPerSecond = (double)sysconf(_SC_CLK_TCK); // Get the clock ticks per second
	setdataStructSol(num_pop);
//	IniSol(CurrentS);
//	CopySolution(CurrentS, BestS);
	BestS->cbmp=num_vex/2+1;
	times(&glo_start);
	startTimeSeconds = glo_start.tms_utime/clockTicksPerSecond;


	// Initialise the population
	for(i=0;i<num_pop;i++){
		temp=&PopS[i];
		IniSol(temp);
		inten_Search(temp);
//		cout<<"finish the "<<i<< " initialisation"<<endl;
	}
	cout<<"finish all initial"<<endl;
	//LS and memetic
	while(BestS->cbmp>value_alb){
		generate_child(num_pop, ChS1, ChS2, num_vex);
		if(count%200==0){
			ofstream Saout(file_S,ios::out|ios::app);
			if (Saout.is_open()){
				Saout<<ChS1->cbmp<<endl;
				Saout.close();
			}

		}

		inten_Search(ChS1);							//Intensification
		quality_update_pop(ChS1, num_pop);

//		cout<<"generation="<<count<<endl;

		if(count%200==0){
			ofstream Gaout(file_G,ios::out|ios::app);
			if (Gaout.is_open()){

				Gaout<<get_bestCB_pop(value_alb,num_pop)<<endl;
				Gaout.close();
			}

			ofstream Daout(file_D,ios::out|ios::app);
			if (Daout.is_open()){
				Daout<<get_averagedistance_pop(value_alb,num_pop)<<endl;
				Daout.close();
			}

			ofstream Eaout(file_E,ios::out|ios::app);
			if (Eaout.is_open()){
				Eaout<<get_entropy_pop(value_alb,num_pop)<<endl;
				Eaout.close();
			}
		}
		count++;
		if(count>=20000) break;

		time_total=get_time();
//		if(time_total>=600) break;
	}
//	for(i=0;i<num_pop;i++) cout<<i<<" "<<PopS[i].cbmp<<endl;
	//verify the result and return the final cbmp
	EvaSol(BestS);
	bc_final=BestS->cbmp;
	freedataStructSol(num_pop);
	return bc_final;
}

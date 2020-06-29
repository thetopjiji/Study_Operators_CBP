/*
 * BCP_X.cpp
 *
 *  Created on: 7 mars 2019
 *      Author: ren
 */

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
using namespace std;
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))


//Declaration of function
int judge_one(int *oldwc, int *newwc, int num_wc);
int judge_or(int *newwc, int *otherwc,int num_wc);
void PMX(int *Father, int *Mother, int length, int *c1, int *c2);
void CX(int *Father, int *Mother, int length, int *c1, int *c2);
void OX(int *Father, int *Mother, int length, int *c1, int *c2);
void OX2(int *Father, int *Mother, int length, int *c1, int *c2);
int DPX(int *Father, int *Mother, int length, int *c1, int *c2, int **mat_s);
int get_dis_tsp(int *Father, int *Mother, int length, int **matrix_sum);


//local function
int func_2connection(int **mat_s, int *temp_0, int *temp_2, int *c1, int *st_pos, int &num_set, int length, int pt_set);
int func_1connection(int **mat_s, int *temp_0, int *temp_2, int *c1, int *st_pos, int &num_set, int length, int pt_set, int num_0);
int func_0connection(int **mat_s, int *temp_0, int *temp_2, int *c1, int *st_pos, int &num_set, int length, int pt_set, int num_0);
int get_newpos(int nowpos,int length, int *c1, int *st_pos);

//The function
int judge_or(int *otherwc,int *newwc, int num_wc){
	int o;
	int cb_new;
	int cb_old;

	for(o=num_wc-1;newwc[o]==0;o--);
	cb_new=o;

	for(o=num_wc-1;otherwc[o]==0;o--);
	cb_old=o;

	return (cb_new-cb_old);
}


int judge_one(int *oldwc, int *newwc, int num_wc){

	//If return value <0, means get improving; return value=0, equal solution; get worse otherwise
	int o,m;
	int h;
	int key_c;
	int deltawc[num_wc];


	for (o=0;o<num_wc;o++) deltawc[o]=newwc[o]-oldwc[o];
	for (h=num_wc-1;oldwc[h]==0;h--);
	key_c=h;
	for (m=num_wc-1;deltawc[m]==0 && m>=key_c;m--);
	return deltawc[m];
}

void PMX(int *Father, int *Mother, int length, int *c1, int *c2){
	//partialle mapped crossover
	//Input the father permutation, mother permutation and the length of the permutation
	//Return child permutation 1 and child permutation 2
	int left,right,temp;
	int ll[length];
	int lf,lm;

	for(int i=0;i<length;i++) ll[i]=-1;

	left=rand()%length;
	right=rand()%length;
	//ensure left<right
	if(left>right){
		temp=left;
		right=left;
		right=temp;
	}
	//copy from parent to children
	for(int i=left;i<=right;i++){
		c1[i]=Father[i];
		c2[i]=Mother[i];
	}
	for(int i=0;i<left;i++){
		c1[i]=Mother[i];
		c2[i]=Father[i];
	}
	for(int i=right+1;i<length;i++){
		c1[i]=Mother[i];
		c2[i]=Father[i];
	}
	//fill the ll
	for(int i=left;i<=right;i++){
		lf=Father[i];
		lm=Mother[i];
		ll[lf]=lm;
	}

	//exchange to get a feasible solution
	for(int i=0;i<left;i++){
		temp=ll[c1[i]];
		if(temp==-1) continue;
		while(ll[temp]!=-1){
			temp=ll[temp];
		}
		c1[i]=temp;
	}
	for(int i=right+1;i<length;i++){
		temp=ll[c1[i]];
		if(temp==-1) continue;
		while(ll[temp]!=-1){
			temp=ll[temp];
		}
		c1[i]=temp;
	}

}

void CX(int *Father, int *Mother, int length, int *c1, int *c2){
	//Cycle crossover
	//Input the father permutation, mother permutation and the length of the permutation
	//Return child permutation 1 and child permutation 2
	int start_p=rand()%length;
	int permu_done[length];
	int index;
	bool done;

	for(int i=0;i<length;i++) permu_done[i]=0;
	index=start_p;

	c1[index]=Father[index];
	c2[index]=Mother[index];
	while(permu_done[index]==0){
		permu_done[index]=1;
		for(int i=0;i<length;i++){
			if((done=(Father[index]==Mother[i]))){
				c1[i]=Father[i];
				c2[i]=Mother[i];
				index=i;
			}
			if(done) break;
		}
	}
	for(int i=0;i<length;i++) if(permu_done[i]==0) {
		c1[i]=Mother[i];
		c2[i]=Father[i];
	}
}


void OX(int *Father, int *Mother, int length, int *c1, int *c2){
	//Order crossover
	//Input the father permutation, mother permutation and the length of the permutation
	//Return child permutation 1 and child permutation 2
	int left, right, temp;
	int temp_per[length];

	left=rand()%length;
	right=rand()%length;
	//ensure left<right
	if(left>right){
		temp=left;
		right=left;
		right=temp;
	}
	//copy from parent to children
	for(int i=left;i<right+1;i++){
		c1[i]=Father[i];
		c2[i]=Mother[i];
	}
	//exchange and turn around for c1
	int index=0;
	int k=left;
	//generate for new permutation
	for(int i=right+1;i<length;i++) temp_per[i-right-1]=Mother[i];
	for(int i=0;i<right+1;i++)		temp_per[i+length-right-1]=Mother[i];
	//tour
	for(int i=right+1;i<length;i++){
		 k=left;
		while(k<right+1){
			if(temp_per[index]==c1[k]) {
				index++;
				k=left;
			}
			else k++;
		}
		c1[i]=temp_per[index];
		index++;
	}
	for(int i=0;i<left;i++){
		 k=left;
		while(k<right+1){
			if(temp_per[index]==c1[k]) {
				index++;
				k=left;
			}
			else k++;
		}
		c1[i]=temp_per[index];
		index++;
	}
	//exchange and turn around for c2
		index=0;
		for(int i=right+1;i<length;i++) temp_per[i-right-1]=Father[i];
		for(int i=0;i<right+1;i++)		temp_per[i+length-right-1]=Father[i];

		for(int i=right+1;i<length;i++){
			 k=left;
			while(k<right+1){
				if(temp_per[index]==c2[k]) {
					index++;
					k=left;
				}
				else k++;
			}
			c2[i]=temp_per[index];
			index++;
		}
		for(int i=0;i<left;i++){
			k=left;
			while(k<right+1){
				if(temp_per[index]==c2[k]) {
					index++;
					k=left;
				}
				else k++;
			}
			c2[i]=temp_per[index];
			index++;
		}
}


void OX2(int *Father, int *Mother, int length, int *c1, int *c2){
	//Order based crossover
	//Input the father permutation, mother permutation and the length of the permutation
	//Return child permutation 1 and child permutation 2
	int i,index,j,count=0;
	int num_sel=rand()%length;	//num_sel from 0 to length-1
	int pos[length];
	int pos2[length];
	int seq[length];
	int seq2[length];

	for(i=0;i<length;i++) pos[i]=0;
	for(i=0;i<length;i++) pos2[i]=0;
	for(i=0;i<num_sel;i++){
		do{
			index=rand()%length;
		}while(pos[index]==1);
		pos[index]=1;
	}
	count=0;
	for(i=0;i<length;i++) if(pos[i]) {
		seq[count]=i;
		count++;
	}
	//copy Father to c1, copy Mother to c2
	for(i=0;i<length;i++) c1[i]=Father[i];
	for(i=0;i<length;i++) c2[i]=Mother[i];
	//find the corresponding position for c2 from Mother
	for(i=0;i<num_sel;i++){
		for(j=0;j<length;j++) if(Mother[j]==Father[seq[i]]){
			pos2[j]=1;
			break;
		}
	}
	for(i=0;i<num_sel;i++){
		for(j=0;j<length;j++) if(pos2[j]){
			c2[j]=Father[seq[i]];
			pos2[j]=0;
			break;
		}
	}
	//find the corresponding position for c1 from Father
	for(i=0;i<num_sel;i++){
		for(j=0;j<length;j++) if(Father[j]==Mother[seq[i]]){
			pos2[j]=1;
			break;
		}
	}
	for(i=0;i<num_sel;i++){
		for(j=0;j<length;j++) if(pos2[j]){
			c1[j]=Mother[seq[i]];
			pos2[j]=0;
			break;
		}
	}
}

int DPX(int *Father, int *Mother, int length, int *c1, int *c2,int **mat_s){
	//Distance Preserving crossover
	//Input the father permutation, mother permutation and the length of the permutation
	//Return child permutation 1 and child permutation 2
	//use a matrix "mat_s" to realize

	int i,j;

	for(i=0;i<length;i++) for(j=0;j<length;j++)	mat_s[i][j]=0;
	//get the all connections of adjacent nodes
	for(i=0;i<length-1;i++) {
		mat_s[Father[i]][Father[i+1]]++;
		mat_s[Father[i+1]][Father[i]]++;
	}
	mat_s[Father[length-1]][Father[0]]++;
	mat_s[Father[0]][Father[length-1]]++;

	for(i=0;i<length-1;i++) {
		mat_s[Mother[i]][Mother[i+1]]++;
		mat_s[Mother[i+1]][Mother[i]]++;
	}
	mat_s[Mother[length-1]][Mother[0]]++;
	mat_s[Mother[0]][Mother[length-1]]++;



	//if there is no 1 in some column or rang, fill it with -1
	bool done=true;
	for(i=0;i<length;i++){
		done=true;
		for(j=0;j<length;j++) if(mat_s[i][j]==1) {done=false;break;}
		if(done) for(int k=0;k<length;k++) if(mat_s[i][k]==0) mat_s[i][k]=-1;
		if(done) for(int k=0;k<length;k++) if(mat_s[k][i]==0) mat_s[k][i]=-1;
	}

	//fill the diagram with -1
	for(i=0;i<length;i++) mat_s[i][i]=-1;


	//Regroup the tour
	int st_pos[length];
	int temp_0[length];
	int temp_2[length];
	int num_0=0;
	int	pt_set=0;
	int num_set=0;
	int num_1;
	int num_2;


	for(i=0;i<length;i++) st_pos[i]=-1;
	for(i=0;i<length;i++) c1[i]=-1;

	pt_set=rand()%length;
	c1[0]=pt_set;
	st_pos[pt_set]=0;						//pt_set is the number of point, st_pos is its position

	//first point cannot be chosen
	for(j=0;j<length;j++) if(mat_s[j][pt_set]==0) mat_s[j][pt_set]=-1;

	num_set=1;
	while(num_set<length){
		num_0=0;
		num_1=0;
		num_2=0;

//		cout<<pt_set<<" "<<st_pos[pt_set]<<endl;
//		sortmat_s(mat_s,length);
		// to judge the state of pt_set, one commune edge or two
		// get the candidate connecting node
		for(i=0;i<length;i++) {
			if(mat_s[pt_set][i]==1) num_1++;
			if(mat_s[pt_set][i]==2) {
				temp_2[num_2]=i;
				num_2++;
			}
			if(mat_s[pt_set][i]==0) {
				temp_0[num_0]=i;
				num_0++;
			}
		}
		// num_2=2 has two commune edge
		if(num_2==2){
			pt_set=func_2connection(mat_s, temp_0,temp_2,c1,st_pos,num_set,length, pt_set);
			if(pt_set==-1) {
				return 0;
//				cout<<"pt_set =-1 "<<endl;
//				exit(-1);
			}
		}
		// num_2=1 has two commune edge
		if(num_2==1) pt_set=func_1connection(mat_s, temp_0,temp_2,c1,st_pos,num_set,length, pt_set, num_0);
		if(pt_set==-1) return 0;

		// num_2=0 has two commune edge
		if(num_2==0) pt_set=func_0connection(mat_s, temp_0,temp_2,c1,st_pos,num_set,length, pt_set, num_0);
		if(pt_set==-1) return 0;

	}
	if(mat_s[pt_set][c1[0]]==1) return 0;

	return 1;


}

int get_newpos(int nowpos,int length, int *c1, int *st_pos){
	int newpos;

	if(nowpos+1>length-1)	newpos=0;
	else newpos=nowpos+1;

	if(c1[newpos]==-1)	return newpos;	//-1: no vertex set here

	if(nowpos-1<0)	newpos=length-1;
	else newpos=nowpos-1;

	if(c1[newpos]==-1)	return newpos;	//-1: no vertex set here

	cout<<"error in 2con, the 2 voisin are already set"<<endl;
	exit(-1);
}

int func_2connection(int **mat_s, int *temp_0, int *temp_2, int *c1, int *st_pos, int &num_set, int length, int pt_set){
	int tvex;
	int i,j;
	int nowpos=st_pos[pt_set];
	int pt_return=-1;
	int newpos;

	for(i=0;i<2;i++){
		tvex=temp_2[i];
		//decide the connection node
		if(st_pos[tvex]!=-1) continue;	//already set
		//get the node to set and decide the pos to set
		newpos=get_newpos( nowpos, length, c1, st_pos);
		// set the tvex
		c1[newpos]=tvex;
		st_pos[tvex]=newpos;
		pt_return=tvex;
		num_set++;
		// update the mat_s
		for(j=0;j<length;j++) if(mat_s[j][tvex]==0) mat_s[j][tvex]=-1;
	}

	return pt_return;
}
int func_1connection(int **mat_s, int *temp_0, int *temp_2, int *c1, int *st_pos, int &num_set, int length, int pt_set, int num_0){
	int tvex;
	int j;
	int nowpos=st_pos[pt_set];
	int pt_return=-1;
	int newpos;
	/* its adjacent vertex is not set*/
	tvex=temp_2[0];
	if(st_pos[tvex]==-1){
		//get the node to set and decide the pos to set
		newpos=get_newpos( nowpos, length, c1, st_pos);
		// set the tvex
		c1[newpos]=tvex;
		st_pos[tvex]=newpos;
		pt_return=tvex;
		num_set++;
		// update the mat_s
//		for(j=0;j<length;j++) if(mat_s[tvex][j]==0) mat_s[tvex][j]=-1;
		for(j=0;j<length;j++) if(mat_s[j][tvex]==0) mat_s[j][tvex]=-1;
		return pt_return;
	}
	/* its adjacent vertex is already set*/
	//choose random one from temp_0[];
	if(num_0==0) return (-1);
	tvex=temp_0[rand()%num_0];

//	cout<<tvex<<endl;
	if(st_pos[tvex]!=-1) {
		cout<<"error in mat_s update in 1con"<<endl;
		exit(-1);
	}
	newpos=get_newpos( nowpos, length, c1, st_pos);
	// set the tvex
	c1[newpos]=tvex;
	st_pos[tvex]=newpos;
	pt_return=tvex;
	num_set++;
	mat_s[pt_set][tvex]=2;
	mat_s[tvex][pt_set]=2;
	// update the mat_s only if num_2=2
//	for(j=0;j<length;j++) if(mat_s[pt_set][j]==0) mat_s[pt_set][j]=-1;
	for(j=0;j<length;j++) if(mat_s[j][tvex]==0) mat_s[j][tvex]=-1;
	return pt_return;

}
int func_0connection(int **mat_s, int *temp_0, int *temp_2, int *c1, int *st_pos, int &num_set, int length, int pt_set, int num_0){
	int tvex;
	int i,j;
	int nowpos=st_pos[pt_set];
	int pt_return=-1;
	int newpos;
	/*it has no adjacent vertex*/
	if(num_0==0) return (-1);
	tvex=temp_0[rand()%num_0];
//	cout<<tvex<<endl;
	if(st_pos[tvex]!=-1) {
		cout<<"error in mat_s update in 0con"<<endl;
		exit(-1);
	}
	newpos=get_newpos(nowpos, length, c1, st_pos);
	// set the tvex
	c1[newpos]=tvex;
	st_pos[tvex]=newpos;
	pt_return=tvex;
	num_set++;
	mat_s[pt_set][tvex]=2;
	mat_s[tvex][pt_set]=2;
	// update the mat_s
//	for(j=0;j<length;j++) if(mat_s[tvex][j]==0) mat_s[tvex][j]=-1;
	for(j=0;j<length;j++) if(mat_s[j][tvex]==0) mat_s[j][tvex]=-1;

	return pt_return;
}

int get_dis_tsp(int *Father, int *Mother, int length, int **matrix_sum){

	int dis=0;

//	int matrix_sum[length][length];

	int i,j;
	int num_2=0;

	for(i=0;i<length;i++) {
		for(j=0;j<length;j++)	{
			matrix_sum[i][j]=0;
		}
	}
	//get the all connections of adjacent nodes
	for(i=0;i<length-1;i++) {
		matrix_sum[Father[i]][Father[i+1]]++;
		matrix_sum[Father[i+1]][Father[i]]++;
	}
	matrix_sum[Father[length-1]][Father[0]]++;
	matrix_sum[Father[0]][Father[length-1]]++;

	for(i=0;i<length-1;i++) {
		matrix_sum[Mother[i]][Mother[i+1]]++;
		matrix_sum[Mother[i+1]][Mother[i]]++;
	}
	matrix_sum[Mother[length-1]][Mother[0]]++;
	matrix_sum[Mother[0]][Mother[length-1]]++;
	// get the number of element whose value is 2 in mat_s, then n-this number
	for(i=0;i<length;i++) for(j=0;j<length;j++) if(matrix_sum[i][j]==2) num_2++;

//	cout<<"num_2="<<num_2<<endl;

	if(num_2%2!=0){
		cout<<"error in mat_s"<<endl;
		for(i=0;i<length;i++){
			for(j=0;j<length;j++) cout<<matrix_sum[i][j];

			cout<<endl;

		}

		exit(-1);
	}
	num_2=num_2/2;

	return (length-num_2);
}


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include "mt19937ar2.h"
#include "MemoryOperation.h"

using namespace std;

//** [PATH] **//
// location file (input)
#define LOC_FILE "../data/Foursquare_TIST15_MH/all_records.csv"
// census file (input)
#define CEN_FILE "../data/USCensus1990/all_records.csv"
// prefix of result file (output)
#define RES_FILE "../results/Res_Lap"

//** [PARAMETER] **//
// #runs
#define RUN_NUM 100
// normalize counts to probabilities (1:yes, 0:no)
//#define NORMALIZE 1
#define NORMALIZE 0

// #vertical regions
#define V_NUM 25
// #horizontal regions
#define H_NUM 25
// minimum of latitude
#define V_MIN 40.7
// maximum of latitude
#define V_MAX 40.8
// minimum of longitude
#define H_MIN -74.04
// maximum of longitude
#define H_MAX -73.92

#define MAX_CHAR_NUM		100000

// #categories
int CatNum;
// maximum #events
int MaxEventNum;

// event
int *Event;
// #users
int UserNum;
// event distribution
double *EventDist;
// noisy event distribution
double *NEventDist;
// estimated event distribution
double *EstEventDist;

FILE *FileOpen(string filename, const char *mode) {
	FILE *fp;

	if ((fp = fopen(filename.c_str(), mode)) == NULL) {
		cout << "cannot open " << filename << endl;
		exit(-1);
	}
	return fp;
}

// Read location file
void ReadLocFile(){
	double bou_v[V_NUM], bou_h[H_NUM];
	double v, h;
	int v_no, h_no, r_no;
	char *poi;
	char s[MAX_CHAR_NUM];
	FILE *fp;
	char *tok;
	int i, j;

	UserNum = 0;

	// calculate the boundary of regions --> bou_v, bou_h
	for (i = 0; i < V_NUM; i++){
		bou_v[i] = V_MIN + (V_MAX - V_MIN) * (double)i / (double)V_NUM;
	}
	for (i = 0; i < H_NUM; i++){
		bou_h[i] = H_MIN + (H_MAX - H_MIN) * (double)i / (double)H_NUM;
	}
//	for(i=0;i<V_NUM;i++) cout << bou_v[i] << endl;
//	for(i=0;i<H_NUM;i++) cout << bou_h[i] << endl;

	fp = FileOpen(LOC_FILE, "r");

	// read location file
	i = 0;
	while (fgets(s, MAX_CHAR_NUM, fp) != NULL){
		i++;

		// longitude --> h
		if ((tok = strtok(s, ",")) == NULL){
			printf("Error: File format is incorrect.\n");
			exit(-1);
		}
		for (j = 0; j < 2; j++){
			if ((tok = strtok(NULL, ",")) == NULL){
				printf("Error: File format is incorrect.\n");
				exit(-1);
			}
		}
		h = atof(tok);

		// latitude --> v
		if ((tok = strtok(NULL, ",")) == NULL){
			printf("Error: File format is incorrect.\n");
			exit(-1);
		}
		v = atof(tok);

		// latitude ID --> h_no
		if (h < H_MIN || h > H_MAX){
			printf("Warning: There is a location out of range (line:%d, longitude:%f, latitude:%f). We skip this location.\n", i, h, v);
			continue;
		}
		for (j = 0; j < H_NUM - 1; j++){
			if (h >= bou_h[j] && h < bou_h[j + 1]) break;
		}
		h_no = j;

		// longitude ID --> v_no
		if (v < V_MIN || v > V_MAX){
			printf("Warning: There is a location out of range (line:%d, longitude:%f, latitude:%f). We skip this location.\n", i, h, v);
			continue;
		}
		for (j = 0; j < V_NUM - 1; j++){
			if (v >= bou_v[j] && v < bou_v[j + 1]) break;
		}
		v_no = j;

		// region ID --> r_no
		r_no = v_no * H_NUM + h_no;
//		cout << h << " " << v << " " << h_no << " " << v_no << " " << r_no << endl;

		if (UserNum == MaxEventNum){
			printf("Error: The total number of events exceeds the maximum limit (=%d).\n", MaxEventNum);
			exit(-1);
		}

		// original region ID --> Event[UserNum]
		Event[UserNum] = v_no * H_NUM + h_no;
		UserNum++;
	}

	fclose(fp);

	printf("Total number of events: %d\n", UserNum);
}

// Read census file
void ReadCenFile() {
	char s[MAX_CHAR_NUM];
	FILE *fp;
	char *tok;
	int user_num;
	int i;

	UserNum = 0;

	fp = FileOpen(CEN_FILE, "r");

	// read category file
	while (fgets(s, MAX_CHAR_NUM, fp) != NULL) {
		if ((tok = strtok(s, ",")) == NULL) {
			printf("Error: File format is incorrect.\n");
			exit(-1);
		}

		if ((tok = strtok(NULL, ",")) == NULL) {
			printf("Error: File format is incorrect.\n");
			exit(-1);
		}

		// original category ID --> Event[UserNum]
		Event[UserNum] = atoi(tok);
		UserNum++;
	}

	fclose(fp);

	printf("Total number of events: %d\n", UserNum);
}

// Calculate (absolute frequency) distribution from events
void CalcDist(int *ev, int ev_num, double *dist){
	int i;

	// initialization
	for (i = 0; i < CatNum; i++) dist[i] = 0.0;

	// increase the population by one for each event
	for (i = 0; i < ev_num; i++) dist[ev[i]] += 1.0;

	// normalize counts to probabilities
	if(NORMALIZE){
		for (i = 0; i < CatNum; i++) dist[i] /= (double)ev_num;
	}
}

// Calculate noisy distribution from events by adding Lap
void CalcDistLap(int* ev, int ev_num, double* dist, double eps, int obf_mechanism) {
	int i, j;

	int ev_now;
	double exp_rnd1, exp_rnd2, lap_rnd;

	// initialization
	for (i = 0; i < CatNum; i++) dist[i] = 0.0;

	// increase the population by one for each event
	for (i = 0; i < ev_num; i++) {
		// current event --> ev_now
		ev_now = ev[i];

		dist[ev_now] += 1.0;
	}

	// add Lap (obf_mechanism times)
	if (eps != -1) {
//		double scale = 1.0 / eps;	// sensitivity = 1 in unbounded DP
		double scale = 2.0 / eps;	// sensitivity = 2 in bounded DP
		std::mt19937 rand_src(12345);
		std::exponential_distribution<> exp_dist(1.0 / scale);

		for (i = 0; i < CatNum; i++) {
			for (j = 0; j < obf_mechanism; j++) {
				// generate Lap(scale) by the difference between two exponential random valuables Exp(lambda), where lambda = 1/scale
				// c.f. https://www.johndcook.com/blog/2018/03/13/generating-laplace-random-variables/
				exp_rnd1 = exp_dist(rand_src);
				exp_rnd2 = exp_dist(rand_src);

				lap_rnd = exp_rnd1 - exp_rnd2;
				dist[i] += lap_rnd;
//				cout << i << " " << j << " " << lap_rnd << " " << dist[i] << endl;
			}
		}
	}

	// normalize counts to probabilities
	if(NORMALIZE){
		for (i = 0; i < CatNum; i++) dist[i] /= (double)ev_num;
	}
}

// Calculate MAE between two distributions
double MAE(double *dist1, double *dist2){
	double mae = 0.0;
	double l1_loss;
	int i;

	for (i = 0; i < CatNum; i++){
		l1_loss = fabs(dist1[i] - dist2[i]);
		mae += l1_loss;
	}
	mae /= (double)CatNum;

	return mae;
}

// Calculate MSE between two distributions
double MSE(double *dist1, double *dist2){
	double mse = 0.0;
	double l2_loss;
	int i;

	for (i = 0; i < CatNum; i++){
		l2_loss = (dist1[i] - dist2[i]) * (dist1[i] - dist2[i]);
		mse += l2_loss;
	}
	mse /= (double)CatNum;

	return mse;
}

int main(int argc, char *argv[])
{
	// distribution estimation method
	// (0:Empirical Estimation (w/o normalization), 1::Empirical Estimation (w/ normalization))
	int est_method_num = 2;
	int est_method[est_method_num] = {0,1};
	// parameter epsilon (-1: do not add noise)
	int eps_num = 24;
	double eps[eps_num] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,3,4,5,6,7,8,9,10,-1};

	int data_type;
	int obf_mechanism;
	char resfile[1025];
	int ep, em, rn;
	double mae_est_avg[est_method_num];
	double mse_est_avg[est_method_num];
	double mae_est_std[est_method_num];
	double mse_est_std[est_method_num];
	int zero_cat_num;
	double est_dist_sum;
	double mae_est, mse_est;
	FILE *fp;
	double dtmp;
	int i;

	// Initialization of Mersennne Twister
	unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
	init_by_array(init, length);
//	init_genrand(1);

	if (argc < 3) {
	    printf("Usage: %s [DATA_TYPE] [OBF_MECHANISM]\n\n", argv[0]);
		printf("[DATA_TYPE]: Data type (1: Foursquare_TIST15_MH, 2: USCensus1990)\n");
		printf("[OBF_MECHANISM]: Obfuscation mechanism (1:Lap, k(>=2): Lap x k times)\n");
		return -1;
	}

	data_type = atoi(argv[1]);
	obf_mechanism = atoi(argv[2]);

	// Set CatNum and MaxEventNum
	switch(data_type){
		// Foursquare_TIST15_MH
		case 1:
			CatNum = V_NUM * H_NUM;
			MaxEventNum = 359054;
			break;
		// USCensus1990
		case 2:
			CatNum = 400;
			MaxEventNum = 2458285;
			break;
		default:
			printf("Incorrect [DATA_TYPE]\n");
			exit(-1);
	}

	// malloc
	malloc1D(&EventDist, CatNum);
	malloc1D(&NEventDist, CatNum);
	malloc1D(&EstEventDist, CatNum);
	malloc1D(&Event, MaxEventNum);

	// read input file --> Event, UserNum
	switch(data_type){
		// Foursquare_TIST15_MH
		case 1:
			ReadLocFile();
			break;
		// USCensus1990
		case 2:
			ReadCenFile();
			break;
	}
//	for(i=0;i<UserNum;i++) cout << i << " " << Event[i] << endl;

	// calculate true distribution from events --> EventDist
	CalcDist(Event, UserNum, EventDist);
//	for(i=0;i<CatNum;i++) cout << i << " " << EventDist[i] << endl;

	// output header to result file
	for (em = 0; em < est_method_num; em++){
		sprintf(resfile, "%s_dt%d_om%d_em%d.csv", RES_FILE, data_type, obf_mechanism, est_method[em]);
		fp = FileOpen(resfile, "w");
		fprintf(fp, "epsilon,MAE_avg,MSE_avg,MAE_std,MSE_std\n");
		fclose(fp);
	}

	// for each parameter epsilon
	for (ep = 0; ep < eps_num; ep++){
		// initialization
		for (em = 0; em < est_method_num; em++) {
			mae_est_avg[em] = mse_est_avg[em] = 0.0;
			mae_est_std[em] = mse_est_std[em] = 0.0;
		}

		// for each run
		for (rn = 0; rn < RUN_NUM; rn++){
			// calculate noisy distribution from events by adding Lap --> NEventDist
			CalcDistLap(Event, UserNum, NEventDist, eps[ep], obf_mechanism);
//			for(i=0;i<CatNum;i++) cout << i << " " << EventDist[i] << " " << NEventDist[i] << endl;

			// for each distribution estimation method
			for (em = 0; em < est_method_num; em++) {

				// distribution estimation --> EstEventDist
				for (i = 0; i < CatNum; i++) EstEventDist[i] = NEventDist[i];

				// w/ normalization
				if (est_method[em] == 1 && eps[em] != -1) {
					zero_cat_num = 0;
					// set 0 for negative probabilities
					for (i = 0; i < CatNum; i++) {
						if (EstEventDist[i] < 0.0){
							EstEventDist[i] = 0.0;
							zero_cat_num++;
						}
					}

					/*
					// sum of estimated distribution --> est_dist_sum
					est_dist_sum = 0.0;
					for (i = 0; i < CatNum; i++) est_dist_sum += EstEventDist[i];

					// assign uniform distribution for zeros
					for (i = 0; i < CatNum; i++){
						if (est_dist_sum <= (double)UserNum && EstEventDist[i] == 0) EstEventDist[i] = ((double)UserNum - est_dist_sum) / zero_cat_num;
					}
					*/

					// normalize so that the sum is 1
					if(NORMALIZE){
						est_dist_sum = 0.0;
						for (i = 0; i < CatNum; i++) est_dist_sum += EstEventDist[i];
						if (est_dist_sum != 0) {
							for (i = 0; i < CatNum; i++) EstEventDist[i] = EstEventDist[i] / est_dist_sum;
						}
						else {
							for (i = 0; i < CatNum; i++) EstEventDist[i] = 1.0 / (double)CatNum;
						}
					}

				}

//				for(i=0;i<CatNum;i++) cout << i << " " << EventDist[i] << " " << NEventDist[i] << " " << EstEventDist[i] << endl;
//				est_dist_sum = 0.0;
//				for (i = 0; i < CatNum; i++) est_dist_sum += EstEventDist[i];
//				cout << est_dist_sum << endl;

				// MAE between true event distribution and estimated distribution --> mae_est
				mae_est = MAE(EventDist, EstEventDist);
				mae_est_avg[em] += mae_est;
				mae_est_std[em] += (mae_est * mae_est);

				// MSE between true event distribution and estimated distribution --> mse_est
				mse_est = MSE(EventDist, EstEventDist);
				mse_est_avg[em] += mse_est;
				mse_est_std[em] += (mse_est * mse_est);
			}
		}

		// calculate average MAE and MSE
		for (em = 0; em < est_method_num; em++){
			mae_est_avg[em] /= (double)RUN_NUM;
			mse_est_avg[em] /= (double)RUN_NUM;
		}

		// calculate standard deviation of MAE and MSE
		for (em = 0; em < est_method_num; em++){
			dtmp = mae_est_std[em] / (double)RUN_NUM - mae_est_avg[em] * mae_est_avg[em];
			if (dtmp <= 0.0) mae_est_std[em] = 0.0;
			else mae_est_std[em] = sqrt(dtmp);
			dtmp = mse_est_std[em] / (double)RUN_NUM - mse_est_avg[em] * mse_est_avg[em];
			if (dtmp <= 0.0) mse_est_std[em] = 0.0;
			else mse_est_std[em] = sqrt(dtmp);
		}

		// output results to result file
		for (em = 0; em < est_method_num; em++){
			printf("EPSILON:%f, EST_METHOD:%d, MAE: %e, MSE: %e\n", eps[ep], est_method[em], mae_est_avg[em], mse_est_avg[em]);

			sprintf(resfile, "%s_dt%d_om%d_em%d.csv", RES_FILE, data_type, obf_mechanism, est_method[em]);
			fp = FileOpen(resfile, "a");
			fprintf(fp, "%f,%e,%e,%e,%e\n", eps[ep], mae_est_avg[em], mse_est_avg[em], mae_est_std[em], mse_est_std[em]);
			fclose(fp);
		}
	}

	// free
	free1D(EventDist);
	free1D(NEventDist);
	free1D(EstEventDist);
	free1D(Event);

	return 0;
}


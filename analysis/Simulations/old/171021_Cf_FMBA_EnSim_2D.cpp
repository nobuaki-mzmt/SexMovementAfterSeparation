// Separation Search

// Encounter sim for Brownian walk with various pauses
// 170924 N Mizumoto created

// summary
// 2D condition without boundary (NBC)
// a male and a female search for each other
// movement from experimental data of termites (before/after - male/female)

// 1 step = 0.2 sec
// CRW (resample every 0.2 sec)
// pause/move cut in 0.2*n sec

// libraries
#include "stdafx.h"
#include <iostream>	// cout; cin
#include <sstream>	// to_string
#include <fstream>	// writing
using namespace std;
#include "MovementSimulation.hpp"

///*// opencv
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;
//*/

//// parameters ////

double rho[4] = { 0.8229967, 0.8541201, 0.7921616, 0.8388241 };
double speed[4] = { 12.5609, 22.37913, 19.22387, 20.86983 };

// space
double LocX[2], LocY[2], Angle[2];

int Pnokori[2];
int Wnokori[2];

vector<int> PP(4);

// result
std::string FileName;

int steps, nsteps;
int detection;
int iter, OK, search_mode;
int i, j, simtime, endtime;
double xdis, ydis, sepdis, L, sep_dis;

// functions


double Rn;
// initialization (1:f-after, 2:f-before, 3:m-after, 4:m-before)
int initial(int id, int pattern) {
	Rn = rnd();
	Angle[id] = allocate();
	switch (pattern) {
	case 1:
		if (Rn <= 0.5) { Pnokori[id] = rbeki(2.0186); Wnokori[id] = 0; }
		else { Pnokori[id] = 0; Wnokori[id] = rbeki(1.8508); }
		break;
	case 2:
	case 3:
	case 4:
		break;
	}
	return(0);
}

// change position (1:f-before, 2:f-after, 3:m-before, 4:m-after)
int change_position(int &j, int pattern) {
	if (pattern < 2) {
		if (Pnokori[j] > 0) {
			Pnokori[j] -= 1;
			if (Pnokori[j] <= 0) {
				Wnokori[j] = rbeki(1.8508);
				Pnokori[j] = 0;
			}
		}
		else {
			Angle[j] += WrappedCauchy(rho[pattern - 1]);
			LocX[j] += cos(Angle[j]) * 0.2 * speed[pattern - 1];
			LocY[j] += sin(Angle[j]) * 0.2 * speed[pattern - 1];
			Wnokori[j] -= 1;
			if (Wnokori[j] <= 0) {
				Pnokori[j] = rbeki(2.0186);
				Wnokori[j] = 0;
			}
		}
	}
	else {
		Angle[j] += WrappedCauchy(rho[pattern - 1]);
		LocX[j] += cos(Angle[j]) * 0.2 * speed[pattern - 1];
		LocY[j] += sin(Angle[j]) * 0.2 * speed[pattern - 1];
	}
	return(0);
}


// Sep search
double Sep_search(vector<int>& PP, int& iter, int& nsteps, vector<int>& res, double& sep_dis) {
	for (i = 0; i < iter; i++) {
		printf("female:%d, male:%d, rep:%d        \r", PP[0], PP[1], i + 1);
		// initialization
		LocX[0] = 0;
		LocY[0] = 0;
		LocX[1] = sep_dis;
		LocY[1] = 0;

		// begins from pause or walk randomly
		for (j = 0; j < 2; j++) {
			initial(j, PP[j]);
		}
		endtime = nsteps;
		for (simtime = 0; simtime < nsteps + 1; simtime++) {
			// change in position
			for (j = 0; j < 2; j++) {
				change_position(j, PP[j]);
			}
			// encounter determination
			xdis = LocX[1] - LocX[0];
			ydis = LocY[1] - LocY[0];
			sepdis = sqrt(xdis*xdis + ydis*ydis);
			if (sepdis <= detection) {
				endtime = simtime + 1;
				break;
			}
			else if (sepdis > (nsteps - simtime) / 2) {
				break;
			}
			/*///*
			cv::Mat Img(cv::Size(500, 500), CV_8UC3, cv::Scalar(255, 255, 255));
			cv::namedWindow("search image", cv::WINDOW_AUTOSIZE);
			cv::circle(Img, cv::Point(LocX[0] + (500 / 2), LocY[0] + (500 / 2)), 7, cv::Scalar(200, 0, 0), -1, CV_AA);
			cv::circle(Img, cv::Point(LocX[1] + (500 / 2), LocY[1] + (500 / 2)), 7, cv::Scalar(0, 0, 200), -1, CV_AA);
			cv::imshow("search image", Img);
			cv::waitKey(1);
			*///*/
		}

		for (j = endtime; j < nsteps; j++) {
			res[j] ++;
		}
	}
	return(0);
}

double Ran_search(vector<int>& PP, int& iter, int& nsteps, vector<int>& res, double& L) {
	for (i = 0; i < iter; i++) {
		printf("female:%d, male:%d, rep:%d        \r", PP[0], PP[1], i + 1);
		// initialization
		LocX[0] = 0;
		LocY[0] = 0;
		LocX[1] = rnd()*L;
		LocY[1] = rnd()*L;

		// begins from pause or walk randomly
		for (j = 0; j < 2; j++) {
			initial(j, PP[j]);
		}
		endtime = nsteps;
		for (simtime = 0; simtime < nsteps + 1; simtime++) {
			// change in position
			for (j = 0; j < 2; j++) {
				change_position(j, PP[j]);
			}
			LocX[0] = aPBC(LocX[0], L);
			LocX[1] = aPBC(LocX[1], L);
			LocY[0] = aPBC(LocY[0], L);
			LocY[1] = aPBC(LocY[1], L);


			// encounter determination
			xdis = min(abs(LocX[1] - LocX[0]), abs(L - abs(LocX[1] - LocX[0])));
			ydis = min(abs(LocY[1] - LocY[0]), abs(L - abs(LocY[1] - LocY[0])));
			sepdis = sqrt(xdis*xdis + ydis*ydis);
			if (sepdis <= detection) {
				endtime = simtime + 1;
				break;
			}
			/*///*
			cv::Mat Img(cv::Size(L, L), CV_8UC3, cv::Scalar(255, 255, 255));
			cv::namedWindow("search image", cv::WINDOW_AUTOSIZE);
			cv::circle(Img, cv::Point(LocX[0], LocY[0]), 7, cv::Scalar(0, 0, 200), -1, CV_AA);
			cv::circle(Img, cv::Point(LocX[1], LocY[1]), 7, cv::Scalar(200, 0, 0), -1, CV_AA);
			cv::imshow("search image", Img);
			cv::waitKey(1);
			*///*///*/
		}

		for (j = endtime; j < nsteps; j++) {
			res[j] ++;
		}
	}
	return(0);
}

// main function
int main() {

	cout << "Encounter simulation for Female (exp-surr) with Brownian Walers" << endl;

	cout << "Analysisng steps" << endl;
	cout << "1 step = 0.2 sec" << endl;
	cout << "60 sec = 300 steps" << endl;
	steps = 300;

	cout << "Detection = 10 mm" << endl;
	detection = 10;

	cout << "Enter the num of rep (enter number)" << endl;
	cin >> iter;

	cout << "Enter the search condition" << endl;
	cout << "0: random search, 1: separation search" << endl;
	cin >> search_mode;

	if (search_mode == 0) {
		L = 223.606;
		cout << "You selected random search" << endl;
		cout << "Enter the size of system (PBC)" << endl;
		cout << "Default: 100*5^0.5 = 223.606 mm, OK?" << endl;
		cin >> OK;
		if (OK != 1) {
			cout << "Enter (e.g., High: 50*5^0.5 = 111.803 mm)" << endl;
			cin >> L;
		}
	}
	else {
		sep_dis = 22.97364;
		cout << "You selected separation search" << endl;
		cout << "The separated distance (NBC)" << endl;
		cout << "Default: 22.97364 mm, OK?" << endl;
		cin >> OK;
		if (OK != 1) {
			cout << "Enter the separation distance (mm)" << endl;
			cin >> sep_dis;
		}
	}

	cout << "Are you OK?" << endl;
	cout << "Yes: 1, No: 2" << endl;
	cin >> OK;
	if (OK == 2) { return(2); }

	// output data
	nsteps = steps;
	vector<int> res(nsteps);

	FileName = "EnSim_FMBA";
	FileName += "_rep";
	FileName += std::to_string(iter);
	if (search_mode == 0) {
		FileName += "_RanSearch_L";
		FileName += std::to_string(int(L));
	}
	else {
		FileName += "_SepSearch_SepDis";
		FileName += std::to_string(int(sep_dis));
	}
	FileName += ".csv";
	std::ofstream ofs(FileName);

	ofs << "Female,Male";
	for (i = 0; i < nsteps; i++) {
		ofs << "," << i + 1;
	}
	ofs << endl;

	int MPP[6] = { 3,3,4,4,1,3 };
	int FPP[6] = { 1,2,1,2,1,3 };
	int p;
	for (p = 0; p < 6; p++) {		// 3,3,4,4,1,3
		PP[1] = MPP[p];
		PP[0] = FPP[p];
		fill(res.begin(), res.end(), 0);
		if (search_mode == 0) {
			Ran_search(PP, iter, nsteps, res, L);
		}
		else if (search_mode == 1) {
			Sep_search(PP, iter, nsteps, res, sep_dis);
		}
		cout << endl;
		for (j = 0; j < 2; j++) {
			switch (PP[j]) {
			case 1:
				ofs << "female_after" << ",";
				break;
			case 2:
				ofs << "female_before" << ",";
				break;
			case 3:
				ofs << "male_after" << ",";
				break;
			case 4:
				ofs << "male_before" << ",";
				break;
			}
		}
		for (i = 0; i < nsteps - 1; i++) {
			ofs << res[i] << ",";
		}
		ofs << res[nsteps - 1] << endl;
	}
}


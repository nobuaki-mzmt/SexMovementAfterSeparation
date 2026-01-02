// Separation Search

// Encounter sim for Brownian walk with various pauses
// 170924 N Mizumoto created

// summary
// 2D condition without boundary (NBC)
// a male and a female search for each other
// real female vs surrogate female

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

// space
double LocX[2], LocY[2], Angle[2];

int Pnokori[2];
int Wnokori[2];

vector<int> PP(2);

// result
std::string FileName;

int steps, nsteps;
int detection;
int iter, OK;
int i, j, simtime, endtime;
double xdis, ydis, sepdis, sep_dis;

// functions

double Rn;

// initialization (1:PL_data, 2:EXP_data, 0:WO_pause, 3:Male_CRW, 4:Only_pause)
int initial(int id, int pattern) {
	Rn = rnd();
	switch (pattern) {
	case 1:
		if (Rn <= 0.5) { Pnokori[id] = rbeki(2.0186); Wnokori[id] = 0; }
		else { Pnokori[id] = 0; Wnokori[id] = rbeki(1.8508); }
		break;
	case 2:
		if (Rn <= 0.5) { Pnokori[id] = rexp(0.6620324); Wnokori[id] = 0; }
		else { Pnokori[id] = 0; Wnokori[id] = rexp(0.534403); }
		break;
	case 0:
	case 3:
	case 4:
		break;
	}
	return(0);
}

// change position (1:PL_data, 2:EXP_data, 0:WO_pause, 3:Male_CRW, 4:Only_pause)
int change_position(int &j, int pattern) {
	if (pattern == 1 || pattern == 2) {
		if (Pnokori[j] > 0) {
			Pnokori[j] -= 1;
			if (Pnokori[j] <= 0) {
				if (pattern == 1) {
					Wnokori[j] = rbeki(1.8508);
				}
				else {
					Wnokori[j] = rexp(0.534403);
				}
				Pnokori[j] = 0;
			}
		}
		else {
			Angle[j] += WrappedCauchy(0.7576250);
			LocX[j] += cos(Angle[j]) * 0.2 * 12.6;
			LocY[j] += sin(Angle[j]) * 0.2 * 12.6;
			Wnokori[j] -= 1;
			if (Wnokori[j] <= 0) {
				if (pattern == 1) {
					Pnokori[j] = rbeki(2.0186);
				}
				else {
					Pnokori[j] = rexp(0.6620324);
				}
				Wnokori[j] = 0;
			}
		}
	}
	else if (pattern == 3) {
		Angle[j] += WrappedCauchy(0.7939391);
		LocX[j] += cos(Angle[j]) * 0.2 * 19.2;
		LocY[j] += sin(Angle[j]) * 0.2 * 19.2;
	}
	else if (pattern == 0) {
		Angle[j] += WrappedCauchy(0.7576250);
		LocX[j] += cos(Angle[j]) * 0.2 * 12.6;
		LocY[j] += sin(Angle[j]) * 0.2 * 12.6;
	}
	return(0);
}


// SepSearch (1:PL_data, 2:EXP_data, 0:WO_pause, 3:Male_CRW, 4:Only_pause)
double Bwalk(vector<int>& PP, int& iter, int& nsteps, vector<int>& res) {
	for (i = 0; i<iter; i++) {
		printf("female:%d, male:%d, rep:%d        \r", PP[0], PP[1], i + 1);
		// initialization
		LocX[0] = 0;
		LocY[0] = 0;
		Angle[0] = allocate();

		LocX[1] = sep_dis;
		LocY[1] = 0;
		Angle[1] = allocate();

		// begins from pause or walk randomly
		for (j = 0; j<2; j++) {
			initial(j, PP[j]);
		}
		endtime = nsteps;
		for (simtime = 0; simtime < nsteps + 1; simtime++) {
			// change in position
			for (j = 0; j<2; j++) {
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

		for (j = endtime; j<nsteps; j++) {
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

	sep_dis = 22.97364;
	cout << "The separation distance (NBC)" << endl;
	cout << "Default: 22.97364 mm, OK?" << endl;
	cin >> OK;
	if (OK != 1) {
		cout << "Enter the separation distance (mm)" << endl;
		cin >> sep_dis;
	}

	cout << "Are you OK?" << endl;
	cout << "Yes: 1, No: 2" << endl;
	cin >> OK;
	if (OK == 2) { return(2); }

	// output data
	nsteps = steps;
	vector<int> res(nsteps);

	FileName = "EnSim_Cf_F-Surr";
	FileName += "_rep";
	FileName += std::to_string(iter);
	FileName += "_SepSearch_SepDis";
	FileName += std::to_string(int(sep_dis));
	FileName += ".csv";
	std::ofstream ofs(FileName);

	ofs << "Female,Male";
	for (i = 0; i < nsteps; i++) {
		ofs << "," << i + 1;
	}
	ofs << endl;

	PP[1] = 3;
	for (PP[0] = 0; PP[0] < 5; PP[0]++) {
		fill(res.begin(), res.end(), 0);
		Bwalk(PP, iter, nsteps, res);
		cout << endl;
		for (j = 0; j<2; j++) {
			switch (PP[j]) {
			case 0:
				ofs << "WO_pause" << ",";
				break;
			case 1:
				ofs << "Exp_data" << ",";
				break;
			case 2:
				ofs << "Sur_data" << ",";
				break;
			case 3:
				ofs << "Male_CRW" << ",";
				break;
			case 4:
				ofs << "only_pause" << ",";
				break;
			}
		}
		for (i = 0; i < nsteps - 1; i++) {
			ofs << res[i] << ",";
		}
		ofs << res[nsteps - 1] << endl;
	}
}


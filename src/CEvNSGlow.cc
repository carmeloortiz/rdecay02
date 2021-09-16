#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <random>
#include <math.h>
#include <iomanip>
#include "CEvNSGlow.hh"

double FM = 2600; //tons
double FV = 2160*(1/0.001); //liters, i.e. 2160 cubic meters times 1 liter/cubic meter
double Ar39_decay_rate = 1.0000/1500; //decay/L/s

int write_CEvNS_events(){
	//couts the process
	std::cout << "Writing CEvNS events to events.dat..." << std::endl;

	//gets the total number of CEvNS events to be written
	int total_no_events;
	total_no_events = int (return_hist_sum("CEvNSGlow/CEvNS_events/CEvNS_event_dist.dat")*FM);
	std::cout << "Total number of CEvNS events: " << total_no_events << std::endl;

	//loop over the total and for each iteration, select a time point from the histogram given in the file
	std::vector<double> t; //time points vector
	for(int i=0;i<total_no_events;i++){
		t.push_back(sample_point("CEvNSGlow/CEvNS_events/CEvNS_event_dist_timesdt.dat"));
	}

	//takes two adjacent .out files, and computes an energy spectrum for each point in time
	std::fstream outfileA;
	std::fstream outfileB;
	std::fstream outfileC;
	std::fstream out_events;
	std::fstream timekey;
	std::string line;
	std::string lineA;
	std::string lineB;
	double time;
	double test_time;
	double time_interval;
	double E_A;
	double E_B;
	double E_C;
	double N_A;
	double N_B;
	double N_C;
	double slope;
	int index;
	int i_1;
	int i_2;

	for(int i=0;i<total_no_events;i++){
		time = t.at(i);
		
		//open timekey file
		timekey.open("CEvNSGlow/CEvNS_events/outfiles_timekey.dat", std::fstream::in);

		//find the two indices such that time is in between the values corresponding to the indices
		index=0;
		while (std::getline(timekey,line)){
			i_1 = line.find(' ');
			i_2 = line.find(' ',i_1+1);
			test_time = stod(line.substr(i_1+1, i_2-i_1));
			if (test_time > time){
				i_1 = i_2;
				i_2 = line.find(' ',i_1+1);
				time_interval =  stod(line.substr(i_1+1, i_2-i_1));
				break;
			}
			index++;
		}

		outfileA.open(("CEvNSGlow/CEvNS_events/out/"+std::to_string(index-1)+".out"), std::fstream::in);
		outfileB.open(("CEvNSGlow/CEvNS_events/out/"+std::to_string(index)+".out"), std::fstream::in);
		outfileC.open("CEvNSGlow/CEvNS_events/temp.out", std::fstream::out);

		//loop over indices
		while (std::getline(outfileA,lineA)){
			std::getline(outfileB,lineB);
			//find recoil energy corresponding to that index
			i_1 = lineA.find(' ');
			E_A=stod(lineA.substr(0, i_1));
			N_A=stod(lineA.substr(i_1+1, 10));
			
			i_1 = lineB.find(' ');
			E_B=stod(lineB.substr(0, i_1));
			N_B=stod(lineB.substr(i_1+1, 10));
			
			//take the point on the line defined by points N_A, N_B
			slope = (N_B-N_A) / (time_interval);
			N_C = slope * (time-(test_time-time_interval)) + N_A;
			
			//save result to a new histogram in outfileC
			outfileC << E_A << " " << N_C << std::endl;
		}

		outfileA.close();
		outfileB.close();
		outfileC.close();

		//sample an energy from outfileC
		E_C = sample_point("CEvNSGlow/CEvNS_events/temp.out");

		//write to events.dat the event tag, time, and recoil energy
		out_events.open("CEvNSGlow/events.dat",std::fstream::app);
		out_events << "0 " << time << " " << std::setprecision(10) << E_C << std::endl;
		//std::cout << "CEvNS event (time (s), recoil energy (MeV): " << time << ", " << E_C << std::endl;
		out_events.close();
		timekey.close();
	}

	return total_no_events;
}

int write_CC_events(){
	//couts the process
	std::cout << "Writing CC decay events to events.dat..." << std::endl;
	
	//gets the total number of CC events to be written
	int total_no_events;
	total_no_events = int (return_hist_sum("CEvNSGlow/CC_events/CC_event_dist.dat")*FM);
	std::cout << "Total number of CC events: " << total_no_events << std::endl;

	//loop over the total and for each iteration, select a time point from the histogram given in the file
	std::vector<double> t; //time points vector
	for(int i=0;i<total_no_events;i++){
		t.push_back(sample_point("CEvNSGlow/CC_events/CC_event_dist.dat"));
		//std::cout << "CC event at: " << t.at(i) << " s" << std::endl;
	}

	//computes the beta value and <E_nue> corresponding to this time point
	std::vector<double> beta;
	std::vector<double> E_ave;
	for(int i=0;i<total_no_events;i++){
		beta.push_back(return_interpolated_y(t.at(i),"CEvNSGlow/CC_events/beta_dist.dat"));
		E_ave.push_back(return_interpolated_y(t.at(i),"CEvNSGlow/CC_events/E-ave_dist.dat"));
	}

	//writes to events.dat the event tag, the time, beta value, and E_ave
	std::fstream file;
	file.open("CEvNSGlow/events.dat", std::fstream::app);
	for(int i=0;i<total_no_events;i++){
		file << "1 " << t.at(i) << " " << beta.at(i) << " " << E_ave.at(i) << std::endl;
		std::cout << "CC event (time (s),beta,E_ave (MeV)): " << t.at(i) << " " << beta.at(i) << " " << E_ave.at(i) << std::endl;
	}
	file.close();
	
	return total_no_events;
}

int write_Ar39_events(){
	//initialize randomizer
	std::random_device rd;
	std::mt19937 gen(rd());

	//couts the process
	std::cout << "Writing Ar39 decay events to events.dat..." << std::endl;

	//opens file
	std::string line;
	std::string last_line;
	std::fstream file;
	file.open("CEvNSGlow/CEvNS_events/outfiles_timekey.dat", std::fstream::in);

	//finds maximum and minimum time points
	double t_min;
	double t_max;
	int i_1;
	int i_2;

	std::getline(file,line);
	i_1 = line.find(' ');
	i_2 = line.find(' ',i_1+1);
	t_min = stod(line.substr(i_1+1, i_2-i_1));

	while (std::getline(file,line)){
		last_line = line;
	}
	
	i_1 = last_line.find(' ');
	i_2 = last_line.find(' ',i_1+1);
	t_max = stod(last_line.substr(i_1+1, i_2-i_1));

	std::cout << "t_min: " << t_min << std::endl;
	std::cout << "t_max: " << t_max << std::endl;
	file.close();

	//compute event rate for Ar39 decay based on fiducial volume
	double event_rate = Ar39_decay_rate*FV; //events/second

	//for now, set event times to be evenly spaced and write to a file
	double time;
	file.open("CEvNSGlow/events.dat",std::fstream::app);

	int i;
	for(i=0;i*(1/event_rate)<t_max;i++){
		time = i*(1/event_rate);
		file << "2 " << std::setprecision(10) << time << std::endl;
		//std::cout << "created Ar39 decay event at: " << time <<  " s" << std::endl;
	}
	file.close();

	return i;
}	


double return_hist_sum(std::string filename){
	//two vectors, one for bins and the other for weights
	std::vector<double> b;
	std::vector<double> w;

	//opens file
	std::string line;
	std::fstream file;
	file.open(filename, std::fstream::in);

	//assigns first column elements to bin vector, second column to weight vector
	int i_1;
	int i_2;
	while (std::getline(file,line)){
		i_1 = line.find(' ');
		i_2 = line.find(' ',i_1+1);
		b.push_back(stod(line.substr(0, i_1)));
		w.push_back(stod(line.substr(i_1+1, i_2-i_1)));
	}

	
	//computes sum of (frequencies*bin_width)
	double sum = 0;
	for (int i=0;i<b.size()-1;i++){
		sum = sum + w.at(i)*(b.at(i+1)-b.at(i));
	}

	return sum;
}

double sample_point(std::string filename){
	//initialize randomizer
	std::random_device rd;
	std::mt19937 gen(rd());

	//two vectors, one for bins and the other for weights
	std::vector<double> b;
	std::vector<double> w;

	//opens file
	std::string line;
	std::fstream file;
	file.open(filename, std::fstream::in);

	//assigns first column elements to bin vector, second column to weight vector
	int i_1;
	int i_2;
	while (std::getline(file,line)){
		i_1 = line.find(' ');
		i_2 = line.find(' ',i_1+1);
		b.push_back(stod(line.substr(0, i_1)));
		w.push_back(stod(line.substr(i_1+1, i_2-i_1)));
	}

	//creates the probability distribution
	std::piecewise_constant_distribution<> d(b.begin(),b.end(),w.begin());
	
	file.close();
	return d(gen);
}

double return_interpolated_y(double x0, std::string filename){
	//two vectors, x and y
	std::vector<double> x;
	std::vector<double> y;

	//opens file
	std::string line;
	std::fstream file;
	file.open(filename, std::fstream::in);

	//assigns first column elements to bin vector, second column to weight vector
	int i_1;
	int i_2;
	while (std::getline(file,line)){
		i_1 = line.find(' ');
		i_2 = line.find(' ',i_1+1);
		x.push_back(stod(line.substr(0, i_1)));
		y.push_back(stod(line.substr(i_1+1, i_2-i_1)));
	}

	//performs a simple piecewise linear interpolation
	double slope;
	double y0;
	for (int i = 1; i < x.size(); i++) {
		if (x.at(i) > x0) {
			slope = (y.at(i) - y.at(i-1)) / (x.at(i)  - x.at(i-1) );
			y0 = slope * (x0 - x.at(i-1)) + y.at(i-1);
			if (y0 > 0.0000000){
				return y0;
			}
			else{
				return 0.0;
			}
		}
	}
	return 0.0;
}

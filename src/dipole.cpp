//
//  correlation.cpp
//
//  code to calculate correlation function
//  from LAMMPS trajectory
//
//  Created by Michelle Liu (2014)
//  Chemical and Biomolecular Engineering, UC Berkeley
//

#include <fstream>
#include <iostream>
#include <cmath>
#include <array>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include "common.h"

#define PI 3.14159265359

void read_dipole(std::string file_name, double* Cm,
        double* Mnet, double* plot_time, double* full_time,
        int num_particles, int num_timesteps, int max_meas_time,
        std::string suffix) {
    printf("Reading %s\ntotal time %d, meas time %d, %d particles\nsuffix: %s\n",
            &file_name[0],num_timesteps,max_meas_time,num_particles,&suffix[0]);
    using namespace std;
    filebuf file_buffer;
    string str;
    if (file_buffer.open(file_name,std::ios::in))
    {
        istream file_stream(&file_buffer);
        vector<double> line_split;

        // skip/parse beginning lines
        //for (int i=0; i<5; ++i) {
        //    file_stream.ignore(1000,'\n');
        //}

        //getline(file_stream,str);
        //split(str,' ',line_split);
        //int num_atoms=line_split[0];
        int num_atoms=3000;

        //for (int i=0; i<72; ++i) {
        //    file_stream.ignore(1000,'\n');
        //}

        plot_time=(double*) malloc(max_meas_time*sizeof(double));
        full_time=(double*) malloc(num_timesteps*sizeof(double));
        Mnet=(double*) malloc(3*num_timesteps*sizeof(double)); // net dipole of confined region
        Cm=(double*) malloc(max_meas_time*sizeof(double)); // dipole-dipole correlation

        double q_hydrogen=0.4238;
        double q_oxygen=-0.8476;
        double Mave [3];
        double Msqave [3];
        double Mfluc [3];
        double Msum [3];
        double Msqsum [3];
        double exp [3];
        double Msqtot=0;
        double Mexptot=0;
        double Mtotfluc=0;
        double Vol=0;
        double Volsq=0;
        double Temp=0;
        for (int i=0; i<3; ++i) {
            Msum[i]=0;
            Msqsum[i]=0;
            exp[i]=0;
            Mfluc[i]=0;
        }

        for (int step=0; step<num_timesteps;++step) {
            getline(file_stream,str);
            split(str,' ',line_split);
            //printf("%f %f %f %f %f %f %f\n",line_split[0],line_split[1],line_split[2],line_split[3],line_split[4],line_split[5],line_split[6]);
            for (int k=0; k<3; ++k) {
                double ox_contribution = q_oxygen*line_split[1+k];
                double hy_contribution = q_hydrogen*line_split[4+k];
                //printf("Net dipole: %f\n",ox_contribution+hy_contribution);
                Msum[k]=Msum[k]+ox_contribution+hy_contribution;
                Msqsum[k]=Msqsum[k]+(ox_contribution+hy_contribution)*(ox_contribution+hy_contribution);
                Msqtot=Msqtot+(ox_contribution+hy_contribution)*(ox_contribution+hy_contribution);
                Mnet[3*step+k]+=ox_contribution;
                Mnet[3*step+k]+=hy_contribution;
            }
            if (line_split.size() > 7) {
                Vol+=line_split[7];
                Volsq+=line_split[7]*line_split[7];
                Temp+=line_split[8];
            } else {
                printf("Error! No volume or temperature output\n");
                exit(1);
            }
            line_split.clear();
            // calculating correlation
            for (int step_meas=0; step_meas < min(step,max_meas_time); ++step_meas) {
                int occu_tmp=0;
                int t0=step-step_meas;
                for (int k=0; k<3; ++k) {
                    Cm[step_meas]+=Mnet[3*step+k]*Mnet[3*t0+k];
                }
            }
        }

        // correction to get averages
        Msqtot=Msqtot/num_timesteps;
        for (int i=0; i<3; ++i) {
            Mave[i]=Msum[i]/num_timesteps;
            Msqave[i]=Msqsum[i]/num_timesteps;
            printf("Mave[%d] = %f\n",i,Mave[i]);
            printf("Msqave[%d] = %f\n",i,Msqave[i]);
        }

        Vol=1.0*Vol/num_timesteps;
        Volsq=1.0*Volsq/num_timesteps;
        Temp=1.0*Temp/num_timesteps;

        Mexptot=Mave[0]*Mave[0]+Mave[1]*Mave[1]+Mave[2]*Mave[2];

        //double delta=0;
        //for (int i=0; i<num_timesteps; ++i) {
        //    for (int k=0; k<3; ++k) {
        //        delta=Mnet[3*i+k]-Mave[i];
        //        exp[k]+=delta/num_timesteps; // E[x-xave]
        //        Mfluc[k]+=delta*delta/num_timesteps; // E[(x-xave)^2]
        //    }
        //}

        //for (int i=0; i<3; ++i) {
        //    Mfluc[i]=Mfluc[i]-exp[i]*exp[i];
        //    Mfluc[i]=Mfluc[i]*1.0e-10; // in q*m
        //}

        for (int i=0; i<3; ++i) {
            Mfluc[i]=Msqave[i]-Mave[i];
            //Mfluc[i]=Mfluc[i]*(1.0e-20)*(1.60218e-19)*(1.60218e-19); // in (q*m)^2
            printf("Mfluc[%d] = %f\n",i,Mfluc[i]);
        }
        Mtotfluc=Msqtot-Mexptot;

        long double epsilontot=0;
        long double epsilon [3];
        // FIXME: UNITS !!!!!!
        Vol=Vol*1.0e-24; // from AA^3 to cm^3
        Volsq=Volsq*1.0e-24*1.0e-24;
        //long double V=5.12940e-20; // in cm^3 // this was just a placeholder
        long double kB=1.38065e-16; // in cgs
        //double T=298; // K
        for (int i=0; i<3; ++i) {
            epsilon[i]=1+(4.0*PI*Mfluc[i]*(1.0e-16)*(3*1.60218e-10)*(3*1.60218e-10))/(Vol*kB*Temp);
        }
        epsilontot=1+(4.0*PI*Mtotfluc*(1.0e-16)*(3*1.60218e-10)*(3*1.60218e-10))/(3*Vol*kB*Temp);

        kB=1.38065e-23; // now in SI units
        Vol=Vol*1.0e-6; // in SI (m^3)
        Volsq=Volsq*1.0e-6*1.0e-6;
        long double kappa = (Volsq - Vol*Vol)/(kB*Temp*Vol);
        kappa=kappa/(9.86923e-6); // converting from Pa^-1 to atm^-1


        printf("dielectric tensors:\n  %Lf\n  %Lf\n  %Lf\n",epsilon[0],epsilon[1],epsilon[2]);
        printf("total dielectric:   %Lf\n",epsilontot);
        printf("Volsq:  %e    Vol:  %e\n",Volsq,Vol);
        //printf("Volsq - Vol*Vol:  %e\n",Volsq - Vol*Vol);
        printf("Temp:   %f\n",Temp);
        printf("isothermal compressibility:   %Le\n",kappa);

        // normalizing
        int t0_max;
        int dump_time=20;
        for (int i=0; i<max_meas_time; ++i) {
            t0_max=num_timesteps-i;
            Cm[i] = (Cm[i])/(t0_max);
            plot_time[i] = (i*dump_time)/1000.0; // now in ps
        }
        for (int i=0; i<num_timesteps; ++i) {
            full_time[i] = (i*dump_time)/1000.0; // now in ps
        }
        ofstream CMout;
        CMout.open("results/C_M-"+suffix+"-w"+to_string(num_particles)+"_m"+to_string(max_meas_time)
                +"_t"+to_string(num_timesteps)+"_d"+to_string(dump_time)+"fs");
        ofstream MNETout;
        MNETout.open("results/Mnet-"+suffix+"-w"+to_string(num_particles)+"_m"+to_string(max_meas_time)
                +"_t"+to_string(num_timesteps)+"_d"+to_string(dump_time)+"fs");
        CMout << "#t Cm" << endl;
        for (int i=0; i< max_meas_time; ++i) {
            CMout << plot_time[i] << " " << Cm[i] << endl;
        }
        for (int i=0; i< num_timesteps; ++i) {
            MNETout << full_time[i] << " " << Mnet[3*i+0] << " " << Mnet[3*i+1] << " " << Mnet[3*i+2] << endl;
        }
        CMout.close();
        MNETout.close();
        file_buffer.close();
    }
    else {
        cout << "ERROR opening file!" << endl;
    }

}

int main(int argc, const char * argv[]) {
    using namespace std;
    int* theta;

    double* MSD;
    double* Cm;
    double* Mnet;
    double* mu;

    double* plot_time;
    double* full_time;

    double x0=-6.05;
    double x1=6.05;
    double y0=-5.25;
    double y1=5.25;
    double z0=16.1;
    double z1=22.90;

    int skip_steps=0;
    int num_particles=1000;
    int num_timesteps=200;
    int max_meas_time=100;
    int O_id=2;
    clock_t start, time;
    string suffix;

    string file_name=argv[argc-1];
    for (int i = 0; i < (argc-1); ++i) {
        if (strcmp("-t", argv[i])==0) {
            num_timesteps = atoi(argv[i+1]);
        }
        if (strcmp("-d", argv[i])==0) {
            double z_dist = atof(argv[i+1]);
            double z_center = 19.5;
            z0=z_center-z_dist/2;
            z1=z_center+z_dist/2;
        }
        if (strcmp("-m", argv[i])==0) {
            max_meas_time = atoi(argv[i+1]);
        }
        if (strcmp("-s", argv[i])==0) {
            skip_steps = atoi(argv[i+1]);
        }
        if (strcmp("-e", argv[i])==0) {
            suffix = string(argv[i+1]);
        }
        if (strcmp("-o", argv[i])==0) {
            O_id = atoi(argv[i+1]);
        }
        if (strcmp("-w", argv[i])==0) {
            num_particles = atoi(argv[i+1]);
        }
        if (strcmp("-C", argv[i])==0) {
            x0+=15;
            x1+=15;
            y0+=15;
            y1+=15;
        }
    }
    start=clock();

    read_dipole(file_name,Cm,Mnet,plot_time,full_time,
            num_particles,num_timesteps,max_meas_time,suffix);

    time=clock()-start;
    time=time/CLOCKS_PER_SEC;
    cout << "Took " << time << " sec to run" << endl;
    cout << "--------------------" << endl;
    return 0;
}

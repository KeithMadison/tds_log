#include <cstdio>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include "getopt.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraph.h"

#define MAX_CNT 1000

const int ibottom_sample = 34200; // sample of bottom reflection
const int nsample_reflect = 256;
const int pre_reflection_i0 = 16000;
const int post_reflection_i0 = 35000;
const int nsample_prepost = 16384;
const int MVbin = 2048;
const int Mtseg = 100;
const int iVmax = 2048;
const float Nsighit = 5.5;

bool SigHit = false;
int sideband_i0;
int sideband_len;
int signal_i0;
int signal_len;
int nptsid_5us;

float abssidsum_5us;
float Vavghist[Mtseg];
float Vhist[MVbin][Mtseg], timehist[Mtseg];
int iVbin, itbin;
int iavghist[Mtseg];

using namespace std;

void read_bin(char *basename);
void read_dat(char *basename);
long convert_utc(char **stamp);
void do_root_1ch(char *stamp, char *comment, int numBins, double minX, double maxX, double *x, double *y);
void do_root_4ch(char *stamp, char *comment, int numBins, double minX, double maxX, double *x, double *y1, double *y2, double *y3, double *y4);

int debug = 0;   // debug mode
int gifs = 0;    // create gif files
int epscripts = 0;  // create encapsulated postscript files
int noroots = 0;  // do not create ROOT files

TNtuple *ntuple_wave;
TFile *results;

void usage(char *execname) {
    printf("usage: %s [-dghprv] filename\n", execname);
    printf("\t-d\tdebug mode\n");
    printf("\t-g\tcreate gif files\n");
    printf("\t-p\tcreate encapsulated postscript files\n");
    printf("\t-r\tdo not create ROOT files; default action is to create ROOT files\n");
    printf("\t-h\tthis help menu\n");
    printf("\t-z\tPerform time-domain sideband subtraction\n");
    printf("\t-v\tversion\n");
    exit(1);
}

int main(int argc, char *argv[]) {
    char out_filename[500];
    sprintf(out_filename, "tcal.root");

    results = new TFile(out_filename, "RECREATE", "Title", 9);
    results->cd();

    ntuple_wave = new TNtuple("tcal", "tcal", "ch:v_max:t_max");

    signal_i0 = 36720;
    signal_len = 550;
    sideband_i0 = 34000;
    sideband_len = 2000;

    char clswitch;
    if (argc > 1) {
        while ((clswitch = getopt(argc, argv, "dghprv")) != EOF) {
            switch (clswitch) {
                case 'd':
                    debug = 1;
                    cout << "DEBUG MODE" << endl << endl;
                    break;
                case 'g':
                    gifs = 1;
                    break;
                case 'p':
                    epscripts = 1;
                    break;
                case 'r':
                    noroots = 1;
                    break;
                case 'h':
                    usage(argv[0]);
                    break;
                case 'v':
                    if (sideband_i0 == 0.) {
                        signal_i0 = 36720;
                        signal_len = 550;
                        sideband_i0 = 34000;
                        sideband_len = 2000;
                    }
                    break;
            }
        }
    } else {
        usage(argv[0]);
        exit(1);
    }

    if (optind == argc)
        usage(argv[0]);

    char *basename, *extension;
    if (!strchr(argv[argc - 1], '.')) {
        cout << "error: file must be of type bin, hdr, dat" << endl;
        exit(1);
    }
    basename = strtok(argv[argc - 1], ".");
    extension = strtok(NULL, "\n");

    if (strstr(extension, "bin") || strstr(extension, "hdr")) {
        read_bin(basename);
    } else {
        if (strstr(extension, "dat")) {
            read_dat(basename);
        } else {
            cout << "error: file must be of type bin, hdr, dat" << endl;
            exit(1);
        }
    }

    results->Close();
}

void read_bin(char *basename) {
    ifstream data;
    char filename[256], name[256];
    int nsample, nchan, type, nevent;
    float *buffer = NULL;

    sprintf(filename, "%s.bin", basename);
    data.open(filename, ios::in | ios::binary);

    if (!data) {
        cout << "error: unable to open file " << filename << endl;
        exit(1);
    }

    data.read((char *)&nsample, sizeof(int));
    data.read((char *)&nchan, sizeof(int));
    data.read((char *)&type, sizeof(int));
    data.read((char *)&nevent, sizeof(int));

    if (debug)
        cout << "nsample: " << nsample << ", nchan: " << nchan << ", type: " << type << ", nevent: " << nevent << endl;

    int n = 0;
    while (1) {
        if (n % 1000 == 0 && n != 0)
            cout << n << " ";
        else if (n % 100 == 0 && n != 0)
            cout << ".";
        else if (n % 10 == 0 && n != 0)
            cout << ",";
        else if (n == 0)
            cout << "reading " << filename << " data: ";

        buffer = new float[nsample];
        data.read((char *)buffer, sizeof(float) * nsample);
        if (!data)
            break;

        if (n % 1000 == 0 && n != 0)
            cout << n << " ";
        else if (n % 100 == 0 && n != 0)
            cout << ".";
        else if (n % 10 == 0 && n != 0)
            cout << ",";
        else if (n == 0)
            cout << "reading " << filename << " data: ";

        n++;
    }

    cout << endl << "read " << n << " waveforms" << endl;
    delete buffer;
    data.close();
}

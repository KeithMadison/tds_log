#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>
#include <string>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#define MAX_CNT 1000

TString options;

void do_root(char* stamp, char* comment, int num_bins, double min_x, double max_x, double* x, double* y, int num_ch);
void read_bin(const char* basename);
void read_dat(const char* basename);
long convert_utc(char* stamp, int* ms);

void usage() {

   printf("usage: '[-dghprv]', 'filename'\n");
   printf("\t-g\tcreate gif files\n");
   printf("\t-h\tthis help menu\n");
   printf("\t-p\tcreate encapsulated postscript files\n");
   printf("\t-r\tdo not create ROOT files; default action is to create ROOT files\n");
   exit(1);

}

void tdslog4_display(TString args = " ", TString filename = " ") {

   if (args.Contains('h') || args == " " || filename == " ")
      usage();

   auto basename = ((TObjString*)(filename.Tokenize(".")->At(0)))->String();
   auto extension = ((TObjString*)(filename.Tokenize(".")->At(1)))->String();

   options = args;

   if (extension == "bin" || extension == "hdr") {

      read_bin(basename.Data());

   } else if (extension == "dat") {

      read_dat(basename.Data());

   } else {

      std::cout << "error: file must be of type bin, hdr, dat" << std::endl;
      exit(1);

   }

}

void read_dat(const char* basename) {

   char filename[80], stamp[160], line[160], comment[260];
   sprintf(filename, "%s.dat", basename);
   std::ifstream fin(filename);

   if (!fin) {

      std::cout << "error: unable to open " << filename << std::endl;
      exit(1);

   }

   int num_ch = 1;
   double x[20000], y[4][20000];
   char* val;

   fin.getline(stamp, 160, '\n');

   while (fin.good()) {

      if (line[0] == '#')
         strcpy(stamp, line);

      fin.getline(line, 160, '\n');
      int bins = 0;
      x[bins] = atof(strtok(line, " "));
      y[0][bins] = atof(strtok(NULL, " "));
      val = strtok(NULL, " ");

      if (val != " ") {

         num_ch = 4;
         x[bins] = atof(val);
         y[1][bins] = atof(strtok(NULL, " "));

         for (int i = 2; i < 4; ++i) {

            x[bins] = atof(strtok(NULL, " "));
            y[i][bins] = atof(strtok(NULL, " "));

         }

      }

      ++bins;

      while (1) {

         if (!(fin.getline(line, 160, '\n'))) break;
         if (line[0] == '#') break;

         x[bins] = atof(strtok(line, " "));
         y[0][bins] = atof(strtok(NULL, " "));

         if (num_ch != 1) {

            for (int i = 1; i < 4; ++i) {

               x[bins] = atof(strtok(NULL, " "));
               y[i][bins] = atof(strtok(NULL, " "));

            }

         }

         ++bins;
      }

      do_root(stamp, comment, bins, x[0], x[bins - 1], x, (double*)y, num_ch);
   }

   fin.close();

}

void read_bin(const char* basename) {

   char filename[80], of1[80], of0[80];
   char hold_line[1000];
   char stamps[MAX_CNT][160];
   char comment[260];
   char* wave[MAX_CNT];
   long filesize;
   int iter;

   FILE* hdrfile, * binfile, * logfile, * ch1out, * ch0out;

   sprintf(of1, "ch1_%s.txt", basename);
   sprintf(of0, "ch0_%s.txt", basename);

   ch1out = fopen(of1, "w");
   ch0out = fopen(of0, "w");

   sprintf(filename, "%s.hdr", basename);

   if (!(hdrfile = fopen(filename, "rb"))) {

      std::cout << "error: unable to open " << filename << std::endl;
      exit(1);

   }

   sprintf(filename, "%s.bin", basename);

   if (!(binfile = fopen(filename, "rb"))) {

      std::cout << "error: unable to open " << filename << std::endl;
      exit(1);

   }

   sprintf(filename, "%s.log", basename);

   if (!(logfile = fopen(filename, "rb"))) {

      std::cout << "error: unable to open " << filename << std::endl;
      exit(1);

   }

   fgets(hold_line, 1000, logfile);
   iter = 0;

   while (hold_line[iter] != '\n') {

      comment[iter] = hold_line[iter];
      ++iter;

   }

   std::cout << comment << std::endl;

   fgets(hold_line, 1000, hdrfile);
   iter = 0;

   while (hold_line[iter] != '\n') {

      stamps[0][iter] = hold_line[iter];
      ++iter;

   }

   char curr_ch = stamps[0][iter - 1];

   fgets(hold_line, 1000, hdrfile);
   strtok(strstr(hold_line, "#"), " ");

   int length = atoi(strtok(NULL, "\n"));
   double yzero[MAX_CNT], yoff[MAX_CNT], ymult[MAX_CNT], pt_off[MAX_CNT], xincr[MAX_CNT];

   int width, num_ch = 1, waveiter = 0;

   do {

      strtok(strstr(hold_line, "YZERO"), " ");
      yzero[waveiter] = atof(strtok(NULL, ";"));

      strtok(strstr(hold_line, "YOFF"), " ");
      yoff[waveiter] = atof(strtok(NULL, ";"));

      strtok(strstr(hold_line, "YMULT"), " ");
      ymult[waveiter] = atof(strtok(NULL, ";"));

      strtok(strstr(hold_line, "PT_OFF"), " ");
      pt_off[waveiter] = atoi(strtok(NULL, ";"));

      strtok(strstr(hold_line, "XINCR"), " ");
      xincr[waveiter] = atof(strtok(NULL, ";"));

      if (!waveiter) {

         strtok(strstr(hold_line, "BYT_NR"), " ");
         width = atoi(strtok(NULL, ";"));

      }

      ++waveiter;
      fgets(hold_line, 1000, hdrfile); // Carriage return

      if (fgets(hold_line, 1000, hdrfile)) {
         iter = 0;

         while (hold_line[iter] != '\n') {

            stamps[waveiter][iter] = hold_line[iter];
            ++iter;

         }

         if (stamps[waveiter][iter - 1] != curr_ch) {

            ++num_ch;
            curr_ch = stamps[waveiter][iter - 1];

         }

         fgets(hold_line, 1000, hdrfile);
      }

   } while (!feof(hdrfile));

   int num_acq = waveiter;
   fseek(binfile, 0, SEEK_END);
   rewind(binfile);

   for (iter = 0; iter < num_acq; iter++) {

      if ((wave[iter] = (char*)malloc((size_t)(length))) == NULL) {

         std::cout << "error: out of memory" << std::endl;
         exit(1);

      }

   }

   iter = 0;
   size_t size;

   while (!feof(binfile)) {

      size = fread(wave[iter], sizeof(char), length, binfile);
      iter++;

   }

   std::cout << "done: " << wave[0][0] << " " << wave[1][0] << std::endl;

   fclose(binfile);

   std::cout << "length = " << length << ", width = " << width << std::endl;

   double x[length / width], y[4][length / width];
   float twob;

   for (int acq = 0; acq < num_acq; acq += num_ch) {

      for (int ch = 0; ch < num_ch; ch++) {

         for (int pt = 0; pt < length; pt += width) {

            if (width == 1) {

               x[pt] = (pt - pt_off[acq + ch]) * xincr[acq + ch];
               y[ch][pt] = ((wave[acq + ch][pt] - yoff[acq + ch]) * ymult[acq + ch]) + yzero[acq + ch];
               float temp = wave[acq + ch][pt];

               if (ch == 0)
                  fprintf(ch0out, "1 %g %g\n", x[pt], y[ch][pt]);
               else
                  fprintf(ch1out, "2 %g %g\n", x[pt], y[ch][pt]);

            } else {

               short* ptr;
               ptr = (short*)& wave[acq + ch][pt];
               twob = (float)(*ptr);
               x[pt / 2] = ((pt / 2) - pt_off[acq + ch]) * xincr[acq + ch];
               y[ch][pt / 2] = ((twob - yoff[acq + ch]) * ymult[acq + ch]) + yzero[acq + ch];
               float temp = wave[acq + ch][pt];

               if (ch == 0)
                  fprintf(ch0out, "3 %g %g\n", x[pt / 2], y[ch][pt / 2]);
               else
                  fprintf(ch1out, "4 %g %g\n", x[pt / 2], y[ch][pt / 2]);

            }

         }

      }
      do_root(stamps[acq], comment, (int)length / width, x[0], x[length / width - 1], x, (double*)y, num_ch);

   }

}

void do_root(char* stamp, char* comment, int num_bins, double min_x, double max_x, double* x, double* y, int num_ch) {

   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0);

   TCanvas* c1 = new TCanvas();
   TLegend* leg = new TLegend(0.88, 0.85, 0.99, 0.99);

   double time, voltage_ch[4];

   TTree* scope_data = new TTree("scope_data", "scope_data");
   scope_data->Branch("time", &time, "time/D");

   TH1D** histo = new TH1D*[num_ch];

   for (int i = 0; i < num_ch; ++i) {

      histo[i] = new TH1D(Form("h%d", i), Form("%s;Time [seconds];Amplitude [volts]", stamp), num_bins, min_x, max_x);
      scope_data->Branch(Form("voltage_ch%d", i + 1), &voltage_ch[i], Form("voltage_ch%d/D", i + 1));

      for (int j = 0; j < num_bins; ++j) {

         voltage_ch[i] = *(y + num_ch * i + j);
         time = x[j];

         scope_data->Fill();
         histo[i]->Fill(x[j], *(y + num_ch * i + j));

      }

      leg->AddEntry(histo[i], Form("Channel %d", i + 1));

      histo[i]->SetLineColor(i + 1);
      histo[i]->Draw("SAME");

   }

   leg->Draw();
   int ms;
   long timestamp = convert_utc(stamp, &ms);

   if (options.Contains('g'))
      c1->Print(Form("%ld+%d.gif", timestamp, ms));

   if (options.Contains('p'))
      c1->Print(Form("%ld+%d.eps", timestamp, ms));

   if (!options.Contains('r')) {

      TFile* fout = new TFile(Form("%ld+%d.root", timestamp, ms), "RECREATE");

      for (int i = 0; i < num_ch; ++i)
         histo[i]->Write();

      leg->Write();
      c1->Write();
      scope_data->Write();

      fout->Close();

   }

}

long convert_utc(char* stamp, int* ms) {

   struct tm t_utc;
   strtok(strstr(stamp, "UTC"), " ");
   strtok(NULL, " ");
   std::string month = strtok(NULL, " ");
   std::map<std::string, int> months {{"Jan", 0}, {"Feb", 1}, {"Mar", 2}, {"Apr", 3}, {"May", 4}, {"Jun", 5},
                                      {"Jul", 6}, {"Aug", 7}, {"Sep", 8}, {"Oct", 9}, {"Nov", 10}, {"Dec", 11}};
   t_utc.tm_mon = months.find(month)->second;
   t_utc.tm_mday = atoi(strtok(NULL, " "));
   t_utc.tm_hour = atoi(strtok(NULL, ":"));
   t_utc.tm_min = atoi(strtok(NULL, ":"));
   t_utc.tm_sec = atoi(strtok(NULL, " "));
   t_utc.tm_year = atoi(strtok(NULL, " ")) - 1900;

   strtok(NULL, " ");

   *ms = atoi(strtok(NULL, " "));
   time_t utc = mktime(&t_utc);

   return (long)utc;

}

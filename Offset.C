#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <memory>
#include <string>
#include <cmath>         // For math functions (with _USE_MATH_DEFINES)
#include <sys/stat.h>    // For file system operations
#include <sys/types.h>   // For file system types
#include <cassert>       // For debugging assertions
#include <thread>        // For sleeping

// ROOT Core Libraries
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TKey.h"
#include "TParameter.h"
#include <TString.h>      // For ROOT string handling
#include <TObjArray.h>    // For working with ROOT arrays

// ROOT Graphical Libraries
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include <TPaletteAxis.h>
#include <TText.h>
#include <TLatex.h>
#include <TLegend.h>
#include "TMultiGraph.h"

// ROOT Graphing and Histogram Libraries
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>
#include <TGraphErrors.h>

// ROOT Random Number Generator
#include <TRandom3.h>

// ROOT Fitting Libraries
#include <TF1.h>

// Additional ROOT Utilities
#include <TPaveText.h>
#include <TArrayF.h>
#include <unordered_map>
#include <TVirtualFFT.h>
#include <complex>
#include <TComplex.h>
#include <map>
#include <vector>
#include <cmath>
#include <THStack.h>

#include "TSpectrum.h"

#include <cstdlib> // for rand()
// Macro Definition
#define _USE_MATH_DEFINES

using namespace std;

const short EMPTY_BIN = -1;
//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Functions //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

// Calculate the OM number from the calo udd variables
int calculate_om_num(std::vector<int> *calo_type, std::vector<int> *calo_side, 
                     std::vector<int> *calo_wall, std::vector<int> *calo_column, 
                     std::vector<int> *calo_row, int k) {
    int om_num = -1;
    if (calo_type->at(k) == 0) {
        om_num = calo_side->at(k) * 260 + calo_column->at(k) * 13 + calo_row->at(k);
    } else if (calo_type->at(k) == 1) {
        om_num = 520 + calo_side->at(k) * 64 + calo_wall->at(k) * 32 + calo_column->at(k) * 16 + calo_row->at(k);
    } else if (calo_type->at(k) == 2) {
        om_num = 648 + calo_side->at(k) * 32 + calo_wall->at(k) * 16 + calo_column->at(k);
    }
    return om_num;
}

////////////////////////////// FUNCTIONS from the WAVEFORMS ////////////////////////////////////
// Function to calculate baseline
double calculate_baseline976(const std::vector<short>& waveform) {
    double baseline = 0;
    for (int i = 0; i < 976; ++i) {
        baseline += waveform[i];
    }
    return baseline / 976;
}
// Function to calculate baseline
double calculate_baseline96(const std::vector<short>& waveform) {
    double baseline = 0;
    for (int i = 0; i < 96; ++i) {
        baseline += waveform[i];
    }
    return baseline / 96;
}
double calculate_stddev976(const std::vector<short>& waveform, double baseline) {
    double sum_sq_diff = 0;
    for (int i = 0; i < 976; ++i) {
        sum_sq_diff += pow(waveform[i] - baseline, 2);
    }
    return sqrt(sum_sq_diff / 976);
}

//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Main ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    // --- Get run number from command line ---
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run_number>" << std::endl;
        return 1;
    }
    int run_number = std::stoi(argv[1]);

    // --- Build filenames dynamically ---
    std::string input_filename  = "snemo_run-" + std::to_string(run_number) + "_udd.root";
    std::string output_filename = "waveform_type_cut_" + std::to_string(run_number) + ".root";

    std::cout << "Reading: "  << input_filename  << std::endl;
    std::cout << "Writing: "  << output_filename << std::endl;

    // Open the ROOT file
    TFile *file = new TFile(input_filename.c_str(), "READ");
    TTree *tree = (TTree *)file->Get("SimData");

    // Variables for branches
    std::vector<std::vector<short>> *wave = nullptr;
    std::vector<int> *calo_wall = nullptr;
    std::vector<int> *calo_side = nullptr;
    std::vector<int> *calo_column = nullptr;
    std::vector<int> *calo_row = nullptr;
    std::vector<int> *calo_type = nullptr;
    std::vector<int> *fcr = nullptr;
    std::vector<int> *timestamp = nullptr;
    int calo_nohits = 0;

    // Set branch addresses
    tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
    tree->SetBranchAddress("digicalo.waveform", &wave);
    tree->SetBranchAddress("digicalo.wall", &calo_wall);
    tree->SetBranchAddress("digicalo.side", &calo_side);
    tree->SetBranchAddress("digicalo.column", &calo_column);
    tree->SetBranchAddress("digicalo.row", &calo_row);
    tree->SetBranchAddress("digicalo.type", &calo_type);
    tree->SetBranchAddress("digicalo.fcr", &fcr);
    tree->SetBranchAddress("digicalo.timestamp", &timestamp);


    // Get the number of entries in the tree
    int max_entries = tree->GetEntries();
    std::cout << "Total entries in tree: " << max_entries << std::endl;

    // Create output ROOT file and TTree for baseline data
    TFile *outfile = new TFile(output_filename.c_str(), "RECREATE");
    TTree *baseline_tree = new TTree("baseline_tree", "OM baseline data");


    // Initialise output variables
    int event_num = -1;
    int om_num = -1;
    double width = 0;
    int maxConsec = 0;
    double baseline_delta = 0;

    // Create branches in the output tree
    baseline_tree->Branch("event_num", &event_num, "event_num/I");
    baseline_tree->Branch("om_num", &om_num, "om_num/I");
    baseline_tree->Branch("width", &width, "width/D");
    baseline_tree->Branch("maxConsec", &maxConsec, "maxConsec/I");
    baseline_tree->Branch("baseline_delta", &baseline_delta, "baseline_delta/D");

    for (int i = 0; i < max_entries; ++i){
        tree->GetEntry(i);
        event_num = i;
        
        for (int j = 0; j <calo_nohits; ++j){
            if (wave->at(j).size() < 1024) continue; // Skip if waveform is too short

            // --- calculate the om number
            int om_num_out = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, j);
            om_num = om_num_out;

            double baseline_orig = calculate_baseline976(wave->at(j));
            double baseline_96  = calculate_baseline96(wave->at(j));

            double stddev_orig = calculate_stddev976(wave->at(j), baseline_orig);

            // --- compute the difference between the baseline of 96 and baseline of 976
            baseline_delta = baseline_orig - baseline_96;
            
            // --- compute min, max, and noisy over bins 0–975 ---
            const std::vector<short>& wave_j = wave->at(j);
            auto [min_it, max_it] = std::minmax_element(wave_j.begin(), wave_j.begin() + 976);

            double min_val = static_cast<double>(*min_it);
            double max_val = static_cast<double>(*max_it);
            double noisy_val = max_val - min_val;
            width = noisy_val;

            // --- compute the number of consecutive bins outside the 2sigma band
            size_t consecutiveCountPlot = 0;
            size_t maxConsecutiveLocal = 0;

            for (size_t b = 0; b < 976; ++b) {

                if (std::abs(baseline_orig - wave_j[b]) > 2 * stddev_orig) {
                    consecutiveCountPlot++;
                    if (consecutiveCountPlot > maxConsecutiveLocal)
                        maxConsecutiveLocal = consecutiveCountPlot;
                } else {
                    consecutiveCountPlot = 0;
                }
            }

            maxConsec = static_cast<int>(maxConsecutiveLocal);


            baseline_tree->Fill();
        }
    }
    outfile->cd();
    baseline_tree->Write();
    outfile->Close();


    return 0;
}

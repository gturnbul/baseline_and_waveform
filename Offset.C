#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <memory>
#include <string>
#include <cmath>         // For math functions (with _USE_MATH_DEFINES)
#include <cstdlib>       // for rand()
#include <cassert>
#include <thread>
#include <map>
#include <unordered_map>
#include <tuple>
#include <complex>

// System / OS
#include <sys/stat.h>
#include <sys/types.h>

// ROOT Core
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TKey.h"
#include "TParameter.h"
#include <TString.h>
#include <TObjArray.h>

// ROOT Graphics / Canvas
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include <TPaletteAxis.h>
#include <TText.h>
#include <TLatex.h>
#include <TLegend.h>
#include "TMultiGraph.h"
#include <THStack.h>

// ROOT Histograms / Graphs / Fitting
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TSpectrum.h"

// ROOT Utilities
#include <TRandom3.h>
#include <TPaveText.h>
#include <TArrayF.h>
#include <TVirtualFFT.h>
#include <TComplex.h>

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

double calculate_baseline48(const std::vector<short>& waveform) {
    double baseline = 0;
    for (int i = 976; i < 1024; ++i) {
        baseline += waveform[i];
    }
    return baseline / 48;
}

// Function to calculate chi2df
double calculate_stddev976(const std::vector<short>& waveform, double baseline) {
    double sum_sq_diff = 0;
    for (int i = 0; i < 976; ++i) {
        sum_sq_diff += pow(waveform[i] - baseline, 2);
    }
    return sqrt(sum_sq_diff / 976);
}

double calculate_stddev48(const std::vector<short>& waveform, double baseline) {
    double sum_sq_diff = 0;
    for (int i = 976; i < 1024; ++i) {
        sum_sq_diff += pow(waveform[i] - baseline, 2);
    }
    return sqrt(sum_sq_diff / 48);
}
// // Function to calculate error on the mean (EOM), with sigma calculated inside
// double calculate_eom976(const std::vector<short>& waveform, double baseline) {
//     double sum_sq_diff = 0.0;

//     for (int i = 0; i < 976; ++i) {
//         double diff = waveform[i] - baseline;
//         sum_sq_diff += diff * diff;
//     }

//     double variance = sum_sq_diff / 976;      // Population variance
//     double sigma = std::sqrt(variance);     // Standard deviation
//     double eom = sigma / std::sqrt(976);      // Error on the mean

//     return eom;
// }
// double calculate_eom48(const std::vector<short>& waveform, double baseline) {
//     double sum_sq_diff = 0.0;

//     for (int i = 976; i < 1024; ++i) {
//         double diff = waveform[i] - baseline;
//         sum_sq_diff += diff * diff;
//     }

//     double variance = sum_sq_diff / 48;          // Population variance
//     double sigma = std::sqrt(variance);          // Standard deviation
//     double eom = sigma / std::sqrt(48);          // Error on the mean

//     return eom;
// }

struct Mod16Stats {
    std::vector<double> offsets;
    std::vector<double> sigmas;
    std::vector<double> chi2ndfs;
};

Mod16Stats compute_mod16_from_fits(const std::vector<double>& means,
                                    const std::vector<double>& sigmas,
                                    const std::vector<double>& chi2ndfs) {
    std::vector<double> mean_acc(16, 0.0);
    std::vector<double> sigma_acc(16, 0.0);
    std::vector<double> chi2ndf_acc(16, 0.0);
    std::vector<int> counts(16, 0);

    for (int i = 0; i < 976; ++i) {
        int mod16_index = i % 16;
        mean_acc[mod16_index] += means[i];
        sigma_acc[mod16_index] += sigmas[i];
        if (chi2ndfs[i] >= 0) {
            chi2ndf_acc[mod16_index] += chi2ndfs[i];
            counts[mod16_index]++;
            } else {
                static std::ofstream badfit_log("bad_fits_mod16_1278_width.csv", std::ios::app);
                badfit_log << "OM_UNKNOWN," << i << ",ndf=0" << std::endl;  
                // If you want OM numbers, pass it in as an argument
        }
    }

    for (int i = 0; i < 16; ++i) {
        if (counts[i] > 0) {
            mean_acc[i] /= counts[i];
            sigma_acc[i] /= counts[i];
            chi2ndf_acc[i] /= counts[i];
        }
    }

    return { mean_acc, sigma_acc, chi2ndf_acc };
}


// Get memory cell from bin number and FCR
int get_cell(int bin_no, int fcr) {
    return (bin_no + fcr) % 1024; // Now ranges from 0 to 1023
}

std::vector<short> reorder_waveform(const std::vector<short>& waveform, int fcr) {
    std::vector<short> reordered(1024, EMPTY_BIN);  // initialize all bins as empty

    int max_bins = std::min((int)waveform.size(), 976);  // only use bins 0 to 975

    for (int i = 0; i < max_bins; ++i) {
        int reordered_index = get_cell(i, fcr);
        reordered[reordered_index] = waveform[i];
    }

    // bins 976 to 1023 remain EMPTY_BIN, indicating gaps

    return reordered;
}
std::vector<short> inverse_reorder_waveform(const std::vector<short>& reordered, int fcr) {
    int max_bins = 976;  // original waveform length
    
    std::vector<short> original(max_bins, EMPTY_BIN);
    
    for (int i = 0; i < max_bins; ++i) {
        int reordered_index = get_cell(i, fcr);
        // Safety check: reordered_index should be within bounds
        if (reordered_index >= 0 && reordered_index < (int)reordered.size()) {
            original[i] = reordered[reordered_index];
        } else {
            // handle unexpected index (optional)
            original[i] = EMPTY_BIN;
        }
    }
    
    return original;
}
void saveWaveformAsPng(const std::vector<short>& wave,
                       const std::string& filename,
                       double baseline,
                       double stddev,
                       int om_num,
                       int event_num)
{
    // Create canvas
    TCanvas* c = new TCanvas("c", "Waveform", 800, 600);

    // Draw waveform
    TGraph* g = new TGraph(wave.size());
    for (size_t i = 0; i < wave.size(); ++i) {
        g->SetPoint(i, i, wave[i]);
    }
    g->SetTitle(Form("OM %d, Event %d;Bin;Amplitude", om_num, event_num));
    g->GetXaxis()->SetLimits(0, 1024); // full 1024 bins
    g->Draw("AL");

    double xMin = 0;
    double xMax = wave.size();

    // --- ±2σ band (green, semi-transparent) ---
    TBox* box2 = new TBox(xMin,
                          baseline - 2 * stddev,
                          xMax,
                          baseline + 2 * stddev);
    box2->SetFillColorAlpha(kGreen, 0.3);
    box2->Draw("same");

    // --- ±1σ band (yellow, more opaque) ---
    TBox* box1 = new TBox(xMin,
                          baseline - stddev,
                          xMax,
                          baseline + stddev);
    box1->SetFillColorAlpha(kYellow, 0.6);
    box1->Draw("same");

    // --- Baseline (red line) ---
    TLine* line = new TLine(xMin, baseline, xMax, baseline);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("same");

    // Redraw waveform on top
    g->Draw("L same");

    // Save to file
    c->SaveAs(filename.c_str());

    // Cleanup
    delete line;
    delete box1;
    delete box2;
    delete g;
    delete c;
}

// 2D histogram: X = OM number, Y = Bin number, Z = counts
int max_om = 712;        // number of OMs (from your QA loop)
int max_bins = 1024;     // waveform bins
TH2D *h2_bin_om = new TH2D("h2_bin_om", "Bin vs OM;Bin number;OM number;Counts",
                           max_bins, 0, max_bins,
                           max_om, 0, max_om);


//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Main ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <run_number>" << endl;
        return 1;
    }

    int run_number = stoi(argv[1]);

    string input_filename  = "snemo_run-" + to_string(run_number) + "_udd.root";
    string output_filename = "baseline_offset_calc_V3-" + to_string(run_number) + ".root";

    cout << "Reading:  " << input_filename  << endl;
    cout << "Writing:  " << output_filename << endl;

    // Open input file
    TFile* file = new TFile(input_filename.c_str(), "READ");
    TTree* tree = (TTree*)file->Get("SimData");

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
    TFile* outfile = new TFile(output_filename.c_str(), "RECREATE");
    TTree *baseline_tree = new TTree("baseline_tree", "OM baseline data");

    TTree *qaTree = new TTree("qaTree", "FEB/Wavecatcher QA");

    // Initialise output variables
    int event_num = -1;
    int om_num = -1;
    int calo_fcr = -1;
    double width = 0;
    // int maxConsec = 0;
    double baseline_orig = 0;
    double baseline_delta = 0;
    double baseline_diff = 0;
    double stddev_orig = 0;
    // double eom_orig = 0;
    double end_baseline = 0;
    double end_stddev = 0;
    // double end_eom = 0;
    
    double gradient = 0;
    double intercept = 0;
    double chi2ndf = 0;
    int peak_flag = -1;

    // Variables for the 1024 FEB offset correction
    double baseline_feb = 0;
    double stddev_feb = 0;
    // double eom_feb = 0;
  
    //first 976 variables after wavecatcher (mod16) correction
    double baseline_976 = 0;
    double stddev_976 = 0;
    // double eom_976 = 0;

    //end spike variables after end of waveform correction
    double baseline_end = 0;
    double stddev_end = 0;
    // double eom_end = 0;

    //Quality assurance variables
    int qa_om_num = -1;
    std::vector<double> mean_mem;
    std::vector<double> sigma_mem;
    std::vector<double> chi2ndf_mem; 
    
    std::vector<double> mean_timeo;
    std::vector<double> sigma_timeo;
    std::vector<double> chi2ndf_timeo;

    // Create branches in the output tree
    baseline_tree->Branch("event_num", &event_num, "event_num/I");
    baseline_tree->Branch("om_num", &om_num, "om_num/I");
    baseline_tree->Branch("calo_fcr", &calo_fcr, "calo_fcr/I");
    baseline_tree->Branch("width", &width, "width/D");
    // baseline_tree->Branch("maxConsec", &maxConsec, "maxConsec/I");
    baseline_tree->Branch("baseline", &baseline_orig, "baseline/D");
    baseline_tree->Branch("baseline_delta", &baseline_delta, "baseline_delta/D");
    baseline_tree->Branch("baseline_diff", &baseline_diff, "baseline_diff/D");
    baseline_tree->Branch("stddev", &stddev_orig, "stddev/D");
    // baseline_tree->Branch("eom", &eom_orig, "eom/D");
    baseline_tree->Branch("end_baseline", &end_baseline, "end_baseline/D");
    baseline_tree->Branch("end_stddev", &end_stddev, "end_stddev/D");
    // baseline_tree->Branch("end_eom", &end_eom, "end_eom/D");

    baseline_tree->Branch("gradient", &gradient, "gradient/D");
    baseline_tree->Branch("chi2ndf", &chi2ndf, "chi2ndf/D");
    baseline_tree->Branch("peak_flag", &peak_flag, "peak_flag/I");

    baseline_tree->Branch("baseline_feb", &baseline_feb, "baseline_feb/D");
    baseline_tree->Branch("stddev_feb", &stddev_feb, "stddev_feb/D");
    // baseline_tree->Branch("eom_feb", &eom_feb, "eom_feb/D");

    baseline_tree->Branch("baseline_976", &baseline_976, "baseline_976/D");
    baseline_tree->Branch("stddev_976", &stddev_976, "stddev_976/D");
    // baseline_tree->Branch("eom_976", &eom_976, "eom_976/D");

    baseline_tree->Branch("baseline_end", &baseline_end, "baseline_end/D");
    baseline_tree->Branch("stddev_end", &stddev_end, "stddev_end/D");
    // baseline_tree->Branch("eom_end", &eom_end, "eom_end/D");

    qaTree->Branch("om_num", &qa_om_num, "om_num/I");
    qaTree->Branch("mean_mem", &mean_mem);
    qaTree->Branch("sigma_mem", &sigma_mem);
    qaTree->Branch("chi2ndf_mem", &chi2ndf_mem);    
    
    qaTree->Branch("mean_timeo", &mean_timeo);
    qaTree->Branch("sigma_timeo", &sigma_timeo);
    qaTree->Branch("chi2ndf_timeo", &chi2ndf_timeo);

    //////////////////////////// MAPS ////////////////////////////////////////
    std::map<int, std::vector<std::vector<float>>> adc_values;
    std::map<int, std::vector<double>> mod16_offsets_per_om;
    std::map<int, std::vector<double>> mod16_sigmas_per_om;
    std::map<int, std::vector<double>> mod16_chi2ndf_per_om;
    std::map<int, std::vector<double>> mod16_data_per_om;

    std::map<int, std::vector<std::vector<float>>> feb_adc_values;
    std::map<int, std::vector<double>> feb_offsets_per_om;
    std::map<int, std::vector<double>> feb_sigmas_per_om;
    std::map<int, std::vector<double>> feb_chi2ndf_per_om;

    std::map<int, std::vector<double>> eow_offsets_per_om;
    std::map<int, std::vector<double>> eow_sigmas_per_om;
    std::map<int, std::vector<double>> eow_chi2ndf_per_om;
    std::map<int, std::vector<double>> eow_gausses_per_om;

    static std::map<int, std::vector<std::vector<double>>> om_waveforms;
    std:map<int, std::vector<double>> baseline_values_per_om;
    std::map<int, double> avg_baseline_per_om;
    // kept events per OM per FEB memory cell (0..1023).
    // We'll fill this by counting non-EMPTY_BIN entries in feb_adc_values
    std::map<int, std::vector<int>> kept_events_per_om_bin;


    //////////////////////////////////////////////////////////////////////////////////////
    ///////////////////// Finding and saving ADC_VAL range limit per OM //////////////////
    //////////////////////////////////////////////////////////////////////////////////////

    // Initialize a random generator
    TRandom3 *gRandom = new TRandom3(0); // Seed with 0 for different sequences each run
    int plot_count = 0;
    // int count_eom_cut = 0;
    int count_noisy_cut = 0;
    int savedKeptPngCount = 0;


    // --- Check for 10 consecutive bins > 2σ within first 976 bins ---
    size_t savedPngCount = 0;
    int count_peak_flag = 0;
    int count_chi2 = 0;
    int count_bld = 0;
    int cut_gradient = 0;

    for (int i = 0; i < max_entries; ++i){
            tree->GetEntry(i);
            event_num = i;
            //if (i < 10) continue; // Skip first 10 events due to cleanliness of early data

            for (int j = 0; j <calo_nohits; ++j){
                if (wave->at(j).size() < 1024) continue; // Skip if waveform is too short
                const vector<short>& w = wave->at(j);

                //calculate the om number
                int om_num_out = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, j);

                // Calculate original baseline, stddev, and eom
                baseline_orig = calculate_baseline976(wave->at(j));
                stddev_orig = calculate_stddev976(wave->at(j), baseline_orig);
                // eom_orig = calculate_eom976(wave->at(j), baseline_orig);
                baseline_delta = baseline_orig - calculate_baseline96(wave->at(j));

                // Linear fit
                int npts = 976;
                vector<double> x(npts), y(npts);
                for (int ii = 0; ii < npts; ++ii) { x[ii] = ii; y[ii] = w[ii]; }
                TGraph* gr_fit = new TGraph(npts, x.data(), y.data());
                TFitResultPtr fit_result = gr_fit->Fit("pol1", "QS");
                if (fit_result.Get() != nullptr) {
                    intercept = fit_result->Parameter(0);
                    gradient  = fit_result->Parameter(1);
                    chi2ndf   = fit_result->Chi2() / fit_result->Ndf();
                } else {
                    intercept = y[0];
                    gradient = 0;
                    chi2ndf = 0;
                }
                baseline_diff = baseline_orig - intercept;

                // Smooth waveform
                vector<double> smooth(976);
                for (int k = 0; k < 976; ++k) {
                    double sum = 0; int count = 0;
                    for (int s = -3; s <= 3; ++s) {
                        int idx = k + s;
                        if (idx >= 0 && idx < 976) { sum += w[idx]; count++; }
                    }
                    smooth[k] = sum / count;
                }

                // Pulse detection
                peak_flag = 0;
                int min_bins = 3, consec = 0, start_bin = -1, end_bin = -1;
                for (int k = 0; k < 976; ++k) {
                    double fit_k = intercept + gradient * k;
                    double deviation = smooth[k] - fit_k;
                    double slope = (k > 0 ? smooth[k] - smooth[k-1] : 0);
                    bool above = fabs(deviation) > 2 * stddev_orig;
                    bool nonflat = fabs(slope) > 0.25 * stddev_orig / 5;

                    if (above && nonflat) { consec++; if (consec >= min_bins && start_bin<0) start_bin = k - min_bins +1; }
                    else { if (start_bin>=0) { end_bin = k-1; break; } consec =0; }
                }
                if (start_bin>=0 && end_bin>start_bin) peak_flag = 1;


                const std::vector<short>& wave_j = wave->at(j);
                // Peak flag must be 0
                if (peak_flag != 0) continue;

                /// Fit must be good
                if (chi2ndf > 5.7) continue;

                // Baseline delta within ±0.6
                if (baseline_delta < -0.6 || baseline_delta > 0.6) continue;

                // Gradient essentially flat
                if (gradient < -0.001 || gradient > 0.001) continue;

            
            // --- Save a PNG for waveforms that are being kept and meet quality criteria ---
            // if (eom_orig> 0.078) {
                if (savedKeptPngCount < 10) {  // optional: limit number of saved PNGs
                TCanvas* c = new TCanvas("c_kept", "Kept Waveform", 800, 600);

                // Create the waveform graph
                TGraph* g = new TGraph(wave_j.size());
                for (size_t k = 0; k < wave_j.size(); ++k)
                    g->SetPoint(k, k, wave_j[k]);

                g->SetTitle(Form("OM %d, Event %d (Kept);Bin;Amplitude", om_num_out, event_num));
                g->SetLineColor(kBlack);
                g->SetLineWidth(1);
                g->GetXaxis()->SetLimits(0, 1024);
                g->GetXaxis()->SetNdivisions(510);
                g->GetXaxis()->SetNoExponent();
                g->Draw("AL");

                // ±1σ and ±2σ bands
                double xMin = 0, xMax = 1024;
                TBox* box3 = new TBox(xMin, baseline_orig - 3*stddev_orig, xMax, baseline_orig + 3*stddev_orig);
                box3->SetFillColorAlpha(kRed, 0.3);
                box3->Draw("same");

                TBox* box2 = new TBox(xMin, baseline_orig - 2*stddev_orig, xMax, baseline_orig + 2*stddev_orig);
                box2->SetFillColorAlpha(kGreen, 0.3);
                box2->Draw("same");

                TBox* box1 = new TBox(xMin, baseline_orig - stddev_orig, xMax, baseline_orig + stddev_orig);
                box1->SetFillColorAlpha(kYellow, 0.6);
                box1->Draw("same");

                // Baseline line
                TLine* line = new TLine(xMin, baseline_orig, xMax, baseline_orig);
                line->SetLineColor(kRed);
                line->SetLineWidth(2);
                line->Draw("same");

                // Redraw waveform on top
                g->Draw("L SAME");

                // Add annotation
                TLatex latex;
                latex.SetTextSize(0.04);
                // latex.DrawLatexNDC(0.65, 0.85, Form("EOM = %.4f", eom_orig));
                // latex.DrawLatexNDC(0.65, 0.78, Form("Width = %.2f", noisy_val));
                // latex.DrawLatexNDC(0.65, 0.78, Form("Consecutive = %zu", maxConsecutive));

                // Save PNG 
                c->SaveAs(Form("run_%d_kept_waveform_event%d_om%d_%d.png", run_number, event_num, om_num_out, savedKeptPngCount));

                // Clean up
                delete line;
                delete box1;
                delete box2;
                delete box3;
                delete g;
                delete c;

                savedKeptPngCount++;
                }
            // }


                int om_counter = 0;

                baseline_values_per_om[om_num_out].push_back(baseline_orig);
                //////////////////////////// Wavecatcher data /////////////////////////////////////
                // Initialise storage for the adc counts of bins per optical module
                if (adc_values.find(om_num_out) == adc_values.end()) {
                    adc_values[om_num_out] = std::vector<std::vector<float>>(1024, std::vector<float>());
                }
                
                // loop each bin in the waveform and collect the ADC values per bin
                for (int bin = 0; bin < 1024; ++bin) {
                    adc_values[om_num_out][bin].push_back(wave->at(j)[bin]);
                }

                ///////////////////////////// FEB data /////////////////////////////////////////////
                // initialise storage for the ADC counts of bins per om for FEB
                if (feb_adc_values.find(om_num_out) == feb_adc_values.end()) {
                    feb_adc_values[om_num_out] = std::vector<std::vector<float>>(1024, std::vector<float>());

                }
                std::vector<short> reordered_waveform = reorder_waveform(wave->at(j), fcr->at(j));
                for (int bin = 0; bin < 1024; ++bin) {
                    feb_adc_values[om_num_out][bin].push_back(reordered_waveform[bin]);
                }

            }
        
    }
  
    for (auto &pair : baseline_values_per_om) {
        double sum = 0;
        for (double v : pair.second) sum += v;
        avg_baseline_per_om[pair.first] = sum / pair.second.size();
    }
    
    // ----------------- compute kept-event counts from feb_adc_values -----------------
    for (const auto &om_pair : feb_adc_values) {
        int om = om_pair.first;
        const auto &bins = om_pair.second; // vector<vector<float>> size 1024

        // ensure a vector sized 1024, initialized to zero
        kept_events_per_om_bin[om] = std::vector<int>(1024, 0);

        for (int bin = 0; bin < 1024; ++bin) {
            int cnt = 0;
            for (float val : bins[bin]) {
                // only count real values (skip EMPTY_BIN sentinel)
                if (static_cast<short>(val) != EMPTY_BIN) ++cnt;
            }
            kept_events_per_om_bin[om][bin] = cnt;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////// FEB Calclulations /////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    std::ofstream csv_file_feb("run_" + std::to_string(run_number) + "_MemCell_offsets.csv");
    csv_file_feb << "om,bin,data_mean,fit_mean,sigma,chi2ndf,baseline_orig\n";

    // Process FEB_ADC_VALUES for histogramming
    for (auto& om_pair : feb_adc_values) {
        int om_num = om_pair.first;

        // if (om_num != 1) continue; // Skip all OMs except OM 1 for now -----------------
        // TCanvas *c = new TCanvas("c", "OM 1 Histograms", 800, 600); //------------------

        auto& adc_bins = om_pair.second;

        std::vector<double> fit_means_feb_om; // filled with 976 Gaussian means
        std::vector<double> fit_sigmas_feb_om; // filled with 976 Gaussian sigmas
        std::vector<double> fit_chi2_ndf_feb_om; // filled with 976 Gaussian chi2/ndf values
        std::vector<double> data_means_feb_om; // filled with 976 data means

        for (int bin = 0; bin < 1024; ++bin){
            std::vector<float>& adc_values_bin = adc_bins[bin];

            // Filter out EMPTY_BIN values
            std::vector<float> valid_adc_values;
            for (float val : adc_values_bin) {
                if (static_cast<short>(val) != EMPTY_BIN) {
                    valid_adc_values.push_back(val);
                }
            }

            if (valid_adc_values.empty()) {
                // No valid ADC values in this bin; skip histogramming/fitting
                continue;
            }

            // Calculate histogram range from valid ADC values only
            float min_adc = *std::min_element(valid_adc_values.begin(), valid_adc_values.end()) - 2;
            float max_adc = *std::max_element(valid_adc_values.begin(), valid_adc_values.end()) + 2;
            int nbins = static_cast<int>(std::ceil(max_adc - min_adc + 1));

            TString hist_name = Form("feb_om%d_bin%d", om_num, bin);
            TString hist_title = Form("FEB OM %d - Bin %d;ADC Value;Counts", om_num, bin);

            TH1D *hist = new TH1D(hist_name, hist_title, nbins, min_adc, max_adc);

            for (float adc_value : valid_adc_values) {
                hist->Fill(adc_value);
            }

            // extract data mean
            double data_mean = hist->GetMean();
            data_means_feb_om.push_back(data_mean);

            // Fit gaussian, extract mean/sigma
            TF1 *fit_func = new TF1("fit_func", "gaus", min_adc, max_adc);
            hist->Fit(fit_func, "Q", "", min_adc, max_adc);

            double mean = fit_func->GetParameter(1);
            double sigma = fit_func->GetParameter(2);
            double chi2 = fit_func->GetChisquare();
            int ndf = fit_func->GetNDF();
            double chi2_ndf = -1;
            if (ndf !=0){
                chi2_ndf = chi2 / ndf;
            } else {
                std::cout << "[BAD FIT] ndf=0 for OM=" << om_num
                        << " Bin=" << bin << std::endl;
            }
            

            fit_means_feb_om.push_back(mean);
            fit_sigmas_feb_om.push_back(sigma);
            fit_chi2_ndf_feb_om.push_back(chi2_ndf);

            csv_file_feb << om_num << ","
                        << bin << ","
                        << data_mean << ","
                        << mean << ","
                        << sigma << ","
                        << chi2_ndf << ","
                        << avg_baseline_per_om[om_num] << "\n";

            hist->Draw();
            // TString filename = Form("feb_bin%d_histogram.png", bin);
            // c->SaveAs(filename);

            delete hist;
            delete fit_func;
        }
        // Store the results for this OM
        feb_offsets_per_om[om_num] = fit_means_feb_om;
        feb_sigmas_per_om[om_num] = fit_sigmas_feb_om;
        feb_chi2ndf_per_om[om_num] = fit_chi2_ndf_feb_om;

        // delete c; // Delete the canvas after use
 
    }
    csv_file_feb.close();

    // //////////////////////////////////////////////////////////////////////////////////////
    // ////////////////////////////// Wavecater Calclulations ///////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////////

    // Process ADC_VALUES for histogramming, fitting and mod16 calculations
    for (auto& om_pair : adc_values) {
        int om_num = om_pair.first;
        auto& adc_bins = om_pair.second;

        std::vector<double> fit_means_om; // filled with 976 Gaussian means
        std::vector<double> fit_sigmas_om; // filled with 976 Gaussian sigmas
        std::vector<double> fit_chi2_ndf_om; // filled with 976 Gaussian chi2/ndf values
        std::vector<double> data_means_om; // filled with 976 data means

        for (int bin = 0; bin < 976; ++bin){
            std::vector<float>& adc_values_bin = adc_bins[bin];

            // Create and fill histogram
            TString hist_name = Form("om%d_bin%d", om_num, bin);
            TString hist_title = Form("OM %d - Bin %d;ADC Value;Counts", om_num, bin);
            // set histogram range based on the ADC values per om
            float min_adc = *std::min_element(adc_values_bin.begin(), adc_values_bin.end())-2;
            float max_adc = *std::max_element(adc_values_bin.begin(), adc_values_bin.end())+2;
            // Ensure bin count is an integer:
            int nbins = static_cast<int>(std::ceil(max_adc - min_adc + 1));
            // create histogram with nbins between min and max ADC values
            TH1D *hist = new TH1D(hist_name, hist_title, nbins, min_adc, max_adc);
            for (float adc_value : adc_values_bin) {
                hist->Fill(adc_value);
            }
            // extract data mean
            double data_mean = hist->GetMean();
            data_means_om.push_back(data_mean);

            // Fit gaussian, extract mean/sigma
            TF1 *fit_func = new TF1("fit_func", "gaus", min_adc, max_adc);
            hist->Fit(fit_func, "Q", "", min_adc, max_adc);
            double mean = fit_func->GetParameter(1);
            double sigma = fit_func->GetParameter(2);
            double chi2 = fit_func->GetChisquare();
            int ndf = fit_func->GetNDF();
            double chi2_ndf = -1;
            if (ndf !=0){
                chi2_ndf = chi2 / ndf;
            } else {
                std::cout << "[BAD FIT] ndf=0 for OM=" << om_num
                        << " Bin=" << bin << std::endl;
            }
            
            fit_means_om.push_back(mean);
            fit_sigmas_om.push_back(sigma);
            fit_chi2_ndf_om.push_back(chi2_ndf);

            delete hist;
            delete fit_func;

        }

        // Calculate mod16 offsets for this OM
        Mod16Stats stats = compute_mod16_from_fits(fit_means_om, fit_sigmas_om, fit_chi2_ndf_om);
        mod16_offsets_per_om[om_num] = stats.offsets;

        // calculate mod 16 data means for this OM
        mod16_data_per_om[om_num] = std::vector<double>(16, 0.0);   
        std::vector<int> counts(16, 0);
        for (int i = 0; i < 976; ++i) {
            int mod16_index = i % 16;
            mod16_data_per_om[om_num][mod16_index] += data_means_om[i];
            counts[mod16_index]++;
        }
        for (int i = 0; i < 16; ++i) {
            if (counts[i] > 0) {
                mod16_data_per_om[om_num][i] /= counts[i];
            }
        }

        // optionally store sigmas and chi2ndf too:
        mod16_sigmas_per_om[om_num] = stats.sigmas;
        mod16_chi2ndf_per_om[om_num] = stats.chi2ndfs;
    }
    // Open CSV file for writing
    std::ofstream csv_file("run_" + std::to_string(run_number) + "_Mod16_offsets.csv");

    // Write header
    csv_file << "om,offset_index,data_mean,offset_value,sigma_value,chi2ndf_value,baseline_orig\n";


    // Loop through the map and write data
    for (const auto& pair : mod16_offsets_per_om) {
        int om = pair.first;
        const std::vector<double>& offsets = pair.second;
        for (int i = 0; i < (int)offsets.size(); ++i) {
            csv_file << om << ","
                 << i << ","
                 << mod16_data_per_om[om][i] << ","
                 << mod16_offsets_per_om[om][i] << ","
                 << mod16_sigmas_per_om[om][i] << ","
                 << mod16_chi2ndf_per_om[om][i] << ","
                 << avg_baseline_per_om[om] << "\n";
        }
    }

    csv_file.close();
     // fit_file->Close();

    //////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  End spike offsets  /////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    std::ofstream spike_csv_file("run_" + std::to_string(run_number) + "_EoW_offsets.csv");
    spike_csv_file << "om,bin,offset_value,sigma_value,chi2ndf_value,baseline_orig\n";

    for (auto& om_pair : adc_values) {
    int om_num = om_pair.first;
    auto& adc_bins = om_pair.second;

        // vectors to store spike offsets for this OM
        std::vector<double> spike_means_om;
        std::vector<double> spike_sigmas_om;
        std::vector<double> spike_chi2ndf_om;
        std::vector<double> spike_gaus_om;

        auto it_mod16 = mod16_offsets_per_om.find(om_num);
        const std::vector<double>& mod16_offsets = (it_mod16 != mod16_offsets_per_om.end()) ? it_mod16->second : std::vector<double>{};


        // Loop bins 976-1023 
        for (int bin = 976; bin <= 1023; ++bin) {
            std::vector<float>& adc_values_bin = adc_bins[bin];

            TString hist_name = Form("om%d_bin%d", om_num, bin);
            TString hist_title = Form("OM %d - Bin %d;ADC Value;Counts", om_num, bin);
            float min_adc = *std::min_element(adc_values_bin.begin(), adc_values_bin.end()) - 5;
            float max_adc = *std::max_element(adc_values_bin.begin(), adc_values_bin.end()) + 5;
            int nbins = static_cast<int>(std::ceil(max_adc - min_adc + 1));
            TH1D* hist = new TH1D(hist_name, hist_title, nbins, min_adc, max_adc);
            for (float adc_value : adc_values_bin) {
                hist->Fill(adc_value);
            }

            TF1* fit_func = new TF1("fit_func", "gaus", min_adc, max_adc);
            hist->Fit(fit_func, "Q", "", min_adc, max_adc);
            double mean = fit_func->GetParameter(1);
            double sigma = fit_func->GetParameter(2);
            double chi2 = fit_func->GetChisquare();
            int ndf = fit_func->GetNDF();
            double chi2_ndf = -1;
              if (ndf !=0){
                chi2_ndf = chi2 / ndf;
            } else {
                std::cout << "[BAD FIT] ndf=0 for OM=" << om_num
                        << " Bin=" << bin << std::endl;
            }

            // Calculate difference: spike mean - mod16 offset for (bin % 16)
            int mod16_index = bin % 16;
            double mod16_value = (mod16_index < (int)mod16_offsets.size()) ? mod16_offsets[mod16_index] : 0.0;
            double corrected_spike_offset = mean - mod16_value;
            

            // Save results directly instead of computing mod16 offsets
            spike_gaus_om.push_back(mean);
            spike_means_om.push_back(corrected_spike_offset);
            spike_sigmas_om.push_back(sigma);
            spike_chi2ndf_om.push_back(chi2_ndf);

            delete hist;
            delete fit_func;
        }
        // Store the results for this OM
        eow_gausses_per_om[om_num] = spike_gaus_om;
        eow_offsets_per_om[om_num] = spike_means_om;
        eow_sigmas_per_om[om_num] = spike_sigmas_om;
        eow_chi2ndf_per_om[om_num] = spike_chi2ndf_om;


    }

    // Write to CSV

    for (const auto& om_pair : eow_offsets_per_om) {
        int om_num = om_pair.first;
        const auto& offsets = om_pair.second;
        const auto& sigmas = eow_sigmas_per_om[om_num];
        const auto& chi2s = eow_chi2ndf_per_om[om_num];

        for (int i = 0; i < (int)offsets.size(); ++i) {
            int bin = 976 + i;
            spike_csv_file << om_num << ","
                        << bin << ","
                        << offsets[i] << ","
                        << sigmas[i] << ","
                        << chi2s[i] << ","
                        << avg_baseline_per_om[om_num] << "\n";
        }
    }
    spike_csv_file.close();
    

    //////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Baseline Tree Fill /////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////

    int cut_count = 0;
    bool foundNonZero = false;

    for (int i = 0; i < max_entries; ++i){
        tree->GetEntry(i);
        event_num = i;
        //if (i < 10) continue; // Skip first 10 events due to cleanliness of early data

        for (int j = 0; j <calo_nohits; ++j){
            if (wave->at(j).size() < 1024) continue; // Skip if waveform is too short
            const vector<short>& w = wave->at(j);

            //calculate the om number
            om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, j);

            // Calculate original baseline, stddev, and eom
            baseline_orig = calculate_baseline976(wave->at(j));
            baseline_delta = baseline_orig - calculate_baseline96(wave->at(j));
            stddev_orig = calculate_stddev976(wave->at(j), baseline_orig);
            // eom_orig = calculate_eom976(wave->at(j), baseline_orig);
            end_baseline = calculate_baseline48(wave->at(j));
            end_stddev = calculate_stddev48(wave->at(j), end_baseline);
            // end_eom = calculate_eom48(wave->at(j), end_baseline);

            // Linear fit
            int npts = 976;
            vector<double> x(npts), y(npts);
            for (int ii = 0; ii < npts; ++ii) { x[ii] = ii; y[ii] = w[ii]; }
            TGraph* gr_fit = new TGraph(npts, x.data(), y.data());
            TFitResultPtr fit_result = gr_fit->Fit("pol1", "QS");
            if (fit_result.Get() != nullptr) {
                intercept = fit_result->Parameter(0);
                gradient  = fit_result->Parameter(1);
                chi2ndf   = fit_result->Chi2() / fit_result->Ndf();
            } else {
                intercept = y[0];
                gradient = 0;
                chi2ndf = 0;
            }
            baseline_diff = baseline_orig - intercept;

            // Smooth waveform
            vector<double> smooth(976);
            for (int k = 0; k < 976; ++k) {
                double sum = 0; int count = 0;
                for (int s = -3; s <= 3; ++s) {
                    int idx = k + s;
                    if (idx >= 0 && idx < 976) { sum += w[idx]; count++; }
                }
                smooth[k] = sum / count;
            }

            // Pulse detection
            peak_flag = 0;
            int min_bins = 3, consec = 0, startbin = -1, endbin = -1;
            for (int k = 0; k < 976; ++k) {
                double fit_k = intercept + gradient * k;
                double deviation = smooth[k] - fit_k;
                double slope = (k > 0 ? smooth[k] - smooth[k-1] : 0);
                bool above = fabs(deviation) > 2 * stddev_orig;
                bool nonflat = fabs(slope) > 0.25 * stddev_orig / 5;

                if (above && nonflat) { consec++; if (consec >= min_bins && startbin<0) startbin = k - min_bins +1; }
                else { if (startbin>=0) { endbin = k-1; break; } consec =0; }
            }
            if (startbin>=0 && endbin>startbin) peak_flag = 1;

            const std::vector<short>& wave_j = wave->at(j);

            // Peak flag must be 0
            if (peak_flag != 0) {
                count_peak_flag++;   // Count it
                continue;         // Skip it
            }

            // Fit must be good
            if (chi2ndf > 5.7) {
                count_chi2++;   // Count it
                continue;         // Skip it
            }
        
            if (baseline_delta < -0.6 || baseline_delta > 0.6) {
                count_bld++;   // Count this event as rejected
                continue;         
            }

            if (gradient < -0.001 || gradient > 0.001) {
                cut_gradient++;      // event rejected
                continue;
            }


            // feb waveform correction //////////////////////////////////////////////////////////////
            // Step 1: reorder the original waveform
            std::vector<short> reordered_wave = reorder_waveform(wave->at(j), fcr->at(j));
            
            
            // Step 2: find the FEB offsets for this OM and subtract them from reordered waveform
            auto it_feb = feb_offsets_per_om.find(om_num);
            if (it_feb != feb_offsets_per_om.end()) {
                const std::vector<double>& feb_offsets = it_feb->second;
                for (int bin = 0; bin < 1024; ++bin) {
                    reordered_wave[bin] -= feb_offsets[bin];
                }
            } else {
                std::cerr << "Warning: No FEB offsets found for OM " << om_num << std::endl;
            }

            // Step 3: invert reorder to get corrected waveform back in original order
            std::vector<short> corrected_wave = inverse_reorder_waveform(reordered_wave, fcr->at(j));

            for (int bin = 0; bin < 1024; ++bin) {
                if (corrected_wave[bin] != EMPTY_BIN) {
                    h2_bin_om->Fill(bin, om_num);
                    h2_bin_om->SetMinimum(300);  // lower edge of color scale
                    h2_bin_om->SetMaximum(1000); // upper edge of color scale
                }
            }

            //////////////////////////////// print comparison of before/after FEB correction
            static int saved_feb_examples = 0;
            if (saved_feb_examples < 5) {
                std::vector<double> corrected_with_baseline(977);
                for (int bin = 0; bin <= 976; ++bin)
                    corrected_with_baseline[bin] = corrected_wave[bin] + baseline_orig;

                TCanvas* c = new TCanvas("c_feb", "FEB Correction Comparison", 1000, 500);

                // Original waveform
                TGraph* g_before = new TGraph(977);
                for (int bin = 0; bin <= 976; ++bin)
                    g_before->SetPoint(bin, bin, wave->at(j)[bin]);
                g_before->SetTitle(Form("OM %d Event %d: Before (Blue) & After (Red) FEB Correction;Bin;ADC", om_num, event_num));
                g_before->GetXaxis()->SetRangeUser(0, 976);
                g_before->SetLineColor(kBlue);
                g_before->Draw("AL");

                // Corrected waveform
                TGraph* g_after = new TGraph(977);
                for (int bin = 0; bin <= 976; ++bin)
                    g_after->SetPoint(bin, bin, corrected_with_baseline[bin]);
                g_after->SetLineColor(kRed);
                g_after->Draw("L SAME");

                c->SaveAs(Form("run_%d_before_after_event%d_om%d_%d.png", run_number, event_num, om_num, saved_feb_examples));
                
                delete c;
                saved_feb_examples++;
            }

            //////////////////////////////////// print zoomed in correction
            static int saved_feb_zoom_examples = 0;
            if (saved_feb_zoom_examples < 5) {
                const int n_bins = 49; // bins 0-48 inclusive
                std::vector<double> corrected_with_baseline(n_bins);
                for (int bin = 0; bin < n_bins; ++bin)
                    corrected_with_baseline[bin] = corrected_wave[bin] + baseline_orig;

                TCanvas* c = new TCanvas("c_feb", "FEB Correction Comparison (Zoomed 0-48)", 1000, 500);

                // Original waveform (bins 0-48)
                TGraph* g_before = new TGraph(n_bins);
                for (int bin = 0; bin < n_bins; ++bin)
                    g_before->SetPoint(bin, bin, wave->at(j)[bin]);
                g_before->SetTitle(Form("OM %d Event %d: Zoomed 0-48 Bins;Bin;ADC", om_num, event_num));
                g_before->SetLineColor(kBlue);
                g_before->Draw("AL");

                // Corrected waveform overlay
                TGraph* g_after = new TGraph(n_bins);
                for (int bin = 0; bin < n_bins; ++bin)
                    g_after->SetPoint(bin, bin, corrected_with_baseline[bin]);
                g_after->SetLineColor(kRed);
                g_after->Draw("L SAME");

                // Set X-axis range explicitly 0-48
                g_before->GetXaxis()->SetLimits(0, n_bins-1);

                // Optional: add legend
                TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
                leg->AddEntry(g_before, "Original", "l");
                leg->AddEntry(g_after, "Corrected", "l");
                leg->Draw();

                // Save
                c->SaveAs(Form("run_%d_MemCell_zoomed(0-48)%d_OM%d.png", run_number, event_num, om_num));
              
                delete c;
                saved_feb_zoom_examples++;
            }


            ////////////////////////////////////////////////////////////////////////

            // Step 4: calculate baseline, stddev, and eom on corrected waveform in original order
            baseline_feb = calculate_baseline976(corrected_wave);
            stddev_feb   = calculate_stddev976(corrected_wave, baseline_feb);
            // eom_feb      = calculate_eom976(corrected_wave, baseline_feb);

            // 976 waveform correction ///////////////////////////////////////////////////////////////

            // --- Save original waveform BEFORE applying correction ---
            std::vector<double> original_waveform_bins(wave->at(j).begin(), wave->at(j).begin() + 976);

            // --- Apply Mod16 correction ---
            auto it = mod16_offsets_per_om.find(om_num);
            if (it != mod16_offsets_per_om.end()) {
                const std::vector<double>& offsets = it->second;

                for (int bin = 0; bin < 976; ++bin) {
                    wave->at(j)[bin] -= offsets[bin % 16]; // Apply mod16 offset
                }
            } else {
                // OM not found in offsets map, handle error or skip correction
                std::cerr << "Warning: No mod16 offsets found for OM " << om_num << std::endl;
            }

            ////////////////////////////////////////// Mod16 comparison ///////////////////////////////////////////

            static int saved_mod16_examples = 0;
            if (saved_mod16_examples < 5) {
                std::vector<double> corrected_with_baseline(976);
                for (int bin = 0; bin < 976; ++bin)
                    corrected_with_baseline[bin] = wave->at(j)[bin] + baseline_orig;

                // Create single canvas
                TCanvas* c = new TCanvas("c_mod16", "Mod16 Correction Comparison", 800, 600);

                // Create blue (original) waveform
                TGraph* g_before = new TGraph(976);
                for (int bin = 0; bin < 976; ++bin)
                    g_before->SetPoint(bin, bin, original_waveform_bins[bin]);
                g_before->SetTitle(Form("OM %d Event %d Before/After Mod16;Bin;ADC", om_num, event_num));
                g_before->SetLineColor(kBlue);
                g_before->Draw("AL");  // Draw axes + blue line

                // Force x-axis range 0–976
                g_before->GetXaxis()->SetLimits(0, 976);   // sets logical x range
                g_before->GetHistogram()->GetXaxis()->SetRangeUser(0, 976);  // alternative way

                // Create red (corrected) waveform
                TGraph* g_after = new TGraph(976);
                for (int bin = 0; bin < 976; ++bin)
                    g_after->SetPoint(bin, bin, corrected_with_baseline[bin]);
                g_after->SetLineColor(kRed);
                g_after->Draw("L SAME");  // Draw red line on top of blue

                // Add legend
                auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
                legend->AddEntry(g_before, "Before Correction", "l");
                legend->AddEntry(g_after, "After Correction", "l");
                legend->Draw();

                // Save and clean up
                c->SaveAs(Form("run_%d_Mod16_overlay_event%d_OM%d.png", run_number, event_num, om_num));
                delete c;
                saved_mod16_examples++;
            }

            static int saved_mod16_zoom_examples = 0;
            if (saved_mod16_zoom_examples < 5) {
                std::vector<double> corrected_with_baseline(976);
                for (int bin = 0; bin < 976; ++bin)
                    corrected_with_baseline[bin] = wave->at(j)[bin] + baseline_orig;

                // Create single canvas
                TCanvas* c = new TCanvas("c_mod16", "Mod16 Correction Comparison", 800, 600);

                // Create blue (original) waveform
                TGraph* g_before = new TGraph(976);
                for (int bin = 0; bin < 976; ++bin)
                    g_before->SetPoint(bin, bin, original_waveform_bins[bin]);
                g_before->SetTitle(Form("OM %d Event %d Before/After Mod16;Bin;ADC", om_num, event_num));
                g_before->SetLineColor(kBlue);
                g_before->Draw("AL");  // Draw axes + blue line

                // Force x-axis range 0–976
                g_before->GetXaxis()->SetLimits(0, 48);   // sets logical x range
                g_before->GetHistogram()->GetXaxis()->SetRangeUser(0, 48);  // alternative way

                // Create red (corrected) waveform
                TGraph* g_after = new TGraph(976);
                for (int bin = 0; bin < 976; ++bin)
                    g_after->SetPoint(bin, bin, corrected_with_baseline[bin]);
                g_after->SetLineColor(kRed);
                g_after->Draw("L SAME");  // Draw red line on top of blue

                // Add legend
                auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
                legend->AddEntry(g_before, "Before Correction", "l");
                legend->AddEntry(g_after, "After Correction", "l");
                legend->Draw();

                // Save and clean up
                c->SaveAs(Form("run_%d_Mod16_overlay_zoom_event%d_OM%d.png", run_number, event_num, om_num));
                delete c;
                saved_mod16_zoom_examples++;
            }

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////


            // --- Calculate adjusted baseline, stddev, and eom ---
            baseline_976 = calculate_baseline976(wave->at(j));
            stddev_976 = calculate_stddev976(wave->at(j), baseline_976);
            // eom_976 = calculate_eom976(wave->at(j), baseline_976);

         
            // Spike correction for bins 976–1023 /////////////////////////////////////////////////
            static int saved_wave_snippets = 0;

            int start_bin = 976;
            int end_bin = 1023;
            int n_bins = end_bin - start_bin + 1;

            // Copy original bins BEFORE correction
            std::vector<double> original_spike_bins(wave->at(j).begin() + start_bin, wave->at(j).begin() + end_bin + 1);

            // Apply spike + mod16 correction
            auto it_spike = eow_offsets_per_om.find(om_num);

            if (it_spike != eow_offsets_per_om.end()) {
                const std::vector<double>& spike_offsets = it_spike->second;

                for (int bin = start_bin; bin <= end_bin; ++bin) {
                    int spike_bin_index = bin - start_bin;
                    if (spike_bin_index < (int)spike_offsets.size()) {
                        wave->at(j)[bin] -= spike_offsets[spike_bin_index]; // remove spike offset

                        // Apply mod16 offset to spike-corrected bin
                        auto it_mod16 = mod16_offsets_per_om.find(om_num);  
                        if (it_mod16 != mod16_offsets_per_om.end()) {
                            const std::vector<double>& offsets = it_mod16->second;
                            wave->at(j)[bin] -= offsets[bin % 16]; // apply mod16 offset
                        }
                    }
                }
            } else {
                std::cerr << "Warning: No spike offsets found for OM " << om_num << std::endl;
            }

            ///////////////////////////////////////////// Spike correction before/after plot ////////////////////////////////////////////
            static int saved_spike_examples = 0;
            if (saved_spike_examples < 5) {
                int start_bin = 976;
                int end_bin = 1023;
                int n_bins = end_bin - start_bin + 1;

                // Add baseline back to corrected waveform so it's on the same scale
                std::vector<double> corrected_with_baseline(n_bins);
                for (int i = 0; i < n_bins; ++i)
                    corrected_with_baseline[i] = wave->at(j)[start_bin + i] + baseline_orig;

                // Create one canvas (no divide)
                TCanvas* c = new TCanvas("c_spike", "End Spike Correction Overlay", 800, 600);

                // Create both graphs
                TGraph* g_before = new TGraph(n_bins);
                TGraph* g_after  = new TGraph(n_bins);

                // Fill points
                for (int i = 0; i < n_bins; ++i) {
                    g_before->SetPoint(i, start_bin + i, original_spike_bins[i]);
                    g_after->SetPoint(i, start_bin + i, corrected_with_baseline[i]);
                }

                // Original waveform (blue)
                g_before->SetLineColor(kBlue);
                g_before->SetLineWidth(2);
                g_before->GetXaxis()->SetRangeUser(976, 1023);
                g_before->SetTitle(Form("End Spike Correction OM %d, Event %d (bins 976-1023);Bin;ADC", om_num, event_num));

                // Draw first
                g_before->Draw("AL");

                // Corrected waveform (red) overlayed on same axes
                g_after->SetLineColor(kRed);
                g_after->SetLineWidth(2);
                g_after->Draw("L SAME");

                // Optional: keep consistent y-axis range
                g_before->GetYaxis()->SetRangeUser(baseline_orig - 40, baseline_orig + 40);

                // Add a legend
                auto legend = new TLegend(0.65, 0.8, 0.9, 0.9);
                legend->AddEntry(g_before, "Original", "l");
                legend->AddEntry(g_after, "Corrected (End Spike)", "l");
                legend->Draw();

                // Save and clean up
                c->SaveAs(Form("run_%d_Spike_overlay_event%d_OM%d.png", run_number, event_num, om_num));
                delete c;

                saved_spike_examples++;
            }


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Recalculate baseline, stddev, and eom
            baseline_end = calculate_baseline48(wave->at(j));
            stddev_end   = calculate_stddev48(wave->at(j), baseline_end);
            // eom_end      = calculate_eom48(wave->at(j), baseline_end);

            calo_fcr = fcr->at(j);
            
            ////////////////////////// Fill the Baseline Tree ////////////////////////////
            baseline_tree->Fill();
        } 
    }
    std::cout << "Number of waveforms cut from peak flag: " << count_peak_flag << std::endl;
    std::cout << "Number of waveforms cut from goodness of fit: " << count_chi2 << std::endl;
    std::cout << "Number of waveforms cut from baseline differences: " << count_bld << std::endl;
    std::cout << "Number of waveforms cut from the gradient: " << cut_gradient << std::endl;

    // --- Mark empty bins with sentinel value (-1) ---
    for (int x = 1; x <= h2_bin_om->GetNbinsX(); ++x) {
        for (int y = 1; y <= h2_bin_om->GetNbinsY(); ++y) {
            if (h2_bin_om->GetBinContent(x, y) == 0) {
                h2_bin_om->SetBinContent(x, y, -1);
            }
        }
    }

    // --- Now draw ---
   TCanvas *c2 = new TCanvas("c2", "OM vs Bin", 1024, 712);

    gStyle->SetPalette(kRainBow);  // set rainbow first
    int ncolors = TColor::GetNumberOfColors();
    std::vector<int> revColors(ncolors);
    for (int i = 0; i < ncolors; ++i) {
        revColors[i] = TColor::GetColorPalette(ncolors - i - 1);
    }
    gStyle->SetPalette(ncolors, revColors.data());  // reversed rainbow

    h2_bin_om->SetStats(0);
    h2_bin_om->SetMinimum(-1);   // include sentinel in color scale
    h2_bin_om->Draw("COLZ");

    c2->SaveAs(Form("run_%d_OM_vs_clock_tick.png", run_number));
    // if (!foundNonZero) {
    //     std::cout << "No non-zero mod64 values were found." << std::endl;
    // }
    // storing quality of fit assurance data ////////////////////////////////////////////////////
    for (int omnum = 0; omnum < 712; ++omnum) {
        qa_om_num = omnum;
        // Clear previous event's data
        mean_timeo.clear();
        sigma_timeo.clear();
        chi2ndf_timeo.clear();

        auto it_offset = mod16_offsets_per_om.find(omnum);
        auto it_sigma  = mod16_sigmas_per_om.find(omnum);
        auto it_chi2   = mod16_chi2ndf_per_om.find(omnum);

        if (it_offset != mod16_offsets_per_om.end() &&
            it_sigma != mod16_sigmas_per_om.end() &&
            it_chi2 != mod16_chi2ndf_per_om.end()) {

            mean_timeo = it_offset->second;
            sigma_timeo = it_sigma->second;
            chi2ndf_timeo = it_chi2->second;

        } else {
            std::cerr << "Warning: Missing QA data for OM " << om_num << std::endl;

            // Optionally fill with NaNs or -999 if QA data missing
            mean_timeo.assign(16, -999);
            sigma_timeo.assign(16, -999);
            chi2ndf_timeo.assign(16, -999);
        }

        // feb QA ///////////////////////////////////////////////////////////////
        mean_mem.clear();
        sigma_mem.clear();
        chi2ndf_mem.clear();

        auto it_feb_offset = feb_offsets_per_om.find(omnum);
        auto it_feb_sigma  = feb_sigmas_per_om.find(omnum);
        auto it_feb_chi2   = feb_chi2ndf_per_om.find(omnum);

        if (it_feb_offset != feb_offsets_per_om.end() &&
            it_feb_sigma != feb_sigmas_per_om.end() &&
            it_feb_chi2 != feb_chi2ndf_per_om.end()) {

            mean_mem = it_feb_offset->second;
            sigma_mem = it_feb_sigma->second;
            chi2ndf_mem = it_feb_chi2->second;

        } else {
            std::cerr << "Warning: Missing FEB QA data for OM " << omnum << std::endl;

            // Optionally fill with NaNs or -999 if QA data missing
            mean_mem.assign(1024, -999);
            sigma_mem.assign(1024, -999);
            chi2ndf_mem.assign(1024, -999);

        }
        qaTree->Fill();
    }   

    // ===================== SAVE KEPT EVENTS PER OM + MEMORY CELL =====================
    std::ofstream kept_feb_csv( "run_" + std::to_string(run_number) + "_Kept_events_per_om_bin.csv");
    kept_feb_csv << "om,bin,kept_events\n";

    for (const auto &om_pair : kept_events_per_om_bin) {
        int om = om_pair.first;
        const std::vector<int> &vec = om_pair.second;
        for (int bin = 0; bin < (int)vec.size(); ++bin) {
            if (vec[bin] == 0) continue; // optional: skip zeros
            kept_feb_csv << om << "," << bin << "," << vec[bin] << "\n";
        }
    }
    kept_feb_csv.close();
    std::cout << "Saved kept-event counts to kept_events_per_om_bin.csv" << std::endl;


    //////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// End of code, writes and closes ////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // write and close
    outfile->cd();
    baseline_tree->Write();
    qaTree->Write();
    outfile->Close();
    file->Close();

    return 0;
}

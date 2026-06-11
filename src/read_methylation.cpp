#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <zlib.h>
using namespace Rcpp;

// Helper function to read gzipped files
std::vector<std::string> readGzippedLines(const char* filename) {
  std::vector<std::string> lines;
  gzFile file = gzopen(filename, "rb");
  
  if (!file) {
    Rcpp::stop("Could not open file: %s", filename);
  }
  
  char buffer[4096];
  while (gzgets(file, buffer, sizeof(buffer)) != NULL) {
    lines.push_back(std::string(buffer));
  }
  
  gzclose(file);
  return lines;
}

// [[Rcpp::export]]
DataFrame readMethylationFile(std::string filename, int min_coverage = 0) {
  Rcout << "Reading file: " << filename << "\n";
  // Read all lines from the gzipped file
  std::vector<std::string> lines = readGzippedLines(filename.c_str());

  // Use std::vector with push_back so we only store rows that pass the filter.
  std::vector<std::string> chr_vec, strand_vec, site_vec;
  std::vector<int>         pos_vec, mc_vec, cov_vec;

  // Reserve capacity based on total lines to avoid repeated reallocation.
  // Most rows will pass when min_coverage is 0; this is a safe upper bound.
  int n_lines = lines.size();
  chr_vec.reserve(n_lines);
  pos_vec.reserve(n_lines);
  strand_vec.reserve(n_lines);
  site_vec.reserve(n_lines);
  mc_vec.reserve(n_lines);
  cov_vec.reserve(n_lines);
  
  // Parse each line
  for (int i = 0; i < n_lines; i++) {
    std::string line = lines[i];
    std::istringstream iss(line);
    std::string chr_val, strand_val, site_val;
    int pos_val, mc_val, cov_val;
    
    // Read tab-separated values
    if (std::getline(iss, chr_val, '\t') && 
        (iss >> pos_val) && iss.ignore(1, '\t') &&
        std::getline(iss, strand_val, '\t') &&
        std::getline(iss, site_val, '\t') &&
        (iss >> mc_val) && iss.ignore() &&
        (iss >> cov_val)) {
        
      // Apply minimum coverage filter before storing
      if (cov_val < min_coverage) continue;
        chr_vec.push_back(chr_val);
        pos_vec.push_back(pos_val);
        strand_vec.push_back(strand_val);
        site_vec.push_back(site_val);
        mc_vec.push_back(mc_val);
        cov_vec.push_back(cov_val);
      }
    }
    if (min_coverage > 0) {
    Rcout << "  Retained " << chr_vec.size() << " / " << n_lines
          << " sites (cov >= " << min_coverage << ")\n";
  }
  
  DataFrame df = DataFrame::create(
    _["chr"]    = chr_vec,
    _["pos"]    = pos_vec,
    _["strand"] = strand_vec,
    _["site"]   = site_vec,
    _["mc"]     = mc_vec,
    _["cov"]    = cov_vec
  );
  
  return df;
}


// [[Rcpp::export]]
List readMethylationFiles(CharacterVector filenames, int min_coverage = 0) {
  Rcout << "Entering readMethylationFiles with " << filenames.size() << " files\n";
  int n_files = filenames.size();
  List results(n_files);
  
  for (int i = 0; i < n_files; i++) {
    std::string filename = as<std::string>(filenames[i]);
    Rcout << "Processing file: " << filename << "\n";
    Rcpp::checkUserInterrupt(); // Allow user to cancel
    results[i] = readMethylationFile(filename, min_coverage);
  }
  
  Rcout << "Exiting readMethylationFiles\n";
  return results;
} 

// [[Rcpp::export]]
std::string testCpp() {
  Rcout << "Test C++ function called successfully\nwhy this so painfull";
  return "Success";
}



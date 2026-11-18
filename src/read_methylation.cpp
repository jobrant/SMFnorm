#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <zlib.h>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <vector>
using namespace Rcpp;

// Helper function to read gzipped (or plain text) files. Stops if the file
// cannot be opened or the stream ends in an error (e.g. a truncated gzip file).
std::vector<std::string> readGzippedLines(const char* filename) {
  std::vector<std::string> lines;
  gzFile file = gzopen(filename, "rb");

  if (!file) {
    Rcpp::stop("Could not open file: %s", filename);
  }

  // A line longer than the buffer arrives in pieces; join them.
  char buffer[65536];
  std::string current;
  bool pending = false;
  while (gzgets(file, buffer, sizeof(buffer)) != NULL) {
    current += buffer;
    pending = true;
    if (!current.empty() && current[current.size() - 1] == '\n') {
      lines.push_back(current);
      current.clear();
      pending = false;
    }
  }
  if (pending) lines.push_back(current);   // last line without a newline

  int errnum = Z_OK;
  const char* msg = gzerror(file, &errnum);
  if (errnum != Z_OK && errnum != Z_STREAM_END) {
    std::string detail(msg ? msg : "unknown error");
    gzclose(file);
    Rcpp::stop("Error reading %s: %s (the file may be truncated or corrupt).",
               filename, detail.c_str());
  }

  gzclose(file);
  return lines;
}

static bool isBlank(const std::string& s) {
  return s.find_first_not_of(" \t\r\n") == std::string::npos;
}

// Parse a whole field as an integer; false if it is empty or has other text.
static bool parseInt(const std::string& field, int& out) {
  if (field.empty()) return false;
  size_t i = (field[0] == '+' || field[0] == '-') ? 1 : 0;
  if (i == field.size()) return false;
  for (size_t k = i; k < field.size(); k++) {
    if (field[k] < '0' || field[k] > '9') return false;
  }
  char* end = NULL;
  errno = 0;
  long v = std::strtol(field.c_str(), &end, 10);
  if (errno != 0 || v > INT_MAX || v < INT_MIN) return false;
  out = static_cast<int>(v);
  return true;
}

// Reads records of chr, pos, strand, site, mc, cov (tab-separated; further
// columns are ignored). A first line that is not a record is taken as a
// header and skipped, and blank lines are skipped. Any other line that is not
// a valid record stops with an error, as does a file with no records: a
// damaged file must not load as a smaller or empty sample.
// The R reader (.load_single_file_r) follows the same rules.
// [[Rcpp::export]]
DataFrame readMethylationFile(std::string filename, int min_coverage = 0) {
  std::vector<std::string> lines = readGzippedLines(filename.c_str());

  // Use std::vector with push_back so we only store rows that pass the filter.
  std::vector<std::string> chr_vec, strand_vec, site_vec;
  std::vector<int>         pos_vec, mc_vec, cov_vec;

  // Reserve capacity based on total lines to avoid repeated reallocation.
  size_t n_lines = lines.size();
  chr_vec.reserve(n_lines);
  pos_vec.reserve(n_lines);
  strand_vec.reserve(n_lines);
  site_vec.reserve(n_lines);
  mc_vec.reserve(n_lines);
  cov_vec.reserve(n_lines);

  bool seen_content = false;   // has a non-blank line been seen yet?
  size_t n_records = 0;        // valid records, before the coverage filter

  for (size_t i = 0; i < n_lines; i++) {
    std::string& line = lines[i];
    while (!line.empty() && (line[line.size() - 1] == '\n' ||
                             line[line.size() - 1] == '\r')) {
      line.erase(line.size() - 1);
    }
    if (isBlank(line)) continue;
    bool first = !seen_content;
    seen_content = true;

    // Split the first six tab-separated fields
    std::string f[6];
    size_t start = 0;
    int n_fields = 0;
    while (n_fields < 6) {
      size_t tab = line.find('\t', start);
      if (tab == std::string::npos) {
        f[n_fields++] = line.substr(start);
        break;
      }
      f[n_fields++] = line.substr(start, tab - start);
      start = tab + 1;
    }

    int pos_val = 0, mc_val = 0, cov_val = 0;
    bool ok = n_fields == 6 && parseInt(f[1], pos_val) &&
              parseInt(f[4], mc_val) && parseInt(f[5], cov_val);
    if (!ok) {
      if (first) continue;   // header line
      Rcpp::stop("%s: line %d is not a valid record (expected chr, pos, "
                 "strand, site, mc, cov separated by tabs).",
                 filename.c_str(), static_cast<int>(i + 1));
    }
    if (mc_val < 0 || cov_val < 0 || mc_val > cov_val) {
      Rcpp::stop("%s: line %d has impossible counts (mc = %d, cov = %d).",
                 filename.c_str(), static_cast<int>(i + 1), mc_val, cov_val);
    }
    n_records++;

    // Apply minimum coverage filter before storing
    if (cov_val < min_coverage) continue;
    chr_vec.push_back(f[0]);
    pos_vec.push_back(pos_val);
    strand_vec.push_back(f[2]);
    site_vec.push_back(f[3]);
    mc_vec.push_back(mc_val);
    cov_vec.push_back(cov_val);

    if ((i & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
  }

  if (n_records == 0) {
    Rcpp::stop("No records found in %s.", filename.c_str());
  }

  DataFrame df = DataFrame::create(
    _["chr"]    = chr_vec,
    _["pos"]    = pos_vec,
    _["strand"] = strand_vec,
    _["site"]   = site_vec,
    _["mc"]     = mc_vec,
    _["cov"]    = cov_vec,
    _["stringsAsFactors"] = false
  );

  return df;
}


// [[Rcpp::export]]
List readMethylationFiles(CharacterVector filenames, int min_coverage = 0) {
  int n_files = filenames.size();
  List results(n_files);

  for (int i = 0; i < n_files; i++) {
    std::string filename = as<std::string>(filenames[i]);
    Rcpp::checkUserInterrupt(); // Allow user to cancel
    results[i] = readMethylationFile(filename, min_coverage);
  }

  return results;
}

// [[Rcpp::export]]
std::string testCpp() {
  Rcout << "Test C++ function called successfully\nwhy this so painfull";
  return "Success";
}



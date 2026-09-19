// Kplex: convenience wrapper for FastK + Histex to produce k-plexity curves
// Usage:
//   Kplex -i genome.fa -k 5:151:1 -h 1:1000 [-T threads] [-P jobs] [-M mem_gb] [-o prefix] [-c out.csv] [-K]
//
// Produces CSV with columns: k,unique_kmers,total_kmers,fraction_unique
//
// Each k is an independent FastK+Histex computation. Kplex auto-detects the
// number of CPU cores and available RAM and runs several k values concurrently
// (a process pool), each FastK using a share of the threads. Output is identical
// to the sequential version; only wall-clock time changes.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>

static void usage(const char *prog)
{
  fprintf(stderr,
          "Usage: %s -i genome.fa -k start:end:step -h h1:h2 [-T threads] [-P jobs] [-M mem_gb] [-o prefix] [-c out.csv] [-K]\n"
          "  -i  input FASTA\n"
          "  -k  k range, e.g. 5:151:1\n"
          "  -h  Histex range, e.g. 1:1000 (passed as -h1:1000)\n"
          "  -T  total thread budget (default: auto = number of cores)\n"
          "  -P  parallel FastK jobs to run at once (default: auto from cores & RAM)\n"
          "  -M  memory budget in GB (default: auto = available RAM)\n"
          "  -o  prefix for FastK output (default: basename of fasta)\n"
          "  -c  output CSV (default: <prefix>.kplex.csv)\n"
          "  -K  keep intermediate .hist and histex output files\n", prog);
}

static int parse_k_range(const char *s, int *k0, int *k1, int *step)
{
  int a=0,b=0,c=0;
  if (sscanf(s, "%d:%d:%d", &a, &b, &c) == 3 ||
      sscanf(s, "%d-%d:%d", &a, &b, &c) == 3)
  {
    if (a <= 0 || b <= 0 || c <= 0 || b < a) return 0;
    *k0 = a; *k1 = b; *step = c;
    return 1;
  }
  return 0;
}

static char *basename_noext(const char *path)
{
  const char *base = path;
  const char *slash = strrchr(path, '/');
  if (slash) base = slash + 1;
  size_t len = strlen(base);
  char *out = (char *) malloc(len + 1);
  if (out == NULL) return NULL;
  strcpy(out, base);
  char *dot = strrchr(out, '.');
  if (dot) *dot = '\0';
  return out;
}

static int parse_histex_counts(const char *path, double *unique, double *total)
{
  FILE *fp = fopen(path, "r");
  if (!fp) return 0;
  double u = 0.0, t = 0.0;
  char line[256];
  while (fgets(line, sizeof(line), fp))
  {
    char *p = line;
    while (isspace((unsigned char)*p)) p++;
    if (!isdigit((unsigned char)*p)) continue;
    long long mult = 0, count = 0;
    if (sscanf(p, "%lld %lld", &mult, &count) == 2)
    {
      u += (double) count;
      t += (double) mult * (double) count;
    }
  }
  fclose(fp);
  *unique = u;
  *total = t;
  return 1;
}

// ---- resource detection --------------------------------------------------

static long detect_cores(void)
{
  long n = sysconf(_SC_NPROCESSORS_ONLN);
  return (n > 0) ? n : 1;
}

// Available RAM in bytes. Prefer /proc/meminfo MemAvailable; fall back to
// total physical RAM from sysconf.
static double detect_avail_ram_bytes(void)
{
  FILE *fp = fopen("/proc/meminfo", "r");
  if (fp)
  {
    char key[64]; long val; char unit[16];
    while (fscanf(fp, "%63s %ld %15s", key, &val, unit) == 3)
    {
      if (strcmp(key, "MemAvailable:") == 0)
      {
        fclose(fp);
        return (double) val * 1024.0;   // kB -> bytes
      }
    }
    fclose(fp);
  }
  long pages = sysconf(_SC_PHYS_PAGES);
  long psize = sysconf(_SC_PAGE_SIZE);
  if (pages > 0 && psize > 0) return (double) pages * (double) psize;
  return 8.0e9;  // last-resort guess: 8 GB
}

static double file_size_bytes(const char *path)
{
  struct stat st;
  if (stat(path, &st) == 0) return (double) st.st_size;
  return 0.0;
}

// ---- run one k: FastK + Histex, write a CSV row to rowfile -----------------

static void run_one_k(int k, const char *input, const char *prefix,
                      const char *h_range, int threads, int keep,
                      const char *rowfile)
{
  char name[600];       // unique per-k FastK output name so jobs don't collide
  snprintf(name, sizeof(name), "%s.k%d", prefix, k);

  char fastk_log[700];
  snprintf(fastk_log, sizeof(fastk_log), "%s.fastk.log", name);

  char cmd[2048];
  snprintf(cmd, sizeof(cmd), "FastK -k%d -T%d \"%s\" -N\"%s\" > \"%s\" 2>&1",
           k, threads, input, name, fastk_log);
  int rc = system(cmd);

  char hist_path[700];
  snprintf(hist_path, sizeof(hist_path), "%s.hist", name);
  char histex_path[700];
  snprintf(histex_path, sizeof(histex_path), "%s.histex.txt", name);

  FILE *rf = fopen(rowfile, "w");
  if (!rf) return;

  if (rc != 0 || access(hist_path, F_OK) != 0)
  {
    fprintf(rf, "%d,ERROR_FASTK,ERROR_FASTK,ERROR_FASTK\n", k);
    fclose(rf);
    return;   // keep fastk_log for debugging
  }

  snprintf(cmd, sizeof(cmd), "Histex -h%s -A \"%s\" > \"%s\"",
           h_range, hist_path, histex_path);
  int rc2 = system(cmd);

  double unique = 0.0, total = 0.0;
  if (rc2 != 0 || !parse_histex_counts(histex_path, &unique, &total))
    fprintf(rf, "%d,ERROR_HISTEX,ERROR_HISTEX,ERROR_HISTEX\n", k);
  else
  {
    double frac = (total > 0.0) ? (unique / total) : 0.0;
    fprintf(rf, "%d,%.0f,%.0f,%.6f\n", k, unique, total, frac);
  }
  fclose(rf);

  if (!keep)
  {
    // Fastrm removes every FastK product for this name (.hist/.ktab/hidden parts)
    char rmcmd[2048];
    snprintf(rmcmd, sizeof(rmcmd),
             "Fastrm -f \"%s\" >/dev/null 2>&1; rm -f \"%s\" \"%s\"",
             name, histex_path, fastk_log);
    system(rmcmd);
  }
}

int main(int argc, char *argv[])
{
  const char *input = NULL, *k_range = NULL, *h_range = NULL;
  const char *prefix = NULL, *outcsv = NULL;
  int thread_budget = 0;   // 0 => auto
  int jobs_opt = 0;        // 0 => auto
  double mem_budget_gb = 0.0; // 0 => auto
  int keep = 0;

  for (int i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) input = argv[++i];
    else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) k_range = argv[++i];
    else if (strcmp(argv[i], "-h") == 0 && i + 1 < argc) h_range = argv[++i];
    else if (strcmp(argv[i], "-T") == 0 && i + 1 < argc) thread_budget = atoi(argv[++i]);
    else if (strcmp(argv[i], "-P") == 0 && i + 1 < argc) jobs_opt = atoi(argv[++i]);
    else if (strcmp(argv[i], "-M") == 0 && i + 1 < argc) mem_budget_gb = atof(argv[++i]);
    else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) prefix = argv[++i];
    else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) outcsv = argv[++i];
    else if (strcmp(argv[i], "-K") == 0) keep = 1;
    else { usage(argv[0]); return 1; }
  }

  if (input == NULL || k_range == NULL || h_range == NULL) { usage(argv[0]); return 1; }

  int k0=0, k1=0, step=0;
  if (!parse_k_range(k_range, &k0, &k1, &step))
  {
    fprintf(stderr, "ERROR: invalid -k range: %s\n", k_range);
    return 1;
  }

  char *auto_prefix = NULL;
  if (prefix == NULL)
  {
    auto_prefix = basename_noext(input);
    if (auto_prefix == NULL) { fprintf(stderr, "ERROR: cannot allocate prefix\n"); return 1; }
    prefix = auto_prefix;
  }

  char *auto_out = NULL;
  if (outcsv == NULL)
  {
    size_t n = strlen(prefix) + 12;
    auto_out = (char *) malloc(n);
    if (auto_out == NULL) { fprintf(stderr, "ERROR: cannot allocate outcsv\n"); free(auto_prefix); return 1; }
    snprintf(auto_out, n, "%s.kplex.csv", prefix);
    outcsv = auto_out;
  }

  // ---- build list of k values ----
  int nk = (k1 - k0) / step + 1;
  if (nk < 1) nk = 1;
  int *ks = (int *) malloc((size_t) nk * sizeof(int));
  if (ks == NULL) { fprintf(stderr, "ERROR: cannot allocate k list\n"); free(auto_prefix); free(auto_out); return 1; }
  int idx = 0;
  for (int k = k0; k <= k1 && idx < nk; k += step) ks[idx++] = k;
  nk = idx;

  // ---- resource-aware parallel plan ----
  long cores = detect_cores();
  double avail_ram = detect_avail_ram_bytes();
  double avail_gb  = (mem_budget_gb > 0.0) ? mem_budget_gb : avail_ram / 1e9;

  int budget = (thread_budget > 0) ? thread_budget : (int) cores;
  if (budget < 1) budget = 1;

  // Per-job memory estimate for FastK: a large fixed working-set component plus
  // a term that scales with genome size. Calibrated from measured single-FastK
  // peak RSS: ~5 GB for a 140 Mb genome and ~11 GB for a 3.1 Gb genome
  // (=> ~5.0 + 2.0*genome_gb), with a small safety margin.
  double genome_gb = file_size_bytes(input) / 1e9;
  double est_job_gb = 5.5 + 2.2 * genome_gb;
  if (est_job_gb < 1.0) est_job_gb = 1.0;

  int jobs_by_ram = (int) ((avail_gb * 0.90) / est_job_gb);
  if (jobs_by_ram < 1) jobs_by_ram = 1;

  int jobs;
  if (jobs_opt > 0) jobs = jobs_opt;                 // user override
  else {
    jobs = budget;                                   // start from thread budget
    if (jobs > jobs_by_ram) jobs = jobs_by_ram;      // cap by RAM
  }
  if (jobs > nk) jobs = nk;                          // no more jobs than k values
  if (jobs < 1) jobs = 1;

  int tpj = budget / jobs;                            // threads per FastK job
  if (tpj < 1) tpj = 1;

  fprintf(stderr,
          "Kplex: %ld cores, %.1f GB available; genome %.2f GB; %d k values\n"
          "       plan: %d parallel job(s) x %d thread(s)  (est %.1f GB/job, %.1f GB total)\n",
          cores, avail_gb, genome_gb, nk, jobs, tpj, est_job_gb, est_job_gb * jobs);

  // ---- per-k temporary row files ----
  char **rowfiles = (char **) malloc((size_t) nk * sizeof(char *));
  for (int i = 0; i < nk; i++)
  {
    size_t n = strlen(prefix) + 32;
    rowfiles[i] = (char *) malloc(n);
    snprintf(rowfiles[i], n, "%s.k%d.row", prefix, ks[i]);
  }

  // ---- process pool: keep up to `jobs` children running ----
  int next = 0, running = 0;
  while (next < nk || running > 0)
  {
    while (running < jobs && next < nk)
    {
      int i = next++;
      pid_t pid = fork();
      if (pid < 0)
      {
        // fork failed: run inline so we never drop a k value
        run_one_k(ks[i], input, prefix, h_range, tpj, keep, rowfiles[i]);
      }
      else if (pid == 0)
      {
        run_one_k(ks[i], input, prefix, h_range, tpj, keep, rowfiles[i]);
        _exit(0);
      }
      else running++;
    }
    if (running > 0) { int st; wait(&st); running--; }
  }

  // ---- collect rows in k order → CSV ----
  FILE *out = fopen(outcsv, "w");
  if (!out)
  {
    fprintf(stderr, "ERROR: cannot write %s\n", outcsv);
    free(auto_prefix); free(auto_out); free(ks);
    return 1;
  }
  fprintf(out, "k,unique_kmers,total_kmers,fraction_unique\n");
  for (int i = 0; i < nk; i++)
  {
    FILE *rf = fopen(rowfiles[i], "r");
    if (rf)
    {
      char line[256];
      if (fgets(line, sizeof(line), rf)) fputs(line, out);
      else fprintf(out, "%d,ERROR_FASTK,ERROR_FASTK,ERROR_FASTK\n", ks[i]);
      fclose(rf);
      remove(rowfiles[i]);
    }
    else
      fprintf(out, "%d,ERROR_FASTK,ERROR_FASTK,ERROR_FASTK\n", ks[i]);
    free(rowfiles[i]);
  }
  fclose(out);

  free(rowfiles);
  free(ks);
  free(auto_prefix);
  free(auto_out);
  return 0;
}

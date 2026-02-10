// Kplex: convenience wrapper for FastK + Histex to produce k-plexity curves
// Usage:
//   Kplex -i genome.fa -k 5:151:1 -h 1:1000 [-T 4] [-o prefix] [-c out.csv] [-K]
//
// Produces CSV with columns: k,unique_kmers,total_kmers,fraction_unique

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

static void usage(const char *prog)
{
  fprintf(stderr,
          "Usage: %s -i genome.fa -k start:end:step -h h1:h2 [-T threads] [-o prefix] [-c out.csv] [-K]\n"
          "  -i  input FASTA\n"
          "  -k  k range, e.g. 5:151:1\n"
          "  -h  Histex range, e.g. 1:1000 (passed as -h1:1000)\n"
          "  -T  FastK threads (default 4)\n"
          "  -o  prefix for FastK output (default: basename of fasta)\n"
          "  -c  output CSV (default: <prefix>.kplex.csv)\n"
          "  -K  keep intermediate .hist and histex output files\n", prog);
}

static int parse_k_range(const char *s, int *k0, int *k1, int *step)
{
  // accept start:end:step or start-end:step
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
    // Skip non-numeric lines
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

int main(int argc, char *argv[])
{
  const char *input = NULL;
  const char *k_range = NULL;
  const char *h_range = NULL;
  const char *prefix = NULL;
  const char *outcsv = NULL;
  int threads = 4;
  int keep = 0;

  for (int i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "-i") == 0 && i + 1 < argc)
      input = argv[++i];
    else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc)
      k_range = argv[++i];
    else if (strcmp(argv[i], "-h") == 0 && i + 1 < argc)
      h_range = argv[++i];
    else if (strcmp(argv[i], "-T") == 0 && i + 1 < argc)
      threads = atoi(argv[++i]);
    else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc)
      prefix = argv[++i];
    else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc)
      outcsv = argv[++i];
    else if (strcmp(argv[i], "-K") == 0)
      keep = 1;
    else
    {
      usage(argv[0]);
      return 1;
    }
  }

  if (input == NULL || k_range == NULL || h_range == NULL)
  {
    usage(argv[0]);
    return 1;
  }

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
    if (auto_prefix == NULL)
    {
      fprintf(stderr, "ERROR: cannot allocate prefix\n");
      return 1;
    }
    prefix = auto_prefix;
  }

  char *auto_out = NULL;
  if (outcsv == NULL)
  {
    size_t n = strlen(prefix) + 12;
    auto_out = (char *) malloc(n);
    if (auto_out == NULL)
    {
      fprintf(stderr, "ERROR: cannot allocate outcsv\n");
      free(auto_prefix);
      return 1;
    }
    snprintf(auto_out, n, "%s.kplex.csv", prefix);
    outcsv = auto_out;
  }

  FILE *out = fopen(outcsv, "w");
  if (!out)
  {
    fprintf(stderr, "ERROR: cannot write %s\n", outcsv);
    free(auto_prefix);
    free(auto_out);
    return 1;
  }
  fprintf(out, "k,unique_kmers,total_kmers,fraction_unique\n");

  for (int k = k0; k <= k1; k += step)
  {
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "FastK -k%d -T%d \"%s\" -N\"%s\" > /dev/null 2>&1",
             k, threads, input, prefix);
    system(cmd);

    char hist_path[512];
    snprintf(hist_path, sizeof(hist_path), "%s.hist", prefix);

    char histex_path[512];
    snprintf(histex_path, sizeof(histex_path), "%s.k%d.histex.txt", prefix, k);

    snprintf(cmd, sizeof(cmd), "Histex -h%s -A \"%s\" > \"%s\"",
             h_range, hist_path, histex_path);
    system(cmd);

    double unique = 0.0, total = 0.0;
    if (!parse_histex_counts(histex_path, &unique, &total))
    {
      fprintf(out, "%d,ERROR,ERROR,ERROR\n", k);
    }
    else
    {
      double frac = (total > 0.0) ? (unique / total) : 0.0;
      fprintf(out, "%d,%.0f,%.0f,%.6f\n", k, unique, total, frac);
    }

    if (!keep)
    {
      char rmcmd[1024];
      snprintf(rmcmd, sizeof(rmcmd), "rm -f \"%s\" \"%s\"", hist_path, histex_path);
      system(rmcmd);
    }
  }

  fclose(out);
  free(auto_prefix);
  free(auto_out);
  return 0;
}


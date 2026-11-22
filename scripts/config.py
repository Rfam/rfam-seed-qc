"""
Configuration parameters for Stockholm validation and fixing.

Modify these values to adjust the behavior of the validation script.
"""

import os

# =============================================================================
# BLAST Search Parameters
# =============================================================================

# Minimum identity percentage to accept a BLAST hit (0-100)
BLAST_MIN_IDENTITY = 95

# Minimum query coverage percentage to accept a BLAST hit (0-100)
BLAST_MIN_COVERAGE = 90

# Maximum e-value to accept a BLAST hit
BLAST_MAX_EVALUE = 1e-10

# Number of BLAST hits to retrieve (top N hits). BLAST also uses this while
# searching (-max_target_seqs), and small values can miss the best hit, so
# keep it generous: nt holds many near-identical copies of a sequence.
BLAST_HITLIST_SIZE = 250

# BLAST task: standard blastn (NOT megablast), with its scoring parameters
# stated explicitly so results do not change silently if NCBI changes its
# server defaults (the NCBI website and BLAST+ default to megablast).
# blastn is slower than megablast but more sensitive: it also finds
# divergent hits and short queries that megablast's 28-nt seed can miss.
BLAST_WORD_SIZE = 11
BLAST_MATCH_REWARD = 2
BLAST_MISMATCH_PENALTY = -3
BLAST_GAP_COSTS = "5 2"   # gap open, gap extend (affine)

# Local BLAST database. Do not edit here: it is set at run time from
# --local --blast-db BLAST_LOCAL_DB (default "nt"). BLAST searches and
# sequence lookups then use this database (NCBI only for accessions it
# lacks). None means BLAST uses NCBI's web API and sequences come from NCBI.
BLAST_LOCAL_DB = None

# CPU threads for local blastn (override with --threads). Default: the CPUs
# this process may run on (on Linux this follows CPU pinning, e.g. SLURM's),
# capped by SLURM_CPUS_PER_TASK when set, for clusters that do not pin CPUs.
# CPU-time quotas (e.g. docker --cpus) are not detected: use --threads there.
def _default_threads():
    n = (len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity')
         else (os.cpu_count() or 1))
    slurm = os.environ.get('SLURM_CPUS_PER_TASK', '')
    if slurm.isdigit() and int(slurm) > 0:
        n = min(n, int(slurm))
    return n


BLAST_NUM_THREADS = _default_threads()


# =============================================================================
# Rfam CM Filtering Parameters (cmscan)
# =============================================================================

# Maximum e-value to consider a cmscan hit significant (sequences with
# E-value below this threshold are flagged as warnings)
CMSCAN_MAX_EVALUE = 1e-3


# =============================================================================
# Base-pair support (SS_cons checks; warnings only, nothing is removed)
# =============================================================================

# Warn about a base pair if fewer than this fraction of the sequences with
# both bases can form it (Watson-Crick or G-U)
STRUCTURE_MIN_PAIR_SUPPORT = 0.5

# Warn about a base pair if more than this fraction of sequences have a gap
# at one or both of its positions
STRUCTURE_MAX_PAIR_GAPS = 0.5

# Warn about a sequence if it cannot form more than this fraction of the
# base pairs it has both bases for (often a misaligned sequence)
STRUCTURE_MAX_SEQ_MISPAIRED = 0.3


# =============================================================================
# NCBI Request Parameters
# =============================================================================

# Delay between NCBI requests in seconds (to avoid rate limiting)
NCBI_REQUEST_DELAY = 0.5

# Timeout in seconds for a single NCBI FASTA download
NCBI_TIMEOUT = 60

# Identify this tool to NCBI (BLAST and E-utilities), as NCBI's usage
# guidelines request, so NCBI can contact the user before blocking an IP.
# Set NCBI_EMAIL to your own address; leave as None to send no email.
NCBI_TOOL = 'rfam-seed-qc'
NCBI_EMAIL = None

# Optional NCBI API key, read from the NCBI_API_KEY environment variable (never
# store it here: the repository is public). Get one from your NCBI account
# settings. With a key, NCBI E-utilities allow 10 requests/s per key instead of
# 3 per IP, which matters on clusters where many users share one IP. Not used
# for web BLAST (NCBI's BLAST service does not take E-utilities keys).
NCBI_API_KEY = os.environ.get('NCBI_API_KEY', '').strip() or None
NCBI_REQUEST_DELAY_WITH_KEY = 0.1

# Seconds to wait before retrying a web BLAST request that NCBI refused
# (429/5xx); multiplied by the attempt number (30 s, then 60 s).
NCBI_BLAST_RETRY_WAIT = 30


# =============================================================================
# RNAcentral Request Parameters (for URS identifiers)
# =============================================================================

# Delay between RNAcentral /api/v1/rna/ requests in seconds. The enforced limit
# is 30 requests/minute per IP; RNAcentral recommends one request every 2 s.
RNACENTRAL_REQUEST_DELAY = 2.0

# Retries on HTTP 429 / 5xx / network errors (Retry-After is honoured)
RNACENTRAL_MAX_RETRIES = 3

# Per-request timeout in seconds
RNACENTRAL_TIMEOUT = 30

# The species-specific endpoint (URS_taxid) used to verify the taxid can be
# very slow or return HTTP 500 for URS found in many species. Allow a longer
# timeout but few retries; if it still fails the taxid is reported as
# unverified and the (sequence-validated) entry is kept.
RNACENTRAL_TAXID_TIMEOUT = 60
RNACENTRAL_TAXID_MAX_RETRIES = 1

# Descriptive User-Agent, as requested by RNAcentral's access guidelines
RNACENTRAL_USER_AGENT = 'rfam-seed-qc (https://github.com/Rfam/rfam-seed-qc)'

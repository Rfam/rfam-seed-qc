"""
Configuration parameters for Stockholm validation and fixing.

Modify these values to adjust the behavior of the validation script.
"""

# =============================================================================
# BLAST Search Parameters
# =============================================================================

# Minimum identity percentage to accept a BLAST hit (0-100)
BLAST_MIN_IDENTITY = 95

# Minimum query coverage percentage to accept a BLAST hit (0-100)
BLAST_MIN_COVERAGE = 90

# Maximum e-value to accept a BLAST hit
BLAST_MAX_EVALUE = 1e-10

# Number of BLAST hits to retrieve (top N hits)
BLAST_HITLIST_SIZE = 5

# BLAST task: standard blastn (NOT megablast), with its scoring parameters
# stated explicitly so results do not change silently if NCBI changes its
# server defaults (the NCBI website and BLAST+ default to megablast).
# blastn is slower than megablast but more sensitive: it also finds
# divergent hits and short queries that megablast's 28-nt seed can miss.
BLAST_WORD_SIZE = 11
BLAST_MATCH_REWARD = 2
BLAST_MISMATCH_PENALTY = -3
BLAST_GAP_COSTS = "5 2"   # gap open, gap extend (affine)


# =============================================================================
# Rfam CM Filtering Parameters (cmscan)
# =============================================================================

# Maximum e-value to consider a cmscan hit significant (sequences with
# E-value below this threshold are flagged as warnings)
CMSCAN_MAX_EVALUE = 1e-3


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

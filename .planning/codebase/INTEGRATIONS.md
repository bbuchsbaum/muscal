# External Integrations

**Analysis Date:** 2026-01-22

## APIs & External Services

**Not Detected:**
- No external HTTP APIs are consumed
- No authentication tokens or API keys required
- No remote service integrations

## Data Storage

**Databases:**
- Not applicable - This is a statistical analysis package
- Data is passed as R data structures (matrices, lists, multiblock objects)
- No persistent database storage

**File Storage:**
- Local filesystem only - data read/written through standard R file I/O
- No cloud storage integration (S3, GCS, etc.)

**Caching:**
- None - Computations are performed in-memory
- Results stored in R objects returned to user

## Authentication & Identity

**Auth Provider:**
- Not applicable - No user authentication system
- CRAN package distribution uses GPG signing for integrity (standard R package practice)

## Monitoring & Observability

**Error Tracking:**
- None - Uses standard R error handling via `tryCatch()` blocks in `R/mfa.R:44-47`
- Validation errors thrown via `chk::` functions

**Logs:**
- Console output only via `message()`, `cat()`, and crayon-colored output
- CLI alerts via `cli::cli_alert_*()` functions in output messages
- Verbose mode available in bootstrap functions (see `verbose=TRUE` parameter in `R/bada.R:82`)

**Debugging:**
- No structured logging framework
- Messages displayed via R console

## CI/CD & Deployment

**Hosting:**
- CRAN (The Comprehensive R Archive Network) - standard R package repository
- GitHub repository at `https://github.com/bbuchsbaum/muscal`
- Installable via `install.packages("muscal")` once on CRAN

**CI Pipeline:**
- Not detected - No CI configuration files found (.github/workflows, .travis.yml, etc.)
- Manual testing via testthat test suite located in `tests/testthat/`

**Build Process:**
- R CMD CHECK - Standard R package validation
- Roxygen2 documentation generation (handles documentation, namespace, exports)

## Environment Configuration

**Required env vars:**
- None - No environment variables required for package functionality

**Secrets location:**
- Not applicable - No secrets or credentials needed

## Webhooks & Callbacks

**Incoming:**
- None - Not applicable for statistical analysis package

**Outgoing:**
- None - No outbound callbacks or webhook notifications

## Inter-Package Dependencies

**Critical Package Interactions:**
- Imports from `multivarious` package: projectors, bootstrap utilities, matrix operations
- Imports from `multidesign` package: multiblock object creation and design handling
- Imports from `genpca` package: underlying generalized PCA algorithm

**Data Flow with External Packages:**
1. User provides data blocks to functions like `mfa()` in `R/mfa.R:127`
2. Data preprocessed using preprocessing utilities from `multivarious`
3. Passed to `genpca::genpca()` at `R/mfa.R:181` for computation
4. Results wrapped in `multivarious::multiblock_biprojector()` at `R/mfa.R:199`

---

*Integration audit: 2026-01-22*

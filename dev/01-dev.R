###########################
###########################
#### 01-dev.R

#### Aims
# 1) Record project development

#### Prerequisites
# 1) NA


###########################
###########################
#### Project set up

#### Use git
# usethis::use_git()
# usethis::use_github()
usethis::git_vaccinate()

#### Use dependency management
renv::init()

#### Install package(s)
# Project management 
# (NB: requireNamespace() is necessary for renv to register the packages
if (!requireNamespace("proj.verse", quietly = TRUE)) {
  renv::install("edwardlavender/proj.verse", prompt = FALSE)
}
# Modelling
if (!requireNamespace("patter", quietly = TRUE)) {
  renv::install("edwardlavender/patter", prompt = FALSE)
}
if (!requireNamespace("patter.workflows", quietly = TRUE)) {
  renv::install("edwardlavender/patter.workflows", prompt = FALSE)
}
if (!requireNamespace("Rfast", quietly = TRUE)) {
  renv::install("Rfast", prompt = FALSE)
}
# Utilities
if (!requireNamespace("utils.add", quietly = TRUE)) {
  renv::install("edwardlavender/utils.add", prompt = FALSE)
}
if (!requireNamespace("usethis", quietly = TRUE)) {
  renv::install("usethis", prompt = FALSE)
}
if (!requireNamespace("lobstr", quietly = TRUE)) {
  renv::install("lobstr", prompt = FALSE)
}
# commonmark/(r)markdown packages (for README documentation)
if (!requireNamespace("commonmark", quietly = TRUE)) {
  renv::install("commonmark", prompt = FALSE)
}
if (!requireNamespace("markdown", quietly = TRUE)) {
  renv::install("markdown", prompt = FALSE)
}
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  renv::install("rmarkdown", prompt = FALSE)
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  renv::install("yaml", prompt = FALSE)
}

#### Use proj.templates templates
library(proj.verse)

# Set up template project structure
use_template_proj()

# Update .gitignore
use_template_gitignore()

# Add a README and associated files
usethis::use_code_of_conduct("edward.lavender@eawag.ch")
use_template_readme(title = "Estimating residency in acoustic telemetry systems",
                        author = "Edward Lavender",
                        email = "edward.lavender@eawag.ch")

# Add template scripts
if (!requireNamespace("pacman", quietly = TRUE))
  renv::install("pacman", prompt = FALSE)
use_template_script(here_r("insert_script_name_1.R"))
use_template_script(here_r("insert_script_name_2.R"))


###########################
###########################
#### Project maintenance

#### Enforce consistent syntax
# usethis::use_tidy_style()
# Check code is syntactically valid
lapply(list.files(here_r(), full.names = TRUE, pattern = ".R"), parse)

#### Check project spelling
spelling::spell_check_files("README.Rmd", lang = "en-GB")
spelling::spell_check_files(
  list.files(here_r(), full.names = TRUE, pattern = ".R"),
  lang = "en-GB")

#### List project dependencies
# Define helper functions
pkg_version <- function(pkg) {
  packageDescription(pkg)$Version
}
pkg_github_user <- function(pkg) {
  usr <- packageDescription(pkg)$GithubUsername
  ifelse(is.null(usr), "", paste0(usr, "/"))
}
pkg_github_commit <- function(pkg) {
  com <- packageDescription(pkg)$RemoteSha
  ifelse(is.null(com), "", com)
}
# Get dependencies
pkg <- renv::dependencies()
unique(pkg$Require); unique(pkg$Version); unique(pkg$Dev)
pkg <- data.frame(package = sort(unique(pkg$Package)),
                  version = NA,
                  github_user = NA,
                  github_commit = NA)
# Get versions and github source (if necessary)
for (i in seq_len(nrow(pkg))) {
  pkg$version[i] <- pkg_version(pkg$package[i])
  pkg$github_user[i]   <- pkg_github_user(pkg$package[i])
  pkg$github_commit[i] <- pkg_github_commit(pkg$package[i])
}
pkg$version[pkg$github_user != ""] <- pkg$github_commit[pkg$github_user != ""]
# Define install code
pkg$install <- paste0("renv::install('",
                      pkg$github_user,
                      pkg$package,
                      "@", pkg$version,
                      "', prompt = FALSE)")
pkg <- pkg[, c("package", "install")]
# Save dataframe
# View(pkg)
saveRDS(pkg, here_data("inst", "dependencies.rds"))

#### Update renv
## Take snapshot
renv::snapshot()
## Clean snapshot
# Note that this may attempt to drop 'suggested packages'
# ... that are required by (some) functions from other packages
# ... but which are not used directly. To guard against this,
# ... make an arbitrary call to the required packages
# ... where they are needed.
renv::clean()

#### Save sessionInfo
saveRDS(sessionInfo(), here_data("inst", "session-info.rds"))

#### Save the project directory 'tree'
# ... This enables the project directory tree to be rebuilt on another machine
# ... This function should be re-run when the directory tree is updated
use_template_tree(save = here_data("inst", "tree.rds"))


#### End of code.
###########################
###########################

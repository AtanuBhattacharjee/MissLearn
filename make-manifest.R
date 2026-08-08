# =============================================================================
# make-manifest.R
#
# Run this once from the MissLearn folder, before the first push to GitHub,
# and again any time you add or remove a package in app.R.
#
#   Rscript make-manifest.R
#
# or, in RStudio, open this file and click Source.
# =============================================================================

needed <- c("shiny", "mice", "rsconnect")

missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

# --- check the app at least parses and loads its packages -------------------
message("\nParsing app.R ...")
parsed <- tryCatch(parse("app.R"), error = function(e) {
  stop("app.R does not parse:\n", conditionMessage(e), call. = FALSE)
})
message("  ", length(parsed), " top-level expressions, no syntax errors.")

suppressPackageStartupMessages({
  library(shiny)
  library(mice)
})
message("  shiny ", as.character(utils::packageVersion("shiny")),
        " and mice ", as.character(utils::packageVersion("mice")), " loaded.")

# --- write the manifest -----------------------------------------------------
message("\nWriting manifest.json ...")
rsconnect::writeManifest(appDir = ".", appPrimaryDoc = "app.R")

mf <- jsonlite::fromJSON("manifest.json")
message("  R version recorded: ", mf$platform)
message("  Packages locked:    ", length(mf$packages))

message("\nDone. Commit manifest.json alongside app.R, then push.")
message("If the Connect Cloud build fails on a missing package, re-run this ",
        "script and push again.")

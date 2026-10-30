# Run archived R scripts from the repository root.
plain_path <- function(root, relative = "") {
  repo <- normalizePath(Sys.getenv("PLAIN_REPO_ROOT", "."), mustWork = TRUE)
  if (root == "PYTHON") return(Sys.getenv("PLAIN_PYTHON", Sys.which("python3")))
  defaults <- c(REPO = repo,
    PROJECT = file.path(repo, "private_workspace/project"),
    HEATMAP = file.path(repo, "private_workspace/heatmap"),
    LEGACY_WORK = file.path(repo, "private_workspace/legacy_work"),
    WORKSPACE = file.path(repo, "private_workspace/workspace"),
    DESKTOP_OUTPUTS = file.path(repo, "private_workspace/outputs"),
    R_LIBRARY = .libPaths()[1])
  if (!root %in% names(defaults)) stop("Unknown path root: ", root)
  base <- path.expand(Sys.getenv(paste0("PLAIN_", root, "_ROOT"), defaults[[root]]))
  if (!nzchar(relative)) return(base)
  paste0(sub("/$", "", base), "/", relative)
}

# path wrapper to access file from specific reference/commit
blob <- function(file_path, ref = 'HEAD') {
    # check file_path
    if (file.info(file_path)$isdir) {
        stop('argument "file_path" must be a file.')
    }
    # get parent directory
    dir_path <- normalizePath(dirname(file_path))
    # get file name
    file_name <- basename(file_path)
    # get top level
    git_top <- system(paste0('cd ', dir_path, ' && git rev-parse --show-toplevel'), intern = TRUE)
    # update file path
    git_rel_path <- sub(paste0('^', git_top, '/?'), '', dir_path)
    git_file_path <- file.path(git_rel_path, file_name)
    # use git show to copy file
    # TODO: capture path not valid (file didn't exist)
	gitshow <- paste0("git -C ", git_top, " show ", ref, ":", git_file_path)
    path_out <- tempfile('blob')
    system(paste0("git -C ", git_top, " show ", ref, ":", git_file_path, ' > ', path_out), intern = TRUE)
    # return path to temporary file
    path_out
}


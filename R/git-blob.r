# path wrapper to access file from specific reference/commit
blob <- function(file_path, ref = 'HEAD') {
    # check file_path
    if (isTRUE(isdir <- file.info(file_path)$isdir)) {
        stop('argument "file_path" must be a file.')
    }
    # get parent directory
    dir_path <- normalizePath(dirname(file_path))
    # get top level
    git_top <- system(paste0('cd ', dir_path, ' && git rev-parse --show-toplevel'), intern = TRUE)
    if (length(attr(git_top, 'status')) > 0) {
        stop('file not inside a git repository')
    }
    # get file name
    file_name <- basename(file_path)
    # update file path
    git_rel_path <- sub(paste0('^', git_top, '/?'), '', dir_path)
    git_file_path <- file.path(git_rel_path, file_name)
    # fallthrough if file is actually present as is
    # get current blob hash
    current_hash <- suppressWarnings(try(system(paste0('git -C ', git_top, ' hash-object ', git_file_path), intern = TRUE), silent = TRUE))
    # reference blob hash
    ref_hash <- suppressWarnings(try(system(paste0('git -C ', git_top, ' rev-parse ', ref, ':', git_file_path), intern = TRUE), silent = TRUE))
    if (length(attr(ref_hash, 'status')) > 0) {
        stop('file does not exist in specified reference under specified path!\nProvide a valid path under reference "', ref, '"')
    } else if (is.na(isdir) || current_hash != ref_hash) {
        # use git show to copy file
        gitshow <- paste0("git -C ", git_top, " show ", ref, ":", git_file_path)
        path_out <- tempfile('blob')
        system(paste0("git -C ", git_top, " show ", ref, ":", git_file_path, ' > ', path_out), intern = TRUE)
        # fix lfs pointers
        file_head <- suppressWarnings(readLines(path_out, n = 1))
        if (file_head == 'version https://git-lfs.github.com/spec/v1') {
            path_pointer <- path_out
            path_out <- tempfile('blob')
            system(paste0('git -C ', git_top, ' lfs smudge < ', path_pointer, ' > ', path_out), intern = TRUE)
        }
        # return temporary path
        return(path_out)
    }
    # file in worktree is ok (lfs should be checked out)
    return(file_path)
}


# issue with git due to large files => store only latest version locally,
# other versions are stored on a remote file system
# save file locally with correct name
# save with hash as name on "remote" drive
# add 2 functions to read (incl. check) and write
# add sha1 as attribute to object

# main function to write "hashed" files
write_hashed <- function(x, file_path, remote_path = getOption('remote_path'), 
    save_local = TRUE, gitignore = TRUE, write_fun = qd_save, ...) {
    # get local file name
    local_name <- basename(file_path)
    # TODO: add split to several files option
    # get local path
    local_path <- dirname(file_path)
    # check directories
    local_path <- check_path(local_path)
    if (is.null(remote_path)) {
        stop('argument "remote_path" is missing!')
    }
    remote_path <- check_path(remote_path)
    # get write function
    if (is.function(write_fun)) {
        fun_name <- deparse(substitute(write_fun))
    } else {
        fun_name <- write_fun
        write_fun <- get(write_fun, mode = 'function')
    }
    # be verbose
    fun_call <- sub('^list', fun_name, deparse(list(...)))
    cat('writing files using ', fun_call, '\n', sep = '')
    # save file on remote
    cat('writing remote file...\n')
    tmp_remote <- file.path(remote_path, local_name)
    write_fun(x, tmp_remote, ...)
    # get file hash
    file_hash <- as.character(openssl::sha1(file(tmp_remote)))
    # move file
    remote_file <- file.path(remote_path, 
        paste(local_name, file_hash, sep = '-'))
    cat('-> remote file name:', basename(remote_file), '\n')
    if (file.exists(remote_file)) {
        file.remove(tmp_remote)
    } else {
        rename_ok <- file.rename(tmp_remote, remote_file)
        if (!rename_ok) {
            stop('renaming file from ', tmp_remote, ' to ', 
                remote_file, ' failed!')
        }
    }
    # save hash locally
    cat('updating local hash file...\n')
    # add hash directory to local
    dir.create(hash_path <- file.path(local_path, '.hash'), 
        showWarnings = FALSE)
    # add .hashed to filename
    # TODO: add remote path too?
    local_hash <- paste0(local_name, '.hash')
    writeLines(c(file_hash, fun_name, remote_file), file.path(hash_path, local_hash))
    # save file locally
    if (save_local) {
        cat('copying remote file to local path...\n')
        suppressWarnings(file.copy(remote_file, file_path, recursive = TRUE))
        if (gitignore) {
            gitignore_file(file_path)
        }
    }
    invisible(TRUE)
}

# main function to read "hashed" files
read_hashed <- function(file_path, remote_path = getOption('remote_path'), 
    update_local = TRUE, gitignore = FALSE, read_fun = qd_read, ...) {
    # get local path
    local_path <- dirname(file_path)
    # check local directory
    local_path <- check_path(local_path)
    # get file name
    local_name <- basename(file_path)
    # check if .hash directory exists
    if (dir.exists(hash_path <- file.path(local_path, '.hash'))) {
        cat('checking hash...\n')
        # get repo hash
        local_hash <- paste0(local_name, '.hash')
        hash_commit <- readLines(file.path(hash_path, local_hash))
        # get hash from local file
        local_hash <- suppressWarnings(
            try(as.character(openssl::sha1(file(file_path))), silent = TRUE)
        )
        if (inherits(local_hash, 'try-error')) {
            cat('-> local file does not exist\n')
            local_hash <- NULL
        } else if (hash_commit[1] == local_hash) {
            cat('-> local file is up-to-date\n')
        } else {
            cat('-> local file needs to be updated\n')
        }
        # check file status
        read_remote <- !update_local || is.null(local_hash) || hash_commit[1] != local_hash
        # get read function
        if (missing(read_fun) && grepl('save|write', hash_commit[2])) {
            read_fun <- get(sub('save|write', 'read', hash_commit[2]))
            fun_name <- sub('save|write', 'read', hash_commit[2])
        } else if (is.function(read_fun)) {
            fun_name <- deparse(substitute(read_fun))
        } else {
            fun_name <- read_fun
            read_fun <- get(read_fun, mode = 'function')
        }
    } else {
        if (!file.exists(file_path)) {
            stop('neither .hash directory nor local file exist!')
        } else if (update_local) {
            warning('.hash directory is missing. Cannot check local file status.')
        } else {
            stop('.hash directory is missing. Cannot read remote file.')
        }
        read_remote <- FALSE
    }
    # update local file
    if (read_remote) {
        if (is.null(remote_path)) {
            stop('argument "remote_path" is missing!')
        }
        remote_path <- check_path(remote_path)
        # be verbose
        fun_call <- sub('^list', fun_name, deparse(list(...)))
        cat('reading from remote calling ', fun_call, '...\n', sep = '')
        remote_file <- paste(local_name, hash_commit[1], sep = '-')
        out <- read_fun(file.path(remote_path, remote_file), ...)
        # save to local system
        if (update_local) {
            cat('updating local file from remote...\n')
            suppressWarnings(file.copy(file.path(remote_path, remote_file),
                file_path, recursive = TRUE))
            if (gitignore) {
                gitignore_file(file_path)
            }
        }
        # return remote
        out
    } else {
        # return local
        cat('reading local file...\n')
        read_fun(file_path, ...)
    }
}

check_path <- function(path) {
    path <- normalizePath(path, mustWork = FALSE)
    if (!dir.exists(path)) {
        stop("Can't access directory \"", path, '"')
    }
    path
}

gitignore_file <- function(file_path) {
    dname <- dirname(file_path)
    repo_root <- suppressWarnings(system(paste0(
            'cd ', dname, ' && git rev-parse --show-toplevel '
    ), intern = TRUE, ignore.stderr = TRUE))
    if (length(repo_root)) {
        # file path relative to repo
        repo_file <- sub(paste0(repo_root, '/'), '', normalizePath(file_path), fixed = TRUE)
        # get git status
        git_status <- system(paste0(
                'git -C ', repo_root, ' status --porcelain'
            ), intern = TRUE)
        if (length(git_status)) {
            # check if file is tracked
            if (repo_file %in% sub('.{3}', '', git_status)) {
                cat('adding "', repo_file, '" to .gitignore\n', sep = '')
                gitignore_path <- file.path(repo_root, '.gitignore')
                con <- file(gitignore_path, open = 'at')
                writeLines(
                    c(
                    '# ignore large file created by R function write_hashed()',
                    repo_file
                    )
                    , con = con
                )
                close(con)
            }
        }
    } else {
        warning(normalizePath(file_path), ' not in a git repository!')
    }
    invisible(NULL)
}

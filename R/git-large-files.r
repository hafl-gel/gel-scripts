# issue with git due to large files => store only latest version locally,
# other versions are stored on a remote file system
# save file locally with correct name
# save with hash as name on "remote" drive
# add 2 functions to read (incl. check) and write
# add sha1 as attribute to object

# main function to write "hashed" files
write_hashed <- function(x, file_path, remote_path = getOption('remote_path'), 
    save_local = TRUE, ...) {
    # get file name
    file_name <- basename(file_path)
    # add .hashed to filename
    file_name <- sub('([.]hashed)?$', '.hashed', file_name)
    # get local path
    local_path <- dirname(file_path)
    # check directories
    local_path <- check_path(local_path)
    if (is.null(remote_path)) {
        stop('argument "remote_path" is missing!')
    }
    remote_path <- check_path(remote_path)
    # add hash directory to local
    dir.create(file.path(local_path, '.hash'), showWarnings = FALSE)
    # add file hash
    add_hash(x)
    # save file on remote
    cat('writing remote file...\n')
    qs2::qd_save(x, file.path(remote_path, get_hash(x)), preset = 'fast')
    # save hash locally
    cat('updating hash file...\n')
    writeLines(get_hash(x), file.path(local_path, '.hash', file_name))
    # save file locally
    if (save_local) {
        cat('writing local file...\n')
        write_local(x, file.path(local_path, file_name), ...)
    }
    invisible(TRUE)
}

# main function to read "hashed" files
read_hashed <- function(file_path, remote_path = getOption('remote_path'), 
    update_local = TRUE) {
    # add .hashed to filename
    file_path <- sub('([.]hashed)?$', '.hashed', file_path)
    # get local path
    local_path <- dirname(file_path)
    # check local directory
    local_path <- check_path(local_path)
    # get file name
    file_name <- basename(file_path)
    # check if .hash folder exists
    if (dir.exists(file.path(local_path, '.hash'))) {
        cat('checking hash...\n')
        # get repo hash
        hash_repo <- readLines(file.path(local_path, '.hash', file_name), n = 1L)
        # get hash from local file: NULL -> file is missing
        hash_file <- read_local(file_path, hash_only = TRUE)
        if (is.null(hash_file)) {
            cat('-> local file does not exist\n')
        } else if (hash_repo == hash_file) {
            cat('-> local file is up-to-date\n')
        } else {
            cat('-> local file needs to be updated\n')
        }
        # check file status
        read_remote <- !update_local || is.null(hash_file) || hash_repo != hash_file
    } else {
        if (update_local) {
            warning('.hash folder is missing. Cannot check local file status.')
        } else {
            stop('.hash folder is missing. Cannot read remote file.')
        }
        read_remote <- FALSE
    }
    # update local file
    if (read_remote) {
        if (is.null(remote_path)) {
            stop('argument "remote_path" is missing!')
        }
        remote_path <- check_path(remote_path)
        # read with qs2::qd_read from remote
        cat('reading from remote...\n')
        out <- alloc.col(qs2::qd_read(file.path(remote_path, hash_repo)))
        # save to local system
        if (update_local) {
            cat('updating local file from remote...\n')
            write_local(out, file_path)
        }
        # return remote
        out
    } else {
        # return local
        cat('reading local file...\n')
        read_local(file_path)
    }
}

## helper functions
# sha1 from list of first time, last time, .N, names
hash <- function(x) {
    x[, digest::sha1(
        list(
            # start
            Time[1],
            # end
            Time[.N],
            # .N
            .N,
            # names
            names(x)
        )
    )]
}

add_hash <- function(x) setattr(x, 'hash', hash(x))

get_hash <- function(x) attr(x, 'hash')

check_path <- function(path) {
    path <- normalizePath(path, mustWork = FALSE)
    if (!dir.exists(path)) {
        stop("Can't access directory \"", path, '"')
    }
    path
}

write_local <- function(dat, path, ...) {
    dat_ser <- qs2::qd_serialize(dat, ...)
    con <- file(path, open = 'wb')
    on.exit(close(con))
    # write hash
    writeBin(get_hash(dat), con)
    # write number of bytes
    writeBin(length(dat_ser), con)
    # write serialized dat
    writeBin(dat_ser, con)
}

read_local <- function(path, hash_only = FALSE) {
    # check if file exists
    if (!file.exists(path)) return(NULL)
    # read local file
    con <- file(path, open = 'rb')
    on.exit(close(con))
    # read hash
    hash_file <- readBin(con, 'character')
    if (hash_only) return(hash_file)
    # read number of bytes
    n_bytes <- readBin(con, 'integer')
    # read serialized data
    out_ser <- readBin(con, 'raw', n = n_bytes)
    alloc.col(qs2::qd_deserialize(out_ser))
}



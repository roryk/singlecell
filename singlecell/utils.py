"""Handle file based transactions allowing safe restarts at any point.

To handle interrupts,this defines output files written to temporary
locations during processing and copied to the final location when finished.
This ensures output files will be complete independent of method of
interruption.

lifted from bcbio-nextgen
https://github.com/chapmanb/bcbio-nextgen
"""
import os
import shutil
import tempfile
import contextlib

@contextlib.contextmanager
def file_transaction(*rollback_files):
    """Wrap file generation in a transaction, moving to output if finishes.
    """
    exts = {".vcf": ".idx", ".bam": ".bai", "vcf.gz": ".tbi"}
    safe_names, orig_names = _flatten_plus_safe(rollback_files)
    _remove_files(safe_names)  # remove any half-finished transactions
    try:
        if len(safe_names) == 1:
            yield safe_names[0]
        else:
            yield tuple(safe_names)
    except:  # failure -- delete any temporary files
        _remove_files(safe_names)
        _remove_tmpdirs(safe_names)
        raise
    else:  # worked -- move the temporary files to permanent location
        for safe, orig in zip(safe_names, orig_names):
            if os.path.exists(safe):
                shutil.move(safe, orig)
                for check_ext, check_idx in exts.iteritems():
                    if safe.endswith(check_ext):
                        safe_idx = safe + check_idx
                        if os.path.exists(safe_idx):
                            shutil.move(safe_idx, orig + check_idx)
        _remove_tmpdirs(safe_names)

def _remove_tmpdirs(fnames):
    for x in fnames:
        xdir = os.path.dirname(os.path.abspath(x))
        if xdir and os.path.exists(xdir):
            shutil.rmtree(xdir, ignore_errors=True)

def _remove_files(fnames):
    for x in fnames:
        if x and os.path.exists(x):
            if os.path.isfile(x):
                os.remove(x)
            elif os.path.isdir(x):
                shutil.rmtree(x, ignore_errors=True)

def _flatten_plus_safe(rollback_files):
    """Flatten names of files and create temporary file names.
    """
    tx_files, orig_files = [], []
    for fnames in rollback_files:
        if isinstance(fnames, basestring):
            fnames = [fnames]
        for fname in fnames:
            basedir = safe_makedir(os.path.join(os.path.dirname(fname), "tx"))
            tmpdir = safe_makedir(tempfile.mkdtemp(dir=basedir))
            tx_file = os.path.join(tmpdir, os.path.basename(fname))
            tx_files.append(tx_file)
            orig_files.append(fname)
    return tx_files, orig_files

def safe_makedir(dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not dname:
        return dname
    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dname

def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    return fname and os.path.exists(fname) and os.path.getsize(fname) > 0

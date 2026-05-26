"""Upload FASTQ files to NCBI SRA FTP preload area."""

from ftplib import FTP
from pathlib import Path

import config


def upload_to_sra_ftp(files: list[Path], subfolder: str) -> None:
    """Upload files to SRA FTP under uploads/<username>/<subfolder>/."""
    assert config.SRA_FTP_USERNAME, "Set SRA_FTP_USERNAME in config.py"
    assert config.SRA_FTP_PASSWORD, "Set SRA_FTP_PASSWORD in config.py"

    ftp = FTP(config.SRA_FTP_HOST)
    ftp.login(config.SRA_FTP_USERNAME, config.SRA_FTP_PASSWORD)

    remote_dir = f"{config.SRA_FTP_SUBFOLDER}/{subfolder}"
    _ftp_mkdirs(ftp, remote_dir)
    ftp.cwd(remote_dir)

    for f in files:
        with open(f, "rb") as fh:
            ftp.storbinary(f"STOR {f.name}", fh)

    ftp.quit()


def _ftp_mkdirs(ftp: FTP, path: str) -> None:
    parts = path.strip("/").split("/")
    for part in parts:
        try:
            ftp.cwd(part)
        except Exception:
            ftp.mkd(part)
            ftp.cwd(part)
    # return to root
    ftp.cwd("/")

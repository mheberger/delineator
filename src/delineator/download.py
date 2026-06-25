"""
Utility functions for downloading data files for delineator.
Matthew Heberger, May 2026

"""

import logging
import pooch
from pathlib import Path, PurePosixPath
import requests
from urllib.parse import quote

from delineator.constants import HASHES, DATA_DIR_NAME
from delineator.settings import DelineatorConfig


logger = logging.getLogger(__name__)

BASE_URL = "https://mghydro.com/watersheds/data/"


def _local_path(relative_path: str, config: DelineatorConfig) -> Path:
    r"""
    Gets the user's data directory for delienator's data files.
    This can be a custom path set by the user, if they have set it in the config
    or with an environment variable. Otherwise, it will be a default path.
    I.e. on Windows, it will be something like:
    C:\Users\<USER>\AppData\Local\delineator2.1
    """
    return config.data_dir / DATA_DIR_NAME / relative_path


def _remote_url(relative_path: str) -> str:
    """
    Returns the full download URL for a data file.
    This will be something like:
      https://mghydro.com/watersheds/data/basins12.db
    """

    posix_path = PurePosixPath(Path(relative_path))
    return f"{BASE_URL}{quote(str(posix_path))}"


def _get_remote_file_size(url: str) -> float | None:
    """
    Get the size of a remote file in MB without downloading it.
    Returns None if the server doesn't report a size.
    """
    try:
        response = requests.head(url, timeout=5, allow_redirects=True)
        response.raise_for_status()
        content_length = response.headers.get("Content-Length")
        if content_length:
            return int(content_length) / (1024 * 1024)  # bytes to MB
    except requests.RequestException:
        pass
    return None


def _download_file(relative_path: str, config: DelineatorConfig) -> Path | None:
    """
    Download a data file to the configured data directory.
    If auto-download is disabled or the download fails, return None.
    Uses the `pooch` library for downloading, and

    Parameters
    ----------
    relative_path : str
        The path of the data file relative to the data directory.
        For example: "rivers23.db" or "accum73.tif"

    config : DelineatorConfig dataclass object
        Specifically looking for config.data_dir, where the user
        can set a custom path to the data directory

    Returns
    -------
    a pathlib.Path object pointing to the downloaded file,
    or None if download failed
    """
    if not config.auto_download:
        logger.info("Auto-download is disabled; not downloading %s", relative_path)
        return None

    url = _remote_url(relative_path)
    dest = _local_path(relative_path, config)
    dest.parent.mkdir(parents=True, exist_ok=True)

    file_size_mb = _get_remote_file_size(url)

    if file_size_mb is not None and file_size_mb > 100:
        logger.warning(
            f"Downloading large file ({file_size_mb:.0f} MB): {relative_path}\n"
            f"Destination: {dest}\n"
            f"This may take a few minutes. Set DELINEATOR_AUTO_DOWNLOAD=0 to disable."
        )
    elif file_size_mb is not None:
        logger.info(f"Downloading {relative_path} ({file_size_mb:.0f} MB)...")
    else:
        logger.info(f"Downloading {relative_path}...")

    known_hash = HASHES[relative_path]

    try:
        downloaded = pooch.retrieve(
            url=url,
            known_hash=known_hash,
            fname=dest.name,
            path=dest.parent,
            progressbar=True,
        )
    except Exception as exc:
        logger.warning("Failed to download %s from %s: %s", relative_path, url, exc)
        return None

    return Path(downloaded)

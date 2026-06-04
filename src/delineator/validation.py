"""
Validation function for user-supplied CSV files.
Ensures that the CSV has the required fields and that the id field is unique.
"""

import logging
import pandas as pd
import click

logger = logging.getLogger(__name__)

# Maps common synonyms/variants → canonical field name
COLUMN_SYNONYMS: dict[str, str] = {
    # id
    "id": "id",
    "site_id": "id",
    "station_id": "id",
    "gage_id": "id",
    "gauge_id": "id",
    "reach_id": "id",
    "comid": "id",
    "identifier": "id",
    # lat
    "lat": "lat",
    "latitude": "lat",
    "y": "lat",
    "lat_dd": "lat",
    "latitude_dd": "lat",
    # lon
    "lon": "lon",
    "long": "lon",
    "longitude": "lon",
    "x": "lon",
    "lon_dd": "lon",
    "lng": "lon",
    "longitude_dd": "lon",
}

REQUIRED_FIELDS = {"id", "lat", "lon"}


def _validate_and_normalize_df(df: pd.DataFrame) -> pd.DataFrame | None:
    """Validate a user-supplied DataFrame, normalizing column names to canonical forms.

    Renames columns to lowercase and resolves common synonyms for `id`, `lat`,
    and `lon`. Logs an error and returns None if any required field is missing
    after normalization.

    Args:
        df: Input DataFrame loaded from a user-supplied CSV.

    Returns:
        DataFrame with normalized column names, or None if validation fails.
    """
    # Normalize all column names to lowercase, stripped of whitespace
    df = df.copy()
    df.columns = [c.strip().lower() for c in df.columns]

    # Map recognized synonyms to canonical names (skip columns already canonical)
    rename_map = {
        col: COLUMN_SYNONYMS[col]
        for col in df.columns
        if col in COLUMN_SYNONYMS and COLUMN_SYNONYMS[col] != col
    }
    if rename_map:
        df = df.rename(columns=rename_map)
        logger.debug("Renamed columns: %s", rename_map)

    # Check that all required fields are now present
    missing = REQUIRED_FIELDS - set(df.columns)
    if missing:
        missing_str = ", ".join(sorted(missing))
        present_str = ", ".join(df.columns.tolist())
        msg = f"CSV is missing required field(s): {missing_str}"
        logger.error(msg)
        click.echo(f"Error: {msg}", err=True)
        click.echo(f"  Required : {', '.join(sorted(REQUIRED_FIELDS))}", err=True)
        click.echo(f"  Found    : {present_str}", err=True)
        click.echo(
            "  Tip: column names are case-insensitive. Common synonyms such as\n"
            "  'latitude'/'lat_dd' (→ lat), 'longitude'/'lng' (→ lon), and\n"
            "  'site_id'/'station_id' (→ id) are accepted automatically.",
            err=True,
        )
        return None

    # Check that id values are unique
    duplicated_ids = df.loc[df["id"].duplicated(keep=False), "id"].dropna().unique()
    if len(duplicated_ids) > 0:
        duplicated_str = ", ".join(map(str, duplicated_ids))
        msg = f"CSV field 'id' must contain unique values. Duplicate id value(s): {duplicated_str}"
        logger.error(msg)
        click.echo(f"Error: {msg}", err=True)
        return None

    return df

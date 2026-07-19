"""Build a VRT mosaic over a directory of GRIT raster tiles.

GRIT distributes its 30 m rasters (drainage direction, drainage area, ...)
as 10,000 x 10,000-pixel tiles in EPSG:8857 (Equal Earth) on a common
aligned grid: tile origins are exact multiples of 300 km. That means a VRT
is a lossless mosaic (no resampling), and windowed reads through it are
seamless across tile boundaries — which is what delineator.grit_detailed
relies on. No GDAL CLI is needed; the VRT XML is generated directly.

Usage:
    python build_grit_vrt.py TILE_DIR OUTPUT.vrt [--extent-from OTHER.vrt]

--extent-from forces the output VRT to cover the same extent as another
VRT. Use it to give the drainage-area mosaic the same extent as the flow
direction mosaic: GRIT ships fewer drainage-area tiles than direction tiles
(tiles with no mainstem channels are omitted), and matching extents ensure
that any window that can be read from one mosaic can be read from the other.

After writing, the script verifies the VRT by reading a small window at the
center of every source tile and comparing against the tile itself.
"""
import argparse
import glob
import math
import os
import sys

import numpy as np
import rasterio

GDAL_DTYPES = {
    "uint8": "Byte", "int16": "Int16", "uint16": "UInt16",
    "int32": "Int32", "uint32": "UInt32", "float32": "Float32",
    "float64": "Float64",
}


def build_vrt(tile_dir: str, out_path: str, extent_from: str | None = None) -> None:
    files = sorted(glob.glob(os.path.join(tile_dir, "*.tif")))
    if not files:
        sys.exit(f"No .tif files found in {tile_dir}")
    print(f"{len(files)} tiles in {tile_dir}")

    tiles = []
    for f in files:
        with rasterio.open(f) as src:
            tiles.append((f, src.bounds, src.width, src.height))
    with rasterio.open(files[0]) as src:
        px = src.transform.a
        dtype = src.dtypes[0]
        nodata = src.nodata
        srs = src.crs.wkt
        block_x, block_y = src.block_shapes[0][1], src.block_shapes[0][0]
    gdal_dtype = GDAL_DTYPES[dtype]

    if extent_from:
        with rasterio.open(extent_from) as src:
            xmin, ymin, xmax, ymax = src.bounds
    else:
        xmin = min(b.left for _, b, _, _ in tiles)
        ymax = max(b.top for _, b, _, _ in tiles)
        xmax = max(b.right for _, b, _, _ in tiles)
        ymin = min(b.bottom for _, b, _, _ in tiles)
    W = round((xmax - xmin) / px)
    H = round((ymax - ymin) / px)
    print(f"mosaic extent: {xmin} {ymin} {xmax} {ymax} -> {W} x {H} pixels, "
          f"{gdal_dtype}, nodata={nodata}")

    # GDAL parses "nan" for floating-point nodata
    nodata_str = "nan" if nodata is not None and math.isnan(nodata) else str(nodata)

    # Use paths relative to the VRT when the tiles sit next to it (keeps the
    # mosaic portable if the whole folder moves), absolute paths otherwise.
    same_dir = os.path.abspath(os.path.dirname(os.path.abspath(out_path))) == os.path.abspath(tile_dir)

    parts = [
        f'<VRTDataset rasterXSize="{W}" rasterYSize="{H}">',
        f"  <SRS>{srs}</SRS>",
        f"  <GeoTransform>{xmin}, {px}, 0, {ymax}, 0, {-px}</GeoTransform>",
        f'  <VRTRasterBand dataType="{gdal_dtype}" band="1">',
    ]
    if nodata is not None:
        parts.append(f"    <NoDataValue>{nodata_str}</NoDataValue>")
    for f, b, w, h in tiles:
        xoff = round((b.left - xmin) / px)
        yoff = round((ymax - b.top) / px)
        fname = os.path.basename(f) if same_dir else os.path.abspath(f)
        rel = 1 if same_dir else 0
        source = [
            "    <ComplexSource>",
            f'      <SourceFilename relativeToVRT="{rel}">{fname}</SourceFilename>',
            "      <SourceBand>1</SourceBand>",
            f'      <SourceProperties RasterXSize="{w}" RasterYSize="{h}" '
            f'DataType="{gdal_dtype}" BlockXSize="{block_x}" BlockYSize="{block_y}"/>',
            f'      <SrcRect xOff="0" yOff="0" xSize="{w}" ySize="{h}"/>',
            f'      <DstRect xOff="{xoff}" yOff="{yoff}" xSize="{w}" ySize="{h}"/>',
        ]
        if nodata is not None:
            source.append(f"      <NODATA>{nodata_str}</NODATA>")
        source.append("    </ComplexSource>")
        parts.append("\n".join(source))
    parts += ["  </VRTRasterBand>", "</VRTDataset>"]

    with open(out_path, "w") as fh:
        fh.write("\n".join(parts))
    print(f"wrote {out_path} ({os.path.getsize(out_path) / 1024:.0f} KB)")

    # Verify: a small window at the center of each tile must match the tile.
    with rasterio.open(out_path) as vrt:
        for f, b, w, h in tiles:
            cx, cy = (b.left + b.right) / 2, (b.bottom + b.top) / 2
            bounds = (cx - 16 * px, cy - 16 * px, cx + 16 * px, cy + 16 * px)
            got = vrt.read(1, window=rasterio.windows.from_bounds(*bounds, vrt.transform))
            with rasterio.open(f) as src:
                want = src.read(1, window=rasterio.windows.from_bounds(*bounds, src.transform))
            floating = np.issubdtype(got.dtype, np.floating)
            if not np.array_equal(got, want, equal_nan=floating):
                sys.exit(f"VERIFICATION FAILED at tile {os.path.basename(f)}")
    print(f"verified: center window of all {len(tiles)} tiles reads identically through the VRT")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("tile_dir", help="directory containing the GRIT .tif tiles")
    parser.add_argument("output", help="path of the .vrt file to write")
    parser.add_argument("--extent-from", metavar="OTHER_VRT", default=None,
                        help="force the mosaic extent to match this raster")
    args = parser.parse_args()
    build_vrt(args.tile_dir, args.output, args.extent_from)

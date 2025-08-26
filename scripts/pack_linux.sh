#!/usr/bin/env bash
set -euo pipefail

# Location of the built extension
EXT="build/release/extension/raster/raster.duckdb_extension"
# Output directory for packaged artifacts
OUTDIR="dist/linux"
mkdir -p "$OUTDIR"
cp "$EXT" "$OUTDIR/"

# Directory containing dynamic libraries installed by vcpkg
VCPKG_LIB_DIR="$(pwd)/build/release/vcpkg_installed/x64-linux-dynamic/lib"

# Ensure patchelf is available
if ! command -v patchelf >/dev/null 2>&1; then
  echo "patchelf not found; please install it (apt-get install patchelf or similar)."
  exit 1
fi

# Copy likely dependencies (GDAL and friends) beside the extension.  This list may
# need to be updated when adding new dependencies.
copy_if_present() {
  for n in "$@"; do
    for f in "$VCPKG_LIB_DIR"/$n; do
      [ -f "$f" ] && cp -n "$f" "$OUTDIR/" || true
    done
  done
}
copy_if_present "libgdal*.so*" "libproj*.so*" "libgeos*.so*" "libcurl*.so*" "libtiff*.so*" \
               "libjpeg*.so*" "libpng*.so*" "libwebp*.so*" "libexpat*.so*" "libxml2*.so*" "libz*.so*"

# Set RPATH on extension and copied libs to $ORIGIN so they find each other
patchelf --set-rpath '$ORIGIN' "$OUTDIR/raster.duckdb_extension"
for so in "$OUTDIR"/*.so*; do
  patchelf --set-rpath '$ORIGIN' "$so" || true
done

echo "ldd check:"
ldd "$OUTDIR/raster.duckdb_extension" || true
echo "Linux package at: $OUTDIR"
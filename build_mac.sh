#!/bin/bash
# Build / refresh the Mac app from this source tree.
#   ./build_mac.sh            -> updates ../Clone Bench.app in place (created from macos/ the first time, with the
#                                .command launchers and README-FIRST.txt next to it) and writes a versioned zip there
#   APP=... OUT_DIR=... ./build_mac.sh   -> another bundle location / zip folder
# After building, double-click the app (or run launch.sh) — the launcher notices the new version and swaps it in.
set -e
cd "$(dirname "$0")"
V=$(cat VERSION)
APP="${APP:-../Clone Bench.app}"
if [ ! -d "$APP/Contents" ]; then
  echo "creating $APP from macos/"
  cp -R "macos/Clone Bench.app" "$APP"
  for f in macos/*.command macos/README-FIRST.txt; do cp "$f" "$(dirname "$APP")/"; done
fi
OUT_DIR="${OUT_DIR:-$(dirname "$APP")}"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"   # absolute: the zip step runs in a subshell that changes directory
sed -i.bak "s|<key>CFBundleVersion</key><string>[^<]*</string>|<key>CFBundleVersion</key><string>$V</string>|; s|<key>CFBundleShortVersionString</key><string>[^<]*</string>|<key>CFBundleShortVersionString</key><string>$V</string>|" "$APP/Contents/Info.plist" && rm -f "$APP/Contents/Info.plist.bak"
mkdir -p "$APP/Contents/Resources/app"
rm -rf "$APP/Contents/Resources/app/server" "$APP/Contents/Resources/app/web" "$APP/Contents/Resources/app/examples"
cp -R server web examples requirements.txt README.md VERSION CHANGELOG.md launch.sh "$APP/Contents/Resources/app/"
# keep the launcher icon in step with macos/make_icon.py, so a rebuilt app never keeps an old mark
for ic in AppIcon.icns AppIcon.png; do
  [ -f "macos/Clone Bench.app/Contents/Resources/$ic" ] && cp "macos/Clone Bench.app/Contents/Resources/$ic" "$APP/Contents/Resources/$ic"
done
find "$APP" -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
chmod +x "$APP/Contents/MacOS/"* "$APP/Contents/Resources/app/launch.sh"
echo "app updated: $APP (v$V)"
if command -v zip >/dev/null 2>&1; then
  OUT="$OUT_DIR/Clone Bench v$V.zip"; rm -f "$OUT"
  ( cd "$(dirname "$APP")"
    items=("$(basename "$APP")")
    for f in *.command README-FIRST.txt; do [ -e "$f" ] && items+=("$f"); done
    zip -qr "$OUT" "${items[@]}" )
  echo "zip: $OUT"
fi

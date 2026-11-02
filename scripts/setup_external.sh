#!/usr/bin/env bash
# Clone the external tools into external/ at the commits pinned in DATA.md ("Code dependencies").
# external/ is gitignored; rerunning is safe (existing clones are just checked out at the pinned commit).
set -euo pipefail
P="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"  # repo root
mkdir -p "$P/external"

clone() {  # name url commit
  local dir="$P/external/$1"
  [[ -d $dir/.git ]] || git clone --quiet "$2" "$dir"
  git -C "$dir" fetch --quiet origin "$3" 2>/dev/null || true
  git -C "$dir" checkout --quiet "$3"
  echo "external/$1 @ $(git -C "$dir" rev-parse --short HEAD)"
}

clone burrow              https://github.com/zoidpuff/burrow.git               b2da7abfbbc56d4be4850dfd44e02e55c7b5a648
clone pocketcompare       https://github.com/zoidpuff/pocketcompare.git        d82357e977a08006506ebec59420eb58fb469dd0
clone ancestra-mvp-ai     https://github.com/evocellnet/ancestra-mvp-ai.git    ffb0bbef28bc31d05755264bbe0094e4c0d4a5fa

#!/usr/bin/env bash
set -euo pipefail

export PATH="$HOME/.cargo/bin:$PATH"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

cargo build --workspace --release --offline

mkdir -p "$HOME/.local/bin"
mkdir -p "$HOME/.local/share/applications"
mkdir -p "$HOME/.local/share/icons/hicolor/scalable/apps"

install -m 755 "target/release/chemcraft" "$HOME/.local/bin/chemcraft"
install -m 755 "target/release/chemcraft-cli" "$HOME/.local/bin/chemcraft-cli"
install -m 644 "packaging/chemcraft.desktop" "$HOME/.local/share/applications/chemcraft.desktop"
install -m 644 "packaging/chemcraft.svg" "$HOME/.local/share/icons/hicolor/scalable/apps/chemcraft.svg"

echo "Installed chemcraft and chemcraft-cli into $HOME/.local/bin"

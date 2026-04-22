#!/usr/bin/env bash
set -euo pipefail

echo "=== Codex/VSCode/GitHub repo check ==="

# 1) Czy jesteśmy w repo git?
if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "❌ Nie jesteś w repozytorium Git."
  echo "   Przejdź do katalogu projektu i uruchom ponownie."
  exit 1
fi
echo "✅ Repozytorium Git wykryte."

# 2) Podstawowe dane git
name="$(git config --get user.name || true)"
email="$(git config --get user.email || true)"
if [[ -z "${name}" || -z "${email}" ]]; then
  echo "⚠️ Brak user.name lub user.email w git config."
  echo "   Ustaw:"
  echo "   git config --global user.name \"Twoje Imię\""
  echo "   git config --global user.email \"twoj@email.com\""
else
  echo "✅ git user: ${name} <${email}>"
fi

# 3) Remote origin
if git remote get-url origin >/dev/null 2>&1; then
  origin_url="$(git remote get-url origin)"
  echo "✅ origin: ${origin_url}"
else
  echo "❌ Brak remote 'origin'."
  echo "   Dodaj np.: git remote add origin https://github.com/ORG/REPO.git"
  exit 1
fi

# 4) Aktualna gałąź
branch="$(git branch --show-current || true)"
if [[ -z "${branch}" ]]; then
  echo "⚠️ HEAD detached. Przełącz się na gałąź roboczą."
else
  echo "✅ Aktualna gałąź: ${branch}"
fi

# 5) Stan repo (czystość)
if [[ -n "$(git status --porcelain)" ]]; then
  echo "⚠️ Masz niezacommitowane zmiany."
  git status --short
else
  echo "✅ Working tree clean."
fi

# 6) GitHub CLI auth (opcjonalnie, ale zalecane)
if command -v gh >/dev/null 2>&1; then
  if gh auth status >/dev/null 2>&1; then
    echo "✅ gh auth OK."
  else
    echo "⚠️ gh zainstalowane, ale brak logowania."
    echo "   Uruchom: gh auth login"
  fi
else
  echo "⚠️ gh (GitHub CLI) nie jest zainstalowane."
fi

# 7) Dostęp do remote (fetch)
if git ls-remote --heads origin >/dev/null 2>&1; then
  echo "✅ Dostęp do origin działa (ls-remote OK)."
else
  echo "❌ Brak dostępu do origin (uprawnienia/sieć/URL)."
  exit 1
fi

# 8) Push test (dry-run)
if [[ -n "${branch}" ]]; then
  if git push --dry-run -u origin "${branch}" >/dev/null 2>&1; then
    echo "✅ Push dry-run OK dla gałęzi '${branch}'."
  else
    echo "⚠️ Push dry-run nie przeszedł (możliwy brak uprawnień lub branch protection)."
  fi
fi

echo "=== Check zakończony ==="
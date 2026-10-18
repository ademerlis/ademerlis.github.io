#!/usr/bin/env bash

POST_DATE=$(date '+%Y-%m-%d')

read -r -p "Title: " TITLE
read -r -p "Categories (comma-separated, e.g. Coding, Stress-Hardening): " CATS

SLUG=$(echo "$TITLE" | tr '[:upper:]' '[:lower:]' | sed -E 's/[^a-z0-9]+/-/g; s/^-+|-+$//g')
FILE="posts/${POST_DATE}-${SLUG}.md"

CAT_LIST=$(echo "$CATS" | awk -F',' '{for (i=1; i<=NF; i++) {gsub(/^ +| +$/, "", $i); printf "%s\"%s\"", (i>1 ? ", " : ""), $i}}')

cat > "$FILE" <<EOF
---
title: "${TITLE}"
date: ${POST_DATE}
categories: [${CAT_LIST}]
---

EOF

echo "Created $FILE"
${EDITOR:-nano} "$FILE"

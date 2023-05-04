---
layout: post
title: Optimizing command line on new laptop for bash scripts for downloading raw sequences
date: '2023-05-04'
categories: coding
tags: [coding]
---

Today I received the sequences from UT Austin GSAF via Illumina Basespace. Now I need to download them. 

Michael Studivan provided me with this script (bash not zsh) to run in terminal:

```{bash}
# download BaseSpaceCLI
wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O $HOME/bin/bs

chmod +x ~/bin/bs

# goes to the website to confirm authorization by logging in to your BaseSpace acct.
bs auth

#------------------------------
## Download and concatenate reads with a launcher_creator script

echo '#!/bin/bash' > downloadReads.sh
echo 'bs download project --concurrency=high -q -n ####### -o .' >> downloadReads.sh
# -n is the BaseSpace project name and -o is the output directory

echo "find . -name '*.gz' -exec mv {} . \;" >> downloadReads.sh
echo 'rmdir SA*' >>downloadReads.sh
echo 'mkdir ../concatReads' >> downloadReads.sh
echo 'cp *.gz ../concatReads' >> downloadReads.sh
echo 'cd ../concatReads' >> downloadReads.sh
echo 'for file in *.gz; do mv "${file}" "${file/-2/}"; done' >> downloadReads.sh
# this removes the '-2' in filenames from duplicate samples for downstream merging
```

But first I need to optimize my terminal.
1) download homebrew
2) download wget
3) 

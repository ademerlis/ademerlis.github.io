---
layout: post
title: Optimizing command line on new laptop for bash scripts for downloading raw sequences
date: '2023-05-04'
categories: [Bioinformatics]
tags: [bash, Basespace, Sequencing]
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

End-of-day update: I couldn't get homebrew downloaded onto my computer, every time I tried using the command to install from the homebrew website, it got to the step of needing to update Xcode Command Line Tools and then an error would come up. I think the issue is that the folder hierarchy of the old Xcode doesn't match the file hierarchy of the new Command Line Tools update, so an error comes up before it can finish. I did some googling and maybe the best thing is to uninstall and then re-install Xcode. 

But the whole point of this was to be able to download sequences from BaseSpace, and I found a GUI online to download it from the Project space on Illumina BaseSpace. But then I get an error and all the files fail to download on the external hard drive I'm using. So I'm not sure if maybe using the code Michael provided me is better and will get around the errors from the GUI. However I will need to make some changes to Michael's code because even if I switch wget to curl, the URL in Michael's code doesn't work. Read this link for details on how to download Basespace from terminal: https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview

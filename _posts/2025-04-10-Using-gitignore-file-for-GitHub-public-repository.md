---
layout: post
title: Using gitignore file for GitHub public repository
date: '2025-04-10'
categories: Processing
tags: [Coding, git, GitHub]
---

I recently had to make all of my files for the (temperaturevariabilityAcerPcli)[https://github.com/ademerlis/temperaturevariabilityAcerPcli/tree/main]
publicly available because I was making my repository public for publication, and I needed to transfer a bunch of the files from one computer to the other.

To do this, I deleted my .gitignore file on the old computer which had all the files saved locally, and then pushed everything to GitHub and downloaded it all to my new computer. But then I was having trouble removing those files listed in my .gitignore file again.

The code I ended up running (which worked) in the terminal was:

```{bash}
nano .gitignore
git ls-files -i -c --exclude-from=.gitignore -z | xargs -0 git rm --cached
git status
git commit -m "Stop tracking files listed in .gitignore"
git push
```



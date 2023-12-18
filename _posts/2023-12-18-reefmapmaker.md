---
layout: post
title: installing and using reefmapmaker to create coral reef maps
date: '2023-12-18'
categories: Coding
tags:
  - Coding
published: true
---

I discovered a github repo for a conda package that helps produce high-quality coral reef maps using online data. The package is called [reefmapmaker](https://github.com/didillysquat/reefMapMaker).

Since my Pegasus environment is now conda, I just ran this in the command line after I logged in:

```{bash}
conda install -c didillysquat -c conda-forge reefmapmaker

conda create --name reefmapmaker_env -c didillysquat -c conda-forge reefmapmaker

conda activate reefmapmaker_env
```

I downloaded locally the 2021 UNEP WCMC global distrubiton of warm-water coral reefs dataset from [here](https://data.unep-wcmc.org/datasets/1), then unzipped the file. You need this directory to be in the current working directory for the reefMapMaker commands. I'll create a directory on Pegasus for this and then upload the dataset using scp.

In the installation instructions it says I should be able to run "reefmapmaker" in the directory and it will automatically do it. But i keep getting errors that basically say it can't find the shape file in the dataset I downloaded.

I tried running this and got "Could not automatically find the reference reef dataset. Please specify the directory of the dataset on the command line using --ref-reef-dir"

```{bash}
reefmapmaker --ref-reef-dir 14_001_WCMC008_CoralReefs2021_v4
```

Which one is the shape file?

I think it's this: 01_Data/WCMC008_CoralReef2021_Py_v4_1.shp 

I'm going to download their [config_sheet.tsv](https://github.com/didillysquat/reefMapMaker/blob/main/config_sheet.tsv) and add my sites to a [site_sheet.tsv](https://github.com/didillysquat/reefMapMaker/blob/main/site_sheet.tsv). 

I updated the bounds to fit my data: -84,-80,24,28
and added these sites to my tsv file: KBNursery:	25.6763, -80.0987
PcliCollection	25.77044	-80.15235


attempt #1 to run full code:
```{bash}
reefmapmaker --ref-reef-dir /scratch/projects/and_transcriptomics/programs/reefMapMaker/14_001_WCMC008_CoralReefs2021_v4/01_Data --config-sheet /scratch/projects/and_transcriptomics/programs/reefMapMaker/map/config_sheet.tsv --site-sheet /scratch/projects/and_transcriptomics/programs/reefMapMaker/map/site_sheet.tsv --fig-out-dir /scratch/projects/and_transcriptomics/programs/reefMapMaker/14_001_WCMC008_CoralReefs2021_v4 --bounds=-84,-80,24,28
```

I think I'm running into issues with the config sheet. I'm going to try running the code without it.

```{bash}
reefmapmaker --ref-reef-dir /scratch/projects/and_transcriptomics/programs/reefMapMaker/14_001_WCMC008_CoralReefs2018_v4/01_Data --site-sheet /scratch/projects/and_transcriptomics/programs/reefMapMaker/map/site_sheet.tsv --fig-out-dir /scratch/projects/and_transcriptomics/programs/reefMapMaker/14_001_WCMC008_CoralReefs2021_v4 --bounds=-84,-80,24,28
```

Now it's back to not being able to find the shapefile. I'm going to try to move the zip file onto pegasus first and then unzip it there.

```{bash}
reefmapmaker --ref-reef-dir /scratch/projects/and_transcriptomics/programs/reefMapMaker/14_001_WCMC008_CoralReefs2021_v4_1 --site-sheet /scratch/projects/and_transcriptomics/programs/reefMapMaker/map/site_sheet.tsv --fig-out-dir /scratch/projects/and_transcriptomics/programs/reefMapMaker --bounds=-84,-80,24,28
```

Ok I found the issue. In the second of code of reefmapmaker.py (lines 439-452), it's specifically looking for a shape file that starts with the string "WCMC008_CoralReef2018_Py". However, when you download the dataset from the website link they provide, it's a 2021 version, so all the file names start with "WCMC008_CoralReef2021". ugh. I guess I need to manually edit the reefmapmaker.py code.

Can i do this in the command line? Yes, using sed.

First, navigate to where reefmapmaker was installed: /nethome/and128/anaconda3/envs/reefmapmaker_env/lib/python3.10/site-packages/reefmapmaker/

Then, run:
```{bash}
sed -i 's/WCMC008_CoralReef2018_Py/WCMC008_CoralReef2021_Py/g' reefmapmaker.py

#make executable
chmod +x reefmapmaker.py 
```

Now, navigate to directory which has the WCMC008_CoralReefs2021 data and run:
```{bash}
reefmapmaker --ref-reef-dir /scratch/projects/and_transcriptomics/programs/reefMapMaker/14_001_WCMC008_CoralReefs2021_v4_1 --site-sheet /scratch/projects/and_transcriptomics/programs/reefMapMaker/map/site_sheet.tsv --fig-out-dir /scratch/projects/and_transcriptomics/programs/reefMapMaker --bounds=-84,-80,24,28
```
This didn't work, and neither did running the baseline code: 
```{bash}
reefmapmaker --ref-reef-dir /scratch/projects/and_transcriptomics/programs/reefMapMaker/14_001_WCMC008_CoralReefs2021_v4_1
```

this was the log and the errors that came up:
```{bash}
Shape file found: /scratch/projects/and_transcriptomics/programs/reefMapMaker/14_001_WCMC008_CoralReefs2021_v4_1/01_Data/WCMC008_CoralReef2021_Py_v4_1.shp
Drawing annotations on map

Annotations complete

Annotating reference reefs

reading in reference reefs
done
Traceback (most recent call last):
  File "/nethome/and128/anaconda3/envs/reefmapmaker_env/lib/python3.10/site-packages/reefmapmaker/reefmapmaker.py", line 1002, in _add_reference_reefs
    self._handle_multipolygon(r)
  File "/nethome/and128/anaconda3/envs/reefmapmaker_env/lib/python3.10/site-packages/reefmapmaker/reefmapmaker.py", line 1069, in _handle_multipolygon
    for polygon in r.geometry:
TypeError: 'MultiPolygon' object is not iterable

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/nethome/and128/anaconda3/envs/reefmapmaker_env/bin/reefmapmaker", line 22, in <module>
    mwif.draw_map()
  File "/nethome/and128/anaconda3/envs/reefmapmaker_env/lib/python3.10/site-packages/reefmapmaker/reefmapmaker.py", line 911, in draw_map
    self._add_reference_reefs()
  File "/nethome/and128/anaconda3/envs/reefmapmaker_env/lib/python3.10/site-packages/reefmapmaker/reefmapmaker.py", line 1009, in _add_reference_reefs
    error_count += 1
UnboundLocalError: local variable 'error_count' referenced before assignment
```

I'm going to give up on this for now because I give up. 



# Getting Started
To prevent any conflicts between pip and condas, it is recommended to use MTEX in a fresh virtualized environment (venv). I recommend using python 3.13.x for MTEX for compatability. Once you have Python 3.13 downloaded, you can begin.
```
python3.13 -m venv MTEX
source MTEX/bin/activate
```
You should now have (MTEX) at the beginning of your username in the terminal prompt. Install all packages as follows:
```
pip install "mtex[all]" numpy matplotlib pytictoc
```

# Using the workflow
This process mirrors the MATLAB EBSD processing procedure used by the Taleff group. Images are stored in the folder that the EBSD data resides in, although it is a trivial fix to have a directory specified for saving. No images are shown to the user during the process except the ones that require user interaction. This saves resources by forcing matplotlib to close all figures after being used.
┌─────────────────────────────────────────┐
│        EBSD data is loaded              │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  User specifies material name,          │
│  crystal frame, and plotting convention │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│     Data is gridified for speed         │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  IPF key & map generated from raw EBSD  │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  IQ plot normalized & plotted (raw)     │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Confidence Index plotted (raw)         │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Field Signal Detection plotted (raw)   │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  IQ plot colorized using raw IPF key    │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Low-confidence points forced to        │
│  single known phase                     │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Grains calculated (min angle & size)   │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  User selects mis-identified grains     │
│  (one per parent, iteratively)          │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Mis-identified grains absorbed into    │
│  their parent grains                    │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Copy of raw EBSD made for direct edit  │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  EBSD smoothed (half-quadratic filter)  │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Missing points filled from neighbors   │
│  near boundaries                        │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Grain boundaries smoothed              │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  New IPF map generated from EBSD copy   │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Grain ref. orientation deviation plot  │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Kernel average misorientation plot     │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Final grain boundaries plotted         │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  ODF optimized (optimal kernel          │
│  halfwidth found)                       │
└──────────────────┬──────────────────────┘
                   ▼
┌─────────────────────────────────────────┐
│  Miller indices selected for pole       │
│  figure generation                      │
└─────────────────────────────────────────┘   

This process is set in stone as of now. Future implementations will allow the user to select which processes should be done, allowing for workflow customization.
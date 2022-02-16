## R CMD check results

### Local check
`devtools::check()` result:

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

### Online check
`devtools::check_rhub()` result:

  > checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Jake Hughey <jakejhughey@gmail.com>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    Hughey (8:72)
    Transcriptome (2:42)
    cosinor (7:47)
    transcriptome (8:37)
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.1177/0748730418813785
      From: DESCRIPTION
      Status: Service Unavailable
      Message: 503

  0 errors ✓ | 0 warnings ✓ | 1 note x

Notes:

  - This is the first time this package is being submitted to CRAN, and the words are names and technical terms, so ignore these parts of the first note.
  
  - The second part of the first note (concerning "invalid DOIs") is wrong. The DOI provided is valid, and the link generated correctly navigates to the relevant [doi.](https://doi.org/10.1177/0748730418813785)
  
  - In addition, when checking the package on rhub, the Windows environment produces an error about not being able to install dependencies. This only occurs in rhub, and is due to an error within `install.packages`. See [here.](https://github.com/r-hub/rhub-linux-builders/issues/58#issuecomment-1035506262)

You can also see the results of R CMD check on Windows, Linux, and MacOS by going to the GitHub Actions run labeled `check-deploy` [here](https://github.com/hugheylab/limorhyde/actions).

## Downstream dependencies
There are no downstream dependencies for limorhyde.


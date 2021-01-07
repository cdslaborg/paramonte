
Thank you for your interest in contributing to ParaMonte! Your help is very much appreciated. 
Below are some tips and guidelines to get started. 

## External contributors  

Here is a checklist to help you get started contributing to ParaMonte and walk you through the process,  

+   Take a look at the [issues](https://github.com/cdslaborg/paramonte/issues) page. Make sure that you're not about to duplicate someone else's work.  
+   Post a [new issue](https://github.com/cdslaborg/paramonte/issues/new/choose) discussing the changes you're proposing to implement, 
    whether bug fix(es) or enhancement(s)/feature request(s), 
    or give the rest of developers a heads up that you are going to start work on [an open issue](https://github.com/cdslaborg/paramonte/issues).  
+   [Fork the ParaMonte git project](https://help.github.com/articles/fork-a-repo/) to your private account.  
+   By contributing to the ParaMonte project you are automatically guaranteeing that,  
    1.  The contribution was created in whole or in part by you and you have the right to submit it under the open source license indicated in the file; **or**  
    1.  The contribution is based upon previous work that, to the best of your knowledge, 
    is covered under an appropriate open source license and you have the right under that license to submit that work with modifications, 
    whether created in whole or in part by you, under the same open source license as indicated in the GitHub repository; **or**  
    1.  The contribution was provided directly to me by some other person who guaranteed one of the aforementioned criteria and you have not modified it.  
    1.  You understand and agree that this project and the contribution are public and that a record of the contribution 
        (including all personal information you submit with it, including my sign-off) is maintained indefinitely 
        and may be redistributed consistent with this project or the open source license(s) involved.  
+   Follow the guidelines for [all contributors](#all-contributors) listed below.  

## All contributors  

+   [Create a branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/) and make sure to include the issue number(s) in the branch name, for example: `Provide-binary-OpenMPI-#5` or `fix-issue-#5`.  
+   Make your changes and commit them to your local repository, following these guidelines:  
    +   Each commit should be a logically atomic, self-consistent, cohesive set of changes.  
    +   The code should compile and pass all tests after each commit.  
    +   The code should be legible and any non-obvious features commented appropriately.  
    +   All unit tests should be run locally and pass (see the language-specific guidelines on how to run tests and generate code coverage report). 
    +   Tests should be added for new features and significant new code, 
        steps should be taken to ensure that the total coverage remains the same or increases.
    +   The [commit message](https://robots.thoughtbot.com/5-useful-tips-for-a-better-commit-message) should follow [these guidelines](https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html):  
        +   The first line is the directive phrase, starting with a capitalized imperative verb, and is no longer than 50 characters summarizing your commit.  
        +   The next line, if necessary, is blank.  
        +   The following lines are all wrapped at 72 characters and can include additional paragraphs, bulleted lists, etc.  
        +   Use [Github keywords](https://help.github.com/articles/closing-issues-via-commit-messages/#closing-an-issue-in-a-different-repository), where appropriate, to indicate the commit resolves an open issue.  
        +   Do your best to keep a [clean and coherent history](https://www.notion.so/reviewboard/Keeping-Commit-Histories-Clean-0f717c4e802c4a0ebd852cf9337ce5d2). 
            The commands `git add -p ...`, `git commit --amend` and `git rebase --interactive <root-ref>` 
            can be helpful to rework your commits into a cleaner, clearer state.  
+   Next, [open up a pull request](https://github.com/cdslaborg/paramonte/compare) against the appropriate base branch, [`main` (formerly `master`)](https://github.com/cdslaborg/paramonte/tree/main) of [cdslaborg/paramonte](https://github.com/cdslaborg/paramonte).  
    +   In the title, please include the text `issue-<#>`, where `<#>` is replaced by 
        the issue number of the feature request or bug report corresponding to this pull request (PR).  
    +   If the PR is a work in progress, please add `WIP: ...` to the title, and rename it deleting that text once the PR is ready to be merged.  
    +   If the PR is problematic for any reason please add `DO NOT MERGE` to the title, until it is either abandoned or fixed.  
+   Please be patient and responsive to requests and comments from the PaaraMonte core team members.  
    You may be asked to amend or otherwise alter commits, or push new commits to your branch.  

## Contributors with Write Access  

The ParaMonte core developers and collaborators with push access must wait at least 24 hours before self-approving 
pull requests so that someone else has the chance to review the proposed changes and provide a formal code review. 
Due to the currently small size of the ParaMonte development team, it is unrealistic to *require* a code review under 
all circumstances. This policy ensures that there is at least an opportunity for a formal code review by another developer.  

## The ParaMonte Branches  

The ParaMonte project on GitHub uses the [Github flow](https://guides.github.com/introduction/flow/) workflow. If you are not familiar with [Github flow](https://guides.github.com/introduction/flow/), [this video](https://www.youtube.com/watch?v=EwWZbyjDs9c&feature=youtu.be&list=PLg7s6cbtAD17uAwaZwiykDci_q3te3CTY) might be a good start. The gist of it is that the `main` branch is always deploy-able and deployed. The means at anytime, a new tagged release could be shipped using the `main` branch.

### The main branch  

The `main` branch should remain pristine and stable at all times. Any changes should be applied atomically and exclusively via pull requests. It should be assumed that users are using the code on this branch,
and great care should be taken to ensure its stability. Most bug fixes and incremental improvements will get merged into the `main` branch as soon as they are deemed ready for production.

## General coding style conventions  

The following coding style are enforced within all ParaMonte source files in any programming language. If you do not follows these rules in your contribution, please provide a minimal explanation and justification of the alternative approach that you have taken in your contribution. Most importantly, strict naming conventions are enforced within the entire ParaMonte library.  

+   All names and variables and statements in the library must be self-explanatory to the highest level possible such that minimal comments would be necessary to explain the code behavior.  
    > **WARNING**  
    > Avoid short vague names that hard to decipher for variables and other objects. In particular, avoid single letter variable names, like `x`, `y`, `i`, ... .  
+   [**camelCase**](https://en.wikipedia.org/wiki/Camel_case) writing style is enforced in the entire ParaMonte library (except for constants), like, `sampleSize`, `domainUpperLimitVec`, ... .  Sometimes this convention may be against the common convention used within a particular language, for example, Python or Fortran. However, this deviation from the common practice is needed to bring homogeneity to the ParaMonte library across all programming languages. There are two advantages with using the `camelCase` naming convention:  
    +   The `camelCase` style naturally distinguishes some programming languages' intrinsic entities (for example, Python and Fortran) from the ParaMonte developer’s.  
    +   The `camelCase` style allows extremely long multi-segment variable names within the 63 character limits of many of the programming languages supported in ParaMonte.  
    > **NOTE**  
    > It is understandable that occasionally the `camelCase` style may be hard to follow and enforce in isolated locations in the library. In such cases, it is a good idea to briefly explain the reason for the deviation from the syntax rules of the library where the entity is defined for the first time in the code.  
+   Functions / subroutines / procedures in any programming language always begin with a verb. Example: `getCovarianceMatrix()` or, `getCorMatFromCovMat()`.  
+   All static functions or methods of classes begin with a lowercase verb.
+   Logical functions always begin with `is`. Example: `isDigit()`.  
+   All variables begin with a lower-case character. 
    > **TIP**  
    > An exception to this rule is the ParaMonte kernel routines where non-scalar objects used to begin with an upper-case letter. This old however, is now abandoned in favor of making the first character of all variables lower-case. This is to bring consistency with all other interfaces to the ParaMonte library.  
+   All logical variables must be English propositions that evaluate to either `true` or `false`. Example: `inputFileHasPriority`.  
+   All constants (parameters) or variable that are supposed to not change at runtime must be written upper-case, separated by underscore. Example: `FILE_EXT = ".txt"`.  
    > **NOTE**  
    > Exceptions to this rule sometimes happen in isolated scenarios. In such cases, we recommend that you provide a minimal comment next to the first appearance of the entity explaining why the deviation from the syntax rules of the library was necessary.  
+   The name of any variable that represents a vector of values, that is also anticipated to always represent a vector, is normally suffixed with `Vec`, for example: `startPointVec`, ...
+   The name of any variable that represents a matrix of values, that is also anticipated to always represent a matrix, is normally suffixed with `Mat`, for example: `proposalStartCorMat`, ...
+   The name of any variable that represents a list of varying-size values is normally suffixed with `List`, like: `variableNameList`, ...
+   Separate function arguments with a single space on the left side of the argument. For example,  
    ```python  
    def getCorMatFromCovMat(self, covMat):
    ```  


Thank you for your interest in contributing to ParaMonte! Your help is very much appreciated. 
Below are some tips and guidelines to get started. 

## External contributors  

Here is a checklist to help you get started contributing to ParaMonte and walk you through the process,  

+   Take a look at the [issues](https://github.com/cdslaborg/paramonte/issues) page. Make sure that you're not about to duplicate someone else's work.  
+   Post a [new issue](https://github.com/cdslaborg/paramonte/issues/new/choose) discussing the changes you're proposing to implement, 
    whether bug fix(es) or enhancement(s)/feature request(s) 
    or give the rest of the developers a heads-up that you will start work on [an open issue](https://github.com/cdslaborg/paramonte/issues).  
+   [Fork the ParaMonte git project](https://help.github.com/articles/fork-a-repo/) to your private account.  
+   By contributing to the ParaMonte project, you are automatically guaranteeing that,  
    1.  The contribution was created in whole or in part by you, and you have the right to submit it under the open source license indicated in the file; **or**  
    1.  The contribution is based upon previous work that, to the best of your knowledge, 
    is covered under an appropriate open-source license, and you have the right under that license to submit that work with modifications, 
    whether created in whole or in part by you, under the same open source license as indicated in the GitHub repository; **or**  
    1.  The contribution was provided directly to me by some other person who guaranteed one of the criteria above, and you have not modified it.  
    1.  You understand and agree that this project and the contribution are public and that a record of the contribution 
        (including all personal information you submit with it, including my sign-off) is maintained indefinitely 
        and may be redistributed consistent with this project or the relevant open-source license(s).  
+   Follow the guidelines for [all contributors](#all-contributors) listed below.  

## All contributors  

+   [Create a branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/) and make sure to include the issue number(s) in the branch name, for example: `Provide-binary-OpenMPI-#5` or `fix-issue-#5`.  
+   Make your changes and commit them to your local repository, following these guidelines:  
    +   Each commit should be a logically atomic, self-consistent, cohesive set of changes.  
    +   The code should compile and pass all tests after each commit.  
    +   The code should be legible, and any non-obvious features should be commented appropriately.  
    +   All unit tests should be run locally and pass (see the language-specific guidelines on how to run tests and generate code coverage reports). 
    +   Tests should be added for new features and significant new code; steps should be taken to ensure the total coverage remains the same or increases.
    +   The [commit message](https://robots.thoughtbot.com/5-useful-tips-for-a-better-commit-message) should follow [these guidelines](https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html):  
        +   The first line is the directive phrase, starting with a capitalized imperative verb, and is no longer than 50 characters summarizing your project commits.  
        +   The next line, if necessary, is blank.  
        +   The following lines are all wrapped at 72 characters and can include additional paragraphs, bulleted lists, etc.  
        +   Use [Github keywords](https://help.github.com/articles/closing-issues-via-commit-messages/#closing-an-issue-in-a-different-repository), where appropriate, to indicate the commit resolves an open issue.  
        +   Do your best to keep a [clean and coherent history](https://www.notion.so/reviewboard/Keeping-Commit-Histories-Clean-0f717c4e802c4a0ebd852cf9337ce5d2). 
            The commands `git add -p ...`, `git commit --amend`, and `git rebase --interactive <root-ref>` 
            can be helpful to rework your commits into a cleaner, more transparent state.  
+   Next, [open up a pull request](https://github.com/cdslaborg/paramonte/compare) against the appropriate base branch, [`main` (formerly `master`)](https://github.com/cdslaborg/paramonte/tree/main) of [cdslaborg/paramonte](https://github.com/cdslaborg/paramonte).  
    +   In the title, please include the text `issue-<#>`, where `<#>` is replaced by 
        the issue number of the feature request or bug report corresponding to this pull request (PR).  
    +   If the PR is a work in progress, please add `WIP: ...` to the title and rename it, deleting that text once the PR is ready to be merged.  
    +   If the PR is problematic, please add `DO NOT MERGE` to the title until it is either abandoned or fixed.  
+   Please be patient and responsive to requests and comments from the PaaraMonte core team members.  
    You may be asked to amend, alter, or push new commits to your branch.  

## Contributors with Write Access  

The ParaMonte core developers and collaborators with push access must wait at least 24 hours before self-approving 
pull requests so someone else can review the proposed changes and provide a formal code review. 
Due to the currently small size of the ParaMonte development team, it is unrealistic to *require* a code review under 
all circumstances. This policy ensures at least an opportunity for a formal code review by another developer.  

## The ParaMonte Branches  

The ParaMonte project on GitHub uses the [Github flow](https://guides.github.com/introduction/flow/) workflow. If you are not familiar with [Github flow](https://guides.github.com/introduction/flow/), [this video](https://www.youtube.com/watch?v=EwWZbyjDs9c&feature=youtu.be&list=PLg7s6cbtAD17uAwaZwiykDci_q3te3CTY) might be a good start. The gist is that the `main` branch is always deployable and deployed. A new tagged release could be shipped using the `main` branch at any time.

### The main branch  

The `main` branch should remain pristine and stable at all times. Any changes should be applied atomically and exclusively via pull requests. It should be assumed that users are using the code on this branch,
and great care should be taken to ensure its stability. Most bug fixes and incremental improvements will merge into the `main` branch as soon as they are ready for production.

## General coding style conventions  

The following coding styles are enforced within all ParaMonte source files in any programming language. If you do not follow these rules in your contribution, please provide a minimal explanation and justification of the alternative approach you have taken. Most importantly, strict naming conventions are enforced within the entire ParaMonte library.  

+   All names, variables, and statements in the library must be self-explanatory to the highest level possible, so minimal comments would be necessary to explain the code behavior.  
    > **WARNING**  
    > Avoid short, vague names that are hard to decipher for variables and other objects. In particular, avoid single-letter variable names, like `x`, `y`, `i`, etc.  
+   [**camelCase**](https://en.wikipedia.org/wiki/Camel_case) writing style is enforced in the entire ParaMonte library (except for constants), like, `outputSampleSize`, etc.  Sometimes this convention may be against the standard convention used within a particular language, for example, Python or Fortran. However, this deviation from the standard practice is needed to bring homogeneity to the ParaMonte library across all programming languages. There are two advantages to using the `camelCase` naming convention:  
    +   The `camelCase` style naturally distinguishes some programming languages' intrinsic entities (Python and Fortran) from the ParaMonte developers’.  
    +   The `camelCase` style allows extremely long multi-segment variable names within the 63-character limits of many of the programming languages supported by the ParaMonte library.  
        > **NOTE**  
        > It is understandable that occasionally, the `camelCase` style may be hard to follow and enforce in isolated locations in the library. 
        > In such cases, briefly explain the reason for the deviation from the syntax rules of the library where the entity is defined for the first time in the code.  
+ Functions/subroutines/procedures in any programming language preferably begin with a verb. Example: `getCovarianceMatrix()` or, `getCorMatFromCovMat()`.  
+   All static functions or methods of classes preferably begin with a lowercase verb.
+   Logical functions preferably begin with `is`. Example: `isDigit()`.  
+   All variables begin with a lowercase character. 
+   All logical variables must be English propositions that evaluate either `true` or `false`. Example: `inputFileHasPriority`.  
+   All constants (parameters) or variables that should not change at runtime must be written in upper-case, separated by underscore. Example: `FILE_EXT = ".txt"`.  
    > **NOTE**  
    > Exceptions to this rule sometimes happen in isolated scenarios. In such cases, we recommend that you provide a minimal comment next to the first appearance of the entity explaining why the deviation from the syntax rules of the library was necessary.  
+   Function arguments must be separated with a single space on the left side of the argument. For example,  
    ```python  
    def getCorMatFromCovMat(self, covMat):
    ```  
+   No space should be added immediately after opening brackets or parentheses or before closing them. For example,  
    ```python  
    def getCorMatFromCovMat(self, covMat):
    ```  

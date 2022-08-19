# How to contribute

The way to make contributions is through making pull requests(PR for short). After your PR is merged, the changes you make can be applied by other users.

<div align=center><href="https://sm.ms/image/qt1MWOfYbQzKGJC" target="_blank"><img src="https://s2.loli.net/2022/08/16/qt1MWOfYbQzKGJC.png" width="60%"></div>

Firstly, fork in DP-GEN repository. Then you can clone the repository, build a new branch, make changes and then make a pull request.

---

## How to contribute to DP-GEN

DP-GEN adopts the same convention as other software in DeepModeling Community. 

You can first refer to DeePMD-kit's
[Contributing guide](https://github.com/deepmodeling/deepmd-kit/blob/master/CONTRIBUTING.md)
and [Developer guide](https://docs.deepmodeling.com/projects/deepmd/en/master/development/coding-conventions.html).

You can also read relative chapters on [Github Docs](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request).

Welcome to the repository of [DP-GEN](https://github.com/deepmodeling/dpgen) 

If you have no idea how to fix your problem or where to find the relative code, see [OVERVIEW/Overview](./overview/DP-GEN-structure.html) of the structure of the DP-GEN repository on this website.
  
### Use command line
You can use git with the command line, or open the repository on Github Desktop. Here is a video as a demo of making changes to DP-GEN and publishing it with command line.

<div class="bilitube" data-youtube="DPApc1_eNS8" data-bvid="BV1ZT411c7wp">Loading...</div>
<script src="https://unpkg.com/bilitube@0/dist/bilitube.min.js" defer></script>


> If you have never used Github before, remember to generate your ssh key and configure the public key in Github Settings.
> If you can't configure your username and password, please use token. 
> The explanation from Github see [Github Blog: token authentication requirements for git operations](https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/). [A discussion on StaskOverflow](https://stackoverflow.com/questions/68775869/message-support-for-password-authentication-was-removed-please-use-a-personal) can solve this problem.  

### Use Github Desktop
Also, you can use Github Desktop to make PR.
The following shows the steps to clone the repository and add your doc to tutorials. If it is your first time using Github, Open with Github Desktop is recommended. Github Desktop is a software, which can make your operations on branches visually.

<div align=center><href="https://sm.ms/image/ShdQXosaRM51Jqv" target="_blank"><img src="https://s2.loli.net/2022/08/16/ShdQXosaRM51Jqv.png" width="40%"></div>


After you clone it to your PC, you can open it with Github Desktop.

<div align=center><href="https://sm.ms/image/NMApYxjaqS4DGEz" target="_blank"><img src="https://s2.loli.net/2022/08/16/NMApYxjaqS4DGEz.png" width="40%"></div>

Firstly, create your new branch based on devel branch.

<div align=center><href="https://sm.ms/image/3Eqm162oQ8Lcg9P" target="_blank"><img src="https://s2.loli.net/2022/08/16/3Eqm162oQ8Lcg9P.png" width="40%"></div>

Secondly, add your doc to the certain directory in your local repository, and add its name into index. 

Here is an [example](https://github.com/deepmodeling/tutorials/pull/43). 
Remember to add the filename of your doc into index! 

Thirdly, select the changes that you what to push, and commit to it. Press "Publish branch" to push your origin repository to the remote branch.

<div align=center><href="https://sm.ms/image/3dyQAKplTnR2tX6" target="_blank"><img src="https://s2.loli.net/2022/08/16/3dyQAKplTnR2tX6.png" width="40%"></div>

Finally, you can check it on github and make a pull request. Press "Compare & pull request" to make a PR.

(Note: please commit pr to the devel branch)    

<div align=center><href="https://sm.ms/image/Uj9m6zGtXRh1L3a" target="_blank"><img src="https://s2.loli.net/2022/08/16/Uj9m6zGtXRh1L3a.png" width="80%"></div>

## How to contribute to DP-GEN tutorials and documents
Welcome to [the documents of DP-GEN](https://github.com/deepmodeling/dpgen/tree/master/doc)
- If you want to add the documentation of a toy model, simply put your file in the directory doc/toymodels/ and push;
- If you want to add a new directory for a new category of instructions, make a new directory and add it in doc/index.rst.

Also welcome to [Tutorials repository](https://github.com/deepmodeling/tutorials)
You can find the structure of tutorials and preparations before writing a document in [Writing Tips](https://tutorials.deepmodeling.com/en/devel/Resources/writingTips.html#). 

The latest page of DP-GEN Docs

<div align=center><href="https://sm.ms/image/zEPKuj3TdaHI57b" target="_blank"><img src="https://s2.loli.net/2022/08/16/zEPKuj3TdaHI57b.png" width="60%"></div>

### Examples of contributions
- [Example 1](https://github.com/deepmodeling/dpgen/pull/758) 
- [Example 2](https://github.com/deepmodeling/dpgen/pull/844) (a simple one for beginner)

### 1. Push your doc

<div align=center><href="https://sm.ms/image/T4Zb8uiDXGeyYvc" target="_blank"><img src="https://s2.loli.net/2022/08/16/T4Zb8uiDXGeyYvc.png" width="60%"></div>

### 2. Add the directory in index.rst

<div align=center><href="https://sm.ms/image/q3iKvzQ8oRmfVLt" target="_blank"><img src="https://s2.loli.net/2022/08/16/q3iKvzQ8oRmfVLt.png" width="60%"></div>

### 3. Build and check it

As mentioned in "How to build the website to check if the modification works".

### 4. Make [pull request](https://github.com/deepmodeling/dpgen/pulls) to dpgen

----
  
## Tips
  
1. Please try to submit a PR after finishing all the changes
  
2. Please briefly describe what you do with `git commit -m "what you do"`! "No description provided." will make the maintainer feel confused.

3. It is not recommended to make changes directly in the `devel` branch. It is recommended to pull a branch from devel: `git checkout -b <new-branch-name>`
  
4. When switching branches, remember to check if you want to bring the changes to the next branch!

5. Please fix the errors reported by the unit test. You can firstly test on your local machine before pushing commits. Hint: The way to test the code is to go from the main directory to the tests directory, and use the command `python3 -m unittest`. You can watch the demo video for review. Sometimes you may fail unit tests due to your local circumstance. You can check whether the error reported is related to the part you modified to eliminate this problem. After submitting, as long as there is a green check mark after the PR title on the webpage, it means that the test has been passed.
6. Pay attention to whether there are comments under your PR. If there is a change request, you need to check and modify the code. If there are conflicts, you need to solve them manually.  
---
After successfully making a PR, developers will check it and give comments. It will be merged after everything done. Then CONGRATULATIONS! You become a first-time contributor to DP-GEN!  

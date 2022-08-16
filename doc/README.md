## How to contribute

The way to make contributions is through making pull requests(PR for short). After your PR is merged, the changes you make can be applied by other users.

<div align=center><href="https://sm.ms/image/qt1MWOfYbQzKGJC" target="_blank"><img src="https://s2.loli.net/2022/08/16/qt1MWOfYbQzKGJC.png" width="60%"></div>

Firstly, fork in DP-GEN repository. Then you can clone the repository, build a new branch, make changes and then make a pull request.

---

### How to contribute to DP-GEN
The repository of DP-GEN https://github.com/deepmodeling/dpgen
If you have no idea how to fix your problem or where to find the relative code, see OVERVIEW/Overview of the structure of the DP-GEN repository on this website.
#### Use command line
You can use git with the command line, or open the repository on Github Desktop. Here is a video as a demo of making changes to DP-GEN and publishing it with command line.

https://www.youtube.com/watch?v=DPApc1_eNS8
  
If you have never used Github before, remember to generate your ssh key and configure the public key in Github Settings.
If you can't configure your username and password, please use token. 
The explanation from Github see https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/
Chinese tutorial see https://blog.csdn.net/Saintmm/article/details/119835900
#### Use Github Desktop
Also, you can use Github Desktop to make PR.
The following shows the steps to clone the repository and add your doc to tutorials. If it is your first time using Github, Open with Github Desktop is recommended. Github Desktop is a software, which can make your operations on branches visually.

<div align=center><href="https://sm.ms/image/ShdQXosaRM51Jqv" target="_blank"><img src="https://s2.loli.net/2022/08/16/ShdQXosaRM51Jqv.png" width="40%"></div>


After you clone it to your PC, you can open it with Github Desktop.

<div align=center><href="https://sm.ms/image/NMApYxjaqS4DGEz" target="_blank"><img src="https://s2.loli.net/2022/08/16/NMApYxjaqS4DGEz.png" width="40%"></div>

Firstly, create your new branch based on devel branch.

<div align=center><href="https://sm.ms/image/3Eqm162oQ8Lcg9P" target="_blank"><img src="https://s2.loli.net/2022/08/16/3Eqm162oQ8Lcg9P.png" width="40%"></div>

Secondly, add your doc to the certain directory in your local repository, and add its name into index. 

For example, see https://github.com/deepmodeling/tutorials/pull/43. 
Remember to add the filename of your doc into index! 

This is a case of failed contribution https://github.com/deepmodeling/tutorials/pull/49. 
Without being listed in index, the document will not be shown on the website.

(And here is how it has been fixed https://github.com/deepmodeling/tutorials/pull/50. ) 

Thirdly, select the changes that you what to push, and commit to it. Press "Publish branch" to push your origin repository to the remote branch.

<div align=center><href="https://sm.ms/image/3dyQAKplTnR2tX6" target="_blank"><img src="https://s2.loli.net/2022/08/16/3dyQAKplTnR2tX6.png" width="40%"></div>

Finally, you can check it on github and make a pull request. Press "Compare & pull request" to make a PR.

(Note: please commit pr to the devel branch)    

<div align=center><href="https://sm.ms/image/Uj9m6zGtXRh1L3a" target="_blank"><img src="https://s2.loli.net/2022/08/16/Uj9m6zGtXRh1L3a.png" width="80%"></div>

### How to contribute to DP-GEN tutorials and documents
The documents of DP-GEN https://github.com/deepmodeling/dpgen/tree/master/doc
- If you want to add the documentation of a toy model, simply put your file in the directory doc/toymodels/ and push;
- If you want to add a new directory for a new category of instructions, make a new directory and add it in doc/index.rst.

Tutorials repository: https://github.com/deepmodeling/tutorials
The structure of tutorials and the preparation before writing a document see https://tutorials.deepmodeling.com/en/devel/Resources/writingTips.html#

The latest page of DP-GEN Docs

<div align=center><href="https://sm.ms/image/zEPKuj3TdaHI57b" target="_blank"><img src="https://s2.loli.net/2022/08/16/zEPKuj3TdaHI57b.png" width="60%"></div>

#### Examples of contributions
Example 1(a merged one): https://github.com/deepmodeling/dpgen/pull/758
Example 2(a simple one for beginner): https://github.com/deepmodeling/dpgen/pull/844
#### 1. Push your doc

<div align=center><href="https://sm.ms/image/T4Zb8uiDXGeyYvc" target="_blank"><img src="https://s2.loli.net/2022/08/16/T4Zb8uiDXGeyYvc.png" width="60%"></div>

#### 2. Add the directory in index.rst

<div align=center><href="https://sm.ms/image/q3iKvzQ8oRmfVLt" target="_blank"><img src="https://s2.loli.net/2022/08/16/q3iKvzQ8oRmfVLt.png" width="60%"></div>

#### 3. Build and check it

As mentioned in "How to build the website to check if the modification works".

#### 4. Make pull request to dpgen

https://github.com/deepmodeling/dpgen/pulls

#### How to build the website to check if it works

1. Fork https://github.com/deepmodeling/dpgen
2. Create account on readthedocs
readthedocs.org
https://readthedocs.org/
3. Import https://github.com/<your-github-username>/dpgen.git 
Remember to set Project homepage as https://github.com/<your-github-username>/dpgen , and set the Default branch as the one you do modification in Advanced Settings.

<div align=center><href="https://sm.ms/image/4cVRb7ytT1h5l9n" target="_blank"><img src="https://s2.loli.net/2022/08/16/4cVRb7ytT1h5l9n.png" width="60%"></div>

(“doc_overview” is just an example, please change it into your branch's name.)

4. Build Version 

<div align=center><href="https://sm.ms/image/CsJ7S5xeTYh9EWp" target="_blank"><img src="https://s2.loli.net/2022/08/16/CsJ7S5xeTYh9EWp.png" width="60%"></div>

---
After successfully making a PR, developers will check it and give comments. It will be merged after everything done. Then CONGRATULATIONS! You become a first-time contributor to DP-GEN!  

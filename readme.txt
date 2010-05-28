Population Inference
5/28/2010
tamachado, Columbia University

***

look at run_pi.m for details on how to use these files. 

this code relies upon josh vogelstein's spike inference algorithm (http://github.com/jovo/fast-oopsi) as well as the glmnet/l1 regularization path fitting code by friedman, hastie, and tibshirani (http://www-stat.stanford.edu/~tibs/glmnet-matlab/).

note that the glmnet code has precompiled mex functions, but i've run everything successfully on windows, macos, and linux and so you probably won't have problems.
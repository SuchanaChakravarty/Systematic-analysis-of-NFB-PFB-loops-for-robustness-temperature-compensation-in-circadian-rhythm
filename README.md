# Temperature-compensation-in-negative-and-positive-feedback-oscillators
The codes can be useful to reproduce the results of the research article 'Temperature compensation in negative and positive feedback oscillators'. For our research we have performed the simulations, for 4 different circadian oscillators - 1. Two-Variable-Goodwin-NFB; 2. cyano-KaiABC; 3. cPNFB; and 4. Selkov-PFB. 
The codes are written in Matlab. 
Each of the previously stated network motifs has five folders with estimates for the following scenarios: Arithmatic-Mean calculation; BIC calculation; Temperature compensation with single reaction fixed; Time Course calculation; and Temperature compensation with two reactions fixed. On the other hand, the folder CV-Calculations measures and compares how robust these 4 different circadian oscillators are.
To validate of our findings, we have tested the models, at different parameter regimes. Our investigation has been expanded to include systematic analysis (folder: Selkov-PFB-plus-Extra-NFB-Code) in order to achieve our goal of better understanding the effects of both positive and negative feedback in a Selkov-like positive feedback oscillator with additional negative feedback. The folder titled "Selkov-PFB-plus-Extra-NFB-Code" has two primary files, one for the computation of Arithmatic-Mean and the other for the Time Course.  In the course of our investigation, we looked at three separate scenarios: Strong positive feedback (file: Large-k3-Small-k5); strong negative feedback (file: Small-k3-Large-k5); and a combination of positive and negative feedback (file - Large-k3-Largel-k5).

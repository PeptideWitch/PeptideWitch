# PeptideWitch
For producing high stringency protein identification data from label-free shotgun proteomics results

*****SETUP:*****

1) I highly recommend you download Anaconda: https://www.anaconda.com/download/ and make sure to grab the 3.6 version!!!

Anaconda is a cool python environment tailor-made for scientific application, and makes running packages/modules/scripts like mine a lot smoother. After you download anaconda, install it and make sure to add it to your %path% (it's an option in the install menus). 

2) After that's all done, go ahead and download the master .zip file of PeptideWitch. Unzip it to a spot you like. 

3) Next, go and find something called the anaconda prompt from your menus, open it (it looks like a command prompt), and enter:

<cd (the path that you unzipped PeptideWitch to)>

eg 

<cd C:\David\PepWitch>

4) Next, enter the following:

<pip install -r requirements.txt>

This will read the requirements text file and pull in some specific updates and modules that PeptideWitch requires. If you don't do this, PeptideWitch won't run.

And you're done!

*****RUNNING PEPWITCH:*****

Using the same anaconda command prompt, run PepWitch by navigating to the folder PeptideWitch lives in and entering:

PeptideWitch1.2.py

And voila, PepWitch will run for you :) At the moment, it's a headless little command prompt based program, so to get it to analyse files, you will need to dump in all of your .csv peptide table information into the same folder that PepWitch lives in. Then, simply type the name of the files you want to cluster together.

Let's say you want to analyse:

NCL-Control-R1
NCL-Control-R2
NCL-Control-R3
NCL-Treat8-R1
NCL-Treat8-R2
NCL-Treat8-R3
NCL-Treat16-R1
NCL-Treat16-R2
NCL-Treat16-R3

You would type in:

NCL-Control NCL-Treat8 NCL-Treat16

Or, even simpler:

Control1
Control2
Control3
Wild1
Wild2
Wild3

You can type in:

C W

And PepWitch will know which files you're talking about.

More updates to come :)

17-Nov-17

# PIPE-chipSAD


***	Requirements. ***

The PIPE-chipSAD pipeline is composed by tree Python scripts: chipSAD.py, anno-chipSAD.py and align-chipSAD.py.  Required Python libraries are: wx and igraph.



***  Installation. ***

Just simply untar the package in any destination folder:


         tar zxvf PIPE-chipSAD_v1.0.tgz ./


You will find two folders:


./scripts:	containing the PIPE-chipSAD scripts

./example_data:    containing a simple input dataset



*** Using the PIPE-chipSAD pipeline. ***

* _chipSAD provides the segmentation of the hybridization signal. 


chipSAD could be run either via command line (shell) or using the user-friendly GUI. To run the GUI easily type:


         python ChipSAD_GUI.py


Before running the GUI, the users may want to check if the library wx is already installed.


Or you can use the command line script chipSAD.py.


Input files: 


signal_file - Mvalues for each position of the genome

probes_file - microarrays probes used to perform the experiment.

ttest_table - is included in the “scripts” folder


Examples of “signal” and “probe” files are provided in the example_data folder. 


         Usage: python chipSAD.py signal_file probes_file ttest_table [options]


Options:


      -h, --help            show this help message and exit
  
      -o OUTDIR, --outdir=OUTDIR
  
                        directory for output file. Default is current
                        
                        directory.
                        
      -w WIDTH, --width=WIDTH
  
                        sliding window width. Default 40.
                        
      -s STEP, --step=STEP  step of computation. Default 10.
  
      -t THRESHOLD, --threshold=THRESHOLD
  
                        t-test threshold value. Default SAD compute it from
                        
                        ttest table.
                        
      -p PERCENTILE, --percentile=PERCENTILE
  
                        percentile threshold. Default 8.
                        
      -v PVALUE, --pvalue=PVALUE
  
                        pvalue threshold. Default 0.05. If you set the t
                        
                        parameter different from defaul this options do not
                        
                        will be considered.
                        
      -g GAP, --gap=GAP     the gap widht. Default 240.
  
      -b PSEUDOMEDIAN, --pseudomedian=PSEUDOMEDIAN
  
                        boolean var: 0 for no presmoothing. Default 1.
                        
      -m MINIMUM, --minimum=MINIMUM
  
                        minimun probe number for CPRs. Default 3.
                        
      -j STRAND, --strand=STRAND
  
                        strand. Default F.
                        



Example procedure:


      python chipSAD.py ./example_data/signal_example.txt ./example_data/probes_example.txt ./script/ttest_table.txt –w 4 –s 1 –o ./output


For the example dataset we suggest to use a smaller window size (w) and computational step (s) than the default parameters. Anyway we suggest to try different sizes to find the most suitable ones.


Output, chipSAD provides 6 outputs files:


chipSAD_results*.txt: the complete results of chipSAD application

chipSAD_ttest*.txt: the t value for each position

chipSAD_artemis*.txt: the SAS boundaries to load on a Artemis genome browser

chipSAD_comparison.txt: another file to be loaded on any genome browser showing the SAS

chipSAD_probes*.txt: the chipSAD results for each probe

chipSAD_excel*.txt: the chipSAD results to be loaded on excel for easy screening of the results






***	anno-chipSAD (annotation of signal areas), a method to map the detected transcripts onto the genome. ***



Input files: 


chipsad_results_file - file with the complete results of chipSAD application (chipSAD_results*.txt)

probes_file - microarrays probes used to perform the experiment (provided in example_data folder)

 genome.gbk - the file providing the genome annotation in GenBank format.  
 

To access the program without running chipSAD before, the user must format the input files as required by the program. Example files are provided in the folder “example_data”. 


         Usage: python anno-chipSAD.py chipsad_results_file probes_file genome.gbk [options]



Options:


          -h, --help            show this help message and exit
    
          -o OUTDIR, --outdir=OUTDIR
  
                        directory for output file. Default is current
                        
                        directory.
                        
          -u THUP, --thup=THUP  up threshold for DE transcripts. Default is 1.00.
  
          -d THDOWN, --thdown=THDOWN
  
                        down threshold for DE tramscripts. Default is -1.00.
                        
                        


Output, anno-chipSAD provides 21 outputs files:


tot_results*.txt: genomic location assigned to all chipSAD SAS for both strands

* results *.txt: for each category and for each strand a file with all the belonging SAS is created.

* artemis *.txt: same as results but in a special format to be loaded on Artemis genome browser.







***	align-chipSAD analyses multiple experiments from different chip layouts at the same. ***



Input files: 


a  tot_results*.txt file, output of anno-chipSAD, for each experiments is required. 


These files should be in folder and align-chipSAD should be run in this folder.  To access the program without running the two previous programs before, the user must prepare the input files as tab delimited file with the following information: strand, start position of probe, end position of probe, unique name of transcriptional unit, annotation, start position of transcriptional unit, end position of transcriptional unit, pseudomedian value of transcriptional unit, t value, p value.  Example files are in the folder: example_data/data_for_align-chipSAD.



To run align-chipSAD the library igraph must be installed.


Command: 


         cd folder_input_files

         python align-chipSAD.py [options]



Options:


         -h, --help            show this help message and exit
  
         -o OUTDIR, --outdir=OUTDIR
  
                        directory for output file. Default is current
                        
                        directory.
                        
         -u THUP, --thup=THUP  up threshold for DE tramscripts. Default is 0.00.
  
         -d THDOWN, --thdown=THDOWN
  
                        down threshold for DE tramscripts. Default is 0.00.
                        
         -s STRAND, --strand=STRAND
  
                        set the strand. default F.
                        
                        


Example procedure:


         cd example_data/data_for_align-chipSAD

         python ../../scripts/align-chpSAD.py –o ../



Output, align-chipSAD provides 4 outputs files:


tot_results*.txt: one file for each strand is created. The total number of SAS aligned from each input file is reported.

tot_results*.txt: one file for each strand is created. It bears the boundaries and the intensity of the aligned SAS.



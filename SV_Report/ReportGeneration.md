I wrote a R script that will take the OpenCRAVAT tsv output and generate an report of the SVs of Clinical Significance for each individual. The final output is a multipage pdf with a paragraph explaining the genome position, disease associated, and significance of the SV identified in the individual.It is currently just a quick write up that I jotted down. If there is more information in the OpenCRAVAT info fields that you think should be included, please let me know so I can add it.

In the command line run:
```
Rscript SVeedy_GenerateSVReport.R --tsv CravatReportTest.tsv
```
The CravatReportTest.tsv is the first 10 SVs of the chr1 test with a few obviously made up ClinVar data since there were no SVs linked to diseases in our test. I may have to adjust some of the code when we find true disease causing SVs. 

If you run this with the CravatReportTest.tsv file, you should get a pdf that looks like "Genetic Report of Diseases caused by SVs.pdf".


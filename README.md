# CytoBars

This script will generate one (or more) bar charts of the intensity of each parameter / channel from a flow or mass cytometry dataset (FCS file).

e.g.

<img src="https://raw.githubusercontent.com/JimboMahoney/CytoBars/master/2019-10-11 13_02_58-Window.png"
  align="center" />

<b>Requirements:</b>
 - [R](https://cran.r-project.org/) 
 - An interface such as [RStudio](https://www.rstudio.com/) 
 - some flow or mass FCS data!
 
 
~~In addition, you will need to [install the relevant packages](https://www.datacamp.com/community/tutorials/r-packages-guide) (FlowCore, ggplot2, svdialogs, tidyverse) using the install.packages("packagename") command, which is not included in the script.~~

UPDATE - Added installing (if required) and loading libraries to top of script.

<b>This script will:</b>

1) Read in a specified FCS file opened with a dialogue window.
2) Determine an approximate "cutoff" value for mass cytometry data - this is very roughly the point at which we consider values to be potentially lost in the noise of irrlevant markers - e.g. common contaminants such as I 127 / Pb 208. This is considered to be around 20k DC (dual counts) in solution mode (e.g. 1 second acquisition).
3) Plot a single graph of all markers. 
4) If there are markers below the cutoff, it will show a second. If there are markers above the cutoff, it will plot a third.
5) Colour-code the bars according to the CV of the intensity. e.g. if the intensity is highly varied (spread around from low to high), the bars will show as purple / red.
6) Place a frequency above the bar. e.g. if all events include that parameter / channel, this will be 100%. 
7) Title the chart with the filename, the cutoff value used and the number of markers.

<b> How I learned to do this: </b>

I'm totally new to programming in R and very new to flow / mass cytometry.

The following links were incredibly useful for getting the code up and running:

- A very simple [script](http://rforbiochemists.blogspot.com/2015/07/opening-and-plotting-some-flow.html) to import some FCS data and plot it
- An incredibly complex CyTOF [workflow](https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day2/Workshops/CyTOF/doc/cytofWorkflow_BioC2017workshop.html) with some snippets I could understand and implement
- A <b>lot</b> of [this!](https://www.google.com/)





Feedback / suggestions appreciated.

### See also 

https://github.com/JimboMahoney/CytobankGraphs

https://github.com/JimboMahoney/CyTOF-PlotViewer




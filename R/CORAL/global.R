# setwd(paste0(getwd(),"/R/CORAL"))
#---------------------- SOURCE R FILES ----------------------#

source("R/CORAL/R/colorby.R")
source("R/CORAL/R/readinput.R")
source("R/CORAL/R/writekinasetree.R")
source("R/CORAL/R/legendfunctions.R")
source("R/CORAL/R/map2color.R")
source("R/CORAL/R/convertID.R")
source("R/CORAL/R/makejson.R")
source("R/CORAL/R/colors.R")
source("R/CORAL/R/radiobuttonswithimages.R")

#---------------------- READ IN AND ORGANIZE DATA ----------------------#

# read RDS
orig_svginfo = readRDS("R/CORAL/Data/kintree.RDS")

# if you need to up date the data frame read in the new info here
orig_svginfo$dataframe = data.frame(suppressMessages(readr::read_tsv("R/CORAL/Data/coral_dataframe.tsv")))

# remove NAs from subfamilies
NAs = which(is.na(orig_svginfo$dataframe$kinase.subfamily))
orig_svginfo$dataframe$kinase.subfamily[NAs] = ""

# remove NAs from HGNCs
NAs = which(is.na(orig_svginfo$dataframe$id.HGNC))
orig_svginfo$dataframe$id.HGNC[NAs] = ""

# add correct header
#orig_svginfo$header = "<svg viewBox=\"50 -10 800 640\"  preserveAspectRatio=\"xMidYMid meet\"\n
#orig_svginfo$header = "<svg width=\"940\" height=\"940\"\n
orig_svginfo$header = "<svg width=\"850\" height=\"630\"\n

xmlns=\"http://www.w3.org/2000/svg\"\n
xmlns:xlink=\"http://www.w3.org/1999/xlink\" >\n"

# <defs>\n    <style type=\"text/css\">\n
# @import url('https://fonts.googleapis.com/css?family=Roboto:700');\n
# text {font-family: \"Roboto-Bold\";\n  }\n  </style>\n    </defs>"

# intitialize title
orig_svginfo$title = ""

# initilize legend
orig_svginfo$legend = c()

# add node opacity
orig_svginfo$dataframe$node.opacity = 1

# add node order
orig_svginfo$dataframe$node.selected = -1

# make svginfo (leaving the original intact)
svginfo = orig_svginfo

# assign node and branch orders
svginfo$dataframe$nodeorder = 1:nrow(svginfo$dataframe)
svginfo$dataframe$branchorder = 1:nrow(svginfo$dataframe)

# get example RNA data
rna_data     = paste(readLines("R/CORAL/Data/RNAdata.txt"),collapse="\n")
rna_abs_data = paste(readLines("R/CORAL/Data/RNAdata_pluripotent.txt"),collapse="\n")

#---------------------- DEFAULT COLORS ----------------------#

# Default tree branch color
BG_col1 = "#D3D3D3"

# Default selected color
Cor_col = "#FA6958"

# 2-color heatmap colors
HM2_low = "#1EA0F7"
HM2_hi = "#FACE1E"

# default heatmap colors
HM_low = "#1b8ed1"
HM_med = "#e0e0e0"
HM_hi = "#FA6958"




### Qualtative Palettes ###

qualitative_palette_choices <- c('<img src="images/Erika.png">' = 'Erika',
                               '<img src="images/Accent.png">' = 'Accent',
                               '<img src="images/Dark2.png">' = 'Dark2',
                               '<img src="images/Paired.png">' = 'Paired',
                               '<img src="images/Pastel1.png">' = 'Pastel1',
                               '<img src="images/Pastel2.png">' = 'Pastel2',
                               '<img src="images/Set1.png">' = 'Set1',
                               '<img src="images/Set2.png">' = 'Set2',
                               '<img src="images/Set3.png">' = 'Set3')

# my qualitative palettes
Erika = c("#FA6958","#3F9FFC","#FAD53F","#B0E6C2","#B348A1","#2CD1E0","#BEC956","#7C64FF","#C2374A","#70BD93","#FFBB99","#BA97F2")
Accent = RColorBrewer::brewer.pal(8,"Accent")
Dark2 = RColorBrewer::brewer.pal(8,"Dark2")
Paired = RColorBrewer::brewer.pal(12,"Paired")
Pastel1 = RColorBrewer::brewer.pal(9,"Pastel1")
Pastel2 = RColorBrewer::brewer.pal(8,"Pastel2")
Set1 = RColorBrewer::brewer.pal(9,"Set1")
Set2 = RColorBrewer::brewer.pal(8,"Set2")
Set3 = RColorBrewer::brewer.pal(12,"Set3")

qualpalettes = list(Erika,Accent,Dark2,Paired,Pastel1,Pastel2,Set1,Set2,Set3)
names(qualpalettes) = c("Erika","Accent","Dark2","Paired","Pastel1","Pastel2","Set1","Set2","Set3")

if (! dir.exists('www/images')) {
 dir.create('www/images', showWarnings = F) 
}

# drawmypalettes("Erika",Erika,"www/images",boxes =5)
# drawmypalettes("Accent",Accent,"www/images",boxes =5)
# drawmypalettes("Dark2",Dark2,"www/images",boxes =5)
# drawmypalettes("Paired",Paired,"www/images",boxes =5)
# drawmypalettes("Pastel1",Pastel1,"www/images",boxes =5)
# drawmypalettes("Pastel2",Pastel2,"www/images",boxes =5)
# drawmypalettes("Set1",Set1,"www/images",boxes =5)
# drawmypalettes("Set2",Set2,"www/images",boxes =5)
# drawmypalettes("Set3",Set3,"www/images",boxes =5)

### Sequential Palettes ###

sequential_palette_choices <- c(
 '<img src="images/Greys.png">' = 'Greys',
 '<img src="images/Reds.png">' = 'Reds',
 '<img src="images/Oranges.png">' = 'Oranges',
 '<img src="images/Greens.png">' = 'Greens',
  '<img src="images/Blues.png">' = 'Blues',
 '<img src="images/Purples.png">' = 'Purples')


# my sequential palettes
Greys = RColorBrewer::brewer.pal(3,"Greys")
Reds = RColorBrewer::brewer.pal(3,"Reds")
Oranges = RColorBrewer::brewer.pal(3,"Oranges")
Greens = RColorBrewer::brewer.pal(3,"Greens")
Blues = RColorBrewer::brewer.pal(3,"Blues")
Purples = RColorBrewer::brewer.pal(3,"Purples")

seqpalettes = list(Greys,Reds,Oranges,Greens,Blues,Purples)
names(seqpalettes) = c("Greys","Reds","Oranges","Greens","Blues","Purples")

# drawmypalettes("Greys",Greys,"www/images")
# drawmypalettes("Reds",Reds,"www/images")
# drawmypalettes("Oranges",Oranges,"www/images")
# drawmypalettes("Greens",Greens,"www/images")
# drawmypalettes("Blues",Blues,"www/images")
# drawmypalettes("Purples",Purples,"www/images")

### Divergent Palettes ###

divergent_palette_choices <- c('<img src="images/Blue_Grey_Coral.png">' = 'Blue_Grey_Coral',
                               '<img src="images/Bro_Grey_Tur.png">' = 'Bro_Grey_Tur',
                               '<img src="images/Pink_Grey_Gre.png">' = 'Pink_Grey_Gre',
                               '<img src="images/Pur_Grey_Gre.png">' = 'Pur_Grey_Gre',
                               '<img src="images/Pur_Grey_Or.png">' = 'Pur_Grey_Or',
                               '<img src="images/Red_Grey_Gre.png">' = 'Red_Grey_Gre')

# my divergent palettes
Blue_Grey_Coral = c("#1B8ED1","#e5e5e5","#FA6958")
Bro_Grey_Tur = c("#A6611A","#e5e5e5", "#018571")
Pink_Grey_Gre = c("#D01C8B","#e5e5e5", "#4DAC26")
Pur_Grey_Gre = c("#7B3294","#e5e5e5", "#008837")
Pur_Grey_Or = c("#E66101","#e5e5e5", "#5E3C99")
Red_Grey_Gre = c("#CA0020","#e5e5e5", "#404040")

divpalettes = list(Blue_Grey_Coral,Bro_Grey_Tur,Pink_Grey_Gre,Pur_Grey_Gre,Pur_Grey_Or,Red_Grey_Gre)
names(divpalettes) = c("Blue_Grey_Coral","Bro_Grey_Tur","Pink_Grey_Gre","Pur_Grey_Gre","Pur_Grey_Or","Red_Grey_Gre")

# drawmypalettes("Blue_Grey_Coral",Blue_Grey_Coral,"www/images")
# drawmypalettes("Bro_Grey_Tur",Bro_Grey_Tur,"www/images")
# drawmypalettes("Pink_Grey_Gre",Pink_Grey_Gre,"www/images")
# drawmypalettes("Pur_Grey_Gre",Pur_Grey_Gre,"www/images")
# drawmypalettes("Pur_Grey_Or",Pur_Grey_Or,"www/images")
# drawmypalettes("Red_Grey_Gre",Red_Grey_Gre,"www/images")


# Default group color palette
defaultpalette = colorRampPalette( c(
 "#FA6958",
 "#3F9FFC",
 "#FAD53F",
 "#B0E6C2",
 "#B348A1",
 "#2CD1E0",
 "#BEC956",
 "#7C64FF",
 "#C2374A",
 "#70BD93",
 "#FFBB99",
 "#BA97F2"
))(12)

# # Default group color palette
# defaultpalette = colorRampPalette( RColorBrewer::brewer.pal(11,"Paired"))(11)

# pal_divs = row.names(RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info[,2]=="div",])
# pal_seqs = row.names(RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info[,2]=="seq",])
# pal_quals = row.names(RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info[,2]=="qual",])

CurrentInfoPage = "About"

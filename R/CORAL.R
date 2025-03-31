# CORAL merging file
# All function of CORAL package are merged
# in a single file.
# Link: https://github.com/dphansti/CORAL/


#---------------------- READ IN AND ORGANIZE DATA ----------------------#----

# read RDS
orig_svginfo = readRDS(system.file("extdata", "kintree.RDS", package = "proTN"))

# if you need to up date the data frame read in the new info here
orig_svginfo$dataframe = data.frame(suppressMessages(readr::read_tsv(system.file("extdata", "coral_dataframe.tsv", package = "proTN"))))

# remove NAs from subfamilies
NAs = which(is.na(orig_svginfo$dataframe$kinase.subfamily))
orig_svginfo$dataframe$kinase.subfamily[NAs] = ""

# remove NAs from HGNCs
NAs = which(is.na(orig_svginfo$dataframe$id.HGNC))
orig_svginfo$dataframe$id.HGNC[NAs] = ""

# add correct header
orig_svginfo$header = "<svg width=\"850\" height=\"630\"\n

xmlns=\"http://www.w3.org/2000/svg\"\n
xmlns:xlink=\"http://www.w3.org/1999/xlink\" >\n"

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

#---------------------- DEFAULT COLORS ----------------------#----

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

CurrentInfoPage = "About"

#---------------------- FUNCTIONS----

# Define a function to create a vector of colors based on selected kinases
color.by.selected <- function(df,sel,bg.col,sel.col)
{
  
  # set background color
  color.vector = rep(bg.col,nrow(df))
  
  # recolor selected kinases
  if (length(sel) > 0)
  {
    color.vector[which(df$id.coral %in% sel)] = sel.col
  }
  
  # return color vector
  return (color.vector)
}

# Define a function creates color vector from group
color.by.group <- function(df,recolordf,colors,bg.col="#D3D3D3",categories=NULL)
{
  # set background color
  color.vector = rep(bg.col,nrow(df))
  
  # keep track of group labels
  group.vector = rep("none",nrow(df))
  
  # get group names
  group.names = names(table(recolordf[,2]))
  
  # used user supplied groups
  if (is.null(categories) == FALSE)
  {
    group.names = categories
  }
  
  # Determine the number of groups
  numgroups = length(group.names)
  
  # set palette
  if (numgroups > length(colors))
  {
    colors = colorRampPalette(colors)(numgroups)
  }
  pal = colors[1:numgroups]
  
  groupcolormapping = c()
  for (i in 1:numgroups)
  {
    # get group and color
    group.name  = group.names[i]
    group.color = pal[i]
    
    # record the group/color mapping
    groupcolormapping = c(groupcolormapping,group.color)
    names(groupcolormapping)[length(groupcolormapping)] = group.name
    
    # get kinases from this group
    kinsase.from.this.group = recolordf[which(recolordf[,2]==group.name),1]
    
    # update vector of colors
    color.vector[which(df$id.coral %in% kinsase.from.this.group)] = group.color
    group.vector[which(df$id.coral %in% kinsase.from.this.group)] = group.name
  }
  
  return( list(color.vector,group.vector,groupcolormapping))
}


# Define a function creates color vector from values
color.by.value <- function(df ,recolordf ,colors  ,heatrange , bg.col="#D3D3D3")
{
  # set background color
  color.vector = rep(bg.col,nrow(df))
  
  # kep track of group labels
  value.vector = rep(NA,nrow(df))
  
  # convert to numeric
  recolordf[,2] = as.numeric(recolordf[,2])
  
  # convert to numeric
  recolordf$color = map2color(recolordf[,2],pal = colors, limits = heatrange)
  
  # find indices to recolor
  dflookup = match(recolordf[,1],df[,1])
  
  # update colors and values
  color.vector[dflookup] = recolordf$color
  value.vector[dflookup] = recolordf[,2]
  
  return (list(color.vector, value.vector))  
}


# Define a function creates radius vector from values
resizes.by.value <- function(df, resizedf, sizerange, controlledrange = FALSE, minvalue=0, maxvalue = 5,showall="show")
{
  # set values for non supplied kinases
  radius.vector = rep(0,nrow(df))
  if (showall == "show"){radius.vector = rep(sizerange[1],nrow(df))}
  
  # keep track of group labels
  value.vector = rep(NA,nrow(df))
  
  # convert to numeric
  resizedf[,2] = as.numeric(resizedf[,2])
  
  if (controlledrange == FALSE)
  {
    # (1) get range
    rangesize = sizerange[2] - sizerange[1]
    
    minvalue = min(resizedf[,2])
    maxvalue = max(resizedf[,2])
  }
  
  # if controlledrange == TRUE
  if (controlledrange == TRUE)
  {
    # truncate values beyond the range
    resizedf[,2][which(resizedf[,2] < minvalue)] = minvalue
    resizedf[,2][which(resizedf[,2] > maxvalue)] = maxvalue
    
    # (1) get range
    rangesize = sizerange[2] - sizerange[1]
  }
  
  # (2) shift values such that they start at zero
  radii = resizedf[,2] - minvalue
  
  # (3) scale so max = 1
  radii[which(radii !=0)] = radii[which(radii !=0)] / maxvalue
  
  # (3) multiply to fit range
  radii = radii * rangesize
  
  # (4) increase all values to be within range
  radii = radii+ sizerange[1]
  
  resizedf$radii = radii
  
  # find indices to resize
  dflookup = match(resizedf[,1],df[,1])
  
  # update colors and values
  radius.vector[dflookup] = resizedf$radii
  value.vector[dflookup] = resizedf[,2]
  
  return (list(radius.vector, value.vector))  
}

# Define a function that reads in 2 column text area input
read.text.input <- function(text.input)
{
  # extract user info
  recolordf = data.frame(matrix(unlist(strsplit(x=text.input,split="\\s+")),ncol=2,byrow = T),stringsAsFactors = F)
  colnames(recolordf) = c("kinase","userinfo")
  return (recolordf) 
}


# Define a function that writes the group names
build.group.labels <- function(l,font,groupcolor)
{
  
  # change color
  colortag = paste(" fill=\"",groupcolor,"\" font-family",sep="")
  l = gsub(pattern = "font-family",colortag,l)
  
  # use correct font
  grouplabel = gsub(pattern = "'Roboto-Bold'",font,l)
  
  
  # make bold
  grouplabel = gsub(pattern = "letter-spacing","font-weight=\"700\" letter-spacing",grouplabel)
  
  
  
  
  return(grouplabel)
}

# Define a function that make a branch
build.branch <- function(l)
{
  branch = paste("<path id=\"b_x5F_",l["id.coral"],
                 "\" fill=\"",l["branch.col"],
                 "\" d=\"",l["branch.coords"],"\"/>",sep = "")
  return(branch)
}

# Define a function that make a label
build.text <- function(l,labelselect)
{
  
  # choose label type
  label = ""
  if (labelselect == "Default"){label = l["text.label"]}
  if (labelselect == "coralID"){label = l["id.coral"]}
  if (labelselect == "uniprot"){label = l["id.uniprot"]}
  if (labelselect == "ensembl"){label = l["id.ensembl"]}
  if (labelselect == "entrez"){label = l["id.entrez"]}
  if (labelselect == "HGNC"){label = l["id.HGNC"]}
  
  label = paste("<a xlink:href=\"http://www.uniprot.org/uniprot/",l["id.uniprot"],"\" target=\"_blank\">",
                "<text id=\"t_x5F_",l["id.coral"],"\" ",
                "x=\"",  l["text.x"],"\" ",
                "y=\"", trimws(l["text.y"]),"\" ",
                "font-weight=\"700\" ",
                " font-size=\"",l["text.size"],"px\" ",
                "fill=\"",l["text.col"],"\" ",
                "font-family=\"", l["text.font"], "\" ",
                "letter-spacing=","\".035\"",
                
                # # mouse over effects
                # " \nonmouseover=\"evt.target.setAttribute('font-size', '10');\"",
                # " \nonmouseout=\"evt.target.setAttribute('font-size','",origfontsize,"');\"",
                
                ">",label,"</text>","</a>",sep = "")
  return(label)
}


# Define a function that makes a node
build.node <- function(l)
{
  if (l["node.col"] == "none")
  {
    return()
  }
  
  circle = paste("<circle id=\"n_x5F_",l["id.coral"],"\" ",
                 "cx=\"",l["node.x"],
                 "\" cy=\"",gsub(" ","",l["node.y"]),
                 "\" r=\"",l["node.radius"],
                 "\" opacity=\"",l["node.opacity"],
                 "\" stroke=\"",l["node.strokecol"],
                 "\" stroke-width=\"",l["node.strokewidth"],
                 "\" fill=\"",l["node.col"],"\"/>",sep="")
  return(circle)
}

# Define a function that writes an kinase tree svg file
writekinasetree <- function(svginfo,destination,font,labelselect,groupcolor)
{
  outputlines = c()
  
  # add header
  outputlines = c(outputlines,svginfo$header)
  
  # add title
  outputlines = c(outputlines,"<g id=\"TITLE\">")
  outputlines = c(outputlines,paste("<text x=\"425\" y=\"10\" text-anchor=\"middle\" font-weight=\"700\" font-family=\"",font,"\"  font-size=\"15px\">",svginfo$title,"</text>",sep=""))
  outputlines = c(outputlines,"</g>")
  
  # add legend
  outputlines = c(outputlines,"<g id=\"LEGEND\">")
  outputlines = c(outputlines,svginfo$legend)
  outputlines = c(outputlines,"</g>")
  
  # reorder branches by branch order
  branchorderesDF = svginfo$dataframe[svginfo$dataframe$branchorder,]
  
  # add branches
  outputlines = c(outputlines,"<g id=\"BRANCHES\">")
  outputlines = c(outputlines,unlist(     apply(branchorderesDF,1, build.branch)       ))
  outputlines = c(outputlines,"</g>")
  
  # add circles
  outputlines = c(outputlines,"<g id=\"CIRCLES\">")
  
  # reorder by node order
  nodeorderesDF = svginfo$dataframe[svginfo$dataframe$nodeorder,]
  
  # reorder circles by size
  nodeorderesDF = nodeorderesDF[order(nodeorderesDF$node.radius,decreasing = TRUE),]
  nodeorderesDF = nodeorderesDF[order(nodeorderesDF$node.selected,decreasing = FALSE),]
  
  outputlines = c(outputlines,unlist(apply(nodeorderesDF,1, build.node )))
  outputlines = c(outputlines,"</g>")
  
  # add labels
  outputlines = c(outputlines,"<g id=\"LABELS\">")
  outputlines = c(outputlines,unlist(apply(svginfo$dataframe,1, build.text,labelselect=labelselect)))
  outputlines = c(outputlines,"</g>")
  
  # add tail
  outputlines = c(outputlines,"<g id=\"GROUPS\">")
  outputlines = c(outputlines,unlist(lapply(svginfo$groups, build.group.labels, font=font,groupcolor=groupcolor)))
  outputlines = c(outputlines,"</g>")
  outputlines = c(outputlines,"</svg>")
  
  writeLines(outputlines,destination)
  return(outputlines)
}


# Define a function that makes boxes and labels for group type legends
build.group.legend.elements <- function(x=99.208,y,color,width=6.584,height=6.584,label="group",fontsize=5,elementtype="Branch",fontfamily="'AvenirNext-Bold'")
{
  # build the square
  square = paste("<rect x=\"", x,"\"",
                 " y=\"", y, "\"",
                 " fill=\"", color, "\"",
                 " width=\"", width, "\"",
                 " height=\"", height,"\"/>",
                 sep="")
  
  # build the circle
  circle = paste("<circle cx=\"", x + width/2,"\"",
                 " cy=\"", y+ width/2 , "\"",
                 " fill=\"", color, "\"",
                 " r=\"", width/2, "\"/>",
                 sep="")
  
  
  # build the text
  textx = 110.8889
  texty = y + 4.5
  text = paste("<text x=\"", textx,"\"",
               " y=\"", texty, "\"",
               " font-weight=\"700\" ",
               " font-size=\"", fontsize, "\"",
               " font-family=\"", fontfamily,"\">",
               label,"</text>",
               sep="")
  
  if (elementtype == "Branch")
  {
    return(c(square,text))
  }
  
  if (elementtype == "Node")
  {
    return(c(circle,text))
  }
  
}


# Define a function that builds a legend for group color
build.group.legend  <- function(yoffset=0,groupslabels,groupcolors,elementtype = "Branch",fontfamily="'AvenirNext-Bold'")
{
  # write the header
  header = paste("<text x=\"98.8075\"",
                 " y=\"", yoffset + 8.8451, "\"",
                 " font-family=\"", fontfamily, "\" ",
                 " font-weight=\"700\" ",
                 " font-size=\"9px\">", elementtype," Color</text>",
                 sep="")
  
  # write the grey line
  greylineheight = 14 * length(groupslabels) + 14
  greyline       = paste("<rect x=\"", 89.807,"\"",
                         " y=\"", yoffset, "\"",
                         " fill=\"", "#D3D3D3", "\"",
                         " width=\"", 2.333, "\"",
                         " height=\"", greylineheight,"\"/>",
                         sep="")
  
  # add all of the labels
  legendstuff = c()
  yoffset = yoffset + 19
  for (i in 1:length(groupslabels))
  {
    legendstuff = c(legendstuff,
                    build.group.legend.elements(
                      x=99.208,
                      y=yoffset,
                      color = groupcolors[i],
                      width=6.584,
                      height=6.584,
                      label=groupslabels[i],
                      fontsize=5,
                      fontfamily=fontfamily,
                      elementtype)
    )
    yoffset = yoffset + 14
  }
  
  return (list(c(header,greyline,legendstuff),yoffset))
}


# Define a function that draws a rect
drawrect <- function(x,y,fill,width=6.584,height=11.27)
{
  rectline = paste("<rect x=\"",x,"\"",
                   " y=\"",y,"\"",
                   " fill=\"",fill,"\"",
                   " width=\"",width,"\"",
                   " height=\"",height, "\"/>", 
                   sep = "")
  
  return (rectline)
}


# Define a function that builds a legend for values
build.value.legend  <- function(yoffset=0,minval,maxval, palette,elementtype = "Branch",fontfamily="'AvenirNext-Bold'",subtitle="test")
{
  # write the header
  header = paste("<text x=\"98.8075\"",
                 " y=\"", yoffset + 8.8451, "\"",
                 " font-family=\"", fontfamily, "\" ",
                 " font-weight=\"700\" ",
                 " font-size=\"9px\">", elementtype," Color</text>",
                 sep="")
  
  subtitle.height = 0
  if (subtitle != ""){subtitle.height = 8.8451}
  
  # write the grey line
  greylineheight = 41.58 + subtitle.height
  greyline       = paste("<rect x=\"", 89.807,"\"",
                         " y=\"", yoffset, "\"",
                         " fill=\"", "#D3D3D3", "\"",
                         " width=\"", 2.333, "\"",
                         " height=\"", greylineheight,"\"/>",
                         sep="")
  
  # add the subtitle
  if (subtitle != ""){
    
    yoffset = yoffset + 8.8451
    
    subtitleline = paste("<text x=\"98.8075\"",
                         " y=\"", yoffset + 8.8451*1.5, "\"",
                         " font-family=\"", fontfamily, "\" ",
                         " font-weight=\"700\" ",
                         " font-size=\"7px\">", subtitle,"</text>",
                         sep="")
  }
  
  # add the gradient
  heatrange = seq(minval,maxval,length.out = 11)
  legcols = map2color(x=heatrange,pal=palette,limits=NULL)
  
  # Draw the rectangle
  rects = c()
  for (i in 1:11)
  {
    rects = c(rects, drawrect (x=92.632 + (6.576 * i), y=yoffset + 26.51 ,fill=legcols[i],width=6.584,height=11.27))
  }
  
  text.min = paste("<text x=\"", 98.8075,"\"",
                   " y=\"", yoffset + 23.1251, "\"",
                   " font-size=\"", "5px", "\"",
                   " font-weight=\"700\" ",
                   " font-family=\"", fontfamily, "\">",
                   minval,"</text>",
                   sep="")
  
  text.mid = paste("<text x=\"", 133.8944,"\"",
                   " y=\"", yoffset + 23.1251, "\"",
                   " font-size=\"", "5px", "\"",
                   " font-weight=\"700\" ",
                   " font-family=\"", fontfamily, "\">",
                   mean(c(minval , maxval)),"</text>",
                   sep="")
  
  text.max = paste("<text x=\"", 166.7776,"\"",
                   " y=\"", yoffset + 23.1251, "\"",
                   " font-size=\"", "5px", "\"",
                   " font-weight=\"700\" ",
                   " font-family=\"", fontfamily,"\">",
                   maxval,"</text>",
                   sep="")
  
  # assemble output
  output = c(header, greyline, rects, text.min, text.mid, text.max)
  if (subtitle != "")
  {
    output = c(header, subtitleline, greyline, rects, text.min, text.mid, text.max)
  }
  
  yoffset = yoffset + 41.58
  return(list(output,yoffset))
}


# Define a function that builds a legend for values
build.nodesize.legend  <- function(yoffset=0,minval,maxval,minsize ,maxsize,fontfamily="'AvenirNext-Bold'",subtitle="test")
{
  extrayoff = 0
  
  if (maxsize > 6)
  {
    extrayoff = maxsize - 6
  }
  
  # write the header
  header = paste("<text x=\"98.8075\"",
                 " y=\"", yoffset + 8.8451, "\"",
                 " font-weight=\"700\" ",
                 " font-family=\"", fontfamily, "\" ",
                 " font-size=\"9px\">","Node Size</text>",
                 sep="")
  
  subtitle.height = 0
  if (subtitle != ""){subtitle.height = 8.8451}
  
  # write the grey line
  greylineheight = 41.58 + 2 * extrayoff + subtitle.height
  greyline       = paste("<rect x=\"", 89.807,"\"",
                         " y=\"", yoffset, "\"",
                         " fill=\"", "#D3D3D3", "\"",
                         " width=\"", 2.333, "\"",
                         " height=\"", greylineheight,"\"/>",
                         sep="")
  
  # add the subtitle
  if (subtitle != ""){
    
    yoffset = yoffset + 8.8451
    
    subtitleline = paste("<text x=\"98.8075\"",
                         " y=\"", yoffset + 8.8451*1.5, "\"",
                         " font-family=\"", fontfamily, "\" ",
                         " font-weight=\"700\" ",
                         " font-size=\"7px\">", subtitle,"</text>",
                         sep="")
  }
  
  # Make circles
  circles = c()
  
  xs = c(100.266,109.45,120.846,134.454,150.273,168.303)
  
  sizes = seq(minsize,maxsize,length.out = length(xs))
  
  for (i in 1:length(xs))
  {
    circle = paste("<circle cx=\"", xs[i] ,"\"",
                   " cy=\"", yoffset + 33.932 + extrayoff, "\"",
                   " fill=\"", "#D3D3D3", "\"",
                   " stroke=\"white\"",
                   " r=\"", sizes[i], "\"/>",
                   sep="")
    circles = c(circles,circle)
  }
  
  # add text labels
  text.min = paste("<text x=\"",  min(xs),"\"",
                   " y=\"", yoffset + 23.1251, "\"",
                   " font-size=\"", "5px", "\"",
                   " text-anchor=\"middle\"",
                   " font-weight=\"700\" ",
                   " font-family=\"", fontfamily, "\">",
                   minval,"</text>",
                   sep="")
  text.max = paste("<text x=\"",  max(xs),"\"",
                   " y=\"", yoffset + 23.1251, "\"",
                   " font-size=\"", "5px", "\"",
                   " text-anchor=\"middle\"",
                   " font-weight=\"700\" ",
                   " font-family=\"", fontfamily, "\">",
                   maxval,"</text>",
                   sep="")
  
  # asssemble output
  output = c(header, greyline, circles, text.min, text.max)
  if (subtitle != "")
  {
    output = c(header, greyline,subtitleline,  circles, text.min, text.max)
  }
  yoffset = yoffset + 41.58 + extrayoff
  return(list(output,yoffset))
}


# Define a function that maps numbers to colors
map2color <- function(x=NULL,pal=colorRampPalette(c("deepskyblue2","black","gold"))(100),limits=NULL)
{
  
  # set the limits to the min and max
  if(is.null(limits)) 
  {
    # get just the finite values
    finitevalues = x[which(is.finite(x)==TRUE)]
    
    # set the limits
    limits = range(x)
  }
  
  # get the new colors
  newcolors = pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
  
  # return the new colors
  return(newcolors)
}


# Define a function the converts input identifiers to coral IDs
convertID <- function(df,recolordf,inputtype)
{
  if (inputtype == "coralID")
  {
    inputtype = "id.coral"
  }
  if (inputtype != "id.coral")
  {
    inputtype = paste("id.",inputtype,sep="")
  }
  
  if (length(which(recolordf[,1] %in% df[,inputtype])) == 0)
  {
    message("no matches found")
    return(data.frame())
  }
  
  message(paste(length(which(recolordf[,1] %in% df[,inputtype])),"matches found"))
  
  # filter input df for those found in table
  recolordf = recolordf[which(recolordf[,1] %in% df[,inputtype]),]
  
  # Make a new recolordf accounting for the fact that one other ID could correspond
  # to multilpe coral IDs
  newids = c()
  newvalues = c()
  for (i in 1:nrow(recolordf))
  {
    otherid  = recolordf[i,1]
    quantval = recolordf[i,2]
    coralids = as.character(df[which(as.character(df[,inputtype]) == otherid),1])
    for (coralid in coralids)
    {
      newids = c(newids,coralid)
      newvalues = c(newvalues,quantval)
    }
  }
  recolordfconverted = data.frame(ids=newids,values=newvalues,stringsAsFactors = F)
  
  return(recolordfconverted)
}

# Define a function that shifts a leged in the x direction
shiftlegend <- function (legend,xshift=0,yshift=0,noderadiusexpansion)
{
  NodeSizeSection = FALSE
  
  newlegend = c()
  for (legendtext in legend)
  {
    if (grepl(">Node Size<",legendtext)==TRUE)
    {
      NodeSizeSection = TRUE
    }
    
    # split legend by space
    legendsplit = unlist(strsplit(x=legendtext,split=" "))
    
    # perform the shift
    for (i in 1:length(legendsplit))
    {
      entry = legendsplit[i]
      
      if (startsWith(entry,"x="))
      {
        origvalue = as.numeric(substr(entry,4,nchar(entry)-1))
        newvalue  = origvalue - xshift
        newentry  = paste("x=\"",newvalue,"\"",sep="")
        legendsplit[i] = newentry
      }
      if (startsWith(entry,"y="))
      {
        origvalue = as.numeric(substr(entry,4,nchar(entry)-1))
        newvalue  = origvalue - yshift
        newentry  = paste("y=\"",newvalue,"\"",sep="")
        legendsplit[i] = newentry
      }
      if (startsWith(entry,"cx="))
      {
        origvalue = as.numeric(substr(entry,5,nchar(entry)-1))
        newvalue  = origvalue - xshift
        newentry  = paste("cx=\"",newvalue,"\"",sep="")
        legendsplit[i] = newentry
      }
      if (startsWith(entry,"cy="))
      {
        origvalue = as.numeric(substr(entry,5,nchar(entry)-1))
        newvalue  = origvalue - yshift
        newentry  = paste("cy=\"",newvalue,"\"",sep="")
        legendsplit[i] = newentry
      }
      
      if (  NodeSizeSection == TRUE)
      {
        # scale the node (for force plots)
        if (startsWith(entry,"r="))
        {
          origvalue = substr(entry,4,nchar(entry)-3)
          newvalue  = as.numeric(origvalue) * noderadiusexpansion
          newentry  = paste("r=\"",newvalue,"\"/>",sep="")
          legendsplit[i] = newentry
        }
        
        # move the cirle down just a little
        if (startsWith(legendsplit[i],"cy="))
        {
          extraspace = (noderadiusexpansion-1) * 16
          origvalue = substr(legendsplit[i],5,nchar(legendsplit[i])-1)
          newvalue  = as.numeric(origvalue) + extraspace
          newentry  = paste("cy=\"",newvalue,"\"",sep="")
          legendsplit[i] = newentry
        }
        
        # increase the height of the grey bar
        if (startsWith(legendsplit[i],"height="))
        {
          extraspace = (noderadiusexpansion-1) * 16
          origvalue = substr(legendsplit[i],9,nchar(legendsplit[i])-3)
          newvalue  = as.numeric(origvalue) + extraspace
          newentry  = paste("height=\"",newvalue,"\"/>",sep="")
          legendsplit[i] = newentry
        }
        
      }
      
    }
    
    # recombine
    newlegendtext = paste(legendsplit,collapse=" ")
    
    newlegend = c(newlegend,newlegendtext)
  }
  
  # return shifted legend
  return(newlegend)
}

makejson <- function(df,tmp="www/subdf.txt",output="www/kinome_tree.json",BGcol="#D3D3D3",
                     BGstrolecol="#ffffff",colsubnodes=FALSE,legend="",
                     labelselect,defaultnoderadius,
                     xshift=85,yshift=50,noderadiusexpansion=1)
{
  
  # shift legend to the left
  legend = shiftlegend(legend,xshift=xshift,yshift=yshift,noderadiusexpansion=noderadiusexpansion)
  
  # reorder so so that groups, families, and subfamilies are properly colored
  df<- df[seq(dim(df)[1],1),]
  
  label = ""
  if (labelselect == "Default"){label = "text.label"}
  if (labelselect == "coralID"){label = "id.coral"}
  if (labelselect == "uniprot"){label = "id.uniprot"}
  if (labelselect == "ensembl"){label = "id.ensembl"}
  if (labelselect == "entrez"){label = "id.entrez"}
  if (labelselect == "HGNC"){label = "id.HGNC"}
  
  # filter df
  df = df[,c(label,"kinase.group","kinase.family","kinase.subfamily","branch.col","node.col","node.radius","text.size","node.strokecol","node.opacity")]
  
  # write df to file
  write_tsv(df,tmp,col_names = T)
  
  # read json file
  data<-read.delim(tmp, stringsAsFactors=F)
  root<-list("name"=list("    "), "nodecol"=list(BGcol),"noderadius"=list(defaultnoderadius), nodestrokecol= BGcol,"children"=list(),"legend"=legend)
  i = 1
  
  # print (legend)
  
  for(i in 1:nrow(data)) {
    
    # grab current row
    row<-data[i, ]
    
    # add white spaces
    group<-row$kinase.group
    family<-paste0(" ", row$kinase.family)
    family<-row$kinase.family
    subfamily<-paste0("  ", row$kinase.subfamily)
    kinase<-paste0("   ", row[label])
    branchcol<-row$branch.col 
    nodecol<-row$node.col 
    noderadius<-row$node.radius
    nodestrokecol<-row$node.strokecol
    subnodestrokecol = row$node.strokecol
    nodeopacity<-row$node.opacity
    textsize<-row$text.size
    subnodecol=nodecol
    if (colsubnodes == FALSE)
    {
      subnodecol = BGcol
      subnodestrokecol = BGstrolecol # only highlight the outer nodes
    }
    
    # Add Group if not already there
    g<-match(group, unlist(unlist(root$children, F)[names(unlist(root$children, F))=="name"]))
    if(is.na(g)) {
      root$children[[length(root$children)+1]]<-list("name"=list(group),"branchcol"=list(branchcol) ,"nodecol"=list(subnodecol),"noderadius"=list(defaultnoderadius),"nodestrokecol"=list(subnodestrokecol) , "nodeopacity"=list(nodeopacity) ,"textsize"=list(textsize),"children"=list())
      g<-length(root$children)
    }
    
    # Add Group
    f<-match(family, unlist(unlist(root$children[[g]]$children, F)[names(unlist(root$children[[g]]$children, F))=="name"]))
    if(is.na(f)) {
      root$children[[g]]$children[[length(root$children[[g]]$children)+1]]<-list("name"=list(family),"branchcol"=list(branchcol) ,"nodecol"=list(subnodecol),"noderadius"=list(defaultnoderadius),"nodestrokecol"=list(subnodestrokecol) , "nodeopacity"=list(nodeopacity) , "textsize"=list(textsize),"children"=list())
      f<-length(root$children[[g]]$children)
    }
    
    # Determine whether to skip subfamily or not and add kinase
    if(subfamily == "  ") {
      root$children[[g]]$children[[f]]$children[[length(root$children[[g]]$children[[f]]$children)+1]]<-list("name"=list(kinase),"branchcol"=list(branchcol) ,"nodecol"=list(nodecol),"noderadius"=list(noderadius),"nodestrokecol"=list(nodestrokecol), "nodeopacity"=list(nodeopacity),"textsize"=list(textsize) )
    } else {
      sf<-match(subfamily, unlist(unlist(root$children[[g]]$children[[f]]$children, F)[names(unlist(root$children[[g]]$children[[f]]$children, F))=="name"]))
      if(is.na(sf)) {
        root$children[[g]]$children[[f]]$children[[length(root$children[[g]]$children[[f]]$children)+1]]<-list("name"=list(subfamily),"branchcol"=list(branchcol) ,"nodecol"=list(subnodecol),"noderadius"=list(defaultnoderadius), "nodestrokecol"=list(subnodestrokecol) , "nodeopacity"=list(nodeopacity),"textsize"=list(textsize) , "children"=list())
        sf<-length(root$children[[g]]$children[[f]]$children)
      }
      root$children[[g]]$children[[f]]$children[[sf]]$children[[length(root$children[[g]]$children[[f]]$children[[sf]]$children)+1]]<-list("name"=list(kinase),"branchcol"=list(branchcol) ,"nodecol"=list(nodecol),"noderadius"=list(noderadius),"nodestrokecol"=list(nodestrokecol), "nodeopacity"=list(nodeopacity),"textsize"=list(textsize)  )
    }
  }
  
  # write out json file
  write(toJSON(root, pretty = TRUE), file=output)
  
}

#---------------------- SERVER FUNCTIONS----
# ----------------- MAIN REACTIVE FUNCTION ---------------- #

newdf <- function(comp, kinase_df) {
  
  # get current values
  tempdf = svginfo$dataframe
  
  # set font family
  tempdf$text.font = "'Helvetica'"
  
  # establish legend
  legend = c()
  # Set initial yoffset
  yoffset = 79.125
  
  # get current values
  tempdf$text.size = 4
  
  {
    # read in text area input
    recolordf = data.frame(
      "kinase" = names(kinase_df),
      "userinfo" = kinase_df,
      stringsAsFactors = F
    )
    
    # convert to coral id
    recolordf = convertID (tempdf, recolordf, inputtype = "HGNC")
    if (nrow(recolordf) > 0)
    {
      {
        branchcolorpalette_div <- "Blue_Grey_Coral"
        palettecolors = unlist(divpalettes[branchcolorpalette_div])
        branchcolpalette = colorRampPalette(palettecolors)(11)
        
        bg.col = unlist(divpalettes[branchcolorpalette_div])[2]
      }
      
      # set colors based on group
      newcolors_and_colormapping = color.by.value(
        df = tempdf,
        recolordf = recolordf,
        colors  = branchcolpalette,
        heatrange = c(as.integer(min(
          recolordf$values
        )) - 1, as.integer(max(
          recolordf$values
        )) + 1),
        bg.col = bg.col
      )
      tempdf$branch.col = newcolors_and_colormapping[[1]]
      tempdf$branch.val = newcolors_and_colormapping[[2]]
      
      # reorder based on branch value
      tempdf$branchorder = order(abs(tempdf$branch.val),
                                 decreasing = FALSE,
                                 na.last = FALSE)
      
      # add legend info
      lines_and_offset = build.value.legend(
        yoffset = yoffset,
        minval = as.integer(min(recolordf$values)) - 1,
        maxval = as.integer(max(recolordf$values)) + 1,
        palette = branchcolpalette,
        elementtype = "Branch",
        fontfamily = tempdf$text.font,
        subtitle = "Kinase Activity Score"
      )
      lines = lines_and_offset[[1]]
      yoffset = lines_and_offset[[2]] + 14
      legend = c(legend, lines)
    }
  }
  
  # ------------------ NODE COLOR ------------------ #
  {
    # read in text area input
    recolordf = data.frame(
      "kinase" = names(kinase_df),
      "userinfo" = kinase_df,
      stringsAsFactors = F
    )
    
    # convert to coral id
    recolordf = convertID (tempdf, recolordf, inputtype = "HGNC")
    
    if (nrow(recolordf) > 0)
    {
      {
        nodecolorpalette_div <- "Blue_Grey_Coral"
        palettecolors = unlist(divpalettes[nodecolorpalette_div])
        nodecolpalette = colorRampPalette(palettecolors)(11)
        bg.col = unlist(divpalettes[nodecolorpalette_div])[2]
      }
      {
        bg.col = "#FFFFFF"
      }
      
      # set colors based on value
      newcolors_and_colormapping = color.by.value(
        df = tempdf,
        recolordf = recolordf,
        colors  = nodecolpalette,
        heatrange = c(as.integer(min(
          recolordf$values
        )) - 1, as.integer(max(
          recolordf$values
        )) + 1),
        bg.col = bg.col
      )
      tempdf$node.col = newcolors_and_colormapping[[1]]
      tempdf$node.val = newcolors_and_colormapping[[2]]
      
      # reorder based on branch color
      tempdf$nodeorder = order(abs(tempdf$node.val),
                               decreasing = FALSE,
                               na.last = FALSE)
    }
  }
  
  # ------------------ NODE SIZE ------------------ #
  
  # color nodes by single color
  # if (input$nodesizetype == "One Size")
  {
    tempdf$node.radius = 5
  }
  # ------------------ ADVANCED OPTIONS ------------------ #
  
  tempdf$node.opacity = 1
  
  # Change Color and Size of Font for Selected kinases
  # if (input$fontcolorselect == "Manual")
  {
    selkinases = names(kinase_df)
    
    selkinasescoral = ""
    if (length(selkinases) > 0)
    {
      # convert selected to coral ids
      kinasestoconvert = data.frame(kin1 = selkinases, kin2 = selkinases)
      selkinasesconverted = convertID (tempdf, kinasestoconvert, inputtype =
                                         "HGNC")
      if (nrow(selkinasesconverted) > 0)
      {
        selkinasescoral = selkinasesconverted[, 1]
      }
    }
    
    # set background color and font size
    tempdf$text.col = "#FFFFFF"
    tempdf$text.size = 0.15
    
    # set selected font color and size
    tempdf$text.col[which(tempdf$id.coral %in% selkinasescoral)]  = "#000000"
    tempdf$text.size[which(tempdf$id.coral %in% selkinasescoral)] = 6
  }
  
  # if (input$nodestrokecolselect == "Same as Node")
  {
    tempdf$node.strokecol = tempdf$node.col
  }
  
  return(list(tempdf, legend))
} # end reactive

# ----------------- PLOTS ---------------- #

# build the manning tree
renderSvg <- function(comp, kinase_Act, dirOutput_kinase) {
  # recolor the official matrix
  dfandlegend = newdf(comp, kinase_Act)
  svginfo$dataframe = dfandlegend[[1]]
  svginfo$legend = dfandlegend[[2]]
  # set title
  svginfo$title = paste0("Kinase Tree - Comparison: ", comp)
  
  # Write SVG file
  svgoutfile <- writekinasetree(
    svginfo,
    destination = file.path(dirOutput_kinase, paste0(comp, "_kinase_Tree_CORAL.svg")),
    font = svginfo$dataframe$text.font[1],
    labelselect = "HGNC",
    groupcolor = "#000000"
  )
  # Render SVG
  return(file.path(dirOutput_kinase, paste0(comp, "_kinase_Tree_CORAL.svg")))
  # svgPanZoom(svgoutfile,viewBox = T,controlIconsEnabled=F)
}


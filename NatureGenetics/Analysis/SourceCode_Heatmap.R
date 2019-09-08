library(ggplot2)
library(reshape)
library(plyr)
library(gridExtra)
library(grid)
library(lattice)
library(gplots)

getFrac = function(gff,annotation)
{
    
    frac = function(gs,annot_sites)
    {
        sapply(annot_sites,function(x) {
            max( sapply(gs,function(y) mean(x%in%y)) )
        } )
    }
    
    gff_sites = lapply(gff,function(y) apply(y,1,function(z) z[2]:z[3] ) )
    annot_sites = apply(annotation,1,function(z) z[2]:z[3] )
    fracs = do.call(rbind,lapply( gff_sites , function(gs) frac(gs,annot_sites) ) )
    colnames(fracs) = as.character(annotation$Exon)
    rownames(fracs) = NULL
    fracs = data.frame(fracs)
    
    fracs$exon3b[ fracs$exon3a==1 & fracs$exon3b==1 ] = 0
    fracs$exon7a[ fracs$exon7a==1 & fracs$exon7b==1 ] = 0
    fracs$exon23a[ fracs$exon23a==1 & fracs$exon23b==1 ] = 0
    
    fracs
    
}

heatmapMannualColor2 = function(X,case_col,shared_col)
{
    X.m <- melt(X)
    X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
    gg <- ggplot(X.s, aes(x=variable, y=Name))
    gg <- gg + geom_tile(aes(fill = factor(value)), colour = "white")
    gg <- gg + scale_fill_manual(values=c(case_col,"white",shared_col))
    gg <- gg + labs(x="", y="")
    gg <- gg + theme_bw()
    gg <- gg + theme(panel.grid=element_blank(), panel.border=element_blank())
    base_size <- 9
    gg <- gg + theme(axis.ticks=element_blank(), 
                     axis.text.x=element_text(size=base_size*0.8, angle=300, 
                                              hjust = 0, colour="grey50"))
    print(gg)
}


heatmapMannualColor3 = function(X, case_col="darkred", cont_col="darkgrey", shared_col="darkorange" )
{
    X.m <- melt(X)
    X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
    gg <- ggplot(X.s, aes(x=variable, y=Name))
    gg <- gg + geom_tile(aes(fill = factor(value)), colour = "white")
    gg <- gg + scale_fill_manual(values=c(case_col,"white",shared_col,cont_col))
    gg <- gg + labs(x="", y="")
    gg <- gg + theme_bw()
    gg <- gg + theme(panel.grid=element_blank(), panel.border=element_blank())
    base_size <- 9
    gg <- gg + theme(axis.ticks=element_blank(), 
                     axis.text.x=element_text(size=base_size*0.8, angle=300, 
                                              hjust = 0, colour="grey50"))
    print(gg)
}


heatmapMannualValidateColor = function(X)
{
    X.m <- melt(X)
    X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
    gg <- ggplot(X.s, aes(x=variable, y=Name))
    gg <- gg + geom_tile(aes(fill = factor(value)), colour = "lightgrey")
    gg <- gg + scale_fill_manual(values=c("white","black"))
    gg <- gg + labs(x="", y="")
    gg <- gg + theme_bw()
    gg <- gg + theme(panel.grid=element_blank(), panel.border=element_blank())
    base_size <- 9
    gg <- gg + theme(axis.ticks=element_blank(), 
                     axis.text.x=element_text(size=base_size*0.8, angle=300, 
                                              hjust = 0, colour="grey50"))
    print(gg)
}


heatmapGradientColor = function(X)
{
    X.m <- melt(X)
    X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
    gg <- ggplot(X.s, aes(x=variable, y=Name))
    gg <- gg + geom_tile(aes(fill = value), colour = "white")
    gg <- gg + scale_fill_gradient2(low = "white", high = "darkblue")
    gg <- gg + labs(x="", y="")
    gg <- gg + theme_bw()
    gg <- gg + theme(panel.grid=element_blank(), panel.border=element_blank())
    base_size <- 9
    gg <- gg + theme(axis.ticks=element_blank(), 
                     axis.text.x=element_text(size=base_size*0.8, angle=300, 
                                              hjust = 0, colour="grey50"))
    print(gg)
}


heatmapGradientColor = function(X,Y,title)
{
    X.m <- melt(X)
    X.s <- ddply(X.m, .(variable), transform,rescale = scale(value))
    X.s$Name <- ordered( X.s$Name , rev( levels(X.s$Name) ) )
    gg <- ggplot(X.s, aes(x=variable, y=Name))
    gg <- gg + geom_tile(aes(fill = value), colour = "white")
    gg <- gg + scale_fill_gradient2(low = "white", high = "darkblue")
    gg <- gg + labs(x="", y="")
    gg <- gg + theme_bw()
    gg <- gg + theme(panel.border=element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks.x = element_blank()) #
    base_size <- 9
    gg <- gg + theme(axis.ticks=element_blank(), 
                     axis.text.x=element_text(size=base_size*0.8, angle=300, 
                                              hjust = 0, colour="grey50"),
                     legend.position = "bottom", 
                     legend.direction = "horizontal")
    
    p <- ggplot(data=Y, aes(x=Name, y=y)) + geom_bar(stat="identity") +
        ggtitle(title)
    p.nox <- p + theme(legend.position = "none",
                       axis.title = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks.x = element_blank())
    
    grid.arrange(p.nox, gg, nrow=2, heights=0.7:2.5)
    
}


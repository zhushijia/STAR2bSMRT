parseCigar = function( cigar )
{
	#cigar =  "7742M1I1954M1D1139M1I69M1I552M1I76M1I1964M1I4228M1D2023M1I2288M2D1326M1I1486M1I393M1I427M1I162M2D542M"
	cigar = as.character(cigar)
	pattern = "(\\d+M)|(\\d+I)|(\\d+D)|(\\d+N)|(\\d+S)|(\\d+H)|(\\d+P)|(\\d+=)|(\\d+X)"
	ind = gregexpr( pattern , cigar )[[1]]
	numL = as.numeric(ind)
	numR = numL + attr(ind, "match.length")- 2
	OpI  = numL + attr(ind, "match.length")- 1
	num  = mapply( function(s,e) as.integer(substr( cigar , s , e )) , numL , numR )
	Op   = sapply( OpI , function(s) substr( cigar , s , s ) )
	
	data.frame(num,Op)
	
}

thres = 20
sam <- file("SRR1184043.sam", "r", blocking = FALSE)

num=0
read = readLines(sam,n=1)
while( length(read)>0 )
{
	num = num+1
	cat(read,"\n")
	sep = strsplit(read,"\t")[[1]]
	if(length(sep)<10)
	{
		read = readLines(sam,n=1)
	} else {
		break
	}
}


juncs = list()

while( length(read)>0 )
{
	sep = strsplit(read,"\t")[[1]]
	start = as.integer(sep[4])
	cigar = parseCigar(sep[6])
	cigar = subset( cigar , !Op%in%c("I","S") )
	index = which( cigar[,1]>thres & as.character(cigar[,2]) %in% c("N","D") )
	
	juncL = c()
	juncR = c()

	for(i in 1:length(index))
	{
		juncL[i] = start + sum(cigar$num[1:(index[i]-1)])
		juncR[i] = start + sum(cigar$num[1:index[i]]) - 1
	}
	
	junc = paste(apply(data.frame(juncL,juncR),1,function(x)paste(x,collapse=",")) , collapse="," )
	juncs = c(juncs,junc)
	
	if(0)
	{
		junc1 = gsub( "^jI:B:i,","",sep[19] )
		junc2 = paste(apply(data.frame(juncL,juncR),1,function(x)paste(x,collapse=",")) , collapse="," )
		
		num = num+1
		cat(num,junc1==junc2,"\n")
		
		if( junc1!=junc2 )
		break
	}
	read = readLines(sam,n=1)
	read
}

close(sam)

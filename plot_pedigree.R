library(shape)

###############################################################
#get pedigree file name:
args=commandArgs(TRUE)

if(length(args)==0) {
	stop(error="\nUsage:\tRscript plot_pedigree.R <input pedigree file>\n",call. = FALSE)
}else{
	file=args[1]
}

###############################################################
#read funcitons:
source("functions.R")

###############################################################
#output name:
description=gsub(".*/|_ped.txt","",file)
out=gsub("txt","pdf",file)

###############################################################
#read pedigree file:
ped = read.table(file, stringsAsFactors=F)
header=scan(	file, 
							nlines=1,
							what=character(),
							quiet=T,
							sep="\t",
							skip = 0
						)
											
colnames(ped) =header[1:dim(ped)[2]]
colnames(ped)[1] = gsub("#","",colnames(ped)[1] )

#################################################################
#define all parents:
ped=get_all_parents(ped$ID,ped)

###############################################################
#fill up pedigree (every indiviual has a mother and a father):
all_ids=unique(c(ped$ID,ped$MOTHER,ped$FATHER))
coords=data.frame(ID=all_ids,x=rep(NA,length(all_ids)),y=rep(NA,length(all_ids)),Gender=rep(NA,length(all_ids)),FAM_ID =rep(NA,length(all_ids)),stringsAsFactors=F)

#get Gender:
for(i in 1:dim(coords)[1]){
	if(coords$ID[i] %in% ped$ID)	{
		if(ped$GENDER[which(ped$ID == coords$ID[i])] != "male" & ped$GENDER[which(ped$ID == coords$ID[i])] != "female"){
			coords$Gender[i] = "undefined"
		}else{
			coords$Gender[i] = ped$GENDER[which(ped$ID == coords$ID[i])]
		}
	}else{
		coords$Gender[i] = gsub(".*_","",coords$ID[i] )
	}
}

#print(ped)
#go through the graph:
fam = 1	#family ID
while(is.na(min(coords$x))){
	ids=intersect(ped$ID,coords$ID[is.na(coords$x)])
	
	ped_tmp=ped[intersect(c(1:length(ped$ID)),which(is.na(coords$x))),]
	this_id=ids[which(ped_tmp$CHILDREN == "-" & ped_tmp$PARENT_CHILD_RELATIONSHIP!= "-")][1]
	if(is.na(this_id)){
		this_id=ids[1]
	}
	coords$x[which(coords$ID == this_id)] = 0
	coords$y[which(coords$ID == this_id)] = 0
	coords=get_coords(coords,ped,ids,this_id)
	
	coords$FAM_ID[which(!is.na(coords$x) & is.na(coords$FAM_ID)) ] = fam
	fam = fam+1
}
final_coords=coords
print(coords)
# #############################
#plot:
#############################
#setEPS()
#postscript(file="FIGURES/Precision_vs_VariantCount.eps", width =20,height=8)
pdf(out,8,8)

#############################
#plot edges:

for(f in 1:(fam-1)){
	coords = final_coords[which(final_coords$FAM_ID == f), ]

	dist=0.2

	coords$x=as.numeric(coords$x)
	coords$y=as.numeric(coords$y)

	x_left=min(coords$x)-dist
	x_right=max(coords$x)+dist
	y_bottom=min(coords$y)-dist
	y_top=max(coords$y)+dist

	plot(x=coords$x,y=coords$y,type="n",xlab="",ylab="",axes=F,xlim=c(x_left,x_right),ylim=c(y_bottom,y_top),main=description)

	for(i in 1:dim(coords)[1]){
		index=which(final_coords$ID == coords$ID[i])
		x1=coords$x[i]
		y1=coords$y[i]
		y2=coords$y[which(coords$ID == ped$MOTHER[index])]		#equals y2 of Father
	
		if(coords$ID[i] %in% ped$ID){		
			siblings=c()
			if(ped$SIBLINGS[index] != "-"){
				siblings=unlist(strsplit(as.character(ped$SIBLINGS[index]),","))
			}
		
			x2_mother=coords$x[which(coords$ID == ped$MOTHER[index])]
			x2_father=coords$x[which(coords$ID == ped$FATHER[index])]

			y_mid=y2
			lines(x=c(x2_mother,x2_father),y=c(y2,y2))
		
			if(length(siblings)== 0){
				lines(x=c(x1,x1),y=c(y1,y_mid))
			}else{	
				y_mid=y2-0.75	
			
				x_sib_coords=c(x1)
				for(j in 1:length(siblings)){
					x_sib=coords$x[which(coords$ID == siblings[j])]		
					x_sib_coords = c(x_sib_coords, x_sib)
					lines(x=c(x1,x1),y=c(y_mid,y1))
					lines(x=c(x1,x_sib),y=c(y_mid,y_mid))
					lines(x=c(x_sib,x_sib),y=c(y_mid,y1))
				}
				lines(x=c(mean(x_sib_coords),mean(x_sib_coords)),y=c(y_mid,y2))		
			}	
		}
	}

	#############################
	#plot plot circles/rectangles:

	for(i in 1:dim(coords)[1]){
			index=which(final_coords$ID == coords$ID[i])

			 x_left=coords$x[i]-dist
			 x_right=coords$x[i]+dist
			 y_bottom=coords$y[i]-dist
			 y_top=coords$y[i]+dist
		 
			 if(coords$Gender[i] == "male"){
			 		rect(x_left,y_bottom,x_right,y_top,col="white",border="dimgray")
			 }else if(coords$Gender[i] == "female"){
					plotcircle(mid = c(coords$x[i], coords$y[i]), r = dist,col="white", lcol ="dimgray",lwd=1)
			}else{
					polygon(x=c(x_left,(x_left+(x_right-x_left)/2),x_right,(x_left+(x_right-x_left)/2)),
									y=c((y_bottom+(y_top-y_bottom)/2), y_bottom,(y_bottom+(y_top-y_bottom)/2), y_top),col="white")	
			}
	
		if(coords$ID[i] %in% ped$ID){
			text(x=(coords$x[i]),y=coords$y[i],labels=coords$ID[i],cex=0.8)
		}
	}
}
dev.off()

 ###############################################################
#fill up pedigree (every indiviual has a mother and a father):
get_all_parents = function(ids, ped){
	
	for(i in 1:length(ids)){
		siblings=c(ids[i])
		if(ped$SIBLINGS[which(ped$ID == ids[i])] != "-"){
			siblings=c(siblings,unlist(strsplit(as.character(ped$SIBLINGS[which(ped$ID == ids[i])]),",")))
		}
	
		for(j in 1:length(siblings)){
			if(ped$MOTHER[which(ped$ID == siblings[j])] == "-"){
				m=paste(i,"_female",sep="")
				ped$MOTHER[which(ped$ID == siblings[j])]  = m
			}
		
			if(ped$FATHER[which(ped$ID == siblings[j])] == "-"){
				f=paste(i,"_male",sep="")	
				ped$FATHER[which(ped$ID == siblings[j])]  = f
			}
		}
	}
	return(ped)
}

###############################################################
#check if relation exist:
check_if_exist=function(ped, id, relation){
	index=which(colnames(ped) ==relation)
	ret = FALSE
	
	if(length(index) >0 & length(which(ped$ID ==id))>0 ){
		if(ped[which(ped$ID ==id),index] != "-"){
			ret = TRUE
		}else{
			ret = FALSE
		}
	}
	return(ret)
}

###############################################################
#get coordinates for a family:
get_coords=function(coords,ped, ids, this_id){

	all_ids = ids
	while(length(which(ids!="-")) != 0){
	  
		if(!this_id %in% ids | this_id=="-"){
				this_id=ids[which(ids != "-")[1]]
		}
	
		if(is.na(coords$x[which(coords$ID ==this_id)])){
			break
		}
	  
	  message("-----------")
	  message(this_id)
	 # message(coords[,1:4])
	  
#print(this_id)
		#children?
		children=c()
		if(check_if_exist(ped,this_id,"CHILDREN")){
				children=unlist(strsplit(as.character(ped$CHILDREN[which(coords$ID == this_id)]),","))
				
				dist_x_coords_children = 0.5
				if(length(children) >1){
					dist_x_coords_children = seq(0,1,length.out=length(children))
				}
				
				#adapt coordinates, if some children are already defined
				for(i in 1:length(children)){
						if(!is.na(coords$x[which(coords$ID == children[i])])){
							dist_x_coords_children = dist_x_coords_children-coords$x[which(coords$ID == children[i])]
							break
						}
				}

				for(i in 1:length(children)){
					if(is.na(coords$x[which(coords$ID == children[i])])){
						coords$x[which(coords$ID == children[i])] = as.numeric(coords$x[which(coords$ID == this_id)])+dist_x_coords_children[i]
						coords$y[which(coords$ID == children[i])] = as.numeric(coords$y[which(coords$ID == this_id)])-1
			
						#proceed with children:
						ids[which(ids==children[i])] = "-"
						ids = c(children[i],ids)	
					}
				}
		}
	
		#siblings?
		siblings=c()
		siblings_idx=c(which(coords$ID ==this_id))	
		if(check_if_exist(ped,this_id,"SIBLINGS")){
				siblings=unlist(strsplit(as.character(ped$SIBLINGS[which(coords$ID == this_id)]),","))
			
				for(i in 1:length(siblings)){
					sib=which(coords$ID ==siblings[i])
					siblings_idx=c(siblings_idx,sib)

					if(is.na(coords$x[which(coords$ID == siblings[i])])){
						dist_x_coords_siblings=1
					
						while(is.na(coords$x[which(coords$ID == siblings[i])])){
						
							#define new coordinates:
							new_x= as.numeric(coords$x[which(coords$ID == this_id)])+dist_x_coords_siblings
							new_y= as.numeric(coords$y[which(coords$ID == this_id)])
					
							#check if cocordinates are already taken:
							if(length(which(coords$x == new_x & coords$y ==new_y )) >0){
								dist_x_coords_siblings = dist_x_coords_siblings+1
							}else{
								coords$x[which(coords$ID == siblings[i])] = new_x
								coords$y[which(coords$ID == siblings[i])] = new_y
							}
						}
				
						#proceed with siblings:
						ids[which(ids==siblings[i])] = "-"
						ids = c(siblings[i],ids)	
					}
				}
		}
	
		#define parents coordinates:
		mother=as.character(ped$MOTHER[which(coords$ID == this_id)])
		father=as.character(ped$FATHER[which(coords$ID == this_id)])

		mean_x_coords=mean(as.numeric(coords$x[siblings_idx]))
		new_y= as.numeric(coords$y[which(coords$ID == this_id)])+1
		new_x= mean_x_coords

 		dist_x_coords_parent = 0.5	
 		#change distance if parents have both parents:
 		if(check_if_exist(ped,mother,"MOTHER") & check_if_exist(ped,father,"MOTHER")){
 			dist_x_coords_parent = 1
 		}
 		
		#define new x-coordinates:
		 if(is.na(coords$x[which(coords$ID == mother)])){
			coords$y[which(coords$ID == mother)] = new_y

		 	if(!is.na(coords$x[which(coords$ID == father)])){		 		
				coords$x[which(coords$ID == mother)] = coords$x[which(coords$ID == father)]+dist_x_coords_parent*2
		 	}else{
				coords$y[which(coords$ID == father)] = new_y
				
				coords$x[which(coords$ID == mother)] = new_x-dist_x_coords_parent			
				coords$x[which(coords$ID == father)] = (new_x-dist_x_coords_parent+dist_x_coords_parent*2)
			}
		 }else{
		 		if(is.na(coords$x[which(coords$ID == father)])){
					coords$y[which(coords$ID == father)] = new_y
					coords$x[which(coords$ID == father)] = coords$x[which(coords$ID == mother)]+dist_x_coords_parent*2
				}
		 }
	
		#check if cocordinates are already taken:
		# if(length(which(coords$x == coords$x[which(coords$ID == mother)] & coords$y == new_y)) >1){
				# this_already_taken = intersect(all_ids,coords$ID[which(coords$x == coords$x[which(coords$ID == mother)] & coords$y == new_y)])
				# already_taken_siblings = unlist(strsplit(as.character(ped$SIBLINGS[which(coords$ID == this_already_taken)]),","))
				# already_taken_mother = unlist(strsplit(as.character(ped$MOTHER[which(coords$ID == this_already_taken)]),","))
				# already_taken_father = unlist(strsplit(as.character(ped$FATHER[which(coords$ID == this_already_taken)]),","))

				# coords$x[which(coords$ID == this_already_taken)] = NA
				# coords$y[which(coords$ID == this_already_taken)] = NA
				# coords$x[which(coords$ID == already_taken_mother)] = NA
				# coords$y[which(coords$ID == already_taken_mother)] = NA
				# coords$x[which(coords$ID == already_taken_father)] = NA
				# coords$y[which(coords$ID == already_taken_father)] = NA
				
				# ids = c(this_already_taken,already_taken_siblings,ids)	
		# }
		
		# if(length(which(coords$x == coords$x[which(coords$ID == father)] & coords$y == new_y)) >1){
				# this_already_taken = intersect(all_ids,coords$ID[which(coords$x == coords$x[which(coords$ID == father)] & coords$y == new_y)])
				# already_taken_siblings = unlist(strsplit(as.character(ped$SIBLINGS[which(coords$ID == this_already_taken)]),","))
				# already_taken_mother = unlist(strsplit(as.character(ped$MOTHER[which(coords$ID == this_already_taken)]),","))
				# already_taken_father = unlist(strsplit(as.character(ped$FATHER[which(coords$ID == this_already_taken)]),","))

				# coords$x[which(coords$ID == this_already_taken)] = NA
				# coords$y[which(coords$ID == this_already_taken)] = NA
				# coords$x[which(coords$ID == already_taken_mother)] = NA
				# coords$y[which(coords$ID == already_taken_mother)] = NA
				# coords$x[which(coords$ID == already_taken_father)] = NA
				# coords$y[which(coords$ID == already_taken_father)] = NA
				
				# ids = c(this_already_taken,already_taken_siblings,ids)	
		# }
		
		#check if another node is between the two parent nodes:
		coords_x_father = coords$x[which(coords$ID == father)]
		coords_x_mother = coords$x[which(coords$ID == mother)]
		coords_y_parents = coords$y[which(coords$ID == father)]

		if(length(coords$ID[which(coords$x>coords_x_father & coords$x<coords_x_mother & coords$y == coords_y_parents)]) >= 1){
			inbetween_id = coords$ID[which(coords$x>coords_x_father & coords$x<coords_x_mother & coords$y == coords_y_parents)]
			
			coords$x[which(coords$ID == inbetween_id)] = coords$x[which(coords$ID == inbetween_id)]-1
			while(length(which(coords$x == coords$x[which(coords$ID == inbetween_id)] & coords$y == coords_y_parents)) >1){
				coords$x[which(coords$ID == inbetween_id)] = coords$x[which(coords$ID == inbetween_id)]-1
				
				already_taken_mother = unlist(strsplit(as.character(ped$MOTHER[which(coords$ID == inbetween_id)]),","))
				already_taken_father = unlist(strsplit(as.character(ped$FATHER[which(coords$ID == inbetween_id)]),","))
				coords$x[which(coords$ID == already_taken_mother)] = NA
				coords$x[which(coords$ID == already_taken_father)] = NA
				
				ids = c(inbetween_id,ids)	
			}
		}
		
		if(length(coords$ID[which(coords$x>coords_x_mother & coords$x<coords_x_father & coords$y == coords_y_parents)]) >= 1){
			inbetween_id = coords$ID[which(coords$x>coords_x_father & coords$x<coords_x_mother & coords$y == coords_y_parents)]
			
			coords$x[which(coords$ID == inbetween_id)] = coords$x[which(coords$ID == inbetween_id)]-1
			while(length(which(coords$x == coords$x[which(coords$ID == inbetween_id)] & coords$y == coords_y_parents)) >1){
				coords$x[which(coords$ID == inbetween_id)] = coords$x[which(coords$ID == inbetween_id)]-1
				
				already_taken_mother = unlist(strsplit(as.character(ped$MOTHER[which(coords$ID == inbetween_id)]),","))
				already_taken_father = unlist(strsplit(as.character(ped$FATHER[which(coords$ID == inbetween_id)]),","))
				coords$x[which(coords$ID == already_taken_mother)] = NA
				coords$x[which(coords$ID == already_taken_father)] = NA
				
				ids = c(inbetween_id,ids)	
			}
		}

		#check if parents are too far away:
		if(coords_x_mother > coords$x[which(coords$ID == this_id)] & coords_x_father > coords$x[which(coords$ID == this_id)] & length(siblings) == 0){
			coords$x[which(coords$ID == this_id)] = coords$x[which(coords$ID == this_id)]+1
			for(i in 1:length(children)){
							coords$x[which(coords$ID == children[i])] = coords$x[which(coords$ID == children[i])] +1
			}
		}
		
		if(coords_x_mother < coords$x[which(coords$ID == this_id)] & coords_x_father < coords$x[which(coords$ID == this_id)] & length(siblings) == 0){
			coords$x[which(coords$ID == this_id)] = coords$x[which(coords$ID == this_id)]-1
			for(i in 1:length(children)){
							coords$x[which(coords$ID == children[i])] = coords$x[which(coords$ID == children[i])] -1
			}
		}
		
		#mark current id as already processed:
		ids[which(ids==this_id)] = "-"
	
		if(mother %in% ids){
				ids[which(ids==mother)] = "-"
				ids = c(mother,ids)	

		}
		if(father %in% ids){
				ids[which(ids==father)] = "-"
				ids = c(father,ids)	
		}
		ids=setdiff(ids,"-")

		}
		return(coords)
}